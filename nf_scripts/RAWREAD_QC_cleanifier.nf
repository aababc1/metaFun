#!/usr/bin/env nextflow
/*
 * RAWREAD_QC_cleanifier.nf
 * -----------------------------------------------------------------------------
 * Raw-read QC for whole-metagenome shotgun runs:
 *   fastp (adapter/quality trimming)  ->  Cleanifier (human-read removal, gapped k-mer index)
 *
 * Designed for an HPC node that is billed per node and cannot run containers
 * (e.g. KISTI K-BDS): one `nextflow run` per node with the local executor
 * (-c config/kbds.config -profile kbds). Compared with RAWREAD_QC_apptainer.nf:
 *   - input is a per-cohort samplesheet (CSV) instead of a directory glob
 *   - fastp and Cleanifier are connected by two named pipes, so trimmed reads
 *     are never gzip-compressed and re-read between the two tools
 *   - the final FASTQ pair is placed by the task itself (rename within the
 *     destination directory) and a `<run>.done` marker is written last, so an
 *     external monitor can rely on marker files instead of publishDir timing
 *   - no FastQC/MultiQC (fastp JSON/HTML per run + one cohort summary TSV)
 *
 * Written in Nextflow's strict syntax: requires Nextflow >= 25.04 (with
 * NXF_SYNTAX_PARSER=v2 on 25.x; default parser from 26.x). Tested with 26.10.3.
 *
 * Samplesheet columns (header required; base_count optional = SRA bases, used to detect truncated input):
 *   run_accession,sample_accession,cohort_accession,fastq_1,fastq_2[,base_count]
 *
 * Output layout (params.outdir = <QC root>):
 *   <outdir>/<cohort>/<run>_1.fastq.gz, <run>_2.fastq.gz     final reads
 *   <outdir>/<cohort>/<run>.done                               marker (written last)
 *   <outdir>/<cohort>/qc_logs/<run>.{fastp.json,fastp.html,fastp.log,cleanifier.log,qc.tsv}
 *   <outdir>/<cohort>/<cohort>.qc_summary.tsv                  one row per run
 */

// ---------------------------------------------------------------------------
// Parameters (defaults; -profile kbds overrides some of them, CLI overrides all)
// ---------------------------------------------------------------------------
params {
    samplesheet          = ''        // required
    outdir               = "${launchDir}/results/metagenome/RAWREAD_QC"
    cleanifier_index     = ''        // required: prefix of <index>.filter + <index>.info (e.g. ~/metaCRC/db/cleanifier/human)
    cleanifier_shared    = false     // load the index once per node into /dev/shm (cleanifier load / --shared)
    cleanifier_sensitive = false     // --sensitive: query every k-mer instead of sampling (slower)
    cleanifier_threshold = 0.5       // fraction of read bases covered by indexed k-mers -> read pair removed
    cpus                 = 4         // Nextflow cpus per task (slot accounting on the node)
    fastp_threads        = 2         // fastp -w (worker threads; fastp adds its own reader/writer threads)
    cleanifier_threads   = 3         // cleanifier --threads (= 1 reader + 1 writer + N-2 classifier threads)
    gz_threads           = 1         // cleanifier --compression-threads -> bgzip -@ per output file (2 files)
    gz_level             = 4         // cleanifier --compression-level -> bgzip -l (1 fastest .. 9)
    fastp_args           = '--detect_adapter_for_pe --trim_poly_g --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 --length_required 50 --n_base_limit 5 --dont_eval_duplication'
    task_memory          = '12 GB'   // informational for the local executor (Cleanifier filter index ~7 GB + buffers)
    max_retries          = 1         // a failing run is retried once, then ignored (workflow keeps going)
    skip_done            = true      // runs with an existing <run>.done are not re-processed
    bases_tolerance      = 0.05      // with a base_count column: fail a run whose raw bases are below (1 - tolerance) x base_count
    fifo_dir             = '/tmp'    // where the named pipes are created (node-local; falls back to the task dir)
    relay_buffer         = '64M'     // mbuffer RAM per mate between fastp and cleanifier (must exceed fastp's 4 MB output buffer)
    scripts_dir          = "${projectDir}/../scripts/rawreadqc"
    help                 = false
}

def helpMessage() {
    return """
    RAWREAD_QC_cleanifier  (fastp -> Cleanifier human-read removal)

    nextflow run RAWREAD_QC_cleanifier.nf -c config/kbds.config -profile kbds \\
        --samplesheet samplesheet.csv --outdir ~/metaCRC/qc \\
        --cleanifier_index ~/metaCRC/db/cleanifier/<index-prefix> [-resume]

    Required
      --samplesheet        CSV: run_accession,sample_accession,cohort_accession,fastq_1,fastq_2[,base_count]
      --cleanifier_index   index prefix (files <prefix>.filter and <prefix>.info must exist)
    Output
      --outdir             QC root; results go to <outdir>/<cohort_accession>/   [${params.outdir}]
    Threads per task (keep the sum near --cpus)
      --cpus ${params.cpus}  --fastp_threads ${params.fastp_threads}  --cleanifier_threads ${params.cleanifier_threads}  --gz_threads ${params.gz_threads}  --gz_level ${params.gz_level}
    Cleanifier
      --cleanifier_shared ${params.cleanifier_shared}  --cleanifier_sensitive ${params.cleanifier_sensitive}  --cleanifier_threshold ${params.cleanifier_threshold}
    fastp
      --fastp_args "${params.fastp_args}"
    Other
      --skip_done ${params.skip_done}  --max_retries ${params.max_retries}  --bases_tolerance ${params.bases_tolerance}  --fifo_dir ${params.fifo_dir}  --relay_buffer ${params.relay_buffer}  --task_memory '${params.task_memory}'
    """.stripIndent()
}

// absolute, normalised path as a string (Cleanifier strips leading '.' characters from its
// output prefix, so a relative prefix such as ./x would be written to /x)
def absPath(p) {
    return file(p.toString()).toAbsolutePath().normalize().toString()
}

def indexPrefix() {
    return absPath(params.cleanifier_index.toString().replaceAll(/\.(filter|info|hash)$/, ''))
}

/*
 * Load the Cleanifier index into POSIX shared memory once per node
 * (/dev/shm/<basename of index>), so that concurrent filter tasks share one copy
 * instead of each reading ~7 GB from the parallel file system.
 * `cleanifier load` (without --keep-running) returns once the object is filled; the
 * object outlives the loader. A small detached guard process removes it when the
 * batch job ends (SIGTERM from Slurm) and the workflow onComplete handler removes it as well.
 * (cleanifier 1.3.2's own --keep-running mode does not survive SIGTERM: it calls
 * sigwait() without blocking the signal first, so the default action kills it.)
 * `cache false`: must run on every (re)start, because the node changes.
 */
process load_index {
    tag "${index_prefix.tokenize('/')[-1]}"
    cache false
    cpus 1
    memory { params.task_memory }
    errorStrategy 'terminate'

    input:
    val index_prefix

    output:
    val index_prefix

    script:
    def name = index_prefix.tokenize('/')[-1]
    """
    set -euo pipefail
    if [ -e /dev/shm/${name} ]; then
        echo "shared memory object /dev/shm/${name} already present (\$(stat -c %s /dev/shm/${name}) bytes), reusing it"
    else
        cleanifier load --name '${index_prefix}' > cleanifier_load.log 2>&1
        cat cleanifier_load.log
        grep -qE "was created|already loaded" cleanifier_load.log || { echo "cleanifier load failed"; exit 1; }
        [ -e /dev/shm/${name} ] || { echo "/dev/shm/${name} missing after cleanifier load"; exit 1; }
    fi
    ls -l /dev/shm/${name}
    # guard: a detached shell that removes the shared object on SIGTERM/SIGINT (end of the Slurm job).
    # It closes every inherited descriptor > 2, otherwise Nextflow's output pipes stay open and this
    # task never completes. The guard writes its own PID to <launchDir>/.cleanifier_shm_guard.pid.
    pidfile='${launchDir}/.cleanifier_shm_guard.pid'
    if [ -s "\$pidfile" ] && kill -0 "\$(cat "\$pidfile")" 2>/dev/null; then
        echo "guard already running (pid \$(cat "\$pidfile"))"
    else
        rm -f "\$pidfile"
        (
            for fd in \$(seq 3 63); do eval "exec \$fd>&-"; done
            exec setsid bash -c 'echo \$\$ > "\$1"; trap "kill \$! 2>/dev/null; cleanifier remove --name \\"\$0\\" >/dev/null 2>&1; exit 0" TERM INT; while :; do sleep 60 & wait \$!; done' '${index_prefix}' "\$pidfile" > shm_guard.log 2>&1 < /dev/null
        ) &
        for i in \$(seq 1 20); do [ -s "\$pidfile" ] && break; sleep 0.5; done
        [ -s "\$pidfile" ] || { echo "shared-memory guard did not start"; cat shm_guard.log; exit 1; }
        echo "guard pid \$(cat "\$pidfile")"
    fi
    """
}

/*
 * One run (one FASTQ pair): fastp -> named pipes -> mbuffer relays -> named pipes -> cleanifier filter -> bgzip.
 * Everything that must be visible to the outside world is written by the task
 * itself, in this order:
 *   1. logs + per-run TSV copied to <dest>/qc_logs/
 *   2. temporary FASTQ (<dest>/.<run>.part_{1,2}_keep.fastq.gz) renamed to <dest>/<run>_{1,2}.fastq.gz
 *   3. <dest>/<run>.done written (contains the per-run TSV row)
 * Nothing is left in a half-written state on failure (temp files removed by the trap).
 */
process fastp_cleanifier {
    tag "${meta.run}"
    cpus { params.cpus }
    memory { params.task_memory }
    maxRetries { params.max_retries }
    errorStrategy { task.attempt <= params.max_retries ? 'retry' : 'ignore' }

    input:
    tuple val(meta), path(fq1), path(fq2)
    val index_prefix

    output:
    tuple val(meta), path("${meta.run}.qc.tsv"), emit: tsv
    path "${meta.run}.done",                     emit: done   // symlink to <dest>/<run>.done: if the marker is deleted, -resume re-runs the task
    path "${meta.run}.fastp.json",               emit: json
    path "${meta.run}.fastp.html"
    path "${meta.run}.fastp.log"
    path "${meta.run}.cleanifier.log"

    script:
    def run       = meta.run
    def shared    = params.cleanifier_shared    ? '--shared'    : ''
    def sensitive = params.cleanifier_sensitive ? '--sensitive' : ''
    """
    set -euo pipefail
    dest='${meta.dest}'
    part="\$dest/.${run}.part"
    tmp1="\${part}_1_keep.fastq.gz"
    tmp2="\${part}_2_keep.fastq.gz"
    fifo_dir=\$(mktemp -d -p '${params.fifo_dir}' "qc_${run}.XXXXXX" 2>/dev/null || mktemp -d "qc_${run}.XXXXXX")
    raw1="\$fifo_dir/${run}_1.raw"          # fastp -> mbuffer (no .gz suffix: fastp writes plain text)
    raw2="\$fifo_dir/${run}_2.raw"
    fifo1="\$fifo_dir/${run}_1.fastq"        # mbuffer -> cleanifier (extension must be .fastq = uncompressed for Cleanifier)
    fifo2="\$fifo_dir/${run}_2.fastq"
    cleanup() {
        rc=\$?
        if [ "\$rc" -ne 0 ]; then
            for p in .fastp.pid .cleanifier.pid .mbuf1.pid .mbuf2.pid; do [ -s \$p ] && kill "\$(cat \$p)" 2>/dev/null || true; done
            rm -f "\$tmp1" "\$tmp2"
        fi
        rm -rf "\$fifo_dir" .fastp.pid .cleanifier.pid .mbuf1.pid .mbuf2.pid .fastp.rc .cleanifier.rc .mbuf1.rc .mbuf2.rc
    }
    trap cleanup EXIT

    mkdir -p "\$dest/qc_logs"
    rm -f "\$tmp1" "\$tmp2" .fastp.pid .cleanifier.pid .mbuf1.pid .mbuf2.pid .fastp.rc .cleanifier.rc .mbuf1.rc .mbuf2.rc
    mkfifo "\$raw1" "\$raw2" "\$fifo1" "\$fifo2"
    t0=\$(date +%s)

    # fastp -> [raw FIFO] -> mbuffer (${params.relay_buffer} RAM) -> [FIFO] -> cleanifier, one relay per mate.
    # The relays decouple the two streams: Cleanifier's reader fills ~50 kB from mate 1, then ~50 kB from
    # mate 2, whereas fastp (>= 0.23) buffers 4 MB per output and flushes mate 1 completely before mate 2
    # when it finishes -- without the relays the two block each other on a 64 kB kernel pipe (verified
    # deadlock with fastp 1.4). mbuffer absorbs the 4 MB flush, so any open/write order works.
    # Each tool runs in a subshell that records its PID and exit status in a file, so the task can react
    # as soon as one of them fails (bash 4.2 on CentOS 7 has no `wait -n`); partners left blocked on a
    # pipe are killed instead of hanging the task.
    (
        mbuffer -q -f -m ${params.relay_buffer} -s 1M -i "\$raw1" -o "\$fifo1" 2> .mbuf1.log &
        echo \$! > .mbuf1.pid; wait \$! && echo 0 > .mbuf1.rc || echo \$? > .mbuf1.rc
    ) &
    (
        mbuffer -q -f -m ${params.relay_buffer} -s 1M -i "\$raw2" -o "\$fifo2" 2> .mbuf2.log &
        echo \$! > .mbuf2.pid; wait \$! && echo 0 > .mbuf2.rc || echo \$? > .mbuf2.rc
    ) &
    (
        fastp --in1 ${fq1} --in2 ${fq2} --out1 "\$raw1" --out2 "\$raw2" \\
            --thread ${params.fastp_threads} ${params.fastp_args} \\
            --json ${run}.fastp.json --html ${run}.fastp.html --report_title '${run}' \\
            2> ${run}.fastp.log &
        echo \$! > .fastp.pid; wait \$! && echo 0 > .fastp.rc || echo \$? > .fastp.rc
    ) &
    (
        cleanifier filter --index '${index_prefix}' ${shared} ${sensitive} \\
            --fastq "\$fifo1" --pairs "\$fifo2" --out "\$part" \\
            --threads ${params.cleanifier_threads} --threshold ${params.cleanifier_threshold} \\
            --compression gz --compression-threads ${params.gz_threads} --compression-level ${params.gz_level} \\
            > ${run}.cleanifier.log 2>&1 &
        echo \$! > .cleanifier.pid; wait \$! && echo 0 > .cleanifier.rc || echo \$? > .cleanifier.rc
    ) &
    failed=""
    while [ ! -s .fastp.rc ] || [ ! -s .cleanifier.rc ] || [ ! -s .mbuf1.rc ] || [ ! -s .mbuf2.rc ]; do
        for t in fastp cleanifier mbuf1 mbuf2; do
            if [ -s .\$t.rc ] && [ "\$(cat .\$t.rc)" != 0 ]; then failed="\$failed \$t=\$(cat .\$t.rc)"; fi
        done
        [ -n "\$failed" ] && break
        sleep 1
    done
    for t in fastp cleanifier mbuf1 mbuf2; do
        if [ -s .\$t.rc ] && [ "\$(cat .\$t.rc)" != 0 ]; then case "\$failed" in *"\$t="*) ;; *) failed="\$failed \$t=\$(cat .\$t.rc)";; esac; fi
    done
    if [ -n "\$failed" ]; then
        echo "ERROR: tool failure:\$failed (see ${run}.fastp.log / ${run}.cleanifier.log; mbuffer: .mbuf*.log)" >&2
        cat .mbuf1.log .mbuf2.log 2>/dev/null >&2 || true
        exit 1
    fi
    wait
    t1=\$(date +%s)

    # per-run TSV + integrity checks (fastp pair count == keep+filter, BGZF EOF block present)
    python3 '${params.scripts_dir}/run_summary.py' \\
        --run '${run}' --sample '${meta.sample}' --cohort '${meta.cohort}' \\
        --fastp-json ${run}.fastp.json --cleanifier-log ${run}.cleanifier.log \\
        --fastq "\$tmp1" "\$tmp2" --wall-sec \$((t1 - t0)) --cpus ${task.cpus} --attempt ${task.attempt} \\
        --expected-bases ${meta.bases} --bases-tolerance ${params.bases_tolerance} \\
        --cleanifier-version "\$(cleanifier --version 2>/dev/null | tail -n 1)" \\
        > ${run}.qc.tsv

    # final placement: logs first, then the two FASTQ files (rename inside the same directory), marker last
    cp -f ${run}.fastp.json ${run}.fastp.html ${run}.fastp.log ${run}.cleanifier.log ${run}.qc.tsv "\$dest/qc_logs/"
    mv -f "\$tmp1" "\$dest/${run}_1.fastq.gz"
    mv -f "\$tmp2" "\$dest/${run}_2.fastq.gz"
    tail -n 1 ${run}.qc.tsv > "\$dest/${run}.done"
    ln -sfn "\$dest/${run}.done" ${run}.done
    """
}

/*
 * Cohort summary: merge every <dest>/qc_logs/*.qc.tsv (all chunks, all attempts so far)
 * into <dest>/<cohort>.qc_summary.tsv (written atomically). Runs after all tasks of this
 * workflow run have finished; the staged files only establish that dependency.
 */
process qc_summary {
    tag "${cohort}"
    cpus 1
    memory '1 GB'
    errorStrategy 'ignore'

    input:
    tuple val(cohort), val(dest), path(tsvs, stageAs: 'staged/*')

    output:
    path "${cohort}.qc_summary.tsv"

    script:
    """
    set -euo pipefail
    python3 '${params.scripts_dir}/qc_summary.py' --qc-logs '${dest}/qc_logs' --cohort '${cohort}' \\
        --out '${dest}/${cohort}.qc_summary.tsv'
    cp -f '${dest}/${cohort}.qc_summary.tsv' ${cohort}.qc_summary.tsv
    """
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------
workflow {
    main:
    if (params.help) {
        log.info helpMessage()
        exit 0
    }
    if (!params.samplesheet)      error "Missing --samplesheet (CSV with run_accession,sample_accession,cohort_accession,fastq_1,fastq_2)"
    if (!params.cleanifier_index) error "Missing --cleanifier_index (prefix of <index>.filter / <index>.info)"

    def outdir_abs = absPath(params.outdir)
    def index_abs  = indexPrefix()
    if (!file("${index_abs}.filter").exists() && !file("${index_abs}.hash").exists())
        error "Cleanifier index not found: ${index_abs}.filter (or .hash)"
    if (!file("${index_abs}.info").exists())
        error "Cleanifier index info not found: ${index_abs}.info"
    def required_cols = ['run_accession', 'sample_accession', 'cohort_accession', 'fastq_1', 'fastq_2']

    log.info """
    =============================================================
      RAWREAD_QC_cleanifier   fastp -> Cleanifier (human removal)
    =============================================================
    samplesheet       : ${params.samplesheet}
    outdir            : ${outdir_abs}
    cleanifier index  : ${index_abs}  (shared memory: ${params.cleanifier_shared}, sensitive: ${params.cleanifier_sensitive}, threshold: ${params.cleanifier_threshold})
    threads per task  : cpus=${params.cpus}  fastp=${params.fastp_threads}  cleanifier=${params.cleanifier_threads}  bgzip=${params.gz_threads}x2 (level ${params.gz_level})
    fastp args        : ${params.fastp_args}
    retries / skip    : max_retries=${params.max_retries}  skip_done=${params.skip_done}
    workDir           : ${workflow.workDir}
    =============================================================
    """.stripIndent()

    runs_ch = Channel.fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row ->
            def missing = required_cols.findAll { !row.containsKey(it) }
            if (missing) error "Samplesheet ${params.samplesheet} lacks column(s): ${missing.join(', ')} (required: ${required_cols.join(',')})"
            if (!row.run_accession || !row.cohort_accession) error "Samplesheet row without run_accession/cohort_accession: ${row}"
            def meta = [
                run    : row.run_accession,
                sample : row.sample_accession ?: row.run_accession,
                cohort : row.cohort_accession,
                bases  : (row.base_count ?: '0').toString(),
                dest   : "${outdir_abs}/${row.cohort_accession}".toString()
            ]
            tuple(meta, file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true))
        }

    // runs whose marker already exists are not processed again (unless --skip_done false)
    branched = runs_ch.branch { meta, fq1, fq2 ->
        done: params.skip_done && file("${meta.dest}/${meta.run}.done").exists()
        todo: true
    }
    branched.done.count().subscribe { n -> if (n > 0) log.info "skip_done: ${n} run(s) already have a .done marker and are skipped" }
    branched.todo.count().subscribe { n -> log.info "runs to process: ${n}" }

    index_ch = Channel.value(index_abs)
    if (params.cleanifier_shared) {
        load_index(index_ch)
        ready_ch = load_index.out
    } else {
        ready_ch = index_ch
    }

    fastp_cleanifier(branched.todo, ready_ch)

    // one summary per cohort present in the samplesheet, after every run of this invocation finished
    cohorts_ch = runs_ch.map { meta, fq1, fq2 -> tuple(meta.cohort, meta.dest) }.unique()
    tsv_list   = fastp_cleanifier.out.tsv.map { meta, tsv -> tsv }.collect().ifEmpty([]).map { [it] }
    qc_summary(cohorts_ch.combine(tsv_list))

    onComplete:
    if (params.cleanifier_shared && params.cleanifier_index) {
        // remove the shared-memory index: ask the guard (SIGTERM -> cleanifier remove), then check ourselves.
        // The sbatch wrapper should also run `cleanifier remove --name <index>` after nextflow exits.
        def idx  = indexPrefix()
        def name = idx.tokenize('/')[-1]
        def pidfile = "${launchDir}/.cleanifier_shm_guard.pid"
        def cmd = """
            if [ -s '${pidfile}' ]; then kill -TERM "\$(cat '${pidfile}')" 2>/dev/null || true; fi
            for i in 1 2 3 4 5; do [ -e /dev/shm/${name} ] || break; sleep 1; done
            [ -e /dev/shm/${name} ] && cleanifier remove --name '${idx}' || true
            rm -f '${pidfile}'
        """
        try { ["bash", "-c", cmd].execute().waitFor() } catch (Exception e) { log.warn "shared-memory cleanup: ${e.message}" }
    }
    log.info "RAWREAD_QC_cleanifier ${workflow.success ? 'finished' : 'FAILED'} after ${workflow.duration} (exit ${workflow.exitStatus})"
}
