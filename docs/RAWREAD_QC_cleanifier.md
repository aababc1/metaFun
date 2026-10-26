# RAWREAD_QC_cleanifier — fastp + Cleanifier host removal (samplesheet, marker files, K-BDS profile)

`nf_scripts/RAWREAD_QC_cleanifier.nf` is a variant of `RAWREAD_QC` for large run collections on
HPC nodes that are billed per node and may not run containers. It replaces Bowtie2 host removal
with [Cleanifier](https://gitlab.com/rahmannlab/cleanifier) (Zentgraf, Schmitz & Rahmann,
*Bioinformatics* 2025, doi:10.1093/bioinformatics/btaf632): a gapped k-mer (k=29, w=33) Cuckoo-filter
index of the T2T reference + 1000G variants (AF ≥ 0.01) + IPD-IMGT/HLA + 47 HPRC assemblies (< 7 GB).
A read pair is removed when either mate has ≥ `--threshold` (default 0.5) of its bases covered by
indexed k-mers; by default k-mers are *sampled* (every k/2 positions), `--cleanifier_sensitive`
queries all of them.

Requires Nextflow >= 26 (the script and `config/kbds.config` use the strict syntax, including the
`onComplete:` workflow section); verified with Nextflow 26.10.3, fastp 1.4.0, Cleanifier 1.3.2,
htslib 1.24, mbuffer 20160228 (the `scripts/rawreadqc/kbds_env.yml` environment).

## Pipeline per run

```
<run>_1.fastq.gz ─┐            ┌─ <run>_1.raw (FIFO) ─ mbuffer 64M ─ <run>_1.fastq (FIFO) ─┐
                  ├─ fastp ────┤                                                           ├─ cleanifier filter ─ bgzip ─► <dest>/.<run>.part_{1,2}_keep.fastq.gz
<run>_2.fastq.gz ─┘  (-w 2)   └─ <run>_2.raw (FIFO) ─ mbuffer 64M ─ <run>_2.fastq (FIFO) ─┘   (--threads 3, 1 classifier)   (-@ 1 -l 4, one per file)
```

The trimmed reads never touch the disk between the two tools (no gzip round trip). The two
`mbuffer` relays (RAM buffers, 64 MB each, `--relay_buffer`) are required: Cleanifier's reader
alternates between the mates in ~50 kB quanta, while fastp >= 0.23 buffers 4 MB per output and
flushes mate 1 completely before mate 2 when it finishes — connected directly through 64 kB
kernel pipes the two block each other (reproducible deadlock with fastp 1.4.0; fastp 0.22 wrote
both mates incrementally and did not need the relays). Then, inside the same task and in this order:

1. `run_summary.py` writes `<run>.qc.tsv` and checks: fastp pairs after filtering == Cleanifier
   `keep + filter`; both outputs end with the BGZF EOF block; with a `base_count` column, raw bases
   seen by fastp ≥ 95 % of the manifest value (fastp accepts a truncated gzip stream silently).
2. logs are copied to `<dest>/qc_logs/` (`<run>.fastp.json/.html/.log`, `<run>.cleanifier.log`, `<run>.qc.tsv`)
3. the two FASTQ files are renamed to `<dest>/<run>_1.fastq.gz`, `<run>_2.fastq.gz` (same directory → atomic)
4. `<dest>/<run>.done` is written last (one TSV row: counts, host fraction, wall time, versions)

A failing run is retried once, then ignored (`trace.txt` keeps the FAILED rows; temporary files are
removed, no marker is written). `<dest>` = `--outdir/<cohort_accession>`. Nothing uses `publishDir`.

After all runs of an invocation, `qc_summary` merges every `<dest>/qc_logs/*.qc.tsv` into
`<dest>/<cohort>.qc_summary.tsv` (so chunks and re-runs of the same cohort add up).

## Samplesheet

CSV with header `run_accession,sample_accession,cohort_accession,fastq_1,fastq_2[,base_count]`.
`scripts/rawreadqc/make_samplesheet.py` builds it from a run manifest and can split a cohort into
`samplesheet.<k>.csv` chunks balanced by `base_count` (`--chunks N` or `--max-bases B`). Runs whose
`<run>.done` already exists are skipped (`--skip_done false` to force).

## Threads per task (`--cpus 4`)

| tool | option | threads | note |
|---|---|---|---|
| fastp | `--fastp_threads 2` | 2 workers (+ reader/writer threads) | usually the bottleneck |
| cleanifier | `--cleanifier_threads 3` | 1 reader + 1 writer + 1 classifier | classification is cheap |
| bgzip × 2 | `--gz_threads 1`, `--gz_level 4` | 1 per output file | level 1 ≈ 3× faster than 6, ~15 % larger |

`trace.txt` column `%cpu` shows what a task really used; adjust `--fastp_threads` first.

## fastp parameters (`--fastp_args`, defaults and why)

| option | reason |
|---|---|
| `--detect_adapter_for_pe` | 25+ public cohorts with unknown library preps: sequence-based adapter detection in addition to overlap analysis |
| `--trim_poly_g` | forced on: SRA-dumped read names lose the instrument tag fastp uses to auto-enable polyG trimming for NovaSeq/NextSeq data |
| `--cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20` | sliding-window end trimming (Trimmomatic `SLIDINGWINDOW:4:20` equivalent); removes low-quality tails that create spurious k-mers in assembly |
| `--length_required 50` | reads < 50 bp are useless for MEGAHIT (k ≥ 21) and k = 31 profilers; default 15 keeps junk |
| `--n_base_limit 5` | fastp default, kept explicit |
| `--dont_eval_duplication` | duplication estimate costs time/memory and is not used downstream |
| (not set) `--dedup`, `--correction`, `--low_complexity_filter` | not standard for metagenome assembly; dedup would remove real high-abundance reads |

## K-BDS profile (`-c config/kbds.config -profile kbds`)

`config/kbds.config` is a standalone, strict-syntax config (the shared `config/nextflow.config` is
legacy syntax and does not parse under Nextflow >= 26; it is untouched). One `sbatch` = one exclusive
node = one `nextflow run`, local executor whose capacity defaults to nproc, `process.cpus = 4`
(8 concurrent runs on cpu32-only, 16 on cpu64). No containers; tools come from the conda environment
`scripts/rawreadqc/kbds_env.yml` (activate it before `nextflow run`, or pass `--conda_env <name>` so
every task activates it). Fixed in the profile:

```
trace  { enabled = true; overwrite = true; raw = true; file = "${params.run_dir}/trace.txt"
         fields = 'task_id,hash,native_id,process,tag,name,status,exit,attempt,cpus,memory,submit,start,complete,duration,realtime,%cpu,peak_rss,rchar,wchar,workdir' }
report { enabled = true; overwrite = true; file = "${params.run_dir}/report.html" }
workDir = params.work_dir
```

`params.run_dir` defaults to the launch directory, `params.work_dir` to `<launch>/work`.
`params.cleanifier_shared = true` in this profile: `load_index` loads the index once per node into
`/dev/shm/<index name>` (`cleanifier load`) and every task uses `--shared`; a detached guard process
removes the object when the job ends (SIGTERM) and `workflow.onComplete` removes it too. Set
`--cleanifier_shared false` to let each task read the index file instead (~7 GB RSS per task).

```
cd ~/metaCRC/runs/qc/<BioProject>[.<k>]
nextflow run ~/metaFun/nf_scripts/RAWREAD_QC_cleanifier.nf -c ~/metaFun/config/kbds.config -profile kbds \
    --samplesheet $PWD/samplesheet.csv --outdir ~/metaCRC/qc \
    --run_dir $PWD --work_dir ~/metaCRC/work/qc/<BioProject>[.<k>] \
    --cleanifier_index ~/metaCRC/db/cleanifier/<prefix> -resume
```

Process names in `trace.txt`: `load_index` (once per run, only with shared memory), `fastp_cleanifier`
(one per run, `tag` = run accession), `qc_summary` (one per cohort).

## Checking the dumped FASTQ first

`scripts/rawreadqc/check_fastq_integrity.sh <cohort|all>` (run it on a compute node, e.g.
`sbatch -p cpu32-only -N 1 --exclusive -t 4:00:00 --wrap "bash .../check_fastq_integrity.sh PRJNA531273"`)
decompresses every pair once and reports per run: gzip stream complete (a job killed during
compression leaves a truncated file), record count a multiple of 4, `_1`/`_2` in sync, and
pairs ≈ manifest `base_count / (len1 + len2)` (≥ 97 %). It writes `<cohort>.fastq_check.tsv` and
`<cohort>.bad_runs.txt` under `~/metaCRC/manifest/fastq_check/`; pass the latter to
`make_samplesheet.py --exclude`. Inside the pipeline, fastp 1.4 itself fails on a truncated gzip
stream and `run_summary.py` additionally compares raw bases with the manifest (95 %).

## Fallbacks if Cleanifier does not run on glibc 2.17

(Not needed on K-BDS: verified 2026-10-26 — glibc 2.17, Python 3.15.2, numba 0.68 JIT OK.) Cleanifier is pure Python ≥ 3.12 + numba (+ `bgzip` from htslib); conda-forge builds of python,
numba and llvmlite target glibc ≥ 2.17, so a solver error mentioning `__glibc >=2.28` would point at
a single package that can be pinned to an older build. If it still fails: (a) **Deacon**
(Rust single binary, minimizer-based, also no alignment) with the same pipe layout (`deacon filter`
reads FASTQ pairs; output naming differs, adjust the process); (b) **Apptainer** image with the
conda env if `apptainer` exists on compute nodes; (c) last resort, the Bowtie2 process of
`RAWREAD_QC_apptainer.nf` (slow: alignment).

## Notes / known behaviour

- Cleanifier 1.3.2 strips leading `.` characters from `--out`: a relative prefix such as `./x`
  is written to `/x`. The workflow always passes an absolute prefix.
- `cleanifier load --keep-running` does not survive SIGTERM (sigwait without blocking the signal),
  hence the separate guard process.
- `NUMBA_CACHE_DIR` is set in the profile but Cleanifier 1.3.2 does not use numba's on-disk cache
  (`cache=True` is never passed), so each invocation pays the JIT compilation (≈ 5–15 s).
- `report.html` is overwritten on `-resume`; `trace.txt` too (`overwrite = true`), with CACHED rows.
- Under the strict parser, process directives that read `params` must be closures
  (`cpus { params.cpus }`); a plain `cpus params.cpus` evaluated to null and logged
  "Access to undefined parameter". No `memory` directive is set: with the shared-memory index a
  task's private RSS is small (fastp buffers + mbuffer 2 x 64 MB) and the local executor would only
  use it as a scheduling constraint.
- The shared-memory object `/dev/shm/<index>` is removed by the workflow's `onComplete`, by the
  guard process on SIGTERM, and (belt and braces) by the sbatch wrapper after `nextflow` exits:
  `cleanifier remove --name <index prefix>`.
