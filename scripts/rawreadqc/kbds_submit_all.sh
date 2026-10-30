#!/bin/bash
# kbds_submit_all.sh — submit RAWREAD_QC_cleanifier for every cohort that is ready, one node job per cohort.
#
# Idempotent: run it as often as you like (e.g. once an hour from tmux). Each pass:
#   1. reads the manifest (PAIRED runs per cohort, bases per cohort)
#   2. skips cohorts that are finished (every run has qc/<cohort>/<run>.done) or that have nothing new to do
#      (every dumped run already has its .done; the rest is not dumped yet)
#   3. skips cohorts that already have a QC job pending/running (squeue name qc_<cohort>)
#   4. skips cohorts whose dump is incomplete (fewer FASTQ pairs than PAIRED runs), unless ALLOW_PARTIAL=1
#      (the pipeline's skip_done + -resume make a later re-submission pick up the stragglers)
#   5. builds runs/qc/<cohort>/samplesheet.csv (excluding manifest/fastq_check/<cohort>.bad_runs.txt if present)
#   6. picks a partition with a free slot (cpu32-only first, then cpu64-only; limits = node counts)
#   7. submits kbds_test_qc.sbatch with -J qc_<cohort>, -t from the cohort size (7.9 s/Gb wall on cpu32, x2 margin);
#      a cohort whose estimate exceeds CHUNK_H (24 h) is split into k balanced chunks -> runs/qc/<cohort>.<k>/,
#      one job each (-J qc_<cohort>.<k>), all writing into qc/<cohort>/
# Dry run by default: prints what it would do. SUBMIT=1 actually submits.
#
# Usage:  ~/metaFun/scripts/rawreadqc/kbds_submit_all.sh                        # dry run
#         ALLOW_PARTIAL=1 SUBMIT=1 ~/metaFun/scripts/rawreadqc/kbds_submit_all.sh  # submit
# (run it directly, not via `bash file`: on K-BDS `bash` is a login-shell function in the user's rc)
# Environment overrides: MANIFEST FASTQ_DIR QC_ROOT RUNS_DIR INDEX METAFUN ALLOW_PARTIAL=1 MAX_SUBMIT=n
#                        PARTITIONS="cpu32-only:5 cpu64-only:5" ONLY="PRJEB6070 PRJNA447983" EXTRA_SBATCH="-A acct"
set -uo pipefail

MANIFEST=${MANIFEST:-$HOME/metaCRC/manifest/run_manifest.tsv}
FASTQ_DIR=${FASTQ_DIR:-$HOME/metaCRC/fastq}
QC_ROOT=${QC_ROOT:-$HOME/metaCRC/qc}
RUNS_DIR=${RUNS_DIR:-$HOME/metaCRC/runs/qc}
CHECK_DIR=${CHECK_DIR:-$HOME/metaCRC/manifest/fastq_check}
METAFUN=${METAFUN:-$HOME/metaFun}
INDEX=${INDEX:-$(find $HOME/metaCRC/db/cleanifier -name '*.info' 2>/dev/null | head -n 1 | sed 's/\.info$//')}
PARTITIONS=${PARTITIONS:-"cpu32-only:5 cpu64-only:5"}     # partition:max concurrent qc jobs (= node count)
SEC_PER_GB=${SEC_PER_GB:-7.9}                               # node wall seconds per Gbase on cpu32 (measured 2026-10-27)
MARGIN=${MARGIN:-2}                                         # -t = estimate x MARGIN + 1 h, capped at 120 h
CHUNK_H=${CHUNK_H:-24}                                      # split a cohort into ceil(estimate / CHUNK_H) chunks
MAX_SUBMIT=${MAX_SUBMIT:-100}
ALLOW_PARTIAL=${ALLOW_PARTIAL:-0}
SUBMIT=${SUBMIT:-0}
ONLY=${ONLY:-}
EXTRA_SBATCH=${EXTRA_SBATCH:-}
WRAPPER=$METAFUN/scripts/rawreadqc/kbds_test_qc.sbatch

[ -s "$MANIFEST" ] || { echo "manifest not found: $MANIFEST"; exit 1; }
[ -n "$INDEX" ] && [ -f "$INDEX.info" ] || { echo "Cleanifier index not found (INDEX=$INDEX)"; exit 1; }
[ -f "$WRAPPER" ] || { echo "wrapper not found: $WRAPPER"; exit 1; }

hdr=$(head -n 1 "$MANIFEST")
col() { echo "$hdr" | tr '\t' '\n' | grep -n -x "$1" | cut -d: -f1; }
c_run=$(col run_accession); c_coh=$(col cohort_accession); c_lay=$(col library_layout); c_bas=$(col base_count)
[ -n "$c_run" ] && [ -n "$c_coh" ] && [ -n "$c_lay" ] || { echo "manifest lacks run_accession/cohort_accession/library_layout"; exit 1; }

# current qc jobs per partition (pending + running)
declare -A running
while read -r part name; do
    [ -n "$part" ] || continue
    running[$part]=$(( ${running[$part]:-0} + 1 ))
    running[name:$name]=1
done < <(squeue -u "${USER:-$(id -un)}" -h -o "%P %j" 2>/dev/null | awk '$2 ~ /^qc_/')

declare -A cap
for pc in $PARTITIONS; do cap[${pc%%:*}]=${pc##*:}; done
free_partition() {
    for pc in $PARTITIONS; do
        local p=${pc%%:*}
        if [ "${running[$p]:-0}" -lt "${cap[$p]}" ]; then echo "$p"; return 0; fi
    done
    return 1
}

printf "%-14s %6s %6s %6s %6s %8s %-9s %s\n" cohort runs dumped done todo Tb action detail
n_sub=0
while read -r coh nruns bases; do
    [ -n "$coh" ] || continue
    if [ -n "$ONLY" ]; then case " $ONLY " in *" $coh "*) ;; *) continue;; esac; fi
    tb=$(awk -v b="$bases" 'BEGIN{printf "%.2f", b/1e12}')
    runs=$(awk -F'\t' -v c="$c_coh" -v l="$c_lay" -v r="$c_run" -v coh="$coh" 'NR>1 && $c==coh && toupper($l)=="PAIRED"{print $r}' "$MANIFEST")
    dumped=0; done_n=0; todo=0
    for r in $runs; do
        d=0; [ -s "$FASTQ_DIR/$coh/${r}_1.fastq.gz" ] && [ -s "$FASTQ_DIR/$coh/${r}_2.fastq.gz" ] && { dumped=$((dumped+1)); d=1; }
        if [ -e "$QC_ROOT/$coh/$r.done" ]; then done_n=$((done_n+1)); elif [ "$d" = 1 ]; then todo=$((todo+1)); fi
    done
    act=""; detail=""
    if [ "$done_n" -ge "$nruns" ]; then
        act=finished
    elif [ -n "${running[name:qc_$coh]:-}" ] || [ -n "${running[name:qc_$coh.1]:-}" ]; then
        act=queued; detail="qc_$coh already pending/running"
    elif [ "$dumped" -eq 0 ]; then
        act=wait; detail="no FASTQ yet"
    elif [ "$todo" -eq 0 ]; then
        act=wait; detail="every dumped run is done; $((nruns-dumped)) run(s) not dumped yet"
    elif [ "$dumped" -lt "$nruns" ] && [ "$ALLOW_PARTIAL" != 1 ]; then
        act=wait; detail="dump $dumped/$nruns (ALLOW_PARTIAL=1 to submit now)"
    elif [ "$n_sub" -ge "$MAX_SUBMIT" ]; then
        act=later; detail="MAX_SUBMIT reached"
    else
        part=$(free_partition) || { act=later; detail="no free partition slot (${PARTITIONS})"; }
    fi
    if [ -z "$act" ]; then
        est_h=$(awk -v b="$bases" -v s="$SEC_PER_GB" 'BEGIN{printf "%.1f", b/1e9*s/3600}')
        nchunk=$(awk -v e="$est_h" -v c="$CHUNK_H" 'BEGIN{k=int(e/c); if(k*c<e)k++; print (k<1)?1:k}')
        excl=""; [ -s "$CHECK_DIR/$coh.bad_runs.txt" ] && excl="--exclude $CHECK_DIR/$coh.bad_runs.txt"
        mkdir -p "$RUNS_DIR/$coh"
        if ! python3 "$METAFUN/scripts/rawreadqc/make_samplesheet.py" --manifest "$MANIFEST" --fastq-dir "$FASTQ_DIR" \
                --cohort "$coh" $excl --chunks "$nchunk" --out "$RUNS_DIR/$coh/samplesheet.csv" > "$RUNS_DIR/$coh/samplesheet.log" 2>&1; then
            act=error; detail="make_samplesheet failed (see $RUNS_DIR/$coh/samplesheet.log)"
        else
            # chunk k -> its own run dir runs/qc/<coh>.<k> (k = 1..nchunk); unchunked -> runs/qc/<coh>
            parts=""
            for k in $(seq 1 "$nchunk"); do
                if [ "$nchunk" -eq 1 ]; then rd="$RUNS_DIR/$coh"; tag="$coh"; wd="$HOME/metaCRC/work/qc/$coh"
                else
                    rd="$RUNS_DIR/$coh.$k"; tag="$coh.$k"; wd="$HOME/metaCRC/work/qc/$coh.$k"; mkdir -p "$rd"
                    mv -f "$RUNS_DIR/$coh/samplesheet.$k.csv" "$rd/samplesheet.csv"
                fi
                [ -n "${running[name:qc_$tag]:-}" ] && { parts="$parts $tag:queued"; continue; }
                part=$(free_partition) || { parts="$parts $tag:no-slot"; continue; }
                nss=$(($(wc -l < "$rd/samplesheet.csv") - 1))
                [ "$nss" -gt 0 ] || { parts="$parts $tag:empty"; continue; }
                wall_h=$(awk -v e="$est_h" -v n="$nchunk" -v m="$MARGIN" 'BEGIN{h=e/n*m+1; if(h>120)h=120; printf "%d", (h<1)?1:h}')
                cmd="sbatch -J qc_$tag -p $part -t ${wall_h}:00:00 $EXTRA_SBATCH $WRAPPER"
                if [ "$SUBMIT" = 1 ]; then
                    out=$(cd "$rd" && COHORT=$coh RUN_DIR=$rd WORK_DIR=$wd INDEX=$INDEX QC_ROOT=$QC_ROOT METAFUN=$METAFUN $cmd 2>&1)
                    if echo "$out" | grep -q "Submitted batch job"; then
                        parts="$parts $tag:$part/${wall_h}h/${nss}runs/job${out##* }"; running[$part]=$(( ${running[$part]:-0} + 1 )); n_sub=$((n_sub+1))
                    else
                        parts="$parts $tag:ERROR($out)"
                    fi
                else
                    parts="$parts $tag:$part/${wall_h}h/${nss}runs"; running[$part]=$(( ${running[$part]:-0} + 1 )); n_sub=$((n_sub+1))
                fi
            done
            act=$([ "$SUBMIT" = 1 ] && echo submitted || echo would-run); detail="est ${est_h}h, ${nchunk} job(s):$parts"
        fi
    fi
    printf "%-14s %6d %6d %6d %6d %8s %-9s %s\n" "$coh" "$nruns" "$dumped" "$done_n" "$todo" "$tb" "$act" "$detail"
done < <(awk -F'\t' -v c="$c_coh" -v l="$c_lay" -v b="$c_bas" 'NR>1 && toupper($l)=="PAIRED"{n[$c]++; s[$c]+=(b?$b:0)} END{for(k in n) print k, n[k], s[k]}' "$MANIFEST" | sort -k3 -n -r)
echo "qc jobs now: $(for pc in $PARTITIONS; do p=${pc%%:*}; printf "%s=%s/%s " "$p" "${running[$p]:-0}" "${cap[$p]}"; done)"
