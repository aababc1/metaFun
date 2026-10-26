#!/bin/bash
# check_fastq_integrity.sh — verify dumped FASTQ pairs before QC (run on a compute node, not the login node).
#
# For every PAIRED run of a cohort in the manifest:
#   - both <run>_1.fastq.gz / <run>_2.fastq.gz exist and are non-empty
#   - gzip stream decodes to the end (gzip -dc exit 0; a job killed mid-compression leaves a truncated file)
#   - line count is a multiple of 4 and identical in _1 and _2 (pairs in sync)
#   - exact bases in _1 + _2 vs manifest base_count: < 97 % -> SHORT (incomplete dump; SRA may count a few
#     technical bases, hence the margin), > 103 % -> MORE_THAN_MANIFEST (manifest discrepancy, file kept)
# Output: <out-dir>/<cohort>.fastq_check.tsv (one row per run) and <out-dir>/<cohort>.bad_runs.txt
# (run ids that are not OK — feed it to make_samplesheet.py --exclude).
#
# Usage:
#   bash check_fastq_integrity.sh <cohort|all> [threads]
#   sbatch -p cpu32-only -N 1 --exclusive -t 4:00:00 -J fqcheck_PRJNA531273 -o fqcheck_%j.out \
#          --wrap "bash ~/metaFun/scripts/rawreadqc/check_fastq_integrity.sh PRJNA531273"
# Environment overrides: MANIFEST, FASTQ_DIR, OUT_DIR
set -uo pipefail

COHORT=${1:?cohort accession or 'all'}
THREADS=${2:-$(nproc)}
MANIFEST=${MANIFEST:-$HOME/metaCRC/manifest/run_manifest.tsv}
FASTQ_DIR=${FASTQ_DIR:-$HOME/metaCRC/fastq}
OUT_DIR=${OUT_DIR:-$HOME/metaCRC/manifest/fastq_check}
mkdir -p "$OUT_DIR"

check_run() {
    # args: run sample cohort base_count   -> prints one TSV row
    local run=$1 sample=$2 cohort=$3 bases=$4
    local f1="$FASTQ_DIR/$cohort/${run}_1.fastq.gz" f2="$FASTQ_DIR/$cohort/${run}_2.fastq.gz"
    local size1=0 size2=0 lines1=NA lines2=NA bases1=NA bases2=NA gz1=NA gz2=NA pairs=NA pct=NA verdict out
    if [ ! -s "$f1" ] || [ ! -s "$f2" ]; then
        verdict=MISSING
    else
        size1=$(stat -c %s "$f1"); size2=$(stat -c %s "$f2")
        # one pass per file: line count + exact base count (sequence lines) + gzip exit status (bash 4.2 compatible)
        out=$(gzip -dc "$f1" 2>/dev/null | awk 'NR%4==2{b+=length($0)} END{print NR, b+0}'; echo "${PIPESTATUS[0]}")
        lines1=$(echo "$out" | sed -n 1p | cut -d' ' -f1); bases1=$(echo "$out" | sed -n 1p | cut -d' ' -f2); gz1=$(echo "$out" | sed -n 2p)
        out=$(gzip -dc "$f2" 2>/dev/null | awk 'NR%4==2{b+=length($0)} END{print NR, b+0}'; echo "${PIPESTATUS[0]}")
        lines2=$(echo "$out" | sed -n 1p | cut -d' ' -f1); bases2=$(echo "$out" | sed -n 1p | cut -d' ' -f2); gz2=$(echo "$out" | sed -n 2p)
        if [ "$gz1" != 0 ] || [ "$gz2" != 0 ]; then
            verdict=TRUNCATED_GZ
        elif [ $((lines1 % 4)) -ne 0 ] || [ $((lines2 % 4)) -ne 0 ]; then
            verdict=BAD_RECORD_COUNT
        elif [ "$lines1" != "$lines2" ]; then
            verdict=PAIR_MISMATCH
        else
            pairs=$((lines1 / 4))
            if [ -n "$bases" ] && [ "$bases" != NA ] && [ "$bases" -gt 0 ] 2>/dev/null; then
                # exact bases in both files vs SRA base_count (SRA may count a few technical bases -> lower bound 97 %)
                pct=$(awk -v a="$bases1" -v b="$bases2" -v e="$bases" 'BEGIN{printf "%.2f", 100*(a+b)/e}')
                if awk -v p="$pct" 'BEGIN{exit !(p < 97)}'; then verdict=SHORT
                elif awk -v p="$pct" 'BEGIN{exit !(p > 103)}'; then verdict=MORE_THAN_MANIFEST
                else verdict=OK; fi
            else
                verdict=OK_NO_MANIFEST_BASES
            fi
        fi
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$run" "$sample" "$cohort" "$verdict" "$size1" "$size2" "$gz1" "$gz2" "$bases1" "$bases2" "$pairs" "$bases" "$pct" "$(date +%FT%T)"
}
export -f check_run
export FASTQ_DIR

# manifest columns (header): run_accession sample_accession cohort_accession library_layout base_count ...
hdr=$(head -n 1 "$MANIFEST")
col() { echo "$hdr" | tr '\t' '\n' | grep -n -x "$1" | cut -d: -f1; }
c_run=$(col run_accession); c_sam=$(col sample_accession); c_coh=$(col cohort_accession); c_lay=$(col library_layout); c_bas=$(col base_count)
[ -n "$c_run" ] && [ -n "$c_coh" ] && [ -n "$c_lay" ] || { echo "manifest lacks run_accession/cohort_accession/library_layout"; exit 1; }

if [ "$COHORT" = all ]; then
    cohorts=$(awk -F'\t' -v c="$c_coh" -v l="$c_lay" 'NR>1 && toupper($l)=="PAIRED"{print $c}' "$MANIFEST" | sort -u)
else
    cohorts=$COHORT
fi

for coh in $cohorts; do
    out="$OUT_DIR/$coh.fastq_check.tsv"
    echo "== $coh  ($(date))  threads=$THREADS"
    printf 'run\tsample\tcohort\tverdict\tsize_1\tsize_2\tgzip_rc_1\tgzip_rc_2\tbases_1\tbases_2\tpairs\tmanifest_bases\tpct_of_manifest\tchecked\n' > "$out.tmp"
    awk -F'\t' -v c="$c_coh" -v l="$c_lay" -v r="$c_run" -v s="$c_sam" -v b="$c_bas" -v coh="$coh" \
        'NR>1 && $c==coh && toupper($l)=="PAIRED"{print $r"\t"($s==""?$r:$s)"\t"$c"\t"(b==""?"NA":$b)}' "$MANIFEST" \
      | xargs -P "$THREADS" -L 1 bash -c 'check_run "$0" "$1" "$2" "$3"' >> "$out.tmp"
    { head -n 1 "$out.tmp"; tail -n +2 "$out.tmp" | sort; } > "$out" && rm -f "$out.tmp"
    awk -F'\t' 'NR>1 && $4!="OK" && $4!="OK_NO_MANIFEST_BASES" && $4!="MORE_THAN_MANIFEST"{print $1}' "$out" > "$OUT_DIR/$coh.bad_runs.txt"
    echo "-- $coh: $(tail -n +2 "$out" | wc -l) runs; verdicts:"; tail -n +2 "$out" | cut -f4 | sort | uniq -c | sed 's/^/   /'
    echo "   table: $out   bad runs: $OUT_DIR/$coh.bad_runs.txt ($(wc -l < "$OUT_DIR/$coh.bad_runs.txt"))"
done
