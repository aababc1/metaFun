#!/bin/bash
# kbds_login_check.sh — smoke test for RAWREAD_QC_cleanifier on the K-BDS *login node*
# (seconds of CPU, a few MB of data; no real analysis). Checks, in order:
#   1. glibc / conda env / tool versions (cleanifier needs Python>=3.12 + numba on glibc 2.17)
#   2. tiny custom Cleanifier index (so the real 7 GB human index is not needed for the test)
#   3. fastp -> FIFO -> mbuffer -> FIFO -> cleanifier filter -> bgzip, exactly as the Nextflow task does it
#   4. shared-memory index: cleanifier load / filter --shared / cleanifier remove
#   5. Nextflow starts offline and the workflow parses (--help)
# Usage:  bash kbds_login_check.sh [conda_env_name_or_path]   (default: metacrc-qc)
set -uo pipefail

ENV_NAME=${1:-metacrc-qc}
# conda.sh: personal miniconda first, then the system one (override with CONDA_SH=...)
if [ -z "${CONDA_SH:-}" ]; then
    for c in "$HOME/miniconda3/etc/profile.d/conda.sh" "$HOME/miniforge3/etc/profile.d/conda.sh" /apps/application/miniconda3/etc/profile.d/conda.sh; do
        [ -f "$c" ] && CONDA_SH=$c && break
    done
fi
CONDA_SH=${CONDA_SH:-/apps/application/miniconda3/etc/profile.d/conda.sh}
METAFUN=${METAFUN:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}
WORK=${TMPDIR:-$HOME/tmp}/metacrc_qc_check.$$
mkdir -p "$WORK" && cd "$WORK"
export NXF_OFFLINE=true
pass=0; fail=0
ok()   { echo "  [PASS] $*"; pass=$((pass+1)); }
bad()  { echo "  [FAIL] $*"; fail=$((fail+1)); }
step() { echo; echo "== $*"; }

step "1. environment"
echo "  host: $(hostname)   glibc: $(ldd --version 2>&1 | head -n 1)"
if [ -f "$CONDA_SH" ]; then source "$CONDA_SH"; fi
if conda activate "$ENV_NAME" 2>/dev/null; then ok "conda activate $ENV_NAME"; else bad "conda activate $ENV_NAME (CONDA_SH=$CONDA_SH)"; fi
for t in python3 cleanifier fastp bgzip mbuffer nextflow java; do
    if command -v $t >/dev/null 2>&1; then echo "  $t: $(command -v $t)"; else bad "$t not on PATH"; fi
done
python3 --version 2>&1 | sed 's/^/  /'
if v=$(cleanifier --version 2>&1 | tail -n 1); then ok "cleanifier --version -> $v"; else bad "cleanifier --version failed: $v"; fi
fastp --version 2>&1 | sed 's/^/  /'
bgzip --version 2>&1 | head -n 1 | sed 's/^/  /'
java -version 2>&1 | head -n 1 | sed 's/^/  /'
cleanifier filter --help > cleanifier_filter_help.txt 2>&1 && ok "cleanifier filter --help saved to $WORK/cleanifier_filter_help.txt" || bad "cleanifier filter --help"
fastp --help > fastp_help.txt 2>&1
missing=""; for o in detect_adapter_for_pe trim_poly_g cut_front cut_tail cut_window_size cut_mean_quality length_required n_base_limit dont_eval_duplication report_title; do
    grep -q -- "--$o " fastp_help.txt || missing="$missing --$o"
done
[ -z "$missing" ] && ok "fastp accepts every option used in params.fastp_args" || bad "fastp lacks:$missing (adjust --fastp_args)"
grep -qi "dedup" fastp_help.txt && echo "  fastp dedup: $(grep -m1 -- '--dedup ' fastp_help.txt | sed 's/^ *//' | cut -c1-80)"
python3 -c "import numba, numpy; print('  numba', numba.__version__, 'numpy', numpy.__version__)" || bad "numba import"

step "2. tiny index + synthetic reads"
python3 - <<'EOF'
import random, gzip
random.seed(7)
rs = lambda n: ''.join(random.choice('ACGT') for _ in range(n))
host = rs(200000)
with open('host.fa', 'w') as f:
    f.write('>host\n' + '\n'.join(host[i:i+80] for i in range(0, len(host), 80)) + '\n')
rc = lambda s: s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
def pair(src):
    p = random.randint(0, len(src) - 300); f = src[p:p+300]
    return f[:150], rc(f)[:150]
with gzip.open('t_1.fastq.gz', 'wt') as o1, gzip.open('t_2.fastq.gz', 'wt') as o2:
    for i in range(2000):
        r1, r2 = pair(host) if random.random() < 0.3 else pair(rs(400))
        o1.write(f'@r{i} 1\n{r1}\n+\n{"I"*150}\n'); o2.write(f'@r{i} 2\n{r2}\n+\n{"I"*150}\n')
print('  synthetic reads: 2000 pairs, ~30% host')
EOF
if cleanifier index --index tinyhost --files host.fa -n 500000 -b 4 --mask '######_#######_###_#######_######' \
      --fill 0.8 --subtables 3 --threads-read 1 --threads-split 1 --maxwalk 2000 > index.log 2>&1; then
    ok "cleanifier index (tiny Cuckoo filter): $(ls tinyhost.*)"
else
    bad "cleanifier index failed (see $WORK/index.log)"; tail -n 5 index.log
fi

step "3. fastp -> FIFO -> mbuffer -> FIFO -> cleanifier filter -> bgzip (as in the Nextflow task)"
mkfifo t_1.raw t_2.raw t_fp_1.fastq t_fp_2.fastq
t0=$(date +%s)
mbuffer -q -f -m 64M -s 1M -i t_1.raw -o t_fp_1.fastq 2> mbuf1.log &
mbuffer -q -f -m 64M -s 1M -i t_2.raw -o t_fp_2.fastq 2> mbuf2.log &
fastp --in1 t_1.fastq.gz --in2 t_2.fastq.gz --out1 t_1.raw --out2 t_2.raw --thread 2 \
    --detect_adapter_for_pe --trim_poly_g --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
    --length_required 50 --n_base_limit 5 --dont_eval_duplication --json t.fastp.json --html t.fastp.html 2> t.fastp.log &
fp=$!
timeout 600 cleanifier filter --index tinyhost --fastq t_fp_1.fastq --pairs t_fp_2.fastq --out "$WORK/.t.part" \
    --threads 3 --compression gz --compression-threads 1 --compression-level 4 > t.cleanifier.log 2>&1
rc_cl=$?; wait $fp; rc_fp=$?; wait
echo "  cleanifier rc=$rc_cl fastp rc=$rc_fp wall=$(( $(date +%s) - t0 )) s (includes numba JIT compile)"
rm -f t_1.raw t_2.raw t_fp_1.fastq t_fp_2.fastq
grep -A1 -E "^prefix\s+keep\s+filter" t.cleanifier.log | sed 's/^/  /'
if [ $rc_cl -eq 0 ] && [ $rc_fp -eq 0 ] && [ -s .t.part_1_keep.fastq.gz ] && [ -s .t.part_2_keep.fastq.gz ]; then
    ok "pipe chain produced .t.part_{1,2}_keep.fastq.gz"
else
    bad "pipe chain (see $WORK/t.cleanifier.log, t.fastp.log, mbuf*.log)"
fi
if python3 "$METAFUN/scripts/rawreadqc/run_summary.py" --run t --sample s --cohort c --fastp-json t.fastp.json \
      --cleanifier-log t.cleanifier.log --fastq .t.part_1_keep.fastq.gz .t.part_2_keep.fastq.gz --expected-bases 600000 > t.qc.tsv; then
    ok "run_summary.py consistency checks (fastp pairs == keep+filter, BGZF EOF, base_count)"; cut -f1-12 t.qc.tsv | sed 's/^/  /'
else
    bad "run_summary.py"
fi

step "4. shared-memory index (cleanifier load / filter --shared / cleanifier remove)"
if cleanifier load --name "$WORK/tinyhost" > load.log 2>&1 && [ -e /dev/shm/tinyhost ]; then
    ok "cleanifier load -> /dev/shm/tinyhost ($(stat -c %s /dev/shm/tinyhost) bytes)"
else
    bad "cleanifier load"; cat load.log
fi
zcat .t.part_1_keep.fastq.gz | head -n 400 > s_1.fastq; zcat .t.part_2_keep.fastq.gz | head -n 400 > s_2.fastq
if cleanifier filter --index "$WORK/tinyhost" --shared --threads 3 --fastq s_1.fastq --pairs s_2.fastq --count > shared.log 2>&1; then
    ok "cleanifier filter --shared (count mode) works"
else
    bad "cleanifier filter --shared"; tail -n 5 shared.log
fi
cleanifier remove --name "$WORK/tinyhost" > /dev/null 2>&1
if [ ! -e /dev/shm/tinyhost ]; then ok "cleanifier remove deleted /dev/shm/tinyhost"; else bad "/dev/shm/tinyhost still present after cleanifier remove"; rm -f /dev/shm/tinyhost; fi

step "5. Nextflow offline start + workflow parse"
if nextflow -version > nf_version.txt 2>&1; then ok "nextflow -version: $(grep -m1 version nf_version.txt | xargs)"; else bad "nextflow -version"; tail -n 3 nf_version.txt; fi
if nextflow run "$METAFUN/nf_scripts/RAWREAD_QC_cleanifier.nf" -c "$METAFUN/config/kbds.config" -profile kbds --help > nf_help.txt 2>&1; then
    ok "workflow parses with -c config/kbds.config -profile kbds (--help)"
else
    bad "workflow --help failed (see $WORK/nf_help.txt)"; tail -n 8 nf_help.txt
fi

step "summary"
echo "  PASS=$pass FAIL=$fail   (files kept in $WORK)"
[ $fail -eq 0 ]
