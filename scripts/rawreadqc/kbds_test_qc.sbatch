#!/bin/bash
#SBATCH -J qc_PRJNA531273
#SBATCH -p cpu32-only
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -t 12:00:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#
# First-test wrapper for RAWREAD_QC_cleanifier on K-BDS (one node, one Nextflow run).
# This is only for the PRJNA531273 timing test; the production wrapper is produced by the
# HeadQuarters runbook. Submit from the run directory:
#   cd ~/metaCRC/runs/qc/PRJNA531273 && sbatch ~/metaFun/scripts/rawreadqc/kbds_test_qc.sbatch
# Override anything with environment variables, e.g.  COHORT=PRJNA531273 CPUS=4 sbatch ...
set -euo pipefail

COHORT=${COHORT:-PRJNA531273}
METAFUN=${METAFUN:-$HOME/metaFun}                 # git checkout (branch kbds-rawreadqc)
QC_ROOT=${QC_ROOT:-$HOME/metaCRC/qc}
RUN_DIR=${RUN_DIR:-$HOME/metaCRC/runs/qc/$COHORT}
WORK_DIR=${WORK_DIR:-$HOME/metaCRC/work/qc/$COHORT}
INDEX=${INDEX:-$(ls $HOME/metaCRC/db/cleanifier/*.info 2>/dev/null | head -n 1 | sed 's/\.info$//')}
CONDA_SH=${CONDA_SH:-/apps/application/miniconda3/etc/profile.d/conda.sh}
CONDA_ENV=${CONDA_ENV:-metacrc-qc}
CPUS=${CPUS:-4}

source "$CONDA_SH"
conda activate "$CONDA_ENV"
export NXF_OFFLINE=true                 # compute nodes may have no outbound network
export NXF_ANSI_LOG=false
export NXF_OPTS='-Xms512m -Xmx2g'
export TMPDIR=${TMPDIR:-/tmp}

echo "node=$(hostname) nproc=$(nproc) mem=$(free -g | awk '/Mem/{print $2}')G shm=$(df -h /dev/shm | awk 'NR==2{print $2}') date=$(date)"
echo "index=$INDEX  run_dir=$RUN_DIR  work_dir=$WORK_DIR  qc_root=$QC_ROOT"
command -v nextflow fastp cleanifier bgzip
echo "--- facts still unknown before this test (for the runbook):"
echo "outbound: $(timeout 5 curl -sS -o /dev/null -w '%{http_code}' https://zenodo.org 2>&1 || echo none)"
echo "apptainer: $(command -v apptainer || echo none)   sacct: $(command -v sacct || echo none)"
echo "account file: $(ls /scratch/account/ 2>/dev/null | head -n 3 | tr '\n' ' ')"

mkdir -p "$RUN_DIR"
cd "$RUN_DIR"
[ -s samplesheet.csv ] || { echo "samplesheet.csv missing in $RUN_DIR"; exit 1; }

nextflow run "$METAFUN/nf_scripts/RAWREAD_QC_cleanifier.nf" -c "$METAFUN/config/nextflow.config" -profile kbds \
    --samplesheet "$RUN_DIR/samplesheet.csv" \
    --outdir "$QC_ROOT" \
    --run_dir "$RUN_DIR" --work_dir "$WORK_DIR" \
    --cleanifier_index "$INDEX" \
    --cpus "$CPUS" \
    -resume

echo "--- process names in trace.txt:"; cut -f4 "$RUN_DIR/trace.txt" | sort -u
echo "--- done: $(ls "$QC_ROOT/$COHORT"/*.done 2>/dev/null | wc -l) markers"
python3 "$METAFUN/scripts/rawreadqc/trace_report.py" --trace "$RUN_DIR/trace.txt" \
    --qc-summary "$QC_ROOT/$COHORT/$COHORT.qc_summary.tsv" --node-rate 0.5 --total-runs 8420 || true
