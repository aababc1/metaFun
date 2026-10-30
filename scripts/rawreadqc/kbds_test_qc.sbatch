#!/bin/bash
#SBATCH -J qc_PRJNA531273
#SBATCH -p cpu32-only
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -t 12:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err
#
# First-test wrapper for RAWREAD_QC_cleanifier on K-BDS (one node, one Nextflow run).
# This is only for the PRJNA531273 timing test; the production wrapper is produced by the
# HeadQuarters runbook. Submit from the run directory:
#   cd ~/metaCRC/runs/qc/PRJNA531273 && sbatch ~/metaFun/scripts/rawreadqc/kbds_test_qc.sbatch
# Override anything with environment variables, e.g.  COHORT=PRJNA531273 CPUS=4 sbatch ...
# Extra workflow arguments: NF_EXTRA="--bases_tolerance 0.10" (sbatch exports the environment to the job)
set -euo pipefail

COHORT=${COHORT:-PRJNA531273}
METAFUN=${METAFUN:-$HOME/metaFun}                 # git checkout (branch kbds-rawreadqc)
QC_ROOT=${QC_ROOT:-$HOME/metaCRC/qc}
RUN_DIR=${RUN_DIR:-$HOME/metaCRC/runs/qc/$COHORT}
WORK_DIR=${WORK_DIR:-$HOME/metaCRC/work/qc/$COHORT}
INDEX=${INDEX:-$(find $HOME/metaCRC/db/cleanifier -name '*.info' 2>/dev/null | head -n 1 | sed 's/\.info$//')}
CONDA_ENV=${CONDA_ENV:-metacrc-qc}
if [ -z "${CONDA_SH:-}" ]; then
    for c in "$HOME/miniconda3/etc/profile.d/conda.sh" "$HOME/miniforge3/etc/profile.d/conda.sh" /apps/application/miniconda3/etc/profile.d/conda.sh; do
        [ -f "$c" ] && CONDA_SH=$c && break
    done
fi
CPUS=${CPUS:-4}
NF_EXTRA=${NF_EXTRA:-}                   # extra `nextflow run` arguments, e.g. NF_EXTRA="--bases_tolerance 0.10"

source "$CONDA_SH"
conda activate "$CONDA_ENV"
export NXF_OFFLINE=true                 # compute nodes may have no outbound network
export NXF_ANSI_LOG=false
export NXF_OPTS='-Xms512m -Xmx2g'
export TMPDIR=${TMPDIR:-/tmp}

echo "node=$(hostname) nproc=$(nproc) mem=$(free -g | awk '/Mem/{print $2}')G shm=$(df -h /dev/shm | awk 'NR==2{print $2}') date=$(date)"
echo "index=$INDEX  run_dir=$RUN_DIR  work_dir=$WORK_DIR  qc_root=$QC_ROOT"
command -v nextflow fastp cleanifier bgzip mbuffer
[ -n "$INDEX" ] && [ -f "$INDEX.info" ] || { echo "Cleanifier index not found (INDEX=$INDEX)"; exit 1; }
echo "--- facts still unknown before this test (for the runbook):"
echo "outbound: $(timeout 5 curl -sS -o /dev/null -w '%{http_code}' https://zenodo.org 2>&1 || echo none)"
echo "apptainer: $(command -v apptainer || echo none)   sacct: $(command -v sacct || echo none)"
echo "account file: $(ls /scratch/account/ 2>/dev/null | head -n 3 | tr '\n' ' ')"

mkdir -p "$RUN_DIR"
cd "$RUN_DIR"
[ -s samplesheet.csv ] || { echo "samplesheet.csv missing in $RUN_DIR"; exit 1; }

nextflow run "$METAFUN/nf_scripts/RAWREAD_QC_cleanifier.nf" -c "$METAFUN/config/kbds.config" -profile kbds \
    --samplesheet "$RUN_DIR/samplesheet.csv" \
    --outdir "$QC_ROOT" \
    --run_dir "$RUN_DIR" --work_dir "$WORK_DIR" \
    --cleanifier_index "$INDEX" \
    --cpus "$CPUS" $NF_EXTRA \
    -resume && nf_rc=0 || nf_rc=$?

# the workflow removes the shared-memory index itself; this is the safety net for an aborted run
cleanifier remove --name "$INDEX" 2>/dev/null || true
[ -e "/dev/shm/$(basename "$INDEX")" ] && rm -f "/dev/shm/$(basename "$INDEX")" || true

echo "--- nextflow exit: $nf_rc"
echo "--- process names in trace.txt:"; cut -f4 "$RUN_DIR/trace.txt" | sort -u
echo "--- done: $(ls "$QC_ROOT/$COHORT"/*.done 2>/dev/null | wc -l) markers"
python3 "$METAFUN/scripts/rawreadqc/trace_report.py" --trace "$RUN_DIR/trace.txt" \
    --qc-summary "$QC_ROOT/$COHORT/$COHORT.qc_summary.tsv" --node-rate 0.5 --total-runs 8420 || true
