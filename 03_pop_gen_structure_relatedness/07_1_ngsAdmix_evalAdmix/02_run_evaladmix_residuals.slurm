###
# 02_run_evaladmix_residuals.slurm does this:
# Runs evalAdmix as a SLURM array over completed NGSadmix runs,
# using the selected whole-genome beagle file and writing one evalAdmix
# output prefix per qopt run.
###
#!/usr/bin/env bash
#SBATCH -J evalAdmix
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -t 1-00:00:00
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail
source ~/.bashrc

# ---- EDIT THESE ----
WORKDIR="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/popstruct_thin/05_ngsadmix_global"
BEAGLE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/popstruct_thin/04_beagle_wholegenome_majmin/catfish.wg.byRF.thin_500.beagle.gz"
# -------------------

cd "${WORKDIR}"
THREADS=${SLURM_CPUS_PER_TASK:-16}

tag="thin500"
RUNROOT="runs_${tag}"
OUTROOT="eval_${tag}"
mkdir -p "${OUTROOT}"

# Build a list of all NGSadmix outputs (OUTPREFIX)
# NGSadmix produces e.g. runs_thin500/thin500_K02_seed1.qopt
mapfile -t QOPTS < <(ls -1 ${RUNROOT}/${tag}_K*_seed*.qopt 2>/dev/null | sort)
if (( ${#QOPTS[@]} == 0 )); then
  echo "[ERROR] No .qopt found under ${RUNROOT}. Did NGSadmix finish?" >&2
  exit 1
fi

task=${SLURM_ARRAY_TASK_ID:?need array id}
if (( task < 1 || task > ${#QOPTS[@]} )); then
  echo "[ERROR] task=$task out of range 1..${#QOPTS[@]}" >&2
  exit 1
fi

QOPT=${QOPTS[$((task-1))]}
base=$(basename "$QOPT" .qopt)
OUTP="${OUTROOT}/${base}"

echo "[INFO] QOPT=$QOPT"
echo "[INFO] OUTP=$OUTP"
echo "[INFO] threads=$THREADS"

# evalAdmix typical usage:
# evalAdmix -beagle <beagle.gz> -fname <qopt> -o <outprefix> -P <threads>
evalAdmix \
  -beagle "$BEAGLE" \
  -fname "$QOPT" \
  -o "$OUTP" \
  -P "$THREADS"
