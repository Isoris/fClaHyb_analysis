###
# 01_run_ngsadmix_Ksweep_seeds.slurm does this:
# Runs global NGSadmix as a SLURM array on one selected whole-genome beagle file,
# sweeping K=2..12 across multiple seeds and writing standard NGSadmix outputs.
###
#!/usr/bin/env bash
#SBATCH -J ngsadmix_Ksweep
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH -t 2-00:00:00
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail
source ~/.bashrc
NGSADMIX_BIN="$(which NGSadmix)"


# ---- EDIT THESE ----
WORKDIR="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/popstruct_thin/05_ngsadmix_global"
BEAGLE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/popstruct_thin/04_beagle_wholegenome_majmin/catfish.wg.byRF.thin_500.beagle.gz"
MINMAF="0.05"
K_MIN=2
K_MAX=12
SEEDS=(1 2 3)
THREADS=${SLURM_CPUS_PER_TASK:-32}
# -------------------

mkdir -p "${WORKDIR}/logs"
cd "${WORKDIR}"

# total tasks = (#K) * (#seeds)
NK=$((K_MAX - K_MIN + 1))
NS=${#SEEDS[@]}
TOTAL=$((NK * NS))

task=${SLURM_ARRAY_TASK_ID:?need array id}
if (( task < 1 || task > TOTAL )); then
  echo "[ERROR] SLURM_ARRAY_TASK_ID=$task out of range 1..$TOTAL" >&2
  exit 1
fi
idx=$((task-1))

seed_i=$(( idx / NK ))
k_i=$(( idx % NK ))

SEED=${SEEDS[$seed_i]}
K=$((K_MIN + k_i))

tag="thin500"   # purely for filenames (edit if you run other thins)
OUTROOT="runs_${tag}"
mkdir -p "${OUTROOT}"

OUTPREFIX="${OUTROOT}/${tag}_K$(printf "%02d" "$K")_seed${SEED}"

echo "[INFO] BEAGLE=$BEAGLE"
echo "[INFO] K=$K  SEED=$SEED  THREADS=$THREADS"
echo "[INFO] OUTPREFIX=$OUTPREFIX"

# NGSadmix outputs: .qopt (Q), .fopt, .log etc.
NGSadmix \
  -likes "$BEAGLE" \
  -K "$K" \
  -minMaf "$MINMAF" \
  -P "$THREADS" \
  -seed "$SEED" \
  -o "$OUTPREFIX"
