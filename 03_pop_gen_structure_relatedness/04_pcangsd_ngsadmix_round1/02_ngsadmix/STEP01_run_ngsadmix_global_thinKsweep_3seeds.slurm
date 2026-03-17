###
# 01_run_ngsadmix_GLOBAL_thinKsweep_seeds.slurm does this:
# Runs global NGSadmix as a SLURM array across multiple thinned whole-genome
# beagle files, sweeping K=2..12 and multiple seeds, and writes outputs
# grouped by thinning level.
###
#!/usr/bin/env bash
#SBATCH -J ngsadmix_GLOBAL
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH -t 2-00:00:00
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail
source ~/.bashrc

# >>> IMPORTANT: ensure NGSadmix is available inside the job <<<
# If you installed it in your conda env, activate it here:
mamba activate assembly

WORKDIR="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/popstruct_thin/05_ngsadmix_global"
BEAGLEDIR="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/popstruct_thin/04_beagle_byRF_majmin"

THINS=(200 500 1000)

MINMAF="0.05"
K_MIN=2
K_MAX=12
SEEDS=(1 2 3)
THREADS=${SLURM_CPUS_PER_TASK:-32}

mkdir -p "${WORKDIR}/logs"
cd "${WORKDIR}"

NGSADMIX_BIN="$(command -v NGSadmix || true)"
if [[ -z "${NGSADMIX_BIN}" ]]; then
  echo "[ERROR] NGSadmix not found in PATH inside SLURM job." >&2
  echo "[DEBUG] PATH=$PATH" >&2
  exit 1
fi
echo "[INFO] NGSADMIX_BIN=${NGSADMIX_BIN}"

# total tasks = (#thins) * (#K) * (#seeds)
NT=${#THINS[@]}
NK=$((K_MAX - K_MIN + 1))
NS=${#SEEDS[@]}
TOTAL=$((NT * NK * NS))

task=${SLURM_ARRAY_TASK_ID:?need array id}
if (( task < 1 || task > TOTAL )); then
  echo "[ERROR] SLURM_ARRAY_TASK_ID=$task out of range 1..$TOTAL" >&2
  exit 1
fi
idx=$((task-1))

thin_i=$(( idx / (NK*NS) ))
rem=$(( idx % (NK*NS) ))
seed_i=$(( rem / NK ))
k_i=$(( rem % NK ))

THIN=${THINS[$thin_i]}
SEED=${SEEDS[$seed_i]}
K=$((K_MIN + k_i))

BEAGLE="${BEAGLEDIR}/catfish.wholegenome.byRF.thin_${THIN}.beagle.gz"
[[ -s "$BEAGLE" ]] || { echo "[ERROR] Missing BEAGLE: $BEAGLE" >&2; exit 1; }

tag="thin${THIN}"
OUTROOT="runs_${tag}"
mkdir -p "${OUTROOT}"

OUTPREFIX="${OUTROOT}/${tag}_K$(printf "%02d" "$K")_seed${SEED}"

echo "[INFO] BEAGLE=$BEAGLE"
echo "[INFO] THIN=$THIN  K=$K  SEED=$SEED  THREADS=$THREADS"
echo "[INFO] OUTPREFIX=$OUTPREFIX"

"${NGSADMIX_BIN}" \
  -likes "$BEAGLE" \
  -K "$K" \
  -minMaf "$MINMAF" \
  -P "$THREADS" \
  -seed "$SEED" \
  -o "$OUTPREFIX"
