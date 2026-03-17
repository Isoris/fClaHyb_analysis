###
# 01b_run_ngsadmix_LOCAL_byRF_list_Ksweep_seeds.slurm does this:
# Runs NGSadmix as a SLURM array across a list of local per-RF beagle files,
# sweeping K=2..12 and multiple seeds, and writes outputs grouped by input tag.
###
#!/usr/bin/env bash
#SBATCH -J ngsadmix_LOCAL
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH -t 2-00:00:00
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

WORKDIR="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/popstruct_thin/05_ngsadmix_local_byRF"
LIST="${1:?ERROR: provide beagle_local_byRF.list as arg1}"

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

NFILES=$(wc -l < "$LIST")
[[ "$NFILES" -gt 0 ]] || { echo "[ERROR] Empty list: $LIST" >&2; exit 1; }

NK=$((K_MAX - K_MIN + 1))
NS=${#SEEDS[@]}
TOTAL=$((NFILES * NK * NS))

task=${SLURM_ARRAY_TASK_ID:?need array id}
if (( task < 1 || task > TOTAL )); then
  echo "[ERROR] SLURM_ARRAY_TASK_ID=$task out of range 1..$TOTAL" >&2
  exit 1
fi
idx=$((task-1))

file_i=$(( idx / (NK*NS) ))
rem=$(( idx % (NK*NS) ))
seed_i=$(( rem / NK ))
k_i=$(( rem % NK ))

SEED=${SEEDS[$seed_i]}
K=$((K_MIN + k_i))

BEAGLE=$(sed -n "$((file_i+1))p" "$LIST")
[[ -n "$BEAGLE" ]] || { echo "[ERROR] Empty BEAGLE for file_i=$file_i" >&2; exit 1; }
[[ -s "$BEAGLE" ]] || { echo "[ERROR] Missing input: $BEAGLE" >&2; exit 1; }

bn=$(basename "$BEAGLE")
tag="${bn%.beagle.gz}"   # keeps LG + thin in the tag automatically

OUTROOT="runs_${tag}"
mkdir -p "${OUTROOT}"

OUTPREFIX="${OUTROOT}/${tag}_K$(printf "%02d" "$K")_seed${SEED}"

echo "[INFO] BEAGLE=$BEAGLE"
echo "[INFO] K=$K  SEED=$SEED  THREADS=$THREADS"
echo "[INFO] OUTPREFIX=$OUTPREFIX"

"${NGSADMIX_BIN}" \
  -likes "$BEAGLE" \
  -K "$K" \
  -minMaf "$MINMAF" \
  -P "$THREADS" \
  -seed "$SEED" \
  -o "$OUTPREFIX"
