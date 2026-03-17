###
# 01b_run_ngsadmix_LOCAL_byRF_arrayFiles_loopKsSeeds.slurm does this:
# Runs NGSadmix as a SLURM array over local per-RF beagle files listed one per line,
# then loops K=2..12 and multiple seeds for each file, grouping outputs by LG and thinning.
###
#!/usr/bin/env bash
#SBATCH -J ngsadmix_local_byRF
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH -t 2-00:00:00
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail

# ---- ENV ----
source ~/.bashrc || true

# If mamba is not fully initialized in batch, uncomment one of these if needed:
# source ~/mambaforge/etc/profile.d/conda.sh 2>/dev/null || true
# source ~/miniforge3/etc/profile.d/conda.sh 2>/dev/null || true
# source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null || true

mamba activate assembly  # or: conda activate assembly

# ---- INPUTS ----
LIST="${1:?ERROR: provide beagle_local_byRF.list as arg1}"
WORKDIR="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/popstruct_thin/05_ngsadmix_local_byRF"

# Base directory used only when paths in LIST are relative
BASE_FOR_BEAGLES="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/popstruct_thin"

MINMAF="0.05"
K_MIN=2
K_MAX=12
SEEDS=(1 2 3)

THREADS=${SLURM_CPUS_PER_TASK:-32}

mkdir -p "${WORKDIR}/logs"
cd "${WORKDIR}"

echo "[INFO] PWD=$(pwd)"
echo "[INFO] LIST=$LIST"
echo "[INFO] BASE_FOR_BEAGLES=$BASE_FOR_BEAGLES"

# Make sure NGSadmix exists in batch context
if ! command -v NGSadmix >/dev/null 2>&1; then
  echo "[ERROR] NGSadmix not found in PATH inside batch job." >&2
  echo "[HINT] Put it in the conda env, or module load it, or export PATH in ~/.bashrc." >&2
  which NGSadmix || true
  exit 1
fi

[[ -f "$LIST" ]] || { echo "[ERROR] List file not found: $LIST" >&2; exit 1; }

NFILES=$(wc -l < "$LIST")

task=${SLURM_ARRAY_TASK_ID:?need array id}
if (( task < 1 || task > NFILES )); then
  echo "[ERROR] SLURM_ARRAY_TASK_ID=$task out of range 1..$NFILES" >&2
  exit 1
fi

IN=$(sed -n "${task}p" "$LIST")
[[ -n "$IN" ]] || { echo "[ERROR] Empty IN at line $task" >&2; exit 1; }

# If the list contains relative paths, convert them to absolute paths
if [[ "$IN" != /* ]]; then
  IN="${BASE_FOR_BEAGLES%/}/${IN}"
fi

# Clean accidental duplicate slashes
IN=$(echo "$IN" | sed 's#//*#/#g')

[[ -s "$IN" ]] || { echo "[ERROR] Missing input: $IN" >&2; exit 1; }

bn=$(basename "$IN")              # e.g. catfish.C_gar_LG27.rf.txt.thin_500.beagle.gz
tag="${bn%.beagle.gz}"            # e.g. catfish.C_gar_LG27.rf.txt.thin_500
LG=$(echo "$tag" | sed -n 's/.*\(LG[0-9]\+\).*/\1/p')
thin=$(echo "$tag" | sed -n 's/.*thin_\([0-9]\+\).*/\1/p')

[[ -n "$LG" ]] || LG="LGNA"
[[ -n "$thin" ]] || thin="thinNA"

OUTROOT="runs_${LG}_thin${thin}"
mkdir -p "$OUTROOT"

echo "[INFO] task=$task/$NFILES"
echo "[INFO] IN=$IN"
echo "[INFO] LG=$LG thin=$thin"
echo "[INFO] threads=$THREADS"
echo "[INFO] OUTROOT=$OUTROOT"

for SEED in "${SEEDS[@]}"; do
  for K in $(seq "$K_MIN" "$K_MAX"); do
    OUTPREFIX="${OUTROOT}/${tag}_K$(printf "%02d" "$K")_seed${SEED}"
    echo "[INFO] Running K=$K seed=$SEED -> $OUTPREFIX"

    NGSadmix \
      -likes "$IN" \
      -K "$K" \
      -minMaf "$MINMAF" \
      -P "$THREADS" \
      -seed "$SEED" \
      -o "$OUTPREFIX"
  done
done
