###
# 06_run_pcangsd_byLG_thin200_500_1000_K2_12.slurm does this:
# Runs PCAngsd as a SLURM array across per-LG beagle files from the
# 200 bp, 500 bp, and 1000 bp thinned datasets, testing admixture K=2..12
# and writing outputs grouped by thinning, linkage group, and K.
###
#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH -J pcangsd_byLG
#SBATCH -o 05_pcangsd_byLG/logs/%x_%A_%a.out
#SBATCH -e 05_pcangsd_byLG/logs/%x_%A_%a.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/popstruct_thin"
OUTBASE="${BASE}/05_pcangsd_byLG"
mkdir -p "${OUTBASE}/logs"

THREADS=${SLURM_CPUS_PER_TASK:-32}

# sample names in SAME ORDER as bamlist used by ANGSD for this beagle
SAMPLES="${BASE}/list_of_samples_one_per_line_same_bamfile_list.tsv"
[[ -s "$SAMPLES" ]] || { echo "[ERROR] missing samples file: $SAMPLES" >&2; exit 1; }

# precomputed lists
L200="${OUTBASE}/beagle_LG_thin_200.list"
L500="${OUTBASE}/beagle_LG_thin_500.list"
L1000="${OUTBASE}/beagle_LG_thin_1000.list"
[[ -s "$L200" && -s "$L500" && -s "$L1000" ]] || {
  echo "[ERROR] missing list files: $L200 / $L500 / $L1000" >&2
  exit 1
}

THINS=(200 500 1000)
KS=($(seq 2 12))
NK=${#KS[@]}

N200=$(wc -l < "$L200")
N500=$(wc -l < "$L500")
N1000=$(wc -l < "$L1000")

TOTAL=$(((N200 + N500 + N1000) * NK))

task=${SLURM_ARRAY_TASK_ID:?need array id}
if (( task < 1 || task > TOTAL )); then
  echo "[ERROR] SLURM_ARRAY_TASK_ID=$task out of range 1..$TOTAL" >&2
  exit 1
fi

idx=$((task-1))
file_global_i=$(( idx / NK ))
k_i=$(( idx % NK ))
K=${KS[$k_i]}

if (( file_global_i < N200 )); then
  thin=200
  file_i=$file_global_i
  LIST="$L200"
elif (( file_global_i < N200 + N500 )); then
  thin=500
  file_i=$(( file_global_i - N200 ))
  LIST="$L500"
else
  thin=1000
  file_i=$(( file_global_i - N200 - N500 ))
  LIST="$L1000"
fi

IN=$(sed -n "$((file_i+1))p" "$LIST")
[[ -n "$IN" && -s "$IN" ]] || { echo "[ERROR] Bad IN for file_i=$file_i LIST=$LIST" >&2; exit 1; }

bn=$(basename "$IN")
tag="${bn%.beagle.gz}"

# extract LG name for output grouping
LG=$(echo "$bn" | grep -oE 'LG[0-9]+' | head -n1)
[[ -n "${LG:-}" ]] || LG="LGNA"

OUTDIR="${OUTBASE}/thin_${thin}/${LG}/K${K}"
mkdir -p "$OUTDIR"
OUT="${OUTDIR}/${tag}.pcangsd"

echo "[INFO] task=$task / $TOTAL"
echo "[INFO] thin=$thin  LG=$LG  K=$K"
echo "[INFO] IN=$IN"
echo "[INFO] OUT=$OUT"
echo "[INFO] threads=$THREADS"

pcangsd \
  --beagle "$IN" \
  --eig 6 \
  --threads "$THREADS" \
  --maf 0.05 \
  --admix \
  --admix-K "$K" \
  --admix-seed 1 \
  --tree \
  --tree-samples "$SAMPLES" \
  --out "$OUT"
