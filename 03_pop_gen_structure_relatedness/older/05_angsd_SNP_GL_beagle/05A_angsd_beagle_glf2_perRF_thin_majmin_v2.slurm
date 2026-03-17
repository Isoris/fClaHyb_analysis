###
# 05A_angsd_beagle_glf2_perRF_thin_majmin_v2.slurm does this:
# Runs ANGSD per RF chunk to generate beagle genotype likelihood files
# using thinned site lists that preserve major/minor alleles, writing
# outputs into a thinning-specific per-RF beagle directory.
###
#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=80
#SBATCH --ntasks=1
#SBATCH --mem=237GB
#SBATCH -t 2-00:00:00
#SBATCH -J angsd_beagle_rf
#SBATCH -o popstruct_thin/logs/%x.%A_%a.out
#SBATCH -e popstruct_thin/logs/%x.%A_%a.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

###############################################################################
# STEP 5A) BEAGLE GLs per RF chunk using thinned sites WITH major/minor columns
#
# Usage:
#   RFLIST=chunk_rf.list
#   N=$(wc -l < "$RFLIST")
#   sbatch --array=0-$((N-1))%8 05A_angsd_beagle_glf2_perRF_thin_majmin.slurm "$RFLIST" 200
#   sbatch --array=0-$((N-1))%8 05A_angsd_beagle_glf2_perRF_thin_majmin.slurm "$RFLIST" 500
#   sbatch --array=0-$((N-1))%8 05A_angsd_beagle_glf2_perRF_thin_majmin.slurm "$RFLIST" 1000
#
# Args:
#   $1 = chunk_rf.list  (one rf file per line)
#   $2 = thinning bp    (e.g. 200 / 500 / 1000 / 5000 / 10000 / 25000)
###############################################################################

RFLIST="${1:?ERROR: provide chunk_rf.list as arg1}"
W="${2:?ERROR: provide thinning bp as arg2 (e.g. 200/500/1000)}"

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
REF="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
BAMLIST="${BASE}/bamlist.pp.samechr.tlenP99.filtered.txt"

OUTDIR="${BASE}/popstruct_thin"
SITES_FINAL="${OUTDIR}/03_sites/sites.thin_${W}.majmin.tsv"

# Put beagles under a thin-specific folder so runs never overwrite each other
BEAGLEDIR="${OUTDIR}/04_beagle_byRF_majmin/thin_${W}"
LOGDIR="${OUTDIR}/logs"
THREADS="${SLURM_CPUS_PER_TASK:-80}"

mkdir -p "${BEAGLEDIR}" "${LOGDIR}"

# Filters (same as your current)
MINQ=25
MINMAPQ=25
MINDEPTHIND=3
MAXDEPTHIND=57
MININD=200

[[ -s "${SITES_FINAL}" ]] || { echo "[ERROR] Missing sites file: ${SITES_FINAL}" >&2; exit 1; }

# Optional: ensure index exists (recommended)
if [[ ! -s "${SITES_FINAL}.idx" || ! -s "${SITES_FINAL}.bin" ]]; then
  echo "[WARN] sites index missing for ${SITES_FINAL} (.idx/.bin). Indexing now..."
  angsd sites index "${SITES_FINAL}"
fi

RF_FILE="$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$RFLIST")"
[[ -n "${RF_FILE}" ]] || { echo "[ERROR] Empty RF_FILE for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
[[ -s "${RF_FILE}" ]] || { echo "[ERROR] RF file not found: ${RF_FILE}" >&2; exit 1; }

TAG="$(basename "$RF_FILE")"
TAG="${TAG//[^A-Za-z0-9_.-]/_}"

PREFIX="${BEAGLEDIR}/catfish.${TAG}.thin_${W}"
LOG="${LOGDIR}/05A_angsd_beagle.${TAG}.thin_${W}.log"

echo "[INFO] task=${SLURM_ARRAY_TASK_ID}"
echo "[INFO] RF_FILE=${RF_FILE}"
echo "[INFO] SITES_FINAL=${SITES_FINAL}"
echo "[INFO] OUT=${PREFIX}.beagle.gz"
echo "[INFO] THREADS=${THREADS}"

set -x
angsd -b "${BAMLIST}" -ref "${REF}" \
  -GL 1 \
  -doMajorMinor 3 \
  -doGlf 2 \
  -doCounts 1 -doDepth 1 \
  -minQ "${MINQ}" -minMapQ "${MINMAPQ}" -baq 1 -C 50 \
  -setMinDepthInd "${MINDEPTHIND}" -setMaxDepthInd "${MAXDEPTHIND}" \
  -minInd "${MININD}" \
  -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 \
  -rf "${RF_FILE}" \
  -sites "${SITES_FINAL}" \
  -P "${THREADS}" \
  -out "${PREFIX}" \
  > "${LOG}" 2>&1
set +x

echo "[DONE] ${PREFIX}.beagle.gz"
