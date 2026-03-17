###
# run_angsd_saf_chunks_v3.slurm does this:
# Runs ANGSD SAF generation and realSFS bootstrap estimation for one chromosome/region
# chunk selected from a chunk list via SLURM array index.
###
#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=80
#SBATCH --mem=64GB
#SBATCH -t 1-00:00:00
#SBATCH -J angsd_saf_chunks
#SBATCH -o angsd_saf_chunks.%A_%a.out
#SBATCH -e angsd_saf_chunks.%A_%a.err
set -euo pipefail

source ~/.bashrc
mamba activate assembly

# ---- args ----
CHUNKLIST="${1:?ERROR: provide chunk_rf.list as argv1}"

# array index 0-based like your pattern
IDX="${SLURM_ARRAY_TASK_ID}"
LINE=$((IDX+1))

RF_FILE=$(sed -n "${LINE}p" "$CHUNKLIST")
if [[ -z "${RF_FILE}" || ! -s "${RF_FILE}" ]]; then
  echo "[ERROR] Could not read RF_FILE from ${CHUNKLIST} at line ${LINE}" >&2
  exit 1
fi

# derive a tag for naming (basename without extension)
TAG=$(basename "$RF_FILE")
TAG=${TAG%.rf.txt}

BASE="/scratch/lt200308-agbsci/Quentin_project"
REF="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
BAMLIST="${BASE}/bamlist.pp.samechr.tlenP99.filtered.txt"
OUTDIR="${BASE}/popstruct_global"
P="${SLURM_CPUS_PER_TASK}"

# callable mask in ANGSD native format (1-based) AND indexed
CALLABLE_SITES="${BASE}/fClaHyb_Gar_LG.mask_regions.normalACGT.renamed.1-based.angsd"

MINQ=25
MINMAPQ=25
MININD=200
MINDEPTHIND=3
MAXDEPTHIND=57

# sanity for sites index
if [[ ! -s "${CALLABLE_SITES}.idx" || ! -s "${CALLABLE_SITES}.bin" ]]; then
  echo "[ERROR] Missing sites index: expected ${CALLABLE_SITES}.idx and ${CALLABLE_SITES}.bin" >&2
  exit 1
fi

# per-chunk folder
CHUNKDIR="${OUTDIR}/chunks/${TAG}"
mkdir -p "${CHUNKDIR}/01_saf" "${CHUNKDIR}/02_sfs" "${CHUNKDIR}/logs"

SAF_PREFIX="${CHUNKDIR}/01_saf/catfish.${TAG}"
SFS_OUT="${CHUNKDIR}/02_sfs/catfish.${TAG}.sfs"

echo "[INFO] IDX=${IDX} line=${LINE}"
echo "[INFO] RF_FILE=${RF_FILE}"
echo "[INFO] TAG=${TAG}"
echo "[INFO] SAF_PREFIX=${SAF_PREFIX}"
echo "[INFO] SFS_OUT=${SFS_OUT}"

# 1) SAF for this chunk
angsd -b "${BAMLIST}" -ref "${REF}" \
  -GL 1 \
  -doSaf 5 -doMajorMinor 1 -doMaf 1 \
  -doCounts 1 -doDepth 1 \
  -minQ "${MINQ}" -minMapQ "${MINMAPQ}" -baq 1 -C 50 \
  -setMinDepthInd "${MINDEPTHIND}" -setMaxDepthInd "${MAXDEPTHIND}" \
  -minInd "${MININD}" \
  -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 \
  -rf "${RF_FILE}" \
  -sites "${CALLABLE_SITES}" \
  -P "${P}" \
  -out "${SAF_PREFIX}" \
  2> "${CHUNKDIR}/logs/01_angsd_saf.${TAG}.log"

# 2) realSFS for this chunk
realSFS "${SAF_PREFIX}.saf.idx" -bootstrap 100 -P "${P}" \
  > "${SFS_OUT}" \
  2> "${CHUNKDIR}/logs/02_realSFS.${TAG}.log"

echo "[DONE] ${TAG}"


# You can also merge saf files. Then you can run angsd on the separate chromosomes and merge afterwars
# SEE THIS https://github.com/ANGSD/angsd/issues/66
# SEE THIS https://www.popgen.dk/angsd/index.php/RealSFS#Merge_SAF_files

