#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=80
#SBATCH --mem=64GB
#SBATCH -t 1-00:00:00
#SBATCH -J inv_saf_chunks
#SBATCH -o inv_saf_chunks.%A_%a.out
#SBATCH -e inv_saf_chunks.%A_%a.err
set -euo pipefail

source ~/.bashrc
mamba activate assembly

CHUNKLIST="${1:?ERROR: provide chunk_rf.list as argv1}"

IDX="${SLURM_ARRAY_TASK_ID}"
LINE=$((IDX+1))

RF_FILE=$(sed -n "${LINE}p" "$CHUNKLIST")
if [[ -z "${RF_FILE}" || ! -s "${RF_FILE}" ]]; then
  echo "[ERROR] Could not read RF_FILE from ${CHUNKLIST} at line ${LINE}" >&2
  exit 1
fi

TAG=$(basename "$RF_FILE")
TAG=${TAG%.rf.txt}

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
REF="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
BAMLIST="${BASE}/bamlist.pp.samechr.tlenP99.filtered.txt"
OUTDIR="${BASE}/inversion_localpca"
P="${SLURM_CPUS_PER_TASK}"

CALLABLE_SITES="${BASE}/fClaHyb_Gar_LG.mask_regions.inversion_acgt_allcase.1-based.angsd"

MINQ=20
MINMAPQ=60
MININD=203
SETMAXDEPTH=2850

if [[ ! -s "${CALLABLE_SITES}.idx" || ! -s "${CALLABLE_SITES}.bin" ]]; then
  echo "[ERROR] Missing sites index: expected ${CALLABLE_SITES}.idx and ${CALLABLE_SITES}.bin" >&2
  exit 1
fi

CHUNKDIR="${OUTDIR}/chunks/${TAG}"
mkdir -p "${CHUNKDIR}/01_saf" "${CHUNKDIR}/02_sfs" "${CHUNKDIR}/logs"

SAF_PREFIX="${CHUNKDIR}/01_saf/catfish.${TAG}"
SFS_OUT="${CHUNKDIR}/02_sfs/catfish.${TAG}.sfs"

angsd -b "${BAMLIST}" -ref "${REF}" \
  -GL 1 \
  -doSaf 5 -doMajorMinor 2 \
  -doCounts 1 -doDepth 1 \
  -minQ "${MINQ}" -minMapQ "${MINMAPQ}" -baq 1 -C 50 \
  -setMaxDepth "${SETMAXDEPTH}" \
  -minInd "${MININD}" \
  -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 \
  -rf "${RF_FILE}" \
  -sites "${CALLABLE_SITES}" \
  -P "${P}" \
  -out "${SAF_PREFIX}" \
  2> "${CHUNKDIR}/logs/01_angsd_saf.${TAG}.log"

realSFS "${SAF_PREFIX}.saf.idx" -bootstrap 100 -fold 1 -P "${P}" \
  > "${SFS_OUT}" \
  2> "${CHUNKDIR}/logs/02_realSFS.${TAG}.log"

echo "[DONE] ${TAG}"
