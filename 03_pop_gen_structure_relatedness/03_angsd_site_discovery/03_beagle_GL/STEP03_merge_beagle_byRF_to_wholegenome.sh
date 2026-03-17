#!/usr/bin/env bash
set -euo pipefail

###
# STEP03_merge_beagle_byRF_to_wholegenome.sh
#
# Does this:
# - Merges per-RF BEAGLE files into one whole-genome BEAGLE file
# - Keeps the header only once
# - Appends all data rows in sorted RF order
# - Writes and keeps one .list file per thinning level for provenance/reuse
# - Validates gzip output
#
# Expected input layout:
#   ${BASE}/popstruct_beagle_GL/01_beagle_byRF/thin_${W}/catfish.*.thin_${W}.beagle.gz
#
# Output layout:
#   ${BASE}/popstruct_beagle_GL/03_merged_beagle/
#     - beagle_thin_${W}.list
#     - catfish.wholegenome.byRF.thin_${W}.beagle.gz
#
# Usage:
#   bash STEP03_merge_beagle_byRF_to_wholegenome.sh
#   bash STEP03_merge_beagle_byRF_to_wholegenome.sh 200 500 1000
#   bash STEP03_merge_beagle_byRF_to_wholegenome.sh 5000 10000 25000
###

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
OUTDIR="${BASE}/popstruct_beagle_GL"

BYRF_DIR="${OUTDIR}/01_beagle_byRF"
MERGED_DIR="${OUTDIR}/03_merged_beagle"
LOGDIR="${OUTDIR}/logs"

mkdir -p "${MERGED_DIR}" "${LOGDIR}"

# Default thinning values if none provided
if [[ "$#" -eq 0 ]]; then
  THIN_LIST=(200 500 1000 5000 10000 25000)
else
  THIN_LIST=("$@")
fi

for W in "${THIN_LIST[@]}"; do
  INDIR="${BYRF_DIR}/thin_${W}"
  [[ -d "${INDIR}" ]] || {
    echo "[ERROR] Missing byRF directory: ${INDIR}" >&2
    exit 1
  }

  LIST="${MERGED_DIR}/beagle_thin_${W}.list"
  TMP="${MERGED_DIR}/catfish.wholegenome.byRF.thin_${W}.tmp.beagle"
  OUT="${MERGED_DIR}/catfish.wholegenome.byRF.thin_${W}.beagle.gz"

  echo "[INFO] thin=${W}"
  echo "[INFO] input dir : ${INDIR}"
  echo "[INFO] list file : ${LIST}"
  echo "[INFO] output    : ${OUT}"

  # Keep the list file as a real output
  find "${INDIR}" -maxdepth 1 -type f -name "*.thin_${W}.beagle.gz" | sort -V > "${LIST}"

  N=$(wc -l < "${LIST}")
  if [[ "${N}" -eq 0 ]]; then
    echo "[ERROR] No beagle.gz files found for thin=${W} in ${INDIR}" >&2
    exit 1
  fi

  FIRST=$(head -n 1 "${LIST}")
  [[ -n "${FIRST}" ]] || {
    echo "[ERROR] Empty first file for thin=${W}" >&2
    exit 1
  }

  # Write header once
  zcat "${FIRST}" | head -n 1 > "${TMP}"

  # Append all data lines without header
  while read -r F; do
    [[ -s "${F}" ]] || {
      echo "[ERROR] Missing or empty beagle file listed in ${LIST}: ${F}" >&2
      exit 1
    }
    zcat "${F}" | tail -n +2 >> "${TMP}"
  done < "${LIST}"

  gzip -c "${TMP}" > "${OUT}"
  rm -f "${TMP}"

  gzip -t "${OUT}"

  NLINES=$(zcat "${OUT}" | wc -l)

  echo "[OK] wrote ${OUT}"
  echo "[OK] kept  ${LIST}"
  echo "[OK] files merged: ${N}"
  echo "[OK] total lines : ${NLINES}"
  echo
done

echo "[DONE] All requested thinning levels merged successfully."
