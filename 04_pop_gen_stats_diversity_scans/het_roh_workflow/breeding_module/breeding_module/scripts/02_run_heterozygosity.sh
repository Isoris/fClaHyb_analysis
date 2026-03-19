#!/usr/bin/env bash
# =============================================================================
# 02_run_heterozygosity.sh — Per-sample genome-wide het + local theta tracks
# =============================================================================
# For each QC-passing sample:
#   1) ANGSD -GL 1 -doSaf 1 -fold 1 -sites callable  → SAF
#   2) realSFS                                         → est.ml
#   3) Genome-wide het = SFS[1] / sum(SFS)
#   4) ANGSD rerun with -pest est.ml -doThetas 1
#   5) thetaStat make_bed + do_stat                    → local theta tracks
#   6) Optional: multiscale theta windows (5kb/1kb, 10kb/2kb, 50kb/10kb)
#
# NOTE: Local theta tracks are diversity proxies, NOT literal per-site Hobs.
#
# Usage:
#   bash 02_run_heterozygosity.sh              # run all samples
#   bash 02_run_heterozygosity.sh CGA009       # run one sample only
# =============================================================================
set -euo pipefail
SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPTDIR}/00_config.sh"
hr_init_dirs

SINGLE_SAMPLE="${1:-}"

hr_log "=== Step 02: Per-sample heterozygosity + local theta ==="
hr_check_file "${REF}" "reference"
hr_check_file "${CALLABLE_SITES}" "ANGSD callable sites"
hr_check_file "${CALLABLE_SITES}.idx" "ANGSD callable sites index"
hr_check_file "${CALLABLE_SITES}.bin" "ANGSD callable sites bin"
hr_check_file "${BAMLIST_QCPASS}" "BAM list"
hr_check_file "${SAMPLE_LIST}" "sample list"

hr_check_cmd angsd
hr_check_cmd realSFS
hr_check_cmd thetaStat

HET_SUMMARY="${DIR_HET}/04_summary/genomewide_heterozygosity.tsv"
if [[ ! -f "${HET_SUMMARY}" ]]; then
  echo -e "sample\thet_genomewide\tn_sites_sfs\tsfs_bin0\tsfs_bin1" > "${HET_SUMMARY}"
fi

if [[ "${RUN_EXTRA_THETA_SCALES:-0}" -eq 1 ]]; then
  mkdir -p "${DIR_HET}/03_theta/multiscale"
fi

# ── Build sample-to-BAM lookup ──────────────────────────────────────────────
declare -A SAMPLE_BAM
while IFS=$'\t' read -r SAMPLE BAM_MINI BAM_FILT; do
  [[ "${SAMPLE}" == "Sample" ]] && continue
  SAMPLE_BAM["${SAMPLE}"]="${BAM_FILT}"
done < "${SAMPLE_MANIFEST}"

# ── Process each sample ─────────────────────────────────────────────────────
while read -r SAMPLE; do
  [[ -z "${SAMPLE}" ]] && continue
  [[ "${SAMPLE}" =~ ^# ]] && continue

  if [[ -n "${SINGLE_SAMPLE}" && "${SAMPLE}" != "${SINGLE_SAMPLE}" ]]; then
    continue
  fi

  BAM="${SAMPLE_BAM[${SAMPLE}]:-}"
  if [[ -z "${BAM}" ]]; then
    hr_err "No BAM found for sample ${SAMPLE} in manifest, skipping"
    continue
  fi

  hr_log "Processing ${SAMPLE}..."

  if [[ -f "${DIR_HET}/03_theta/${SAMPLE}.thetas.idx" ]] && \
     grep -q "^${SAMPLE}	" "${HET_SUMMARY}" 2>/dev/null; then
    hr_log "  ${SAMPLE} already complete, skipping (delete outputs to rerun)"
    continue
  fi

  # ── Step 1: SAF ────────────────────────────────────────────────────────
  hr_log "  [1/5] SAF estimation..."
  angsd \
    -i "${BAM}" \
    -ref "${REF}" \
    -anc "${REF}" \
    -GL 1 \
    -doSaf 1 \
    -fold 1 \
    -minQ "${MINQ}" \
    -minMapQ "${MINMAPQ}" \
    -C "${CLIP}" \
    -sites "${CALLABLE_SITES}" \
    -nThreads "${THREADS}" \
    -out "${DIR_HET}/01_saf/${SAMPLE}" \
    2> "${DIR_LOGS}/${SAMPLE}.saf.log"

  # ── Step 2: realSFS → est.ml ───────────────────────────────────────────
  hr_log "  [2/5] realSFS..."
  realSFS \
    "${DIR_HET}/01_saf/${SAMPLE}.saf.idx" \
    -maxIter "${REALSFS_MAXITER}" \
    -tole "${REALSFS_TOLE}" \
    -P "${THREADS}" \
    > "${DIR_HET}/02_sfs/${SAMPLE}.est.ml" \
    2> "${DIR_LOGS}/${SAMPLE}.realsfs.log"

  # ── Step 3: Genome-wide het ────────────────────────────────────────────
  hr_log "  [3/5] Computing genome-wide heterozygosity..."
  HET=$(awk '{
    sum=0; for(i=1;i<=NF;i++) sum+=$i;
    if(sum>0) printf "%.12g\t%d\t%.6f\t%.6f", $2/sum, sum, $1, $2;
    else printf "NA\tNA\tNA\tNA";
  }' "${DIR_HET}/02_sfs/${SAMPLE}.est.ml")

  (
    flock -x 200
    if ! grep -q "^${SAMPLE}	" "${HET_SUMMARY}" 2>/dev/null; then
      echo -e "${SAMPLE}\t${HET}" >> "${HET_SUMMARY}"
    fi
  ) 200>"${HET_SUMMARY}.lock"

  # ── Step 4: Local theta with -pest ─────────────────────────────────────
  hr_log "  [4/5] Local theta estimation (diversity proxy)..."
  angsd \
    -i "${BAM}" \
    -ref "${REF}" \
    -anc "${REF}" \
    -GL 1 \
    -doSaf 1 \
    -fold 1 \
    -pest "${DIR_HET}/02_sfs/${SAMPLE}.est.ml" \
    -doThetas 1 \
    -minQ "${MINQ}" \
    -minMapQ "${MINMAPQ}" \
    -C "${CLIP}" \
    -sites "${CALLABLE_SITES}" \
    -nThreads "${THREADS}" \
    -out "${DIR_HET}/03_theta/${SAMPLE}" \
    2> "${DIR_LOGS}/${SAMPLE}.theta.log"

  # ── Step 5: thetaStat (main windows) ──────────────────────────────────
  hr_log "  [5/5] thetaStat..."

  thetaStat make_bed \
    "${DIR_HET}/03_theta/${SAMPLE}.thetas.idx" \
    -outnames "${DIR_HET}/03_theta/${SAMPLE}"

  thetaStat do_stat \
    "${DIR_HET}/03_theta/${SAMPLE}.thetas.idx" \
    -win "${WIN}" -step "${STEP}" \
    -outnames "${DIR_HET}/03_theta/${SAMPLE}.win${WIN}.step${STEP}"

  # ── Step 6: Multiscale theta windows ──────────────────────────────────
  if [[ "${RUN_EXTRA_THETA_SCALES:-0}" -eq 1 ]]; then
    hr_log "  [6/5] Extra theta scales..."
    MS_DIR="${DIR_HET}/03_theta/multiscale"

    for i in "${!THETA_SCALES[@]}"; do
      SCALE="${THETA_SCALES[$i]}"
      W="${SCALE%%_*}"
      S="${SCALE##*_}"
      thetaStat do_stat \
        "${DIR_HET}/03_theta/${SAMPLE}.thetas.idx" \
        -win "${W}" -step "${S}" \
        -outnames "${MS_DIR}/${SAMPLE}.win${W}.step${S}"
    done
  fi

  hr_log "  ${SAMPLE} done."

done < "${SAMPLE_LIST}"

# ── Post: Sort summary table ───────────────────────────────────────────────
if [[ -f "${HET_SUMMARY}" ]]; then
  (head -1 "${HET_SUMMARY}" && tail -n +2 "${HET_SUMMARY}" | sort -k1,1) \
    > "${HET_SUMMARY}.tmp" && mv "${HET_SUMMARY}.tmp" "${HET_SUMMARY}"
  rm -f "${HET_SUMMARY}.lock"
fi

hr_log "=== Step 02 complete ==="
hr_log "Summary: ${HET_SUMMARY}"
hr_log "Theta windows: ${DIR_HET}/03_theta/*.win${WIN}.step${STEP}.*"
if [[ "${RUN_EXTRA_THETA_SCALES:-0}" -eq 1 ]]; then
  hr_log "Extra scales:  ${DIR_HET}/03_theta/multiscale/"
fi
