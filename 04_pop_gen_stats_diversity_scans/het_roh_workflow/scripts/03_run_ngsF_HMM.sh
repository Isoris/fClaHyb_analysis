#!/usr/bin/env bash
# =============================================================================
# 03_run_ngsF_HMM.sh — Run ngsF-HMM with multiple random starts, keep best
# =============================================================================
# Uses explicit binary path from config; loads required modules.
#
# Usage:
#   bash 03_run_ngsF_HMM.sh
#   bash 03_run_ngsF_HMM.sh 20    # override number of replicates
# =============================================================================
set -euo pipefail
SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPTDIR}/00_config.sh"
hr_init_dirs

# ── Load required modules ──────────────────────────────────────────────────
for mod in "${NGSFHMM_MODULES[@]}"; do
  module load "${mod}" 2>/dev/null || hr_log "WARNING: could not load module ${mod}"
done

NREPS="${1:-${NGSFHMM_REPS}}"

hr_log "=== Step 03: ngsF-HMM (${NREPS} replicates) ==="

# Validate inputs
hr_check_file "${BEAGLE_GL}" "BEAGLE GL"
hr_check_file "${POS_FILE}" "position file"
hr_check_file "${SAMPLES_IND}" "samples.ind"
hr_check_file "${NGSFHMM_BIN}" "ngsF-HMM binary"

# Count individuals and sites
N_IND=$(wc -l < "${SAMPLES_IND}")
N_SITES=$(wc -l < "${POS_FILE}")
hr_log "N individuals: ${N_IND}"
hr_log "N sites: ${N_SITES}"
hr_log "ngsF-HMM binary: ${NGSFHMM_BIN}"

# Create output directories
mkdir -p "${DIR_NGSFHMM}/replicates" "${DIR_NGSFHMM}/best"

# ── Run replicates ──────────────────────────────────────────────────────────
LL_SUMMARY="${DIR_NGSFHMM}/replicate_loglikelihoods.tsv"
echo -e "replicate\tseed\tlog_likelihood\tconverged" > "${LL_SUMMARY}"

for REP in $(seq 1 "${NREPS}"); do
  SEED=$(( NGSFHMM_SEED_BASE + REP ))
  PREFIX="${DIR_NGSFHMM}/replicates/rep${REP}"

  hr_log "  Replicate ${REP}/${NREPS} (seed=${SEED})..."

  # Check if already done
  if [[ -f "${PREFIX}.indF" && -f "${PREFIX}.ibd" ]]; then
    hr_log "    Already complete, skipping"
  else
    "${NGSFHMM_BIN}" \
      --geno "${BEAGLE_GL}" \
      --pos "${POS_FILE}" \
      --n_ind "${N_IND}" \
      --n_sites "${N_SITES}" \
      --lkl \
      --seed "${SEED}" \
      --out "${PREFIX}" \
      2> "${DIR_LOGS}/ngsF-HMM_rep${REP}.log" || {
        hr_err "  Replicate ${REP} failed — see ${DIR_LOGS}/ngsF-HMM_rep${REP}.log"
        echo -e "${REP}\t${SEED}\tNA\tFAILED" >> "${LL_SUMMARY}"
        continue
      }
  fi

  # Extract log-likelihood from log
  LL=$(grep -i "likelihood\|lnL\|loglik" "${DIR_LOGS}/ngsF-HMM_rep${REP}.log" \
       | tail -1 | grep -oP '[-]?[\d.]+(?:e[+-]?\d+)?' | tail -1 || echo "NA")
  CONV="yes"
  if [[ "${LL}" == "NA" ]]; then
    LL=$(tail -1 "${DIR_LOGS}/ngsF-HMM_rep${REP}.log" \
         | grep -oP '[-]?[\d.]+(?:e[+-]?\d+)?' | tail -1 || echo "NA")
    CONV="unknown"
  fi

  echo -e "${REP}\t${SEED}\t${LL}\t${CONV}" >> "${LL_SUMMARY}"
  hr_log "    LL = ${LL}"

done

# ── Select best replicate ──────────────────────────────────────────────────
hr_log "Selecting best replicate..."

BEST_REP=$(awk -F'\t' 'NR>1 && $3!="NA" {
  if (best=="" || $3+0 > best_ll+0) { best=$1; best_ll=$3 }
} END { print best }' "${LL_SUMMARY}")

if [[ -z "${BEST_REP}" ]]; then
  hr_die "No successful replicates found!"
fi

BEST_LL=$(awk -F'\t' -v r="${BEST_REP}" 'NR>1 && $1==r {print $3}' "${LL_SUMMARY}")
hr_log "Best replicate: ${BEST_REP} (LL = ${BEST_LL})"

# Copy best results
for EXT in indF ibd geno; do
  SRC="${DIR_NGSFHMM}/replicates/rep${BEST_REP}.${EXT}"
  if [[ -f "${SRC}" ]]; then
    cp "${SRC}" "${DIR_NGSFHMM}/best/best.${EXT}"
    hr_log "  Copied best .${EXT}"
  else
    hr_err "  Missing: ${SRC}"
  fi
done

# Record best info
echo -e "best_replicate\t${BEST_REP}" > "${DIR_NGSFHMM}/best/best_info.txt"
echo -e "best_seed\t$(( NGSFHMM_SEED_BASE + BEST_REP ))" >> "${DIR_NGSFHMM}/best/best_info.txt"
echo -e "best_loglikelihood\t${BEST_LL}" >> "${DIR_NGSFHMM}/best/best_info.txt"

hr_log "=== Step 03 complete ==="
hr_log "Best results: ${DIR_NGSFHMM}/best/"
hr_log "LL summary: ${LL_SUMMARY}"
