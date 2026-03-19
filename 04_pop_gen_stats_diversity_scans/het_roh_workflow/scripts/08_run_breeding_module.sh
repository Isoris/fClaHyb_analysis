#!/usr/bin/env bash
# v2 note: fixed array argument passing — uses "${ARRAY[@]}" directly in CMD
# =============================================================================
# 08_run_breeding_module.sh — Build all breeding tables + plots
# =============================================================================
set -euo pipefail
SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPTDIR}/00_config.sh"
hr_init_dirs

MAX_PAIRS="${1:-0}"

hr_log "=== Step 08: Build breeding tables + plots ==="

TRACTS="${DIR_TABLES}/roh_tracts_all.bed"
CHROM_SIZES="${DIR_INPUTS}/chrom_sizes.tsv"

hr_check_file "${TRACTS}" "ROH tracts BED"
hr_check_file "${CALLABLE_BED}" "callable BED"
hr_check_file "${CHROM_SIZES}" "chrom sizes"
hr_check_file "${SAMPLE_LIST}" "sample list"

# ── 1. Build breeding tables ────────────────────────────────────────────
hr_log "Building breeding tables..."
CMD=(
  python3 "${SCRIPTDIR}/07_build_breeding_tables.py"
  --tracts-bed "${TRACTS}"
  --callable-bed "${CALLABLE_BED}"
  --chrom-sizes "${CHROM_SIZES}"
  --theta-dir "${DIR_HET}/03_theta/multiscale"
  --theta-main-dir "${DIR_HET}/03_theta"
  --sample-list "${SAMPLE_LIST}"
  --out-dir "${DIR_BREEDING}"
  --max-pairs "${MAX_PAIRS}"
  --theta-scales "${THETA_SCALES[@]}"
  --theta-labels "${THETA_SCALE_LABELS[@]}"
  --segment-sizes "${SEGMENT_SIZES[@]}"
  --segment-labels "${SEGMENT_LABELS[@]}"
  --recurrence-sizes "${RECURRENCE_WINDOWS[@]}"
  --recurrence-labels "${RECURRENCE_LABELS[@]}"
)

if [[ -n "${ANCESTRY_LABELS}" && -f "${ANCESTRY_LABELS}" ]]; then
  CMD+=(--ancestry-labels "${ANCESTRY_LABELS}")
fi
if [[ -n "${COVTREE_ORDER}" && -f "${COVTREE_ORDER}" ]]; then
  CMD+=(--order-file "${COVTREE_ORDER}")
fi
if [[ -n "${PRUNED81_SAMPLES}" && -f "${PRUNED81_SAMPLES}" ]]; then
  CMD+=(--pruned81 "${PRUNED81_SAMPLES}")
fi

"${CMD[@]}" 2>&1 | tee "${DIR_LOGS}/07_breeding_tables.log"

# ── 2. ROH genome map + scatter ─────────────────────────────────────────
hr_log "Plotting ROH genome map..."
FIXED_BINS_FILE="${DIR_BREEDING}/per_sample_roh_bins_wide_fixedBins.tsv"
ADAPT_BINS_FILE="${DIR_BREEDING}/per_sample_roh_bins_wide_adaptiveBins.tsv"
SEG_FILE="${DIR_BREEDING}/per_sample_genome_roh_segments.tsv"

if [[ -f "${SEG_FILE}" && -f "${FIXED_BINS_FILE}" ]]; then
  Rscript "${SCRIPTDIR}/plot_roh_genome_map.R" \
    "${SEG_FILE}" \
    "${FIXED_BINS_FILE}" \
    "${ADAPT_BINS_FILE}" \
    "${CHROM_SIZES}" \
    "${DIR_PLOTS_BREEDING}" \
    ${COVTREE_ORDER:+"${COVTREE_ORDER}"} \
    ${ANCESTRY_LABELS:+"${ANCESTRY_LABELS}"} \
    2>&1 | tee "${DIR_LOGS}/plot_genome_map.log"
fi

# ── 3. ROH recurrence plots ─────────────────────────────────────────────
hr_log "Plotting ROH recurrence..."
Rscript "${SCRIPTDIR}/plot_roh_recurrence.R" \
  "${DIR_BREEDING}" \
  "${DIR_PLOTS_BREEDING}" \
  2>&1 | tee "${DIR_LOGS}/plot_recurrence.log"

# ── 4. Breeding heatmaps ────────────────────────────────────────────────
hr_log "Plotting breeding heatmaps..."
Rscript "${SCRIPTDIR}/plot_breeding_heatmaps.R" \
  "${DIR_BREEDING}" \
  "${DIR_PLOTS_BREEDING}" \
  ${COVTREE_ORDER:+"${COVTREE_ORDER}"} \
  2>&1 | tee "${DIR_LOGS}/plot_breeding_heatmaps.log"

hr_log "=== Step 08 complete ==="
hr_log "Breeding tables: ${DIR_BREEDING}"
hr_log "Breeding plots:  ${DIR_PLOTS_BREEDING}"
