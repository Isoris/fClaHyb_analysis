#!/usr/bin/env bash
# =============================================================================
# 08_run_breeding_module.sh — Build all breeding tables + plots
# =============================================================================
# Runs after steps 01-05 are complete. Reads existing outputs and produces
# multiscale derived tables and breeding-ready plots.
#
# Usage:  bash 08_run_breeding_module.sh [--max-pairs N]
# =============================================================================
set -euo pipefail
SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPTDIR}/00_config.sh"
hr_init_dirs

MAX_PAIRS="${1:-0}"  # 0 = all pairs; use e.g. 1000 for testing

hr_log "=== Step 08: Build breeding tables + plots ==="

# ── Required inputs ──────────────────────────────────────────────────────
TRACTS="${DIR_TABLES}/roh_tracts_all.bed"
CHROM_SIZES="${DIR_INPUTS}/chrom_sizes.tsv"

hr_check_file "${TRACTS}" "ROH tracts BED"
hr_check_file "${CALLABLE_BED}" "callable BED"
hr_check_file "${CHROM_SIZES}" "chrom sizes"
hr_check_file "${SAMPLE_LIST}" "sample list"

# ── Build theta scale args ──────────────────────────────────────────────
THETA_SCALE_ARGS=""
for s in "${THETA_SCALES[@]}"; do
  THETA_SCALE_ARGS="${THETA_SCALE_ARGS} ${s}"
done
THETA_LABEL_ARGS=""
for s in "${THETA_SCALE_LABELS[@]}"; do
  THETA_LABEL_ARGS="${THETA_LABEL_ARGS} ${s}"
done
SEG_SIZE_ARGS=""
for s in "${SEGMENT_SIZES[@]}"; do
  SEG_SIZE_ARGS="${SEG_SIZE_ARGS} ${s}"
done
SEG_LABEL_ARGS=""
for s in "${SEGMENT_LABELS[@]}"; do
  SEG_LABEL_ARGS="${SEG_LABEL_ARGS} ${s}"
done
REC_SIZE_ARGS=""
for s in "${RECURRENCE_WINDOWS[@]}"; do
  REC_SIZE_ARGS="${REC_SIZE_ARGS} ${s}"
done
REC_LABEL_ARGS=""
for s in "${RECURRENCE_LABELS[@]}"; do
  REC_LABEL_ARGS="${REC_LABEL_ARGS} ${s}"
done

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
  --theta-labels "${THETA_LABELS[@]}"
  --segment-sizes "${SEGMENT_SIZES[@]}"
  --segment-labels "${SEGMENT_LABELS[@]}"
  --recurrence-sizes "${RECURRENCE_WINDOWS[@]}"
  --recurrence-labels "${RECURRENCE_LABELS[@]}"
)

#CMD=(
#  python3 "${SCRIPTDIR}/07_build_breeding_tables.py"
#  --tracts-bed "${TRACTS}"
#  --callable-bed "${CALLABLE_BED}"
#  --chrom-sizes "${CHROM_SIZES}"
#  --theta-dir "${DIR_HET}/03_theta/multiscale"
#  --theta-main-dir "${DIR_HET}/03_theta"
#  --sample-list "${SAMPLE_LIST}"
#  --out-dir "${DIR_BREEDING}"
#  --max-pairs "${MAX_PAIRS}"
#  --theta-scales ${THETA_SCALE_ARGS}
#  --theta-labels ${THETA_LABEL_ARGS}
#  --segment-sizes ${SEG_SIZE_ARGS}
#  --segment-labels ${SEG_LABEL_ARGS}
#  --recurrence-sizes ${REC_SIZE_ARGS}
#  --recurrence-labels ${REC_LABEL_ARGS}
#)

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
