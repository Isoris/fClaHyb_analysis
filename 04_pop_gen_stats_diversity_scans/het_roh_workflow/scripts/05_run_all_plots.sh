#!/usr/bin/env bash
# v2 note: passes per_chr + breeding_dir to run_stats.R; adds multiscale theta ideogram loop
# =============================================================================
# 05_run_all_plots.sh — Run all plotting + statistics scripts
# =============================================================================
set -euo pipefail
SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPTDIR}/00_config.sh"
hr_init_dirs

hr_log "=== Step 05: Run all plots and statistics ==="

Rscript -e 'for(p in c("data.table","ggplot2")) if(!require(p, character.only=TRUE)) stop(paste("Missing R package:", p))' || {
  hr_die "Missing required R packages: data.table, ggplot2"
}

MASTER="${DIR_TABLES}/master_summary.tsv"
PER_CHR="${DIR_TABLES}/per_chr_roh_summary.tsv"
BINS="${DIR_TABLES}/catfish_roh.per_sample_roh_bins_long.tsv"
HET_TSV="${DIR_HET}/04_summary/genomewide_heterozygosity.tsv"
CHROM_SIZES="${DIR_INPUTS}/chrom_sizes.tsv"

hr_check_file "${MASTER}" "master_summary.tsv"

# ── 1. Core heterozygosity plots ────────────────────────────────────────
hr_log "Running plot_heterozygosity_core.R..."
if [[ -f "${HET_TSV}" ]]; then
  Rscript "${SCRIPTDIR}/plot_heterozygosity_core.R" \
    "${HET_TSV}" "${DIR_PLOTS_CORE}" \
    ${ANCESTRY_LABELS:+"${ANCESTRY_LABELS}"} \
    2>&1 | tee "${DIR_LOGS}/plot_het_core.log"
fi

# ── 2. Core ROH plots ──────────────────────────────────────────────────
hr_log "Running plot_roh_core.R..."
if [[ -f "${BINS}" ]]; then
  Rscript "${SCRIPTDIR}/plot_roh_core.R" \
    "${MASTER}" "${BINS}" "${DIR_PLOTS_CORE}" \
    2>&1 | tee "${DIR_LOGS}/plot_roh_core.log"
else
  Rscript "${SCRIPTDIR}/plot_roh_core.R" \
    "${MASTER}" "/dev/null" "${DIR_PLOTS_CORE}" \
    2>&1 | tee "${DIR_LOGS}/plot_roh_core.log" || true
fi

# ── 3. Scatter plots ───────────────────────────────────────────────────
hr_log "Running plot_scatter_stats.R..."
Rscript "${SCRIPTDIR}/plot_scatter_stats.R" \
  "${MASTER}" "${DIR_PLOTS_CORE}" \
  ${SAMPLE_QC_TABLE:+""} \
  ${ANCESTRY_LABELS:+"${ANCESTRY_LABELS}"} \
  2>&1 | tee "${DIR_LOGS}/plot_scatter.log"

# ── 4. Chromosome-level plots ─────────────────────────────────────────
hr_log "Running plot_roh_by_chromosome.R..."
if [[ -f "${PER_CHR}" ]]; then
  Rscript "${SCRIPTDIR}/plot_roh_by_chromosome.R" \
    "${PER_CHR}" "${DIR_PLOTS_CORE}" \
    ${CHROM_SIZES:+"${CHROM_SIZES}"} \
    ${ANCESTRY_LABELS:+"${ANCESTRY_LABELS}"} \
    2>&1 | tee "${DIR_LOGS}/plot_chr.log"
fi

# ── 5. Theta ideogram plots (main + multiscale) ─────────────────────
hr_log "Running plot_theta_ideogram.R (main scale)..."
if [[ -d "${DIR_HET}/03_theta" && -f "${CHROM_SIZES}" ]]; then
  Rscript "${SCRIPTDIR}/plot_theta_ideogram.R" \
    "${DIR_HET}/03_theta" "${CHROM_SIZES}" "${DIR_PLOTS_CORE}" \
    "${SAMPLE_LIST}" 10 \
    2>&1 | tee "${DIR_LOGS}/plot_ideogram_main.log"

  # Multiscale theta ideograms
  if [[ "${RUN_EXTRA_THETA_SCALES:-0}" -eq 1 && -d "${DIR_HET}/03_theta/multiscale" ]]; then
    for i in "${!THETA_SCALE_LABELS[@]}"; do
      LABEL="${THETA_SCALE_LABELS[$i]}"
      hr_log "  Theta ideogram: ${LABEL}..."
      MS_OUT="${DIR_PLOTS_CORE}/theta_${LABEL}"
      mkdir -p "${MS_OUT}"
      Rscript "${SCRIPTDIR}/plot_theta_ideogram.R" \
        "${DIR_HET}/03_theta/multiscale" "${CHROM_SIZES}" "${MS_OUT}" \
        "${SAMPLE_LIST}" 10 \
        2>&1 | tee "${DIR_LOGS}/plot_ideogram_${LABEL}.log" || true
    done
  fi
fi

# ── 6. Metadata overlay plots ────────────────────────────────────────
hr_log "Running plot_roh_metadata_overlays.R..."
META_ARGS=()
if [[ -n "${ANCESTRY_LABELS}" && -f "${ANCESTRY_LABELS}" ]]; then
  META_ARGS+=(--ancestry "${ANCESTRY_LABELS}")
fi
if [[ -n "${COVTREE_ORDER}" && -f "${COVTREE_ORDER}" ]]; then
  META_ARGS+=(--order "${COVTREE_ORDER}")
fi
if [[ -n "${PRUNED81_SAMPLES}" && -f "${PRUNED81_SAMPLES}" ]]; then
  META_ARGS+=(--pruned81 "${PRUNED81_SAMPLES}")
fi
if [[ -n "${NGSRELATE_PAIRS}" && -f "${NGSRELATE_PAIRS}" ]]; then
  META_ARGS+=(--kinship "${NGSRELATE_PAIRS}")
fi

if [[ -f "${PER_CHR}" ]]; then
  Rscript "${SCRIPTDIR}/plot_roh_metadata_overlays.R" \
    "${MASTER}" "${PER_CHR}" "${BINS:-/dev/null}" "${DIR_PLOTS_META}" \
    "${META_ARGS[@]}" \
    2>&1 | tee "${DIR_LOGS}/plot_metadata.log"
fi

# ── 7. Statistics (v2: now with chr-level and theta-scale stats) ─────
hr_log "Running run_stats.R..."
STATS_ARGS=("${MASTER}" "${DIR_STATS}")
if [[ -n "${ANCESTRY_LABELS}" && -f "${ANCESTRY_LABELS}" ]]; then
  STATS_ARGS+=("${ANCESTRY_LABELS}")
else
  STATS_ARGS+=("")
fi
if [[ -f "${PER_CHR}" ]]; then
  STATS_ARGS+=("${PER_CHR}")
else
  STATS_ARGS+=("")
fi
if [[ -d "${DIR_BREEDING}" ]]; then
  STATS_ARGS+=("${DIR_BREEDING}")
else
  STATS_ARGS+=("")
fi

Rscript "${SCRIPTDIR}/run_stats.R" "${STATS_ARGS[@]}" \
  2>&1 | tee "${DIR_LOGS}/run_stats.log"

# ── 8. Report ────────────────────────────────────────────────────────
hr_log "Running 06_write_report.py..."
python3 "${SCRIPTDIR}/06_write_report.py" \
  --tables-dir "${DIR_TABLES}" \
  --stats-dir "${DIR_STATS}" \
  --out-dir "${DIR_REPORT}" \
  --ngsf-reps "${NGSFHMM_REPS}" \
  2>&1 | tee "${DIR_LOGS}/write_report.log"

hr_log "=== Step 05 complete ==="
