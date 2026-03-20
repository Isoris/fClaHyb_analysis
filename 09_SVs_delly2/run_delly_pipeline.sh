#!/usr/bin/env bash
# =============================================================================
# run_delly_pipeline.sh — Submit DELLY DEL pipeline to SLURM with dependencies
# =============================================================================
# Usage:
#   cd /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/delly_sv
#   bash run_delly_pipeline.sh
# =============================================================================
set -euo pipefail

# ── Config loading (same pattern as SLURM scripts) ─────────────────────────
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/00_delly_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }

set -a
source "${CONFIG}"
set +a

# ── Create log directory ───────────────────────────────────────────────────
LOG_DIR="${OUTDIR}/logs"
mkdir -p "${LOG_DIR}"

echo "=== DELLY DEL Pipeline Launcher ==="
echo "Script directory: ${SCRIPT_DIR}"
echo "DELLY binary:    ${DELLY_BIN}"
echo "BAM manifest:    ${BAM_MANIFEST}"
echo ""

# Verify DELLY binary
[[ -x "${DELLY_BIN}" ]] || { echo "ERROR: DELLY not found: ${DELLY_BIN}" >&2; exit 1; }

# Verify manifest exists
[[ -f "${BAM_MANIFEST}" ]] || { echo "ERROR: Manifest not found: ${BAM_MANIFEST}" >&2; exit 1; }

echo "Submitting 6-job dependency chain..."
echo ""

JID1=$(sbatch --parsable \
  -o "${LOG_DIR}/01_prep.%j.out" -e "${LOG_DIR}/01_prep.%j.err" \
  "${SCRIPT_DIR}/01_prep_inputs.sh")
echo "  [1/6] Prep inputs:       Job ${JID1}"

JID2=$(sbatch --parsable --dependency=afterok:${JID1} \
  -o "${LOG_DIR}/02_discovery.%j.out" -e "${LOG_DIR}/02_discovery.%j.err" \
  "${SCRIPT_DIR}/02_delly_discovery.sh")
echo "  [2/6] DEL discovery:     Job ${JID2} (after ${JID1})"

JID3=$(sbatch --parsable --dependency=afterok:${JID2} \
  -o "${LOG_DIR}/03_merge_geno.%j.out" -e "${LOG_DIR}/03_merge_geno.%j.err" \
  "${SCRIPT_DIR}/03_delly_merge_genotype.sh")
echo "  [3/6] Merge+genotype:    Job ${JID3} (after ${JID2})"

JID4=$(sbatch --parsable --dependency=afterok:${JID3} \
  -o "${LOG_DIR}/04_annotation.%j.out" -e "${LOG_DIR}/04_annotation.%j.err" \
  "${SCRIPT_DIR}/04_annotation_layers.sh")
echo "  [4/6] Annotation:        Job ${JID4} (after ${JID3})"

JID5=$(sbatch --parsable --dependency=afterok:${JID4} \
  -o "${LOG_DIR}/05_summary.%j.out" -e "${LOG_DIR}/05_summary.%j.err" \
  "${SCRIPT_DIR}/05_summary_report.sh")
echo "  [5/6] Summary:           Job ${JID5} (after ${JID4})"

JID6=$(sbatch --parsable --dependency=afterok:${JID5} \
  -o "${LOG_DIR}/06_plots.%j.out" -e "${LOG_DIR}/06_plots.%j.err" \
  "${SCRIPT_DIR}/06_plot_results.sh")
echo "  [6/6] Plots:             Job ${JID6} (after ${JID5})"

echo ""
echo "=== Pipeline submitted ==="
echo "Chain: ${JID1} → ${JID2} → ${JID3} → ${JID4} → ${JID5} → ${JID6}"
echo ""
echo "Monitor:  squeue -u \$(whoami)"
echo "Logs:     tail -f ${LOG_DIR}/*.out"
