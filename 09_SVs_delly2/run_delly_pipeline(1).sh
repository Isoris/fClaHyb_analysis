#!/usr/bin/env bash
# =============================================================================
# run_delly_pipeline.sh — Submit DELLY DEL pipeline to SLURM with dependencies
# =============================================================================
# Usage:
#   cd /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/delly_sv
#   bash run_delly_pipeline.sh
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${SCRIPT_DIR}/logs"
mkdir -p "${LOG_DIR}"

echo "=== DELLY DEL Pipeline Launcher ==="
echo "Script directory: ${SCRIPT_DIR}"
echo ""

# Verify DELLY binary exists
DELLY_BIN="/project/lt200308-agbsci/01-catfish_assembly/delly/src/delly"
if [[ ! -x "${DELLY_BIN}" ]]; then
  echo "ERROR: DELLY binary not found: ${DELLY_BIN}"
  exit 1
fi
echo "DELLY binary: ${DELLY_BIN}"
echo ""

echo "Submitting 5-job chain..."
echo ""

JID1=$(sbatch --parsable \
  -o "${LOG_DIR}/01_prep.%j.out" -e "${LOG_DIR}/01_prep.%j.err" \
  "${SCRIPT_DIR}/01_prep_inputs.sh")
echo "  [1/5] Prep inputs:      Job ${JID1}"

JID2=$(sbatch --parsable --dependency=afterok:${JID1} \
  -o "${LOG_DIR}/02_discovery.%j.out" -e "${LOG_DIR}/02_discovery.%j.err" \
  "${SCRIPT_DIR}/02_delly_discovery.sh")
echo "  [2/5] DEL discovery:     Job ${JID2} (after ${JID1})"

JID3=$(sbatch --parsable --dependency=afterok:${JID2} \
  -o "${LOG_DIR}/03_merge_geno.%j.out" -e "${LOG_DIR}/03_merge_geno.%j.err" \
  "${SCRIPT_DIR}/03_delly_merge_genotype.sh")
echo "  [3/5] Merge+genotype:    Job ${JID3} (after ${JID2})"

JID4=$(sbatch --parsable --dependency=afterok:${JID3} \
  -o "${LOG_DIR}/04_annotation.%j.out" -e "${LOG_DIR}/04_annotation.%j.err" \
  "${SCRIPT_DIR}/04_annotation_layers.sh")
echo "  [4/5] Annotation:        Job ${JID4} (after ${JID3})"

JID5=$(sbatch --parsable --dependency=afterok:${JID4} \
  -o "${LOG_DIR}/05_summary.%j.out" -e "${LOG_DIR}/05_summary.%j.err" \
  "${SCRIPT_DIR}/05_summary_report.sh")
echo "  [5/5] Summary:           Job ${JID5} (after ${JID4})"

echo ""
echo "=== Pipeline submitted ==="
echo "Chain: ${JID1} → ${JID2} → ${JID3} → ${JID4} → ${JID5}"
echo ""
echo "Monitor:"
echo "  squeue -u \$(whoami)"
echo "  tail -f ${LOG_DIR}/*.out"
echo ""
echo "Estimated wall time:"
echo "  01 prep:       ~15 min"
echo "  02 discovery:  ~6-24h (226 samples × DEL, 30 parallel, 4 threads each)"
echo "  03 merge+geno: ~6-24h (merge + 226 regenotype)"
echo "  04 annotation: ~2-6h"
echo "  05 summary:    ~10 min"
