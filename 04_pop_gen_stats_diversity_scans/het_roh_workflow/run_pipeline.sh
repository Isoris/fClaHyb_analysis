#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh — Master runner for the full het/ROH/FROH workflow
# =============================================================================
# Usage:
#   bash run_pipeline.sh                    # Run everything sequentially
#   bash run_pipeline.sh --step 2           # Run only step 2
#   bash run_pipeline.sh --from 3           # Run from step 3 onward
#   bash run_pipeline.sh --step 2 --slurm   # Submit step 2 as SLURM array
# =============================================================================
set -euo pipefail
SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)/scripts"

# Parse args
STEP=""
FROM=1
USE_SLURM=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --step)  STEP="$2"; shift 2 ;;
    --from)  FROM="$2"; shift 2 ;;
    --slurm) USE_SLURM=true; shift ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done

run_step() {
  local step_num="$1"
  local desc="$2"
  shift 2

  if [[ -n "${STEP}" && "${step_num}" != "${STEP}" ]]; then
    return 0
  fi
  if [[ "${step_num}" -lt "${FROM}" ]]; then
    return 0
  fi

  echo ""
  echo "================================================================"
  echo "  STEP ${step_num}: ${desc}"
  echo "================================================================"
  "$@"
}

# ── Step 1: Prepare inputs ──────────────────────────────────────────────
run_step 1 "Prepare and validate inputs" \
  bash "${SCRIPTDIR}/01_prep_inputs.sh"

# ── Step 2: Per-sample heterozygosity ───────────────────────────────────
if [[ "${USE_SLURM}" == true && ("${STEP}" == "2" || -z "${STEP}") ]]; then
  echo ""
  echo "================================================================"
  echo "  STEP 2: Per-sample heterozygosity (SLURM array)"
  echo "================================================================"
  source "${SCRIPTDIR}/00_config.sh"
  N=$(wc -l < "${SAMPLE_LIST}")
  echo "Submitting SLURM array for ${N} samples..."
  sbatch --array=1-${N} "${SCRIPTDIR}/02_run_heterozygosity_slurm.sh"
  echo "Submitted. Wait for completion before running step 3."
  if [[ -z "${STEP}" ]]; then
    echo "Stopping sequential run. Re-run with --from 3 after SLURM completes."
    exit 0
  fi
else
  run_step 2 "Per-sample heterozygosity (sequential)" \
    bash "${SCRIPTDIR}/02_run_heterozygosity.sh"
fi

# ── Step 3: ngsF-HMM ──────────────────────────────────────────────────
run_step 3 "ngsF-HMM (multi-replicate)" \
  bash "${SCRIPTDIR}/03_run_ngsF_HMM.sh"

# ── Step 4: Parse ROH + het in/out ROH ────────────────────────────────
run_step 4 "Parse ROH, compute FROH, het in/out ROH" \
  bash "${SCRIPTDIR}/04_parse_roh_and_het.sh"

# ── Step 5: Plots + stats + report ───────────────────────────────────
run_step 5 "Generate all plots, statistics, and report" \
  bash "${SCRIPTDIR}/05_run_all_plots.sh"

echo ""
echo "================================================================"
echo "  PIPELINE COMPLETE"
echo "================================================================"
