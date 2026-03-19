#!/usr/bin/env bash
#SBATCH --job-name=het_angsd
#SBATCH --output=logs/het_%A_%a.out
#SBATCH --error=logs/het_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition=compute
# =============================================================================
# 02_run_heterozygosity_slurm.sh — SLURM array wrapper
# =============================================================================
# Submit as:
#   N=$(wc -l < /path/to/het_roh/01_inputs_check/samples_qcpass.txt)
#   sbatch --array=1-${N} 02_run_heterozygosity_slurm.sh
# =============================================================================
set -euo pipefail
SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPTDIR}/00_config.sh"

# Get the sample for this array task
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")

if [[ -z "${SAMPLE}" ]]; then
  echo "No sample at line ${SLURM_ARRAY_TASK_ID}, exiting."
  exit 0
fi

hr_log "SLURM array task ${SLURM_ARRAY_TASK_ID}: processing ${SAMPLE}"

# Run the main heterozygosity script for this single sample
bash "${SCRIPTDIR}/02_run_heterozygosity.sh" "${SAMPLE}"
