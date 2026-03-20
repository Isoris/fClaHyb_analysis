#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=64G
#SBATCH -t 0-04:00:00
#SBATCH -J delly_plot
#SBATCH -o logs/06_plots.%j.out
#SBATCH -e logs/06_plots.%j.err
set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/00_delly_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a
source "${CONFIG}"
set +a

dv_init_dirs
dv_log "=== STEP 6: Publication-style DEL plots ==="

# Check required intermediate tables from step 5
for f in \
  "${DIR_SUMMARY}/per_sample_DEL_counts.tsv" \
  "${DIR_SUMMARY}/DEL_svlen_distribution.tsv" \
  "${DIR_SUMMARY}/per_chromosome_DEL_counts.tsv" \
  "${DIR_SUMMARY}/DEL_window_counts_1Mb.tsv" \
  "${DIR_SUMMARY}/pairwise_shared_DEL.tsv" \
  "${DIR_SUMMARY}/DEL_binary_genotype_matrix.tsv" \
  "${DIR_FINAL}/catalog_226.DEL.GT_matrix.tsv"; do
  dv_check_file "$f" "$(basename "$f")"
done

# Run R plotting script
Rscript "${SCRIPT_DIR}/06_plot_delly_results.R" \
  --summary_dir "${DIR_SUMMARY}" \
  --final_dir   "${DIR_FINAL}" \
  --plot_dir    "${DIR_PLOTS}" \
  --ref_fai     "${REF_FAI}" \
  --samples_unrelated "${SAMPLES_UNRELATED}"

dv_log "=== PLOTS COMPLETE ==="
dv_log "Output: ${DIR_PLOTS}/"
ls -lh "${DIR_PLOTS}"/*.pdf "${DIR_PLOTS}"/*.png 2>/dev/null | while read line; do dv_log "  ${line}"; done
