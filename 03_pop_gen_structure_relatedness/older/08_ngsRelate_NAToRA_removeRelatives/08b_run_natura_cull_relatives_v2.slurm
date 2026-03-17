###
# 08b_run_natura_cull_relatives_v2.slurm does this:
# Converts ngsRelate output into NAToRA input, runs NAToRA across multiple
# theta cutoffs for different relatedness classes, and writes a summary table
# of removed vs retained individuals.
###
#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1 -n 16
#SBATCH --mem=64GB
#SBATCH -t 08:00:00
#SBATCH -J natora_theta_compare
#SBATCH -o natora_theta_compare.%j.out
#SBATCH -e natora_theta_compare.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
OUTDIR="${BASE}/popstruct_thin/05_ngsrelate"
NATORA="${BASE}/NAToRA_Public/NAToRA_Public.py"
RES="${OUTDIR}/catfish_226_relatedness.res"
INPUT="${OUTDIR}/catfish_226_for_natora.txt"
SUMMARY="${OUTDIR}/catfish_theta_cutoff_summary.tsv"

mkdir -p "${OUTDIR}"

if [[ ! -s "${RES}" ]]; then
  echo "[ERROR] Missing ngsRelate result: ${RES}"
  exit 1
fi

# Build NAToRA input from ngsRelate:
# col3 = ida, col4 = idb, col18 = theta
awk 'BEGIN{OFS="\t"} NR>1 {print $3, $4, $18}' "${RES}" > "${INPUT}"

# Count unique individuals from input
N_TOTAL=$(awk '{print $1"\n"$2}' "${INPUT}" | sort -u | wc -l)

echo "[INFO] NAToRA input: ${INPUT}"
echo "[INFO] Total individuals: ${N_TOTAL}"
echo "[INFO] Total pairs: $(wc -l < "${INPUT}")"

# Header for summary
printf "group\tcutoff_label\ttheta_cutoff\tremoved_n\tretained_n\tremoved_pct\toutput_prefix\n" > "${SUMMARY}"

run_cutoff () {
  local GROUP="$1"
  local LABEL="$2"
  local CUTOFF="$3"
  local PREFIX="${OUTDIR}/natora_${GROUP}_${LABEL}_c${CUTOFF}"

  echo
  echo "[INFO] Running ${GROUP} ${LABEL} at theta cutoff ${CUTOFF}"

  # Force heuristic mode for consistency across runs
  python "${NATORA}" \
    -i "${INPUT}" \
    -o "${PREFIX}" \
    -c "${CUTOFF}" \
    -v 0.0 \
    -V 0.5 \
    -e 1

  local REMOVED=0
  if [[ -s "${PREFIX}_toRemove.txt" ]]; then
    REMOVED=$(sort -u "${PREFIX}_toRemove.txt" | wc -l)
  fi

  local RETAINED=$((N_TOTAL - REMOVED))
  local REMOVED_PCT
  REMOVED_PCT=$(awk -v r="${REMOVED}" -v n="${N_TOTAL}" 'BEGIN{printf "%.2f", (100*r)/n}')

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${GROUP}" "${LABEL}" "${CUTOFF}" "${REMOVED}" "${RETAINED}" "${REMOVED_PCT}" "${PREFIX}" >> "${SUMMARY}"
}

# Duplicate / MZ
run_cutoff "duplicate_mz" "low"  "0.3540"
run_cutoff "duplicate_mz" "high" "0.5000"

# First degree
run_cutoff "first_degree" "low"  "0.1770"
run_cutoff "first_degree" "high" "0.3540"

# Second degree
run_cutoff "second_degree" "low"  "0.0884"
run_cutoff "second_degree" "high" "0.1770"

# Third degree
run_cutoff "third_degree" "low"  "0.0442"
run_cutoff "third_degree" "high" "0.0884"

echo
echo "[DONE] Summary written to:"
echo "  ${SUMMARY}"

echo
echo "[INFO] Preview:"
column -t -s $'\t' "${SUMMARY}"
