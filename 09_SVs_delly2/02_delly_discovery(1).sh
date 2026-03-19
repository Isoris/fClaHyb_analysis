#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 127
#SBATCH --mem=237G
#SBATCH -t 2-00:00:00
#SBATCH -J delly_disc
#SBATCH -o logs/02_discovery.%j.out
#SBATCH -e logs/02_discovery.%j.err
set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_delly_config.sh"

dv_init_dirs
dv_log "=== STEP 2: Per-sample DELLY DEL discovery (226 samples) ==="

# Validate
dv_check_file "$REF" "Reference"
dv_check_file "$EXCL_BED" "Exclusion BED"
dv_check_file "$SAMPLES_ALL" "Sample list"
[[ -x "${DELLY_BIN}" ]] || dv_die "DELLY binary not executable: ${DELLY_BIN}"

N_SAMPLES=$(wc -l < "$SAMPLES_ALL")
dv_log "Samples: ${N_SAMPLES}"
dv_log "DELLY threads per call: ${DELLY_THREADS_PER_CALL}"
dv_log "Parallel calls: ${DELLY_PARALLEL}"
dv_log "Exclusion BED: ${EXCL_BED} ($(wc -l < "${EXCL_BED}") regions)"

# ─────────────────────────────────────────────────────────────────────────────
# Per-sample discovery function
# Uses -t DEL to restrict to deletions at discovery time (no post-hoc subset)
# Uses -h <threads> for DELLY internal threading
# ─────────────────────────────────────────────────────────────────────────────
run_discovery() {
  local sid="$1"
  local bam="${BAMDIR}/${sid}${BAM_SUFFIX}"
  local out_bcf="${DIR_DISC}/${sid}.disc.DEL.bcf"
  local logf="${DIR_LOGS}/disc_${sid}.log"

  # Skip if done
  if [[ -f "${out_bcf}" && -f "${out_bcf}.csi" ]]; then
    echo "[$(date '+%T')] ${sid}: already done, skipping"
    return 0
  fi

  echo "[$(date '+%T')] ${sid}: starting DEL discovery..."

  # Index BAM if needed
  if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
    samtools index "${bam}" 2>>"${logf}"
  fi

  # DELLY call: -t DEL restricts to deletions, -h threads
  "${DELLY_BIN}" call \
    -t DEL \
    -g "${REF}" \
    -x "${EXCL_BED}" \
    -h "${DELLY_THREADS_PER_CALL}" \
    -o "${out_bcf}" \
    "${bam}" \
    2>>"${logf}"

  bcftools index "${out_bcf}" 2>>"${logf}"

  local n_del
  n_del=$(bcftools view -H "${out_bcf}" 2>/dev/null | wc -l)
  echo "[$(date '+%T')] ${sid}: ${n_del} DELs discovered"
}

export -f run_discovery
export BAMDIR BAM_SUFFIX REF EXCL_BED DELLY_BIN DELLY_THREADS_PER_CALL DIR_DISC DIR_LOGS

# Run in parallel
parallel -j "${DELLY_PARALLEL}" --joblog "${DIR_LOGS}/discovery_parallel.log" \
  run_discovery {} \
  :::: "${SAMPLES_ALL}"

# ─────────────────────────────────────────────────────────────────────────────
# Verify outputs
# ─────────────────────────────────────────────────────────────────────────────
dv_log "Verifying discovery outputs..."
FAIL=0
while IFS= read -r sid; do
  if [[ ! -f "${DIR_DISC}/${sid}.disc.DEL.bcf" ]]; then
    dv_err "Missing: ${sid}"
    ((FAIL++)) || true
  fi
done < "$SAMPLES_ALL"

[[ ${FAIL} -eq 0 ]] || dv_die "${FAIL} samples failed. Check ${DIR_LOGS}/disc_*.log"

# Summary
dv_log "Discovery complete for ${N_SAMPLES} samples."
dv_log "Per-sample DEL counts:"
{
  echo -e "sample\tn_DEL"
  for bcf in "${DIR_DISC}"/*.disc.DEL.bcf; do
    sid=$(basename "$bcf" .disc.DEL.bcf)
    n=$(bcftools view -H "$bcf" 2>/dev/null | wc -l)
    echo -e "${sid}\t${n}"
  done
} | tee "${DIR_LOGS}/discovery_DEL_counts.tsv"

TOTAL_DEL=$(awk 'NR>1 {s+=$2} END {print s}' "${DIR_LOGS}/discovery_DEL_counts.tsv")
MEDIAN_DEL=$(awk 'NR>1 {print $2}' "${DIR_LOGS}/discovery_DEL_counts.tsv" | sort -n | awk '{a[NR]=$1} END {print a[int(NR/2)+1]}')
dv_log "Total DELs across samples: ${TOTAL_DEL}"
dv_log "Median DELs per sample: ${MEDIAN_DEL}"

dv_log "=== STEP 2 COMPLETE ==="
dv_log "Next: submit 03_delly_merge_genotype.sh"
