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

# ── Config (SLURM-safe) ───────────────────────────────────────────────────
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/00_delly_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a
source "${CONFIG}"
set +a

dv_init_dirs
dv_log "=== STEP 2: Per-sample DELLY DEL discovery ==="

# Validate
dv_check_file "$REF" "Reference"
dv_check_file "$EXCL_BED" "Exclusion BED"
dv_check_file "$SAMPLES_ALL" "Sample list"
[[ -x "${DELLY_BIN}" ]] || dv_die "DELLY not executable: ${DELLY_BIN}"

N_SAMPLES=$(wc -l < "$SAMPLES_ALL")
dv_log "Samples: ${N_SAMPLES}"
dv_log "DELLY threads/call: ${DELLY_THREADS_PER_CALL}, parallel: ${DELLY_PARALLEL}"

# ─────────────────────────────────────────────────────────────────────────────
# Per-sample discovery
# Uses markdup BAMs from DIR_MARKDUP (created by 01_prep)
# Uses -t DEL to restrict to deletions at call time
# Uses -h threads for DELLY internal threading
# ─────────────────────────────────────────────────────────────────────────────
run_discovery() {
  local sid="$1"
  local bam="${DIR_MARKDUP}/${sid}.markdup.bam"
  local out_bcf="${DIR_DISC}/${sid}.disc.DEL.bcf"
  local logf="${DIR_LOGS}/disc_${sid}.log"

  if [[ -f "${out_bcf}" && -f "${out_bcf}.csi" ]]; then
    echo "[$(date '+%T')] ${sid}: already done, skipping"
    return 0
  fi

  if [[ ! -f "${bam}" ]]; then
    echo "[$(date '+%T')] ${sid}: ERROR — markdup BAM missing: ${bam}" >&2
    return 1
  fi

  echo "[$(date '+%T')] ${sid}: starting DEL discovery..."

  "${DELLY_BIN}" call \
    -t DEL \
    -g "${REF}" \
    -x "${EXCL_BED}" \
    -h "${DELLY_THREADS_PER_CALL}" \
    -o "${out_bcf}" \
    "${bam}" \
    2>>"${logf}"

  bcftools index -f "${out_bcf}" 2>>"${logf}"

  local n_del
  n_del=$(bcftools view -H "${out_bcf}" 2>/dev/null | wc -l)
  echo "[$(date '+%T')] ${sid}: ${n_del} DELs"
}

export -f run_discovery
export DIR_MARKDUP REF EXCL_BED DELLY_BIN DELLY_THREADS_PER_CALL DIR_DISC DIR_LOGS

parallel -j "${DELLY_PARALLEL}" --joblog "${DIR_LOGS}/discovery_parallel.log" \
  run_discovery {} \
  :::: "${SAMPLES_ALL}"

# Verify
dv_log "Verifying outputs..."
FAIL=0
while IFS= read -r sid; do
  [[ -f "${DIR_DISC}/${sid}.disc.DEL.bcf" ]] || { dv_err "Missing: ${sid}"; ((FAIL++)) || true; }
done < "$SAMPLES_ALL"
[[ ${FAIL} -eq 0 ]] || dv_die "${FAIL} samples failed. Check ${DIR_LOGS}/disc_*.log"

# Summary
{
  echo -e "sample\tn_DEL"
  for bcf in "${DIR_DISC}"/*.disc.DEL.bcf; do
    sid=$(basename "$bcf" .disc.DEL.bcf)
    n=$(bcftools view -H "$bcf" 2>/dev/null | wc -l)
    echo -e "${sid}\t${n}"
  done
} | tee "${DIR_LOGS}/discovery_DEL_counts.tsv"

dv_log "=== STEP 2 COMPLETE ==="
