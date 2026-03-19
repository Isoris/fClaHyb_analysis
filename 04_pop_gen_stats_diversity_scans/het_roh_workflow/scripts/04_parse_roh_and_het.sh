#!/usr/bin/env bash
# =============================================================================
# 04_parse_roh_and_het.sh — Wrapper to run the ROH parsing pipeline
# =============================================================================
set -euo pipefail
SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPTDIR}/00_config.sh"
hr_init_dirs

hr_log "=== Step 04: Parse ROH + compute het in/out ROH ==="

# Paths to utility scripts (should be in utils/ alongside scripts/)
CONVERT_IBD_PL="${SCRIPTDIR}/../utils/convert_ibd.pl"
SUMMARIZE_ROH_PY="${SCRIPTDIR}/../utils/summarize_ibd_roh.py"

hr_check_file "${CONVERT_IBD_PL}" "convert_ibd.pl"
hr_check_file "${SUMMARIZE_ROH_PY}" "summarize_ibd_roh.py"
hr_check_file "${DIR_NGSFHMM}/best/best.ibd" "best .ibd"

# Build command
CMD=(
  python3 "${SCRIPTDIR}/04_parse_roh_and_het.py"
  --ibd "${DIR_NGSFHMM}/best/best.ibd"
  --pos "${POS_FILE}"
  --ind "${SAMPLES_IND}"
  --callable-bed "${CALLABLE_BED}"
  --het-summary "${DIR_HET}/04_summary/genomewide_heterozygosity.tsv"
  --theta-dir "${DIR_HET}/03_theta"
  --out-dir "${DIR_ROH}"
  --tables-dir "${DIR_TABLES}"
  --convert-ibd-pl "${CONVERT_IBD_PL}"
  --summarize-roh-py "${SUMMARIZE_ROH_PY}"
  --win-size "${WIN}"
)

# Optional inputs
if [[ -n "${REPEAT_BED}" && -f "${REPEAT_BED}" ]]; then
  CMD+=(--repeats-bed "${REPEAT_BED}")
fi
if [[ -n "${Q_MATRIX}" && -f "${Q_MATRIX}" ]]; then
  CMD+=(--q-matrix "${Q_MATRIX}" --q-has-sample-col)
fi
if [[ -n "${COVTREE_ORDER}" && -f "${COVTREE_ORDER}" ]]; then
  CMD+=(--order-file "${COVTREE_ORDER}")
fi
if [[ -n "${NGSRELATE_PAIRS}" && -f "${NGSRELATE_PAIRS}" ]]; then
  CMD+=(--kinship "${NGSRELATE_PAIRS}")
fi

"${CMD[@]}" 2>&1 | tee "${DIR_LOGS}/04_parse_roh.log"

hr_log "=== Step 04 complete ==="
