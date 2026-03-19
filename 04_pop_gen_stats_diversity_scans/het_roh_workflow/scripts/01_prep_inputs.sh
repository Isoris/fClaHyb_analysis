#!/usr/bin/env bash
# =============================================================================
# 01_prep_inputs.sh — Validate inputs and prepare BAM list for het/ROH workflow
# =============================================================================
set -euo pipefail
SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPTDIR}/00_config.sh"
hr_init_dirs

hr_log "=== Step 01: Prepare and validate inputs ==="

# ── 1. Extract BAM list from manifest (column 3 = filtered BAMs) ───────────
hr_check_file "${SAMPLE_MANIFEST}" "sample manifest"

hr_log "Extracting filtered BAM paths from manifest column 3..."
awk -F'\t' 'NR>1 && $3!="" { print $3 }' "${SAMPLE_MANIFEST}" > "${BAMLIST_QCPASS}"
awk -F'\t' 'NR>1 && $3!="" { print $1 }' "${SAMPLE_MANIFEST}" > "${SAMPLE_LIST}"

N_SAMPLES=$(wc -l < "${BAMLIST_QCPASS}")
hr_log "Found ${N_SAMPLES} QC-passing BAMs"

# ── 2. Validate reference ──────────────────────────────────────────────────
hr_check_file "${REF}" "reference FASTA"
hr_check_file "${REF_FAI}" "reference FASTA index"
hr_log "Reference: ${REF}"

# ── 3. Validate callable BED (intervals for span calculation) ─────────────
hr_check_file "${CALLABLE_BED}" "callable BED (interval format)"
CALLABLE_BP=$(awk '{s+=($3-$2)} END{print s}' "${CALLABLE_BED}")
hr_log "Callable genome span: ${CALLABLE_BP} bp"
echo -e "callable_nonrepeat_bp\t${CALLABLE_BP}" > "${DIR_INPUTS}/callable_genome_size.tsv"

# ── 4. Validate ANGSD callable sites (indexed, for -sites) ───────────────
hr_check_file "${CALLABLE_SITES}" "ANGSD callable sites"
hr_check_file "${CALLABLE_SITES}.idx" "ANGSD callable sites index"
hr_check_file "${CALLABLE_SITES}.bin" "ANGSD callable sites bin"
hr_log "ANGSD callable sites: ${CALLABLE_SITES} (+.idx +.bin OK)"

# ── 5. Validate BAM files exist ───────────────────────────────────────────
hr_log "Checking BAM files exist..."
MISSING=0
while read -r BAM; do
  if [[ ! -f "${BAM}" ]]; then
    hr_err "BAM not found: ${BAM}"
    MISSING=$((MISSING + 1))
  fi
done < "${BAMLIST_QCPASS}"

if [[ "${MISSING}" -gt 0 ]]; then
  hr_die "${MISSING} BAM files missing. Fix manifest and re-run."
fi
hr_log "All ${N_SAMPLES} BAM files found."

# ── 6. Generate chromosome sizes from .fai ────────────────────────────────
awk -F'\t' '{print $1"\t"$2}' "${REF_FAI}" > "${DIR_INPUTS}/chrom_sizes.tsv"
hr_log "Wrote chrom_sizes.tsv ($(wc -l < "${DIR_INPUTS}/chrom_sizes.tsv") contigs)"

# ── 7. Create callable_regions.rf from BED (for ANGSD -rf) ───────────────
# ANGSD -rf format: chrom:start-end (1-based, inclusive)
awk '{printf "%s:%d-%d\n", $1, $2+1, $3}' "${CALLABLE_BED}" > "${DIR_INPUTS}/callable_regions.rf"
hr_log "Wrote callable_regions.rf for ANGSD -rf"

# ── 8. Validate ngsF-HMM inputs if present ────────────────────────────────
if [[ -n "${BEAGLE_GL}" && -f "${BEAGLE_GL}" ]]; then
  hr_log "BEAGLE GL: ${BEAGLE_GL} (exists)"
else
  hr_log "WARNING: BEAGLE GL not found at ${BEAGLE_GL} — ngsF-HMM step will need this"
fi

if [[ -n "${POS_FILE}" && -f "${POS_FILE}" ]]; then
  N_POS=$(wc -l < "${POS_FILE}")
  hr_log "Position file: ${POS_FILE} (${N_POS} sites)"
else
  hr_log "WARNING: Position file not found — ngsF-HMM step will need this"
fi

if [[ -n "${SAMPLES_IND}" && -f "${SAMPLES_IND}" ]]; then
  N_IND=$(wc -l < "${SAMPLES_IND}")
  hr_log "Samples.ind: ${SAMPLES_IND} (${N_IND} individuals)"
else
  hr_log "WARNING: samples.ind not found — ngsF-HMM step will need this"
fi

# ── 9. Validate ngsF-HMM binary ──────────────────────────────────────────
if [[ -n "${NGSFHMM_BIN}" && -x "${NGSFHMM_BIN}" ]]; then
  hr_log "ngsF-HMM binary: ${NGSFHMM_BIN} (OK)"
elif [[ -n "${NGSFHMM_BIN}" && -f "${NGSFHMM_BIN}" ]]; then
  hr_log "WARNING: ngsF-HMM binary exists but is not executable: ${NGSFHMM_BIN}"
else
  hr_log "WARNING: ngsF-HMM binary not found at ${NGSFHMM_BIN}"
fi

# ── 10. Summary ────────────────────────────────────────────────────────────
cat > "${DIR_INPUTS}/inputs_summary.txt" <<EOF
=== Het/ROH Workflow Input Summary ===
Date: $(date)
Reference: ${REF}
Callable BED (intervals): ${CALLABLE_BED}
Callable ANGSD sites:     ${CALLABLE_SITES}
Callable span: ${CALLABLE_BP} bp
N QC-pass samples: ${N_SAMPLES}
BAM list: ${BAMLIST_QCPASS}
Sample list: ${SAMPLE_LIST}
BEAGLE GL: ${BEAGLE_GL:-not set}
Position file: ${POS_FILE:-not set}
Samples.ind: ${SAMPLES_IND:-not set}
ngsF-HMM binary: ${NGSFHMM_BIN:-not set}
Repeat BED: ${REPEAT_BED:-not set}
ngsRelate: ${NGSRELATE_PAIRS:-not set}
Pruned 81: ${PRUNED81_SAMPLES:-not set}
Q matrix: ${Q_MATRIX:-not set}
Ancestry labels: ${ANCESTRY_LABELS:-not set}
Covtree order: ${COVTREE_ORDER:-not set}
EOF

hr_log "Input summary written to ${DIR_INPUTS}/inputs_summary.txt"
hr_log "=== Step 01 complete ==="
