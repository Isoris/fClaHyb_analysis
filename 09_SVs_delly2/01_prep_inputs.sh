#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 127
#SBATCH --mem=237G
#SBATCH -t 1-00:00:00
#SBATCH -J delly_prep
#SBATCH -o logs/01_prep.%j.out
#SBATCH -e logs/01_prep.%j.err
# =============================================================================
# 01_prep_inputs.sh — Prepare all inputs for DELLY DEL pipeline
# =============================================================================
# This script does:
#   1. Loads config (SLURM-safe)
#   2. Reads BAM manifest TSV (manifest-driven, no glob)
#   3. Marks duplicates on raw BAMs (samtools markdup)
#   4. Builds exclude BED (callable-based + unconditional chr-end masks)
#   5. Builds sample lists (all 226, unrelated 81)
#   6. Validates all inputs
#   7. Sorts annotation BEDs
# =============================================================================
set -euo pipefail
source ~/.bashrc
mamba activate assembly

# ── Config loading (SLURM-safe) ────────────────────────────────────────────
# IMPORTANT: Do NOT use BASH_SOURCE dirname under SLURM — it resolves to the
# spool directory. Use SLURM_SUBMIT_DIR instead.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/00_delly_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }

# set -a / set +a exports all variables so embedded Python can see them
# via os.environ[...]
set -a
source "${CONFIG}"
set +a

dv_init_dirs
dv_log "=== STEP 1: Prepare inputs for DELLY DEL pipeline ==="

# ─────────────────────────────────────────────────────────────────────────────
# 1A. Read BAM manifest
# ─────────────────────────────────────────────────────────────────────────────
# The manifest TSV has 3 columns:
#   col1 = sample_id
#   col2 = unfiltered BAM path (minimap2 alignment) → USE THIS
#   col3 = P99/TLEN/MAPQ30 filtered BAM → IGNORE for DELLY
#
# DELLY needs the unfiltered BAM (column 2) because filtered BAMs have had
# discordant/split reads removed, which are exactly what DELLY uses for
# SV discovery.
# ─────────────────────────────────────────────────────────────────────────────
dv_log "--- 1A: Reading BAM manifest ---"
dv_read_manifest

dv_log "  First 5 samples:"
for i in 0 1 2 3 4; do
  [[ $i -lt ${MANIFEST_N_SAMPLES} ]] || break
  dv_log "    ${MANIFEST_SIDS[$i]}  →  $(basename "${MANIFEST_RAW_BAMS[$i]}")"
done

# Validate that all raw BAMs exist
dv_log "  Validating raw BAM paths..."
MISSING_BAMS=0
for i in $(seq 0 $((MANIFEST_N_SAMPLES - 1))); do
  if [[ ! -f "${MANIFEST_RAW_BAMS[$i]}" ]]; then
    dv_err "Missing BAM for ${MANIFEST_SIDS[$i]}: ${MANIFEST_RAW_BAMS[$i]}"
    ((MISSING_BAMS++)) || true
  fi
done
[[ ${MISSING_BAMS} -eq 0 ]] || dv_die "${MISSING_BAMS} BAMs missing from manifest. Aborting."
dv_log "  All ${MANIFEST_N_SAMPLES} raw BAMs found."

# ─────────────────────────────────────────────────────────────────────────────
# 1B. Mark duplicates
# ─────────────────────────────────────────────────────────────────────────────
# DELLY REQUIRES duplicate-marked BAMs. The raw minimap2 BAMs from the
# manifest (column 2) do not have duplicates marked.
#
# IMPORTANT:
#   We do NOT run markdup directly on the raw coordinate-sorted BAM because
#   that produced header-only / zero-read BAMs in testing.
#
#   Safe pipeline for paired-end BAMs:
#     raw BAM -> name-sort -> fixmate -m -> coordinate-sort -> markdup -> index
#
#   We only MARK duplicates. We do NOT remove duplicates.
#   We do NOT run clipOverlap here.
# ─────────────────────────────────────────────────────────────────────────────
dv_log "--- 1B: Marking duplicates (name-sort -> fixmate -> coord-sort -> markdup) ---"
dv_log "  Output directory: ${DIR_MARKDUP}"
dv_log "  Parallel jobs: ${DELLY_PARALLEL}"

# Optional command checks
dv_check_cmd samtools
dv_check_cmd parallel

# Conservative per-job threading inside GNU parallel
MARKDUP_THREADS=2

run_markdup() {
  local sid="$1"
  local raw_bam="$2"

  local tmpdir="${DIR_MARKDUP}/tmp_${sid}"
  local ns_bam="${tmpdir}/${sid}.namesort.bam"
  local fx_bam="${tmpdir}/${sid}.fixmate.bam"
  local cs_bam="${tmpdir}/${sid}.coordsort.bam"
  local out_bam="${DIR_MARKDUP}/${sid}.markdup.bam"
  local logf="${DIR_LOGS}/markdup_${sid}.log"

  mkdir -p "${tmpdir}"

  # Skip only if BAM + BAI exist AND BAM has reads
  if [[ -f "${out_bam}" && -f "${out_bam}.bai" ]]; then
    local existing_reads
    existing_reads=$(samtools view -c "${out_bam}" 2>>"${logf}" || echo 0)
    if [[ "${existing_reads}" -gt 0 ]]; then
      echo "[$(date '+%F %T')] ${sid}: markdup BAM exists with ${existing_reads} reads, skipping" >> "${logf}"
      return 0
    else
      echo "[$(date '+%F %T')] ${sid}: existing markdup BAM has 0 reads; rebuilding" >> "${logf}"
      rm -f "${out_bam}" "${out_bam}.bai"
    fi
  fi

  echo "[$(date '+%F %T')] ${sid}: START" > "${logf}"
  echo "[$(date '+%F %T')] ${sid}: raw BAM = ${raw_bam}" >> "${logf}"

  # Sanity: input exists and is non-empty
  if [[ ! -s "${raw_bam}" ]]; then
    echo "[$(date '+%F %T')] ${sid}: ERROR input BAM missing or empty: ${raw_bam}" >> "${logf}"
    return 1
  fi

  # Sanity: input BAM has reads
  local raw_reads
  raw_reads=$(samtools view -c "${raw_bam}" 2>>"${logf}" || echo 0)
  echo "[$(date '+%F %T')] ${sid}: raw read count = ${raw_reads}" >> "${logf}"
  if [[ "${raw_reads}" -eq 0 ]]; then
    echo "[$(date '+%F %T')] ${sid}: ERROR raw BAM has 0 reads" >> "${logf}"
    return 1
  fi

  echo "[$(date '+%F %T')] ${sid}: STEP name-sort" >> "${logf}"
  samtools sort \
    -n \
    -@ "${MARKDUP_THREADS}" \
    -T "${tmpdir}/${sid}.nsort" \
    -o "${ns_bam}" \
    "${raw_bam}" \
    2>>"${logf}"

  echo "[$(date '+%F %T')] ${sid}: STEP fixmate -m" >> "${logf}"
  samtools fixmate \
    -m \
    -@ "${MARKDUP_THREADS}" \
    "${ns_bam}" \
    "${fx_bam}" \
    2>>"${logf}"

  echo "[$(date '+%F %T')] ${sid}: STEP coordinate-sort" >> "${logf}"
  samtools sort \
    -@ "${MARKDUP_THREADS}" \
    -T "${tmpdir}/${sid}.csort" \
    -o "${cs_bam}" \
    "${fx_bam}" \
    2>>"${logf}"

  echo "[$(date '+%F %T')] ${sid}: STEP markdup (mark only)" >> "${logf}"
  samtools markdup \
    -@ "${MARKDUP_THREADS}" \
    -s \
    "${cs_bam}" \
    "${out_bam}" \
    2>>"${logf}"

  echo "[$(date '+%F %T')] ${sid}: STEP index" >> "${logf}"
  samtools index \
    -@ "${MARKDUP_THREADS}" \
    "${out_bam}" \
    2>>"${logf}"

  # Final sanity checks
  local out_reads
  out_reads=$(samtools view -c "${out_bam}" 2>>"${logf}" || echo 0)
  echo "[$(date '+%F %T')] ${sid}: output read count = ${out_reads}" >> "${logf}"

  if [[ "${out_reads}" -eq 0 ]]; then
    echo "[$(date '+%F %T')] ${sid}: ERROR markdup BAM has 0 reads" >> "${logf}"
    return 1
  fi

  samtools quickcheck -v "${out_bam}" >> "${logf}" 2>&1 || {
    echo "[$(date '+%F %T')] ${sid}: ERROR samtools quickcheck failed" >> "${logf}"
    return 1
  }

  echo "[$(date '+%F %T')] ${sid}: DONE" >> "${logf}"

  # Keep intermediates for now for safety/debugging
  # Do NOT delete tmp files yet until pipeline is confirmed working.
  return 0
}

export -f run_markdup
export DIR_MARKDUP DIR_LOGS MARKDUP_THREADS

# Build a two-column TSV for parallel: sample_id <tab> raw_bam_path
MARKDUP_JOBS="${DIR_LOGS}/markdup_jobs.tsv"
for i in $(seq 0 $((MANIFEST_N_SAMPLES - 1))); do
  echo -e "${MANIFEST_SIDS[$i]}\t${MANIFEST_RAW_BAMS[$i]}"
done > "${MARKDUP_JOBS}"

parallel -j "${DELLY_PARALLEL}" --colsep '\t' \
  --joblog "${DIR_LOGS}/markdup_parallel.log" \
  run_markdup {1} {2} \
  :::: "${MARKDUP_JOBS}"

# Verify all markdup BAMs exist and are non-empty in terms of read count
dv_log "  Verifying markdup outputs..."
MARKDUP_FAIL=0
for i in $(seq 0 $((MANIFEST_N_SAMPLES - 1))); do
  sid="${MANIFEST_SIDS[$i]}"
  mbam="${DIR_MARKDUP}/${sid}.markdup.bam"

  if [[ ! -f "${mbam}" ]]; then
    dv_err "Missing markdup BAM: ${mbam}"
    ((MARKDUP_FAIL++)) || true
    continue
  fi

  if [[ ! -f "${mbam}.bai" ]]; then
    dv_err "Missing index for: ${mbam}"
    ((MARKDUP_FAIL++)) || true
    continue
  fi

  nreads=$(samtools view -c "${mbam}" 2>/dev/null || echo 0)
  if [[ "${nreads}" -eq 0 ]]; then
    dv_err "Markdup BAM has 0 reads: ${mbam}"
    ((MARKDUP_FAIL++)) || true
    continue
  fi
done

[[ ${MARKDUP_FAIL} -eq 0 ]] || dv_die "${MARKDUP_FAIL} markdup BAMs missing / unindexed / empty."
dv_log "  All ${MANIFEST_N_SAMPLES} markdup BAMs ready and non-empty."

# ─────────────────────────────────────────────────────────────────────────────
# 1C. Build sample lists
# ─────────────────────────────────────────────────────────────────────────────
dv_log "--- 1C: Building sample lists ---"

# All samples: from manifest column 1 (NOT from ls glob)
printf '%s\n' "${MANIFEST_SIDS[@]}" | sort > "${SAMPLES_ALL}"
N_ALL=$(wc -l < "${SAMPLES_ALL}")
dv_log "  All samples: ${N_ALL} → ${SAMPLES_ALL}"

if [[ ${N_ALL} -ne 226 ]]; then
  dv_log "  WARNING: Expected 226 samples, found ${N_ALL}."
fi

# Unrelated 81: from NAToRA keep list
dv_check_file "${NATORA_KEEP}" "NAToRA keep list"
cp "${NATORA_KEEP}" "${SAMPLES_UNRELATED}"
N_UNREL=$(wc -l < "${SAMPLES_UNRELATED}")
dv_log "  Unrelated: ${N_UNREL} → ${SAMPLES_UNRELATED}"

# ─────────────────────────────────────────────────────────────────────────────
# 1D. Build exclude BED
# ─────────────────────────────────────────────────────────────────────────────
# Two sources, merged:
#   1. Callable-based: 50-kb bins with callable_bp < threshold → uncallable
#   2. Chromosome-end: unconditionally mask first/last CHR_END_MASK_BP of
#      every chromosome (catches telomeric/centromeric junk)
# ─────────────────────────────────────────────────────────────────────────────
dv_log "--- 1D: Building exclude BED ---"
dv_check_file "$PA_CALLABLE_BP" "callable_bp_per_bin.tsv"
dv_check_file "$REF_FAI" "Reference FASTA index"

python3 << 'PYEOF'
import os

callable_bp_file = os.environ['PA_CALLABLE_BP']
ref_fai          = os.environ['REF_FAI']
out_bed          = os.environ['EXCL_BED']
min_callable     = int(os.environ['EXCL_MIN_CALLABLE_BP'])
min_block        = int(os.environ['EXCL_MIN_BLOCK_BP'])
chr_end_mask     = int(os.environ['CHR_END_MASK_BP'])

# --- Read chromosome lengths from .fai ---
chr_lengths = {}
chr_order = []
with open(ref_fai) as f:
    for line in f:
        parts = line.strip().split('\t')
        chr_lengths[parts[0]] = int(parts[1])
        chr_order.append(parts[0])

# --- Source 1: callable-based uncallable bins ---
uncallable = []
with open(callable_bp_file) as f:
    header = f.readline()
    for line in f:
        parts = line.strip().split('\t')
        chrom, start, end, cbp = parts[0], int(parts[1]), int(parts[2]), int(parts[3])
        if cbp < min_callable:
            uncallable.append((chrom, start, end))
print(f"Callable-based uncallable bins: {len(uncallable)}")

# Merge adjacent bins
merged_callable = []
if uncallable:
    cc, cs, ce = uncallable[0]
    for chrom, start, end in uncallable[1:]:
        if chrom == cc and start <= ce:
            ce = max(ce, end)
        else:
            merged_callable.append((cc, cs, ce))
            cc, cs, ce = chrom, start, end
    merged_callable.append((cc, cs, ce))

# Filter to big blocks
big_callable = [(c, s, e) for c, s, e in merged_callable if (e - s) >= min_block]
print(f"Callable-based blocks >= {min_block} bp: {len(big_callable)}")

# --- Source 2: unconditional chromosome-end masking ---
# Always exclude first and last chr_end_mask bp of every chromosome.
# This is NOT conditional on callable_bp.
chr_end_blocks = []
for chrom in chr_order:
    length = chr_lengths[chrom]
    # First chr_end_mask bp
    end1 = min(chr_end_mask, length)
    chr_end_blocks.append((chrom, 0, end1))
    # Last chr_end_mask bp
    start2 = max(0, length - chr_end_mask)
    if start2 < length:
        chr_end_blocks.append((chrom, start2, length))
print(f"Chromosome-end blocks: {len(chr_end_blocks)} ({chr_end_mask} bp each end)")

# --- Combine and merge ---
all_blocks = big_callable + chr_end_blocks
all_blocks.sort(key=lambda x: (chr_order.index(x[0]) if x[0] in chr_order else 999, x[1]))

final = []
for chrom, start, end in all_blocks:
    if final and final[-1][0] == chrom and start <= final[-1][2]:
        final[-1] = (chrom, final[-1][1], max(final[-1][2], end))
    else:
        final.append((chrom, start, end))

# --- Write BED ---
total_bp = 0
with open(out_bed, 'w') as f:
    for chrom, start, end in final:
        f.write(f"{chrom}\t{start}\t{end}\n")
        total_bp += (end - start)

genome_bp = sum(chr_lengths.values())
pct = 100.0 * total_bp / genome_bp if genome_bp > 0 else 0
n_chrs = len(chr_lengths)
print(f"\nFinal exclude BED: {len(final)} regions")
print(f"  Total excluded: {total_bp:,} bp ({pct:.2f}% of {genome_bp:,} bp genome)")
print(f"  Chromosomes: {n_chrs}")
print(f"  Chromosome-end mask: {chr_end_mask} bp (unconditional)")
print(f"  Callable threshold: callable_bp < {min_callable} in 50-kb bins")
print(f"Written to: {out_bed}")
PYEOF

N_EXCL=$(wc -l < "${EXCL_BED}")
dv_log "  Exclude BED: ${N_EXCL} regions"
dv_log "  Preview (first 15 lines):"
head -15 "${EXCL_BED}" | while IFS= read -r line; do dv_log "    ${line}"; done
dv_log "  Preview (last 10 lines):"
tail -10 "${EXCL_BED}" | while IFS= read -r line; do dv_log "    ${line}"; done

# ─────────────────────────────────────────────────────────────────────────────
# 1E. Sort annotation BEDs for bedtools intersect
# ─────────────────────────────────────────────────────────────────────────────
dv_log "--- 1E: Preparing annotation BEDs ---"
ANNOT_OUT="${DIR_ANNOT}/beds"
mkdir -p "${ANNOT_OUT}"

for bed_label in GENE EXON CDS; do
  src="${ANNOT_DIR}/features/${bed_label}.bed"
  if [[ -f "$src" ]]; then
    bedtools sort -i "$src" > "${ANNOT_OUT}/${bed_label}.sorted.bed"
    dv_log "  ${bed_label}: $(wc -l < "${ANNOT_OUT}/${bed_label}.sorted.bed") features"
  else
    dv_log "  WARNING: ${bed_label}.bed not found at ${src}"
  fi
done

if [[ -n "${REPEAT_BED}" && -f "${REPEAT_BED}" ]]; then
  bedtools sort -i "${REPEAT_BED}" > "${ANNOT_OUT}/REPEATS.sorted.bed"
  dv_log "  REPEATS: $(wc -l < "${ANNOT_OUT}/REPEATS.sorted.bed") intervals"
fi

# ─────────────────────────────────────────────────────────────────────────────
# 1F. Verify DELLY binary
# ─────────────────────────────────────────────────────────────────────────────
dv_log "--- 1F: Checking DELLY binary ---"
[[ -x "${DELLY_BIN}" ]] || dv_die "DELLY binary not executable: ${DELLY_BIN}"
DELLY_VER=$("${DELLY_BIN}" 2>&1 | grep -oP 'Version: \K[0-9.]+' || echo "unknown")
dv_log "  DELLY ${DELLY_VER} at ${DELLY_BIN}"

# ─────────────────────────────────────────────────────────────────────────────
# Done
# ─────────────────────────────────────────────────────────────────────────────
dv_log "=== STEP 1 COMPLETE ==="
dv_log "  Markdup BAMs:   ${DIR_MARKDUP}/ (${MANIFEST_N_SAMPLES} samples)"
dv_log "  Exclude BED:    ${EXCL_BED} (${N_EXCL} regions)"
dv_log "  Samples all:    ${SAMPLES_ALL} (${N_ALL})"
dv_log "  Samples unrel:  ${SAMPLES_UNRELATED} (${N_UNREL})"
dv_log "Next: submit 02_delly_discovery.sh"
