#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=16G
#SBATCH -t 0-01:00:00
#SBATCH -J delly_prep
#SBATCH -o logs/01_prep.%j.out
#SBATCH -e logs/01_prep.%j.err
set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_delly_config.sh"

dv_init_dirs
dv_log "=== STEP 1: Prepare exclude BED + sample lists ==="

# ─────────────────────────────────────────────────────────────────────────────
# 1A. Build exclude BED from callable_bp_per_bin.tsv
#     Strategy: 50-kb bins where callable_bp < threshold across ALL 226 samples
#     = regions where nobody has reads = telomeric/centromeric/assembly-gap junk
#     Then merge adjacent uncallable bins into big blocks (>= 50 kb)
# ─────────────────────────────────────────────────────────────────────────────
dv_log "Building exclude BED from callable_bp_per_bin.tsv..."
dv_check_file "$PA_CALLABLE_BP" "callable_bp_per_bin.tsv"
dv_check_file "$REF_FAI" "Reference FASTA index"

python3 << 'PYEOF'
import os, sys

callable_bp_file = os.environ['PA_CALLABLE_BP']
ref_fai = os.environ['REF_FAI']
out_bed = os.environ['EXCL_BED']
min_callable = int(os.environ['EXCL_MIN_CALLABLE_BP'])
min_block = int(os.environ['EXCL_MIN_BLOCK_BP'])

# Read chromosome lengths
chr_lengths = {}
with open(ref_fai) as f:
    for line in f:
        parts = line.strip().split('\t')
        chr_lengths[parts[0]] = int(parts[1])

# Read callable_bp_per_bin.tsv: chrom, agg_start, agg_end, callable_bp
# Find bins where callable_bp < threshold → uncallable
uncallable_bins = []  # (chrom, start, end)
with open(callable_bp_file) as f:
    header = f.readline()
    for line in f:
        parts = line.strip().split('\t')
        chrom = parts[0]
        start = int(parts[1])
        end = int(parts[2])
        cbp = int(parts[3])
        if cbp < min_callable:
            uncallable_bins.append((chrom, start, end))

print(f"Found {len(uncallable_bins)} uncallable 50-kb bins (callable_bp < {min_callable})")

# Merge adjacent bins per chromosome
merged = []
if uncallable_bins:
    curr_chrom, curr_start, curr_end = uncallable_bins[0]
    for chrom, start, end in uncallable_bins[1:]:
        if chrom == curr_chrom and start <= curr_end:
            # Adjacent or overlapping → extend
            curr_end = max(curr_end, end)
        else:
            merged.append((curr_chrom, curr_start, curr_end))
            curr_chrom, curr_start, curr_end = chrom, start, end
    merged.append((curr_chrom, curr_start, curr_end))

print(f"After merging adjacent: {len(merged)} blocks")

# Filter to blocks >= min_block_bp
big_blocks = [(c, s, e) for c, s, e in merged if (e - s) >= min_block]
print(f"After size filter (>= {min_block} bp): {len(big_blocks)} blocks")

# Also add chromosome tips: first bin and last bin if uncallable
# (catch telomeric regions even if they're not in callable_bp)
# Strategy: for each chromosome, if the first 50 kb or last 50 kb have
# callable_bp < threshold, extend the exclude block to the chromosome edge

# Build quick lookup of callable_bp per (chrom, start)
cbp_lookup = {}
with open(callable_bp_file) as f:
    header = f.readline()
    for line in f:
        parts = line.strip().split('\t')
        cbp_lookup[(parts[0], int(parts[1]))] = int(parts[3])

tip_blocks = []
for chrom, length in chr_lengths.items():
    # Check first 50-kb bin
    first_cbp = cbp_lookup.get((chrom, 0), 0)
    if first_cbp < min_callable:
        # Extend: find how far the uncallable region goes
        pos = 0
        while pos < length:
            cbp = cbp_lookup.get((chrom, pos), -1)
            if cbp == -1 or cbp >= min_callable:
                break
            pos += 50000
        if pos >= min_block:
            tip_blocks.append((chrom, 0, pos))

    # Check last 50-kb bin
    last_bin_start = (length // 50000) * 50000
    if last_bin_start >= length:
        last_bin_start -= 50000
    last_cbp = cbp_lookup.get((chrom, last_bin_start), 0)
    if last_cbp < min_callable:
        # Walk backwards
        pos = last_bin_start
        while pos >= 0:
            cbp = cbp_lookup.get((chrom, pos), -1)
            if cbp == -1 or cbp >= min_callable:
                pos += 50000
                break
            pos -= 50000
        if pos < 0:
            pos = 0
        if (length - pos) >= min_block:
            tip_blocks.append((chrom, pos, length))

print(f"Tip blocks (telomeric uncallable): {len(tip_blocks)}")

# Combine big_blocks + tip_blocks, sort, merge
all_blocks = big_blocks + tip_blocks
all_blocks.sort(key=lambda x: (x[0], x[1]))

# Merge overlapping
final = []
for chrom, start, end in all_blocks:
    if final and final[-1][0] == chrom and start <= final[-1][2]:
        final[-1] = (chrom, final[-1][1], max(final[-1][2], end))
    else:
        final.append((chrom, start, end))

# Write BED
total_bp = 0
with open(out_bed, 'w') as f:
    for chrom, start, end in final:
        f.write(f"{chrom}\t{start}\t{end}\n")
        total_bp += (end - start)

genome_bp = sum(chr_lengths.values())
pct = 100.0 * total_bp / genome_bp if genome_bp > 0 else 0
print(f"\nExclude BED: {len(final)} regions, {total_bp:,} bp ({pct:.2f}% of {genome_bp:,} bp genome)")
print(f"Written to: {out_bed}")
PYEOF

N_EXCL=$(wc -l < "${EXCL_BED}")
dv_log "Exclude BED: ${N_EXCL} regions"
dv_log "Preview:"
head -20 "${EXCL_BED}" | while read line; do dv_log "  ${line}"; done

# ─────────────────────────────────────────────────────────────────────────────
# 1B. Build sample lists
# ─────────────────────────────────────────────────────────────────────────────
dv_log "Building sample lists..."

ls "${BAMDIR}"/*"${BAM_SUFFIX}" 2>/dev/null \
  | xargs -I{} basename {} "${BAM_SUFFIX}" \
  | sort \
  > "${SAMPLES_ALL}"

N_ALL=$(wc -l < "${SAMPLES_ALL}")
dv_log "All samples: ${N_ALL} -> ${SAMPLES_ALL}"

if [[ ${N_ALL} -ne 226 ]]; then
  dv_log "WARNING: Expected 226, found ${N_ALL}. Check BAMDIR='${BAMDIR}' BAM_SUFFIX='${BAM_SUFFIX}'"
fi

# Unrelated 81
dv_check_file "${NATORA_KEEP}" "NAToRA keep list"
cp "${NATORA_KEEP}" "${SAMPLES_UNRELATED}"
N_UNREL=$(wc -l < "${SAMPLES_UNRELATED}")
dv_log "Unrelated samples: ${N_UNREL} -> ${SAMPLES_UNRELATED}"

# ─────────────────────────────────────────────────────────────────────────────
# 1C. Validate BAMs + indexes
# ─────────────────────────────────────────────────────────────────────────────
dv_log "Validating BAMs..."
MISSING=0; NO_INDEX=0
while IFS= read -r sid; do
  bam="${BAMDIR}/${sid}${BAM_SUFFIX}"
  if [[ ! -f "$bam" ]]; then
    dv_err "Missing BAM: ${bam}"
    ((MISSING++)) || true
    continue
  fi
  if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
    dv_err "Missing index: ${bam}"
    ((NO_INDEX++)) || true
  fi
done < "${SAMPLES_ALL}"

if [[ ${MISSING} -gt 0 ]]; then
  dv_die "${MISSING} BAMs missing. Fix BAMDIR/BAM_SUFFIX in config."
fi
if [[ ${NO_INDEX} -gt 0 ]]; then
  dv_log "WARNING: ${NO_INDEX} BAMs lack index. Will attempt samtools index in discovery step."
fi

# ─────────────────────────────────────────────────────────────────────────────
# 1D. Sort annotation BEDs for bedtools
# ─────────────────────────────────────────────────────────────────────────────
dv_log "Preparing annotation BEDs..."
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
# 1E. Verify DELLY binary
# ─────────────────────────────────────────────────────────────────────────────
dv_log "Checking DELLY binary..."
if [[ ! -x "${DELLY_BIN}" ]]; then
  dv_die "DELLY binary not found or not executable: ${DELLY_BIN}"
fi
DELLY_VER=$("${DELLY_BIN}" 2>&1 | grep -oP 'Version: \K[0-9.]+' || echo "unknown")
dv_log "  DELLY version: ${DELLY_VER}"
dv_log "  Path: ${DELLY_BIN}"

dv_log "=== STEP 1 COMPLETE ==="
dv_log "Exclude BED: ${EXCL_BED}"
dv_log "Samples all: ${SAMPLES_ALL} (${N_ALL})"
dv_log "Samples unrel: ${SAMPLES_UNRELATED} (${N_UNREL})"
dv_log "Next: submit 02_delly_discovery.sh"
