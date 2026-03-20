#!/usr/bin/env bash
set -euo pipefail

###
# STEP02_make_inversion_callable_angsd_mask.sh
#
# Usage:
#   bash STEP02_make_inversion_callable_angsd_mask.sh \
#     --fasta genome.fa \
#     --mask-script STEP01_mask_regions_from_fasta.py \
#     --prefix outprefix
###

usage() {
    cat <<EOF
Usage:
  $0 --fasta genome.fa --mask-script STEP01_mask_regions_from_fasta.py --prefix outprefix
EOF
}

FASTA=""
MASK_SCRIPT=""
PREFIX=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --fasta) FASTA="$2"; shift 2 ;;
        --mask-script) MASK_SCRIPT="$2"; shift 2 ;;
        --prefix) PREFIX="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
    esac
done

[[ -n "$FASTA" && -n "$MASK_SCRIPT" && -n "$PREFIX" ]] || { usage >&2; exit 1; }
[[ -f "$FASTA" ]] || { echo "[ERROR] Missing FASTA: $FASTA" >&2; exit 1; }
[[ -f "$MASK_SCRIPT" ]] || { echo "[ERROR] Missing script: $MASK_SCRIPT" >&2; exit 1; }

python3 "$MASK_SCRIPT" -i "$FASTA" -p "$PREFIX"

INPUT_BED="${PREFIX}.inversion_acgt_allcase.bed"
OUTPUT_ANGSD="${PREFIX}.inversion_acgt_allcase.1-based.angsd"

[[ -f "$INPUT_BED" ]] || { echo "[ERROR] Expected BED not found: $INPUT_BED" >&2; exit 1; }

awk 'BEGIN{OFS="\t"} {print $1, $2+1, $3}' "$INPUT_BED" > "$OUTPUT_ANGSD"

angsd sites index "$OUTPUT_ANGSD"

echo "[INFO] Done."
echo "[INFO] Callable BED : $INPUT_BED"
echo "[INFO] ANGSD sites  : $OUTPUT_ANGSD"
