#!/usr/bin/env bash
set -euo pipefail

###
# STEP01_make_callable_angsd_mask_from_fasta.sh
#
# Does this:
# 1. Runs mask_regions_from_fasta.py on the reference FASTA
# 2. Uses the generated normalACGT BED as callable regions
# 3. Converts BED start from 0-based to 1-based for ANGSD
# 4. Indexes the ANGSD sites file
#
# Usage:
#   bash STEP01_make_callable_angsd_mask_from_fasta.sh \
#     --fasta /path/to/fClaHyb_Gar_LG.fa \
#     --mask-script /project/lt200308-agbsci/01-catfish_assembly/02-annot/BRAKER/mask_regions_from_fasta.py \
#     --prefix /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/fClaHyb_Gar_LG.mask_regions
#
# Outputs:
#   <prefix>.normalACGT.bed
#   <prefix>.normalACGT.angsd
#   <prefix>.normalACGT.angsd.idx
#   <prefix>.normalACGT.angsd.bin
###

usage() {
    cat <<EOF
Usage:
  $0 --fasta genome.fa --mask-script mask_regions_from_fasta.py --prefix outprefix
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

INPUT_BED="${PREFIX}.normalACGT.bed"
OUTPUT_ANGSD="${PREFIX}.normalACGT.1-based.angsd"

[[ -f "$INPUT_BED" ]] || { echo "[ERROR] Expected BED not found: $INPUT_BED" >&2; exit 1; }

awk 'BEGIN{OFS="\t"} {print $1, $2+1, $3}' "$INPUT_BED" > "$OUTPUT_ANGSD"

angsd sites index "$OUTPUT_ANGSD"

echo "[INFO] Done."
echo "[INFO] Callable BED : $INPUT_BED"
echo "[INFO] ANGSD sites  : $OUTPUT_ANGSD"
