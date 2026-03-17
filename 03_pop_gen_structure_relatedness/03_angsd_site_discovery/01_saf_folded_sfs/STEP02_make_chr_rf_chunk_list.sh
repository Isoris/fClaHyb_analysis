#!/usr/bin/env bash
set -euo pipefail

###
# STEP02_make_chr_rf_chunk_list.sh
#
# Does this:
# - Reads sequence names from genome FASTA
# - Creates one .rf.txt file per sequence
# - Creates chunk_rf.list for SLURM array jobs
#
# Usage:
#   bash STEP02_make_chr_rf_chunk_list.sh --genome genome.fa --output output_dir/
#   bash STEP02_make_chr_rf_chunk_list.sh --genome genome.fa --output output_dir/chunk_rf.list
###

usage() {
    cat <<EOF
Usage:
  $0 --genome genome.fa --output output_dir/
  $0 --genome genome.fa --output chunk_rf.list
EOF
}

GENOME=""
OUTPUT=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --genome) GENOME="$2"; shift 2 ;;
        --output) OUTPUT="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
    esac
done

[[ -n "$GENOME" && -n "$OUTPUT" ]] || { usage >&2; exit 1; }
[[ -f "$GENOME" ]] || { echo "[ERROR] Missing genome: $GENOME" >&2; exit 1; }

if [[ "$OUTPUT" == */ ]] || [[ -d "$OUTPUT" ]]; then
    OUTDIR="${OUTPUT%/}"
    CHUNK_LIST="${OUTDIR}/chunk_rf.list"
else
    OUTDIR="$(dirname "$OUTPUT")"
    CHUNK_LIST="$OUTPUT"
fi

RF_DIR="${OUTDIR}/rf_files"
mkdir -p "$RF_DIR"

awk '
    /^>/ {
        name = substr($1, 2)
        if (name != "") print name
    }
' "$GENOME" | while read -r chr; do
    printf "%s\n" "$chr" > "${RF_DIR}/${chr}.rf.txt"
done

find "$RF_DIR" -maxdepth 1 -type f -name "*.rf.txt" | sort > "$CHUNK_LIST"

echo "[INFO] RF_DIR     : $RF_DIR"
echo "[INFO] CHUNK_LIST : $CHUNK_LIST"
echo "[INFO] N chunks   : $(wc -l < "$CHUNK_LIST")"
head "$CHUNK_LIST"
