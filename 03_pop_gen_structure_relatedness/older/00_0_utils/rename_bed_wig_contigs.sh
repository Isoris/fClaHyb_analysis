###
# rename_bed_wig_contigs.sh does this:
# Renames chromosome/contig names inside BED and WIG files using a generated sed map,
# then checks whether renamed contigs match the allowed contigs from the Mac and Gar FASTA headers.
###
#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# ---- paths ----
WORKDIR="/scratch/lt200308-agbsci/Quentin_project/00-samples"
MAC_FA="${WORKDIR}/fClaHyb_Mac_LG.fa"
GAR_FA="${WORKDIR}/fClaHyb_Gar_LG.fa"

cd "$WORKDIR"

# ---- inputs (edit if you have more files) ----
FILES=(
  "Hybrid_catfish_polished_19_k31_only.bed"
  "Hybrid_catfish_polished_19_k31_only.wig"
  "Hybrid_catfish_polished_19_only.bed"
  "Hybrid_catfish_polished_19_only.wig"
)

# ---- outputs ----
OUTDIR="${WORKDIR}/renamed_contigs"
mkdir -p "$OUTDIR"

# ---- build allowed contig lists from fasta headers ----
grep '^>' "$MAC_FA" | sed 's/^>//' | sort -u > "${OUTDIR}/mac.contigs.txt"
grep '^>' "$GAR_FA" | sed 's/^>//' | sort -u > "${OUTDIR}/gar.contigs.txt"

echo "[INFO] Mac contigs: $(wc -l < "${OUTDIR}/mac.contigs.txt")"
echo "[INFO] Gar contigs: $(wc -l < "${OUTDIR}/gar.contigs.txt")"

# ---- helper: zero-pad to 2 digits ----
pad2() { printf "%02d" "$1"; }

# ---- make sed script ----
SED_SCRIPT="${OUTDIR}/rename.sed"
: > "$SED_SCRIPT"

# bighead -> mac (01..27)
for i in $(seq 1 27); do
  ii=$(pad2 "$i")
  echo "s/\bfClaHyb_bighead_1_Chr_${ii}\b/C_mac_LG${ii}/g" >> "$SED_SCRIPT"
done

# african -> gar (01..28)
for i in $(seq 1 28); do
  ii=$(pad2 "$i")
  echo "s/\bfClaHyb_african_1_Chr_${ii}\b/C_gar_LG${ii}/g" >> "$SED_SCRIPT"
done

echo "[INFO] Wrote sed rules: $SED_SCRIPT"

# ---- process files ----
for f in "${FILES[@]}"; do
  if [[ ! -s "$f" ]]; then
    echo "[WARN] missing: $f (skipping)"
    continue
  fi

  b=$(basename "$f")
  out="${OUTDIR}/${b%.*}.renamed.${b##*.}"
  cp -a "$f" "${OUTDIR}/${b}.bak"

  sed -E -f "$SED_SCRIPT" "$f" > "$out"
  echo "[INFO] wrote: $out"
done

# ---- sanity checks ----
echo
echo "[CHECK] Contigs present after renaming (BED):"
for f in "${OUTDIR}"/*.bed; do
  [[ -e "$f" ]] || continue
  echo "  - $(basename "$f")"
  cut -f1 "$f" | sort -u | head
done

echo
echo "[CHECK] Contigs present after renaming (WIG):"
for f in "${OUTDIR}"/*.wig; do
  [[ -e "$f" ]] || continue
  echo "  - $(basename "$f")"
  grep -E '^variableStep' "$f" | sed -E 's/.*chrom=([^\t ]+).*/\1/' | sort -u | head
done

# ---- find any contigs not in either fasta list ----
ALL_ALLOWED="${OUTDIR}/allowed.contigs.txt"
cat "${OUTDIR}/mac.contigs.txt" "${OUTDIR}/gar.contigs.txt" | sort -u > "$ALL_ALLOWED"

echo
echo "[CHECK] Any contigs in outputs NOT found in either FASTA?"
for f in "${OUTDIR}"/*.bed "${OUTDIR}"/*.wig; do
  [[ -e "$f" ]] || continue
  if [[ "$f" == *.bed ]]; then
    cut -f1 "$f" | sort -u
  else
    grep -E '^variableStep' "$f" | sed -E 's/.*chrom=([^\t ]+).*/\1/' | sort -u
  fi | comm -23 - "$ALL_ALLOWED" | sed 's/^/  UNMATCHED: /' || true
done

echo
echo "[DONE] Renamed files are in: $OUTDIR"
