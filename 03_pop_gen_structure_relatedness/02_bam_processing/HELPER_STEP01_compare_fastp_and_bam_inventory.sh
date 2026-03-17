###
# compare_fastp_and_bam_inventory.sh does this:
# Compares fastp output pairs against generated .sr.bam files and their indexes,
# then reports whether each sample is complete, missing, or inconsistent.
###
FASTP_DIR="/scratch/lt200308-agbsci/Quentin_project/00-samples/fastp"
BAMROOT="/scratch/lt200308-agbsci/Quentin_project/01-minimap2-bams"
OUT="fastp_vs_bam.tsv"

echo -e "id\tfastp_R1\tfastp_R2\tfastp_pair_ok\tbam_path\tbam_ok\tbai_ok\tstatus" > "$OUT"

# 1) collect fastp IDs (must have both R1 and R2)
tmp_fastp=$(mktemp)
find "$FASTP_DIR" -maxdepth 1 -type f -name "*.R1.fastp.fq.gz" | sort \
| sed 's#.*/##; s/\.R1\.fastp\.fq\.gz$//' > "$tmp_fastp"

# 2) collect bam IDs (ignore tmp chunks)
tmp_bam=$(mktemp)
find "$BAMROOT"/batch*/ -type f -name "*.sr.bam" ! -name "*.tmp.*" | sort \
| sed 's#.*/##; s/\.sr\.bam$//' > "$tmp_bam"

# 3) union of IDs
tmp_all=$(mktemp)
cat "$tmp_fastp" "$tmp_bam" | sort -u > "$tmp_all"

while read -r id; do
  r1="${FASTP_DIR}/${id}.R1.fastp.fq.gz"
  r2="${FASTP_DIR}/${id}.R2.fastp.fq.gz"
  r1_ok=0; [[ -s "$r1" ]] && r1_ok=1
  r2_ok=0; [[ -s "$r2" ]] && r2_ok=1
  pair_ok=0; [[ "$r1_ok" -eq 1 && "$r2_ok" -eq 1 ]] && pair_ok=1

  bam="$(find "$BAMROOT"/batch*/ -type f -name "${id}.sr.bam" ! -name "*.tmp.*" | head -n 1)"
  bam_ok=0; [[ -n "$bam" && -s "$bam" ]] && bam_ok=1
  bai_ok=0
  if [[ "$bam_ok" -eq 1 && -s "${bam}.bai" ]]; then
    bai_ok=1
  fi

  status="OK"
  if [[ "$pair_ok" -eq 1 && "$bam_ok" -eq 0 ]]; then
    status="FASTP_OK_BUT_NO_BAM"
  elif [[ "$pair_ok" -eq 0 && "$bam_ok" -eq 1 ]]; then
    status="BAM_OK_BUT_NO_FASTP_PAIR"
  elif [[ "$bam_ok" -eq 1 && "$bai_ok" -eq 0 ]]; then
    status="BAM_OK_BUT_MISSING_BAI"
  elif [[ "$pair_ok" -eq 0 && "$bam_ok" -eq 0 ]]; then
    status="NEITHER_PRESENT"
  fi

  echo -e "${id}\t${r1}\t${r2}\t${pair_ok}\t${bam:-NA}\t${bam_ok}\t${bai_ok}\t${status}"
done < "$tmp_all" >> "$OUT"

rm -f "$tmp_fastp" "$tmp_bam" "$tmp_all"
echo "Wrote: $OUT"
