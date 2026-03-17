###
# summarize_snp_stats_for_manuscript.sh does this:
# Summarizes raw and thinned SNP sets for manuscript writing, including per-chromosome
# counts, SNP density, Ti/Tv ratio, allele composition, and draft copy-paste text.
###
#!/usr/bin/env bash
# snp_stats_for_manuscript.sh
# Run from: /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/popstruct_thin/03_sites
# Or pass the sites directory as $1
#
# Produces a summary table for the Methods/Results text

set -euo pipefail

SITESDIR="${1:-.}"
RAW="${SITESDIR}/sites.raw.ALL.majmin.tsv"

echo "=============================================="
echo " SNP statistics for manuscript"
echo "=============================================="
echo ""

# --- 1. Total biallelic SNPs (raw, before thinning) ---
TOTAL_RAW=$(wc -l < "$RAW")
echo "Total biallelic SNPs (raw, all chromosomes): ${TOTAL_RAW}"
echo ""

# --- 2. Per-chromosome SNP counts ---
echo "--- Per-chromosome SNP counts (raw) ---"
awk '{print $1}' "$RAW" | sort | uniq -c | sort -k2,2V
echo ""

# --- 3. Number of chromosomes with SNPs ---
NCHR=$(awk '{print $1}' "$RAW" | sort -u | wc -l)
echo "Number of chromosomes/scaffolds with SNPs: ${NCHR}"
echo ""

# --- 4. SNP density (SNPs per kb) ---
FAI=""
for candidate in \
  /scratch/lt200308-agbsci/Quentin_project/00-samples/fClaHyb_Gar_LG.fa.fai \
  /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/00-samples/fClaHyb_Gar_LG.fa.fai; do
  if [[ -s "$candidate" ]]; then
    FAI="$candidate"
    break
  fi
done

if [[ -n "$FAI" ]]; then
  echo "--- Reference genome size (from $FAI) ---"
  GENOME_SIZE=$(awk '{s+=$2} END{print s}' "$FAI")
  echo "Total reference size: ${GENOME_SIZE} bp ($(awk -v g="$GENOME_SIZE" 'BEGIN{printf "%.2f Mb", g/1e6}'))"
  SNP_PER_KB=$(awk -v n="$TOTAL_RAW" -v g="$GENOME_SIZE" 'BEGIN{printf "%.2f", n/(g/1000)}')
  echo "SNP density (genome-wide): ${SNP_PER_KB} SNPs/kb"
  echo ""

  echo "--- Per-chromosome SNP density ---"
  echo -e "chr\tchr_length_bp\tn_snps\tsnps_per_kb"
  awk '{print $1}' "$RAW" | sort | uniq -c | awk '{print $2"\t"$1}' | sort -k1,1V > /tmp/_snp_counts.tmp
  awk '{print $1"\t"$2}' "$FAI" | sort -k1,1V > /tmp/_chr_lengths.tmp
  join -t$'\t' -1 1 -2 1 /tmp/_chr_lengths.tmp /tmp/_snp_counts.tmp | \
    awk -F'\t' '{printf "%s\t%s\t%s\t%.2f\n", $1, $2, $3, $3/($2/1000)}'
  rm -f /tmp/_snp_counts.tmp /tmp/_chr_lengths.tmp
  echo ""
else
  echo "[WARN] Could not find .fai file — skipping density calculations"
  echo "       Set FAI path manually if needed"
  echo ""
fi

# --- 5. Thinned site counts ---
echo "--- Thinned site counts ---"
echo -e "thinning_bp\tn_sites\tfile"
for W in 200 500 1000 5000 10000 25000; do
  F="${SITESDIR}/sites.thin_${W}.majmin.tsv"
  if [[ -s "$F" ]]; then
    N=$(wc -l < "$F")
    echo -e "${W}\t${N}\t$(basename $F)"
  else
    echo -e "${W}\tNA\t(not found)"
  fi
done
echo ""

# --- 6. Allele composition check ---
echo "--- Allele composition (major/minor) ---"
echo "Unique major alleles:"
awk '{print $3}' "$RAW" | sort | uniq -c | sort -rn
echo ""
echo "Unique minor alleles:"
awk '{print $4}' "$RAW" | sort | uniq -c | sort -rn
echo ""

# --- 7. Transition/transversion ratio ---
echo "--- Transition / Transversion ---"
awk '{
  a=$3; b=$4
  if (a > b) { tmp=a; a=b; b=tmp }
  pair = a">"b
  if (pair == "A>G" || pair == "C>T") ti++
  else tv++
  counts[pair]++
}
END {
  print "Transitions (Ti): " ti+0
  print "Transversions (Tv): " tv+0
  if (tv > 0) printf "Ti/Tv ratio: %.3f\n", ti/tv
  print ""
  print "Per-pair counts:"
  for (p in counts) print "  " p ": " counts[p]
}' "$RAW"
echo ""

# --- 8. Quick manuscript sentence generator ---
echo "=============================================="
echo " COPY-PASTE for manuscript"
echo "=============================================="

N200=$(wc -l < "${SITESDIR}/sites.thin_200.majmin.tsv" 2>/dev/null || echo "NA")
N500=$(wc -l < "${SITESDIR}/sites.thin_500.majmin.tsv" 2>/dev/null || echo "NA")
N1000=$(wc -l < "${SITESDIR}/sites.thin_1000.majmin.tsv" 2>/dev/null || echo "NA")
N5000=$(wc -l < "${SITESDIR}/sites.thin_5000.majmin.tsv" 2>/dev/null || echo "NA")
N10000=$(wc -l < "${SITESDIR}/sites.thin_10000.majmin.tsv" 2>/dev/null || echo "NA")
N25000=$(wc -l < "${SITESDIR}/sites.thin_25000.majmin.tsv" 2>/dev/null || echo "NA")

echo ""
echo "After quality filtering, we identified ${TOTAL_RAW} biallelic SNPs"
echo "across ${NCHR} chromosomes (MAF >= 0.05, p < 1e-6),"
if [[ -n "$FAI" ]]; then
  echo "corresponding to ${SNP_PER_KB} SNPs per kb."
fi
echo ""
echo "Linkage-thinned sets:"
echo "  200 bp:    ${N200} SNPs"
echo "  500 bp:    ${N500} SNPs"
echo "  1,000 bp:  ${N1000} SNPs"
echo "  5,000 bp:  ${N5000} SNPs"
echo "  10,000 bp: ${N10000} SNPs"
echo "  25,000 bp: ${N25000} SNPs"
echo ""
echo "Draft sentence:"
echo "\"After quality filtering, we identified $(printf "%'d" $TOTAL_RAW) biallelic SNPs"
echo "across ${NCHR} chromosomes (likelihood ratio test p < 10^-6, MAF >= 0.05),"
if [[ -n "$FAI" ]]; then
  echo "corresponding to ${SNP_PER_KB} SNPs per kb."
fi
echo "To reduce the effects of physical linkage, SNP sets were thinned"
echo "at minimum inter-site distances of 200 bp ($(printf "%'d" $N200) SNPs),"
echo "500 bp ($(printf "%'d" $N500)), 1 kb ($(printf "%'d" $N1000)),"
echo "5 kb ($(printf "%'d" $N5000)), 10 kb ($(printf "%'d" $N10000))"
echo "and 25 kb ($(printf "%'d" $N25000)) for structure and per-chromosome analyses.\""
echo ""

# --- 9. BAM list stats (if accessible) ---
for BAMLIST in \
  /scratch/lt200308-agbsci/Quentin_project/bamlist.pp.samechr.tlenP99.filtered.txt \
  /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/bamlist.pp.samechr.tlenP99.filtered.txt; do
  if [[ -s "$BAMLIST" ]]; then
    NBAM=$(wc -l < "$BAMLIST")
    echo "--- BAM list: $BAMLIST ---"
    echo "Total individuals in BAM list: ${NBAM}"
    echo ""
    break
  fi
done

# --- 10. Sample list stats ---
SAMPLES="${SITESDIR}/../list_of_samples_one_per_line_same_bamfile_list.tsv"
if [[ -s "$SAMPLES" ]]; then
  NSAMP=$(wc -l < "$SAMPLES")
  echo "--- Sample list ---"
  echo "Total samples: ${NSAMP}"
  echo "First 5:"
  head -5 "$SAMPLES"
  echo "..."
  echo ""
fi

echo "[DONE]"
