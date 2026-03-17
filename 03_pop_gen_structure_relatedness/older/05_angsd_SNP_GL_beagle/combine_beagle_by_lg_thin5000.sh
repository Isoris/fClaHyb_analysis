###
# combine_beagle_by_lg_thin5000.sh does this:
# Combines per-linkage-group beagle files into one whole-genome beagle file
# for the 5000 bp thinning set.
###
OUT=catfish.wholegenome.byRF.thin_5000.beagle.gz

# sort ensures LG01..LG28 order (important for reproducibility)
ls -1 catfish.C_gar_LG*.rf.txt.thin_5000.beagle.gz | sort -V > beagle_thin_5000.list

# write header once + append data lines
first=$(head -n1 beagle_thin_5000.list)

zcat "$first" | head -n1 > catfish.wholegenome.byRF.thin_5000.beagle

while read -r f; do
  zcat "$f" | tail -n +2 >> catfish.wholegenome.byRF.thin_5000.beagle
done < beagle_thin_5000.list

gzip -c catfish.wholegenome.byRF.thin_5000.beagle > "$OUT"
rm -f catfish.wholegenome.byRF.thin_5000.beagle
