###
# combine_beagle_by_lg_thin10000.sh does this:
# Combines per-linkage-group beagle files into one whole-genome beagle file
# for the 10000 bp thinning set.
###
OUT=catfish.wholegenome.byRF.thin_10000.beagle.gz
ls -1 catfish.C_gar_LG*.rf.txt.thin_10000.beagle.gz | sort -V > beagle_thin_10000.list
first=$(head -n1 beagle_thin_10000.list)

zcat "$first" | head -n1 > catfish.wholegenome.byRF.thin_10000.beagle
while read -r f; do
  zcat "$f" | tail -n +2 >> catfish.wholegenome.byRF.thin_10000.beagle
done < beagle_thin_10000.list

gzip -c catfish.wholegenome.byRF.thin_10000.beagle > "$OUT"
rm -f catfish.wholegenome.byRF.thin_10000.beagle
