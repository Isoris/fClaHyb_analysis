###
# make_pcaangsd_per_lg_lists_thin5000_10000_25000.sh does this:
# Builds validated beagle file lists for per-linkage-group PCAngsd analyses
# at 5000 bp, 10000 bp, and 25000 bp thinning distances.
###
cd /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/popstruct_thin

mkdir -p 05_pcangsd_byLG

for thin in 5000 10000 25000; do
  out="05_pcangsd_byLG/beagle_LG_thin_${thin}.list"
  : > "$out"
  for f in 04_beagle_byRF_majmin/catfish.C_gar_LG*.rf.txt.thin_${thin}.beagle.gz; do
    if gzip -t "$f" 2>/dev/null; then
      echo "$f" >> "$out"
    else
      echo "[BAD] $f" >&2
    fi
  done
  echo "thin_${thin}: $(wc -l < "$out") files -> $out"
done
