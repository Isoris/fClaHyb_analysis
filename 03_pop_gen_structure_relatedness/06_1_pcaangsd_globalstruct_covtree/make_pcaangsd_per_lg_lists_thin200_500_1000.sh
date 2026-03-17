###
# make_pcaangsd_per_lg_lists_thin200_500_1000.sh does this:
# Builds validated beagle file lists for per-linkage-group PCAngsd analyses
# at 200 bp, 500 bp, and 1000 bp thinning distances.
###
cd /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/popstruct_thin

mkdir -p 05_pcangsd_byLG

for thin in 200 500 1000; do
  out="05_pcangsd_byLG/beagle_LG_thin_${thin}.list"
  : > "$out"
  for f in 04_beagle_byRF_majmin/*/catfish.C_gar_LG*.rf.txt.thin_${thin}.beagle.gz; do
    if gzip -t "$f" 2>/dev/null; then
      echo "$f" >> "$out"
    else
      echo "[BAD] $f" >&2
    fi
  done
  echo "thin_${thin}: $(wc -l < "$out") files -> $out"
done
