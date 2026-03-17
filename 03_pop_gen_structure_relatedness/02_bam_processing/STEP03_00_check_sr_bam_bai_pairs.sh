###
# check_sr_bam_bai_pairs.sh does this:
# Checks that each .sr.bam has a matching .sr.bam.bai and vice versa,
# and reports missing BAM or BAI files.
###
find batch*/ -type f -name "*.sr.bam" ! -name "*.tmp.*" -print \
| sed 's/\.sr\.bam$//' \
| while read -r p; do
    bam="${p}.sr.bam"
    bai="${p}.sr.bam.bai"
    if [[ ! -s "$bai" ]]; then
      echo "MISSING_BAI  $bam"
    fi
  done

find batch*/ -type f -name "*.sr.bam.bai" -print \
| sed 's/\.sr\.bam\.bai$//' \
| while read -r p; do
    bam="${p}.sr.bam"
    bai="${p}.sr.bam.bai"
    if [[ ! -s "$bam" ]]; then
      echo "MISSING_BAM  $bai"
    fi
  done
