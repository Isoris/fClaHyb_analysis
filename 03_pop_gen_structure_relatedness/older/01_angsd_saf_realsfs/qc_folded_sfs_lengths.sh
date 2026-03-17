###
# qc_folded_sfs_lengths.sh does this:
# Reports the number of bins (NF) and total sum for each folded SFS file
# to quickly check consistency across chunk-level spectra.
###
BASE="/scratch/lt200308-agbsci/Quentin_project"
SFS_DIR="${BASE}/popstruct_global/02_sfs"

for f in ${SFS_DIR}/catfish.*.sfs; do
  n=$(awk '{print NF}' "$f")
  s=$(awk '{sum=0; for(i=1;i<=NF;i++) sum+=$i; print sum}' "$f")
  echo -e "$(basename "$f")\tNF=$n\tsum=$s"
done | sort -k2,2n | head
