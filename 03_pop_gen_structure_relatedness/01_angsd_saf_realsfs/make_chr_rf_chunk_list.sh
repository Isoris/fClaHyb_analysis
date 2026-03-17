###
# make_chr_rf_chunk_list.sh does this:
# Creates one PCAngsd/ANGSD region file (-rf) per chromosome and builds
# a master chunk list for array job indexing.
###
BASE="/scratch/lt200308-agbsci/Quentin_project"
OUTDIR="${BASE}/popstruct_global"
RF_DIR="${OUTDIR}/rf_files"

mkdir -p "$RF_DIR"

# Make one -rf file per chromosome
while read -r chr; do
  printf "%s\n" "$chr" > "${RF_DIR}/${chr}.rf.txt"
done < "${BASE}/fClaHyb_Gar.chromosomes.names.txt"

# List of rf-files (this is what the array will index into)
ls -1 "${RF_DIR}"/*.rf.txt > "${OUTDIR}/chunk_rf.list"

# Check
wc -l "${OUTDIR}/chunk_rf.list"
head "${OUTDIR}/chunk_rf.list"
