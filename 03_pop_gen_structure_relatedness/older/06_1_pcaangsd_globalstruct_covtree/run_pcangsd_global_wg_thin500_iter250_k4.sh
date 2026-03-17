###
# run_pcangsd_global_wg_thin500_iter250_k4.sh does this:
# Runs whole-genome PCAngsd on the 500 bp thinned beagle file
# with 250 iterations and admixture K=4, then prints a short log tail.
###
source ~/.bashrc
mamba activate assembly

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/popstruct_thin"
IN="${BASE}/04_beagle_byRF_majmin/catfish.wholegenome.byRF.thin_500.beagle.gz"
SAMPLES="${BASE}/list_of_samples_one_per_line_same_bamfile_list.tsv"
OUTDIR="${BASE}/05_pcangsd_global_WGthin/iter_250_thin_500/K4_cmdline"
mkdir -p "$OUTDIR"
OUT="${OUTDIR}/catfish.wg.byRF.thin_500.pcangsd"

pcangsd \
  -b "$IN" \
  -e 4 \
  -t 32 \
  --maf 0.05 \
  --iter 250 \
  --admix \
  --admix-K 4 \
  --admix-seed 1 \
  --tree \
  --tree-samples "$SAMPLES" \
  -o "$OUT" \
  > "${OUTDIR}/pcangsd_run.log" 2>&1

tail -n 30 "${OUTDIR}/pcangsd_run.log"
ls -lh "${OUTDIR}"
