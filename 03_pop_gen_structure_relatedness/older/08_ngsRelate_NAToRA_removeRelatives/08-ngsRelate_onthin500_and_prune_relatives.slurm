#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 80
#SBATCH --mem=237GB
#SBATCH -t 23:00:00
#SBATCH -J ngsRelate
#SBATCH -o ngsRelate.%j.out
#SBATCH -e ngsRelate.%j.err



set -euo pipefail
source ~/.bashrc
mamba activate assembly

module load PrgEnv-gnu
module load gcc/11.2.0
module load zlib/1.2.13
module load bzip2/1.0.8-cpeGNU-23.03
module load XZ/5.4.3-cpeGNU-23.03
module load cURL/8.1.2
module unload GSL 2>/dev/null || true
module load GSL/2.8-cpeGNU-25.03


BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
BEAGLE="${BASE}/popstruct_thin/04_beagle_byRF_majmin/catfish.wholegenome.byRF.thin_500.beagle.gz"
MAF_DIR="${BASE}/popstruct_global/02_snps"
IDFILE="${BASE}/popstruct_thin/list_of_samples_one_per_line_same_bamfile_list.tsv"
OUTDIR="${BASE}/popstruct_thin/05_ngsrelate"
mkdir -p "${OUTDIR}"

# --- 1. Build merged frequency file from all mafs.gz, keyed by chr_pos ---
echo "[INFO] Merging mafs.gz frequencies..."
zcat ${MAF_DIR}/catfish.*.mafs.gz \
  | awk 'BEGIN{OFS="\t"} $1!="chromo" {print $1"_"$2, $6}' \
  | sort -k1,1 \
  > "${OUTDIR}/all_mafs_freq.tmp"

# --- 2. Extract site order from beagle file (marker column) ---
echo "[INFO] Extracting beagle site order..."
zcat "${BEAGLE}" | tail -n +2 | cut -f1 > "${OUTDIR}/beagle_sites.tmp"

NSITES=$(wc -l < "${OUTDIR}/beagle_sites.tmp")
echo "[INFO] Number of sites in beagle: ${NSITES}"

# --- 3. Match frequencies to beagle site order ---
echo "[INFO] Matching frequencies..."
awk 'NR==FNR {freq[$1]=$2; next} {if($1 in freq) print freq[$1]; else print "NA"}' \
  "${OUTDIR}/all_mafs_freq.tmp" \
  "${OUTDIR}/beagle_sites.tmp" \
  > "${OUTDIR}/freq_for_ngsrelate.txt"

NA_COUNT=$(grep -c "^NA$" "${OUTDIR}/freq_for_ngsrelate.txt" || true)
echo "[INFO] Sites with missing frequency: ${NA_COUNT} / ${NSITES}"
if [[ "${NA_COUNT}" -gt 0 ]]; then
  echo "[WARN] Some sites have no matching frequency — check marker format consistency"
fi

# --- 4. Run ngsRelate ---
echo "[INFO] Running ngsRelate on 226 samples, ${NSITES} sites..."
ngsRelate \
  -G "${BEAGLE}" \
  -f "${OUTDIR}/freq_for_ngsrelate.txt" \
  -n 226 \
  -O "${OUTDIR}/catfish_226_relatedness.res" \
  -p 80 \
  -m 1 \
  -z "${IDFILE}"

echo "[DONE] Results in ${OUTDIR}/catfish_226_relatedness.res"
