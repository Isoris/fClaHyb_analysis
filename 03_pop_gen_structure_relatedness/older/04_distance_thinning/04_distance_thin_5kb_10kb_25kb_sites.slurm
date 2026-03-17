###
# 04_distance_thin_5kb_10kb_25kb_sites.slurm does this:
# Merges SNP sites from chunked .mafs.gz files, creates genome-wide site lists,
# applies minimum-distance thinning at 5 kb, 10 kb, and 25 kb, and ANGSD-indexes
# the resulting site files for reuse with -sites.
###
#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mem=16GB
#SBATCH -t 02:00:00
#SBATCH -J make_thin_sites
#SBATCH -o make_thin_sites.%j.out
#SBATCH -e make_thin_sites.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

###############################################################################
# STEP 1: Merge SNP sites from *.mafs.gz, thin at 5/10/25kb, and ANGSD-index them
###############################################################################

# --------- EDIT THESE PATHS ----------
BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
MAF_DIR="${BASE}/popstruct_global/02_snps"     # where catfish.*.mafs.gz live
OUTDIR="${BASE}/popstruct_thin"               # output directory for merged/thinned sites
# ------------------------------------

mkdir -p "${OUTDIR}/03_sites"

RAW_ALL="${OUTDIR}/03_sites/sites.raw.ALL.tsv"

MAF_GLOB="${BASE}/popstruct_global/02_snps/catfish.*.mafs.gz"

# sanity: make sure glob matches files
NFILES=$(ls -1 ${MAF_GLOB} 2>/dev/null | wc -l)
if [[ "$NFILES" -eq 0 ]]; then
  echo "[ERROR] No .mafs.gz matched: ${MAF_GLOB}" >&2
  exit 1
fi
echo "[INFO] Found ${NFILES} mafs.gz files"

# Merge all .mafs.gz -> 2-column (chr, pos) skipping repeated headers
zcat ${MAF_GLOB} \
  | awk 'BEGIN{OFS="\t"} $1!="chromo" {print $1,$2}' \
  | sort -k1,1 -k2,2n \
  > "${RAW_ALL}"

# sanity checks
echo "[INFO] Raw merged site count: $(wc -l < "${RAW_ALL}")"
echo "[INFO] First 3 lines:"
head -n 3 "${RAW_ALL}"
echo "[INFO] Last 3 lines:"
tail -n 3 "${RAW_ALL}"

# Confirm all positions are >= 1
BADPOS=$(awk '$2<1{c++} END{print c+0}' "${RAW_ALL}")
if [[ "${BADPOS}" -ne 0 ]]; then
  echo "[ERROR] Found ${BADPOS} positions < 1. Something is wrong (positions should be 1-based)."
  exit 1
fi

# 2) Thinning function: true minimum distance (>= W) within each chromosome
thin_min_dist () {
  local W="$1"
  local OUT="${OUTDIR}/03_sites/sites.thin_${W}.tsv"

  awk -v W="${W}" 'BEGIN{OFS="\t"}
    {chr=$1; pos=$2;
     if(chr!=pchr){pchr=chr; last=-1e18}
     if(pos-last>=W){print chr,pos; last=pos}
    }' "${RAW_ALL}" > "${OUT}"

  echo "[INFO] Thinned ${W} bp site count: $(wc -l < "${OUT}") -> ${OUT}"
}

thin_min_dist 5000
thin_min_dist 10000
thin_min_dist 25000

# 3) ANGSD index each thinned sites file so it can be used efficiently with -sites
#    This produces: sites.thin_W.tsv.bin and sites.thin_W.tsv.idx (names may vary by ANGSD build)
index_sites () {
  local F="$1"
  echo "[INFO] Indexing sites with ANGSD: ${F}"
  angsd sites index "${F}"
}

index_sites "${OUTDIR}/03_sites/sites.thin_5000.tsv"
index_sites "${OUTDIR}/03_sites/sites.thin_10000.tsv"
index_sites "${OUTDIR}/03_sites/sites.thin_25000.tsv"

echo "[DONE] Step 1 complete."
echo "       Raw:     ${RAW_ALL}"
echo "       Thinned: ${OUTDIR}/03_sites/sites.thin_5000.tsv"
echo "               ${OUTDIR}/03_sites/sites.thin_10000.tsv"
echo "               ${OUTDIR}/03_sites/sites.thin_25000.tsv"
echo "       Indexed: run 'ls -lh ${OUTDIR}/03_sites/sites.thin_*.tsv*' to see .bin/.idx"
