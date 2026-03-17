###
# 04_distance_thin_200b_500b_1kb_sites_v3.slurm does this:
# Merges SNP sites from chunked .mafs.gz files while preserving major/minor alleles,
# applies minimum-distance thinning at 200 bp, 500 bp, and 1 kb, and ANGSD-indexes
# the resulting site files for downstream analyses that reuse the same allele coding.
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
# STEP 1: Merge SNP sites from *.mafs.gz (KEEP major/minor), thin at 5/10/25kb,
#         and ANGSD-index them (so -sites reuses the same allele coding)
###############################################################################

# --------- EDIT THESE PATHS ----------
BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
MAF_DIR="${BASE}/popstruct_global/02_snps"     # where catfish.*.mafs.gz live
OUTDIR="${BASE}/popstruct_thin"               # output directory for merged/thinned sites
# ------------------------------------

mkdir -p "${OUTDIR}/03_sites"

# 4-col sites file: chr  pos  major  minor
RAW_ALL="${OUTDIR}/03_sites/sites.raw.ALL.majmin.tsv"

MAF_GLOB="${MAF_DIR}/catfish.*.mafs.gz"

# sanity: make sure glob matches files
NFILES=$(ls -1 ${MAF_GLOB} 2>/dev/null | wc -l)
if [[ "$NFILES" -eq 0 ]]; then
  echo "[ERROR] No .mafs.gz matched: ${MAF_GLOB}" >&2
  exit 1
fi
echo "[INFO] Found ${NFILES} mafs.gz files"

# Merge all .mafs.gz -> 4-column (chr, pos, major, minor) skipping repeated headers
# Filter: alleles must be single bases, non-N, and major!=minor
zcat ${MAF_GLOB} \
  | awk 'BEGIN{OFS="\t"} $1!="chromo" {print $1,$2,$3,$4}' \
  | awk 'BEGIN{OFS="\t"}
         length($3)==1 && length($4)==1 &&
         $3!="N" && $3!="n" && $4!="N" && $4!="n" &&
         $3!=$4 {print}' \
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
# Confirm alleles look OK
BADAL=$(awk 'length($3)!=1 || length($4)!=1 || $3==$4 {c++} END{print c+0}' "${RAW_ALL}")
if [[ "${BADPOS}" -ne 0 || "${BADAL}" -ne 0 ]]; then
  echo "[ERROR] Found BADPOS=${BADPOS} or BADAL=${BADAL} in ${RAW_ALL}" >&2
  exit 1
fi

# 2) Thinning function: true minimum distance (>= W) within each chromosome
#    Keeps the major/minor columns aligned with the chosen positions.
thin_min_dist () {
  local W="$1"
  local OUT="${OUTDIR}/03_sites/sites.thin_${W}.majmin.tsv"

  awk -v W="${W}" 'BEGIN{OFS="\t"}
    {chr=$1; pos=$2; maj=$3; min=$4;
     if(chr!=pchr){pchr=chr; last=-1e18}
     if(pos-last>=W){print chr,pos,maj,min; last=pos}
    }' "${RAW_ALL}" > "${OUT}"

  echo "[INFO] Thinned ${W} bp site count: $(wc -l < "${OUT}") -> ${OUT}"
}

thin_min_dist 200
thin_min_dist 500
thin_min_dist 1000

# 3) ANGSD index each thinned sites file so it can be used efficiently with -sites
#    This produces: .bin and .idx alongside each .tsv
index_sites () {
  local F="$1"
  echo "[INFO] Indexing sites with ANGSD: ${F}"
  angsd sites index "${F}"
}

index_sites "${OUTDIR}/03_sites/sites.thin_200.majmin.tsv"
index_sites "${OUTDIR}/03_sites/sites.thin_500.majmin.tsv"
index_sites "${OUTDIR}/03_sites/sites.thin_1000.majmin.tsv"

echo "[DONE] Step 1 complete."
echo "       Raw:     ${RAW_ALL}"
echo "       Thinned: ${OUTDIR}/03_sites/sites.thin_200.majmin.tsv"
echo "               ${OUTDIR}/03_sites/sites.thin_500.majmin.tsv"
echo "               ${OUTDIR}/03_sites/sites.thin_1000.majmin.tsv"
echo "       Indexed: run 'ls -lh ${OUTDIR}/03_sites/sites.thin_*.majmin.tsv*' to see .bin/.idx"
