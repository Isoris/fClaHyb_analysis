###
# 05B_angsd_beagle_glf2_wholegenome_thin_majmin.slurm does this:
# Runs ANGSD on the whole genome to generate a beagle genotype likelihood file
# using a thinned site list that preserves major/minor alleles.
###
#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=80
#SBATCH --ntasks=1
#SBATCH --mem=237GB
#SBATCH -t 2-00:00:00
#SBATCH -J angsd_beagle_wg
#SBATCH -o angsd_beagle_wg.%j.out
#SBATCH -e angsd_beagle_wg.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

###############################################################################
# STEP 5B) Whole-genome BEAGLE GLs using thinned sites WITH major/minor columns
#
# Usage:
#   sbatch 05B_angsd_beagle_glf2_wholegenome_thin_majmin.slurm 10000
#
# Args:
#   $1 = thinning bp (5000 / 10000 / 25000)
###############################################################################

W="${1:?ERROR: thinning bp (5000/10000/25000) as arg1}"

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
REF="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
BAMLIST="${BASE}/bamlist.pp.samechr.tlenP99.filtered.txt"

OUTDIR="${BASE}/popstruct_thin"
SITES_FINAL="${OUTDIR}/03_sites/sites.thin_${W}.majmin.tsv"

BEAGLEDIR="${OUTDIR}/04_beagle_wholegenome_majmin"
LOGDIR="${OUTDIR}/logs"
THREADS="${SLURM_CPUS_PER_TASK}"

mkdir -p "${BEAGLEDIR}" "${LOGDIR}"

# Filters
MINQ=25
MINMAPQ=25
MINDEPTHIND=3
MAXDEPTHIND=57
MININD=200

[[ -s "${SITES_FINAL}" ]] || { echo "[ERROR] Missing sites file: ${SITES_FINAL}" >&2; exit 1; }

PREFIX="${BEAGLEDIR}/catfish.wholegenome.thin_${W}"
LOG="${LOGDIR}/05B_angsd_beagle.wholegenome.thin_${W}.log"

echo "[INFO] SITES_FINAL=${SITES_FINAL}"
echo "[INFO] OUT=${PREFIX}.beagle.gz"

set -x
angsd -b "${BAMLIST}" -ref "${REF}" \
  -GL 1 \
  -doMajorMinor 3 \
  -doGlf 2 \
  -doCounts 1 -doDepth 1 \
  -minQ "${MINQ}" -minMapQ "${MINMAPQ}" -baq 1 -C 50 \
  -setMinDepthInd "${MINDEPTHIND}" -setMaxDepthInd "${MAXDEPTHIND}" \
  -minInd "${MININD}" \
  -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 \
  -sites "${SITES_FINAL}" \
  -P "${THREADS}" \
  -out "${PREFIX}" \
  > "${LOG}" 2>&1
set +x

echo "[DONE] ${PREFIX}.beagle.gz"
