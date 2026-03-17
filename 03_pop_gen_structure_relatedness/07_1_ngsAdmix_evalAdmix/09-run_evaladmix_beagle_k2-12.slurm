###
# run_evaladmix_thin500_k2_k12.slurm does this:
# Runs evalAdmix on the thin-500 whole-genome beagle file for K=2 to K=12
# across three seeds, using existing NGSadmix fopt/qopt outputs.
###
#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1 -n 127
#SBATCH --mem=237GB
#SBATCH -t 08:00:00
#SBATCH -J evaladmix
#SBATCH -o evaladmix.%j.out
#SBATCH -e evaladmix.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
EVALADMIX_BIN="${BASE}/evalAdmix"   # edit if needed
RUN_DIR="${BASE}/popstruct_thin/05_ngsadmix_global/runs_thin500"
OUTDIR="${BASE}/popstruct_thin/05_ngsadmix_global/evaladmix_thin500"
BEAGLE="${BASE}/popstruct_thin/04_beagle_byRF_majmin/catfish.wholegenome.byRF.thin_500.beagle.gz"

mkdir -p "${OUTDIR}"

for K in $(seq 2 12); do
  KPAD=$(printf "%02d" "${K}")
  for SEED in 1 2 3; do
    PREFIX="thin500_K${KPAD}_seed${SEED}"

    "${EVALADMIX_BIN}/evalAdmix" \
      -beagle "${BEAGLE}" \
      -fname  "${RUN_DIR}/${PREFIX}.fopt.gz" \
      -qname  "${RUN_DIR}/${PREFIX}.qopt" \
      -P 110 \
      -nIts 5 \
      -misTol 0.05 \
      -minMaf 0.05 \
      -o "${OUTDIR}/${PREFIX}.corres.txt"
  done
done
