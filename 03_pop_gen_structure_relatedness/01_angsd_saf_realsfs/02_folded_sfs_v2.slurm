###
# 02_folded_sfs_v2.slurm does this:
# Runs realSFS on the global SAF index to compute a folded bootstrap SFS
# and writes the folded spectrum plus the realSFS log.
###
#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 124
#SBATCH --mem=237GB
#SBATCH -t 2-00:00:00
#SBATCH -J realSFS_folded
#SBATCH -o realSFS_folded.%j.out
#SBATCH -e realSFS_folded.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

realSFS global_sfs/catfish.global.saf.idx -P 124 -bootstrap 15 -fold 1 \
  > global_sfs/catfish.global.folded.sfs \
  2> global_sfs/realSFS.global.folded.log
