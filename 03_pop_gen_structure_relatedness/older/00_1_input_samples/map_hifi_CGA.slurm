###
# map_hifi_CGA.slurm does this:
# Maps HiFi reads to the reference genome with minimap2, sorts and indexes
# the alignment with samtools, and writes the final BAM/BAI in scratch space.
###
#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 127
#SBATCH --mem=237GB
#SBATCH -t 2-00:00:00
##SBATCH -A lt200308
#SBATCH -J map_hifi
#SBATCH -o map_hifi.%j.out
#SBATCH -e map_hifi.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

# -----------------------
# Inputs (EDIT IF NEEDED)
# -----------------------
READS="/project/lt200308-agbsci/01-catfish_assembly/00-reads/hifi/CG_F_HIFI_None_reads.fastq"

# Put the reference you want "together in scratch" here:
REF="/scratch/lt200308-agbsci/Quentin_project/00-samples/fClaHyb_Gar_LG.fa"

# Output in scratch (as you requested)
OUTDIR="/scratch/lt200308-agbsci/Quentin_project/00-samples/01-minimap2-hifi"
mkdir -p "$OUTDIR"

# -----------------------
# Threads (don’t actually use all 127 for minimap2)
# -----------------------
MM2_T=64
SM_T=32

BAM="${OUTDIR}/CG_F_HIFI_None_vs_ref.sorted.bam"
LOG="${OUTDIR}/CG_F_HIFI_None_vs_ref.minimap2.log"

echo "[INFO] REF=$REF"
echo "[INFO] READS=$READS"
echo "[INFO] OUTDIR=$OUTDIR"
echo "[INFO] minimap2 threads=$MM2_T ; samtools threads=$SM_T"

# -----------------------
# Index reference (safe if already indexed)
# -----------------------
if [[ ! -s "${REF}.fai" ]]; then
  samtools faidx "$REF"
fi

# -----------------------
# Map + sort + index
# -----------------------
minimap2 -ax map-hifi -t "$MM2_T" --secondary=no "$REF" "$READS" 2> "$LOG" \
  | samtools view -@ "$SM_T" -bh \
  | samtools sort -@ "$SM_T" -o "$BAM" -

samtools index -@ "$SM_T" "$BAM"

echo "[DONE] BAM=$BAM"
echo "[DONE] BAI=${BAM}.bai"
