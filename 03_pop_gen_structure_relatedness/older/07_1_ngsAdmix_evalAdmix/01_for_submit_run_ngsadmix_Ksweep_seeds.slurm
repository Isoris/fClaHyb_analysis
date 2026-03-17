###
# 01_for_submit_run_ngsadmix_Ksweep_seeds.slurm does this:
# Creates the logs directory for global NGSadmix runs and submits the
# K-sweep x seed array job for the predefined whole-genome run.
###
mkdir -p /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/popstruct_thin/05_ngsadmix_global/logs
cd /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/popstruct_thin/05_ngsadmix_global

# total = (K_MAX-K_MIN+1) * 3 seeds = 11*3 = 33 if K=2..12
sbatch --array=1-33 01_run_ngsadmix_Ksweep_seeds.slurm
