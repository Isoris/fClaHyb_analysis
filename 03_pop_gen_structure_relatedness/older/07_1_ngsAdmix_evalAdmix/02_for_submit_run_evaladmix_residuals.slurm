###
# 02_for_submit_run_evaladmix_residuals.slurm does this:
# Counts completed NGSadmix qopt files and submits the evalAdmix array job
# with the matching number of array tasks.
###
cd /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/popstruct_thin/05_ngsadmix_global
NRUNS=$(ls -1 runs_thin500/thin500_K*_seed*.qopt | wc -l)
sbatch --array=1-$NRUNS 02_run_evaladmix_residuals.slurm
