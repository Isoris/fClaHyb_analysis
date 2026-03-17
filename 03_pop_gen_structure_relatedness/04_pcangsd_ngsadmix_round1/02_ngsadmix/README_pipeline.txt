# 02_ngsadmix

This module runs first-round NGSadmix ancestry inference from BEAGLE genotype-likelihood files.

## Purpose

This step is used to:
- estimate admixture proportions across K values
- compare multiple thinning levels
- evaluate seed stability
- generate global and local-byRF ancestry solutions for downstream comparison

The main global workflow uses:
- thinning levels: 200 / 500 / 1000 bp
- K = 2..12
- seeds = 1, 2, 3

## Main scripts

### STEP01_run_ngsadmix_global_thinKsweep_3seeds.slurm
Runs global NGSadmix on whole-genome BEAGLE files across:
- multiple thinning levels
- K = 2..12
- 3 seeds

Purpose:
- first-round global ancestry inference
- compare structure signal across thinning levels
- generate the main global NGSadmix run directories

Typical outputs:
- `runs_thin200/`
- `runs_thin500/`
- `runs_thin1000/`

### HELPER_submit_ngsadmix_global_K2_12_3seeds.slurm
Helper submit script for the global NGSadmix array job.

Purpose:
- creates logs directory if needed
- submits the full global K × seed array

### STEP02_run_ngsadmix_local_byRF_K2_12_3seeds.slurm
Runs local NGSadmix by RF / linkage-group BEAGLE files.

Purpose:
- local ancestry exploration by chromosome / RF
- compare local structure patterns across LGs
- evaluate whether local ancestry signals differ from global patterns

Inputs:
- one list file containing per-RF BEAGLE file paths

Typical outputs:
- `runs_LGXX_thinYYY/` style folders

## Run order

### Global route
1. Submit global NGSadmix with the helper script or directly with STEP01

### Local-byRF route
1. Prepare the per-RF BEAGLE list file
2. Run STEP02 as an array across listed BEAGLE files

## Inputs

This module depends on outputs from:
- `../01_pcangsd/` for sample ordering/context
- `../../03_angsd_site_discovery/03_beagle_GL/` for BEAGLE files

Main required inputs:
- whole-genome BEAGLE files
- local-byRF BEAGLE files
- sample ordering consistency across all runs

## Important outputs

Most important outputs:
- global NGSadmix run folders across thin levels
- local-byRF NGSadmix run folders
- `.qopt`
- `.fopt.gz`
- `.log`

## Notes

- This is the first-round NGSadmix module.
- The global runs are the main ancestry inference outputs.
- The local-byRF runs are exploratory and help assess chromosome-level ancestry heterogeneity.
- Seeds are fixed to 1, 2, and 3 for reproducibility and basic stability checks.
