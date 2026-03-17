# 01_pcangsd

This module runs first-round PCAngsd analyses from BEAGLE genotype-likelihood files.

## Purpose

This step is mainly used to:
- estimate covariance / PCA-like structure
- generate tree outputs
- get a reasonable sample ordering for later NGSadmix and evalAdmix steps

The PCAngsd runs here use `--iter 250` so admixture optimization has more time to converge.

## Main scripts

### STEP01_make_pcangsd_byLG_lists_200_500_1000.sh
Creates per-linkage-group BEAGLE list files for:
- 200 bp
- 500 bp
- 1000 bp

Outputs:
- `beagle_LG_thin_200.list`
- `beagle_LG_thin_500.list`
- `beagle_LG_thin_1000.list`

### STEP02_make_pcangsd_byLG_lists_5000_10000_25000.sh
Creates per-linkage-group BEAGLE list files for:
- 5000 bp
- 10000 bp
- 25000 bp

Outputs:
- `beagle_LG_thin_5000.list`
- `beagle_LG_thin_10000.list`
- `beagle_LG_thin_25000.list`

### STEP03_run_pcangsd_byLG_thin200_500_1000_K2_12_iter250.slurm
Runs PCAngsd per linkage group using BEAGLE files from the 200/500/1000 bp thinned panels.

Main settings:
- `K = 2..12`
- `--iter 250`
- admixture seed fixed to `1`
- tree output enabled

Purpose:
- local / per-LG structure exploration
- local covariance/tree patterns
- first-round admixture signal by chromosome

### STEP04_run_pcangsd_global_WG_thin200_500_1000_K2_12_iter250.slurm
Runs PCAngsd on whole-genome merged BEAGLE files for:
- 200 bp
- 500 bp
- 1000 bp

Main settings:
- `K = 2..12`
- `--iter 250`
- tree output enabled

Purpose:
- global covariance / tree inference
- sample ordering for downstream NGSadmix
- first-round whole-genome structure exploration

## Run order

1. Run STEP01 to create by-LG lists for 200/500/1000
2. Run STEP02 if 5000/10000/25000 by-LG lists are also needed
3. Run STEP03 for per-LG PCAngsd
4. Run STEP04 for whole-genome PCAngsd

## Inputs

This module depends on outputs from:
- `../../03_angsd_site_discovery/03_beagle_GL/`

Main required inputs:
- BEAGLE genotype-likelihood files
- per-LG BEAGLE lists
- sample name file in the same order as the BAM/BEAGLE input

## Important outputs

Typical outputs include:
- PCAngsd covariance-related results
- tree files
- admixture outputs across K values
- organized by thinning level and linkage group or whole genome

## Notes

- `--iter 250` is used here because lower iteration settings did not always converge reliably.
- This module is mainly a first-round exploration step.
- In practice, the global PCAngsd tree/covariance outputs are especially useful for ordering samples before running NGSadmix in the next module.
