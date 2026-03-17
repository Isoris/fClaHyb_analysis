# 04_pcangsd_ngsadmix_round1

This module contains the first-round structure and admixture inference workflow based on:
- PCAngsd
- NGSadmix
- evalAdmix

It is designed as the first integrated population-structure exploration block after:
- site discovery
- distance thinning
- BEAGLE genotype-likelihood generation

## Purpose

This module is used to:
- explore first-round global and local population structure
- estimate covariance / tree structure with PCAngsd
- estimate admixture proportions with NGSadmix
- evaluate model fit with evalAdmix
- choose the best seed per K
- generate stable sample-order and ancestry metadata for downstream analysis

## Folder structure

### 01_pcangsd
PCAngsd-based first-round structure analysis.

Main uses:
- covariance / PCA-like structure
- tree generation
- first-pass sample ordering
- early by-LG and whole-genome structure exploration

Important note:
- `--iter 250` is used in PCAngsd runs here because lower iteration settings did not always converge reliably
- in practice, PCAngsd outputs from this module are mainly used to help order samples and guide later NGSadmix interpretation

### 02_ngsadmix
NGSadmix-based ancestry inference.

Main uses:
- global ancestry inference across multiple thinning levels
- local-byRF ancestry inference
- K sweep and seed sweep
- first-round admixture model exploration

Global runs here use:
- thin200
- thin500
- thin1000
- K = 2..12
- seeds = 1, 2, 3

### 03_evaladmix
evalAdmix-based evaluation and best-seed selection.

Main uses:
- residual model-fit evaluation
- compare runs across K and seeds
- select the best seed per K
- build stable metadata tables for downstream plotting and interpretation

The first-round evalAdmix workflow here focuses on:
- thin500
- K = 2..12
- seeds = 1, 2, 3

## Overall run logic

### Step group 1: PCAngsd
1. Create by-LG BEAGLE list files
2. Run PCAngsd by linkage group
3. Run PCAngsd on whole-genome merged BEAGLE files

Main purpose:
- get covariance/tree structure
- guide sample ordering
- provide first-round structure context before NGSadmix

### Step group 2: NGSadmix
1. Run global NGSadmix across thin levels, K values, and seeds
2. Optionally run local-byRF NGSadmix

Main purpose:
- estimate ancestry components
- compare thin levels and seed behavior
- generate the main first-round admixture outputs

### Step group 3: evalAdmix
1. Run evalAdmix on the selected thin500 NGSadmix runs
2. Summarize residual strength across K and seeds
3. Select the best seed per K
4. Write final metadata tables and standardized best-run files

Main purpose:
- evaluate model fit
- stabilize run selection
- prepare clean downstream interpretation objects

## Dependencies

This module depends on outputs from:
- `../03_angsd_site_discovery/`

In particular, it requires:
- BEAGLE genotype-likelihood files
- thinned SNP panels
- sample order consistency across all downstream runs

## Important outputs

Most important outputs from this module include:
- PCAngsd tree/covariance outputs
- NGSadmix `.qopt`, `.fopt.gz`, and logs
- evalAdmix residual outputs
- best seed per K tables
- canonical sample order table
- cluster palette table
- per-sample ancestry metadata

## Notes

- This is a first-round structure workflow, not necessarily the final fully curated structure analysis.
- PCAngsd is mainly used here for covariance/tree structure and practical sample ordering.
- NGSadmix is the main first-round ancestry inference engine.
- evalAdmix is used to evaluate NGSadmix fit and help select the best seed per K.
- The outputs of this module are intended to feed later interpretation, figure-making, inversion analysis, and population-structure refinement steps.
