# 03_evaladmix

This module evaluates first-round NGSadmix runs using evalAdmix and selects the best seed per K.

## Purpose

This step is used to:
- evaluate residual structure after NGSadmix
- compare seeds across K values
- identify the best seed per K
- generate stable metadata for downstream plotting and interpretation

The final evaluation workflow here is focused on:
- thin500
- K = 2..12
- seeds = 1, 2, 3

## Main scripts

### STEP01_run_evaladmix_thin500_K2_12_3seeds.slurm
Runs evalAdmix on the selected thin500 global NGSadmix outputs.

Purpose:
- compute residual correlation structure
- evaluate model fit across K and seeds
- generate `.corres.txt` residual outputs for each run

Inputs:
- thin500 whole-genome BEAGLE file
- corresponding NGSadmix `.fopt.gz` and `.qopt` files

Typical outputs:
- `evaladmix_thin500/`
- one residual output per K/seed combination

### STEP02_summarize_evaladmix.R
Reads evalAdmix residual-like outputs and summarizes residual strength across K and seeds.

Purpose:
- compute a simple residual-strength metric
- write summary tables
- generate K-by-seed plots for quick interpretation

Typical outputs:
- `summary_residuals.tsv`
- residual-strength plots

### STEP03_select_best_seed_by_K.R
Combines NGSadmix and evalAdmix information to select the best seed per K.

Purpose:
- rank seeds within each K
- select the best run per K
- create standardized best-file copies/links
- write sample-order and ancestry metadata
- create stable cluster color palettes

Typical outputs:
- `all_seed_metrics_by_K.tsv`
- `best_seed_by_K.tsv`
- `best_seed_copied_files.tsv`
- `cluster_palette_by_K.tsv`
- `sample_order_reference.tsv`
- `sample_main_ancestry_by_K.tsv`

## Run order

1. Run STEP01 on the chosen thin500 global NGSadmix runs
2. Run STEP02 to summarize evalAdmix residual strength
3. Run STEP03 to select the best seed per K and write metadata tables

## Inputs

This module depends on outputs from:
- `../02_ngsadmix/`

Main required inputs:
- thin500 NGSadmix run directory
- thin500 whole-genome BEAGLE file
- sample file or canonical sample manifest

## Important outputs

Most important final outputs:
- best seed per K
- standardized best `.qopt`, `.fopt.gz`, `.log`, and evalAdmix outputs
- palette table
- canonical sample order table
- per-sample dominant ancestry table

## Notes

- This module is intentionally focused on thin500 for first-round model evaluation.
- evalAdmix is used here as a model-fit check, not as the primary ancestry inference method.
- The best-seed selection step is the most downstream and most metadata-rich part of the first-round structure workflow.
