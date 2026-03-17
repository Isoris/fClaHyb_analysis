# 05_relatedness_pruning

This module handles pairwise relatedness estimation, relative filtering, and relatedness plotting.
It is inspired from 'https://www.science.org/doi/10.1126/sciadv.adm7980 '
## Steps

### STEP01
Current script:
- `08-ngsRelate_onthin500_and_prune_relatives.slurm`

Purpose:
- run ngsRelate on the thin500 whole-genome BEAGLE panel
- generate the main pairwise relatedness result table

Main output:
- `catfish_226_relatedness.res`

### STEP02
Current script:
- `08b_run_natura_cull_relatives_v2.slurm`

Helper:
- `NAToRA_Public.py`

Purpose:
- convert ngsRelate output to NAToRA input
- run NAToRA across multiple theta cutoffs
- summarize removed vs retained individuals

### STEP03
Current script:
- `pairwise_first_degree_curation.py`

Purpose:
- alternate pruning route based on pairwise theta threshold
- generate keep/remove lists with explicit decision rules

### STEP04
Current script:
- `plot_relatedness_3panel_OLD_FROM_HPC.R`

Purpose:
- plot relatedness heatmaps, pair counts, and summary figures

## Run order

1. Run STEP01
2. Run STEP02
3. Optionally run STEP03 as an alternate filtering strategy
4. Run STEP04 for visualization

## Notes

- STEP02 is the main filtering route in this workflow.
- STEP03 is an alternate explicit pairwise pruning method.
- thin500 is the relatedness panel used here.
