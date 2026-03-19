# Catfish Heterozygosity + ROH/FROH Workflow

## Overview

Complete modular workflow for:
1. **Per-sample genome-wide heterozygosity** (ANGSD GL1 + doSaf1 + realSFS)
2. **Local diversity tracks** (theta proxy via -pest + doThetas + thetaStat)
3. **ROH / FROH** (ngsF-HMM multi-replicate → convert_ibd.pl → summarize_ibd_roh.py)
4. **Het inside vs outside ROH** (theta proxy intersected with ROH intervals)
5. **Rich plots and tables** with optional ancestry/relatedness metadata overlays
6. **Manuscript-ready report** (auto-filled markdown template)

## Important Interpretation Notes

- **Local theta tracks** are diversity proxies (NOT literal per-site Hobs)
- **ROH interval length** is a tract span, not direct homozygosity evidence for every base
- **Repeats and masked regions** are uninformative, not homozygous
- **No `-doHWE`** is used — this hatchery population has family structure / Wahlund effect
- All analysis runs on **all QC-passing samples** (not pruned 81)
- Pruned-81 is used only for optional cleaner visualization

## Directory Structure

```
het_roh/
├── 01_inputs_check/     # Validated inputs, BAM list, regions file
├── 02_heterozygosity/   # Per-sample SAF, SFS, theta
│   ├── 01_saf/
│   ├── 02_sfs/
│   ├── 03_theta/
│   └── 04_summary/
├── 03_ngsF_HMM/         # Multi-replicate ngsF-HMM runs
│   ├── replicates/
│   └── best/
├── 04_roh_summary/      # ROH BED tracts, summaries
├── 05_inversion_support/ # (future: local PCA + kmeans)
├── 06_plots_core/       # Core genome-wide plots
├── 07_plots_metadata/   # Ancestry/family overlay plots
├── 08_stats/            # Statistical test results
├── 09_final_tables/     # All clean TSV output tables
├── 10_report/           # Markdown methods/results/limitations
└── logs/                # Per-step log files
```

## Quick Start

```bash
# 1. Edit 00_config.sh paths to match your setup
vim scripts/00_config.sh

# 2. Run everything sequentially
bash run_pipeline.sh

# OR run step-by-step:
bash run_pipeline.sh --step 1       # validate inputs
bash run_pipeline.sh --step 2 --slurm  # submit SLURM array
# ... wait for SLURM ...
bash run_pipeline.sh --from 3       # continue from step 3
```

## Pipeline Steps

| Step | Script | What it does |
|------|--------|-------------|
| 1 | `01_prep_inputs.sh` | Validate BAMs, ref, callable BED; extract sample list from manifest |
| 2 | `02_run_heterozygosity.sh` | Per-sample ANGSD het + local theta (adapted from existing script) |
| 2s | `02_run_heterozygosity_slurm.sh` | SLURM array wrapper for step 2 |
| 3 | `03_run_ngsF_HMM.sh` | ngsF-HMM multi-start, select best replicate |
| 4 | `04_parse_roh_and_het.sh` | Convert .ibd→BED, compute FROH, het in/out ROH, master table |
| 5 | `05_run_all_plots.sh` | All plots + statistics + report generation |

## Output Tables

| Table | Description |
|-------|-------------|
| `genomewide_heterozygosity.tsv` | Per-sample het from 1D-SFS |
| `per_sample_roh.tsv` | ROH total bp, n_tracts, longest, mean/median length, FROH |
| `per_chr_roh_summary.tsv` | Per-sample per-chromosome ROH/FROH |
| `per_sample_het_in_out_roh.tsv` | Theta proxy inside vs outside ROH |
| `roh_tracts_all.bed` | All ROH tracts (chr, start, end, sample, length) |
| `catfish_roh.per_sample_roh.tsv` | Extended ROH summary with ancestry/kinship annotations |
| `catfish_roh.per_sample_roh_bins_long.tsv` | ROH size-class bins (plot-ready long format) |
| `master_summary.tsv` | Merged het + ROH + FROH + theta in/out ROH |

## Output Plots

### Core (06_plots_core/)
- Genome-wide het: boxplot, violin, histogram, ranked bar
- FROH: boxplot, violin
- Longest ROH, ROH count, total ROH: boxplots
- ROH bins: stacked bar, distribution
- Scatterplots: het vs FROH, het vs ROH, depth vs het/FROH, theta in vs out ROH
- Per-chromosome: FROH boxplot, ROH burden, sample×chr heatmaps, ROH frequency
- Theta ideograms: per-sample genome-wide, per-chromosome, multi-sample mean

### Metadata overlays (07_plots_metadata/)
- All core plots grouped/colored by ancestry
- Pruned-81 subset heatmaps and boxplots
- ROH bins stacked bars faceted by ancestry

## Statistical Tests (08_stats/)

- `spearman_correlations.tsv` — het vs FROH, het vs ROH, het vs longest ROH
- `group_comparisons.tsv` — Wilcoxon / Kruskal-Wallis by ancestry (if available)
- `descriptive_summary.tsv` — mean, median, SD, range for all main variables

## Dependencies

- ANGSD (angsd, realSFS, thetaStat)
- ngsF-HMM
- bedtools
- perl (for convert_ibd.pl)
- python3
- R with: data.table, ggplot2

## Key Config Settings (00_config.sh)

```bash
MINQ=20           # Base quality
MINMAPQ=30        # Mapping quality
CLIP=50           # Clip overlapping reads
REALSFS_MAXITER=2000
REALSFS_TOLE=1e-16
WIN=500000         # thetaStat window
STEP=500000        # thetaStat step
NGSFHMM_REPS=10    # Random-start replicates
```

## Existing Scripts Reused

- `per_sample_het_angsdsaf1.sh` → adapted as `02_run_heterozygosity.sh`
- `convert_ibd.pl` → called by `04_parse_roh_and_het.py`
- `summarize_ibd_roh.py` → called by `04_parse_roh_and_het.py`
- Delly config pattern → adapted as `00_config.sh`
