# 03_beagle_GL

This module generates BEAGLE genotype-likelihood (GL) files from the distance-thinned
biallelic SNP panels produced in the previous module. It supports both chromosome-wise
(per-RF) BEAGLE generation and direct whole-genome BEAGLE generation, and can merge
per-RF BEAGLE files into whole-genome BEAGLE files for downstream analyses.

## Purpose

This module starts after `02_biSNP_discovery_distance_thin` is complete.

It uses:
- the thinned SNP site panels with preserved major/minor alleles
- the reference FASTA
- the BAM list
- the chromosome chunk list (`chunk_rf.list`) for per-RF runs

It produces:
- per-RF BEAGLE genotype-likelihood files
- direct whole-genome BEAGLE genotype-likelihood files
- merged whole-genome BEAGLE files created from per-RF outputs
- per-thinning list files documenting which BEAGLE chunks were merged

## Main scripts

### STEP01_run_beagle_GL_byRF_chunks.slurm
Runs ANGSD per chromosome/region (`-rf`) to generate BEAGLE genotype-likelihood files
from a given thinned SNP panel.

Inputs:
- `chunk_rf.list`
- reference FASTA
- BAM list
- thinned SNP panel with major/minor columns:
  - `sites.thin_<W>.majmin.tsv`
  - plus its `.idx` and `.bin` files

Typical usage:
- submit as a SLURM array
- one run per thinning level (`W`)
- one array task per chromosome/region RF file

Typical outputs:
- per-RF `.beagle.gz` files under thinning-specific directories
- logs

Purpose:
- allows BEAGLE GL generation in smaller chromosome-wise jobs
- useful when whole-genome jobs are too large or when chromosome-level provenance is desired

### STEP02_run_beagle_GL_wholegenome.slurm
Runs ANGSD on the whole genome to generate a single BEAGLE genotype-likelihood file
for a selected thinning level.

Inputs:
- reference FASTA
- BAM list
- thinned SNP panel with major/minor columns:
  - `sites.thin_<W>.majmin.tsv`
  - plus its `.idx` and `.bin` files

Typical outputs:
- one whole-genome `.beagle.gz` file per thinning level
- logs

Purpose:
- direct whole-genome BEAGLE generation
- useful as a simpler single-step alternative to chromosome-wise generation + merge

### STEP03_merge_beagle_byRF_to_wholegenome.sh
Merges per-RF BEAGLE files into one whole-genome BEAGLE file for each requested
thinning level.

This step:
- finds all per-RF `.beagle.gz` files for a thinning level
- writes a persistent `.list` file for provenance
- keeps the BEAGLE header only once
- appends all data rows in sorted RF order
- compresses the merged output
- validates the gzip output

Inputs:
- per-RF `.beagle.gz` files from STEP01

Outputs:
- `beagle_thin_<W>.list`
- `catfish.wholegenome.byRF.thin_<W>.beagle.gz`

Purpose:
- creates whole-genome BEAGLE files from chromosome-wise results
- preserves provenance of which chromosome files were included
- useful when per-RF generation is preferred for compute reasons

## Run order

### Minimal route A: direct whole-genome BEAGLE
1. Run STEP02 for a chosen thinning level

### Minimal route B: chromosome-wise BEAGLE then merge
1. Run STEP01 for a chosen thinning level
2. Run STEP03 for that thinning level

### Typical multi-thinning workflow
1. Run STEP01 for each required thinning level
2. Run STEP03 to merge per-RF BEAGLE files across all required thinning levels
3. Optionally run STEP02 as a direct whole-genome comparison or backup output

## Dependencies

This module depends on outputs from:
- `../02_biSNP_discovery_distance_thin/`

Required upstream files include:
- thinned SNP site panels with preserved major/minor alleles
- ANGSD site indices (`.idx` and `.bin`)
- chromosome chunk list for per-RF runs

## Important inputs

Typical upstream thinned SNP files:
- `sites.thin_200.majmin.tsv`
- `sites.thin_500.majmin.tsv`
- `sites.thin_1000.majmin.tsv`
- `sites.thin_5000.majmin.tsv`
- `sites.thin_10000.majmin.tsv`
- `sites.thin_25000.majmin.tsv`

Typical reference-related inputs:
- reference FASTA
- BAM list
- chromosome chunk list (`chunk_rf.list`)

## Important outputs

Reusable outputs from this module include:
- per-RF BEAGLE GL files for each thinning level
- direct whole-genome BEAGLE GL files
- merged whole-genome BEAGLE GL files
- `.list` files documenting the per-RF files included in each merge

Most important final outputs:
- whole-genome `.beagle.gz` files for downstream population-structure analyses

## Notes

- This module uses ANGSD genotype likelihoods, not hard genotype calls.
- Major/minor allele coding is preserved from the upstream thinned site panels.
- Per-RF BEAGLE generation and direct whole-genome BEAGLE generation are both valid routes.
- Merged per-RF whole-genome BEAGLE files and directly generated whole-genome BEAGLE files may both be useful for QC or workflow comparison.
- The `.list` files written during merging are intentionally retained for provenance and reuse.
- Consistent sample order across all BEAGLE files is required.
- Consistent thinning-panel naming across modules is important.

## Final goal of this module

Produce BEAGLE genotype-likelihood files at multiple thinning levels for downstream
population-structure, PCA, admixture, relatedness, and other ANGSD-compatible analyses.
