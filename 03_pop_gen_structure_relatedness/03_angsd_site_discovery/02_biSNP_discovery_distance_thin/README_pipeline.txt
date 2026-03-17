# 02_biSNP_discovery_distance_thin

This module performs candidate bi-allelic SNP discovery from ANGSD likelihood-based calls
using the global folded SFS mean prior generated in the previous module, then builds
distance-thinned SNP site panels for downstream population-genetic analyses.

## Purpose

This module starts after `01_saf_folded_sfs` is complete.

It uses:
- the global folded SFS mean prior (`catfish.global.folded.mean.pest`)
- the callable ANGSD sites mask
- the chromosome chunk list (`chunk_rf.list`)
- the filtered BAM list
- the reference FASTA

It produces:
- per-chromosome candidate biSNP calls (`.mafs.gz`)
- merged raw SNP site tables
- distance-thinned SNP panels
- ANGSD-indexed site files for downstream analyses
- summary statistics for manuscript writing

---

## Step order

### STEP01_run_angsd_biSNP_discovery_chunks.slurm
Runs ANGSD per chromosome/region (`-rf`) using the global folded SFS prior.

Input:
- chunk_rf.list
- reference FASTA
- BAM list
- callable ANGSD sites file
- global folded mean pest

Output:
- per-chromosome `.mafs.gz` files
- logs

Main result directory:
- `/scratch/.../popstruct_biSNP_discovery/01_snps/`

---

### STEP02_make_distance_thinned_sites_200_500_1000bp.slurm
Merges all per-chromosome `.mafs.gz` files, preserves major/minor alleles,
and builds distance-thinned SNP panels at:
- 200 bp
- 500 bp
- 1000 bp

Input:
- `/scratch/.../popstruct_biSNP_discovery/01_snps/*.mafs.gz`

Output:
- `sites.raw.ALL.majmin.tsv`
- `sites.thin_200.majmin.tsv`
- `sites.thin_500.majmin.tsv`
- `sites.thin_1000.majmin.tsv`
- ANGSD `.idx` and `.bin` files for each

Main result directory:
- `/scratch/.../popstruct_biSNP_discovery/02_thinned_sites/`

---

### STEP03_make_distance_thinned_sites_5_10_25kb.slurm
Builds larger-distance thinned SNP panels at:
- 5 kb
- 10 kb
- 25 kb

Input:
- merged raw SNP site table from STEP02

Output:
- `sites.thin_5000.majmin.tsv`
- `sites.thin_10000.majmin.tsv`
- `sites.thin_25000.majmin.tsv`
- ANGSD `.idx` and `.bin` files for each

Main result directory:
- `/scratch/.../popstruct_biSNP_discovery/02_thinned_sites/`

---

### STEP04_summarize_biSNP_stats_for_manuscript.sh
Summarizes raw and thinned SNP site sets for manuscript writing.

Reports:
- total raw SNP count
- per-chromosome SNP counts
- genome-wide SNP density
- per-chromosome SNP density
- thinned site counts
- allele composition
- Ti/Tv ratio
- draft copy-paste manuscript text

Input:
- `/scratch/.../popstruct_biSNP_discovery/02_thinned_sites/`

Output:
- console summary
- copy-paste text for methods/results

---

## Run order

1. Run STEP01
2. Run STEP02
3. Run STEP03
4. Run STEP04

---

## Important dependencies

This module depends on outputs from:
- `../01_saf_folded_sfs/`

Required upstream files include:
- global folded mean pest
- callable ANGSD mask
- chromosome chunk list

---

## Key final outputs

Most important reusable outputs:
- raw merged SNP table
- thinned SNP panels at 200 / 500 / 1000 bp
- thinned SNP panels at 5 / 10 / 25 kb
- indexed ANGSD site files for each panel

These files are used for:
- PCA / structure analyses
- chromosome-specific analyses
- distance-thinned downstream ANGSD runs
- manuscript summaries

---

## Notes

- All SNP sites are discovered from genotype likelihoods with ANGSD, not hard genotype calls.
- Distance thinning is applied separately within each chromosome.
- Major/minor allele columns are retained so downstream ANGSD runs can reuse consistent allele coding.
- This module is specifically for bi-allelic SNP site discovery and thinning, not final genotype calling.
