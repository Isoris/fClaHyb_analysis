# 07_polymorphic_inversions
# Inversion Local PCA Workflow from ANGSD Genotype Likelihoods

This workflow is an adaptation of the local PCA / lostruct framework of **Li and Ralph**, recoded for our data and our pipeline. Instead of going through hard genotype calls in VCF format, this workflow starts from **ANGSD genotype likelihoods (GLs)** and converts them into a **BEAGLE-based dosage representation** for downstream local PCA, MDS, candidate inversion region discovery, and follow-up validation plots.

The goal is to keep the logic of local PCA in genomic windows, while using an input format that is natural for our ANGSD workflow and avoids unnecessary hard-calling.

---

## Main idea

The workflow has six main parts:

1. Build a callable inversion mask from the reference FASTA.
2. Compute a folded SFS prior under the inversion-specific filtering scheme.
3. Use that prior to call candidate biallelic SNPs and export BEAGLE genotype likelihoods.
4. Convert BEAGLE to dosage matrices, run local PCA in windows, compute window-to-window distances, perform MDS, and define candidate outlier regions.
5. Overlap candidate regions with theta / heterozygosity windows from another workflow.
6. For a chosen candidate region, run a regional PCA, assign genotype-like groups, and make panel A/B/C validation figures.

---

## Why this workflow exists

The original Li and Ralph local PCA framework commonly reads genotypes from VCF-like inputs and performs PCA in local windows. Here, we adapt that idea to ANGSD by:

- using **GL model 1**
- estimating the **SFS prior** under the same filters and callable mask used for the inversion scan
- using **`-pest`** so SNP calling is informed by the correct folded SFS prior
- exporting **BEAGLE** genotype likelihoods instead of creating a fake VCF
- retaining the two alleles per SNP from the BEAGLE file, so we preserve exact allele identities (`A/C`, `G/T`, etc.)
- converting genotype likelihoods to expected dosage for local PCA

This keeps the workflow consistent with ANGSD and avoids unnecessary hard genotype conversion.

---

## Why we recompute SAF / SFS / pest

The inversion workflow uses a different callable mask and a different filtering strategy from the global population-structure panel. Because of that, we recompute:

- SAF
- folded SFS
- mean pest prior

under the **inversion-specific settings**.

This is important because the prior used with `-pest` should match the exact filtering and callable-space definition of the downstream SNP-calling step.

---

## Why we use `-doSaf 5` and `-doMajorMinor 2`

In the SAF step we use:

- `-doSaf 5`
- `-doMajorMinor 2`

This is done so that ANGSD estimates the site allele frequency likelihoods and the folded SFS prior in a way that is appropriate for our GL-based workflow and does not depend on an already hard-called SNP panel.

Conceptually:

- `-doSaf 5` gives the SAF values needed to estimate the folded SFS
- `-doMajorMinor 2` lets ANGSD infer alleles appropriately for the SAF/SFS step

Then `realSFS` is used to obtain the folded SFS, and bootstrap rows are averaged to obtain the mean `.pest` prior used later.

---

## Why we do not use VCF

We do **not** need to create a VCF for the inversion discovery step.

Instead:

1. ANGSD calls candidate SNPs with:
   - the inversion-specific callable mask
   - inversion-specific filters
   - the inversion-specific folded SFS prior via `-pest`
2. ANGSD writes a **BEAGLE** file with genotype likelihoods
3. The BEAGLE file already contains:
   - marker
   - allele1
   - allele2
   - genotype likelihoods for each sample

So we already have the two alleles for each SNP and the genotype uncertainty per individual. That means we can stay in ANGSD space and do not need an intermediate fake VCF.

---

## Callable mask strategy

Two masks are conceptually useful:

- **strict uppercase `normalACGT` mask** for conservative global population structure
- **broader `inversion_acgt_allcase` mask** for inversion discovery

For inversion discovery, we allow:

- uppercase `A/C/G/T`
- lowercase `a/c/g/t` (soft-masked sequence)

but exclude:

- `N/n`
- other ambiguous symbols

This broader inversion mask preserves local marker density while still excluding hard-masked sequence and non-ACGT symbols.

---

## Core inversion filters used here

Current inversion settings:

- `-GL 1`
- `-minQ 20`
- `-minMapQ 60`
- `-minInd 203`
- `-setMaxDepth 2850`
- `-minMaf 0.05`
- `-SNP_pval 1e-6`

And:

- `-remove_bads 1`
- `-uniqueOnly 1`
- `-only_proper_pairs 1`

These are intended for the inversion-discovery panel, not the strict global structure panel.

---

## Why BEAGLE is converted to dosage

For each biallelic SNP in the BEAGLE file, we have:

- allele1
- allele2
- `P(AA)`
- `P(AB)`
- `P(BB)`

For local PCA we convert this to **expected dosage** of allele2:

`dosage = P(AB) + 2 * P(BB)`

This gives a numeric matrix suitable for covariance estimation and PCA, while the exact alleles are still preserved in the accompanying site metadata table.

So:

- the **sites file** keeps the exact alleles and genomic coordinates
- the **dosage matrix** is used for PCA math

---

## Route 1 vs Route 2 logic

There are two ways one could think about local PCA from genotype likelihoods.

### Route 1: collapse GLs to dosage, then do PCA
This is the route implemented here.

At each SNP and for each individual, the three genotype likelihoods are collapsed to one expected dosage value:

`dosage = P(AB) + 2 * P(BB)`

Then local covariance and PCA are computed from that numeric matrix.

Why use this route:

- simple
- transparent
- easy to debug
- easy to preserve exact allele identities in metadata
- easy to connect back to SNP-level interpretation

Limitation:

- a highly uncertain genotype and a confident intermediate genotype can have similar dosage values

### Route 2: likelihood-aware covariance estimation
This is closer to what PCAngsd is designed to do.

Instead of first collapsing each genotype to one expected dosage value, the covariance matrix is estimated more directly from genotype-likelihood information.

Why this route exists:

- more faithful to genotype uncertainty
- useful for low-depth data with variable confidence among sites and individuals

Why it is **not** used here:

- more complex to inspect and rewrite
- less transparent for a first custom implementation
- the dosage route is a reasonable first approximation and easier to integrate with this recoded Li-and-Ralph-style workflow

So this pipeline uses:

**Route 1 = GLs → dosage → covariance → PCA**

---

## What the PCA points mean

In the regional PCA plots, the **points are individuals**.

For one SNP, an individual can be:

- `AA`
- `AB`
- `BB`

Across many SNPs in a candidate inversion region, individuals may separate into the classic three inversion-like groups:

- one homokaryotype-like class
- one heterokaryotype-like class
- the other homokaryotype-like class

The PCA itself is done on the individual covariance structure across many SNPs in the region/window.

---

## Why heterozygosity plot C matters

The MDS / outlier windows identify candidate genomic regions with deviant local population structure.

But a panel C–style heterozygosity plot provides a biological validation layer:

- individuals are grouped by the **regional PCA** in the candidate region
- chromosome-wide heterozygosity curves are then plotted for those groups
- if the groups differ most strongly in the candidate region, that supports an inversion-like interpretation

Important:

- heterozygosity is **allele-polarity invariant**
- this means plot C is safe even if allele coding flips across SNPs

---

## Allele polarity and optional SNP flipping

If a candidate region is later summarized by regional haplotype score, allele coding direction can matter.

For example:

- one SNP may be oriented so that allele2 corresponds to one inversion background
- another SNP may be oriented the opposite way

This does **not** affect heterozygosity, but it can affect regional dosage summaries.

So in the candidate-region follow-up PCA step, SNPs can optionally be oriented to the regional PC1 axis:

- if a SNP is negatively associated with regional PC1, flip its dosage:
  - `dosage_flipped = 2 - dosage`

This creates a more consistent regional haplotype score.

This is **optional** and used for interpretation; it is not required for plot C.

---

## Workflow steps

### STEP01 — `STEP01_mask_regions_from_fasta.py`

Reads the FASTA and writes BED tracks for:

- hard-masked regions
- soft-masked regions
- normal uppercase ACGT
- masked soft+hard regions
- inversion-callable all-case ACGT regions

Also writes a per-contig stats table.

### STEP02 — `STEP02_make_inversion_callable_angsd_mask.sh`

Converts the inversion BED file to ANGSD 1-based sites format and indexes it with:

`angsd sites index`

### STEP03 — `STEP03_mask_depth_mapq_stats.py`

Computes manuscript-ready mask statistics:

- total bp
- percent of genome
- number of intervals
- mean interval size
- mean depth in the mask
- mean depth after MAPQ filtering

This can be used to report coverage characteristics of the strict and inversion masks.

### STEP04 — `STEP04_run_angsd_saf_chunks_inversion.slurm`

Runs ANGSD per RF chunk to compute SAF under inversion-specific filters using:

- `-GL 1`
- `-doSaf 5`
- `-doMajorMinor 2`

and then computes a folded bootstrap SFS for each chunk.

### STEP05 — `STEP05_merge_chunk_saf_to_global_inversion.slurm`

Merges chunk SAF files into a single global SAF using `realSFS cat`.

### STEP06 — `STEP06_make_global_folded_sfs_and_mean_pest_inversion.slurm`

Runs `realSFS` on the merged SAF to compute a folded bootstrap SFS and averages the bootstrap rows into a mean `.pest` file.

This `.pest` file is the inversion-specific prior used in the next step.

### STEP07 — `STEP07_angsd_call_snps_and_beagle_inversion.slurm`

Calls candidate biallelic SNPs using:

- the inversion-specific callable mask
- the inversion-specific folded SFS prior via `-pest`
- `-SNP_pval`
- `-minMaf`

and writes:

- `.mafs.gz`
- `.beagle.gz`

This is the key SNP-discovery step for the inversion pipeline.

### STEP08 — `STEP08_beagle_to_dosage_by_chr.py`

Reads the `.beagle.gz` file and converts genotype likelihoods to expected dosage per sample:

`dosage = P(AB) + 2 * P(BB)`

Writes per-chromosome:

- `CHR.sites.tsv.gz`
- `CHR.dosage.tsv.gz`

The sites file preserves exact allele identities and coordinates.

### STEP09 — `STEP09_local_pca_windows_by_chr.R`

Reads one chromosome dosage file and its matching sites file, splits the chromosome into fixed SNP windows (default 100 SNPs), and computes local PCA summaries per window.

This is the local PCA step.

### STEP10 — `STEP10_window_mds_outliers.R`

Takes all STEP09 window PCA summaries, computes lostruct-style distances among windows, performs MDS, and identifies outlier windows.

Nearby outlier windows are clustered into candidate regions.

This is the main window-MDS and candidate-region discovery step.

### STEP11 — `STEP11_overlap_candidate_regions_with_theta_and_het.R`

Overlaps candidate local-PCA outlier regions with external theta / heterozygosity windows from another workflow using genomic coordinates.

This step links candidate inversion regions to independent diversity statistics.

### STEP12 — `STEP12_candidate_region_pca_groups_and_plotC.R`

For one chosen candidate region:

- extract SNPs from the regional dosage matrix
- run a regional PCA on individuals
- assign 3 genotype-like groups using PC1
- compute regional heterozygosity per individual
- optionally orient SNP dosages to the regional PC1 axis
- merge group labels with chromosome-wide heterozygosity windows
- generate:
  - panel A–style regional PCA
  - panel B–style PCA colored by regional heterozygosity
  - panel C–style chromosome heterozygosity curves by inferred group

### STEP13 — `STEP13_make_combined_panels_ABC.R`

Combines the outputs from STEP12 into a single manuscript-style multi-panel figure:

- **(a)** plain regional PCA
- **(b)** regional PCA colored by regional heterozygosity, with representative individuals highlighted
- **(c)** chromosome-wide heterozygosity curves by inferred group, with candidate region shaded

This writes a combined PDF and PNG.

---

## Important conceptual note

The local PCA windows and theta/heterozygosity windows do **not** need to be identical.

For example:

- local PCA may use **100-SNP windows**
- theta may use **10 kb windows**

These are matched later by genomic coordinate overlap, not by assuming they are the same window set.

---

## Minimal run order

```bash
STEP01 -> STEP02 -> STEP04 -> STEP05 -> STEP06 -> STEP07 -> STEP08 -> STEP09 -> STEP10 -> STEP11 -> STEP12 -> STEP13

Here is the fully corrected block, with the nested code fences escaped properly for markdown.

````markdown id="zbg3xo"
## Example usage outline

### 1. Build mask

```bash
bash STEP02_make_inversion_callable_angsd_mask.sh \
  --fasta /path/to/ref.fa \
  --mask-script STEP01_mask_regions_from_fasta.py \
  --prefix /path/to/ref.mask_regions
````

### 2. Run SAF per chunk

```bash
sbatch --array=0-$((N-1))%8 STEP04_run_angsd_saf_chunks_inversion.slurm chunk_rf.list
```

### 3. Merge SAF

```bash
sbatch STEP05_merge_chunk_saf_to_global_inversion.slurm
```

### 4. Make folded SFS and pest

```bash
sbatch STEP06_make_global_folded_sfs_and_mean_pest_inversion.slurm
```

### 5. Call SNPs and write BEAGLE

```bash
sbatch --array=0-$((N-1))%8 STEP07_angsd_call_snps_and_beagle_inversion.slurm chunk_rf.list
```

### 6. Convert BEAGLE to dosage by chromosome

```bash
python3 STEP08_beagle_to_dosage_by_chr.py \
  --beagle /path/to/catfish.LG01.beagle.gz \
  --outdir /path/to/dosage_by_chr
```

### 7. Run local PCA for one chromosome

```bash
Rscript STEP09_local_pca_windows_by_chr.R \
  LG01.sites.tsv.gz \
  LG01.dosage.tsv.gz \
  STEP09_LG01 \
  2 \
  100
```

### 8. Run MDS / outlier windows across all chromosomes

```bash
Rscript STEP10_window_mds_outliers.R \
  /path/to/step09_outputs \
  inversion_localpca \
  2 \
  20 \
  3 \
  1000000 \
  3
```

### 9. Overlap candidates with theta / heterozygosity windows

```bash
Rscript STEP11_overlap_candidate_regions_with_theta_and_het.R \
  inversion_localpca.candidate_regions.tsv.gz \
  theta_het_windows.tsv.gz \
  inversion_localpca \
  candidate_id
```

### 10. Build regional PCA / plot C for one candidate

```bash
Rscript STEP12_candidate_region_pca_groups_and_plotC.R \
  inversion_localpca.candidate_regions.tsv.gz \
  3 \
  LG01.sites.tsv.gz \
  LG01.dosage.tsv.gz \
  sample_het_windows.tsv.gz \
  inversion_localpca \
  sample \
  Hobs
```

### 11. Combine panels A/B/C

```bash
Rscript STEP13_make_combined_panels_ABC.R \
  inversion_localpca \
  3 \
  inversion_localpca
```

## Output logic

Main key outputs:

* inversion callable ANGSD sites file
* inversion global folded `.pest`
* per-chunk `.beagle.gz`
* per-chromosome dosage matrices
* per-window local PCA summaries
* MDS coordinates for windows
* candidate outlier regions
* overlap summaries with theta / heterozygosity
* regional PCA sample groups
* panel A / B / C plots
* combined A / B / C figure

## Summary

This is a GL-based local PCA inversion workflow adapted from Li and Ralph’s local PCA / lostruct logic, but recoded to work directly from ANGSD genotype likelihood outputs.

The main design choices are:

* keep the workflow in ANGSD space
* recompute the SFS prior under the inversion-specific mask and filters
* use BEAGLE output instead of VCF
* preserve exact allele identities
* convert GLs to dosage for local PCA
* use MDS on window PCA summaries to identify candidate inversion regions
* connect candidate regions to theta / heterozygosity in a later overlap step
* validate selected candidates with regional PCA and chromosome-wide heterozygosity curves

```
```
