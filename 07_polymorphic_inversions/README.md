# 07_polymorphic_inversions
# Inversion Local PCA Workflow from ANGSD Genotype Likelihoods

This workflow is an adaptation of the local PCA / lostruct framework of **Li and Ralph**, recoded for our data and our pipeline. Instead of going through hard genotype calls in VCF format, this workflow starts from **ANGSD genotype likelihoods (GLs)** and converts them into a **BEAGLE-based dosage representation** for downstream local PCA, MDS, candidate inversion region discovery, and follow-up validation plots.

---

## Pipeline steps

```
STEP01 → STEP02 → STEP03 (optional QC)
                ↓
STEP04 → STEP05 → STEP06 → STEP07 → STEP08
                                        ↓
STEP09 → STEP10 → STEP11 → STEP12 → STEP13
```

### Mask building
- **STEP01** — Extract BED tracks from reference FASTA (hard-masked, soft-masked, inversion-callable)
- **STEP02** — Convert inversion BED to ANGSD 1-based sites and index
- **STEP03** — Compute mask depth/MAPQ stats for manuscript QC tables

### SFS prior estimation
- **STEP04** — Run ANGSD SAF per RF chunk under inversion-specific filters
- **STEP05** — Merge chunk SAFs into a global SAF (`realSFS cat`)
- **STEP06** — Estimate folded SFS, average bootstrap rows → mean `.pest` prior

### SNP calling and dosage conversion
- **STEP07** — Call biallelic SNPs with `-pest`, export BEAGLE GLs
- **STEP08** — Convert BEAGLE GLs to per-chromosome dosage matrices

### Local PCA and candidate discovery
- **STEP09** — Local PCA in SNP windows per chromosome
- **STEP10** — Window-to-window distances, MDS, outlier detection, candidate region clustering

### Validation with heterozygosity
- **STEP11** — Overlap candidates with theta/heterozygosity from `het_roh` workflow
- **STEP12** — Regional PCA, k-means groups, panel A/B/C plots per candidate
- **STEP13** — Combined multi-panel manuscript figure

---

## Key paths and integration with het_roh workflow

### Project root
```
/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04
```

### Inversion pipeline output
```
${PROJECT}/inversion_localpca/
├── chunks/             # STEP04 per-chunk SAF/SFS
├── global_sfs/         # STEP05-06 merged SAF + pest
├── 02_snps_beagle/     # STEP07 beagle output
├── 03_sites/           # STEP07 site lists
├── 04_dosage_by_chr/   # STEP08 dosage matrices
├── 05_local_pca/       # STEP09 window PCA
├── 06_mds_candidates/  # STEP10 MDS + candidates
├── 07_het_overlap/     # STEP11 theta overlap
├── 08_regional_pca/    # STEP12 per-candidate PCA
├── 09_combined_plots/  # STEP13 final figures
└── logs/
```

### Heterozygosity inputs (from het_roh workflow)
```
${PROJECT}/het_roh/02_heterozygosity/03_theta/multiscale/*.pestPG
${PROJECT}/het_roh/02_heterozygosity/04_summary/genomewide_heterozygosity.tsv
${PROJECT}/het_roh/01_inputs_check/samples_qcpass.txt
${PROJECT}/het_roh/01_inputs_check/samples.ind
```

### Key .pestPG column for this workflow
- **tP** = pairwise theta / pi-like diversity proxy (main per-window diversity track)
- Physical coordinates from **WinStart / WinStop / WinCenter** columns
- NOT firstPos_withData (which is data-dependent, not window-anchored)

---

## Callable mask strategy

- **strict `normalACGT`** mask — conservative global population structure
- **broader `inversion_acgt_allcase`** mask — inversion discovery (includes soft-masked)

---

## Core inversion filters

```
-GL 1  -minQ 20  -minMapQ 60  -minInd 203  -setMaxDepth 2850
-minMaf 0.05  -SNP_pval 1e-6
-remove_bads 1  -uniqueOnly 1  -only_proper_pairs 1
```

---

## Dosage formula

For each biallelic SNP from BEAGLE GLs:
```
dosage = P(AB) + 2 × P(BB)
```
This is the expected count of allele2, used directly for local PCA.

---

## Plot C interpretation

Panel C validates candidate inversions by showing:
1. Group individuals by regional PCA k-means (3 groups on PC1)
2. Plot chromosome-wide tP diversity curves for each group
3. If groups differ most strongly in the candidate region → supports inversion

Heterozygosity (tP) is allele-polarity invariant, so plot C is safe even if allele coding flips.
