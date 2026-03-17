# 01_saf_folded_sfs

This module builds the callable-site mask, computes chromosome-wise ANGSD SAF files,
merges them into a genome-wide SAF, and estimates the global folded site-frequency spectrum (SFS)
to generate a mean folded prior (`.pest`) for downstream ANGSD SNP discovery.

## Purpose

This module prepares the likelihood-based prior used later for bi-allelic SNP discovery.

It does this by:
- extracting callable/non-masked genomic intervals from the reference FASTA
- converting those intervals into ANGSD callable-site format
- creating one chromosome/region file (`.rf.txt`) per reference sequence
- running ANGSD SAF estimation per chromosome/region
- merging all chromosome-level SAFs into one genome-wide SAF
- estimating the global folded bootstrap SFS
- averaging bootstrap SFS replicates into a single mean folded prior (`.pest`)

## Main scripts

### STEP00_mask_regions_from_fasta.py
Extracts interval tracks from the reference FASTA:
- hard-masked regions (`N/n`)
- soft-masked regions (`a/c/g/t`)
- normal callable regions (`A/C/G/T`)
- masked regions combined

Main output used downstream:
- `<prefix>.normalACGT.bed`

Notes:
- output BED is 0-based, half-open
- this step does not run ANGSD

### STEP01_make_callable_angsd_mask.sh
Converts the callable BED produced in STEP00 into ANGSD callable-site format and indexes it.

Typical logic:
- input: normal callable BED
- convert start coordinate from 0-based to 1-based
- write ANGSD interval file
- run `angsd sites index`

Expected outputs:
- callable ANGSD sites file
- `.idx`
- `.bin`

Important:
- the final callable filename produced here must exactly match the filename expected in STEP03
- if a contig renaming step is needed to match BAM/reference naming, it must happen before indexing

### STEP02_make_chr_rf_chunk_list.sh
Creates one region file per chromosome/scaffold from the reference FASTA and writes a master chunk list.

Outputs:
- `rf_files/*.rf.txt`
- `chunk_rf.list`

Purpose:
- allows SLURM array jobs to run one chromosome/region at a time

### STEP03_run_angsd_saf_chunks.slurm
Runs ANGSD per chromosome/region using the chunk list.

For each chunk:
- reads one `.rf.txt`
- runs `angsd -doSaf`
- writes one per-chromosome SAF
- runs `realSFS` bootstrap for that chunk

Typical outputs per chromosome:
- `chunks/<TAG>/01_saf/catfish.<TAG>.saf.idx`
- `chunks/<TAG>/01_saf/catfish.<TAG>.saf.gz`
- `chunks/<TAG>/01_saf/catfish.<TAG>.saf.pos.gz`
- `chunks/<TAG>/02_sfs/catfish.<TAG>.sfs`
- logs

Important:
- STEP03 depends on:
  - reference FASTA
  - BAM list
  - callable ANGSD sites file from STEP01
  - chunk list from STEP02
- the `BASE` path here must match the same project root used in STEP04 and STEP05

### STEP04_merge_chunk_saf_to_global.slurm
Collects all per-chromosome `.saf.idx` files and merges them into one genome-wide SAF.

This step:
- finds all `*.saf.idx`
- writes `global_sfs/saf.idx.list`
- runs `realSFS cat`

Outputs:
- `global_sfs/saf.idx.list`
- `global_sfs/catfish.global.saf.idx`
- `global_sfs/catfish.global.saf.gz`
- `global_sfs/catfish.global.saf.pos.gz`

Important:
- this step merges SAF files, not SFS files

### STEP05_make_global_folded_sfs_and_mean_pest.slurm
Runs `realSFS` on the merged genome-wide SAF to estimate the global folded bootstrap SFS,
then averages bootstrap replicates into one mean folded prior.

Outputs:
- `global_sfs/catfish.global.folded.sfs`
- `global_sfs/realSFS.global.folded.log`
- `global_sfs/catfish.global.folded.mean.pest`

Purpose:
- this mean folded `.pest` file is the main prior used by the downstream SNP discovery module

## Run order

1. Run STEP00 to extract callable/non-masked intervals from the reference FASTA
2. Run STEP01 to convert callable intervals into ANGSD callable-site format and index them
3. Run STEP02 to create chromosome `.rf.txt` files and `chunk_rf.list`
4. Submit STEP03 as a SLURM array across all chunks
5. Run STEP04 to merge chunk SAFs into one genome-wide SAF
6. Run STEP05 to estimate the global folded bootstrap SFS and mean folded prior

## Typical submission pattern

### STEP02
Run once to create the chunk list.

### STEP03
Count the number of chunks and submit as an array:

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
OUTDIR="${BASE}/popstruct_global"
N=$(wc -l < "${OUTDIR}/chunk_rf.list")
echo "N=$N"

sbatch --array=0-$((N-1))%8 STEP03_run_angsd_saf_chunks.slurm "${OUTDIR}/chunk_rf.list"

### STEP04
Run once after all chunk jobs complete.

### STEP05
Run once after STEP04 completes successfully.

## Inputs

Core inputs required by this module:
- reference FASTA
- BAM list
- callable/non-masked BED derived from the FASTA
- callable ANGSD sites file
- chromosome/region chunk list

Typical reference-related files:
- reference FASTA (`.fa`)
- reference index (`.fai`)
- callable BED from STEP00
- callable ANGSD sites file from STEP01

## Important outputs

Most important reusable outputs from this module:
- callable ANGSD sites file
- `chunk_rf.list`
- per-chromosome SAF files
- genome-wide merged SAF
- global folded bootstrap SFS
- `catfish.global.folded.mean.pest`

Most important final file for downstream analyses:
- `catfish.global.folded.mean.pest`

## Notes

- This module was built because computing the full SAF directly was not the preferred final workflow; instead, SAF was estimated per chromosome and merged afterward.
- The genome-wide folded SFS is estimated from the merged genome-wide SAF, not by concatenating per-chromosome SFS files.
- The per-chromosome SFS files are mainly intermediate/QC outputs.
- Bootstrap replicates from the global folded SFS are averaged into one mean prior vector for reuse in later ANGSD runs.
- Coordinate consistency is critical:
  - BED from STEP00 is 0-based half-open
  - ANGSD callable-site intervals in STEP01 are converted to 1-based
- Naming consistency is also critical:
  - callable mask filename used in STEP01 must match the path expected in STEP03
  - project `BASE` path must be consistent across STEP03, STEP04, and STEP05

## Final goal of this module

Produce a robust genome-wide folded prior for downstream ANGSD SNP discovery and related population-genetic analyses.
