# DELLY DEL Pipeline — Catfish Hatchery Cohort (v2)
## 226 samples + 81 unrelated | DELLY 1.7.3 | ~9X lcWGS

### Changes from v1
- **DELLY binary**: Uses compiled `/project/.../delly/src/delly` (not conda)
- **Module loads**: `HTSlib/1.17-cpeGNU-23.03` + `Boost/1.81.0-cpeGNU-23.03` in every script
- **`-t DEL`**: Restricts discovery to deletions at call time (no post-hoc bcftools subset)
- **`-h 4`**: Uses DELLY's built-in threading (4 threads/call × 30 parallel = 120 cores)
- **Exclude BED**: Built from `callable_bp_per_bin.tsv` — identifies universally uncallable 50-kb bins (callable_bp < 500), merges adjacent ones, keeps blocks ≥ 50 kb. Catches telomeric/centromeric dead zones empirically rather than scanning for N-blocks.
- **Depth support**: Reuses PA-Roary mosdepth data if available (no redundant mosdepth runs)
- **Removed intermediate DEL-subset step**: `-t DEL` makes the old `02_discovery_DEL/` directory unnecessary

### Quick Start

```bash
# 1. Copy to project
mkdir -p /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/delly_sv
cp *.sh README.md /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/delly_sv/

# 2. cd and edit config
cd /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/delly_sv/
vi 00_delly_config.sh   # CHECK: BAMDIR, BAM_SUFFIX

# 3. Verify
ls ${BAMDIR}/*${BAM_SUFFIX} | wc -l   # should be 226
/project/lt200308-agbsci/01-catfish_assembly/delly/src/delly   # should print help

# 4. Submit
bash run_delly_pipeline.sh
```

### Config checklist

| Variable | Expected | Verify |
|---|---|---|
| `BAMDIR` | Dir with 226 markdup BAMs | `ls $BAMDIR/*markdup.bam \| wc -l` |
| `BAM_SUFFIX` | `.markdup.bam` or similar | `ls $BAMDIR/ \| head` |
| `DELLY_BIN` | `/project/.../delly/src/delly` | Already set |
| `NATORA_KEEP` | 81-line keep list | `wc -l` |
| `PA_CALLABLE_BP` | `callable_bp_per_bin.tsv` | `head` should show chrom/start/end/callable_bp |

### Pipeline flow

```
01_prep          Build exclude BED from callable data, validate BAMs, sort annotation BEDs
02_discovery     delly call -t DEL per sample (30 parallel × 4 threads)
03_merge_geno    delly merge → delly call -v (regenotype) → bcftools merge → subset 81 → delly filter -f germline
04_annotation    Gene/exon/CDS overlap, repeat flags, depth support, SVLEN QC
05_summary       Collate counts, per-sample stats, private/shared spectrum
```

### Output structure
```
delly_sv/
├── 01_discovery/          Per-sample DEL BCFs
├── 02_merged_sites/       Shared DEL site list
├── 03_genotyped/          Per-sample regenotyped BCFs
├── 04_merged_cohort/      226-sample merged BCF
├── 05_subset_81/          81-sample subset
├── 06_germline_filtered/  Germline-filtered 81 BCF
├── 07_final_catalogs/     Final VCFs + GT matrices + BEDs
├── 08_annotation/         Functional + repeat annotation
├── 09_depth_support/      Depth-based QC
├── 10_mate_distance_qc/   SVLEN-based flags
├── 11_summary/            Report + counts
└── logs/
```
