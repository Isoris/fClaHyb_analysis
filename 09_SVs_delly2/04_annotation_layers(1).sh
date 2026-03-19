#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 127
#SBATCH --mem=237G
#SBATCH -t 1-00:00:00
#SBATCH -J delly_annot
#SBATCH -o logs/04_annotation.%j.out
#SBATCH -e logs/04_annotation.%j.err
set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_delly_config.sh"

dv_init_dirs
dv_log "=== ANNOTATION LAYERS ==="

# ─────────────────────────────────────────────────────────────────────────────
# Locate DEL BEDs from step 8
# ─────────────────────────────────────────────────────────────────────────────
DEL_BED_226="${DIR_FINAL}/catalog_226.DEL.bed"
dv_check_file "${DEL_BED_226}" "226 DEL BED"

DEL_BED_81=$(ls "${DIR_FINAL}"/catalog_81.DEL.*.PASS.bed 2>/dev/null | head -1)
ANNOT_BEDS="${DIR_ANNOT}/beds"

# =============================================================================
# LAYER 1: Gene / Exon / CDS functional overlap
# =============================================================================
dv_log "--- LAYER 1: Functional overlap ---"

annotate_functional() {
  local del_bed="$1" label="$2"
  local outpfx="${DIR_ANNOT}/${label}"

  # bedtools intersect
  for feat in GENE EXON CDS; do
    local feat_bed="${ANNOT_BEDS}/${feat}.sorted.bed"
    if [[ -f "${feat_bed}" ]]; then
      bedtools intersect -a "${del_bed}" -b "${feat_bed}" -wa -wb \
        > "${outpfx}.${feat}_overlap.tsv" 2>/dev/null || true
    fi
  done

  # Classify each DEL
  python3 << PYEOF
from collections import defaultdict

del_bed = "${del_bed}"
outpfx = "${outpfx}"

dels = {}
with open(del_bed) as f:
    for line in f:
        p = line.strip().split('\t')
        key = f"{p[0]}:{p[1]}-{p[2]}"
        dels[key] = {'chrom':p[0], 'start':p[1], 'end':p[2],
                     'id':p[3] if len(p)>3 else '.', 'svlen':p[4] if len(p)>4 else '.',
                     'genes':set(), 'exons':set(), 'cds':set()}

def load_ovl(path, field):
    try:
        with open(path) as f:
            for line in f:
                p = line.strip().split('\t')
                key = f"{p[0]}:{p[1]}-{p[2]}"
                if key in dels:
                    name = p[8] if len(p)>8 else p[7] if len(p)>7 else '.'
                    dels[key][field].add(name)
    except FileNotFoundError:
        pass

load_ovl(f"{outpfx}.GENE_overlap.tsv", 'genes')
load_ovl(f"{outpfx}.EXON_overlap.tsv", 'exons')
load_ovl(f"{outpfx}.CDS_overlap.tsv", 'cds')

with open(f"{outpfx}.functional_class.tsv", 'w') as out:
    out.write("del_key\tchrom\tstart\tend\tid\tsvlen\tn_genes\tn_exons\tn_cds\tclass\tgene_names\n")
    for key, d in dels.items():
        ng, ne, nc = len(d['genes']), len(d['exons']), len(d['cds'])
        cls = "CDS_overlap" if nc>0 else "exon_overlap" if ne>0 else "intronic" if ng>0 else "intergenic"
        gs = ','.join(sorted(d['genes'])) if d['genes'] else '.'
        out.write(f"{key}\t{d['chrom']}\t{d['start']}\t{d['end']}\t{d['id']}\t{d['svlen']}\t{ng}\t{ne}\t{nc}\t{cls}\t{gs}\n")
PYEOF

  dv_log "  ${label} functional:"
  tail -n +2 "${outpfx}.functional_class.tsv" | cut -f10 | sort | uniq -c | sort -rn \
    | while read cnt cls; do dv_log "    ${cls}: ${cnt}"; done
}

annotate_functional "${DEL_BED_226}" "catalog_226"
[[ -n "${DEL_BED_81:-}" && -f "${DEL_BED_81}" ]] && annotate_functional "${DEL_BED_81}" "catalog_81"

# =============================================================================
# LAYER 2: Repeat overlap
# =============================================================================
dv_log "--- LAYER 2: Repeat overlap ---"

annotate_repeats() {
  local del_bed="$1" label="$2"
  local outpfx="${DIR_ANNOT}/${label}"

  if [[ -f "${ANNOT_BEDS}/REPEATS.sorted.bed" ]]; then
    # DELs with >=50% reciprocal repeat overlap
    bedtools intersect -a "${del_bed}" -b "${ANNOT_BEDS}/REPEATS.sorted.bed" \
      -wa -f 0.5 -u > "${outpfx}.DELs_in_repeats.bed" 2>/dev/null || true
    bedtools intersect -a "${del_bed}" -b "${ANNOT_BEDS}/REPEATS.sorted.bed" \
      -wa -f 0.5 -v > "${outpfx}.DELs_not_in_repeats.bed" 2>/dev/null || true

    local n_in n_out
    n_in=$(wc -l < "${outpfx}.DELs_in_repeats.bed")
    n_out=$(wc -l < "${outpfx}.DELs_not_in_repeats.bed")
    dv_log "  ${label}: ${n_in} in repeats, ${n_out} outside"
  else
    dv_log "  ${label}: No repeat BED, skipping"
  fi
}

annotate_repeats "${DEL_BED_226}" "catalog_226"
[[ -n "${DEL_BED_81:-}" && -f "${DEL_BED_81}" ]] && annotate_repeats "${DEL_BED_81}" "catalog_81"

# =============================================================================
# LAYER 3: Depth support using PA-Roary mosdepth data
# =============================================================================
dv_log "--- LAYER 3: Depth support (from PA-Roary mosdepth) ---"

# Reuse existing mosdepth data from PA-Roary instead of re-running
PA_MOSDEPTH_DIR="${DELLY_PROJECT}/pa_roary_results/01_mosdepth"

# Check if we have per-sample mosdepth regions files
if ls "${PA_MOSDEPTH_DIR}"/*.regions.bed.gz &>/dev/null 2>&1; then
  MOSDEPTH_SOURCE="${PA_MOSDEPTH_DIR}"
  dv_log "  Using PA-Roary mosdepth data from: ${MOSDEPTH_SOURCE}"
elif ls "${PA_MOSDEPTH_DIR}"/*/*.regions.bed.gz &>/dev/null 2>&1; then
  # Try subdirectories
  MOSDEPTH_SOURCE="${PA_MOSDEPTH_DIR}"
  dv_log "  Using PA-Roary mosdepth data from: ${MOSDEPTH_SOURCE} (subdirs)"
else
  dv_log "  WARNING: No PA-Roary mosdepth regions found. Running mosdepth fresh."
  MOSDEPTH_SOURCE=""
fi

# Also check if the callable mosdepth was done differently
PA_MOSDEPTH_CALLABLE="${DELLY_PROJECT}/pa_roary_results/01_mosdepth_callable"
if [[ -z "${MOSDEPTH_SOURCE}" ]] && ls "${PA_MOSDEPTH_CALLABLE}"/*.regions.bed.gz &>/dev/null 2>&1; then
  MOSDEPTH_SOURCE="${PA_MOSDEPTH_CALLABLE}"
  dv_log "  Using mosdepth_callable data from: ${MOSDEPTH_SOURCE}"
fi

if [[ -n "${MOSDEPTH_SOURCE}" ]]; then
  # Use existing data
  python3 << 'PYEOF'
import gzip, os, glob, statistics
from collections import defaultdict

del_bed = os.environ['DEL_BED_226']
mosdepth_dir = os.environ.get('MOSDEPTH_SOURCE', '')
samples_file = os.environ['SAMPLES_ALL']
out_file = os.environ['DIR_DEPTH'] + '/depth_support_226.tsv'
flank_bp = 5000

# Read DEL sites
del_sites = []
with open(del_bed) as f:
    for line in f:
        p = line.strip().split('\t')
        del_sites.append((p[0], int(p[1]), int(p[2]), p[3] if len(p)>3 else f"{p[0]}:{p[1]}-{p[2]}"))
print(f"DEL sites: {len(del_sites)}")

# Find mosdepth regions files
region_files = glob.glob(os.path.join(mosdepth_dir, '*.regions.bed.gz'))
if not region_files:
    region_files = glob.glob(os.path.join(mosdepth_dir, '*', '*.regions.bed.gz'))
print(f"Found {len(region_files)} mosdepth regions files")

# Use up to 40 samples for speed
region_files = region_files[:40]

def load_regions(path):
    data = defaultdict(list)
    with gzip.open(path, 'rt') as f:
        for line in f:
            p = line.strip().split('\t')
            data[p[0]].append((int(p[1]), int(p[2]), float(p[3])))
    return data

def mean_depth(regions, chrom, start, end):
    if chrom not in regions:
        return None
    total_bp = total_d = 0
    for ws, we, wd in regions[chrom]:
        os_ = max(ws, start)
        oe = min(we, end)
        if os_ < oe:
            bp = oe - os_
            total_bp += bp
            total_d += wd * bp
    return total_d / total_bp if total_bp > 0 else None

del_ratios = defaultdict(list)

for i, rfile in enumerate(region_files):
    if i % 10 == 0:
        print(f"  Processing sample {i+1}/{len(region_files)}...")
    regions = load_regions(rfile)
    # Chromosome medians
    chr_med = {}
    for c, wins in regions.items():
        ds = [d for _,_,d in wins if d > 0]
        if ds:
            chr_med[c] = statistics.median(ds)

    for chrom, start, end, del_id in del_sites:
        med = chr_med.get(chrom)
        if not med or med == 0:
            continue
        inside = mean_depth(regions, chrom, start, end)
        if inside is None:
            continue
        left_d = mean_depth(regions, chrom, max(0, start - flank_bp), start)
        right_d = mean_depth(regions, chrom, end, end + flank_bp)
        flanks = [d for d in [left_d, right_d] if d is not None]
        if not flanks:
            continue
        flank_m = statistics.mean(flanks)
        if flank_m == 0:
            continue
        del_ratios[del_id].append(inside / flank_m)

with open(out_file, 'w') as out:
    out.write("del_id\tn_samples\tmedian_ratio\tmean_ratio\tdepth_label\n")
    for chrom, start, end, del_id in del_sites:
        ratios = del_ratios.get(del_id, [])
        if not ratios:
            out.write(f"{del_id}\t0\tNA\tNA\tdepth_no_data\n")
            continue
        med_r = statistics.median(ratios)
        mean_r = statistics.mean(ratios)
        if med_r < 0.3:
            label = "strong_DEL_depth_support"
        elif med_r < 0.6:
            label = "moderate_DEL_depth_support"
        elif med_r < 0.85:
            label = "weak_DEL_depth_support"
        else:
            label = "no_depth_support"
        out.write(f"{del_id}\t{len(ratios)}\t{med_r:.4f}\t{mean_r:.4f}\t{label}\n")

print(f"Written: {out_file}")
PYEOF

  export DEL_BED_226 MOSDEPTH_SOURCE SAMPLES_ALL DIR_DEPTH
else
  # Fall back: run mosdepth fresh
  dv_log "  Running mosdepth (${DEPTH_WINDOW} bp windows)..."
  mkdir -p "${DIR_DEPTH}/mosdepth"

  run_depth_sample() {
    local sid="$1"
    local bam="${BAMDIR}/${sid}${BAM_SUFFIX}"
    local outpfx="${DIR_DEPTH}/mosdepth/${sid}"
    [[ -f "${outpfx}.regions.bed.gz" ]] && return 0
    mosdepth --by "${DEPTH_WINDOW}" --mapq "${DEPTH_MAPQ}" \
      --no-per-base --threads "${DEPTH_THREADS}" \
      "${outpfx}" "${bam}" 2>>"${DIR_LOGS}/depth_${sid}.log"
  }

  export -f run_depth_sample
  export BAMDIR BAM_SUFFIX DIR_DEPTH DIR_LOGS DEPTH_WINDOW DEPTH_MAPQ DEPTH_THREADS

  parallel -j 20 --joblog "${DIR_LOGS}/depth_parallel.log" \
    run_depth_sample {} :::: "${SAMPLES_ALL}"

  MOSDEPTH_SOURCE="${DIR_DEPTH}/mosdepth"
  export MOSDEPTH_SOURCE

  # Same python as above
  python3 << 'PYEOF'
import gzip, os, glob, statistics
from collections import defaultdict

del_bed = os.environ['DEL_BED_226']
mosdepth_dir = os.environ['MOSDEPTH_SOURCE']
out_file = os.environ['DIR_DEPTH'] + '/depth_support_226.tsv'
flank_bp = 5000

del_sites = []
with open(del_bed) as f:
    for line in f:
        p = line.strip().split('\t')
        del_sites.append((p[0], int(p[1]), int(p[2]), p[3] if len(p)>3 else f"{p[0]}:{p[1]}-{p[2]}"))

region_files = glob.glob(os.path.join(mosdepth_dir, '*.regions.bed.gz'))[:40]

def load_regions(path):
    data = defaultdict(list)
    with gzip.open(path, 'rt') as f:
        for line in f:
            p = line.strip().split('\t')
            data[p[0]].append((int(p[1]), int(p[2]), float(p[3])))
    return data

def mean_depth(regions, chrom, start, end):
    if chrom not in regions: return None
    total_bp = total_d = 0
    for ws, we, wd in regions[chrom]:
        os_ = max(ws, start); oe = min(we, end)
        if os_ < oe:
            bp = oe - os_; total_bp += bp; total_d += wd * bp
    return total_d / total_bp if total_bp > 0 else None

del_ratios = defaultdict(list)
for rfile in region_files:
    regions = load_regions(rfile)
    chr_med = {}
    for c, wins in regions.items():
        ds = [d for _,_,d in wins if d > 0]
        if ds: chr_med[c] = statistics.median(ds)
    for chrom, start, end, del_id in del_sites:
        med = chr_med.get(chrom)
        if not med: continue
        inside = mean_depth(regions, chrom, start, end)
        if inside is None: continue
        left_d = mean_depth(regions, chrom, max(0, start - flank_bp), start)
        right_d = mean_depth(regions, chrom, end, end + flank_bp)
        flanks = [d for d in [left_d, right_d] if d is not None]
        if not flanks: continue
        flank_m = statistics.mean(flanks)
        if flank_m == 0: continue
        del_ratios[del_id].append(inside / flank_m)

with open(out_file, 'w') as out:
    out.write("del_id\tn_samples\tmedian_ratio\tmean_ratio\tdepth_label\n")
    for chrom, start, end, del_id in del_sites:
        ratios = del_ratios.get(del_id, [])
        if not ratios:
            out.write(f"{del_id}\t0\tNA\tNA\tdepth_no_data\n"); continue
        med_r = statistics.median(ratios); mean_r = statistics.mean(ratios)
        label = "strong_DEL_depth_support" if med_r < 0.3 else \
                "moderate_DEL_depth_support" if med_r < 0.6 else \
                "weak_DEL_depth_support" if med_r < 0.85 else "no_depth_support"
        out.write(f"{del_id}\t{len(ratios)}\t{med_r:.4f}\t{mean_r:.4f}\t{label}\n")
PYEOF

  export DEL_BED_226 DIR_DEPTH
fi

dv_log "  Depth support counts:"
if [[ -f "${DIR_DEPTH}/depth_support_226.tsv" ]]; then
  tail -n +2 "${DIR_DEPTH}/depth_support_226.tsv" | cut -f5 | sort | uniq -c | sort -rn \
    | while read cnt lbl; do dv_log "    ${lbl}: ${cnt}"; done
fi

# =============================================================================
# LAYER 4: Mate distance / SVLEN QC
# =============================================================================
dv_log "--- LAYER 4: SVLEN / mate distance QC ---"

FINAL_226_VCF="${DIR_FINAL}/catalog_226.DEL.vcf.gz"

python3 << 'PYEOF'
import gzip, os

vcf = os.environ['FINAL_226_VCF']
out = os.environ['DIR_MATDIST'] + '/mate_distance_qc_226.tsv'
warn = int(os.environ['MATE_WARN_KB'])
susp = int(os.environ['MATE_SUSPICIOUS_KB'])
extr = int(os.environ['MATE_EXTREME_KB'])

with gzip.open(vcf, 'rt') as v, open(out, 'w') as o:
    o.write("chrom\tpos\tend\tid\tsvlen_bp\tabs_svlen_kb\tmate_dist_flag\n")
    for line in v:
        if line.startswith('#'): continue
        p = line.strip().split('\t')
        info = dict(f.split('=',1) for f in p[7].split(';') if '=' in f)
        end_val = info.get('END', p[1])
        svlen = abs(int(info.get('SVLEN', '0')))
        kb = svlen / 1000
        flag = "extreme_artifact_candidate" if kb >= extr else \
               "very_suspicious" if kb >= susp else \
               "warning_large" if kb >= warn else "normal"
        o.write(f"{p[0]}\t{p[1]}\t{end_val}\t{p[2]}\t{svlen}\t{kb:.2f}\t{flag}\n")
PYEOF

export FINAL_226_VCF DIR_MATDIST MATE_WARN_KB MATE_SUSPICIOUS_KB MATE_EXTREME_KB

dv_log "  Mate distance QC:"
if [[ -f "${DIR_MATDIST}/mate_distance_qc_226.tsv" ]]; then
  tail -n +2 "${DIR_MATDIST}/mate_distance_qc_226.tsv" | cut -f7 | sort | uniq -c | sort -rn \
    | while read cnt flg; do dv_log "    ${flg}: ${cnt}"; done
fi

dv_log "=== ANNOTATION LAYERS COMPLETE ==="
dv_log "Next: submit 05_summary_report.sh"
