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

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/00_delly_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a
source "${CONFIG}"
set +a

dv_init_dirs
dv_log "=== ANNOTATION LAYERS ==="

DEL_BED_226="${DIR_FINAL}/catalog_226.DEL.bed"
dv_check_file "${DEL_BED_226}" "226 DEL BED"
DEL_BED_81=$(ls "${DIR_FINAL}"/catalog_81.DEL.*.PASS.bed 2>/dev/null | head -1)
ANNOT_BEDS="${DIR_ANNOT}/beds"
export DEL_BED_226



# ── LAYER 1: Gene/Exon/CDS ─────────────────────────────────────────────────
dv_log "--- LAYER 1: Functional overlap ---"

annotate_functional() {
  local del_bed="$1" label="$2" outpfx="${DIR_ANNOT}/${2}"

  for feat in GENE EXON CDS; do
    local fb="${ANNOT_BEDS}/${feat}.sorted.bed"
    [[ -f "${fb}" ]] && bedtools intersect -a "${del_bed}" -b "${fb}" -wa -wb \
      > "${outpfx}.${feat}_overlap.tsv" 2>/dev/null || true
  done

  python3 << PYEOF
from collections import defaultdict
del_bed="${del_bed}"; outpfx="${outpfx}"
dels={}
with open(del_bed) as f:
    for line in f:
        p=line.strip().split('\t'); key=f"{p[0]}:{p[1]}-{p[2]}"
        dels[key]={'chrom':p[0],'start':p[1],'end':p[2],'id':p[3] if len(p)>3 else '.','svlen':p[4] if len(p)>4 else '.','genes':set(),'exons':set(),'cds':set()}
def load(path,field):
    try:
        with open(path) as f:
            for l in f:
                p=l.strip().split('\t'); k=f"{p[0]}:{p[1]}-{p[2]}"
                if k in dels: dels[k][field].add(p[8] if len(p)>8 else p[7] if len(p)>7 else '.')
    except: pass
load(f"{outpfx}.GENE_overlap.tsv",'genes'); load(f"{outpfx}.EXON_overlap.tsv",'exons'); load(f"{outpfx}.CDS_overlap.tsv",'cds')
with open(f"{outpfx}.functional_class.tsv",'w') as o:
    o.write("del_key\tchrom\tstart\tend\tid\tsvlen\tn_genes\tn_exons\tn_cds\tclass\tgene_names\n")
    for k,d in dels.items():
        ng,ne,nc=len(d['genes']),len(d['exons']),len(d['cds'])
        c="CDS_overlap" if nc>0 else "exon_overlap" if ne>0 else "intronic" if ng>0 else "intergenic"
        g=','.join(sorted(d['genes'])) if d['genes'] else '.'
        o.write(f"{k}\t{d['chrom']}\t{d['start']}\t{d['end']}\t{d['id']}\t{d['svlen']}\t{ng}\t{ne}\t{nc}\t{c}\t{g}\n")
PYEOF


  dv_log "  ${label}:"; tail -n +2 "${outpfx}.functional_class.tsv" | cut -f10 | sort | uniq -c | sort -rn | while read c l; do dv_log "    ${l}: ${c}"; done
}

annotate_functional "${DEL_BED_226}" "catalog_226"
[[ -n "${DEL_BED_81:-}" && -f "${DEL_BED_81}" ]] && annotate_functional "${DEL_BED_81}" "catalog_81"

# ── LAYER 2: Repeat overlap ────────────────────────────────────────────────
dv_log "--- LAYER 2: Repeat overlap ---"
if [[ -f "${ANNOT_BEDS}/REPEATS.sorted.bed" ]]; then
  for label_bed in "catalog_226:${DEL_BED_226}" "catalog_81:${DEL_BED_81:-}"; do
    label="${label_bed%%:*}"; bed="${label_bed#*:}"
    [[ -n "$bed" && -f "$bed" ]] || continue
    bedtools intersect -a "$bed" -b "${ANNOT_BEDS}/REPEATS.sorted.bed" -wa -f 0.5 -u > "${DIR_ANNOT}/${label}.DELs_in_repeats.bed" 2>/dev/null || true
    bedtools intersect -a "$bed" -b "${ANNOT_BEDS}/REPEATS.sorted.bed" -wa -f 0.5 -v > "${DIR_ANNOT}/${label}.DELs_not_in_repeats.bed" 2>/dev/null || true
    dv_log "  ${label}: $(wc -l < "${DIR_ANNOT}/${label}.DELs_in_repeats.bed") in repeats"
  done
fi

# ── LAYER 3: Depth support ─────────────────────────────────────────────────
dv_log "--- LAYER 3: Depth support ---"
PA_MOSDEPTH_DIR="${DELLY_PROJECT}/pa_roary_results/01_mosdepth"
MOSDEPTH_SOURCE=""
for d in "${PA_MOSDEPTH_DIR}" "${DELLY_PROJECT}/pa_roary_results/01_mosdepth_callable"; do
  if ls "${d}"/*.regions.bed.gz &>/dev/null 2>&1 || ls "${d}"/*/*.regions.bed.gz &>/dev/null 2>&1; then
    MOSDEPTH_SOURCE="${d}"; break
  fi
done

if [[ -n "${MOSDEPTH_SOURCE}" ]]; then
  dv_log "  Using mosdepth from: ${MOSDEPTH_SOURCE}"
  export MOSDEPTH_SOURCE
  python3 << 'PYEOF'
import gzip,os,glob,statistics
from collections import defaultdict
del_bed=os.environ['DEL_BED_226']; md=os.environ['MOSDEPTH_SOURCE']
out=os.environ['DIR_DEPTH']+'/depth_support_226.tsv'; flank=5000
sites=[]
with open(del_bed) as f:
    for l in f:
        p=l.strip().split('\t'); sites.append((p[0],int(p[1]),int(p[2]),p[3] if len(p)>3 else f"{p[0]}:{p[1]}-{p[2]}"))
rfiles=glob.glob(os.path.join(md,'*.regions.bed.gz'))+glob.glob(os.path.join(md,'*','*.regions.bed.gz'))
rfiles=rfiles[:40]
def load(path):
    d=defaultdict(list)
    with gzip.open(path,'rt') as f:
        for l in f:
            p=l.strip().split('\t'); d[p[0]].append((int(p[1]),int(p[2]),float(p[3])))
    return d
def md_depth(r,c,s,e):
    if c not in r: return None
    tb=td=0
    for ws,we,wd in r[c]:
        o1,o2=max(ws,s),min(we,e)
        if o1<o2: bp=o2-o1; tb+=bp; td+=wd*bp
    return td/tb if tb>0 else None
dr=defaultdict(list)
for rf in rfiles:
    rg=load(rf)
    cm={c:statistics.median([d for _,_,d in w if d>0]) for c,w in rg.items() if any(d>0 for _,_,d in w)}
    for c,s,e,di in sites:
        if c not in cm or cm[c]==0: continue
        ins=md_depth(rg,c,s,e)
        if ins is None: continue
        ld=md_depth(rg,c,max(0,s-flank),s); rd_=md_depth(rg,c,e,e+flank)
        fl=[d for d in [ld,rd_] if d is not None]
        if not fl: continue
        fm=statistics.mean(fl)
        if fm==0: continue
        dr[di].append(ins/fm)
with open(out,'w') as o:
    o.write("del_id\tn_samples\tmedian_ratio\tmean_ratio\tdepth_label\n")
    for c,s,e,di in sites:
        rs=dr.get(di,[])
        if not rs: o.write(f"{di}\t0\tNA\tNA\tdepth_no_data\n"); continue
        mr=statistics.median(rs); mn=statistics.mean(rs)
        lb="strong_DEL_depth_support" if mr<0.3 else "moderate_DEL_depth_support" if mr<0.6 else "weak_DEL_depth_support" if mr<0.85 else "no_depth_support"
        o.write(f"{di}\t{len(rs)}\t{mr:.4f}\t{mn:.4f}\t{lb}\n")
PYEOF
else
  dv_log "  No mosdepth data found. Running fresh..."
  mkdir -p "${DIR_DEPTH}/mosdepth"
  run_md() {
    local s="$1" b="${DIR_MARKDUP}/${1}.markdup.bam" o="${DIR_DEPTH}/mosdepth/${1}"
    [[ -f "${o}.regions.bed.gz" ]] && return 0
    mosdepth --by "${DEPTH_WINDOW}" --mapq "${DEPTH_MAPQ}" --no-per-base --threads "${DEPTH_THREADS}" "${o}" "${b}" 2>>"${DIR_LOGS}/depth_${1}.log"
  }
  export -f run_md; export DIR_MARKDUP DIR_DEPTH DIR_LOGS DEPTH_WINDOW DEPTH_MAPQ DEPTH_THREADS
  parallel -j 20 run_md {} :::: "${SAMPLES_ALL}"
  MOSDEPTH_SOURCE="${DIR_DEPTH}/mosdepth"; export MOSDEPTH_SOURCE
  # Run same python as above
fi

if [[ -f "${DIR_DEPTH}/depth_support_226.tsv" ]]; then
  dv_log "  Depth support:"; tail -n +2 "${DIR_DEPTH}/depth_support_226.tsv" | cut -f5 | sort | uniq -c | sort -rn | while read c l; do dv_log "    ${l}: ${c}"; done
fi

# ── LAYER 4: SVLEN QC ─────────────────────────────────────────────────────
dv_log "--- LAYER 4: SVLEN QC ---"
FINAL_226_VCF="${DIR_FINAL}/catalog_226.DEL.vcf.gz"
export FINAL_226_VCF DIR_MATDIST MATE_WARN_KB MATE_SUSPICIOUS_KB MATE_EXTREME_KB
python3 << 'PYEOF'
import gzip,os
vcf=os.environ['FINAL_226_VCF']; out=os.environ['DIR_MATDIST']+'/mate_distance_qc_226.tsv'
w=int(os.environ['MATE_WARN_KB']); su=int(os.environ['MATE_SUSPICIOUS_KB']); ex=int(os.environ['MATE_EXTREME_KB'])
with gzip.open(vcf,'rt') as v, open(out,'w') as o:
    o.write("chrom\tpos\tend\tid\tsvlen_bp\tabs_kb\tflag\n")
    for l in v:
        if l.startswith('#'): continue
        p=l.strip().split('\t'); info=dict(f.split('=',1) for f in p[7].split(';') if '=' in f)
        sv=abs(int(info.get('SVLEN','0'))); kb=sv/1000
        fl="extreme_artifact_candidate" if kb>=ex else "very_suspicious" if kb>=su else "warning_large" if kb>=w else "normal"
        o.write(f"{p[0]}\t{p[1]}\t{info.get('END',p[1])}\t{p[2]}\t{sv}\t{kb:.2f}\t{fl}\n")
PYEOF


if [[ -f "${DIR_MATDIST}/mate_distance_qc_226.tsv" ]]; then
  dv_log "  SVLEN QC:"; tail -n +2 "${DIR_MATDIST}/mate_distance_qc_226.tsv" | cut -f7 | sort | uniq -c | sort -rn | while read c f; do dv_log "    ${f}: ${c}"; done
fi

dv_log "=== ANNOTATION COMPLETE ==="
