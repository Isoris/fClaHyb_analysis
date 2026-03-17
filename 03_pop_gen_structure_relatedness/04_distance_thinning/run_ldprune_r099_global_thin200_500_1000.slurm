###
# run_ldprune_r099_global_thin200_500_1000.slurm does this:
# Runs ngsLD and prune_graph.pl on whole-genome beagle files for 200/500/1000 bp
# thinned datasets, keeps markers with pairwise r2 below the pruning threshold,
# and writes pruned beagle files plus summary/QC outputs.
###
#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 80
#SBATCH --mem=237GB
#SBATCH -t 23:00:00
#SBATCH -J ldprune_r099
#SBATCH -o ldprune_r099.%j.out
#SBATCH -e ldprune_r099.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly
  module load PrgEnv-gnu
  module load gcc/11.2.0
  module load zlib/1.2.13
  module load bzip2/1.0.8-cpeGNU-23.03
  module load XZ/5.4.3-cpeGNU-23.03
  module load cURL/8.1.2
  module unload GSL 2>/dev/null || true
#  module load GSL/2.7-cpeGNU-23.03
module load GSL/2.8-cpeGNU-25.03

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/popstruct_thin"

THREADS=80
MAX_KB=100
THINS=(200 500 1000)

OUTDIR="${BASE}/06_ldprune_global_WGthin"
mkdir -p "${OUTDIR}"

LOGDIR="${OUTDIR}/logs"
mkdir -p "${LOGDIR}"
RUNLOG="${LOGDIR}/ldprune_r099.${SLURM_JOB_ID}.log"
exec > >(tee -a "$RUNLOG") 2>&1

log () {
  echo "[$(date '+%F %T')] $*"
}

PRUNE_PL="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/ngsTools/ngsLD/scripts/prune_graph.pl"
[[ -s "$PRUNE_PL" ]] || { echo "[ERROR] missing prune_graph.pl: $PRUNE_PL" >&2; exit 1; }

MANIFEST="${OUTDIR}/07_files_and_paths.tsv"
{
  echo -e "thin\tlabel\tpath\texists\tbytes"
  for W in "${THINS[@]}"; do
    BEAGLE="${BASE}/04_beagle_byRF_majmin/catfish.wholegenome.byRF.thin_${W}.beagle.gz"
    if [[ -e "$BEAGLE" ]]; then
      BYTES=$(ls -ln "$BEAGLE" | awk '{print $5}')
      EXISTS="yes"
    else
      BYTES="NA"
      EXISTS="no"
    fi
    echo -e "${W}\tbeagle\t${BEAGLE}\t${EXISTS}\t${BYTES}"
  done

  echo -e "all\tprune_graph.pl\t${PRUNE_PL}\t$( [[ -e "$PRUNE_PL" ]] && echo yes || echo no )\t$( [[ -e "$PRUNE_PL" ]] && ls -ln "$PRUNE_PL" | awk '{print $5}' || echo NA )"
  echo -e "all\toutdir\t${OUTDIR}\t$( [[ -d "$OUTDIR" ]] && echo yes || echo no )\tNA"
  echo -e "all\tlogdir\t${LOGDIR}\t$( [[ -d "$LOGDIR" ]] && echo yes || echo no )\tNA"
} > "$MANIFEST"

log "[INFO] Wrote input manifest: $MANIFEST"
cat "$MANIFEST"

for W in "${THINS[@]}"; do
  BEAGLE="${BASE}/04_beagle_byRF_majmin/catfish.wholegenome.byRF.thin_${W}.beagle.gz"
  [[ -s "$BEAGLE" ]] || { echo "[ERROR] missing beagle: $BEAGLE" >&2; exit 1; }
done

for W in "${THINS[@]}"; do
  log "============================================================"
  log "[INFO] Starting W=${W} bp"
  log "============================================================"

  BEAGLE="${BASE}/04_beagle_byRF_majmin/catfish.wholegenome.byRF.thin_${W}.beagle.gz"
  PREFIX="${OUTDIR}/catfish.wholegenome.byRF.thin_${W}.r099"

  [[ -s "$BEAGLE" ]] || { echo "[ERROR] missing beagle: $BEAGLE" >&2; exit 1; }
  gzip -t "$BEAGLE" || { echo "[ERROR] corrupt gzip: $BEAGLE" >&2; exit 2; }

  log "[INFO] Using BEAGLE = $BEAGLE"

  log "[INFO] Header preview for W=${W}:"
  set +o pipefail
  zcat "$BEAGLE" | awk '
  NR==1 { print "[BEAGLE_HEADER]", $1, $2, $3, "..."; next }
  NR==2 { print "[BEAGLE_FIRST_ROW]", $1, $2, $3; exit }
  '
  set -o pipefail

  zcat "$BEAGLE" | awk 'END{print NR-1}' > "${PREFIX}.before.total.txt"

  zcat "$BEAGLE" | awk '
  NR==1{next}
  {
    id=$1
    chr=id
    sub(/_[^_]+$/, "", chr)
    c[chr]++
  }
  END{
    for (k in c) print k "\t" c[k]
  }' | sort -k1,1 > "${PREFIX}.before.per_chr.tsv"

  zcat "$BEAGLE" | awk '
  NR==1{next}
  {
    id=$1
    pos=id
    sub(/^.*_/, "", pos)
    chr=id
    sub(/_[^_]+$/, "", chr)
    print chr "\t" pos
  }' > "${PREFIX}.pos"

  zcat "$BEAGLE" | awk '
  NR==1{next}
  {
    id=$1
    pos=id
    sub(/^.*_/, "", pos)
    chr=id
    sub(/_[^_]+$/, "", chr)
    print id "\t" chr "\t" pos
  }' > "${PREFIX}.marker_map.tsv"

  N_SITES=$(zcat "$BEAGLE" | awk 'END{print NR-1}')

  set +o pipefail
  N_IND=$(zcat "$BEAGLE" | awk 'NR==1{print (NF-3)/3; exit}')
  set -o pipefail

  log "[INFO] W=${W} N_SITES=$N_SITES"
  log "[INFO] W=${W} N_IND=$N_IND"

  log "[INFO] W=${W} starting ngsLD"
  ngsLD \
    --geno "$BEAGLE" \
    --probs \
    --n_ind "$N_IND" \
    --n_sites "$N_SITES" \
    --pos "${PREFIX}.pos" \
    --max_kb_dist "$MAX_KB" \
    --n_threads "$THREADS" \
    --out "${PREFIX}.ld"
  log "[INFO] W=${W} finished ngsLD"

  [[ -s "${PREFIX}.ld" ]] || { echo "[ERROR] ngsLD output missing: ${PREFIX}.ld" >&2; exit 3; }

  tail -n +2 "${PREFIX}.ld" > "${PREFIX}.ld.noheader"

  log "[INFO] W=${W} starting prune_graph.pl"
  /project/lt200308-agbsci/13-programs/mambaforge/bin/perl "$PRUNE_PL" \
    --in_file "${PREFIX}.ld.noheader" \
    --field_dist 3 \
    --max_kb_dist "$MAX_KB" \
    --field_weight 7 \
    --min_weight 0.99 \
    --out "${PREFIX}.keep.ids"
  log "[INFO] W=${W} finished prune_graph.pl"

  [[ -s "${PREFIX}.keep.ids" ]] || { echo "[ERROR] keep list missing: ${PREFIX}.keep.ids" >&2; exit 4; }

  awk '
  NR==FNR{
    keep[$1]=1
    next
  }
  ($1 in keep){
    print $2 "\t" $3
  }' "${PREFIX}.keep.ids" "${PREFIX}.marker_map.tsv" > "${PREFIX}.keep.pos"

  zcat "$BEAGLE" | awk '
  NR==FNR{
    keep[$1]=1
    next
  }
  NR==1{
    print
    next
  }
  {
    if ($1 in keep) print
  }' "${PREFIX}.keep.ids" - | gzip > "${PREFIX}.beagle.gz"

  gzip -t "${PREFIX}.beagle.gz" || { echo "[ERROR] pruned beagle is corrupt for W=${W}" >&2; exit 5; }

  zcat "${PREFIX}.beagle.gz" | awk 'END{print NR-1}' > "${PREFIX}.after.total.txt"

  zcat "${PREFIX}.beagle.gz" | awk '
  NR==1{next}
  {
    id=$1
    chr=id
    sub(/_[^_]+$/, "", chr)
    c[chr]++
  }
  END{
    for (k in c) print k "\t" c[k]
  }' | sort -k1,1 > "${PREFIX}.after.per_chr.tsv"

  BEFORE=$(cat "${PREFIX}.before.total.txt")
  AFTER=$(cat "${PREFIX}.after.total.txt")
  REMOVED=$((BEFORE - AFTER))

  {
    echo -e "metric\tcount"
    echo -e "before\t${BEFORE}"
    echo -e "after\t${AFTER}"
    echo -e "removed\t${REMOVED}"
  } > "${PREFIX}.summary.tsv"

  log "[DONE] LD pruning complete for W=${W}"
  log "[INFO] Summary for W=${W}:"
  cat "${PREFIX}.summary.tsv"

  log "[INFO] Output files for W=${W}:"
  echo "  ${PREFIX}.ld"
  echo "  ${PREFIX}.ld.noheader"
  echo "  ${PREFIX}.keep.ids"
  echo "  ${PREFIX}.keep.pos"
  echo "  ${PREFIX}.marker_map.tsv"
  echo "  ${PREFIX}.beagle.gz"
  echo "  ${PREFIX}.before.per_chr.tsv"
  echo "  ${PREFIX}.after.per_chr.tsv"
  echo "  ${PREFIX}.summary.tsv"
done

log "============================================================"
log "[DONE] All thinnings complete: ${THINS[*]}"
log "============================================================"
