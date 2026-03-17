###
# extract_fastp_qc_tables.sh does this:
# Parses fastp log files into a per-sample QC table and a compact numeric summary
# (min/mean/median/max) for key read-quality metrics.
###
#!/usr/bin/env bash
set -euo pipefail

LOGDIR="/scratch/lt200308-agbsci/Quentin_project/00-samples/fastp/logs"
OUTTABLE="${LOGDIR}/fastp_qc_table.tsv"
OUTSUM="${LOGDIR}/fastp_qc_summary.tsv"

cd "$LOGDIR"

# ---- 1) Build per-log table -------------------------------------------------
awk -v OFS="\t" '
function pct_from_line(line,   m){
  if (match(line, /\(([0-9.]+)%\)/, m)) return m[1];
  return "";
}
BEGIN{
  print "ID","Sample","Run","Lane",
        "Adapter_R1","Adapter_R2",
        "R1_reads_before","R2_reads_before",
        "R1_Q20pct_before","R1_Q30pct_before","R2_Q20pct_before","R2_Q30pct_before",
        "R1_reads_after","R2_reads_after",
        "R1_Q20pct_after","R1_Q30pct_after","R2_Q20pct_after","R2_Q30pct_after",
        "Reads_passed","Fail_lowqual","Fail_toomanyN","Fail_tooshort",
        "Reads_adapter_trimmed","Bases_adapter_trimmed",
        "Reads_polyX","Bases_polyX",
        "Dup_rate_pct","Insert_peak"
}
FNR==1{
  id=FILENAME; sub(/\.fastp\.log$/,"",id)
  sample=run=lane=""
  split(id,a,"."); sample=a[1]; run=a[2]; lane=a[3]

  ad1=ad2=""
  r1rb=r2rb=r1ra=r2ra=""
  r1q20b=r1q30b=r2q20b=r2q30b=""
  r1q20a=r1q30a=r2q20a=r2q30a=""
  pass=flq=fn=fts=0
  r_ad=0; b_ad=0
  r_px=0; b_px=0
  dup=""; ins=""
}
/^No adapter detected for read1/ { ad1="No" }
/^No adapter detected for read2/ { ad2="No" }
/^Detected adapter for read1/     { ad1="Yes" }
/^Detected adapter for read2/     { ad2="Yes" }

/^Read1 before filtering:/ { mode="before"; which="R1"; next }
/^Read2 before filtering:/ { mode="before"; which="R2"; next }
/^Read1 after filtering:/  { mode="after";  which="R1"; next }
/^Read2 after filtering:/  { mode="after";  which="R2"; next }

/^[[:space:]]*total reads:/{
  gsub(/,/,"",$3)
  if (mode=="before" && which=="R1") r1rb=$3
  if (mode=="before" && which=="R2") r2rb=$3
  if (mode=="after"  && which=="R1") r1ra=$3
  if (mode=="after"  && which=="R2") r2ra=$3
}

/^[[:space:]]*Q20 bases:/{
  p=pct_from_line($0)
  if (mode=="before" && which=="R1") r1q20b=p
  if (mode=="before" && which=="R2") r2q20b=p
  if (mode=="after"  && which=="R1") r1q20a=p
  if (mode=="after"  && which=="R2") r2q20a=p
}
/^[[:space:]]*Q30 bases:/{
  p=pct_from_line($0)
  if (mode=="before" && which=="R1") r1q30b=p
  if (mode=="before" && which=="R2") r2q30b=p
  if (mode=="after"  && which=="R1") r1q30a=p
  if (mode=="after"  && which=="R2") r2q30a=p
}

/^reads passed filter:/{
  gsub(/,/,"",$4); pass=$4
}
/^reads failed due to low quality:/{
  gsub(/,/,"",$6); flq=$6
}
/^reads failed due to too many N:/{
  gsub(/,/,"",$7); fn=$7
}
/^reads failed due to too short:/{
  gsub(/,/,"",$6); fts=$6
}
/^reads with adapter trimmed:/{
  gsub(/,/,"",$5); r_ad=$5
}
/^bases trimmed due to adapters:/{
  gsub(/,/,"",$6); b_ad=$6
}
/^reads with polyX in 3'\'' end:/{
  gsub(/,/,"",$7); r_px=$7
}
/^bases trimmed in polyX tail:/{
  gsub(/,/,"",$7); b_px=$7
}

/^Duplication rate:/{
  dup=$3; gsub(/%/,"",dup)
}
/^Insert size peak/{
  ins=$NF
}

ENDFILE{
  if (ad1=="") ad1="NA"
  if (ad2=="") ad2="NA"
  print id,sample,run,lane,
        ad1,ad2,
        r1rb,r2rb,
        r1q20b,r1q30b,r2q20b,r2q30b,
        r1ra,r2ra,
        r1q20a,r1q30a,r2q20a,r2q30a,
        pass,flq,fn,fts,
        r_ad,b_ad,
        r_px,b_px,
        dup,ins
}
' *.fastp.log > "$OUTTABLE"

echo "[OK] Wrote table: $OUTTABLE"

# ---- 2) Make a small summary (min/median/max/mean) --------------------------
python3 - <<'PY'
import pandas as pd

table = pd.read_csv("fastp_qc_table.tsv", sep="\t")

cols = [
  "R1_Q30pct_before","R2_Q30pct_before",
  "R1_Q30pct_after","R2_Q30pct_after",
  "Fail_toomanyN",
  "Reads_adapter_trimmed",
  "Dup_rate_pct",
  "Insert_peak"
]

for c in cols:
  table[c] = pd.to_numeric(table[c], errors="coerce")

def summarize(series):
  s = series.dropna()
  if s.empty:
    return {"n":0,"min":None,"mean":None,"median":None,"max":None}
  return {
    "n": int(s.shape[0]),
    "min": float(s.min()),
    "mean": float(s.mean()),
    "median": float(s.median()),
    "max": float(s.max()),
  }

rows = []
for c in cols:
  d = summarize(table[c])
  rows.append({"metric": c, **d})

out = pd.DataFrame(rows)
out.to_csv("fastp_qc_summary.tsv", sep="\t", index=False)
print("[OK] Wrote summary: fastp_qc_summary.tsv")
PY

echo "[OK] Wrote summary: $OUTSUM"
