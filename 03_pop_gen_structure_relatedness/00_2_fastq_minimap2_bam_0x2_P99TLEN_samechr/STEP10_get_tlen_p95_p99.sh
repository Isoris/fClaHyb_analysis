###
# get_tlen_p95_p99.sh does this:
# Computes approximate p95 and p99 absolute TLEN values from a properly paired,
# primary high-quality BAM using samtools + awk.
###
#!/usr/bin/env bash
set -euo pipefail

BAM="${1:-02-merged_per_sample/CGA109/CGA109.merged.markdup.clip.bam}"

samtools view -f 0x2 -F 0xF0C -q 60 "$BAM" \
| awk '$9 != 0 { t = $9; if (t < 0) t = -t; print t }' \
| sort -n \
| awk '{
    a[NR] = $1
}
END {
    n = NR
    if (n == 0) {
        print "ERROR: no TLEN values found" > "/dev/stderr"
        exit 1
    }
    p95 = int(n * 0.95); if (p95 < 1) p95 = 1
    p99 = int(n * 0.99); if (p99 < 1) p99 = 1
    print "p95=" a[p95], "p99=" a[p99]
}'
