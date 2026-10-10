#!/usr/bin/env bash
set -euo pipefail
[ -x harness/seqc2_align.sh ] || { echo "seqc2_align.sh not executable" >&2; exit 1; }
grep -q "bwa mem -M -t 6" harness/seqc2_align.sh
grep -q "samtools fixmate -m" harness/seqc2_align.sh
grep -q "samtools sort -@1 -m 2G" harness/seqc2_align.sh
grep -q "samtools markdup -@2" harness/seqc2_align.sh
grep -q "samtools index" harness/seqc2_align.sh
echo "PASS test_seqc2_align"
