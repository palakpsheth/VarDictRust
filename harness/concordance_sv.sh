#!/usr/bin/env bash
# Usage: concordance_sv.sh <calls.vcf.gz> <truth.vcf.gz> <regions.bed> <outdir>
# SV Tier-3 gate spike (HG008 V0.6 stvar/CNV via truvari, NOT rtg).
# GATED on the V0.6 representation fix (sizeless <DEL>/<DUP>/<INV> discarded,
# CNV method none — eval BLOCKED upstream): this builds the gate, not the
# integration, so no V0.6 metric semantics are invented here — the follow-up
# owns operating points and summary fields once V0.6 unblocks.
set -euo pipefail
[ $# -eq 4 ] || { echo "Usage: concordance_sv.sh <calls.vcf.gz> <truth.vcf.gz> <regions.bed> <outdir>" >&2; exit 2; }
CALLS=$1; TRUTH=$2; REGIONS=$3; OUT=$4
case "$CALLS" in
  *.vcf.gz) ;;
  *) echo "calls must be bgzipped + indexed VCF (samtools bgzip + tabix), got: $CALLS" >&2; exit 2; ;;
esac
if [ ! -f "$CALLS.tbi" ] && [ ! -f "$CALLS.csi" ]; then
  echo "calls must be bgzipped + indexed VCF (samtools bgzip + tabix), got: $CALLS" >&2; exit 2
fi
if ! command -v truvari >/dev/null 2>&1; then
  echo "SKIP: truvari not installed (conda install -c bioconda truvari)"
  exit 0
fi
mkdir -p "$OUT"
VER=$(truvari version 2>&1 | head -n 1)
echo "$VER"
truvari bench -b "$TRUTH" -c "$CALLS" -o "$OUT/truvari" --includebed "$REGIONS"
python3 - "$OUT" "$VER" <<'PY'
import json, sys
out, ver = sys.argv[1], sys.argv[2]
d = {"tool": "truvari", "tool_version": ver, "status": "ran"}
json.dump(d, open(out + "/summary.json", "w"), indent=2)
PY
cat "$OUT/summary.json"
