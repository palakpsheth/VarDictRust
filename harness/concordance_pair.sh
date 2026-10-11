#!/usr/bin/env bash
# Usage: concordance_pair.sh <name> <java.vcf> <rust.vcf> <truth.vcf.gz> <eval.bed> <ref.sdf> <outdir> <tumor>
# Gate (D-27): identical 'None'-row counts AND identical tp/tp-baseline/fp/fn site lists. Not a dF1 tolerance.
set -euo pipefail
NAME=$1; J=$2; R=$3; TRUTH=$4; BED=$5; SDF=$6; OUT=$7; TS=$8
bash harness/concordance.sh "$J" "$TRUTH" "$BED" "$SDF" "$OUT/java" "$TS" >/dev/null
bash harness/concordance.sh "$R" "$TRUTH" "$BED" "$SDF" "$OUT/rust" "$TS" >/dev/null
res=EQUAL; cmp -s "$OUT/java/summary.json" "$OUT/rust/summary.json" || res=DIFFERENT
for f in tp tp-baseline fp fn; do cmp -s "$OUT/java/$f.sites" "$OUT/rust/$f.sites" || res=DIFFERENT; done
python3 -c 'import json,sys; o=sys.argv[1]; json.dump({"java": json.load(open(o+"/java/summary.json")), "rust": json.load(open(o+"/rust/summary.json"))}, open(o+"/extra.json","w"), indent=2, sort_keys=True)' "$OUT"
python3 harness/record_evidence.py "$NAME" "$res" "$OUT" -- "$0" "$@"
[ "$res" = EQUAL ]
