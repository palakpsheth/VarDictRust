#!/usr/bin/env bash
# Usage: concordance.sh <calls.vcf> <truth.vcf.gz> <eval.bed> <ref.sdf> <outdir> <tumor-sample>
# Calls -> FILTER=PASS & STATUS in {StrongSomatic,LikelySomatic}; truth FILTER stripped (else rtg drops it);
# rtg vcfeval with sites-only baseline (sample ALT). Prints the 'None' row as JSON; fails on empty baseline.
set -euo pipefail
CALLS=$1; TRUTH=$2; BED=$3; SDF=$4; OUT=$5; TS=$6; mkdir -p "$OUT"
bcftools view -f PASS -i 'INFO/STATUS="StrongSomatic" || INFO/STATUS="LikelySomatic"' "$CALLS" \
  | bcftools sort -Oz -o "$OUT/calls.som.vcf.gz"; tabix -f -p vcf "$OUT/calls.som.vcf.gz"
bcftools annotate -x FILTER "$TRUTH" | bcftools sort -Oz -o "$OUT/truth.nofilter.vcf.gz"
tabix -f -p vcf "$OUT/truth.nofilter.vcf.gz"
rm -rf "$OUT/vcfeval"
rtg vcfeval -b "$OUT/truth.nofilter.vcf.gz" -c "$OUT/calls.som.vcf.gz" -t "$SDF" -e "$BED" \
  --squash-ploidy --sample "ALT,$TS" --vcf-score-field QUAL -T 8 -o "$OUT/vcfeval"
python3 - "$OUT/vcfeval/summary.txt" > "$OUT/summary.json" <<'PY'
import json, sys
rows = [l.split() for l in open(sys.argv[1]) if l.strip() and not l.startswith(("Threshold", "-"))]
none = [r for r in rows if r[0] == "None"]
if len(none) != 1: sys.exit("concordance: no unique 'None' row in summary.txt")
r = none[0]
d = dict(tp_baseline=int(r[1]), tp_call=int(r[2]), fp=int(r[3]), fn=int(r[4]),
         precision=float(r[5]), sensitivity=float(r[6]), f_measure=float(r[7]))
if d["tp_baseline"] + d["fn"] == 0: sys.exit("concordance: empty baseline in eval regions (vacuous)")
print(json.dumps(d, indent=2))
PY
for f in tp tp-baseline fp fn; do
  [ -f "$OUT/vcfeval/$f.vcf.gz" ] && bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$OUT/vcfeval/$f.vcf.gz" > "$OUT/$f.sites" || touch "$OUT/$f.sites"
done
cat "$OUT/summary.json"
