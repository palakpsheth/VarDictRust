#!/usr/bin/env bash
# Usage: concordance.sh <calls.vcf> <truth.vcf.gz> <eval.bed> <ref.sdf> <outdir> <tumor-sample>
# Calls -> FILTER=PASS & STATUS in {StrongSomatic,LikelySomatic}; truth FILTER stripped (else rtg drops it);
# rtg vcfeval with sites-only baseline (sample ALT). Prints the 'None' row as JSON; fails on empty baseline.
set -euo pipefail
CALLS=$1; TRUTH=$2; BED=$3; SDF=$4; OUT=$5; TS=$6; mkdir -p "$OUT"

VCFEVAL_STORE="${GOLDEN_STORE_DIR:-scratch/golden-store}/vcfeval"
GMODE="${GOLDEN_MODE:-auto}"
CACHE_KEY=$(python3 - "$CALLS" "$TRUTH" "$BED" "$SDF" "$TS" <<'PY'
import hashlib, os, subprocess, sys
calls, truth, bed, sdf, sample = sys.argv[1:6]
h = hashlib.sha256()
h.update(b"somatic\0" + sample.encode() + b"\0")
try:
    rtg_ver = subprocess.check_output(["rtg", "version"], text=True, stderr=subprocess.DEVNULL).splitlines()[0]
except Exception:
    rtg_ver = "unknown"
h.update(rtg_ver.encode() + b"\0")
for p in (bed, truth, os.path.join(sdf, "mainIndex")):
    if os.path.isfile(p):
        with open(p, "rb") as f:
            while chunk := f.read(1 << 20):
                h.update(chunk)
    h.update(b"\0")
with open(calls, "rb") as f:
    for line in f:
        if line.startswith((b"##fileDate=", b"##source=")):
            continue
        h.update(line)
print(h.hexdigest()[:32])
PY
)
CACHE_ENTRY="$VCFEVAL_STORE/$CACHE_KEY"

if [ "$GMODE" = "auto" ] && [ -s "$CACHE_ENTRY/summary.json" ] \
   && [ -f "$CACHE_ENTRY/tp.sites" ] && [ -f "$CACHE_ENTRY/tp-baseline.sites" ] \
   && [ -f "$CACHE_ENTRY/fp.sites" ] && [ -f "$CACHE_ENTRY/fn.sites" ] \
   && python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); sys.exit(0 if d.get("tp_baseline",0)+d.get("fn",0)>0 else 1)' "$CACHE_ENTRY/summary.json" 2>/dev/null; then
  cp "$CACHE_ENTRY/summary.json" "$OUT/summary.json"
  for f in tp tp-baseline fp fn; do
    cp "$CACHE_ENTRY/$f.sites" "$OUT/$f.sites"
  done
  cat "$OUT/summary.json"
  exit 0
fi

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

if [ "$GMODE" != "off" ] && [ -s "$OUT/summary.json" ]; then
  mkdir -p "$VCFEVAL_STORE"
  stage_dir=$(mktemp -d "$VCFEVAL_STORE/.tmp.XXXXXX")
  cp "$OUT/summary.json" "$stage_dir/summary.json"
  for f in tp tp-baseline fp fn; do
    cp "$OUT/$f.sites" "$stage_dir/$f.sites"
  done
  rm -rf "$CACHE_ENTRY"
  mv "$stage_dir" "$CACHE_ENTRY" 2>/dev/null || rm -rf "$stage_dir"
fi

cat "$OUT/summary.json"
