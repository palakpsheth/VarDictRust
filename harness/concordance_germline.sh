#!/usr/bin/env bash
# Usage: concordance_germline.sh <name> <java.vcf> <rust.vcf> <truth.vcf.gz> <eval.bed> <ref.sdf> <outdir> <sample>
# Gate: identical 'None'-row counts AND identical tp/tp-baseline/fp/fn site lists for germline calls.
set -euo pipefail

if ! command -v rtg &>/dev/null || ! command -v bcftools &>/dev/null; then
  export PATH="/home/palak/miniconda3/envs/g4-tools/bin:$PATH"
fi

NAME=$1; J=$2; R=$3; TRUTH=$4; BED=$5; SDF=$6; OUT=$7; SAMPLE=$8
mkdir -p "$OUT/java" "$OUT/rust"

# 1. Prepare truth VCF (strip filters so rtg doesn't drop truth variants)
TRUTH_NOFILTER="$OUT/truth.nofilter.vcf.gz"
if [ ! -f "$TRUTH_NOFILTER" ] || [ ! -f "$TRUTH_NOFILTER.tbi" ]; then
  bcftools annotate -x FILTER "$TRUTH" | bcftools sort -Oz -o "$TRUTH_NOFILTER"
  tabix -f -p vcf "$TRUTH_NOFILTER"
fi

# 2. Process each engine
run_eval() {
  local eng=$1
  local vcf_in=$2
  local eng_out="$OUT/$eng"

  echo "Filtering $eng calls with -f PASS..."
  bcftools view -f PASS "$vcf_in" | bcftools sort -Oz -o "$eng_out/calls.pass.vcf.gz"
  tabix -f -p vcf "$eng_out/calls.pass.vcf.gz"

  rm -rf "$eng_out/vcfeval"
  echo "Running rtg vcfeval for $eng against truth..."
  rtg vcfeval -b "$TRUTH_NOFILTER" -c "$eng_out/calls.pass.vcf.gz" -t "$SDF" -e "$BED" \
    --sample "$SAMPLE,$SAMPLE" --vcf-score-field QUAL -T 8 -o "$eng_out/vcfeval"

  python3 - "$eng_out/vcfeval/summary.txt" > "$eng_out/summary.json" <<'PY'
import json, sys
rows = [l.split() for l in open(sys.argv[1]) if l.strip() and not l.startswith(("Threshold", "-"))]
none = [r for r in rows if r[0] == "None"]
if len(none) != 1: sys.exit("concordance_germline: no unique 'None' row in summary.txt")
r = none[0]
d = dict(tp_baseline=int(r[1]), tp_call=int(r[2]), fp=int(r[3]), fn=int(r[4]),
         precision=float(r[5]), sensitivity=float(r[6]), f_measure=float(r[7]))
if d["tp_baseline"] + d["fn"] == 0: sys.exit("concordance_germline: empty baseline in eval regions (vacuous)")
print(json.dumps(d, indent=2))
PY

  for f in tp tp-baseline fp fn; do
    [ -f "$eng_out/vcfeval/$f.vcf.gz" ] && bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$eng_out/vcfeval/$f.vcf.gz" > "$eng_out/$f.sites" || touch "$eng_out/$f.sites"
  done
}

run_eval "java" "$J"
run_eval "rust" "$R"

# 3. Compare java vs rust metrics and site lists
res=EQUAL
if ! cmp -s "$OUT/java/summary.json" "$OUT/rust/summary.json"; then
  echo "WARNING: summary.json differs between Java and Rust" >&2
  diff -u "$OUT/java/summary.json" "$OUT/rust/summary.json" >&2 || true
  res=DIFFERENT
fi

for f in tp tp-baseline fp fn; do
  if ! cmp -s "$OUT/java/$f.sites" "$OUT/rust/$f.sites"; then
    echo "WARNING: $f.sites differs between Java and Rust" >&2
    diff -u "$OUT/java/$f.sites" "$OUT/rust/$f.sites" | head -n 20 >&2 || true
    res=DIFFERENT
  fi
done

python3 -c 'import json,sys; o=sys.argv[1]; json.dump({"java": json.load(open(o+"/java/summary.json")), "rust": json.load(open(o+"/rust/summary.json"))}, open(o+"/extra.json","w"), indent=2, sort_keys=True)' "$OUT"
python3 harness/record_evidence.py "$NAME" "$res" "$OUT" -- "$0" "$@"

echo "Germline concordance result: $res"
[ "$res" = EQUAL ]
