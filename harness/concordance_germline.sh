#!/usr/bin/env bash
# Usage: concordance_germline.sh <name> <java.vcf> <rust.vcf> <truth.vcf.gz> <eval.bed> <ref.sdf> <outdir> <sample> [cpp.vcf]
# Gate: identical 'None'-row counts AND identical tp/tp-baseline/fp/fn site lists for germline calls.
# Optional 9th arg adds the advisory vardictcpp third arm (pinned vardictcpp 3
# hc52dbad_0 bioconda): the same run_eval below re-runs rtg vcfeval on it with
# identical truth/BED/sdf. cpp F1 lands in extra.json (cpp_f1 + cpp_version) and
# stdout as a rank line; the gate and exit code stay java==rust only (Task-3
# ruling: header/END-key deltas are known-normalized, cpp rank never gated).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/portable.sh"
source "$SCRIPT_DIR/lib/concordance_cpp.sh"

if ! command -v rtg &>/dev/null || ! command -v bcftools &>/dev/null; then
  export PATH="$G4_TOOLS_BIN:$PATH"
fi

NAME=$1; J=$2; R=$3; TRUTH=$4; BED=$5; SDF=$6; OUT=$7; SAMPLE=$8; CPP_VCF=${9:-}
mkdir -p "$OUT/java" "$OUT/rust"
TRUTH_NOFILTER="$OUT/truth.nofilter.vcf.gz"
TRUTH_SHA=$(portable_sha256 "$TRUTH")

# 1. Process each engine
VCFEVAL_STORE="${GOLDEN_STORE_DIR:-scratch/golden-store}/vcfeval"
GMODE="${GOLDEN_MODE:-auto}"

run_eval() {
  local eng=$1
  local vcf_in=$2
  local eng_out="$OUT/$eng"
  mkdir -p "$eng_out"

  local cache_key
  cache_key=$(python3 - "$vcf_in" "$TRUTH" "$BED" "$SDF" "$SAMPLE" <<'PY'
import hashlib, os, subprocess, sys
calls, truth, bed, sdf, sample = sys.argv[1:6]
h = hashlib.sha256()
h.update(b"germline\0" + sample.encode() + b"\0")
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
  local cache_entry="$VCFEVAL_STORE/$cache_key"

  if [ "$GMODE" = "auto" ] && [ -s "$cache_entry/summary.json" ] \
     && [ -f "$cache_entry/tp.sites" ] && [ -f "$cache_entry/tp-baseline.sites" ] \
     && [ -f "$cache_entry/fp.sites" ] && [ -f "$cache_entry/fn.sites" ] \
     && python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); sys.exit(0 if d.get("tp_baseline",0)+d.get("fn",0)>0 else 1)' "$cache_entry/summary.json" 2>/dev/null; then
    echo "Reusing cached rtg vcfeval for $eng ($cache_key)..."
    cp "$cache_entry/summary.json" "$eng_out/summary.json"
    for f in tp tp-baseline fp fn; do
      cp "$cache_entry/$f.sites" "$eng_out/$f.sites"
    done
    return 0
  fi

  if [ ! -f "$TRUTH_NOFILTER" ] || [ ! -f "$TRUTH_NOFILTER.tbi" ] || [ "$(cat "$TRUTH_NOFILTER.sha256" 2>/dev/null || true)" != "$TRUTH_SHA" ]; then
    bcftools annotate -x FILTER "$TRUTH" | bcftools sort -Oz -o "$TRUTH_NOFILTER"
    tabix -f -p vcf "$TRUTH_NOFILTER"
    printf '%s\n' "$TRUTH_SHA" > "$TRUTH_NOFILTER.sha256"
  fi

  echo "Filtering $eng calls with -f PASS..." &&
  bcftools view -f PASS "$vcf_in" | bcftools sort -Oz -o "$eng_out/calls.pass.vcf.gz" &&
  tabix -f -p vcf "$eng_out/calls.pass.vcf.gz" &&
  rm -rf "$eng_out/vcfeval" &&
  echo "Running rtg vcfeval for $eng against truth..." &&
  rtg vcfeval -b "$TRUTH_NOFILTER" -c "$eng_out/calls.pass.vcf.gz" -t "$SDF" -e "$BED" \
    --sample "$SAMPLE,$SAMPLE" --vcf-score-field QUAL -T 8 -o "$eng_out/vcfeval" &&
  python3 - "$eng_out/vcfeval/summary.txt" > "$eng_out/summary.json" <<'PY' &&
import json, sys
rows = [l.split() for l in open(sys.argv[1]) if l.strip() and not l.startswith(("Threshold", "-"))]
none_rows = [r for r in rows if r[0] == "None"]
if len(none_rows) != 1: sys.exit("concordance_germline: no unique 'None' row in summary.txt")
r_none = none_rows[0]

opt_rows = [r for r in rows if r[0] != "None"]
d_opt = None
if opt_rows:
    r_opt = opt_rows[0]
    d_opt = dict(threshold=float(r_opt[0]), tp_baseline=int(r_opt[1]), tp_call=int(r_opt[2]),
                 fp=int(r_opt[3]), fn=int(r_opt[4]), precision=float(r_opt[5]),
                 sensitivity=float(r_opt[6]), f_measure=float(r_opt[7]))

d = {
    "none": dict(tp_baseline=int(r_none[1]), tp_call=int(r_none[2]), fp=int(r_none[3]), fn=int(r_none[4]),
                 precision=float(r_none[5]), sensitivity=float(r_none[6]), f_measure=float(r_none[7])),
    "optimal": d_opt,
    "tp_baseline": int(r_none[1]), "tp_call": int(r_none[2]), "fp": int(r_none[3]), "fn": int(r_none[4]),
    "precision": float(r_none[5]), "sensitivity": float(r_none[6]), "f_measure": float(r_none[7])
}
if d["none"]["tp_baseline"] + d["none"]["fn"] == 0: sys.exit("concordance_germline: empty baseline in eval regions (vacuous)")
print(json.dumps(d, indent=2))
PY
  [ -s "$eng_out/summary.json" ] &&
  for f in tp tp-baseline fp fn; do
    if [ -f "$eng_out/vcfeval/$f.vcf.gz" ]; then
      bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$eng_out/vcfeval/$f.vcf.gz" > "$eng_out/$f.sites"
    else
      touch "$eng_out/$f.sites"
    fi
  done &&
  if [ "$GMODE" != "off" ] && [ -s "$eng_out/summary.json" ]; then
    mkdir -p "$VCFEVAL_STORE"
    local stage_dir
    stage_dir=$(mktemp -d "$VCFEVAL_STORE/.tmp.XXXXXX")
    cp "$eng_out/summary.json" "$stage_dir/summary.json"
    for f in tp tp-baseline fp fn; do
      cp "$eng_out/$f.sites" "$stage_dir/$f.sites"
    done
    rm -rf "$cache_entry"
    mv "$stage_dir" "$cache_entry" 2>/dev/null || rm -rf "$stage_dir"
  fi
}

run_eval "java" "$J"
run_eval "rust" "$R"

# 2b. Third arm: vardictcpp (advisory only via lib/concordance_cpp.sh; sets CPP_F1 / CPP_F1_NORM / CPP_NOTE)
concordance_cpp_run_third_arm "$CPP_VCF" "$TRUTH" "$J" "$OUT" run_eval

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

# Emits extra.json (java/rust summaries + cpp_f1 / cpp_f1_normalized / cpp_version or cpp_skip_reason)
concordance_cpp_write_extra "$OUT" "$CPP_F1" "$CPP_NOTE" "$CPP_F1_NORM"
python3 harness/record_evidence.py --input "$TRUTH" --input "$BED" "$NAME" "$res" "$OUT" -- "$0" "$@"

echo "Germline concordance result: $res"
[ "$res" = EQUAL ]
