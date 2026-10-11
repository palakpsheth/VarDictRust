#!/usr/bin/env bash
# PLAN.md §15 acceptance, items 1-8 verbatim. Exit 0 only if every item is PASS or user-WAIVED.
set -uo pipefail
OUT=${ACCEPT_DIR:-scratch/acceptance}; mkdir -p "$OUT"; : > "$OUT/log.txt"; declare -A R
FP=$(bash harness/engine_fingerprint.sh) || { echo "acceptance: commit first"; exit 2; }
waived() { grep -qE "^- item $1:" docs/release/waivers.md 2>/dev/null; }
mark() { if [ "$2" = 1 ]; then R[$1]=PASS; elif waived "$1"; then R[$1]=WAIVED; else R[$1]="FAIL: $3"; fi; }
ev() { python3 - "benchmarks/evidence/$1.json" "$FP" <<'PY'
import json, sys
p, fp = sys.argv[1], sys.argv[2]
try: d = json.load(open(p))
except Exception as e: sys.exit(f"{p}: {e}")
if d.get("fingerprint") != fp: sys.exit(f"{p}: stale fingerprint {d.get('fingerprint')} != {fp}")
if d.get("result") not in ("MATCH", "EQUAL", "PASS"): sys.exit(f"{p}: result={d.get('result')}")
PY
}
run() { "$@" >> "$OUT/log.txt" 2>&1; }
# 1 Unit
t=1; run cargo test --workspace || t=0; mark 1 $t "cargo test"
# 2 Tier-1: corpus .var parity vs VarDictJava expected + micro .var golden
bash harness/run_integration.sh > "$OUT/integration.txt" 2>&1
t=0; grep -qx 'PASS: 91  FAIL: 1  SKIP: 0' "$OUT/integration.txt" \
  && [ "$(grep -c '^MISMATCH: ' "$OUT/integration.txt")" -eq 1 ] \
  && grep -qxF 'MISMATCH: Simple;hard_clip_case.fa;hard_clip_next_to_del_test2.bam;test;6434-6564;-f 0.0 -p -r 1.txt' "$OUT/integration.txt" \
  && run cargo build -q --release && run bash harness/smoke.sh && t=1
mark 2 $t "integration/smoke"
# 3 Tier-2: micro + every committed slice + external datasets
t=1; MIN_VARIANTS=1 run bash harness/compare.sh fixtures/micro/micro.ref.fa fixtures/micro/micro.bam fixtures/micro/micro.bed || t=0
for id in $(awk -F'\t' '!/^#/{print $1}' harness/slices.tsv); do run bash harness/run_slice.sh "$id" || { t=0; echo "slice $id failed" >> "$OUT/log.txt"; }; done
for e in hg008-tier2 seqc2-slices-tier2 seqc2-exome-tier2 hg002-tier2 amplicon-tier2 extended-coverage-tier2 multimodal-slices-tier2; do run ev "$e" || t=0; done
mark 3 $t "Tier-2"
# 4 Robustness: crash corpus (#341 family, PR #402, #187) + CRAM twin
t=1; run cargo test -p vardictrust --test robustness || t=0; run ev hg008-cram-twin || t=0; mark 4 $t "robustness/CRAM"
# 5 Scale + Perf Budget (thresholds pre-registered in T21)
t=1; run ev scale || t=0; run ev perf-budget || t=0; mark 5 $t "scale/perf-budget"

# 6 Truth set: rtg exact java==rust equality on both datasets (som.py/truvari via waiver)
t=1; for e in hg008-truth seqc2-truth hg002-truth; do run ev "$e" || t=0; done; mark 6 $t "truth"
# 7 Hygiene
t=1; run cargo fmt --all -- --check || t=0; run cargo clippy --all-targets --all-features -- -D warnings || t=0
run cargo audit || t=0; run cargo deny check || t=0; mark 7 $t "hygiene"
# 8 Distribution: musl CI evidence + --version content
t=1; v=$(target/release/vardictrust --version 2>&1)
for s in "VarDictJava v1.8.3" 0886261 009e017 Lai gkw227 AstraZeneca-NGS "AI disclosure"; do grep -qF "$s" <<<"$v" || t=0; done
run ev musl-ci || t=0; mark 8 $t "distribution"
rc=0; for i in 1 2 3 4 5 6 7 8; do echo "item $i: ${R[$i]}"; [[ "${R[$i]}" == FAIL* ]] && rc=1; done
echo "fingerprint $FP (log: $OUT/log.txt)"; exit $rc
