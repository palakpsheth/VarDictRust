#!/usr/bin/env bash
# AGENTS.md §15 acceptance, items 1-8 verbatim. Exit 0 only if every item is PASS or user-WAIVED.
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
ev_cpp() { python3 - "benchmarks/evidence/$1.json" <<'PY'
import json, sys
p = sys.argv[1]
try: d = json.load(open(p))
except Exception as e: sys.exit(f"{p}: {e}")
x = d.get("extra") or {}
if "cpp_version" not in x and "cpp_skip_reason" not in x:
    sys.exit(f"{p}: missing extra.cpp_version/cpp_skip_reason (cpp third-arm presence)")
PY
}
run() { "$@" >> "$OUT/log.txt" 2>&1; }
# 1 Unit
t=1; run cargo test --workspace || t=0; mark 1 $t "cargo test"
# 2 Tier-1: corpus .var parity vs VarDictJava expected + micro .var golden
status=0; bash harness/run_integration.sh > "$OUT/integration.txt" 2>&1 || status=$?
t=0; [ "$status" -eq 0 ] \
  && grep -qx 'PASS: 91  XFAIL: 1  FAIL: 0  SKIP: 0  XPASS: 0' "$OUT/integration.txt" \
  && [ "$(grep -c '^XFAIL: ' "$OUT/integration.txt")" -eq 1 ] \
  && grep -qxF 'XFAIL: Simple;hard_clip_case.fa;hard_clip_next_to_del_test2.bam;test;6434-6564;-f 0.0 -p -r 1.txt (UB-4 upstream fixture defect, verified against live-oracle signature)' "$OUT/integration.txt" \
  && run cargo build -q --release && run bash harness/smoke.sh && t=1
mark 2 $t "integration/smoke"
# 3 Tier-2: micro + every committed slice + live extended coverage + external datasets
t=1; MIN_VARIANTS=1 run bash harness/compare.sh fixtures/micro/micro.ref.fa fixtures/micro/micro.bam fixtures/micro/micro.bed || t=0
for id in $(awk -F'\t' '!/^#/{print $1}' harness/slices.tsv); do run bash harness/run_slice.sh "$id" || { t=0; echo "slice $id failed" >> "$OUT/log.txt"; }; done
run bash harness/run_extended_coverage.sh || { t=0; echo "run_extended_coverage failed" >> "$OUT/log.txt"; }
if [ ! -d scratch/g4/multimodal ]; then
  echo "acceptance: scratch/g4/multimodal is missing; live multimodal suite required" >> "$OUT/log.txt"
  t=0
else
  run bash harness/run_multimodal_slices.sh || { t=0; echo "run_multimodal_slices failed" >> "$OUT/log.txt"; }
fi
for e in hg008-tier2 seqc2-slices-tier2 seqc2-exome-tier2 hg002-tier2 amplicon-tier2 extended-coverage-tier2 multimodal-slices-tier2; do run ev "$e" || t=0; done
# cpp third-arm presence (advisory arm, loud gate): every Tier-2 evidence extra
# carries cpp_version or cpp_skip_reason; the java==rust verdict above is untouched.
# CPP_SKIP=1 escapes this check (hosts without the vardictcpp pin).
if [ "${CPP_SKIP:-0}" = "1" ] || [ "${CPP_SKIP:-0}" = "true" ]; then
  echo "acceptance: item 3 cpp-presence check skipped via CPP_SKIP=1" >> "$OUT/log.txt"
else
  for e in hg008-tier2 seqc2-slices-tier2 seqc2-exome-tier2 hg002-tier2 amplicon-tier2 extended-coverage-tier2 multimodal-slices-tier2; do run ev_cpp "$e" || { t=0; echo "cpp presence $e failed" >> "$OUT/log.txt"; }; done
fi
mark 3 $t "Tier-2"
# 4 Robustness: crash corpus (#341 family, PR #402, #187) + CRAM twin
t=1; run cargo test -p vardictrust --test robustness || t=0; run ev hg008-cram-twin || t=0; mark 4 $t "robustness/CRAM"
# 5 Scale + Perf Budget (thresholds pre-registered in T21)
t=1; run ev scale || t=0; run ev perf-budget || t=0; mark 5 $t "scale/perf-budget"

# 6 Truth set: rtg exact java==rust equality on all 5 certified cohorts (som.py/truvari via waiver)
t=1; for e in hg008-truth seqc2-truth hg002-truth hg005-truth na24385-b37-truth; do run ev "$e" || t=0; done; mark 6 $t "truth"
# 7 Hygiene
t=1; run cargo fmt --all -- --check || t=0; run cargo clippy --all-targets --all-features -- -D warnings || t=0
run cargo audit || t=0; run cargo deny check || t=0
if command -v gitleaks >/dev/null 2>&1; then
  run gitleaks detect --source . -v || t=0
elif [ -x /home/palak/miniconda3/envs/g4-tools/bin/gitleaks ]; then
  run /home/palak/miniconda3/envs/g4-tools/bin/gitleaks detect --source . -v || t=0
else
  echo "acceptance: gitleaks scanner not found" >> "$OUT/log.txt"; t=0
fi
mark 7 $t "hygiene"
# 8 Distribution: musl CI evidence + --version content
t=1; v=$(target/release/vardictrust --version 2>&1)
for s in "VarDictJava v1.8.3" 0886261 009e017 Lai gkw227 AstraZeneca-NGS "AI disclosure"; do grep -qF "$s" <<<"$v" || t=0; done
run ev musl-ci || t=0; mark 8 $t "distribution"
rc=0; for i in 1 2 3 4 5 6 7 8; do echo "item $i: ${R[$i]}"; [[ "${R[$i]}" == FAIL* ]] && rc=1; done
echo "fingerprint $FP (log: $OUT/log.txt)"; exit $rc
