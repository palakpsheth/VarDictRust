#!/usr/bin/env bash
# AGENTS.md §15 acceptance, items 1-8 verbatim. Exit 0 only if every item is PASS. There is no waiver mechanism.
set -uo pipefail
OUT=${ACCEPT_DIR:-scratch/acceptance}; mkdir -p "$OUT"; : > "$OUT/log.txt"; declare -A R
FP=$(bash harness/engine_fingerprint.sh) || { echo "acceptance: commit first"; exit 2; }
PERF_FP=$(bash harness/engine_fingerprint.sh --perf) || exit 2
ts() { date '+[%H:%M:%S]'; }
mark() {
  local item=$1 ok=$2 desc=$3 elapsed=$((SECONDS - ITEM_T0))
  if [ "$ok" = 1 ]; then
    R[$item]=PASS
    echo "$(ts) [Acceptance $item/8] item $item ($desc): PASS (${elapsed}s)"
  else
    R[$item]="FAIL: $desc"
    echo "$(ts) [Acceptance $item/8] item $item ($desc): FAIL (${elapsed}s)"
  fi
}
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
has_v = isinstance(x.get("cpp_version"), str) and bool(x["cpp_version"].strip())
has_s = isinstance(x.get("cpp_skip_reason"), str) and bool(x["cpp_skip_reason"].strip())
if has_v == has_s:
    sys.exit(f"{p}: extra must have exactly one non-empty cpp_version or cpp_skip_reason (cpp third-arm presence)")
PY
}
# ev_scale: the hours-long scale bench is carried across commits ONLY while the perf-relevant code (engine
# sources, manifests, lockfile) is byte-identical to the run that produced it (extra.perf_fingerprint); a
# perf-relevant change makes the evidence stale with no exception. P>=4 / bar>=2.4 is the pre-registered bar.
ev_scale() { local n=${1:-scale}; python3 - "benchmarks/evidence/$n.json" "$FP" "$PERF_FP" <<'PY'
import json, sys
p, fp, perf = sys.argv[1], sys.argv[2], sys.argv[3]
try: d = json.load(open(p))
except Exception as e: sys.exit(f"{p}: {e}")
x = d.get("extra") or {}
if d.get("fingerprint") != fp and x.get("perf_fingerprint") != perf:
    sys.exit(f"{p}: stale - neither engine fingerprint {fp} nor perf fingerprint {perf} matches (re-run scale_bench.sh under courtesy.sh)")
if d.get("result") not in ("MATCH", "EQUAL", "PASS") or x.get("pass") is not True:
    sys.exit(f"{p}: result={d.get('result')} extra.pass={x.get('pass')}")
if not x.get("certifying"): sys.exit(f"{p}: non-certifying scale run (extra.certifying is not true)")
if x.get("P", 0) < 4 or x.get("bar", 0) < 2.4:
    sys.exit(f"{p}: below the pre-registered bar (needs P>=4 and bar>=2.4): P={x.get('P')} bar={x.get('bar')}")
PY
}
run() { "$@" >> "$OUT/log.txt" 2>&1; }
# 1 Unit
ITEM_T0=$SECONDS; echo "$(ts) [Acceptance 1/8] Starting Unit & Workspace Tests..."
t=1; run cargo test --workspace || t=0; mark 1 $t "cargo test"
# 2 Tier-1: corpus .var parity vs VarDictJava expected + micro .var golden
ITEM_T0=$SECONDS; echo "$(ts) [Acceptance 2/8] Starting Tier-1 Integration Corpus & Smoke..."
status=0; bash harness/run_integration.sh > "$OUT/integration.txt" 2>&1 || status=$?
t=0; [ "$status" -eq 0 ] \
  && grep -qx 'PASS: 91  XFAIL: 1  FAIL: 0  SKIP: 0  XPASS: 0' "$OUT/integration.txt" \
  && [ "$(grep -c '^XFAIL: ' "$OUT/integration.txt")" -eq 1 ] \
  && grep -qxF 'XFAIL: Simple;hard_clip_case.fa;hard_clip_next_to_del_test2.bam;test;6434-6564;-f 0.0 -p -r 1.txt (UB-4 upstream fixture defect, verified against live-oracle signature)' "$OUT/integration.txt" \
  && run cargo build -q --release && run bash harness/smoke.sh && t=1
mark 2 $t "integration/smoke"
# 3 Tier-2: micro + every committed slice + live extended coverage + external datasets
ITEM_T0=$SECONDS; echo "$(ts) [Acceptance 3/8] Starting Tier-2 Slices, Extended Modalities & Dataset Evidence..."
t=1
SLICE_KEEP_TMP=$(mktemp -d)
trap 'rm -rf "$SLICE_KEEP_TMP"' EXIT
export SLICE_KEEP="${SLICE_KEEP:-$SLICE_KEEP_TMP}"
export VARDICT_BIN="${VARDICT_BIN:-$(pwd)/target/release/vardictrust}"
{
  ( set +e; RUST_THREADS=1 MIN_VARIANTS=1 bash harness/compare.sh fixtures/micro/micro.ref.fa fixtures/micro/micro.bam fixtures/micro/micro.bed > "$SLICE_KEEP_TMP/micro.log" 2>&1 ); rc=$?
  printf '%d\n' "$rc" > "$SLICE_KEEP_TMP/micro.rc.tmp" && mv "$SLICE_KEEP_TMP/micro.rc.tmp" "$SLICE_KEEP_TMP/micro.rc"
} &
for id in $(awk -F'\t' '!/^#/{print $1}' harness/slices.tsv); do
  while [ "$(jobs -rp | wc -l)" -ge 8 ]; do
    wait -n || true
  done
  {
    ( set +e; RUST_THREADS=1 bash harness/run_slice.sh "$id" > "$SLICE_KEEP_TMP/$id.log" 2>&1 ); rc=$?
    printf '%d\n' "$rc" > "$SLICE_KEEP_TMP/$id.rc.tmp" && mv "$SLICE_KEEP_TMP/$id.rc.tmp" "$SLICE_KEEP_TMP/$id.rc"
  } &
done
wait
cat "$SLICE_KEEP_TMP/micro.log" >> "$OUT/log.txt" 2>&1 || true
[ "$(cat "$SLICE_KEEP_TMP/micro.rc" 2>/dev/null || echo 1)" -eq 0 ] || t=0
for id in $(awk -F'\t' '!/^#/{print $1}' harness/slices.tsv); do
  cat "$SLICE_KEEP_TMP/$id.log" >> "$OUT/log.txt" 2>&1 || true
  if [ "$(cat "$SLICE_KEEP_TMP/$id.rc" 2>/dev/null || echo 1)" -ne 0 ]; then
    t=0
    echo "slice $id failed" >> "$OUT/log.txt"
  fi
done
if ev extended-coverage-tier2 >> "$OUT/log.txt" 2>&1; then
  echo "acceptance: extended-coverage-tier2 already certified for $FP" >> "$OUT/log.txt"
else
  run bash harness/run_extended_coverage.sh || { t=0; echo "run_extended_coverage failed" >> "$OUT/log.txt"; }
fi
if [ ! -d scratch/g4/multimodal ]; then
  echo "acceptance: scratch/g4/multimodal is missing; live multimodal suite required" >> "$OUT/log.txt"
  t=0
elif ev multimodal-slices-tier2 >> "$OUT/log.txt" 2>&1; then
  echo "acceptance: multimodal-slices-tier2 already certified for $FP" >> "$OUT/log.txt"
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
# pinned golden store integrity (issue #31): every live artifact re-hashes to its manifest sha256
run bash harness/golden_store.sh verify --quiet || { t=0; echo "golden store verify failed" >> "$OUT/log.txt"; }
mark 3 $t "Tier-2"
# 4 Robustness: crash corpus (#341 family, PR #402, #187) + CRAM twin
ITEM_T0=$SECONDS; echo "$(ts) [Acceptance 4/8] Starting Robustness & CRAM Twin..."
t=1; run cargo test --release -p vardictrust-core --test engine_behavior || t=0; run cargo test --release -p vardictrust --test robustness || t=0; run ev hg008-cram-twin || t=0
if [ "${CPP_SKIP:-0}" != "1" ] && [ "${CPP_SKIP:-0}" != "true" ]; then run ev_cpp hg008-cram-twin || t=0; fi
mark 4 $t "robustness/CRAM"
# 5 Scale + Perf Budget (thresholds pre-registered in T21)
ITEM_T0=$SECONDS; echo "$(ts) [Acceptance 5/8] Starting Scale & Performance Budget..."
t=1; run ev_scale || t=0; run ev perf-budget || t=0; mark 5 $t "scale/perf-budget"

# 6 Truth set: rtg exact java==rust equality on all 5 certified cohorts
ITEM_T0=$SECONDS; echo "$(ts) [Acceptance 6/8] Starting Truth Concordance (5 Cohorts)..."
t=1; for e in hg008-truth seqc2-truth hg002-truth hg005-truth na24385-b37-truth; do run ev "$e" || t=0; done; mark 6 $t "truth"
# 7 Hygiene
ITEM_T0=$SECONDS; echo "$(ts) [Acceptance 7/8] Starting Hygiene (fmt, clippy, audit, deny, gitleaks)..."
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
ITEM_T0=$SECONDS; echo "$(ts) [Acceptance 8/8] Starting Distribution & musl Check..."
t=1; v=$(target/release/vardictrust --version 2>&1)
for s in "VarDictJava v1.8.3" 0886261 009e017 Lai gkw227 AstraZeneca-NGS "AI disclosure"; do grep -qF "$s" <<<"$v" || t=0; done
run ev musl-ci || t=0; mark 8 $t "distribution"
rc=0; for i in 1 2 3 4 5 6 7 8; do echo "item $i: ${R[$i]}"; [[ "${R[$i]}" == FAIL* ]] && rc=1; done
echo "fingerprint $FP (log: $OUT/log.txt)"; exit $rc
