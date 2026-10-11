#!/usr/bin/env bash
# Regression gate (plan rev2 §R-GATE). Exit 0 = PASS. Never edit this to make a change pass.
set -uo pipefail
OUT=${RGATE_DIR:-scratch/rgate}; mkdir -p "$OUT"; fail=0
echo "== cargo fmt/clippy/test"
{ cargo fmt --all -- --check && cargo clippy --all-targets --all-features -- -D warnings \
    && cargo test --workspace; } > "$OUT/cargo.txt" 2>&1 || { echo "FAIL cargo ($OUT/cargo.txt)"; fail=1; }
echo "== integration corpus"
status=0; bash harness/run_integration.sh > "$OUT/integration.txt" 2>&1 || status=$?
if [ "$status" -eq 0 ] \
   && grep -qx 'PASS: 91  XFAIL: 1  FAIL: 0  SKIP: 0  XPASS: 0' "$OUT/integration.txt" \
   && [ "$(grep -c '^XFAIL: ' "$OUT/integration.txt")" -eq 1 ] \
   && grep -qxF 'XFAIL: Simple;hard_clip_case.fa;hard_clip_next_to_del_test2.bam;test;6434-6564;-f 0.0 -p -r 1.txt (UB-4 upstream fixture defect, verified against live-oracle signature)' "$OUT/integration.txt"
then :; else echo "FAIL integration ($OUT/integration.txt)"; fail=1; fi
echo "== Tier-2 micro"
MIN_VARIANTS=1 bash harness/compare.sh fixtures/micro/micro.ref.fa fixtures/micro/micro.bam \
  fixtures/micro/micro.bed > "$OUT/micro.txt" 2>&1 || { echo "FAIL micro ($OUT/micro.txt)"; fail=1; }
for id in ${RGATE_SLICES:-s1a s1b s5 s6 s7}; do
  echo "== slice $id"
  bash harness/run_slice.sh "$id" > "$OUT/$id.txt" 2>&1 || { echo "FAIL $id ($OUT/$id.txt)"; fail=1; }
done
[ "$fail" -eq 0 ] && echo "RGATE PASS" || echo "RGATE FAIL"
exit "$fail"
