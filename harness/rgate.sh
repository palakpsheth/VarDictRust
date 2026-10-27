#!/usr/bin/env bash
# Regression gate (plan rev2 §R-GATE). Exit 0 = PASS. Never edit this to make a change pass.
set -uo pipefail
OUT=${RGATE_DIR:-scratch/rgate}; mkdir -p "$OUT"; fail=0
SLICE_KEEP_TMP=$(mktemp -d)
trap 'rm -rf "$SLICE_KEEP_TMP"' EXIT
export SLICE_KEEP="${SLICE_KEEP:-$SLICE_KEEP_TMP}"
export VARDICT_BIN="${VARDICT_BIN:-$(pwd)/target/release/vardictrust}"
export COMPARE_NO_TIME=1
ts() { date '+[%H:%M:%S]'; }
RGATE_T0=$SECONDS

st_t0=$SECONDS
echo "$(ts) [rgate 1/4] == cargo fmt/clippy/test"
{ cargo fmt --all -- --check && cargo clippy --all-targets --all-features -- -D warnings \
    && cargo test --workspace --bins --lib --tests && cargo build -q --release; } > "$OUT/cargo.txt" 2>&1 \
  && echo "$(ts) [rgate 1/4] PASS ($((SECONDS - st_t0))s)" \
  || { echo "FAIL cargo ($OUT/cargo.txt)"; fail=1; }

st_t0=$SECONDS
echo "$(ts) [rgate 2/4] == integration corpus"
status=0
exec 3>&1
INTEGRATION_PROGRESS_FD=3 bash harness/run_integration.sh > "$OUT/integration.txt" 2>&1 || status=$?
exec 3>&-
if [ "$status" -eq 0 ] \
   && grep -qx 'PASS: 91  XFAIL: 1  FAIL: 0  SKIP: 0  XPASS: 0' "$OUT/integration.txt" \
   && [ "$(grep -c '^XFAIL: ' "$OUT/integration.txt")" -eq 1 ] \
   && grep -qxF 'XFAIL: Simple;hard_clip_case.fa;hard_clip_next_to_del_test2.bam;test;6434-6564;-f 0.0 -p -r 1.txt (UB-4 upstream fixture defect, verified against live-oracle signature)' "$OUT/integration.txt"
then
  echo "$(ts) [rgate 2/4] PASS: 91 / XFAIL: 1 ($((SECONDS - st_t0))s)"
else
  echo "FAIL integration ($OUT/integration.txt)"; fail=1
fi

# Launch micro and standard slices concurrently (1 thread per Rust invocation, <= 8 workers total)
st_t0=$SECONDS
(
  {
    s_t0=$(date +%s%3N)
    ( set +e; KEEP_DIR="$SLICE_KEEP_TMP/micro" RUST_THREADS=1 MIN_VARIANTS=1 bash harness/compare.sh fixtures/micro/micro.ref.fa fixtures/micro/micro.bam \
      fixtures/micro/micro.bed > "$OUT/micro.txt" 2>&1 ); rc=$?
    printf '%d\n' "$(( $(date +%s%3N) - s_t0 ))" > "$SLICE_KEEP_TMP/micro.ms"
    printf '%d\n' "$rc" > "$SLICE_KEEP_TMP/micro.rc.tmp" && mv "$SLICE_KEEP_TMP/micro.rc.tmp" "$SLICE_KEEP_TMP/micro.rc"
  } &
  for id in ${RGATE_SLICES:-s1a s1b s5 s6 s7}; do
    while [ "$(jobs -rp | wc -l)" -ge 8 ]; do
      wait -n 2>/dev/null || sleep 0.02
    done
    {
      s_t0=$(date +%s%3N)
      ( set +e; RUST_THREADS=1 bash harness/run_slice.sh "$id" > "$OUT/$id.txt" 2>&1 ); rc=$?
      printf '%d\n' "$(( $(date +%s%3N) - s_t0 ))" > "$SLICE_KEEP_TMP/$id.ms"
      printf '%d\n' "$rc" > "$SLICE_KEEP_TMP/$id.rc.tmp" && mv "$SLICE_KEEP_TMP/$id.rc.tmp" "$SLICE_KEEP_TMP/$id.rc"
    } &
  done
  wait
) &
POOL_PID=$!

while [ ! -f "$SLICE_KEEP_TMP/micro.rc" ] && kill -0 "$POOL_PID" 2>/dev/null; do sleep 0.02; done
echo "$(ts) [rgate 3/4] == Tier-2 micro"
if [ "$(cat "$SLICE_KEEP_TMP/micro.rc" 2>/dev/null || echo 1)" -ne 0 ]; then
  echo "FAIL micro ($OUT/micro.txt)"; fail=1
else
  m_vars=$(grep -vc '^#' "$SLICE_KEEP_TMP/micro/rust.vcf" 2>/dev/null || echo 1)
  m_ms=$(cat "$SLICE_KEEP_TMP/micro.ms" 2>/dev/null || echo 0)
  echo "$(ts) [slice micro] MATCH ($m_vars variants, ${m_ms}ms)"
fi
for id in ${RGATE_SLICES:-s1a s1b s5 s6 s7}; do
  while [ ! -f "$SLICE_KEEP_TMP/$id.rc" ] && kill -0 "$POOL_PID" 2>/dev/null; do sleep 0.02; done
  echo "$(ts) [rgate 4/4] == slice $id"
  if [ "$(cat "$SLICE_KEEP_TMP/$id.rc" 2>/dev/null || echo 1)" -ne 0 ]; then
    echo "FAIL $id ($OUT/$id.txt)"; fail=1
  else
    s_vars=$(grep -vc '^#' "$SLICE_KEEP_TMP/$id/rust.vcf" 2>/dev/null || echo 0)
    s_ms=$(cat "$SLICE_KEEP_TMP/$id.ms" 2>/dev/null || echo 0)
    echo "$(ts) [slice $id] MATCH ($s_vars variants, ${s_ms}ms)"
  fi
done
wait "$POOL_PID" 2>/dev/null || true
[ "$fail" -eq 0 ] && echo "$(ts) RGATE PASS ($((SECONDS - RGATE_T0))s)" || echo "$(ts) RGATE FAIL ($((SECONDS - RGATE_T0))s)"
exit "$fail"
