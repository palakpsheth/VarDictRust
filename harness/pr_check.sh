#!/usr/bin/env bash
# harness/pr_check.sh — Standard PR Regression & Performance Gate
# Enforces full regression testing, code hygiene, and performance metric checks on all PRs.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
COMMON_GIT_DIR="$(git -C "$ROOT_DIR" rev-parse --path-format=absolute --git-common-dir 2>/dev/null || true)"
if [[ -n "$COMMON_GIT_DIR" && -d "$(dirname "$COMMON_GIT_DIR")/VarDictJava" ]]; then
  DATA_ROOT="$(cd "$(dirname "$COMMON_GIT_DIR")" && pwd)"
else
  DATA_ROOT="$ROOT_DIR"
fi

SLICE_KEEP_TMP=$(mktemp -d)
trap 'rm -rf "$SLICE_KEEP_TMP"' EXIT
export SLICE_KEEP="${SLICE_KEEP:-$SLICE_KEEP_TMP}"
export VARDICT_BIN="${VARDICT_BIN:-$ROOT_DIR/target/release/vardictrust}"
export COMPARE_NO_TIME=1
ts() { date '+[%H:%M:%S]'; }
PR_T0=$SECONDS

echo "=== $(ts) VarDictRust PR Regression & Performance Gate ==="

# Stage 1: Hygiene (fmt, clippy, unit + integration + robustness tests, release binary)
st_t0=$SECONDS
echo "$(ts) [Stage 1/6] Running workspace hygiene and unit tests..."
cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --workspace --bins --lib --tests
cargo build -q --release
echo "$(ts) [Stage 1/6] PASS ($((SECONDS - st_t0))s)"

# Stage 2: Crash robustness test suite (executed in Stage 1 workspace test pass)
echo "$(ts) [Stage 2/6] Crash robustness test suite verified in workspace test pass (0s)."

# Stage 3: Tier-1 Upstream Integration Corpus (92 fixtures: 91 PASS / 1 XFAIL)
st_t0=$SECONDS
echo "$(ts) [Stage 3/6] Running Tier-1 integration regression corpus..."
mkdir -p scratch
status=0
exec 3>&1
INTEGRATION_PROGRESS_FD=3 bash harness/run_integration.sh > scratch/pr_integration.log 2>&1 || status=$?
exec 3>&-
[ "$status" -eq 0 ] \
  && grep -qx 'PASS: 91  XFAIL: 1  FAIL: 0  SKIP: 0  XPASS: 0' scratch/pr_integration.log \
  && [ "$(grep -c '^XFAIL: ' scratch/pr_integration.log)" -eq 1 ] \
  && grep -qxF 'XFAIL: Simple;hard_clip_case.fa;hard_clip_next_to_del_test2.bam;test;6434-6564;-f 0.0 -p -r 1.txt (UB-4 upstream fixture defect, verified against live-oracle signature)' scratch/pr_integration.log \
  || {
    echo "PR Regression Failure: Tier-1 corpus did not match expected 91 PASS / 1 XFAIL signature (exit status: $status)!"
    cat scratch/pr_integration.log
    exit 1
  }
echo "$(ts) [Stage 3/6] Tier-1 integration corpus PASS: 91 PASS / 1 XFAIL ($((SECONDS - st_t0))s)"

# Stage 4: Tier-2 Dual-Engine Differential Parity (Pinned Oracle vs Rust)
st_t0=$SECONDS
echo "$(ts) [Stage 4/6] Running Tier-2 differential parity suites against oracle..."
SLICE_IDS=(micro s1a s1b tiny s3 s4 s5 s6 s7 s8 s9)
(
  {
    s_t0=$(date +%s%3N)
    ( set +e; KEEP_DIR="$SLICE_KEEP_TMP/micro" RUST_THREADS=1 MIN_VARIANTS=1 bash harness/compare.sh fixtures/micro/micro.ref.fa fixtures/micro/micro.bam fixtures/micro/micro.bed > "$SLICE_KEEP_TMP/micro.log" 2>&1 ); rc=$?
    printf '%d\n' "$(( $(date +%s%3N) - s_t0 ))" > "$SLICE_KEEP_TMP/micro.ms"
    printf '%d\n' "$rc" > "$SLICE_KEEP_TMP/micro.rc.tmp" && mv "$SLICE_KEEP_TMP/micro.rc.tmp" "$SLICE_KEEP_TMP/micro.rc"
  } &
  for id in s1a s1b tiny s3 s4 s5 s6 s7 s8 s9; do
    while [ "$(jobs -rp | wc -l)" -ge 8 ]; do
      wait -n 2>/dev/null || sleep 0.02
    done
    {
      s_t0=$(date +%s%3N)
      ( set +e; RUST_THREADS=1 bash harness/run_slice.sh "$id" > "$SLICE_KEEP_TMP/$id.log" 2>&1 ); rc=$?
      printf '%d\n' "$(( $(date +%s%3N) - s_t0 ))" > "$SLICE_KEEP_TMP/$id.ms"
      printf '%d\n' "$rc" > "$SLICE_KEEP_TMP/$id.rc.tmp" && mv "$SLICE_KEEP_TMP/$id.rc.tmp" "$SLICE_KEEP_TMP/$id.rc"
    } &
  done
  wait
) &
POOL_PID=$!

for idx in "${!SLICE_IDS[@]}"; do
  id="${SLICE_IDS[$idx]}"
  k=$((idx + 1))
  while [ ! -f "$SLICE_KEEP_TMP/$id.rc" ] && kill -0 "$POOL_PID" 2>/dev/null; do
    sleep 0.02
  done
  if [ "$(cat "$SLICE_KEEP_TMP/$id.rc" 2>/dev/null || echo 1)" -ne 0 ]; then
    echo "FAIL slice $id:"
    cat "$SLICE_KEEP_TMP/$id.log"
    wait "$POOL_PID" 2>/dev/null || true
    exit 1
  fi
  s_vars=$(grep -vc '^#' "$SLICE_KEEP_TMP/$id/rust.vcf" 2>/dev/null || echo 0)
  s_ms=$(cat "$SLICE_KEEP_TMP/$id.ms" 2>/dev/null || echo 0)
  echo "$(ts) [slice $k/${#SLICE_IDS[@]}] $id: MATCH ($s_vars variants, ${s_ms}ms)"
done
wait "$POOL_PID" 2>/dev/null || true
echo "$(ts) All 10 Tier-2 standard slices (s1a..s9) + micro passed."

echo "$(ts) Running 17-modality extended coverage suite (--no-record)..."
bash harness/run_extended_coverage.sh --no-record

if [ -d "$DATA_ROOT/scratch/g4/multimodal" ]; then
  echo "$(ts) Running 18-slice multimodal suite (--no-record)..."
  bash harness/run_multimodal_slices.sh --no-record
fi
echo "$(ts) [Stage 4/6] PASS ($((SECONDS - st_t0))s)"

# Stage 5: Performance & Memory Benchmark Regression Check
st_t0=$SECONDS
echo "$(ts) [Stage 5/6] Checking performance budgets on benchmark targets..."
if [ -f benchmarks/evidence/perf-budget.json ]; then
  python3 harness/perf_budget.py hg008-tier2 seqc2-exome-tier2 hg002-tier2 amplicon-tier2
fi
echo "$(ts) [Stage 5/6] PASS ($((SECONDS - st_t0))s)"

# Stage 6: Security Hygiene (Gitleaks scan)
st_t0=$SECONDS
echo "$(ts) [Stage 6/6] Scanning commits for secrets..."
if command -v gitleaks >/dev/null 2>&1; then
  gitleaks detect --source . -v
elif [ -x /home/palak/miniconda3/envs/g4-tools/bin/gitleaks ]; then
  /home/palak/miniconda3/envs/g4-tools/bin/gitleaks detect --source . -v
else
  echo "ERROR: gitleaks scanner not found. Install gitleaks to satisfy PR security gate." >&2
  exit 1
fi
echo "$(ts) [Stage 6/6] PASS ($((SECONDS - st_t0))s)"

echo "=== $(ts) ALL PR REGRESSION & PERFORMANCE GATES PASSED ($((SECONDS - PR_T0))s) ==="
