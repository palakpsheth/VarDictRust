#!/usr/bin/env bash
# harness/pr_check.sh — Standard PR Regression & Performance Gate
# Enforces full regression testing, code hygiene, and performance metric checks on all PRs.
set -euo pipefail

echo "=== VarDictRust PR Regression & Performance Gate ==="

# Stage 1: Hygiene (fmt, clippy, unit tests)
echo "[1/6] Running workspace hygiene and unit tests..."
cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --workspace

# Stage 2: Crash robustness test suite (Issue #341 family, PR #402, missing BAI #187)
echo "[2/6] Running crash robustness test suite..."
cargo test -p vardictrust --test robustness

# Stage 3: Tier-1 Upstream Integration Corpus (92 fixtures: 91 PASS / 1 XFAIL)
echo "[3/6] Running Tier-1 integration regression corpus..."
mkdir -p scratch
status=0
bash harness/run_integration.sh > scratch/pr_integration.log 2>&1 || status=$?
[ "$status" -eq 0 ] \
  && grep -qx 'PASS: 91  XFAIL: 1  FAIL: 0  SKIP: 0  XPASS: 0' scratch/pr_integration.log \
  && [ "$(grep -c '^XFAIL: ' scratch/pr_integration.log)" -eq 1 ] \
  && grep -qxF 'XFAIL: Simple;hard_clip_case.fa;hard_clip_next_to_del_test2.bam;test;6434-6564;-f 0.0 -p -r 1.txt (UB-4 upstream fixture defect, verified against live-oracle signature)' scratch/pr_integration.log \
  || {
    echo "PR Regression Failure: Tier-1 corpus did not match expected 91 PASS / 1 XFAIL signature (exit status: $status)!"
    cat scratch/pr_integration.log
    exit 1
  }
echo "Tier-1 integration corpus passed (91 PASS / 1 XFAIL expected fixture defect)."

# Stage 4: Tier-2 Dual-Engine Differential Parity (Live Oracle vs Rust)
echo "[4/6] Running Tier-2 differential parity suites against live oracle..."
cargo build -q --release
MIN_VARIANTS=1 bash harness/compare.sh fixtures/micro/micro.ref.fa fixtures/micro/micro.bam fixtures/micro/micro.bed
echo "Running live 17-modality extended coverage suite..."
bash harness/run_extended_coverage.sh

# Stage 5: Performance & Memory Benchmark Regression Check
echo "[5/6] Checking performance budgets on benchmark targets..."
if [ -f benchmarks/evidence/perf-budget.json ]; then
  python3 harness/perf_budget.py hg008-tier2 seqc2-exome-tier2 hg002-tier2 amplicon-tier2
fi

# Stage 6: Security Hygiene (Gitleaks scan)
echo "[6/6] Scanning commits for secrets..."
if command -v gitleaks >/dev/null 2>&1; then
  gitleaks detect --source . -v
elif [ -x /home/palak/miniconda3/envs/g4-tools/bin/gitleaks ]; then
  /home/palak/miniconda3/envs/g4-tools/bin/gitleaks detect --source . -v
else
  echo "ERROR: gitleaks scanner not found. Install gitleaks to satisfy PR security gate." >&2
  exit 1
fi

echo "=== ALL PR REGRESSION & PERFORMANCE GATES PASSED ==="
