#!/usr/bin/env bash
# =============================================================================
# VarDictRust v1.0.4 — Harness Self-Test Runner (harness/run_self_tests.sh)
# Automated Tier 6 Harness Diagnostic Runner for VarDictRust v1.0.4.
#
# Discovers and executes all 25 harness diagnostic scripts (17 shell tests
# and 8 python test suites in harness/tests/), reporting unified status.
# Exits 0 only if all 25 tests pass; exits 1 on any failure.
#
# Adheres to Host Policy A3:
# - Minimal CPU overhead (< 10 seconds total runtime)
# - No persistent artifact generation
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
TESTS_DIR="$SCRIPT_DIR/tests"

# Isolate every self-test from the real pinned store: tests rebuild temp fixtures (edge-case BAMs, VCF-header
# cases) whose cells must never be pinned into scratch/golden-store.
export GOLDEN_STORE_DIR="$(mktemp -d)"
trap 'rm -rf "$GOLDEN_STORE_DIR"' EXIT
unset GOLDEN_MODE

if ! command -v bcftools &>/dev/null; then
  export PATH="/home/palak/miniconda3/envs/g4-tools/bin:$PATH"
fi

# Discover all diagnostic test files
mapfile -t SH_TESTS < <(find "$TESTS_DIR" -maxdepth 1 -type f -name "test_*.sh" | sort)
mapfile -t PY_TESTS < <(find "$TESTS_DIR" -maxdepth 1 -type f -name "test_*.py" | sort)
TOTAL_TESTS=$(( ${#SH_TESTS[@]} + ${#PY_TESTS[@]} ))

echo "================================================================================"
echo " VarDictRust v1.0.4 Harness Diagnostic Self-Test Suite"
echo " Discovered $TOTAL_TESTS Diagnostic Tests in $TESTS_DIR"
echo "================================================================================"

PASSED=0
FAILED=0
declare -a FAILURES=()
CURRENT=0

# 1. Run Shell Diagnostic Tests
for t in "${SH_TESTS[@]}"; do
  CURRENT=$((CURRENT + 1))
  name="$(basename "$t")"
  printf "[%2d/%2d] Running %-32s " "$CURRENT" "$TOTAL_TESTS" "$name..."
  start_ms=$(date +%s%N)

  if out=$(bash "$t" 2>&1); then
    end_ms=$(date +%s%N)
    elapsed_ms=$(( (end_ms - start_ms) / 1000000 ))
    echo "PASS (${elapsed_ms}ms)"
    PASSED=$((PASSED + 1))
  else
    echo "FAIL"
    echo "$out" | sed 's/^/    /'
    FAILED=$((FAILED + 1))
    FAILURES+=("$name")
  fi
done

# 2. Run Python Diagnostic Unit Test Suites
for t in "${PY_TESTS[@]}"; do
  CURRENT=$((CURRENT + 1))
  name="$(basename "$t")"
  printf "[%2d/%2d] Running %-32s " "$CURRENT" "$TOTAL_TESTS" "$name..."
  start_ms=$(date +%s%N)

  if out=$(python3 -m unittest "$t" 2>&1); then
    end_ms=$(date +%s%N)
    elapsed_ms=$(( (end_ms - start_ms) / 1000000 ))
    echo "PASS (${elapsed_ms}ms)"
    PASSED=$((PASSED + 1))
  else
    echo "FAIL"
    echo "$out" | sed 's/^/    /'
    FAILED=$((FAILED + 1))
    FAILURES+=("$name")
  fi
done

echo "================================================================================"
if [[ "$FAILED" -eq 0 ]]; then
  echo " ALL $TOTAL_TESTS HARNESS DIAGNOSTIC TESTS PASSED SUCCESSFULLY"
  echo "================================================================================"
  exit 0
else
  echo " FAILED: $FAILED / $TOTAL_TESTS tests failed: ${FAILURES[*]}"
  echo "================================================================================"
  exit 1
fi
