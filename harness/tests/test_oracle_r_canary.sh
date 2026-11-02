#!/usr/bin/env bash
set -uo pipefail
# Test for Linux x86_64 R 80-bit long double Fisher MLE canaries in vardict-oracle.
fail=0

HOST_OS="$(uname -s)"
HOST_ARCH="$(uname -m)"
if [ "$HOST_OS" != "Linux" ] || [ "$HOST_ARCH" != "x86_64" ]; then
  echo "SKIP: test_oracle_r_canary requires Linux/x86_64 host (saw ${HOST_OS}/${HOST_ARCH})"
  echo "PASS test_oracle_r_canary"
  exit 0
fi

RSCRIPT_CMD=()
if command -v conda >/dev/null 2>&1; then
  CONDA_BASE="$(conda info --base 2>/dev/null || true)"
  if [ -n "$CONDA_BASE" ] && [ -x "$CONDA_BASE/envs/vardict-oracle/bin/Rscript" ]; then
    RSCRIPT_CMD=("$CONDA_BASE/envs/vardict-oracle/bin/Rscript")
  elif conda run -n vardict-oracle Rscript --version >/dev/null 2>&1; then
    RSCRIPT_CMD=(conda run -n vardict-oracle Rscript)
  fi
fi

if [ "${#RSCRIPT_CMD[@]}" -eq 0 ]; then
  echo "SKIP: Rscript not available in conda env vardict-oracle"
  echo "PASS test_oracle_r_canary"
  exit 0
fi

out=$("${RSCRIPT_CMD[@]}" -e '
cat("sizeof.longdouble=", .Machine$sizeof.longdouble, "\n", sep="")
cat("c(9,6,5,1)=", sprintf("%.5f", round(fisher.test(matrix(c(9,6,5,1), nrow=2))$estimate, 5)), "\n", sep="")
cat("c(27,17,3,8)=", sprintf("%.5f", round(fisher.test(matrix(c(27,17,3,8), nrow=2))$estimate, 5)), "\n", sep="")
cat("c(10,8,1,3)=", sprintf("%.5f", round(fisher.test(matrix(c(10,8,1,3), nrow=2))$estimate, 5)), "\n", sep="")
cat("c(2,150,7,183)=", sprintf("%.5f", round(fisher.test(matrix(c(2,150,7,183), nrow=2))$estimate, 5)), "\n", sep="")
' 2>&1) || { echo "FAIL: Rscript execution failed: $out"; exit 1; }

grep -qx 'sizeof.longdouble=16' <<<"$out" || { echo "FAIL: expected sizeof.longdouble=16, got: $out"; fail=1; }
grep -qx 'c(9,6,5,1)=0.31603' <<<"$out" || { echo "FAIL: expected c(9,6,5,1)=0.31603, got: $out"; fail=1; }
grep -qx 'c(27,17,3,8)=4.12268' <<<"$out" || { echo "FAIL: expected c(27,17,3,8)=4.12268, got: $out"; fail=1; }
grep -qx 'c(10,8,1,3)=3.53774' <<<"$out" || { echo "FAIL: expected c(10,8,1,3)=3.53774, got: $out"; fail=1; }
grep -qx 'c(2,150,7,183)=0.34953' <<<"$out" || { echo "FAIL: expected c(2,150,7,183)=0.34953, got: $out"; fail=1; }

[ "$fail" -eq 0 ] && echo "PASS test_oracle_r_canary"
exit "$fail"
