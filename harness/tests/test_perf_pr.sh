#!/usr/bin/env bash
# harness/tests/test_perf_pr.sh
# Tests harness/perf_pr.sh using lightweight stub binaries to verify:
# 1. Slowdown exceeding threshold on both median and min exits non-zero.
# 2. Slowdown within threshold passes (exits 0).
# 3. Non-zero exit code from stub binary fails (exits non-zero).
# 4. Output digest mismatch fails unless --allow-vcf-diff is passed.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"

TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

BASE_STUB="$TMP_DIR/base_vardictrust"
HEAD_SLOW="$TMP_DIR/head_slow"
HEAD_CRASH="$TMP_DIR/head_crash"
HEAD_DIFF="$TMP_DIR/head_diff"

cat > "$BASE_STUB" <<'EOF'
#!/usr/bin/env bash
python3 -c 'import time; time.sleep(0.08)'
echo "##fileformat=VCFv4.2"
echo -e "chr1\t101\t.\tT\tA\t60\tPASS\tDP=15\tGT\t0/1"
EOF

cat > "$HEAD_SLOW" <<'EOF'
#!/usr/bin/env bash
python3 -c 'import time; time.sleep(0.16)'
echo "##fileformat=VCFv4.2"
echo -e "chr1\t101\t.\tT\tA\t60\tPASS\tDP=15\tGT\t0/1"
EOF

cat > "$HEAD_CRASH" <<'EOF'
#!/usr/bin/env bash
exit 101
EOF

cat > "$HEAD_DIFF" <<'EOF'
#!/usr/bin/env bash
python3 -c 'import time; time.sleep(0.08)'
echo "##fileformat=VCFv4.2"
echo -e "chr1\t101\t.\tT\tG\t60\tPASS\tDP=15\tGT\t0/1"
EOF

chmod +x "$BASE_STUB" "$HEAD_SLOW" "$HEAD_CRASH" "$HEAD_DIFF"

fail=0

# 1. 2x slowdown (0.16s vs 0.08s) with --threshold 1.10 must FAIL
if bash harness/perf_pr.sh \
  --base-bin "$BASE_STUB" \
  --head-bin "$HEAD_SLOW" \
  --slices micro \
  --runs 3 \
  --threshold 1.10 \
  --out-json "$TMP_DIR/out1.json" >/dev/null 2>&1; then
  echo "FAIL test_perf_pr: 2x slowdown should fail at threshold 1.10" >&2
  fail=1
fi

# 2. Same 2x slowdown with --threshold 3.00 must PASS
if ! bash harness/perf_pr.sh \
  --base-bin "$BASE_STUB" \
  --head-bin "$HEAD_SLOW" \
  --slices micro \
  --runs 3 \
  --threshold 3.00 \
  --out-json "$TMP_DIR/out2.json" >/dev/null 2>&1; then
  echo "FAIL test_perf_pr: 2x slowdown should pass at threshold 3.00" >&2
  fail=1
fi

# 3. Crashing head binary (exit 101) must FAIL
if bash harness/perf_pr.sh \
  --base-bin "$BASE_STUB" \
  --head-bin "$HEAD_CRASH" \
  --slices micro \
  --runs 2 \
  --threshold 2.00 \
  --out-json "$TMP_DIR/out3.json" >/dev/null 2>&1; then
  echo "FAIL test_perf_pr: crashing head stub should fail" >&2
  fail=1
fi

# 4. Output VCF diff must FAIL by default, and PASS with --allow-vcf-diff
if bash harness/perf_pr.sh \
  --base-bin "$BASE_STUB" \
  --head-bin "$HEAD_DIFF" \
  --slices micro \
  --runs 2 \
  --threshold 2.00 \
  --out-json "$TMP_DIR/out4.json" >/dev/null 2>&1; then
  echo "FAIL test_perf_pr: diverging VCF output should fail without --allow-vcf-diff" >&2
  fail=1
fi

if ! bash harness/perf_pr.sh \
  --base-bin "$BASE_STUB" \
  --head-bin "$HEAD_DIFF" \
  --slices micro \
  --runs 2 \
  --threshold 2.00 \
  --allow-vcf-diff \
  --out-json "$TMP_DIR/out5.json" >/dev/null 2>&1; then
  echo "FAIL test_perf_pr: diverging VCF output should pass with --allow-vcf-diff" >&2
  fail=1
fi

if [ $fail -eq 0 ]; then
  echo "PASS test_perf_pr"
  exit 0
else
  exit 1
fi
