#!/usr/bin/env bash
set -euo pipefail
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT; fail=0

mkdir -p "$T/benchmarks/evidence"

cat > "$T/benchmarks/evidence/case1.json" <<'EOF'
{
  "fingerprint": "abc123",
  "java_time": {"wall": "1:00.00", "max_rss_kb": "2000000", "exit": "0"},
  "rust_time": {"wall": "0:30.00", "max_rss_kb": "1000000", "exit": "0"}
}
EOF

cat > "$T/benchmarks/evidence/case_wall_fail.json" <<'EOF'
{
  "fingerprint": "abc123",
  "java_time": {"wall": "1:00.00", "max_rss_kb": "2000000", "exit": "0"},
  "rust_time": {"wall": "1:01.00", "max_rss_kb": "1000000", "exit": "0"}
}
EOF

cat > "$T/benchmarks/evidence/case_rss_fail.json" <<'EOF'
{
  "fingerprint": "abc123",
  "java_time": {"wall": "1:00.00", "max_rss_kb": "2000000", "exit": "0"},
  "rust_time": {"wall": "0:30.00", "max_rss_kb": "2000001", "exit": "0"}
}
EOF

cat > "$T/benchmarks/evidence/case_missing_time.json" <<'EOF'
{
  "fingerprint": "abc123",
  "java_time": {"wall": "1:00.00", "max_rss_kb": "2000000", "exit": "0"},
  "rust_time": null
}
EOF

HERE=$PWD

# Test 1: rust faster and smaller -> exit 0
if (cd "$T" && python3 "$HERE/harness/perf_budget.py" case1 >/dev/null 2>&1); then
  :
else
  echo "FAIL: case1 should pass"; fail=1
fi

# Test 2: rust wall greater by 1s -> exit 1
if (cd "$T" && python3 "$HERE/harness/perf_budget.py" case_wall_fail >/dev/null 2>&1); then
  echo "FAIL: case_wall_fail should fail"; fail=1
fi

# Test 3: rust RSS greater by 1 kB -> exit 1
if (cd "$T" && python3 "$HERE/harness/perf_budget.py" case_rss_fail >/dev/null 2>&1); then
  echo "FAIL: case_rss_fail should fail"; fail=1
fi

# Test 4: missing rust_time -> non-zero
if (cd "$T" && python3 "$HERE/harness/perf_budget.py" case_missing_time >/dev/null 2>&1); then
  echo "FAIL: case_missing_time should fail"; fail=1
fi

if [ $fail -eq 0 ]; then
  echo "PASS test_perf_budget"
  exit 0
else
  exit 1
fi
