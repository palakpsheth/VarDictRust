#!/usr/bin/env bash
set -euo pipefail

echo "=== Verifying vardictrust installation in conda environment ==="

# 1. Binary presence
command -v vardictrust

# 2. Version check (asserts semantic version regex)
version_output="$(vardictrust --version)"
echo "Output: $version_output"
echo "$version_output" | grep -E "^vardictrust [0-9]+\.[0-9]+\.[0-9]+"

# 3. Help invocation
vardictrust --help | grep -i "Usage:"

# 4. Error-handling contract (invalid flag should exit non-zero without panic)
set +e
vardictrust --nonexistent-flag 2>/dev/null
status=$?
set -e
if [ "$status" -eq 0 ]; then
  echo "FAIL: Expected non-zero exit on invalid flag, got 0" >&2
  exit 1
fi

echo "=== vardictrust conda smoke test passed successfully ==="
