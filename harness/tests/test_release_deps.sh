#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"

FORBIDDEN_REGEX='\<pprof\>|\<dhat\>|\<inferno\>|\<jemalloc_pprof\>|\<pprof_util\>'

# 1. Assert [features] default in crates/vardictrust/Cargo.toml contains no profile-* feature.
# Robust multi-line array parsing handles single-line and multi-line bracketed arrays.
DEFAULT_FEATURES=$(awk '
  /^\[features\]/ { in_features=1; next }
  /^\[/ { in_features=0 }
  in_features && /default *=/ { in_default=1 }
  in_default {
    buf = buf " " $0
    if (/\]/) { in_default=0 }
  }
  END { print buf }
' crates/vardictrust/Cargo.toml)

if echo "$DEFAULT_FEATURES" | grep -q 'profile-'; then
  echo "FAIL test_release_deps: default features contains profile-* feature: $DEFAULT_FEATURES" >&2
  exit 1
fi

# 2. Check default dependency tree (normal and build dependencies)
DEFAULT_TREE=$(cargo tree -p vardictrust -e normal,build)
if echo "$DEFAULT_TREE" | grep -E -q "$FORBIDDEN_REGEX"; then
  echo "FAIL test_release_deps: forbidden profiling dependency found in default dependency tree:" >&2
  echo "$DEFAULT_TREE" | grep -E "$FORBIDDEN_REGEX" >&2
  exit 1
fi

# 3. Check --no-default-features dependency tree
NO_DEFAULT_TREE=$(cargo tree -p vardictrust -e normal,build --no-default-features)
if echo "$NO_DEFAULT_TREE" | grep -E -q "$FORBIDDEN_REGEX"; then
  echo "FAIL test_release_deps: forbidden profiling dependency found in --no-default-features tree:" >&2
  echo "$NO_DEFAULT_TREE" | grep -E "$FORBIDDEN_REGEX" >&2
  exit 1
fi

echo "PASS test_release_deps"
