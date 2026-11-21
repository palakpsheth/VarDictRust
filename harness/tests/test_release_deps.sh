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

if grep -q 'profile-' <<< "$DEFAULT_FEATURES"; then
  echo "FAIL test_release_deps: default features contains profile-* feature: $DEFAULT_FEATURES" >&2
  exit 1
fi

CORE_DEFAULT_FEATURES=$(awk '
  /^\[features\]/ { in_features=1; next }
  /^\[/ { in_features=0 }
  in_features && /default *=/ { in_default=1 }
  in_default {
    buf = buf " " $0
    if (/\]/) { in_default=0 }
  }
  END { print buf }
' crates/vardictrust-core/Cargo.toml)

if grep -q 'profile-' <<< "$CORE_DEFAULT_FEATURES"; then
  echo "FAIL test_release_deps: vardictrust-core default features contains profile-* feature: $CORE_DEFAULT_FEATURES" >&2
  exit 1
fi

# 2. Check default dependency tree (normal and build dependencies)
DEFAULT_TREE=$(cargo tree -p vardictrust -e normal,build)
if grep -E -q "$FORBIDDEN_REGEX" <<< "$DEFAULT_TREE"; then
  echo "FAIL test_release_deps: forbidden profiling dependency found in default dependency tree:" >&2
  grep -E "$FORBIDDEN_REGEX" <<< "$DEFAULT_TREE" >&2
  exit 1
fi

DEFAULT_FEAT_TREE=$(cargo tree -p vardictrust -e features)
if grep -q 'profile-counters' <<< "$DEFAULT_FEAT_TREE"; then
  echo "FAIL test_release_deps: profile-counters enabled in default feature tree" >&2
  exit 1
fi

# 3. Check --no-default-features dependency tree
NO_DEFAULT_TREE=$(cargo tree -p vardictrust -e normal,build --no-default-features)
if grep -E -q "$FORBIDDEN_REGEX" <<< "$NO_DEFAULT_TREE"; then
  echo "FAIL test_release_deps: forbidden profiling dependency found in --no-default-features tree:" >&2
  grep -E "$FORBIDDEN_REGEX" <<< "$NO_DEFAULT_TREE" >&2
  exit 1
fi

NO_DEFAULT_FEAT_TREE=$(cargo tree -p vardictrust -e features --no-default-features)
if grep -q 'profile-counters' <<< "$NO_DEFAULT_FEAT_TREE"; then
  echo "FAIL test_release_deps: profile-counters enabled in --no-default-features feature tree" >&2
  exit 1
fi

echo "PASS test_release_deps"
