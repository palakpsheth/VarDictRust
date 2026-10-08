#!/usr/bin/env bash
set -euo pipefail

ORACLE_YML="benchmarks/env/vardict-oracle.yml"
G4_YML="benchmarks/env/g4-tools.yml"
TOOLCHAIN_TXT="benchmarks/env/toolchain.txt"

# 1. vardict-oracle.yml checks
[ -f "$ORACLE_YML" ] || { echo "FAIL test_env_pinned: $ORACLE_YML missing" >&2; exit 1; }
grep -qE 'vardict-java(=1\.8\.3)?' "$ORACLE_YML" || {
  echo "FAIL test_env_pinned: vardict-java not found in $ORACLE_YML" >&2
  exit 1
}

# 2. g4-tools.yml checks
[ -f "$G4_YML" ] || { echo "FAIL test_env_pinned: $G4_YML missing" >&2; exit 1; }
grep -qE 'rtg-tools' "$G4_YML" || {
  echo "FAIL test_env_pinned: rtg-tools not found in $G4_YML" >&2
  exit 1
}
grep -qE 'samtools' "$G4_YML" || {
  echo "FAIL test_env_pinned: samtools not found in $G4_YML" >&2
  exit 1
}

# 3. toolchain.txt checks
[ -f "$TOOLCHAIN_TXT" ] || { echo "FAIL test_env_pinned: $TOOLCHAIN_TXT missing" >&2; exit 1; }
grep -qi 'rustc' "$TOOLCHAIN_TXT" || {
  echo "FAIL test_env_pinned: rustc not found in $TOOLCHAIN_TXT" >&2
  exit 1
}
grep -qiE 'java|openjdk' "$TOOLCHAIN_TXT" || {
  echo "FAIL test_env_pinned: java/openjdk not found in $TOOLCHAIN_TXT" >&2
  exit 1
}

echo "PASS test_env_pinned"
