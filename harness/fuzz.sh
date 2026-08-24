#!/usr/bin/env bash
# Fuzzing wrapper (see CONTRIBUTING.md "Fuzzing").
#
# Usage: bash harness/fuzz.sh <target> [seconds]
#   <target>  fuzz_cigar | fuzz_regions
#   [seconds] max runtime per invocation, default 300
#
# Requires the nightly toolchain pinned in fuzz/rust-toolchain.toml and
# cargo-fuzz (cargo install cargo-fuzz --locked). Any crash/panic found is a
# real no-panic-rule bug (CLAUDE.md §9) — triage per the trust-batch policy,
# never silence.
set -euo pipefail

TARGET="${1:?usage: harness/fuzz.sh <fuzz_cigar|fuzz_regions> [seconds]}"
SECS="${2:-300}"

cd "$(dirname "$0")/.."
exec cargo +nightly fuzz run "$TARGET" -- -max_total_time="$SECS"
