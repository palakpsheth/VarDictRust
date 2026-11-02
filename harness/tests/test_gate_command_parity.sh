#!/usr/bin/env bash
# harness/tests/test_gate_command_parity.sh
# Verifies that all gate scripts, git hooks, CI workflows, and AGENTS.md
# invoke the exact same cargo fmt & cargo clippy commands without -A drift,
# and that rust-toolchain.toml pins the exact version from benchmarks/env/toolchain.txt.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$ROOT_DIR"

fail=0
ck() {
  local desc="$1"
  shift
  if ! "$@"; then
    echo "FAIL test_gate_command_parity: $desc" >&2
    fail=1
  fi
}

EXPECTED_CLIPPY="cargo clippy --all-targets --all-features -- -D warnings"
EXPECTED_FMT="cargo fmt --all -- --check"

GATE_FILES=(
  "harness/rgate.sh"
  "harness/pr_check.sh"
  "harness/acceptance.sh"
  ".githooks/pre-commit"
  ".githooks/pre-push"
  ".github/workflows/ci.yml"
)

for f in "${GATE_FILES[@]}"; do
  ck "$f contains exact '$EXPECTED_CLIPPY'" grep -Fq "$EXPECTED_CLIPPY" "$f"
  ck "$f contains exact '$EXPECTED_FMT'" grep -Fq "$EXPECTED_FMT" "$f"
  if grep -Eq 'clippy::(needless_late_init|byte_char_slices)' "$f"; then
    echo "FAIL test_gate_command_parity: $f still suppresses clippy lints via -A" >&2
    fail=1
  fi
done

ck "AGENTS.md documents '$EXPECTED_CLIPPY'" grep -Fq "$EXPECTED_CLIPPY" AGENTS.md

if grep -Eq 'clippy::(needless_late_init|byte_char_slices)' .github/workflows/security.yml; then
  echo "FAIL test_gate_command_parity: .github/workflows/security.yml still suppresses clippy lints via -A" >&2
  fail=1
fi
ck "security.yml runs cargo clippy --all-targets --all-features --message-format=json" \
  grep -Fq "cargo clippy --all-targets --all-features --message-format=json" .github/workflows/security.yml

# Verify rust-toolchain.toml is pinned to the exact rustc version in benchmarks/env/toolchain.txt
tc_ver=$(awk '/^rustc /{print $2; exit}' benchmarks/env/toolchain.txt)
[ -n "$tc_ver" ] || { echo "FAIL test_gate_command_parity: could not parse rustc version from benchmarks/env/toolchain.txt" >&2; exit 1; }
ck "rust-toolchain.toml pins channel = \"$tc_ver\"" grep -Fq "channel = \"$tc_ver\"" rust-toolchain.toml
if grep -Fq 'channel = "stable"' rust-toolchain.toml; then
  echo "FAIL test_gate_command_parity: rust-toolchain.toml is unpinned (channel = \"stable\")" >&2
  fail=1
fi

if [ "$fail" -eq 0 ]; then
  echo "PASS test_gate_command_parity"
  exit 0
else
  exit 1
fi
