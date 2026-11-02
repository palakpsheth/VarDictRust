#!/usr/bin/env bash
# harness/lint_ratchet.sh — Enforces AGENTS.md §9 no-panic / byte-safe string-slicing
# lint ratchet on `crates/vardictrust-core`.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

LIB_RS="crates/vardictrust-core/src/lib.rs"

if ! grep -q '#!\[forbid(unsafe_code)\]' "$LIB_RS"; then
    echo "FAIL: $LIB_RS is missing #![forbid(unsafe_code)]" >&2
    exit 1
fi

for lint in \
    "clippy::unwrap_used" \
    "clippy::expect_used" \
    "clippy::panic" \
    "clippy::unreachable" \
    "clippy::todo" \
    "clippy::string_slice"
do
    if ! grep -q "$lint" "$LIB_RS"; then
        echo "FAIL: $LIB_RS does not deny $lint" >&2
        exit 1
    fi
done

# Ensure no unapproved #[allow(...)] bypasses exist in crates/vardictrust-core/src/**
# (only `crates/vardictrust-core/src/util.rs` may allow `clippy::unwrap_used` for `static_regex`).
UNAPPROVED="$(
    grep -rnE '#!?\[allow\([^]]*(unwrap_used|expect_used|clippy::panic|clippy::unreachable|clippy::todo|string_slice)' \
        crates/vardictrust-core/src/ \
        | grep -v '^crates/vardictrust-core/src/util.rs:' \
        || true
)"
if [[ -n "$UNAPPROVED" ]]; then
    echo "FAIL: unapproved #[allow(...)] on no-panic/string_slice lint in crates/vardictrust-core/src/:" >&2
    echo "$UNAPPROVED" >&2
    exit 1
fi

UTIL_ALLOW_COUNT="$(grep -c 'clippy::unwrap_used' crates/vardictrust-core/src/util.rs || true)"
if [[ "$UTIL_ALLOW_COUNT" -ne 1 ]]; then
    echo "FAIL: expected exactly 1 clippy::unwrap_used allow in crates/vardictrust-core/src/util.rs, found $UTIL_ALLOW_COUNT" >&2
    exit 1
fi

cargo clippy --all-targets --all-features -- -D warnings

echo "PASS: lint_ratchet.sh (no-panic + string_slice deny ratchet clean)"
