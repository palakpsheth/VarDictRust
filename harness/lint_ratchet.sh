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

# Ensure no unapproved #[allow(...)] or #[expect(...)] bypasses (including multiline
# attributes or cfg_attr wrappers) exist in crates/vardictrust-core/src/**, and that
# the single approved `clippy::unwrap_used` allowance in `util.rs` is attached to
# `pub(crate) fn static_regex`.
python3 - <<'PY'
import pathlib
import re
import sys

RATCHET_LINTS = (
    "clippy::unwrap_used",
    "clippy::expect_used",
    "clippy::panic",
    "clippy::unreachable",
    "clippy::todo",
    "clippy::string_slice",
    "unwrap_used",
    "expect_used",
    "string_slice",
)

ATTR_RE = re.compile(
    r"#!?\[(?P<body>(?:[^\[\]]|\[(?:[^\[\]]|\[[^\[\]]*\])*\])*)\]",
    re.DOTALL,
)
SUPPRESS_KIND_RE = re.compile(r"\b(?:allow|expect)\s*\(")
LINT_TOKEN_RE = re.compile(
    r"(?:clippy::(?:unwrap_used|expect_used|panic|unreachable|todo|string_slice)|\b(?:unwrap_used|expect_used|string_slice)\b)"
)


def scan_source(rel_path: str, text: str) -> list[tuple[str, int, str, str]]:
    hits: list[tuple[str, int, str, str]] = []
    for m in ATTR_RE.finditer(text):
        body = m.group("body")
        if not SUPPRESS_KIND_RE.search(body):
            continue
        lint_m = LINT_TOKEN_RE.search(body)
        if not lint_m:
            continue
        line_no = text.count("\n", 0, m.start()) + 1
        tail = text[m.end() : m.end() + 240]
        hits.append((rel_path, line_no, m.group(0), tail))
    return hits


# Self-test: verify multiline #[allow(...)] and #[expect(...)] detection
assert len(
    scan_source(
        "synthetic.rs",
        "#[allow(\n    clippy::unwrap_used,\n    reason = \"bad\"\n)]\nfn bad() {}\n",
    )
) == 1, "multiline #[allow(...)] must be detected"
assert len(
    scan_source(
        "synthetic.rs",
        "#[expect(\n    clippy::string_slice\n)]\nfn bad2() {}\n",
    )
) == 1, "multiline #[expect(...)] must be detected"
assert len(
    scan_source(
        "synthetic.rs",
        "#[cfg_attr(not(test), allow(clippy::panic))]\nfn bad3() {}\n",
    )
) == 1, "cfg_attr(..., allow(...)) must be detected"

src_root = pathlib.Path("crates/vardictrust-core/src")
all_hits: list[tuple[str, int, str, str]] = []
for rs_path in sorted(src_root.rglob("*.rs")):
    rel = rs_path.as_posix()
    all_hits.extend(scan_source(rel, rs_path.read_text(encoding="utf-8")))

approved: list[tuple[str, int, str, str]] = []
unapproved: list[tuple[str, int, str, str]] = []

for hit in all_hits:
    rel, line_no, attr_src, tail = hit
    has_only_unwrap = (
        "clippy::unwrap_used" in attr_src
        and not any(
            tok in attr_src
            for tok in (
                "expect_used",
                "clippy::panic",
                "clippy::unreachable",
                "clippy::todo",
                "string_slice",
            )
        )
    )
    if (
        rel == "crates/vardictrust-core/src/util.rs"
        and has_only_unwrap
        and "fn static_regex(" in tail
    ):
        approved.append(hit)
    else:
        unapproved.append(hit)

if unapproved:
    print(
        "FAIL: unapproved #[allow(...)] / #[expect(...)] on no-panic/string_slice lint in crates/vardictrust-core/src/:",
        file=sys.stderr,
    )
    for rel, line_no, attr_src, _ in unapproved:
        compact = " ".join(attr_src.split())
        print(f"  {rel}:{line_no}: {compact}", file=sys.stderr)
    sys.exit(1)

if len(approved) != 1:
    print(
        f"FAIL: expected exactly 1 approved static_regex #[allow(clippy::unwrap_used)] in crates/vardictrust-core/src/util.rs, found {len(approved)}",
        file=sys.stderr,
    )
    sys.exit(1)
PY

cargo clippy --all-targets --all-features -- -D warnings

echo "PASS: lint_ratchet.sh (no-panic + string_slice deny ratchet clean)"
