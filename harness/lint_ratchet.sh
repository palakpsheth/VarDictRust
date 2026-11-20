#!/usr/bin/env bash
# harness/lint_ratchet.sh — Enforces AGENTS.md §9 no-panic / byte-safe string-slicing
# and bounds-checked slice-indexing zero-allow lint ratchet on `crates/vardictrust-core`.
set -euo pipefail

SCAN_ONLY=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --scan-only)
            SCAN_ONLY=1
            shift
            ;;
        *)
            echo "ERROR: unknown argument '$1' (usage: $0 [--scan-only])" >&2
            exit 2
            ;;
    esac
done

REPO_ROOT="${LINT_RATCHET_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$REPO_ROOT"

LIB_RS="crates/vardictrust-core/src/lib.rs"
if [[ ! -f "$LIB_RS" ]]; then
    echo "FAIL: missing $LIB_RS" >&2
    exit 1
fi

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
    "clippy::string_slice" \
    "clippy::indexing_slicing"
do
    if ! grep -q "$lint" "$LIB_RS"; then
        echo "FAIL: $LIB_RS does not deny $lint" >&2
        exit 1
    fi
done

# Ensure zero #[allow(...)] or #[expect(...)] bypasses (including multiline
# attributes or cfg_attr wrappers) exist in non-test source across
# crates/vardictrust-core/src/**, and verify lib.rs header attributes on
# comment/string/cfg(test)-stripped source.
python3 - <<'PY'
import pathlib
import re
import sys

REQUIRED_DENY_LINTS = (
    "clippy::unwrap_used",
    "clippy::expect_used",
    "clippy::panic",
    "clippy::unreachable",
    "clippy::todo",
    "clippy::string_slice",
    "clippy::indexing_slicing",
)

_LEX_RE = re.compile(
    r"""
    (?P<line_comment>//[^\n]*)
    | (?P<block_open>/\*)
    | (?<![A-Za-z0-9_'])(?P<raw_prefix>br|rb|cr|rc|r)(?P<raw_hashes>\#*)"
    | (?<![A-Za-z0-9_'])(?:b|c)?(?P<str_open>")
    | (?P<char_lit>(?<![A-Za-z0-9_'])b?'(?:\\(?:x[0-9a-fA-F]{2}|u\{[0-9a-fA-F_]+\}|[\s\S])|[^'\\\n])')
    """,
    re.VERBOSE,
)
_BLOCK_STEP_RE = re.compile(r"/\*|\*/")
_STR_STEP_RE = re.compile(r'\\[\s\S]|"')
_CFG_TEST_INNER_RE = re.compile(r"#\s*!\s*\[\s*cfg\s*\(\s*test\s*\)\s*\]")
_CFG_TEST_OUTER_RE = re.compile(r"#\s*(?!!)\[\s*cfg\s*\(\s*test\s*\)\s*\]")
_IDENT_RE = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
_DELIM_RE = re.compile(r"[{}\[\]();]")

ATTR_RE = re.compile(
    r"#\s*!?\s*\[(?P<body>(?:[^\[\]]|\[(?:[^\[\]]|\[[^\[\]]*\])*\])*)\]",
    re.DOTALL,
)
SUPPRESS_KIND_RE = re.compile(r"\b(?:allow|expect)\s*\(")
DENY_KIND_RE = re.compile(r"\bdeny\s*\(")
FORBID_UNSAFE_RE = re.compile(r"#\s*!\s*\[\s*forbid\s*\(\s*unsafe_code\s*\)\s*\]")
LINT_TOKEN_RE = re.compile(
    r"(?:clippy\s*::\s*(?:unwrap_used|expect_used|panic|unreachable|todo|string_slice|indexing_slicing|restriction|all)\b"
    r"|\b(?:unwrap_used|expect_used|panic|unreachable|todo|string_slice|indexing_slicing|warnings)\b)"
)


def _blank_segment(segment: str) -> str:
    return "\n".join(" " * len(part) for part in segment.split("\n"))


def _skip_balanced(text: str, open_idx: int, open_ch: str, close_ch: str) -> int:
    depth = 0
    n = len(text)
    i = open_idx
    while i < n:
        ch = text[i]
        if ch == open_ch:
            depth += 1
        elif ch == close_ch:
            depth -= 1
            if depth == 0:
                return i + 1
        i += 1
    return n


def strip_comments_and_strings(text: str) -> str:
    chunks: list[str] = []
    pos = 0
    n = len(text)
    while pos < n:
        m = _LEX_RE.search(text, pos)
        if not m:
            chunks.append(text[pos:])
            break
        start = m.start()
        if start > pos:
            chunks.append(text[pos:start])
        if m.group("line_comment") is not None or m.group("char_lit") is not None:
            end = m.end()
        elif m.group("block_open") is not None:
            depth = 1
            end = n
            for bm in _BLOCK_STEP_RE.finditer(text, m.end()):
                if bm.group(0) == "/*":
                    depth += 1
                else:
                    depth -= 1
                    if depth == 0:
                        end = bm.end()
                        break
        elif m.group("raw_prefix") is not None:
            closing = '"' + (m.group("raw_hashes") or "")
            close_idx = text.find(closing, m.end())
            end = n if close_idx == -1 else close_idx + len(closing)
        else:
            end = n
            for sm in _STR_STEP_RE.finditer(text, m.end()):
                if sm.group(0) == '"':
                    end = sm.end()
                    break
        chunks.append(_blank_segment(text[start:end]))
        pos = end
    return "".join(chunks)


def strip_cfg_test_blocks(code: str) -> str:
    if _CFG_TEST_INNER_RE.search(code):
        return _blank_segment(code)
    chunks: list[str] = []
    cursor = 0
    n = len(code)
    while cursor < n:
        m = _CFG_TEST_OUTER_RE.search(code, cursor)
        if not m:
            chunks.append(code[cursor:])
            break
        cfg_start = m.start()
        pos = m.end()
        while pos < n:
            while pos < n and code[pos].isspace():
                pos += 1
            if pos < n and code[pos] == "#":
                k = pos + 1
                while k < n and code[k].isspace():
                    k += 1
                if k < n and code[k] == "[":
                    pos = _skip_balanced(code, k, "[", "]")
                    continue
            break

        probe = pos
        words: list[str] = []
        while probe < n and len(words) < 6:
            while probe < n and code[probe].isspace():
                probe += 1
            if probe < n and code[probe] == "(":
                probe = _skip_balanced(code, probe, "(", ")")
                continue
            wm = _IDENT_RE.match(code, probe)
            if not wm:
                break
            words.append(wm.group(0))
            probe = wm.end()

        idx_w = 0
        if idx_w < len(words) and words[idx_w] == "pub":
            idx_w += 1
        if idx_w < len(words) and words[idx_w] == "default":
            idx_w += 1
        first_kw = words[idx_w] if idx_w < len(words) else ""
        second_kw = words[idx_w + 1] if idx_w + 1 < len(words) else ""
        semi_only = first_kw in ("use", "type", "static") or (
            first_kw == "const" and second_kw not in ("fn", "unsafe", "async", "extern")
        )

        brace = 0
        paren = 0
        bracket = 0
        item_end = n
        for dm in _DELIM_RE.finditer(code, pos):
            ch = dm.group(0)
            if ch == "[":
                bracket += 1
            elif ch == "]":
                bracket = max(0, bracket - 1)
            elif ch == "(":
                paren += 1
            elif ch == ")":
                paren = max(0, paren - 1)
            elif ch == "{":
                brace += 1
            elif ch == "}":
                if brace > 0:
                    brace -= 1
                    if brace == 0 and paren == 0 and bracket == 0 and not semi_only:
                        item_end = dm.end()
                        break
                else:
                    item_end = dm.start()
                    break
            elif ch == ";":
                if brace == 0 and paren == 0 and bracket == 0:
                    item_end = dm.end()
                    break

        if cfg_start > cursor:
            chunks.append(code[cursor:cfg_start])
        chunks.append(_blank_segment(code[cfg_start:item_end]))
        cursor = item_end
    return "".join(chunks)


def sanitize_non_test_rust(raw_text: str) -> str:
    return strip_cfg_test_blocks(strip_comments_and_strings(raw_text))


def scan_source(rel_path: str, raw_text: str) -> list[tuple[str, int, str]]:
    clean = sanitize_non_test_rust(raw_text)
    hits: list[tuple[str, int, str]] = []
    for m in ATTR_RE.finditer(clean):
        body = m.group("body")
        matched = False
        for sm in SUPPRESS_KIND_RE.finditer(body):
            open_paren = sm.end() - 1
            close_paren = _skip_balanced(body, open_paren, "(", ")")
            suppress_args = body[open_paren + 1 : max(open_paren + 1, close_paren - 1)]
            if LINT_TOKEN_RE.search(suppress_args):
                matched = True
                break
        if not matched:
            continue
        line_no = raw_text.count("\n", 0, m.start()) + 1
        hits.append((rel_path, line_no, raw_text[m.start() : m.end()]))
    return hits


def verify_lib_rs_header(raw_text: str) -> list[str]:
    clean = sanitize_non_test_rust(raw_text)
    errs: list[str] = []
    if not FORBID_UNSAFE_RE.search(clean):
        errs.append("missing active #![forbid(unsafe_code)]")
    denied_text_parts: list[str] = []
    for m in ATTR_RE.finditer(clean):
        attr_src = m.group(0)
        if not re.match(r"#\s*!", attr_src):
            continue
        body = m.group("body")
        for dm in DENY_KIND_RE.finditer(body):
            open_paren = dm.end() - 1
            close_paren = _skip_balanced(body, open_paren, "(", ")")
            denied_text_parts.append(body[open_paren + 1 : max(open_paren + 1, close_paren - 1)])
    denied_blob = "\n".join(denied_text_parts)
    for lint in REQUIRED_DENY_LINTS:
        if not re.search(rf"\b{re.escape(lint)}\b", denied_blob):
            errs.append(f"missing active crate-level deny({lint})")
    return errs


# Inline synthetic self-tests:
# 1) All ratcheted non-test suppressions MUST be flagged
POSITIVE_CASES = (
    ("multiline allow unwrap_used", '#[allow(\n    clippy::unwrap_used,\n    reason = "bad"\n)]\nfn f() {}\n'),
    ("multiline expect string_slice", "#[expect(\n    clippy::string_slice\n)]\nfn f() {}\n"),
    ("cfg_attr allow clippy::panic", "#[cfg_attr(not(test), allow(clippy::panic))]\nfn f() {}\n"),
    ("allow clippy::indexing_slicing", "#[allow(clippy::indexing_slicing)]\nfn f() {}\n"),
    ("expect clippy::indexing_slicing", "#[expect(clippy::indexing_slicing)]\nfn f() {}\n"),
    (
        "cfg_attr allow clippy::indexing_slicing",
        "#[cfg_attr(not(test), allow(clippy::indexing_slicing))]\nfn f() {}\n",
    ),
    ("allow clippy::all", "#[allow(clippy::all)]\nfn f() {}\n"),
    ("allow clippy::restriction", "#[allow(clippy::restriction)]\nfn f() {}\n"),
    ("allow warnings", "#[allow(warnings)]\nfn f() {}\n"),
    ("allow bare panic", "#[allow(panic)]\nfn f() {}\n"),
    ("allow bare unreachable", "#[allow(unreachable)]\nfn f() {}\n"),
    ("allow bare todo", "#[allow(todo)]\nfn f() {}\n"),
)
for label, snippet in POSITIVE_CASES:
    got = scan_source("synthetic.rs", snippet)
    assert len(got) == 1, f"self-test failed ({label}): expected 1 hit, got {got}"

# 2) Occurrences inside comments, string/raw-string literals, reason strings, and #[cfg(test)] blocks MUST be ignored
NEGATIVE_SNIPPET = """
// #[allow(clippy::indexing_slicing)]
/* #[allow(clippy::all)] /* nested #[allow(warnings)] */ #[allow(panic)] */
const S1: &str = "#[allow(clippy::indexing_slicing)] \\" #[allow(todo)]";
const S2: &str = r#"#[expect(clippy::indexing_slicing)]"#;
const S3: &[u8] = br##"#[allow(clippy::restriction)] #[allow(unreachable)]"##;
const Q: char = '"';
#[expect(deprecated, reason = "mentions panic and indexing_slicing inside reason literal")]
#[allow(clippy::too_many_arguments)]
fn prod_clean() {}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
fn test_helper() {}

#[cfg(test)]
use std::collections::{BTreeMap, HashMap};

#[cfg(test)]
mod tests {
    #[allow(clippy::indexing_slicing)]
    #[allow(clippy::all)]
    #[allow(warnings)]
    fn inside_test_mod() {}
}
"""
assert scan_source("synthetic_clean.rs", NEGATIVE_SNIPPET) == [], (
    "self-test failed: comments, strings, raw strings, and #[cfg(test)] blocks must be ignored"
)

# 3) Code AFTER a #[cfg(test)] item must still be scanned!
POST_TEST_SNIPPET = NEGATIVE_SNIPPET + "\n#[allow(clippy::indexing_slicing)]\nfn after_test_mod() {}\n"
assert len(scan_source("synthetic_post.rs", POST_TEST_SNIPPET)) == 1, (
    "self-test failed: non-test code following #[cfg(test)] block must still be scanned"
)

# 4) Self-test header verifier against commented-out forbid/deny attributes
assert verify_lib_rs_header(
    "// #![forbid(unsafe_code)]\n#![deny(clippy::unwrap_used)]\n"
), "self-test failed: commented-out forbid(unsafe_code) and missing deny lints must fail"
assert verify_lib_rs_header(
    "#![forbid(unsafe_code)]\n"
    "// clippy::indexing_slicing\n"
    "#![cfg_attr(not(test), deny(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::unreachable, clippy::todo, clippy::string_slice))]\n"
), "self-test failed: missing active deny(clippy::indexing_slicing) must fail"

# 5) Verify lib.rs active header attributes
lib_path = pathlib.Path("crates/vardictrust-core/src/lib.rs")
header_errs = verify_lib_rs_header(lib_path.read_text(encoding="utf-8"))
if header_errs:
    for err in header_errs:
        print(f"FAIL: {lib_path.as_posix()}: {err}", file=sys.stderr)
    sys.exit(1)

src_root = pathlib.Path("crates/vardictrust-core/src")
all_hits: list[tuple[str, int, str]] = []
for rs_path in sorted(src_root.rglob("*.rs")):
    rel = rs_path.as_posix()
    all_hits.extend(scan_source(rel, rs_path.read_text(encoding="utf-8")))

if all_hits:
    print(
        "FAIL: #[allow(...)] / #[expect(...)] bypass on ratcheted lint in crates/vardictrust-core/src/ (zero-allow policy):",
        file=sys.stderr,
    )
    for rel, line_no, attr_src in all_hits:
        compact = " ".join(attr_src.split())
        print(f"  {rel}:{line_no}: {compact}", file=sys.stderr)
    sys.exit(1)
PY

if [[ "$SCAN_ONLY" -eq 1 ]]; then
    echo "PASS: lint_ratchet.sh --scan-only (zero-allow + indexing_slicing header scan clean)"
    exit 0
fi

cargo clippy --all-targets --all-features -- -D warnings

echo "PASS: lint_ratchet.sh (no-panic + string_slice + indexing_slicing zero-allow ratchet clean)"
