#!/usr/bin/env bash
# harness/lint_ratchet.sh — Enforces AGENTS.md §9 11-lint no-panic / hygiene,
# byte-safe string-slicing, bounds-checked slice-indexing, and unsafe_code
# zero-allow lint ratchet across `crates/vardictrust-core` and `crates/vardictrust`.
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

for req_file in \
    "crates/vardictrust-core/src/lib.rs" \
    "crates/vardictrust/src/main.rs" \
    "crates/vardictrust/src/cli.rs" \
    "crates/vardictrust/src/profile.rs"
do
    if [[ ! -f "$req_file" ]]; then
        echo "FAIL: missing $req_file" >&2
        exit 1
    fi
done

# Ensure zero unauthorized #[allow(...)] or #[expect(...)] bypasses (including
# multiline attributes or cfg_attr wrappers) exist in non-test source across
# crates/vardictrust-core/src/** and crates/vardictrust/src/**, and verify
# crate/module header attributes on comment/string/cfg(test)-stripped source.
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
    "clippy::unimplemented",
    "clippy::dbg_macro",
    "clippy::print_stdout",
    "clippy::print_stderr",
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
    r"(?:clippy\s*::\s*(?:"
    r"unwrap_used|expect_used|panic|unreachable|todo|unimplemented|dbg_macro|"
    r"print_stdout|print_stderr|string_slice|indexing_slicing|"
    r"undocumented_unsafe_blocks|multiple_unsafe_ops_per_block|"
    r"restriction|all"
    r")\b"
    r"|\b(?:"
    r"unwrap_used|expect_used|panic|unreachable|todo|unimplemented|dbg_macro|"
    r"print_stdout|print_stderr|string_slice|indexing_slicing|"
    r"unsafe_code|unsafe_op_in_unsafe_fn|"
    r"undocumented_unsafe_blocks|multiple_unsafe_ops_per_block|"
    r"warnings"
    r")\b)"
)

_MALLOC_CONF_FOLLOW_RE = re.compile(
    r'^\s*#\s*\[\s*export_name\s*=\s*"_rjem_malloc_conf"\s*\]\s*'
    r"pub\s+static\s+MALLOC_CONF\s*:\s*&\s*\[\s*u8\s*;"
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


def _lex_transform(text: str, *, blank_strings: bool) -> str:
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
        is_comment = False
        if m.group("line_comment") is not None:
            end = m.end()
            is_comment = True
        elif m.group("char_lit") is not None:
            end = m.end()
        elif m.group("block_open") is not None:
            is_comment = True
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
        if is_comment or blank_strings:
            chunks.append(_blank_segment(text[start:end]))
        else:
            chunks.append(text[start:end])
        pos = end
    return "".join(chunks)


def strip_comments_and_strings(text: str) -> str:
    return _lex_transform(text, blank_strings=True)


def strip_comments_preserving_strings(text: str) -> str:
    return _lex_transform(text, blank_strings=False)


def _normalize_attr(attr_text: str) -> str:
    no_comments = strip_comments_preserving_strings(attr_text)
    compact = re.sub(r"\s+", "", no_comments)
    while True:
        nxt = re.sub(r",([)\]])", r"\1", compact)
        if nxt == compact:
            break
        compact = nxt
    return compact


def _strip_inner_cfg_test_scopes(code: str) -> str:
    inner_spans = [m.start() for m in _CFG_TEST_INNER_RE.finditer(code)]
    if not inner_spans:
        return code
    chars = list(code)
    brace_stack: list[int] = []
    blank_ranges: list[tuple[int, int]] = []
    span_idx = 0
    i = 0
    n = len(code)
    while i < n and span_idx < len(inner_spans):
        if i == inner_spans[span_idx]:
            if not brace_stack:
                return _blank_segment(code)
            open_brace = brace_stack[-1]
            block_end = _skip_balanced(code, open_brace, "{", "}")
            closed = block_end > open_brace and code[block_end - 1] == "}"
            inner_end = block_end - 1 if closed else block_end
            blank_ranges.append((open_brace + 1, inner_end))
            while span_idx < len(inner_spans) and inner_spans[span_idx] < block_end:
                span_idx += 1
            if closed:
                brace_stack.pop()
            else:
                brace_stack.clear()
            i = block_end
            continue
        ch = code[i]
        if ch == "{":
            brace_stack.append(i)
        elif ch == "}" and brace_stack:
            brace_stack.pop()
        i += 1
    for b_start, b_end in blank_ranges:
        for k in range(b_start, b_end):
            if chars[k] != "\n":
                chars[k] = " "
    return "".join(chars)


def strip_cfg_test_blocks(code: str) -> str:
    code = _strip_inner_cfg_test_scopes(code)
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


def _is_exempt_malloc_conf_expect(
    rel_path: str,
    raw_text: str,
    clean_attr: str,
    start: int,
    end: int,
) -> bool:
    """Permit expect(unsafe_code) solely on the _rjem_malloc_conf #[export_name] statics in main.rs."""
    if rel_path != "crates/vardictrust/src/main.rs":
        return False
    if re.match(r"#\s*!", clean_attr):
        return False
    clean_norm = _normalize_attr(clean_attr)
    if clean_norm not in (
        "#[cfg_attr(not(feature=),expect(unsafe_code,reason=))]",
        "#[cfg_attr(not(feature=),expect(unsafe_code))]",
    ):
        return False
    raw_attr_no_comments = strip_comments_preserving_strings(raw_text[start:end])
    if not re.search(
        r'not\s*\(\s*feature\s*=\s*"profile-cpu"\s*\)', raw_attr_no_comments
    ):
        return False
    tail_no_comments = strip_comments_preserving_strings(raw_text[end : end + 512])
    return bool(_MALLOC_CONF_FOLLOW_RE.match(tail_no_comments))


def scan_source(rel_path: str, raw_text: str) -> list[tuple[str, int, str]]:
    clean = sanitize_non_test_rust(raw_text)
    hits: list[tuple[str, int, str]] = []
    exempt_malloc_conf_hits = 0
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
        if (
            exempt_malloc_conf_hits < 4
            and _is_exempt_malloc_conf_expect(
                rel_path, raw_text, m.group(0), m.start(), m.end()
            )
        ):
            exempt_malloc_conf_hits += 1
            continue
        line_no = raw_text.count("\n", 0, m.start()) + 1
        hits.append((rel_path, line_no, raw_text[m.start() : m.end()]))
    return hits


def _collect_active_inner_attrs_normalized(raw_text: str, clean: str) -> set[str]:
    out: set[str] = set()
    for m in ATTR_RE.finditer(clean):
        if re.match(r"#\s*!", m.group(0)):
            out.add(_normalize_attr(raw_text[m.start() : m.end()]))
    return out


def _collect_active_non_test_denied_blob(clean: str) -> str:
    denied_text_parts: list[str] = []
    for m in ATTR_RE.finditer(clean):
        attr_src = m.group(0)
        if not re.match(r"#\s*!", attr_src):
            continue
        body = m.group("body")
        if not (
            re.match(r"^\s*deny\s*\(", body)
            or re.match(r"^\s*cfg_attr\s*\(\s*not\s*\(\s*test\s*\)\s*,", body)
        ):
            continue
        for dm in DENY_KIND_RE.finditer(body):
            open_paren = dm.end() - 1
            close_paren = _skip_balanced(body, open_paren, "(", ")")
            denied_text_parts.append(
                body[open_paren + 1 : max(open_paren + 1, close_paren - 1)]
            )
    return "\n".join(denied_text_parts)


def verify_lib_rs_header(raw_text: str) -> list[str]:
    clean = sanitize_non_test_rust(raw_text)
    errs: list[str] = []
    if not FORBID_UNSAFE_RE.search(clean):
        errs.append("missing active #![forbid(unsafe_code)]")
    denied_blob = _collect_active_non_test_denied_blob(clean)
    for lint in REQUIRED_DENY_LINTS:
        if not re.search(rf"\b{re.escape(lint)}\b", denied_blob):
            errs.append(f"missing active crate-level deny({lint})")
    if re.search(r"\bunsafe\b", clean):
        errs.append("lib.rs must not contain unsafe items/blocks")
    return errs


def verify_cli_rs_header(raw_text: str) -> list[str]:
    clean = sanitize_non_test_rust(raw_text)
    if not FORBID_UNSAFE_RE.search(clean):
        return ["missing active #![forbid(unsafe_code)]"]
    return []


def verify_profile_rs_header(raw_text: str) -> list[str]:
    clean = sanitize_non_test_rust(raw_text)
    active = _collect_active_inner_attrs_normalized(raw_text, clean)
    errs: list[str] = []
    req_forbid = _normalize_attr(
        '#![cfg_attr(not(feature = "profile-cpu"), forbid(unsafe_code))]'
    )
    req_deny = _normalize_attr(
        '#![cfg_attr(feature = "profile-cpu", deny(unsafe_op_in_unsafe_fn, clippy::undocumented_unsafe_blocks, clippy::multiple_unsafe_ops_per_block))]'
    )
    if req_forbid not in active:
        errs.append(
            'missing active #![cfg_attr(not(feature = "profile-cpu"), forbid(unsafe_code))]'
        )
    if req_deny not in active:
        errs.append(
            'missing active #![cfg_attr(feature = "profile-cpu", deny(unsafe_op_in_unsafe_fn, clippy::undocumented_unsafe_blocks, clippy::multiple_unsafe_ops_per_block))]'
        )
    return errs


def verify_main_rs_header(raw_text: str) -> list[str]:
    clean = sanitize_non_test_rust(raw_text)
    active = _collect_active_inner_attrs_normalized(raw_text, clean)
    errs: list[str] = []
    req_forbid = _normalize_attr(
        '#![cfg_attr(not(any(feature = "profile-cpu", feature = "jemalloc", feature = "profile-heap")), forbid(unsafe_code))]'
    )
    req_deny_unsafe = _normalize_attr(
        '#![cfg_attr(all(any(feature = "jemalloc", feature = "profile-heap"), not(feature = "profile-cpu")), deny(unsafe_code))]'
    )
    req_deny_cpu = _normalize_attr(
        '#![cfg_attr(feature = "profile-cpu", deny(unsafe_op_in_unsafe_fn, clippy::undocumented_unsafe_blocks, clippy::multiple_unsafe_ops_per_block))]'
    )
    if req_forbid not in active:
        errs.append(
            'missing active #![cfg_attr(not(any(feature = "profile-cpu", feature = "jemalloc", feature = "profile-heap")), forbid(unsafe_code))]'
        )
    if req_deny_unsafe not in active:
        errs.append(
            'missing active #![cfg_attr(all(any(feature = "jemalloc", feature = "profile-heap"), not(feature = "profile-cpu")), deny(unsafe_code))]'
        )
    if req_deny_cpu not in active:
        errs.append(
            'missing active #![cfg_attr(feature = "profile-cpu", deny(unsafe_op_in_unsafe_fn, clippy::undocumented_unsafe_blocks, clippy::multiple_unsafe_ops_per_block))]'
        )
    denied_blob = _collect_active_non_test_denied_blob(clean)
    for lint in REQUIRED_DENY_LINTS:
        if not re.search(rf"\b{re.escape(lint)}\b", denied_blob):
            errs.append(f"missing active crate-level deny({lint})")
    if re.search(r"\bunsafe\b", clean):
        errs.append("main.rs must not contain direct unsafe items/blocks (scope profile-cpu unsafe to profile.rs)")
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
    ("allow clippy::unimplemented", "#[allow(clippy::unimplemented)]\nfn f() {}\n"),
    ("allow clippy::dbg_macro", "#[allow(clippy::dbg_macro)]\nfn f() {}\n"),
    ("allow clippy::print_stdout", "#[allow(clippy::print_stdout)]\nfn f() {}\n"),
    ("allow clippy::print_stderr", "#[allow(clippy::print_stderr)]\nfn f() {}\n"),
    ("expect unsafe_code", "#[expect(unsafe_code)]\nfn f() {}\n"),
    ("allow unsafe_op_in_unsafe_fn", "#[allow(unsafe_op_in_unsafe_fn)]\nfn f() {}\n"),
    (
        "allow clippy::undocumented_unsafe_blocks",
        "#[allow(clippy::undocumented_unsafe_blocks)]\nfn f() {}\n",
    ),
    (
        "allow clippy::multiple_unsafe_ops_per_block",
        "#[allow(clippy::multiple_unsafe_ops_per_block)]\nfn f() {}\n",
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

# Unauthorized expect(unsafe_code) inside main.rs on non-_rjem_malloc_conf item MUST be flagged,
# whereas _rjem_malloc_conf static ONLY in crates/vardictrust/src/main.rs is exempted.
MALLOC_CONF_SNIPPET = (
    '#[cfg_attr(not(feature = "profile-cpu"), expect(unsafe_code, reason = "jemalloc C ABI symbol _rjem_malloc_conf requires #[export_name] on thin &[u8; N] static"))]\n'
    '#[export_name = "_rjem_malloc_conf"]\n'
    'pub static MALLOC_CONF: &[u8; 4] = b"ok!\\0";\n'
)
assert scan_source("crates/vardictrust/src/main.rs", MALLOC_CONF_SNIPPET) == [], (
    "self-test failed: _rjem_malloc_conf expect(unsafe_code) in main.rs must be exempted"
)
assert len(scan_source("crates/vardictrust/src/cli.rs", MALLOC_CONF_SNIPPET)) == 1, (
    "self-test failed: _rjem_malloc_conf expect(unsafe_code) outside main.rs must fail"
)
assert len(
    scan_source(
        "crates/vardictrust/src/main.rs",
        '#[cfg_attr(not(feature = "profile-cpu"), expect(unsafe_code, reason = "sneaky"))]\n'
        "static OTHER: u8 = 0;\n",
    )
) == 1, "self-test failed: expect(unsafe_code) on non-_rjem_malloc_conf item in main.rs must fail"

# 2) Occurrences inside comments, string/raw-string literals, reason strings, and #[cfg(test)] blocks MUST be ignored
NEGATIVE_SNIPPET = """
// #[allow(clippy::indexing_slicing)]
/* #[allow(clippy::all)] /* nested #[allow(warnings)] */ #[allow(panic)] */
const S1: &str = "#[allow(clippy::indexing_slicing)] \\" #[allow(todo)]";
const S2: &str = r#"#[expect(clippy::indexing_slicing)]"#;
const S3: &[u8] = br##"#[allow(clippy::restriction)] #[allow(unreachable)]"##;
const Q: char = '"';
#[expect(deprecated, reason = "mentions panic, unsafe_code, and indexing_slicing inside reason literal")]
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

mod inline_inner_cfg_tests {
    #![cfg(test)]
    #[allow(clippy::indexing_slicing)]
    #[allow(clippy::unwrap_used)]
    fn inside_inner_cfg_test_mod() {}
}
"""
assert scan_source("synthetic_clean.rs", NEGATIVE_SNIPPET) == [], (
    "self-test failed: comments, strings, raw strings, and #[cfg(test)] blocks must be ignored"
)
assert scan_source(
    "synthetic_file_cfg_test.rs",
    "#![cfg(test)]\n#[allow(clippy::indexing_slicing)]\nfn file_level_test_only() {}\n",
) == [], "self-test failed: file-level #![cfg(test)] module must be ignored"

# 3) Code BEFORE and AFTER #[cfg(test)] / inline mod { #![cfg(test)] } must still be scanned!
POST_TEST_SNIPPET = (
    "#[allow(clippy::panic)]\nfn before_inner_cfg_mod() {}\n"
    + NEGATIVE_SNIPPET
    + "\n#[allow(clippy::indexing_slicing)]\nfn after_test_mod() {}\n"
)
assert len(scan_source("synthetic_post.rs", POST_TEST_SNIPPET)) == 2, (
    "self-test failed: non-test code surrounding #[cfg(test)] and inline #![cfg(test)] scopes must still be scanned"
)

# 4) Self-test header verifiers against commented-out or inactive-cfg forbid/deny attributes
assert verify_lib_rs_header(
    "// #![forbid(unsafe_code)]\n#![deny(clippy::unwrap_used)]\n"
), "self-test failed: commented-out forbid(unsafe_code) and missing deny lints must fail"
assert verify_lib_rs_header(
    "#![forbid(unsafe_code)]\n"
    "// clippy::print_stderr\n"
    "#![cfg_attr(not(test), deny(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::unreachable, clippy::todo, clippy::unimplemented, clippy::dbg_macro, clippy::print_stdout, clippy::string_slice, clippy::indexing_slicing))]\n"
), "self-test failed: missing active deny(clippy::print_stderr) must fail"
assert verify_lib_rs_header(
    "#![forbid(unsafe_code)]\n"
    "#![cfg_attr(test, deny(" + ", ".join(REQUIRED_DENY_LINTS) + "))]\n"
), "self-test failed: cfg_attr(test, deny(...)) must fail verify_lib_rs_header"
assert verify_cli_rs_header("// #![forbid(unsafe_code)]\n"), (
    "self-test failed: commented-out cli.rs forbid(unsafe_code) must fail"
)
assert verify_profile_rs_header("#![deny(unsafe_op_in_unsafe_fn)]\n"), (
    "self-test failed: incomplete profile.rs header must fail"
)

# 5) Verify active header attributes across lib.rs, main.rs, cli.rs, and profile.rs
header_checks = (
    ("crates/vardictrust-core/src/lib.rs", verify_lib_rs_header),
    ("crates/vardictrust/src/main.rs", verify_main_rs_header),
    ("crates/vardictrust/src/cli.rs", verify_cli_rs_header),
    ("crates/vardictrust/src/profile.rs", verify_profile_rs_header),
)
header_failed = False
for rel_hdr, verifier in header_checks:
    hdr_path = pathlib.Path(rel_hdr)
    errs = verifier(hdr_path.read_text(encoding="utf-8"))
    if errs:
        header_failed = True
        for err in errs:
            print(f"FAIL: {rel_hdr}: {err}", file=sys.stderr)
if header_failed:
    sys.exit(1)

all_hits: list[tuple[str, int, str]] = []
for src_dir in (
    pathlib.Path("crates/vardictrust-core/src"),
    pathlib.Path("crates/vardictrust/src"),
):
    for rs_path in sorted(src_dir.rglob("*.rs")):
        rel = rs_path.as_posix()
        all_hits.extend(scan_source(rel, rs_path.read_text(encoding="utf-8")))

if all_hits:
    print(
        "FAIL: unauthorized #[allow(...)] / #[expect(...)] bypass on ratcheted lint in crates/{vardictrust-core,vardictrust}/src/ (zero-allow policy):",
        file=sys.stderr,
    )
    for rel, line_no, attr_src in all_hits:
        compact = " ".join(attr_src.split())
        print(f"  {rel}:{line_no}: {compact}", file=sys.stderr)
    sys.exit(1)
PY

if [[ "$SCAN_ONLY" -eq 1 ]]; then
    echo "PASS: lint_ratchet.sh --scan-only (11-lint + unsafe_code zero-allow header & source scan clean)"
    exit 0
fi

cargo clippy --workspace --all-targets --no-default-features -- -D warnings
cargo clippy --workspace --all-targets -- -D warnings
cargo clippy --workspace --all-targets --all-features -- -D warnings

echo "PASS: lint_ratchet.sh (11-lint + unsafe_code zero-allow ratchet clean across --no-default-features, default, and --all-features)"
