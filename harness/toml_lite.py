#!/usr/bin/env python3
"""
harness/toml_lite.py — Zero-dependency TOML loader for Python 3.8+.

Uses stdlib `tomllib` on Python 3.11+ and falls back to a pure-stdlib
recursive-descent subset parser on Python 3.8–3.10 covering sections,
arrays of tables, dotted keys, inline tables, arrays, booleans, integers,
floats, single/double-quoted strings, and multiline triple-quoted strings.
"""

from __future__ import annotations

import pathlib
import re
from typing import Any

try:
    import tomllib as _stdlib_tomllib  # Python 3.11+
except ImportError:
    _stdlib_tomllib = None  # type: ignore[assignment]

_TOML_TOKEN_RE = re.compile(
    r"""
    (?P<ws>[ \t\r\n]+)
    | (?P<comment>\#[^\n]*)
    | (?P<ml_basic>\"\"\"[\s\S]*?\"\"\")
    | (?P<ml_lit>'''[\s\S]*?''')
    | (?P<basic>\"(?:\\.|[^\"\\\n])*\")
    | (?P<lit>'[^'\n]*')
    | (?P<lbrack>\[)
    | (?P<rbrack>\])
    | (?P<lbrace>\{)
    | (?P<rbrace>\})
    | (?P<eq>=)
    | (?P<comma>,)
    | (?P<float>[+-]?(?:\d+\.\d+|\d+[eE][+-]?\d+))
    | (?P<dot>\.)
    | (?P<bare>[+-]?[A-Za-z0-9_-]+)
    """,
    re.VERBOSE,
)


def _unescape_toml_basic(raw_inner: str) -> str:
    escapes = {
        "n": "\n",
        "r": "\r",
        "t": "\t",
        "b": "\b",
        "f": "\f",
        '"': '"',
        "\\": "\\",
    }
    out: list[str] = []
    i = 0
    n = len(raw_inner)
    while i < n:
        ch = raw_inner[i]
        if ch == "\\" and i + 1 < n:
            nxt = raw_inner[i + 1]
            if nxt in escapes:
                out.append(escapes[nxt])
                i += 2
                continue
            if nxt == "u" and i + 5 < n:
                out.append(chr(int(raw_inner[i + 2 : i + 6], 16)))
                i += 6
                continue
            if nxt == "U" and i + 9 < n:
                out.append(chr(int(raw_inner[i + 2 : i + 10], 16)))
                i += 10
                continue
        out.append(ch)
        i += 1
    return "".join(out)


def _tokenize_toml(text: str) -> list[tuple[str, str]]:
    tokens: list[tuple[str, str]] = []
    pos = 0
    n = len(text)
    while pos < n:
        m = _TOML_TOKEN_RE.match(text, pos)
        if not m:
            snippet = text[pos : pos + 40]
            raise ValueError(f"Unexpected TOML character at offset {pos}: {snippet!r}")
        kind = m.lastgroup
        assert kind is not None
        val = m.group(kind)
        pos = m.end()
        if kind in ("ws", "comment"):
            continue
        tokens.append((kind, val))
    return tokens


def _fallback_loads(text: str) -> dict[str, Any]:
    """Pure-stdlib TOML subset parser for Python 3.8–3.10."""
    tokens = _tokenize_toml(text)
    idx = 0
    total = len(tokens)
    root: dict[str, Any] = {}
    current_table: dict[str, Any] = root

    def peek() -> tuple[str | None, str | None]:
        return tokens[idx] if idx < total else (None, None)

    def advance() -> tuple[str, str]:
        nonlocal idx
        tok = tokens[idx]
        idx += 1
        return tok

    def parse_key_segment() -> str:
        kind, val = advance()
        if kind == "bare":
            return val
        if kind == "basic":
            return _unescape_toml_basic(val[1:-1])
        if kind == "lit":
            return val[1:-1]
        raise ValueError(f"Expected TOML key segment, got {kind}:{val!r}")

    def parse_dotted_key() -> list[str]:
        parts = [parse_key_segment()]
        while peek()[0] == "dot":
            advance()
            parts.append(parse_key_segment())
        return parts

    def descend_table(base: dict[str, Any], parts: list[str]) -> dict[str, Any]:
        cur = base
        for part in parts:
            if part not in cur:
                cur[part] = {}
            nxt = cur[part]
            if isinstance(nxt, list):
                cur = nxt[-1]
            elif isinstance(nxt, dict):
                cur = nxt
            else:
                raise ValueError(f"Cannot descend into non-table key {part!r}")
        return cur

    def set_dotted(target: dict[str, Any], parts: list[str], value: Any) -> None:
        cur = target
        for part in parts[:-1]:
            if part not in cur:
                cur[part] = {}
            cur = cur[part]
        cur[parts[-1]] = value

    def parse_value() -> Any:
        kind, val = peek()
        assert val is not None or kind is None
        if kind == "ml_basic" and val is not None:
            advance()
            inner = val[3:-3]
            if inner.startswith("\r\n"):
                inner = inner[2:]
            elif inner.startswith("\n"):
                inner = inner[1:]
            return _unescape_toml_basic(inner)
        if kind == "ml_lit" and val is not None:
            advance()
            inner = val[3:-3]
            if inner.startswith("\r\n"):
                inner = inner[2:]
            elif inner.startswith("\n"):
                inner = inner[1:]
            return inner
        if kind == "basic" and val is not None:
            advance()
            return _unescape_toml_basic(val[1:-1])
        if kind == "lit" and val is not None:
            advance()
            return val[1:-1]
        if kind == "float" and val is not None:
            advance()
            return float(val)
        if kind == "bare" and val is not None:
            advance()
            if val == "true":
                return True
            if val == "false":
                return False
            if re.fullmatch(r"[+-]?\d+", val):
                return int(val)
            raise ValueError(f"Unsupported bare TOML value {val!r}")
        if kind == "lbrack":
            advance()
            arr = []
            while peek()[0] != "rbrack":
                if peek()[0] is None:
                    raise ValueError("Unterminated TOML array")
                arr.append(parse_value())
                if peek()[0] == "comma":
                    advance()
            advance()
            return arr
        if kind == "lbrace":
            advance()
            tbl: dict[str, Any] = {}
            while peek()[0] != "rbrace":
                if peek()[0] is None:
                    raise ValueError("Unterminated TOML inline table")
                k_parts = parse_dotted_key()
                eq_k, _ = advance()
                if eq_k != "eq":
                    raise ValueError(f"Expected '=' in inline table, got {eq_k}")
                v = parse_value()
                set_dotted(tbl, k_parts, v)
                if peek()[0] == "comma":
                    advance()
            advance()
            return tbl
        raise ValueError(f"Unexpected token while parsing TOML value: {kind}:{val!r}")

    while idx < total:
        kind, _val = peek()
        if kind == "lbrack":
            advance()
            if peek()[0] == "lbrack":
                advance()
                h_parts = parse_dotted_key()
                for _ in range(2):
                    rk, _ = advance()
                    if rk != "rbrack":
                        raise ValueError(f"Expected ']]' closing array-of-tables header, got {rk}")
                parent = descend_table(root, h_parts[:-1])
                last = h_parts[-1]
                if last not in parent:
                    parent[last] = []
                new_entry: dict[str, Any] = {}
                parent[last].append(new_entry)
                current_table = new_entry
            else:
                h_parts = parse_dotted_key()
                rk, _ = advance()
                if rk != "rbrack":
                    raise ValueError(f"Expected ']' closing section header, got {rk}")
                current_table = descend_table(root, h_parts)
        else:
            k_parts = parse_dotted_key()
            eq_k, _ = advance()
            if eq_k != "eq":
                raise ValueError(f"Expected '=' after key {k_parts}, got {eq_k}")
            val_obj = parse_value()
            set_dotted(current_table, k_parts, val_obj)

    return root


def loads(text: str) -> dict[str, Any]:
    """Parse a TOML string using stdlib `tomllib` (3.11+) or `_fallback_loads` (3.8–3.10)."""
    if _stdlib_tomllib is not None:
        return _stdlib_tomllib.loads(text)
    return _fallback_loads(text)


def load_path(path: str | pathlib.Path) -> dict[str, Any]:
    """Read and parse a UTF-8 TOML file from disk."""
    return loads(pathlib.Path(path).read_text(encoding="utf-8"))


def flatten_nav_entries(nav_list: list[Any]) -> list[tuple[str, str]]:
    """Yield `(title, md_rel_path)` leaf pairs in declaration order from parsed Zensical nav."""
    leaves: list[tuple[str, str]] = []
    for entry in nav_list:
        if isinstance(entry, str):
            leaves.append(("", entry))
        elif isinstance(entry, dict):
            for title, val in entry.items():
                if isinstance(val, str):
                    leaves.append((title, val))
                elif isinstance(val, list):
                    leaves.extend(flatten_nav_entries(val))
                else:
                    raise ValueError(f"Unexpected nav value for {title!r}: {val!r}")
        else:
            raise ValueError(f"Unexpected nav entry: {entry!r}")
    return leaves
