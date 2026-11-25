#!/usr/bin/env python3
"""
harness/tests/test_docs_ui.py

Automated UI, Zensical site configuration, built-HTML, and browser rendering test suite
for VarDictRust documentation (100% Python 3.8+ standard library; optional dynamic
Playwright import for headless Chromium DOM verification).

1. Static Markdown & Zensical Configuration Hygiene (pure stdlib):
   - Native zensical.toml schema, 32-page nav completeness, not_in_nav omitted-files guard,
     strict link/anchor validation, MathJax/Mermaid extension settings, and staging exclusions.
   - ReadTheDocs v2 build.commands single-quoting/schema and GitHub Actions docs.yml invariants.
   - Table column consistency, non-HTML angle-bracket escaping, KaTeX/MathJax hygiene,
     wiki links/anchors, and repository relative Markdown links.

2. Built Zensical Site HTML Verification (pure stdlib html.parser):
   - Verifies built site/ pages, navigation links, in-page and cross-page #anchor targets,
     Mermaid custom fences, staged SVG figures/PDFs, and MathJax configuration scripts.
   - Fails closed when VARDICTRUST_REQUIRE_BUILT_SITE=1 and site/index.html is missing.

3. Browser UI Rendering with Playwright (headless Chromium):
   - Fails closed when VARDICTRUST_REQUIRE_PLAYWRIGHT=1 and Playwright is unavailable.
   - Uses a pure-stdlib local disk-caching CDN route handler under .cache/docs_ui_cdn
     (zero hollow stubs) and a pure-stdlib semantic GFM HTML fallback when gh is unauthenticated.
   - Serves built site/ over loopback ThreadingHTTPServer(("127.0.0.1", 0), ...) and verifies
     /, /ROADMAP/, /paper/manuscript/, and /docs/validation/omnibus-empirical-validation/
     for navigation, DOM tables, MathJax (mjx-container, zero mjx-merror), Mermaid diagrams,
     and SVG figures (naturalWidth > 0).
"""

import fnmatch
import functools
import hashlib
import html
from html.parser import HTMLParser
import http.server
import importlib
import importlib.util
import json
import mimetypes
import os
import pathlib
import re
import subprocess
import sys
import tempfile
import threading
import unittest
import urllib.error
import urllib.parse
import urllib.request

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
if str(REPO_ROOT / "harness") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "harness"))

from toml_lite import (  # noqa: E402
    _fallback_loads,
    flatten_nav_entries as _flatten_nav_entries,
    loads as _parse_zensical_toml,
)

CDN_CACHE_DIR = REPO_ROOT / ".cache" / "docs_ui_cdn"

EXCLUDED_DIRS = (
    ".git/",
    ".superpowers/",
    "docs/superpowers/",
    ".docs_staged/",
    "site/",
    ".cache/",
    "target/",
    "scratch/",
    "VarDictJava/",
    "node_modules/",
    ".venv/",
    ".venv_docs/",
    ".venv_rtd/",
    "venv/",
    ".pytest_cache/",
    ".worktrees/",
    "benchmarks/data/",
)

# Standard permitted HTML tags in Markdown
ALLOWED_HTML_TAGS = frozenset([
    "a", "abbr", "b", "blockquote", "br", "code", "dd", "del", "details",
    "div", "dl", "dt", "em", "h1", "h2", "h3", "h4", "h5", "h6", "hr",
    "i", "img", "ins", "kbd", "li", "ol", "p", "pre", "q", "rp", "rt",
    "ruby", "s", "samp", "small", "span", "strong", "sub", "summary",
    "sup", "table", "tbody", "td", "tfoot", "th", "thead", "tr", "ul",
    "var", "wbr"
])


def find_all_markdown_files():
    """Discover all tracked markdown files in the repository."""
    md_files = []
    for path in REPO_ROOT.rglob("*.md"):
        rel = path.relative_to(REPO_ROOT).as_posix()
        if any(rel.startswith(ignore) for ignore in EXCLUDED_DIRS):
            continue
        md_files.append(path)
    return sorted(md_files)


# ---------------------------------------------------------------------------
# Pure-stdlib GFM-to-HTML fallback renderer (used when gh api is unavailable)
# ---------------------------------------------------------------------------

def _heading_to_slug_gfm(heading_line):
    """Convert a markdown heading line to its GitHub anchor slug."""
    text = re.sub(r"^#{1,6}[\t ]+", "", heading_line.strip())
    text = re.sub(r"[\t ]+#+[\t ]*$", "", text)
    text = re.sub(r"<[^>]+>", "", text)
    text = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", text)
    text = text.replace("`", "")
    text = re.sub(r"[\x00-\x1f]", "", text)
    text = text.lower()
    text = re.sub(r"[^\w\s-]", "", text)
    text = text.replace("\t", " ")
    return text.replace(" ", "-")


def _render_gfm_inline(text):
    """Render inline GFM (code spans, math spans, images, links, bold/italic) cleanly."""
    placeholders = []

    def stash(fragment):
        key = f"\x00PH{len(placeholders)}\x00"
        placeholders.append((key, fragment))
        return key

    # 1. Stash inline code spans (`...`) first so angle brackets, links, and math inside backticks are untouched
    text = re.sub(
        r"`([^`\n]+)`",
        lambda m: stash(f"<code>{html.escape(m.group(1))}</code>"),
        text,
    )

    # 2. Stash display/inline math spans ($$...$$ and $...$) with HTML-safe <, >, &
    text = re.sub(
        r"\$\$(.+?)\$\$",
        lambda m: stash(f"$${html.escape(m.group(1), quote=False)}$$"),
        text,
    )
    text = re.sub(
        r"(?<!\$)\$(?!\s)([^$\n]+?)(?<!\s)\$(?!\$)",
        lambda m: stash(f"${html.escape(m.group(1), quote=False)}$"),
        text,
    )

    # 3. Stash explicit safe raw HTML anchor tags <a id="..."></a>
    text = re.sub(
        r'<a\s+(?:id|name)=["\']([A-Za-z0-9_\-]+)["\']\s*>\s*</a>',
        lambda m: stash(f'<a id="{html.escape(m.group(1), quote=True)}"></a>'),
        text,
    )

    # 4. Stash Markdown images ![alt](src) and links [label](url) linearly
    def convert_links_and_images(src):
        out = []
        i = 0
        n = len(src)
        while i < n:
            is_img = src[i] == "!" and i + 1 < n and src[i + 1] == "["
            if is_img or src[i] == "[":
                lb = i + 1 if is_img else i
                depth = 1
                j = lb + 1
                while j < n and depth > 0:
                    if src[j] == "\\":
                        j += 2
                        continue
                    if src[j] == "[":
                        depth += 1
                    elif src[j] == "]":
                        depth -= 1
                    j += 1
                if depth == 0 and j < n and src[j] == "(":
                    k = j + 1
                    p_depth = 1
                    while k < n and p_depth > 0:
                        if src[k] == "\\":
                            k += 2
                            continue
                        if src[k] == "(":
                            p_depth += 1
                        elif src[k] == ")":
                            p_depth -= 1
                        k += 1
                    if p_depth == 0:
                        label_raw = src[lb + 1 : j - 1]
                        url_full = src[j + 1 : k - 1].strip()
                        parts = url_full.split(None, 1)
                        href = parts[0] if parts else ""
                        href_esc = html.escape(href, quote=True)
                        if is_img:
                            alt_esc = html.escape(label_raw, quote=True)
                            out.append(stash(f'<img alt="{alt_esc}" src="{href_esc}">'))
                        else:
                            label_html = convert_links_and_images(label_raw)
                            label_html = html.escape(label_html, quote=False)
                            out.append(stash(f'<a href="{href_esc}">{label_html}</a>'))
                        i = k
                        continue
            out.append(src[i])
            i += 1
        return "".join(out)

    text = convert_links_and_images(text)

    # 5. Escape remaining plain text HTML entities
    text = html.escape(text, quote=False)

    # 6. Bold / italic emphasis
    text = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", text)
    text = re.sub(r"(?<!\*)\*([^*\n]+)\*(?!\*)", r"<em>\1</em>", text)

    # 7. Restore stashed HTML fragments in reverse order (handles nested code inside links)
    for key, fragment in reversed(placeholders):
        text = text.replace(key, fragment)

    return text


def _split_gfm_table_cells(line_str):
    cells = [c.strip() for c in re.split(r"(?<!\\)\|", line_str.strip())]
    if cells and cells[0] == "":
        cells.pop(0)
    if cells and cells[-1] == "":
        cells.pop(-1)
    return [c.replace(r"\|", "|") for c in cells]


def _is_gfm_delimiter_row(cells):
    return bool(cells) and all(re.fullmatch(r":?-+:?", c.strip()) is not None for c in cells)


def _render_gfm_html_fallback(markdown_text):
    """Pure-stdlib semantic GFM-to-HTML renderer (headings, tables, lists, math, links, code)."""
    lines = markdown_text.splitlines()
    out = []
    i = 0
    n = len(lines)

    while i < n:
        line = lines[i]
        stripped = line.strip()

        # Blank line
        if not stripped:
            i += 1
            continue

        # Fenced code block (``` or ~~~)
        fence_m = re.match(r"^(`{3,}|~{3,})(.*)$", stripped)
        if fence_m:
            fence_marker = fence_m.group(1)
            lang = fence_m.group(2).strip().split()[0] if fence_m.group(2).strip() else ""
            code_lines = []
            i += 1
            while i < n and not lines[i].strip().startswith(fence_marker):
                code_lines.append(lines[i])
                i += 1
            if i < n:
                i += 1
            code_escaped = html.escape("\n".join(code_lines))
            if lang == "mermaid":
                out.append(f'<pre class="mermaid"><code>{code_escaped}</code></pre>')
            elif lang:
                out.append(
                    f'<pre><code class="language-{html.escape(lang, quote=True)}">{code_escaped}</code></pre>'
                )
            else:
                out.append(f"<pre><code>{code_escaped}</code></pre>")
            continue

        # Multi-line or single-line $$ display math block
        if stripped.startswith("$$"):
            if stripped.endswith("$$") and len(stripped) > 4:
                math_inner = html.escape(stripped[2:-2], quote=False)
                out.append(f'<div class="math-display">$${math_inner}$$</div>')
                i += 1
                continue
            math_lines = [line[line.find("$$") + 2 :]]
            i += 1
            while i < n:
                cur = lines[i]
                if "$$" in cur:
                    end_idx = cur.rfind("$$")
                    math_lines.append(cur[:end_idx])
                    i += 1
                    break
                math_lines.append(cur)
                i += 1
            math_inner = html.escape("\n".join(math_lines).strip(), quote=False)
            out.append(f'<div class="math-display">$${math_inner}$$</div>')
            continue

        # Horizontal rule
        if re.fullmatch(r"(?:-\s*){3,}|(?:\*\s*){3,}|(?:_\s*){3,}", stripped):
            out.append("<hr>")
            i += 1
            continue

        # ATX Heading (# .. ######)
        h_m = re.match(r"^(#{1,6})[\t ]+(.*?)(?:[\t ]+#+[\t ]*)?$", stripped)
        if h_m:
            level = len(h_m.group(1))
            slug = _heading_to_slug_gfm(stripped)
            inner_html = _render_gfm_inline(h_m.group(2))
            out.append(f'<h{level} id="{html.escape(slug, quote=True)}">{inner_html}</h{level}>')
            i += 1
            continue

        # GFM Pipe Table
        if stripped.startswith("|") and stripped.endswith("|") and i + 1 < n:
            next_stripped = lines[i + 1].strip()
            if next_stripped.startswith("|") and next_stripped.endswith("|"):
                header_cells = _split_gfm_table_cells(stripped)
                delim_cells = _split_gfm_table_cells(next_stripped)
                if len(header_cells) == len(delim_cells) and _is_gfm_delimiter_row(delim_cells):
                    ths = "".join(f"<th>{_render_gfm_inline(c)}</th>" for c in header_cells)
                    i += 2
                    body_rows = []
                    while i < n:
                        row_s = lines[i].strip()
                        if not (row_s.startswith("|") and row_s.endswith("|")):
                            break
                        row_cells = _split_gfm_table_cells(row_s)
                        tds = "".join(f"<td>{_render_gfm_inline(c)}</td>" for c in row_cells)
                        body_rows.append(f"<tr>{tds}</tr>")
                        i += 1
                    tbody_html = f"<tbody>{''.join(body_rows)}</tbody>" if body_rows else "<tbody></tbody>"
                    out.append(f"<table><thead><tr>{ths}</tr></thead>{tbody_html}</table>")
                    continue

        # Unordered list
        if re.match(r"^[-*+]\s+", stripped):
            items = []
            while i < n and re.match(r"^\s*[-*+]\s+", lines[i]):
                item_txt = re.sub(r"^\s*[-*+]\s+", "", lines[i])
                items.append(f"<li>{_render_gfm_inline(item_txt)}</li>")
                i += 1
            out.append(f"<ul>{''.join(items)}</ul>")
            continue

        # Ordered list
        if re.match(r"^\d+\.\s+", stripped):
            items = []
            while i < n and re.match(r"^\s*\d+\.\s+", lines[i]):
                item_txt = re.sub(r"^\s*\d+\.\s+", "", lines[i])
                items.append(f"<li>{_render_gfm_inline(item_txt)}</li>")
                i += 1
            out.append(f"<ol>{''.join(items)}</ol>")
            continue

        # Blockquote
        if stripped.startswith(">"):
            bq_lines = []
            while i < n and lines[i].strip().startswith(">"):
                bq_lines.append(re.sub(r"^\s*>\s?", "", lines[i]))
                i += 1
            out.append(f"<blockquote>{_render_gfm_inline(' '.join(bq_lines))}</blockquote>")
            continue

        # Paragraph — only stop at block-starter candidates after consuming >= 1 line
        # so unrecognized block candidates (e.g. orphan '| a | b |' lines) always advance i.
        para_lines = []
        while i < n:
            cur_s = lines[i].strip()
            if not cur_s:
                break
            if para_lines and (
                cur_s.startswith(("```", "~~~", "$$", "#", ">"))
                or (cur_s.startswith("|") and cur_s.endswith("|"))
                or re.match(r"^[-*+]\s+", cur_s)
                or re.match(r"^\d+\.\s+", cur_s)
                or re.fullmatch(r"(?:-\s*){3,}|(?:\*\s*){3,}|(?:_\s*){3,}", cur_s)
            ):
                break
            para_lines.append(cur_s)
            i += 1
        if para_lines:
            out.append(f"<p>{_render_gfm_inline(' '.join(para_lines))}</p>")

    return "\n".join(out)


# ---------------------------------------------------------------------------
# Pure-stdlib HTML inspector for Zensical built site/
# ---------------------------------------------------------------------------

class _BuiltHtmlInspector(HTMLParser):
    """Extract element IDs/names, anchor hrefs, script srcs, img srcs, and CSS classes."""

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.ids = set()
        self.hrefs = []
        self.nav_hrefs = []
        self.scripts = []
        self.images = []
        self.mermaid_blocks = 0
        self.arithmatex_nodes = 0
        self._nav_depth = 0

    def handle_starttag(self, tag, attrs):
        attr_map = dict(attrs)
        el_id = attr_map.get("id")
        if el_id:
            self.ids.add(el_id)
        if tag == "a" and attr_map.get("name"):
            self.ids.add(attr_map["name"])

        if tag == "nav":
            self._nav_depth += 1

        classes = set((attr_map.get("class") or "").split())
        if "mermaid" in classes and tag in ("pre", "div"):
            self.mermaid_blocks += 1
        if "arithmatex" in classes:
            self.arithmatex_nodes += 1

        if tag == "a" and "href" in attr_map:
            href = attr_map["href"]
            self.hrefs.append(href)
            if self._nav_depth > 0:
                self.nav_hrefs.append(href)
        elif tag == "script" and attr_map.get("src"):
            self.scripts.append(attr_map["src"])
        elif tag == "img" and attr_map.get("src"):
            self.images.append(attr_map["src"])

    def handle_endtag(self, tag):
        if tag == "nav" and self._nav_depth > 0:
            self._nav_depth -= 1


def _nav_md_to_site_html(site_dir, md_rel):
    """Map a staged Markdown path in zensical.toml nav to its built HTML path under site/."""
    if md_rel in ("README.md", "index.md"):
        return site_dir / "index.html"
    p = pathlib.PurePosixPath(md_rel)
    if p.name in ("README.md", "index.md"):
        return site_dir.joinpath(*p.parent.parts, "index.html")
    return site_dir.joinpath(*p.parent.parts, p.stem, "index.html")


# ---------------------------------------------------------------------------
# 1. Static Markdown & Zensical Configuration Hygiene
# ---------------------------------------------------------------------------

class TestMarkdownStaticHygiene(unittest.TestCase):
    """Test static formatting, Zensical config, table alignment, and escaping across docs."""

    def test_zensical_config_not_in_nav_and_validation(self):
        """Verify native zensical.toml schema, 32-page nav + omitted-files guard, and staging exclusions."""
        self.assertFalse(
            (REPO_ROOT / "mkdocs.yml").exists(),
            "Legacy mkdocs.yml must be removed once zensical.toml cutover is complete",
        )
        zensical_path = REPO_ROOT / "zensical.toml"
        self.assertTrue(zensical_path.is_file(), "Missing native zensical.toml at repository root")

        zensical_raw = zensical_path.read_text(encoding="utf-8")
        cfg = _parse_zensical_toml(zensical_raw)
        self.assertEqual(_fallback_loads(zensical_raw), cfg)
        self.assertIn("project", cfg, "zensical.toml must define [project] table")
        proj = cfg["project"]

        self.assertEqual(proj.get("site_name"), "VarDictRust")
        self.assertEqual(proj.get("docs_dir"), ".docs_staged")
        self.assertEqual(proj.get("site_dir"), "site")
        self.assertEqual(proj.get("site_url"), "https://vardictrust.readthedocs.io/")
        self.assertEqual(proj.get("repo_url"), "https://github.com/palakpsheth/VarDictRust")
        self.assertEqual(proj.get("repo_name"), "palakpsheth/VarDictRust")
        self.assertTrue(
            isinstance(proj.get("copyright"), str) and proj["copyright"].strip(),
            "zensical.toml [project] must set non-empty copyright",
        )

        # Verify strict link & anchor validation
        val_cfg = proj.get("validation", {})
        self.assertIs(val_cfg.get("invalid_links"), True, "zensical.toml must set validation.invalid_links = true")
        self.assertIs(
            val_cfg.get("invalid_link_anchors"),
            True,
            "zensical.toml must set validation.invalid_link_anchors = true",
        )

        # Verify theme features & dual palettes (default + slate)
        theme = proj.get("theme", {})
        features = theme.get("features", [])
        for req_feat in (
            "navigation.tabs",
            "navigation.sections",
            "navigation.top",
            "search.suggest",
            "search.highlight",
            "content.code.copy",
        ):
            self.assertIn(req_feat, features, f"zensical.toml theme.features missing {req_feat!r}")
        palettes = theme.get("palette", [])
        schemes = {entry.get("scheme"): entry for entry in palettes if isinstance(entry, dict)}
        self.assertIn("default", schemes, "zensical.toml missing default light palette")
        self.assertIn("slate", schemes, "zensical.toml missing slate dark palette")
        for scheme_name in ("default", "slate"):
            self.assertEqual(schemes[scheme_name].get("primary"), "teal")
            self.assertEqual(schemes[scheme_name].get("accent"), "deep orange")

        # Verify Markdown extensions (MathJax arithmatex, Mermaid superfence, tables, attr_list, etc.)
        mdx = proj.get("markdown_extensions", {})
        for req_ext in ("tables", "admonition", "attr_list", "md_in_html", "toc"):
            self.assertIn(req_ext, mdx, f"zensical.toml [project.markdown_extensions] missing {req_ext!r}")
        self.assertIs(mdx["toc"].get("permalink"), True)
        pymdownx = mdx.get("pymdownx", {})
        for req_pym in ("details", "highlight", "inlinehilite", "superfences", "tabbed", "arithmatex"):
            self.assertIn(req_pym, pymdownx, f"zensical.toml markdown_extensions.pymdownx missing {req_pym!r}")
        self.assertIs(pymdownx["arithmatex"].get("generic"), True)
        self.assertIs(pymdownx["tabbed"].get("alternate_style"), True)
        fences = pymdownx["superfences"].get("custom_fences", [])
        self.assertTrue(
            any(
                isinstance(f, dict)
                and f.get("name") == "mermaid"
                and f.get("class") == "mermaid"
                and f.get("format") == "pymdownx.superfences.fence_code_format"
                for f in fences
            ),
            "zensical.toml pymdownx.superfences.custom_fences must configure mermaid fence_code_format",
        )

        # Verify MathJax local config + CDN script ordering
        extra_js = proj.get("extra_javascript", [])
        self.assertIn("javascripts/mathjax.js", extra_js)
        self.assertIn("https://unpkg.com/mathjax@3.2.2/es5/tex-mml-chtml.js", extra_js)
        self.assertLess(
            extra_js.index("javascripts/mathjax.js"),
            extra_js.index("https://unpkg.com/mathjax@3.2.2/es5/tex-mml-chtml.js"),
            "javascripts/mathjax.js must load before tex-mml-chtml.js so window.MathJax is configured first",
        )
        mathjax_js_path = REPO_ROOT / "docs" / "javascripts" / "mathjax.js"
        self.assertTrue(mathjax_js_path.is_file(), "Missing docs/javascripts/mathjax.js")
        mathjax_js_text = mathjax_js_path.read_text(encoding="utf-8")
        self.assertIn("window.MathJax", mathjax_js_text)
        self.assertIn("processEscapes: true", mathjax_js_text)
        self.assertIn("processEnvironments: true", mathjax_js_text)
        self.assertIn('typeof document$ !== "undefined"', mathjax_js_text)

        # Verify not_in_nav exclusions & complete nav coverage (omitted-files guard)
        not_in_nav_raw = proj.get("not_in_nav", "")
        self.assertTrue(
            isinstance(not_in_nav_raw, str) and not_in_nav_raw.strip(),
            "zensical.toml [project] must declare multiline not_in_nav string",
        )
        not_in_nav_patterns = [ln.strip() for ln in not_in_nav_raw.splitlines() if ln.strip()]
        for expected_pat in (
            "benchmarks/gcp-prereg.md",
            "docs/validation/gcp-study-protocol.md",
            "docs/validation/scale-recalibration-proposal.md",
            "docs/validation/provenance/task-*.md",
            "paper/tables/*.md",
        ):
            self.assertIn(expected_pat, not_in_nav_patterns)

        nav_leaves = _flatten_nav_entries(proj.get("nav", []))
        self.assertGreaterEqual(len(nav_leaves), 27, "zensical.toml nav must contain all >=27 documentation pages")
        nav_paths = [p for _, p in nav_leaves]
        self.assertIn(
            ("JVM & R Kernel Oracles", "fixtures/jvm-oracles/README.md"),
            nav_leaves,
            "zensical.toml Architecture nav must include fixtures/jvm-oracles/README.md",
        )
        self.assertIn(("Migration Guide", "docs/migration.md"), nav_leaves)

        for title, md_rel in nav_leaves:
            if md_rel.startswith("community/"):
                src_path = REPO_ROOT / ".github" / pathlib.PurePosixPath(md_rel).name
            else:
                src_path = REPO_ROOT / md_rel
            self.assertTrue(
                src_path.is_file(),
                f"zensical.toml nav entry {title!r} -> {md_rel!r} does not exist at {src_path}",
            )

        # Omitted-files guard: every Markdown file staged by prepare_docs_site.py must appear in nav or not_in_nav
        prep_mod_spec = importlib.util.spec_from_file_location(
            "prepare_docs_site", REPO_ROOT / "harness" / "prepare_docs_site.py"
        )
        prep_mod = importlib.util.module_from_spec(prep_mod_spec)
        prep_mod_spec.loader.exec_module(prep_mod)

        expected_staged_rels = set()
        for p in REPO_ROOT.rglob("*.md"):
            rel_posix = p.relative_to(REPO_ROOT).as_posix()
            if not any(rel_posix.startswith(ign) for ign in prep_mod.IGNORED_PREFIXES):
                expected_staged_rels.add(rel_posix)
        gh_dir = REPO_ROOT / ".github"
        if gh_dir.is_dir():
            for gh_md in gh_dir.glob("*.md"):
                expected_staged_rels.add(f"community/{gh_md.name}")

        nav_path_set = set(nav_paths)
        omitted = sorted(
            rel_p
            for rel_p in expected_staged_rels
            if rel_p not in nav_path_set
            and not any(fnmatch.fnmatch(rel_p, pat) for pat in not_in_nav_patterns)
        )
        self.assertEqual(
            omitted,
            [],
            f"Staged Markdown files omitted from both zensical.toml nav and not_in_nav: {omitted}",
        )

        # Verify prepare_docs_site.py exclusions, javascripts staging, and docstring hygiene
        prep_text = (REPO_ROOT / "harness" / "prepare_docs_site.py").read_text(encoding="utf-8")
        for excluded in ('".venv/"', '".venv_docs/"', '".venv_rtd/"', '"venv/"', '".pytest_cache/"', '".cache/"'):
            self.assertIn(
                excluded,
                prep_text,
                f"prepare_docs_site.py IGNORED_PREFIXES must exclude {excluded}",
            )
        self.assertIn('"javascripts"', prep_text, "prepare_docs_site.py must stage docs/javascripts/")
        self.assertNotIn("MkDocs", prep_text, "prepare_docs_site.py contains stale 'MkDocs' docstring/comment")

    @staticmethod
    def _is_valid_zensical_pin(reqs_text):
        if re.search(r"^\s*mkdocs(?:-material)?\b", reqs_text, re.MULTILINE):
            return False
        m = re.search(r"^zensical\s*([><=!~,\d.\s]+?)\s*(?:#.*)?$", reqs_text, re.MULTILINE)
        if not m:
            return False
        clauses = [c.strip() for c in m.group(1).split(",") if c.strip()]
        if len(clauses) != 2:
            return False
        has_lower = False
        has_upper = False
        for clause in clauses:
            low_m = re.fullmatch(r">=\s*0\.1\.(\d+)", clause)
            if low_m:
                if int(low_m.group(1)) >= 1:
                    has_lower = True
                continue
            if re.fullmatch(r"<\s*0\.2(?:\.0)?", clause):
                has_upper = True
        return has_lower and has_upper

    def test_zensical_requirements_and_dependabot_config(self):
        """Verify docs/requirements.txt pins zensical>=0.1.1,<0.2.0 and dependabot.yml ignores zensical >= 0.2.0."""
        reqs_text = (REPO_ROOT / "docs" / "requirements.txt").read_text(encoding="utf-8")
        self.assertTrue(
            self._is_valid_zensical_pin(reqs_text),
            "docs/requirements.txt must pin zensical>=0.1.1,<0.2.0 and remove legacy mkdocs packages",
        )
        for compatible in (
            "zensical>=0.1.1,<0.2.0\n",
            "zensical>=0.1.2,<0.2.0\n",
            "zensical>=0.1.5,<0.2\n",
            "zensical<0.2.0,>=0.1.1\n",
        ):
            self.assertTrue(
                self._is_valid_zensical_pin(compatible),
                f"Validator must accept compatible Zensical 0.1.x pin {compatible.strip()!r}",
            )
        for invalid in (
            "zensical>=0.1.0,<0.2.0\n",
            "zensical>=0.1.1\n",
            "zensical>=0.1.1,<0.3.0\n",
            "zensical>=0.2.0,<0.3.0\n",
            "zensical>=0.1.1,<0.2.0\nmkdocs>=1.6.1,<2.0\n",
        ):
            self.assertFalse(
                self._is_valid_zensical_pin(invalid),
                f"Validator must reject out-of-range or legacy spec {invalid.strip()!r}",
            )

        dependabot_text = (REPO_ROOT / ".github" / "dependabot.yml").read_text(encoding="utf-8")
        pip_m = re.search(
            r'^  - package-ecosystem:\s*"pip"\n(?:(?!^  - package-ecosystem:)[^\n]*\n?)*',
            dependabot_text,
            re.MULTILINE,
        )
        self.assertIsNotNone(pip_m, ".github/dependabot.yml missing pip ecosystem block")
        pip_block = pip_m.group(0)
        self.assertRegex(
            pip_block,
            re.compile(r'^    directory:\s*"/docs"\s*$', re.MULTILINE),
            "pip ecosystem block in .github/dependabot.yml must target directory '/docs'",
        )
        self.assertNotIn("mkdocs", pip_block, ".github/dependabot.yml /docs block must not reference mkdocs")
        ignore_m = re.search(
            r'^    ignore:\s*\n'
            r'      - dependency-name:\s*"zensical"\s*\n'
            r'        versions:\s*\[">= 0\.2\.0"\]\s*$',
            pip_block,
            re.MULTILINE,
        )
        self.assertIsNotNone(
            ignore_m,
            '.github/dependabot.yml /docs pip entry must declare 4-space ignore block for zensical >= 0.2.0',
        )
        groups_m = re.search(r"^    groups:\s*$", pip_block, re.MULTILINE)
        self.assertIsNotNone(groups_m, ".github/dependabot.yml /docs pip entry missing 4-space groups block")
        self.assertLess(
            ignore_m.start(),
            groups_m.start(),
            "ignore block in .github/dependabot.yml /docs pip entry must not be nested under groups",
        )

    def test_readthedocs_and_docs_workflow_schema(self):
        """Validate .readthedocs.yaml v2 build.commands quoting/schema and .github/workflows/docs.yml gates."""
        rtd_path = REPO_ROOT / ".readthedocs.yaml"
        self.assertTrue(rtd_path.is_file(), "Missing .readthedocs.yaml")
        rtd_text = rtd_path.read_text(encoding="utf-8")

        self.assertRegex(rtd_text, re.compile(r"^version:\s*2\s*$", re.MULTILINE))
        self.assertRegex(rtd_text, re.compile(r"^  os:\s*ubuntu-24\.04\s*$", re.MULTILINE))
        self.assertRegex(rtd_text, re.compile(r'^    python:\s*"3\.12"\s*$', re.MULTILINE))

        # Mutually-exclusive ReadTheDocs v2 keys forbidden alongside build.commands
        for forbidden_pat, desc in (
            (r"^mkdocs\s*:", "top-level mkdocs:"),
            (r"^sphinx\s*:", "top-level sphinx:"),
            (r"^python\s*:", "top-level python:"),
            (r"^  jobs\s*:", "build.jobs:"),
        ):
            self.assertNotRegex(
                rtd_text,
                re.compile(forbidden_pat, re.MULTILINE),
                f".readthedocs.yaml must omit {desc} when using build.commands",
            )

        cmds_m = re.search(r"^  commands:\s*\n((?:    - [^\n]+\n?)+)", rtd_text, re.MULTILINE)
        self.assertIsNotNone(cmds_m, ".readthedocs.yaml must define build.commands list")
        cmd_lines = [ln.strip() for ln in cmds_m.group(1).splitlines() if ln.strip()]
        self.assertGreaterEqual(len(cmd_lines), 6)
        unquoted = []
        for cmd_ln in cmd_lines:
            self.assertTrue(cmd_ln.startswith("- "), f"Malformed command entry: {cmd_ln!r}")
            scalar = cmd_ln[2:].strip()
            self.assertTrue(
                len(scalar) >= 2 and scalar.startswith("'") and scalar.endswith("'"),
                f".readthedocs.yaml build.commands entry must be outer single-quoted: {cmd_ln!r}",
            )
            unquoted.append(scalar[1:-1].replace("''", "'"))

        joined_cmds = "\n".join(unquoted)
        self.assertNotIn(".venv_rtd", joined_cmds, ".readthedocs.yaml must use .venv_docs fallback, not .venv_rtd")
        for required_fragment in (
            'python3 -m venv "${READTHEDOCS_VIRTUALENV_PATH:-.venv_docs}"',
            '"${READTHEDOCS_VIRTUALENV_PATH:-.venv_docs}/bin/python" -m pip install --upgrade pip',
            '"${READTHEDOCS_VIRTUALENV_PATH:-.venv_docs}/bin/python" -m pip install -r docs/requirements.txt',
            "harness/prepare_docs_site.py",
            '"${READTHEDOCS_VIRTUALENV_PATH:-.venv_docs}/bin/zensical" build --clean --strict',
            'mkdir -p "$READTHEDOCS_OUTPUT/html"',
            'cp -a site/. "$READTHEDOCS_OUTPUT/html/"',
        ):
            self.assertIn(
                required_fragment,
                joined_cmds,
                f".readthedocs.yaml build.commands missing {required_fragment!r}",
            )

        try:
            yaml_mod = importlib.import_module("yaml")
        except ImportError:
            yaml_mod = None
        if yaml_mod is not None:
            parsed_rtd = yaml_mod.safe_load(rtd_text)
            self.assertEqual(parsed_rtd.get("version"), 2)
            self.assertNotIn("mkdocs", parsed_rtd)
            self.assertNotIn("sphinx", parsed_rtd)
            self.assertNotIn("python", parsed_rtd)
            self.assertNotIn("jobs", parsed_rtd.get("build", {}))
            self.assertEqual(parsed_rtd["build"]["commands"], unquoted)

        # Validate .github/workflows/docs.yml line 18 job name, unfiltered triggers, Playwright, and cleanup order
        docs_wf = REPO_ROOT / ".github" / "workflows" / "docs.yml"
        wf_lines = docs_wf.read_text(encoding="utf-8").splitlines()
        self.assertGreaterEqual(len(wf_lines), 18)
        self.assertEqual(
            wf_lines[17],
            "    name: MkDocs Strict Build Gate",
            "docs.yml line 18 must preserve exact branch-protection status check 'name: MkDocs Strict Build Gate'",
        )
        wf_text = "\n".join(wf_lines)
        self.assertNotIn("paths:", wf_text, "docs.yml must keep unfiltered push/pull_request triggers")
        self.assertNotIn("paths-ignore:", wf_text)
        self.assertIn("branches: [main]", wf_text)
        self.assertIn("pull_request:", wf_text)
        self.assertIn("workflow_dispatch:", wf_text)
        self.assertIn("playwright install --with-deps chromium", wf_text)
        self.assertIn("zensical build --clean --strict", wf_text)
        self.assertIn("GH_TOKEN: ${{ github.token }}", wf_text)

        pos_ui = wf_text.find(
            "VARDICTRUST_REQUIRE_BUILT_SITE=1 VARDICTRUST_REQUIRE_PLAYWRIGHT=1 python3 harness/tests/test_docs_ui.py"
        )
        pos_clean = wf_text.find("rm -rf .docs_staged site .cache")
        pos_status = wf_text.find("git status --porcelain")
        self.assertNotEqual(pos_ui, -1, "docs.yml must run test_docs_ui.py with built-site + Playwright flags")
        self.assertNotEqual(pos_clean, -1, "docs.yml must clean .docs_staged site .cache")
        self.assertNotEqual(pos_status, -1, "docs.yml must verify git status --porcelain after cleanup")
        self.assertLess(pos_ui, pos_clean, "test_docs_ui.py must run BEFORE rm -rf .docs_staged site .cache")
        self.assertLess(pos_clean, pos_status, "rm -rf .docs_staged site .cache must run BEFORE git status --porcelain")

    def test_gfm_fallback_renders_semantic_html_without_pre_wrapper(self):
        """Verify pure-stdlib _render_gfm_html_fallback produces real headings, tables, links, code, and math."""
        sample_md = (
            "# Sample Roadmap\n\n"
            "Inline math $\\varepsilon = 10^{-6}$ and token `<dup-N>` with [Issue #1](https://example.com/1).\n\n"
            "$$\n"
            "F_1 = \\frac{2 \\cdot TP}{2 \\cdot TP + FP + FN}\n"
            "$$\n\n"
            "| ID | Flag |\n"
            "| :--- | :--- |\n"
            "| `#341` | [`--fix a`](docs/upstream-bugs.md#a) |\n"
        )
        rendered = _render_gfm_html_fallback(sample_md)
        self.assertFalse(rendered.strip().startswith("<pre>"), "Fallback must not wrap document in <pre>")
        self.assertIn('<h1 id="sample-roadmap">Sample Roadmap</h1>', rendered)
        self.assertIn("<table><thead><tr><th>ID</th><th>Flag</th></tr></thead><tbody>", rendered)
        self.assertIn("<code>&lt;dup-N&gt;</code>", rendered)
        self.assertIn('<a href="https://example.com/1">Issue #1</a>', rendered)
        self.assertIn('<a href="docs/upstream-bugs.md#a"><code>--fix a</code></a>', rendered)
        self.assertIn("$\\varepsilon = 10^{-6}$", rendered)
        self.assertIn("$$F_1 =", rendered)

        # Orphan pipe-delimited lines that fail table detection must terminate and render as <p>
        orphan_md = "| orphan | pipe |\nFollowing line\n\n| A | B |\n|:--|--:|\n| 1 | 2 |\n"
        orphan_rendered = _render_gfm_html_fallback(orphan_md)
        self.assertIn("<p>| orphan | pipe | Following line</p>", orphan_rendered)
        self.assertIn("<table><thead><tr><th>A</th><th>B</th></tr></thead><tbody><tr><td>1</td><td>2</td></tr></tbody></table>", orphan_rendered)

        # Full omnibus-empirical-validation.md (including :-- 2-hyphen delimiter roster at lines 107-108) must render cleanly
        full_omnibus = (REPO_ROOT / "docs" / "validation" / "omnibus-empirical-validation.md").read_text(
            encoding="utf-8"
        )
        full_omnibus_html = _render_gfm_html_fallback(full_omnibus)
        self.assertIn("<table><thead><tr><th>#</th><th>Cohort / Specimen Family</th>", full_omnibus_html)
        self.assertFalse(full_omnibus_html.strip().startswith("<pre>"))

    def test_toml_parsers_2d_arrays_globs_and_non_cdn_route_abort(self):
        """Verify stdlib TOML parsers handle 2D arrays, [glob] lines in not_in_nav, escapes, and non-CDN route aborts."""
        sample_toml = (
            '[project]\n'
            'site_name = "VarDictRust"\n'
            'inline = [["\\\\(", "\\\\)"], ["\\\\[", "\\\\]"]]\n'
            'single_2d = [["a"]]\n'
            'not_in_nav = """\n'
            '[d]raft/*.md\n'
            'paper/tables/*.md\n'
            '"""\n'
            'nav = [\n'
            '  { "Overview" = "README.md" },\n'
            '  { "Section" = [\n'
            '    { "Quoted \\"Title\\" & \\\\Slash" = "docs/migration.md" },\n'
            '  ] },\n'
            ']\n'
            '\n'
            '[[project.theme.palette]]\n'
            'scheme = "default"\n'
        )
        parsed = _parse_zensical_toml(sample_toml)
        self.assertEqual(_fallback_loads(sample_toml), parsed)
        self.assertEqual(parsed["project"]["inline"], [["\\(", "\\)"], ["\\[", "\\]"]])
        self.assertEqual(parsed["project"]["single_2d"], [["a"]])
        self.assertEqual(parsed["project"]["theme"]["palette"], [{"scheme": "default"}])

        crc_path = (
            REPO_ROOT / ".agents" / "skills" / "vardictrust-release" / "scripts" / "check_release_consistency.py"
        )
        crc_spec = importlib.util.spec_from_file_location("check_release_consistency_test", crc_path)
        crc_mod = importlib.util.module_from_spec(crc_spec)
        crc_spec.loader.exec_module(crc_mod)
        nav_entries = crc_mod.parse_zensical_nav_entries(sample_toml)
        self.assertEqual(
            nav_entries,
            [
                ("Overview", "README.md"),
                ('Quoted "Title" & \\Slash', "docs/migration.md"),
            ],
        )

        # Verify _install_cdn_disk_cache_route aborts non-CDN external hosts (fonts.googleapis.com, shields.io, api.github.com)
        captured_handler = {}

        class _FakeContext:
            def route(self, pattern, fn):
                captured_handler["pattern"] = pattern
                captured_handler["fn"] = fn

        class _FakeRequest:
            def __init__(self, url):
                self.url = url

        class _FakeRoute:
            def __init__(self, url):
                self.request = _FakeRequest(url)
                self.action = None

            def continue_(self):
                self.action = "continue"

            def abort(self):
                self.action = "abort"

            def fulfill(self, **_kwargs):
                self.action = "fulfill"

        pw_case = TestPlaywrightBrowserRendering("test_readme_browser_rendering")
        pw_case._install_cdn_disk_cache_route(_FakeContext())
        handler_fn = captured_handler["fn"]

        for allowed_url in ("http://127.0.0.1:8123/ROADMAP/", "http://localhost:8000/", "data:image/svg+xml,<svg/>"):
            r = _FakeRoute(allowed_url)
            handler_fn(r)
            self.assertEqual(r.action, "continue", f"Expected continue_ for local/inline URL {allowed_url}")

        for blocked_url in (
            "https://fonts.googleapis.com/css?family=Inter",
            "https://fonts.gstatic.com/s/inter/v1/font.woff2",
            "https://img.shields.io/readthedocs/vardictrust",
            "https://api.github.com/repos/palakpsheth/VarDictRust",
        ):
            r = _FakeRoute(blocked_url)
            handler_fn(r)
            self.assertEqual(r.action, "abort", f"Expected abort for non-CDN external URL {blocked_url}")

    def test_table_column_balance(self):
        """Verify that every markdown table has consistent column counts across all rows."""
        md_files = find_all_markdown_files()
        self.assertGreater(len(md_files), 10, "Should discover at least 10 markdown documents.")

        table_errors = []
        for md_path in md_files:
            rel_path = md_path.relative_to(REPO_ROOT)
            with open(md_path, "r", encoding="utf-8") as f:
                lines = f.readlines()

            in_code_block = False
            current_table_lines = []
            current_table_start = 0

            for idx, line in enumerate(lines, 1):
                stripped = line.strip()
                if stripped.startswith("```"):
                    in_code_block = not in_code_block
                    continue
                if in_code_block:
                    continue

                if stripped.startswith("|") and stripped.endswith("|"):
                    if not current_table_lines:
                        current_table_start = idx
                    current_table_lines.append((idx, stripped))
                else:
                    if current_table_lines:
                        self._validate_table(rel_path, current_table_start, current_table_lines, table_errors)
                        current_table_lines = []

            if current_table_lines:
                self._validate_table(rel_path, current_table_start, current_table_lines, table_errors)

        self.assertEqual(
            table_errors, [],
            "Table column count mismatches found:\n" + "\n".join(f"  - {e}" for e in table_errors)
        )

    def _validate_table(self, rel_path, start_line, table_lines, errors):
        """Helper to check column count equality within a single table."""
        if len(table_lines) < 2:
            return

        expected_cols = None
        for line_num, line_str in table_lines:
            cells = [c.strip() for c in re.split(r"(?<!\\)\|", line_str)]
            if cells and cells[0] == "":
                cells.pop(0)
            if cells and cells[-1] == "":
                cells.pop(-1)

            col_count = len(cells)
            if expected_cols is None:
                expected_cols = col_count
            elif col_count != expected_cols:
                errors.append(
                    f"{rel_path}:{line_num}: table starting at line {start_line} has "
                    f"{col_count} columns; expected {expected_cols}"
                )

    def test_angle_bracket_escaping_in_tables(self):
        """Ensure non-HTML angle brackets (e.g. <dup-N>, <DUP>) are backtick-escaped."""
        md_files = find_all_markdown_files()
        escaping_errors = []

        angle_pattern = re.compile(r"<([A-Za-z0-9_\-]+)>")

        for md_path in md_files:
            rel_path = md_path.relative_to(REPO_ROOT)
            with open(md_path, "r", encoding="utf-8") as f:
                content = f.read()

            no_code_blocks = re.sub(r"```[\s\S]*?```", "", content)
            no_inline_code = re.sub(r"`[^`]*`", "", no_code_blocks)

            for match in angle_pattern.finditer(no_inline_code):
                tag = match.group(1).lower()
                if tag not in ALLOWED_HTML_TAGS:
                    if "@" in tag or "." in tag:
                        continue
                    escaping_errors.append(
                        f"{rel_path}: unescaped angle bracket token '<{match.group(1)}>' "
                        "should be enclosed in backticks."
                    )

        self.assertEqual(
            escaping_errors, [],
            "Unescaped angle brackets found (can break HTML/GFM rendering):\n"
            + "\n".join(f"  - {e}" for e in escaping_errors)
        )

    def test_katex_math_hygiene(self):
        """Verify that markdown files do not contain unescaped underscores in \\text{} or malformed math."""
        md_files = find_all_markdown_files()
        math_errors = []

        text_underscore_pattern = re.compile(r"\\text\{[^}]*(?<!\\)_[^}]*\}")

        for md_path in md_files:
            rel_path = md_path.relative_to(REPO_ROOT)
            with open(md_path, "r", encoding="utf-8") as f:
                lines = f.readlines()

            in_code_block = False
            for idx, line in enumerate(lines, 1):
                stripped = line.strip()
                if stripped.startswith("```"):
                    in_code_block = not in_code_block
                    continue
                if in_code_block:
                    continue

                for m in text_underscore_pattern.finditer(line):
                    math_errors.append(
                        f"{rel_path}:{idx}: unescaped underscore in '{m.group(0)}' "
                        "(causes KaTeX ParseError: '_' allowed only in math mode; use hyphens or \\_)"
                    )

        self.assertEqual(
            math_errors, [],
            "KaTeX math formatting errors found:\n" + "\n".join(f"  - {e}" for e in math_errors)
        )

    @staticmethod
    def _heading_to_slug(heading_line):
        return _heading_to_slug_gfm(heading_line)

    @staticmethod
    def _extract_all_anchors(content):
        """Extract all valid anchors (explicit HTML IDs/names and GFM heading slugs) from markdown content."""
        anchors = set()
        for m in re.finditer(r'(?:id|name)=["\']([^"\']+)["\']', content, re.IGNORECASE):
            anchors.add(m.group(1))

        for line in content.splitlines():
            line_str = line.strip()
            if re.match(r"^#{1,6}\s+", line_str):
                slug = TestMarkdownStaticHygiene._heading_to_slug(line_str)
                if slug:
                    anchors.add(slug)
        return anchors

    def test_heading_to_slug_gfm_parity(self):
        """Verify _heading_to_slug produces exact GFM slugs for complex headings."""
        cases = [
            (
                "### 10.3 Release Packaging & SLSA Level 3 (`.github/workflows/release.yml`)",
                "103-release-packaging--slsa-level-3-githubworkflowsreleaseyml",
            ),
            (
                "## Level 1: Unit & Component Tests",
                "level-1-unit--component-tests",
            ),
            (
                "## Pillar 9: Multi-Modal Public Cohorts",
                "pillar-9-multi-modal-public-cohorts",
            ),
            (
                "### Heading with Trailing Hashes ###",
                "heading-with-trailing-hashes",
            ),
            (
                "### Section with <code>inline_code</code> and [link](https://example.com)",
                "section-with-inline_code-and-link",
            ),
            (
                "### Clean Title\r\n",
                "clean-title",
            ),
        ]
        for heading, expected_slug in cases:
            self.assertEqual(
                self._heading_to_slug(heading),
                expected_slug,
                f"Slug mismatch for '{heading}'",
            )

    def test_extract_all_anchors(self):
        """Verify _extract_all_anchors captures explicit HTML IDs and ATX heading slugs."""
        sample_doc = (
            "# Main Heading\n\n"
            '<a id="custom-anchor"></a>\n'
            "### Sub-Section & Details\n\n"
            '<div name="named-div">content</div>\n'
        )
        anchors = self._extract_all_anchors(sample_doc)
        self.assertIn("main-heading", anchors)
        self.assertIn("custom-anchor", anchors)
        self.assertIn("sub-section--details", anchors)
        self.assertIn("named-div", anchors)

    def test_wiki_links_and_anchors(self):
        """Verify wiki links and anchors are valid and compatible with GitHub Wiki."""
        wiki_dir = REPO_ROOT / "wiki"
        if not wiki_dir.exists():
            return

        wiki_page_files = list(wiki_dir.glob("*.md"))
        wiki_page_stems = {p.stem: p for p in wiki_page_files}

        link_pattern = re.compile(r"(?<!\!)\[([^\]]+)\]\(([^)\s]+)(?:\s+[\"'][^\"']*[\"'])?\)")
        errors = []

        for wiki_file in sorted(wiki_page_files):
            rel_path = wiki_file.relative_to(REPO_ROOT)
            with open(wiki_file, "r", encoding="utf-8") as f:
                content = f.read()

            no_code_blocks = re.sub(r"```[\s\S]*?```", "", content)
            no_inline_code = re.sub(r"`[^`\n]*`", "CODE", no_code_blocks)

            for match in link_pattern.finditer(no_inline_code):
                label = match.group(1)
                raw_url = match.group(2)
                url = urllib.parse.unquote(raw_url)

                if raw_url.startswith("../") or url.startswith("../"):
                    errors.append(
                        f"{rel_path}: relative link '{raw_url}' for '{label}' will 404 on GitHub Wiki; "
                        "must use absolute GitHub URL or wiki page slug"
                    )
                    continue

                if raw_url.startswith("#"):
                    anchor = raw_url[1:]
                    file_anchors = self._extract_all_anchors(content)
                    if anchor not in file_anchors and urllib.parse.unquote(anchor) not in file_anchors:
                        errors.append(
                            f"{rel_path}: in-page anchor '{raw_url}' for '{label}' not found in {rel_path}"
                        )
                    continue

                if not (raw_url.startswith("http://") or raw_url.startswith("https://") or raw_url.startswith("mailto:")):
                    target_page, _, page_anchor = raw_url.partition("#")
                    target_page_clean = urllib.parse.unquote(target_page)

                    if target_page_clean.endswith(".md"):
                        errors.append(
                            f"{rel_path}: internal wiki link '{raw_url}' for '{label}' ends with '.md' "
                            "which triggers raw file download; use slug without .md"
                        )
                        continue

                    matched_file = None
                    if target_page_clean in wiki_page_stems:
                        matched_file = wiki_page_stems[target_page_clean]
                    else:
                        norm_target = re.sub(r"[\s_-]+", "-", target_page_clean.lower())
                        for stem, path in wiki_page_stems.items():
                            if re.sub(r"[\s_-]+", "-", stem.lower()) == norm_target:
                                matched_file = path
                                break

                    if not matched_file:
                        errors.append(
                            f"{rel_path}: internal wiki link target page '{target_page}' for '{label}' "
                            "does not exist in wiki/"
                        )
                    elif page_anchor:
                        with open(matched_file, "r", encoding="utf-8") as tf:
                            target_wiki_content = tf.read()
                        target_anchors = self._extract_all_anchors(target_wiki_content)
                        if page_anchor not in target_anchors and urllib.parse.unquote(page_anchor) not in target_anchors:
                            errors.append(
                                f"{rel_path}: anchor '#{page_anchor}' not found in wiki page "
                                f"'{matched_file.name}' (link '{raw_url}')"
                            )
                    continue

                if raw_url.startswith("https://github.com/palakpsheth/VarDictRust/blob/main/") or \
                   raw_url.startswith("https://github.com/palakpsheth/VarDictRust/tree/main/"):
                    prefix = "https://github.com/palakpsheth/VarDictRust/blob/main/" if \
                        raw_url.startswith("https://github.com/palakpsheth/VarDictRust/blob/main/") else \
                        "https://github.com/palakpsheth/VarDictRust/tree/main/"
                    target_ref = raw_url[len(prefix):]
                    if "#" in target_ref:
                        subpath, anchor = target_ref.split("#", 1)
                    else:
                        subpath, anchor = target_ref, None

                    subpath_unquoted = urllib.parse.unquote(subpath)
                    target_file = REPO_ROOT / subpath_unquoted
                    if not target_file.exists():
                        errors.append(
                            f"{rel_path}: target file does not exist: '{subpath}' (link '{raw_url}')"
                        )
                    elif anchor:
                        check_file = target_file if target_file.is_file() else (target_file / "README.md")
                        if not check_file.exists():
                            errors.append(
                                f"{rel_path}: target file for anchor does not exist: '{check_file}' (link '{raw_url}')"
                            )
                        else:
                            with open(check_file, "r", encoding="utf-8") as tf:
                                target_content = tf.read()
                            valid_anchors = self._extract_all_anchors(target_content)
                            if anchor not in valid_anchors and urllib.parse.unquote(anchor) not in valid_anchors:
                                errors.append(
                                    f"{rel_path}: anchor '#{anchor}' not found in '{check_file.relative_to(REPO_ROOT)}' (link '{raw_url}')"
                                )
                    continue

                if raw_url.startswith("https://github.com/palakpsheth/VarDictRust") and "#" in raw_url:
                    anchor = raw_url.split("#", 1)[1]
                    readme_path = REPO_ROOT / "README.md"
                    if not readme_path.exists():
                        errors.append(
                            f"{rel_path}: README.md does not exist for anchor check (link '{raw_url}')"
                        )
                    else:
                        with open(readme_path, "r", encoding="utf-8") as rf:
                            readme_content = rf.read()
                        valid_anchors = self._extract_all_anchors(readme_content)
                        if anchor not in valid_anchors and urllib.parse.unquote(anchor) not in valid_anchors:
                            errors.append(
                                f"{rel_path}: anchor '#{anchor}' not found in README.md (link '{raw_url}')"
                            )
                    continue

        self.assertEqual(
            errors, [],
            "Wiki link / anchor hygiene errors found:\n" + "\n".join(f"  - {e}" for e in errors)
        )

    def test_readme_badges_camo_safe(self):
        """Verify README.md badges use Camo-safe endpoints and avoid raw readthedocs.org badge redirects."""
        readme_path = REPO_ROOT / "README.md"
        with open(readme_path, "r", encoding="utf-8") as f:
            content = f.read()
        self.assertNotIn(
            "readthedocs.org/projects/vardictrust/badge/",
            content,
            "README.md must not use raw readthedocs.org/projects/.../badge/ image URLs "
            "(Cloudflare 307/429 causes GitHub Camo HTTP 502); use img.shields.io/readthedocs/vardictrust.",
        )
        self.assertIn(
            "https://img.shields.io/readthedocs/vardictrust",
            content,
            "README.md should use img.shields.io/readthedocs/vardictrust for the Documentation Status badge.",
        )

    def test_repo_markdown_relative_links(self):
        """Verify relative file links and anchors across all tracked repository markdown files."""
        md_files = find_all_markdown_files()
        link_pattern = re.compile(r"(?<!\!)\[([^\]]+)\]\(([^)\s]+)(?:\s+[\"'][^\"']*[\"'])?\)")
        errors = []
        anchor_cache = {}
        tracked_files = set()
        tracked_dirs = set()
        try:
            tracked_out = subprocess.check_output(
                ["git", "-C", str(REPO_ROOT), "ls-files"],
                stderr=subprocess.DEVNULL,
                text=True,
            )
            for line in tracked_out.splitlines():
                if not line.strip():
                    continue
                tf = (REPO_ROOT / line).resolve()
                tracked_files.add(tf)
                parent = tf.parent
                while parent != REPO_ROOT and REPO_ROOT in parent.parents:
                    tracked_dirs.add(parent)
                    parent = parent.parent
            if tracked_files:
                tracked_dirs.add(REPO_ROOT)
        except (OSError, subprocess.CalledProcessError):
            pass

        for md_path in md_files:
            rel_path = md_path.relative_to(REPO_ROOT)
            if str(rel_path).startswith("wiki/"):
                continue

            with open(md_path, "r", encoding="utf-8") as f:
                content = f.read()

            no_code_blocks = re.sub(r"```[\s\S]*?```", "", content)
            no_inline_code = re.sub(r"`[^`\n]*`", "CODE", no_code_blocks)

            for match in link_pattern.finditer(no_inline_code):
                label = match.group(1)
                raw_url = match.group(2)

                if raw_url.startswith(("http://", "https://", "mailto:", "conversation://", "file://")):
                    continue

                if raw_url.startswith("#"):
                    anchor = raw_url[1:]
                    if not anchor or re.match(r"^L\d+(-L\d+)?$", anchor):
                        continue
                    if md_path not in anchor_cache:
                        anchor_cache[md_path] = self._extract_all_anchors(content)
                    valid_anchors = anchor_cache[md_path]
                    if anchor not in valid_anchors and urllib.parse.unquote(anchor) not in valid_anchors:
                        errors.append(
                            f"{rel_path}: broken in-page anchor '{raw_url}' for '{label}'"
                        )
                    continue

                target_part, _, anchor = raw_url.partition("#")
                target_clean = urllib.parse.unquote(target_part)
                if not target_clean:
                    continue

                resolved = (md_path.parent / target_clean).resolve()
                try:
                    resolved.relative_to(REPO_ROOT)
                except ValueError:
                    errors.append(
                        f"{rel_path}: relative link '{raw_url}' escapes repository root ('{resolved}')"
                    )
                    continue
                if not resolved.exists():
                    errors.append(
                        f"{rel_path}: broken relative link '{raw_url}' for '{label}' "
                        f"(target '{resolved}' does not exist)"
                    )
                    continue
                if tracked_files and resolved not in tracked_files and resolved not in tracked_dirs:
                    errors.append(
                        f"{rel_path}: relative link '{raw_url}' for '{label}' points to "
                        f"untracked path '{resolved.relative_to(REPO_ROOT)}'"
                    )
                    continue

                if anchor:
                    if re.match(r"^L\d+(-L\d+)?$", anchor):
                        continue
                    anchor_file = resolved if resolved.is_file() else (resolved / "README.md")
                    if anchor_file.is_file() and anchor_file.suffix.lower() == ".md":
                        if anchor_file not in anchor_cache:
                            with open(anchor_file, "r", encoding="utf-8") as tf:
                                anchor_cache[anchor_file] = self._extract_all_anchors(tf.read())
                        valid_anchors = anchor_cache[anchor_file]
                        if anchor not in valid_anchors and urllib.parse.unquote(anchor) not in valid_anchors:
                            errors.append(
                                f"{rel_path}: broken anchor '#{anchor}' in "
                                f"'{anchor_file.relative_to(REPO_ROOT)}' (link '{raw_url}')"
                            )

        self.assertEqual(
            errors, [],
            "Broken relative markdown links or anchors found:\n" + "\n".join(f"  - {e}" for e in errors),
        )


# ---------------------------------------------------------------------------
# 2. Pure-Stdlib HTML Inspection of Built Zensical Site (site/)
# ---------------------------------------------------------------------------

class TestZensicalBuiltSiteHtml(unittest.TestCase):
    """Validate built Zensical HTML pages, navigation, anchors, Mermaid, SVGs, and MathJax config."""

    def setUp(self):
        self.site_dir = REPO_ROOT / "site"
        self.require_built = os.environ.get("VARDICTRUST_REQUIRE_BUILT_SITE") == "1"
        index_html = self.site_dir / "index.html"
        if not index_html.is_file():
            if self.require_built:
                self.fail(
                    "VARDICTRUST_REQUIRE_BUILT_SITE=1 is set, but site/index.html does not exist. "
                    "Run 'python3 harness/prepare_docs_site.py && zensical build --clean --strict' first."
                )
            self.skipTest(
                "site/index.html not present; set VARDICTRUST_REQUIRE_BUILT_SITE=1 after building to enforce."
            )

    def test_built_site_nav_anchors_mermaid_mathjax_and_assets(self):
        """Verify all nav pages exist in site/, internal links/anchors resolve, and extensions emitted DOM."""
        zensical_path = REPO_ROOT / "zensical.toml"
        self.assertTrue(zensical_path.is_file(), "zensical.toml must exist when validating built site/")
        zensical_raw = zensical_path.read_text(encoding="utf-8")
        cfg = _parse_zensical_toml(zensical_raw)
        self.assertEqual(_fallback_loads(zensical_raw), cfg)
        nav_leaves = _flatten_nav_entries(cfg["project"]["nav"])

        # 1. Verify staged MathJax config & PDF/SVG assets exist in site/
        self.assertTrue(
            (self.site_dir / "javascripts" / "mathjax.js").is_file(),
            "Built site/ missing site/javascripts/mathjax.js",
        )
        for pdf_name in ("manuscript.pdf", "supplement.pdf"):
            if (REPO_ROOT / "paper" / pdf_name).is_file():
                self.assertTrue(
                    (self.site_dir / "paper" / pdf_name).is_file(),
                    f"Built site/ missing staged paper/{pdf_name}",
                )

        # 2. Verify every nav entry built its expected index.html
        for title, md_rel in nav_leaves:
            built_html = _nav_md_to_site_html(self.site_dir, md_rel)
            self.assertTrue(
                built_html.is_file(),
                f"Missing built HTML for nav entry {title!r} ({md_rel}) at {built_html}",
            )

        # 3. Parse all HTML files under site/ (excluding 404.html relative-base quirks where appropriate)
        html_files = sorted(self.site_dir.rglob("*.html"))
        self.assertGreaterEqual(len(html_files), len(nav_leaves))

        inspected = {}
        for hf in html_files:
            inspector = _BuiltHtmlInspector()
            inspector.feed(hf.read_text(encoding="utf-8"))
            inspected[hf.resolve()] = inspector

        # 4. Check root index.html navigation, scripts, and arithmatex
        root_ins = inspected[(self.site_dir / "index.html").resolve()]
        self.assertGreaterEqual(
            len(root_ins.nav_hrefs),
            20,
            "site/index.html must render full navigation links",
        )
        self.assertTrue(
            any("javascripts/mathjax.js" in s for s in root_ins.scripts),
            "site/index.html must include javascripts/mathjax.js script tag",
        )
        self.assertGreater(root_ins.arithmatex_nodes, 0, "site/index.html must contain .arithmatex math spans")

        # 5. Check Mermaid custom fences & SVG figures on manuscript / testing / omnibus pages
        omnibus_ins = inspected[
            (self.site_dir / "docs" / "validation" / "omnibus-empirical-validation" / "index.html").resolve()
        ]
        matrix_ins = inspected[
            (self.site_dir / "docs" / "testing" / "test-matrix-and-ci-audit" / "index.html").resolve()
        ]
        manuscript_html = (self.site_dir / "paper" / "manuscript" / "index.html").resolve()
        manuscript_ins = inspected[manuscript_html]

        self.assertGreater(omnibus_ins.mermaid_blocks, 0, "Omnibus validation HTML missing Mermaid blocks")
        self.assertGreater(matrix_ins.mermaid_blocks, 0, "Test matrix HTML missing Mermaid blocks")
        self.assertGreater(omnibus_ins.arithmatex_nodes, 5, "Omnibus validation HTML missing arithmatex nodes")
        self.assertGreater(manuscript_ins.arithmatex_nodes, 5, "Manuscript HTML missing arithmatex nodes")

        svg_imgs = [src for src in manuscript_ins.images if src.endswith(".svg")]
        self.assertGreaterEqual(len(svg_imgs), 1, "Manuscript HTML must embed SVG figures")
        for svg_src in svg_imgs:
            resolved_svg = (manuscript_html.parent / svg_src).resolve()
            self.assertTrue(resolved_svg.is_file(), f"Missing built SVG figure referenced by manuscript: {svg_src}")

        # 6. Validate all internal links & #anchor targets across built content pages
        site_root_resolved = self.site_dir.resolve()
        link_errors = []
        for hf_resolved, ins in inspected.items():
            if hf_resolved.name == "404.html":
                continue
            rel_hf = hf_resolved.relative_to(site_root_resolved)
            for raw_href in ins.hrefs:
                if not raw_href or raw_href.startswith(
                    ("http://", "https://", "mailto:", "javascript:", "tel:", "data:")
                ):
                    continue
                if raw_href.startswith("#"):
                    anchor = raw_href[1:]
                    if anchor and anchor not in ins.ids and urllib.parse.unquote(anchor) not in ins.ids:
                        link_errors.append(f"{rel_hf}: broken in-page anchor '#{anchor}'")
                    continue

                path_part, _, anchor = raw_href.partition("#")
                path_part = urllib.parse.unquote(path_part.split("?", 1)[0])
                if not path_part:
                    continue
                if path_part.startswith("/"):
                    target_candidate = (self.site_dir / path_part.lstrip("/")).resolve()
                else:
                    target_candidate = (hf_resolved.parent / path_part).resolve()

                try:
                    target_candidate.relative_to(site_root_resolved)
                except ValueError:
                    link_errors.append(f"{rel_hf}: link '{raw_href}' escapes site_dir")
                    continue
                target_file = (
                    (target_candidate / "index.html").resolve()
                    if target_candidate.is_dir()
                    else target_candidate
                )
                if not target_file.is_file():
                    link_errors.append(f"{rel_hf}: broken internal link '{raw_href}' -> {target_file}")
                    continue
                if anchor and target_file in inspected:
                    target_ids = inspected[target_file].ids
                    if anchor not in target_ids and urllib.parse.unquote(anchor) not in target_ids:
                        link_errors.append(
                            f"{rel_hf}: broken cross-page anchor '#{anchor}' in '{raw_href}'"
                        )

        self.assertEqual(
            link_errors,
            [],
            "Broken internal links or #anchors in built Zensical site/:\n"
            + "\n".join(f"  - {e}" for e in link_errors[:50]),
        )


# ---------------------------------------------------------------------------
# 3. Headless Chromium Browser Rendering (Playwright + Disk-Caching CDN)
# ---------------------------------------------------------------------------

class _QuietStaticHandler(http.server.SimpleHTTPRequestHandler):
    """Serve built site/ over 127.0.0.1 without spamming test stderr."""

    def log_message(self, format, *args):  # noqa: A002
        pass


class TestPlaywrightBrowserRendering(unittest.TestCase):
    """Render documents and built Zensical site/ with Playwright headless Chromium."""

    sync_playwright_fn = None
    has_playwright = False

    @classmethod
    def setUpClass(cls):
        try:
            playwright_mod = importlib.import_module("playwright.sync_api")
            cls.sync_playwright_fn = staticmethod(playwright_mod.sync_playwright)
            cls.has_playwright = True
        except ImportError:
            cls.has_playwright = False

    def setUp(self):
        require_pw = os.environ.get("VARDICTRUST_REQUIRE_PLAYWRIGHT") == "1"
        if not self.has_playwright:
            if require_pw:
                self.fail(
                    "VARDICTRUST_REQUIRE_PLAYWRIGHT=1 is set, but 'playwright.sync_api' is not installed."
                )
            self.skipTest("playwright not installed on host; skipping browser rendering test.")

    def _install_cdn_disk_cache_route(self, context):
        """Attach a pure-stdlib disk-caching CDN route handler under .cache/docs_ui_cdn (no hollow stubs)."""
        CDN_CACHE_DIR.mkdir(parents=True, exist_ok=True)
        cdn_hosts = ("unpkg.com", "jsdelivr.net")

        def handle_route(route):
            url = route.request.url
            parsed = urllib.parse.urlsplit(url)
            host = (parsed.hostname or "").lower()

            # Local loopback / inline data / blob requests bypass external handler
            if host in ("", "127.0.0.1", "localhost") or parsed.scheme in ("data", "blob", "about"):
                route.continue_()
                return

            if any(host == h or host.endswith("." + h) for h in cdn_hosts):
                digest = hashlib.sha256(url.encode("utf-8")).hexdigest()
                body_path = CDN_CACHE_DIR / f"{digest}.bin"
                meta_path = CDN_CACHE_DIR / f"{digest}.meta.json"

                if body_path.is_file() and meta_path.is_file():
                    try:
                        meta = json.loads(meta_path.read_text(encoding="utf-8"))
                        body = body_path.read_bytes()
                        if body:
                            route.fulfill(
                                status=int(meta.get("status", 200)),
                                content_type=meta.get("content_type", "application/octet-stream"),
                                body=body,
                            )
                            return
                    except (OSError, ValueError):
                        pass

                try:
                    req = urllib.request.Request(
                        url,
                        headers={"User-Agent": "VarDictRust-DocsUI-CDN-Cache/1.0"},
                    )
                    with urllib.request.urlopen(req, timeout=20) as resp:
                        body = resp.read()
                        status = getattr(resp, "status", 200) or 200
                        ctype = resp.headers.get("Content-Type")
                        if not ctype:
                            guessed, _ = mimetypes.guess_type(parsed.path)
                            ctype = guessed or "application/octet-stream"
                    if status == 200 and body:
                        tmp_body = CDN_CACHE_DIR / f"{digest}.{os.getpid()}.tmp"
                        tmp_meta = CDN_CACHE_DIR / f"{digest}.{os.getpid()}.meta.tmp"
                        tmp_body.write_bytes(body)
                        tmp_meta.write_text(
                            json.dumps({"url": url, "status": status, "content_type": ctype}),
                            encoding="utf-8",
                        )
                        os.replace(tmp_body, body_path)
                        os.replace(tmp_meta, meta_path)
                    route.fulfill(status=status, content_type=ctype, body=body)
                    return
                except (urllib.error.URLError, OSError):
                    # Fail honestly via network abort — never return hollow JS/CSS stubs
                    route.abort()
                    return

            # Block all non-CDN external URLs (fonts.googleapis.com, shields.io, api.github.com, etc.)
            # so offline / packet-dropping firewalls never stall rendering waiting on TCP timeouts.
            route.abort()

        context.route("**/*", handle_route)

    def _render_gfm_html(self, markdown_text):
        """Render markdown text to GitHub-Flavored HTML using gh api or pure-stdlib semantic fallback."""
        try:
            res = subprocess.run(
                ["gh", "api", "/markdown", "-f", f"text={markdown_text}", "-f", "mode=gfm"],
                capture_output=True,
                text=True,
                check=True,
                timeout=10,
            )
            if res.stdout.strip():
                return res.stdout
        except (subprocess.SubprocessError, OSError):
            pass
        return _render_gfm_html_fallback(markdown_text)

    def test_roadmap_browser_rendering(self):
        """Verify ROADMAP.md renders clean tables and links in headless Chromium."""
        roadmap_path = REPO_ROOT / "ROADMAP.md"
        with open(roadmap_path, "r", encoding="utf-8") as f:
            content = f.read()

        rendered_html = self._render_gfm_html(content)
        page_content = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>VarDictRust Roadmap</title>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif; padding: 24px; }}
table {{ border-collapse: collapse; width: 100%; margin: 16px 0; }}
th, td {{ border: 1px solid #d0d7de; padding: 6px 13px; text-align: left; }}
th {{ background-color: #f6f8fa; }}
code {{ font-family: monospace; background-color: rgba(175, 184, 193, 0.2); padding: 0.2em 0.4em; border-radius: 6px; }}
</style>
</head>
<body>
{rendered_html}
</body>
</html>
"""

        with self.sync_playwright_fn() as p:
            browser = p.chromium.launch(headless=True)
            try:
                context = browser.new_context()
                self._install_cdn_disk_cache_route(context)
                page = context.new_page()
                page.set_content(page_content)

                h1_text = page.inner_text("h1")
                self.assertIn("VarDictRust Roadmap", h1_text)

                tables = page.query_selector_all("table")
                self.assertGreaterEqual(len(tables), 4, "ROADMAP.md should contain at least 4 tables.")

                table1 = tables[0]
                headers1 = [th.inner_text().strip() for th in table1.query_selector_all("th")]
                self.assertEqual(len(headers1), 4)
                rows1 = table1.query_selector_all("tbody tr")
                self.assertGreaterEqual(len(rows1), 21)

                table2 = tables[1]
                headers2 = [th.inner_text().strip() for th in table2.query_selector_all("th")]
                self.assertEqual(len(headers2), 4)
                rows2 = table2.query_selector_all("tbody tr")
                self.assertGreaterEqual(len(rows2), 2)

                table3 = tables[2]
                headers3 = [th.inner_text().strip() for th in table3.query_selector_all("th")]
                self.assertEqual(len(headers3), 4)
                rows3 = table3.query_selector_all("tbody tr")
                self.assertGreaterEqual(len(rows3), 6)

                table4 = tables[3]
                headers4 = [th.inner_text().strip() for th in table4.query_selector_all("th")]
                self.assertEqual(len(headers4), 3)
                rows4 = table4.query_selector_all("tbody tr")
                self.assertGreaterEqual(len(rows4), 4)

                allowed_roadmap_href_prefixes = (
                    "https://github.com/AstraZeneca-NGS/VarDictJava/",
                    "https://github.com/AstraZeneca-NGS/VarDict/",
                    "https://github.com/MHH-Bioinformatics-Hematology/vardictcpp/",
                    "https://github.com/MHH-Bioinformatics-Hematology/galaxytools/",
                    "https://github.com/joachimwolff/VarDictJava/",
                    "https://github.com/clintval/vartovcf/",
                    "https://github.com/PolinaBevad/VarDictJava/",
                    "https://github.com/bioconda/bioconda-recipes/",
                    "https://github.com/nf-core/modules/",
                    "https://github.com/palakpsheth/VarDictRust/",
                    "https://github.com/samtools/hts-specs",
                    "https://github.com/samtools/bcftools/",
                    "https://github.com/bcbio/bcbio-nextgen/",
                    "https://github.com/Clinical-Genomics/BALSAMIC/",
                    "https://github.com/msk-access/postprocessing_variant_calls",
                    "https://github.com/mskcc/vcf2maf/",
                    "https://github.com/bioinform/somaticseq",
                    "https://github.com/sigven/pcgr/",
                    "https://github.com/EBIvariation/vcf-validator/",
                    "https://github.com/Kennedy-Lab-UW/VarDictJava/",
                    "https://github.com/hedera-dx/VarDictJava/",
                    "https://github.com/fengqing-dong/VarDictJava",
                    "https://github.com/Schaudge/VarDictJava",
                    "https://github.com/crutching/VarDictJava",
                    "https://github.com/invitae/VarDictJava",
                    "https://github.com/Belay-Diagnostics/VarDictJava",
                    "https://github.com/Crispy13/VarDictRust2",
                    "https://github.com/LeiHaoa/RabbitVar",
                    "https://github.com/roysomak4/VarDict-Go",
                    "https://pmc.ncbi.nlm.nih.gov/articles/PMC",
                    "https://doi.org/10.",
                    "docs/",
                    "#",
                )
                links = [a for a in page.query_selector_all("a") if a.get_attribute("href")]
                self.assertGreater(len(links), 15, "Should have links for issues, PRs, and docs.")
                for link in links:
                    href = link.get_attribute("href")
                    self.assertIsNotNone(href)
                    assert href is not None
                    self.assertTrue(
                        href.startswith(allowed_roadmap_href_prefixes),
                        f"Unexpected link URL: {href}",
                    )

                with tempfile.TemporaryDirectory() as tmpdir:
                    screenshot_path = os.path.join(tmpdir, "roadmap_browser_rendered.png")
                    page.screenshot(path=screenshot_path, full_page=True)
                    self.assertTrue(os.path.exists(screenshot_path))
                context.close()
            finally:
                browser.close()

    def test_readme_browser_rendering(self):
        """Verify README.md renders cleanly with benchmark tables in headless Chromium."""
        readme_path = REPO_ROOT / "README.md"
        with open(readme_path, "r", encoding="utf-8") as f:
            content = f.read()

        rendered_html = self._render_gfm_html(content[:8000])
        page_content = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>VarDictRust README</title>
</head>
<body>
{rendered_html}
</body>
</html>
"""

        with self.sync_playwright_fn() as p:
            browser = p.chromium.launch(headless=True)
            try:
                context = browser.new_context()
                self._install_cdn_disk_cache_route(context)
                page = context.new_page()
                page.set_content(page_content)

                h1_text = page.inner_text("h1")
                self.assertIn("VarDictRust", h1_text)

                tables = page.query_selector_all("table")
                self.assertGreaterEqual(len(tables), 1, "README.md should render empirical benchmark table.")
                context.close()
            finally:
                browser.close()

    def test_omnibus_validation_rendering_and_katex(self):
        """Verify omnibus-empirical-validation.md renders cleanly with zero KaTeX errors."""
        omnibus_path = REPO_ROOT / "docs" / "validation" / "omnibus-empirical-validation.md"
        with open(omnibus_path, "r", encoding="utf-8") as f:
            content = f.read()

        sample_content = content[:25000]
        rendered_html = self._render_gfm_html(sample_content)

        page_content = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>VarDictRust Omnibus Validation</title>
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.8/dist/katex.min.css">
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif; padding: 24px; }}
table {{ border-collapse: collapse; width: 100%; margin: 16px 0; }}
th, td {{ border: 1px solid #d0d7de; padding: 6px 13px; text-align: left; }}
th {{ background-color: #f6f8fa; }}
.katex-error {{ color: #cf222e; background-color: #ffebe9; padding: 2px 4px; border: 1px solid #cf222e; }}
</style>
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.8/dist/katex.min.js"></script>
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.8/dist/contrib/auto-render.min.js"
        onload="renderMathInElement(document.body, {{delimiters: [{{left: '$$', right: '$$', display: true}}, {{left: '$', right: '$', display: false}}], throwOnError: false}});"></script>
</head>
<body>
{rendered_html}
</body>
</html>
"""

        with self.sync_playwright_fn() as p:
            browser = p.chromium.launch(headless=True)
            try:
                context = browser.new_context()
                self._install_cdn_disk_cache_route(context)
                page = context.new_page()
                page.set_content(page_content)

                page.wait_for_selector(".katex", timeout=15000)

                katex_errors = page.query_selector_all(".katex-error")
                error_texts = [e.inner_text().strip() for e in katex_errors]
                self.assertEqual(
                    error_texts, [],
                    f"KaTeX parse error boxes detected in rendered DOM: {error_texts}",
                )

                content_html = page.content()
                self.assertNotIn("allowed only in math mode", content_html)

                katex_nodes = page.query_selector_all(".katex")
                self.assertGreater(len(katex_nodes), 5, "Should render mathematical expressions via KaTeX.")

                with tempfile.TemporaryDirectory() as tmpdir:
                    screenshot_path = os.path.join(tmpdir, "omnibus_consolidated_browser_rendered.png")
                    page.screenshot(path=screenshot_path, full_page=True)
                    self.assertTrue(os.path.exists(screenshot_path))
                context.close()
            finally:
                browser.close()

    def test_zensical_built_site_loopback_rendering(self):
        """Serve built site/ over loopback HTTP and verify 4 key pages in headless Chromium."""
        site_dir = REPO_ROOT / "site"
        require_built = os.environ.get("VARDICTRUST_REQUIRE_BUILT_SITE") == "1"
        if not (site_dir / "index.html").is_file():
            if require_built:
                self.fail(
                    "VARDICTRUST_REQUIRE_BUILT_SITE=1 is set, but site/index.html is missing for Playwright test."
                )
            self.skipTest("site/index.html not present; skipping Zensical loopback Playwright test.")

        handler = functools.partial(_QuietStaticHandler, directory=str(site_dir))
        server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
        port = server.server_address[1]
        srv_thread = threading.Thread(target=server.serve_forever, daemon=True)
        srv_thread.start()

        try:
            base_url = f"http://127.0.0.1:{port}"
            with self.sync_playwright_fn() as p:
                browser = p.chromium.launch(headless=True)
                try:
                    context = browser.new_context(viewport={"width": 1440, "height": 900})
                    self._install_cdn_disk_cache_route(context)
                    page = context.new_page()

                    # Page 1: / (Overview / README) — navigation, tables, MathJax
                    resp = page.goto(f"{base_url}/", wait_until="domcontentloaded")
                    self.assertEqual(resp.status, 200)
                    self.assertIn("VarDictRust", page.inner_text("h1"))
                    self.assertGreaterEqual(len(page.query_selector_all("nav a")), 10)
                    self.assertGreaterEqual(len(page.query_selector_all("table")), 1)
                    page.wait_for_function(
                        "() => document.querySelectorAll('mjx-container').length > 0",
                        timeout=20000,
                    )
                    self.assertGreater(len(page.query_selector_all("mjx-container")), 0)
                    self.assertEqual(len(page.query_selector_all("mjx-merror, [data-mjx-error]")), 0)

                    # Page 2: /ROADMAP/ — navigation, >= 4 DOM tables
                    resp = page.goto(f"{base_url}/ROADMAP/", wait_until="domcontentloaded")
                    self.assertEqual(resp.status, 200)
                    self.assertIn("VarDictRust Roadmap", page.inner_text("h1"))
                    self.assertGreaterEqual(len(page.query_selector_all("nav a")), 10)
                    self.assertGreaterEqual(len(page.query_selector_all("table")), 4)
                    self.assertEqual(len(page.query_selector_all("mjx-merror, [data-mjx-error]")), 0)

                    # Page 3: /paper/manuscript/ — MathJax equations & SVG figures (naturalWidth > 0)
                    resp = page.goto(f"{base_url}/paper/manuscript/", wait_until="domcontentloaded")
                    self.assertEqual(resp.status, 200)
                    self.assertGreaterEqual(len(page.query_selector_all("nav a")), 10)
                    page.wait_for_function(
                        "() => document.querySelectorAll('mjx-container').length >= 5",
                        timeout=20000,
                    )
                    self.assertGreaterEqual(len(page.query_selector_all("mjx-container")), 5)
                    self.assertEqual(len(page.query_selector_all("mjx-merror, [data-mjx-error]")), 0)
                    page.wait_for_function(
                        "() => { const imgs = Array.from(document.querySelectorAll('img[src$=\".svg\"]')); "
                        "return imgs.length > 0 && imgs.every(i => i.complete && i.naturalWidth > 0); }",
                        timeout=15000,
                    )
                    svg_dims = page.evaluate(
                        "() => Array.from(document.querySelectorAll('img[src$=\".svg\"]')).map("
                        "i => ({src: i.getAttribute('src'), naturalWidth: i.naturalWidth}))"
                    )
                    self.assertGreaterEqual(len(svg_dims), 1)
                    for item in svg_dims:
                        self.assertGreater(
                            item["naturalWidth"],
                            0,
                            f"SVG figure failed to load (naturalWidth == 0): {item}",
                        )

                    # Page 4: /docs/validation/omnibus-empirical-validation/ — MathJax, tables, Mermaid diagrams
                    resp = page.goto(
                        f"{base_url}/docs/validation/omnibus-empirical-validation/",
                        wait_until="domcontentloaded",
                    )
                    self.assertEqual(resp.status, 200)
                    self.assertGreaterEqual(len(page.query_selector_all("nav a")), 10)
                    self.assertGreaterEqual(len(page.query_selector_all("table")), 5)
                    page.wait_for_function(
                        "() => document.querySelectorAll('mjx-container').length >= 5",
                        timeout=20000,
                    )
                    self.assertGreaterEqual(len(page.query_selector_all("mjx-container")), 5)
                    self.assertEqual(len(page.query_selector_all("mjx-merror, [data-mjx-error]")), 0)
                    page.wait_for_function(
                        "() => { const m = Array.from(document.querySelectorAll('div.mermaid')); "
                        "return document.querySelectorAll('pre.mermaid').length === 0 && "
                        "m.length > 0 && m.every(el => el.getBoundingClientRect().height > 0); }",
                        timeout=20000,
                    )
                    self.assertEqual(len(page.query_selector_all("pre.mermaid")), 0)
                    mermaid_boxes = page.evaluate(
                        "() => Array.from(document.querySelectorAll('div.mermaid')).map("
                        "el => ({w: el.getBoundingClientRect().width, h: el.getBoundingClientRect().height}))"
                    )
                    self.assertGreaterEqual(len(mermaid_boxes), 1)
                    for box in mermaid_boxes:
                        self.assertGreater(box["w"], 0)
                        self.assertGreater(box["h"], 0)

                    context.close()
                finally:
                    browser.close()
        finally:
            server.shutdown()
            server.server_close()


if __name__ == "__main__":
    unittest.main()
