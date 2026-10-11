#!/usr/bin/env python3
"""
harness/tests/test_docs_ui.py

Automated UI and document integrity test suite for VarDictRust documentation.

1. Static Markdown Structure Validation (pure stdlib):
   - Table column consistency across header, delimiter, and rows.
   - Tag hygiene: non-HTML angle brackets (e.g. `<dup-N>`) must be backtick-escaped.
   - Markdown link format and anchor validity.

2. Browser UI Rendering with Playwright (dynamic import):
   - Loads GitHub-Flavored Markdown rendering of key documents (ROADMAP.md, README.md, etc.)
     into headless Chromium.
   - Verifies table DOM structure, row/column rendering, and absence of broken HTML nodes.
   - Verifies external link structures.
"""

import html
import importlib
import json
import os
import pathlib
import re
import subprocess
import sys
import unittest

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent.parent

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
        # Exclude git internals, cargo target, conda environments, scratch, and transient dev directories
        rel = str(path.relative_to(REPO_ROOT))
        if any(rel.startswith(ignore) for ignore in [
            ".git/", ".superpowers/", ".agents/", "target/", "scratch/", "VarDictJava/", "node_modules/"
        ]):
            continue
        md_files.append(path)
    return sorted(md_files)


class TestMarkdownStaticHygiene(unittest.TestCase):
    """Test static formatting, table alignment, and escaping across all markdown documents."""

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
                        # Validate completed table
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
            # Split by | excluding escaped \|
            cells = [c.strip() for c in re.split(r"(?<!\\)\|", line_str)]
            # Strip empty leading and trailing cells from outer pipes
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

            # Remove code blocks and inline code spans before inspecting raw text
            no_code_blocks = re.sub(r"```[\s\S]*?```", "", content)
            no_inline_code = re.sub(r"`[^`]*`", "", no_code_blocks)

            for match in angle_pattern.finditer(no_inline_code):
                tag = match.group(1).lower()
                if tag not in ALLOWED_HTML_TAGS:
                    # Check if it looks like an email or autolink
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

        # Pattern: \text{..._...} where _ is not escaped with backslash
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


class TestPlaywrightBrowserRendering(unittest.TestCase):
    """Render documents with Playwright headless Chromium to verify real browser DOM layout."""

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
        if not self.has_playwright:
            self.skipTest("playwright not installed on host; skipping browser rendering test.")

    def _render_gfm_html(self, markdown_text):
        """Render markdown text to GitHub-Flavored HTML using gh api or fallback."""
        try:
            res = subprocess.run(
                ["gh", "api", "/markdown", "-f", f"text={markdown_text}", "-f", "mode=gfm"],
                capture_output=True,
                text=True,
                check=True,
                timeout=10
            )
            return res.stdout
        except (subprocess.SubprocessError, FileNotFoundError):
            # Fallback wrapper for offline/isolated execution
            escaped = html.escape(markdown_text)
            return f"<pre>{escaped}</pre>"

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
            page = browser.new_page()
            page.set_content(page_content)

            # 1. Verify H1 header
            h1_text = page.inner_text("h1")
            self.assertIn("VarDictRust Roadmap", h1_text)

            # 2. Verify all 3 tables exist
            tables = page.query_selector_all("table")
            self.assertGreaterEqual(len(tables), 3, "ROADMAP.md should contain at least 3 tables.")

            # Table 1: v1.1.0 candidates (4 columns, >= 10 rows)
            table1 = tables[0]
            headers1 = [th.inner_text().strip() for th in table1.query_selector_all("th")]
            self.assertEqual(len(headers1), 4)
            rows1 = table1.query_selector_all("tbody tr")
            self.assertGreaterEqual(len(rows1), 10)

            # Table 2: Needs a specification first (4 columns, >= 7 rows)
            table2 = tables[1]
            headers2 = [th.inner_text().strip() for th in table2.query_selector_all("th")]
            self.assertEqual(len(headers2), 4)
            rows2 = table2.query_selector_all("tbody tr")
            self.assertGreaterEqual(len(rows2), 7)

            # Table 3: Open upstream PRs reviewed (3 columns, >= 4 rows)
            table3 = tables[2]
            headers3 = [th.inner_text().strip() for th in table3.query_selector_all("th")]
            self.assertEqual(len(headers3), 3)
            rows3 = table3.query_selector_all("tbody tr")
            self.assertGreaterEqual(len(rows3), 4)

            # 3. Verify links point to AstraZeneca-NGS issues/PRs or docs
            links = page.query_selector_all("a")
            self.assertGreater(len(links), 15, "Should have links for issues, PRs, and docs.")
            for link in links:
                href = link.get_attribute("href")
                self.assertIsNotNone(href)
                self.assertTrue(
                    href.startswith("https://github.com/AstraZeneca-NGS/VarDictJava/")
                    or href.startswith("docs/")
                    or href.startswith("#"),
                    f"Unexpected link URL: {href}"
                )

            # 4. Capture screenshot for verification
            screenshot_path = "/tmp/roadmap_browser_rendered.png"
            page.screenshot(path=screenshot_path, full_page=True)
            self.assertTrue(os.path.exists(screenshot_path))

            browser.close()

    def test_readme_browser_rendering(self):
        """Verify README.md renders cleanly with benchmark tables in headless Chromium."""
        readme_path = REPO_ROOT / "README.md"
        with open(readme_path, "r", encoding="utf-8") as f:
            content = f.read()

        rendered_html = self._render_gfm_html(content[:8000])  # First 8KB covers hook, badges, benchmarks
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
            page = browser.new_page()
            page.set_content(page_content)

            h1_text = page.inner_text("h1")
            self.assertIn("VarDictRust", h1_text)

            tables = page.query_selector_all("table")
            self.assertGreaterEqual(len(tables), 1, "README.md should render empirical benchmark table.")

            browser.close()

    def test_omnibus_validation_rendering_and_katex(self):
        """Verify omnibus-empirical-validation.md renders cleanly with zero KaTeX errors."""
        omnibus_path = REPO_ROOT / "docs" / "validation" / "omnibus-empirical-validation.md"
        with open(omnibus_path, "r", encoding="utf-8") as f:
            content = f.read()

        # Render first 25KB covering executive summary, master roster, NIST ROC math, and derivations
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
            page = browser.new_page()
            page.set_content(page_content)

            # Wait for KaTeX auto-render to process math
            try:
                page.wait_for_selector(".katex", timeout=5000)
            except Exception:
                pass

            # 1. Assert zero KaTeX parse error boxes
            katex_errors = page.query_selector_all(".katex-error")
            error_texts = [e.inner_text().strip() for e in katex_errors]
            self.assertEqual(
                error_texts, [],
                f"KaTeX parse error boxes detected in rendered DOM: {error_texts}"
            )

            # 2. Assert no KaTeX underscore error in page text
            content_html = page.content()
            self.assertNotIn("allowed only in math mode", content_html)

            # 3. Assert math elements were rendered
            katex_nodes = page.query_selector_all(".katex")
            self.assertGreater(len(katex_nodes), 5, "Should render mathematical expressions via KaTeX.")

            # 4. Save verification screenshot
            screenshot_path = "/tmp/omnibus_consolidated_browser_rendered.png"
            page.screenshot(path=screenshot_path, full_page=True)
            self.assertTrue(os.path.exists(screenshot_path))

            browser.close()


if __name__ == "__main__":
    unittest.main()
