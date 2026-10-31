#!/usr/bin/env python3
"""Render paper/manuscript.md and paper/supplement.md into typeset PDFs via headless Chrome."""

from __future__ import annotations

import html
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
import markdown


LATEX_SYMBOLS = {
    r"\varepsilon": "&varepsilon;",
    r"\alpha": "&alpha;",
    r"\beta": "&beta;",
    r"\Theta": "&Theta;",
    r"\Delta": "&Delta;",
    r"\le": "&le;",
    r"\leq": "&le;",
    r"\ge": "&ge;",
    r"\geq": "&ge;",
    r"\in": "&isin;",
    r"\perp": "&perp;",
    r"\approx": "&asymp;",
    r"\times": "&times;",
    r"\cdot": "&middot;",
    r"\pm": "&plusmn;",
    r"\to": "&rarr;",
    r"\rightarrow": "&rarr;",
    r"\quad": "&emsp;",
    r"\,": "&#8201;",
    r"\!": "",
    r"\%": "%",
}


def _convert_math_expr(expr: str) -> str:
    """Convert a LaTeX math expression into clean semantic HTML without raw TeX tokens."""
    s = expr.strip()

    # Replace \left and \right delimiters
    s = re.sub(r"\\left\s*([(\[\{|])", r"\1", s)
    s = re.sub(r"\\right\s*([)\]\}'|])", r"\1", s)

    # Replace \bar{x} -> x with combining macron
    s = re.sub(r"\\bar\{([A-Za-z0-9]+)\}", r"<i>\1</i>&#772;", s)

    # Protect \text{...} blocks before further substitutions
    text_spans: list[str] = []

    def _save_text(m: re.Match[str]) -> str:
        idx = len(text_spans)
        inner = html.escape(m.group(1).strip())
        text_spans.append(f'<span class="math-rm">{inner}</span>')
        return f"@@MATHTEXT{idx}@@"

    s = re.sub(r"\\text\{([^}]*)\}", _save_text, s)

    # Replace \frac{a}{b}
    s = re.sub(
        r"\\frac\{([^{}]+)\}\{([^{}]+)\}",
        r'<span class="math-frac"><span class="num">\1</span>/<span class="den">\2</span></span>',
        s,
    )

    # Replace \sum_{sub}^{sup}
    s = re.sub(
        r"\\sum_\{([^{}]+)\}\^\{([^{}]+)\}",
        r'<span class="math-op">&Sigma;</span><sub>\1</sub><sup>\2</sup>',
        s,
    )

    # Replace braced subscripts and superscripts (loop for nested)
    for _ in range(2):
        s = re.sub(r"_\{([^{}]+)\}", r"<sub>\1</sub>", s)
        s = re.sub(r"\^\{([^{}]+)\}", r"<sup>\1</sup>", s)

    # Replace simple single-token subscripts and superscripts
    s = re.sub(r"_([A-Za-z0-9]+)", r"<sub>\1</sub>", s)
    s = re.sub(r"\^([A-Za-z0-9+-]+)", r"<sup>\1</sup>", s)

    # Replace LaTeX symbol macros
    for tex, ent in sorted(LATEX_SYMBOLS.items(), key=lambda kv: -len(kv[0])):
        s = s.replace(tex, ent)

    # Restore protected \text{...} blocks
    for idx, val in enumerate(text_spans):
        s = s.replace(f"@@MATHTEXT{idx}@@", val)

    return s


def preprocess_markdown_for_html(md_text: str, base_dir: Path) -> str:
    """Protect code blocks/spans, convert LaTeX math to HTML, and resolve figure SVG paths."""
    code_blocks: list[str] = []

    def _save_fenced(m: re.Match[str]) -> str:
        idx = len(code_blocks)
        code_blocks.append(m.group(0))
        return f"@@CODEBLOCK{idx}@@"

    # 1. Protect fenced code blocks ```...```
    text = re.sub(r"```[\s\S]*?```", _save_fenced, md_text)

    # 2. Protect inline code spans `...`
    def _save_inline_code(m: re.Match[str]) -> str:
        idx = len(code_blocks)
        code_blocks.append(m.group(0))
        return f"@@CODEBLOCK{idx}@@"

    text = re.sub(r"`[^`\n]+`", _save_inline_code, text)

    # 3. Convert display math $$...$$
    def _replace_display_math(m: re.Match[str]) -> str:
        inner = _convert_math_expr(m.group(1))
        return f'\n\n<div class="math-display">{inner}</div>\n\n'

    text = re.sub(r"\$\$([\s\S]+?)\$\$", _replace_display_math, text)

    # 4. Convert inline math $...$
    def _replace_inline_math(m: re.Match[str]) -> str:
        inner = _convert_math_expr(m.group(1))
        return f'<span class="math-inline">{inner}</span>'

    text = re.sub(r"\$([^$\n]+)\$", _replace_inline_math, text)

    # 5. Convert Markdown image syntax ![alt](path) into rich <figure> with file:// URI
    def _replace_image(m: re.Match[str]) -> str:
        alt = html.escape(m.group(1).strip())
        rel_path = m.group(2).strip()
        img_path = (base_dir / rel_path).resolve()
        uri = img_path.as_uri()
        return (
            f'\n\n<figure class="paper-figure">\n'
            f'  <img src="{uri}" alt="{alt}" />\n'
            f'</figure>\n\n'
        )

    text = re.sub(r"!\[([^\]]*)\]\(([^)]+)\)", _replace_image, text)

    # 6. Sanitize raw '<' / '<=' / '>=' / '->' in Markdown table rows and ensure blank lines before lists
    processed_lines: list[str] = []
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("|") and stripped.endswith("|"):
            # Do not touch separator rows like | :--- | :--- |
            if not re.match(r"^\|\s*[:\-| ]+\|\s*$", stripped):
                line = line.replace("<=", "&le;").replace(">=", "&ge;").replace("->", "&rarr;")
                # Escape any remaining raw '<' that is not an HTML tag we emitted (<span, </span, <sub, </sub, <sup, </sup, <i, </i)
                line = re.sub(r"<(?!(?:  |/?(?:span|sub|sup|i|code|strong|em|a)\b))", "&lt;", line)
        # Ensure blank line before list items (- ... or 1. ...) so Python-Markdown starts a list block
        if re.match(r"^\s*(?:-\s+|\d+\.\s+)", line):
            if processed_lines and processed_lines[-1].strip() != "":
                prev = processed_lines[-1]
                if not re.match(r"^\s*(?:-\s+|\d+\.\s+|\|)", prev):
                    processed_lines.append("")
        processed_lines.append(line)

    text = "\n".join(processed_lines)

    # 7. Restore code blocks and inline code spans
    for idx, block in enumerate(code_blocks):
        text = text.replace(f"@@CODEBLOCK{idx}@@", block)

    return text


PRINT_CSS = """
@page {
  size: letter;
  margin: 13mm 13mm 15mm 13mm;
  @bottom-left {
    content: "VarDictRust v1.0.6 — Empirical Validation & Architectural Manuscript";
    font-family: "Liberation Sans", "DejaVu Sans", Arial, sans-serif;
    font-size: 7.2pt;
    color: #64748b;
  }
  @bottom-right {
    content: "Page " counter(page) " of " counter(pages);
    font-family: "Liberation Sans", "DejaVu Sans", Arial, sans-serif;
    font-size: 7.2pt;
    color: #64748b;
  }
}
* {
  box-sizing: border-box;
}
html, body {
  margin: 0;
  padding: 0;
  font-family: "Liberation Serif", "Times New Roman", Georgia, serif;
  font-size: 9.2pt;
  line-height: 1.42;
  color: #111827;
  background: #ffffff;
}
h1 {
  font-family: "Liberation Sans", "DejaVu Sans", Arial, sans-serif;
  font-size: 14.5pt;
  line-height: 1.25;
  font-weight: 700;
  color: #0f172a;
  margin: 0 0 8pt 0;
  padding-bottom: 5pt;
  border-bottom: 1.5px solid #1e293b;
  break-after: avoid;
  page-break-after: avoid;
}
h2 {
  font-family: "Liberation Sans", "DejaVu Sans", Arial, sans-serif;
  font-size: 11.2pt;
  line-height: 1.28;
  font-weight: 700;
  color: #0f172a;
  margin: 12pt 0 5pt 0;
  padding-bottom: 2pt;
  border-bottom: 1px solid #cbd5e1;
  break-after: avoid;
  page-break-after: avoid;
}
h3 {
  font-family: "Liberation Sans", "DejaVu Sans", Arial, sans-serif;
  font-size: 9.8pt;
  line-height: 1.28;
  font-weight: 700;
  color: #1e293b;
  margin: 9pt 0 4pt 0;
  break-after: avoid;
  page-break-after: avoid;
}
h4 {
  font-family: "Liberation Sans", "DejaVu Sans", Arial, sans-serif;
  font-size: 9.0pt;
  line-height: 1.25;
  font-weight: 700;
  color: #334155;
  margin: 7pt 0 3pt 0;
  break-after: avoid;
  page-break-after: avoid;
}
p {
  margin: 0 0 6pt 0;
  text-align: justify;
  hyphens: auto;
}
ul, ol {
  margin: 0 0 6pt 0;
  padding-left: 16pt;
}
li {
  margin-bottom: 2.5pt;
}
hr {
  border: 0;
  border-top: 1px solid #cbd5e1;
  margin: 10pt 0;
}
a {
  color: #0369a1;
  text-decoration: none;
}
code {
  font-family: "Liberation Mono", "DejaVu Sans Mono", monospace;
  font-size: 7.8pt;
  background: #f1f5f9;
  padding: 0.5pt 2.5pt;
  border-radius: 2px;
  overflow-wrap: break-word;
  word-break: normal;
}
pre {
  background: #f8fafc;
  border: 1px solid #cbd5e1;
  border-left: 3px solid #0284c7;
  border-radius: 3px;
  padding: 6pt 8pt;
  margin: 5pt 0 8pt 0;
  font-family: "Liberation Mono", "DejaVu Sans Mono", monospace;
  font-size: 7.3pt;
  line-height: 1.28;
  white-space: pre-wrap;
  overflow-wrap: anywhere;
  word-break: break-word;
  break-inside: avoid;
  page-break-inside: avoid;
}
pre code {
  background: transparent;
  padding: 0;
  font-size: inherit;
}
table {
  width: 100%;
  border-collapse: collapse;
  font-family: "Liberation Sans", "DejaVu Sans", Arial, sans-serif;
  font-size: 6.5pt;
  line-height: 1.22;
  margin: 6pt 0 10pt 0;
  table-layout: auto;
}
thead {
  display: table-header-group;
}
tr {
  break-inside: avoid;
  page-break-inside: avoid;
}
th {
  background: #f1f5f9;
  color: #0f172a;
  font-weight: 700;
  border-top: 1.2px solid #334155;
  border-bottom: 1px solid #475569;
  padding: 3pt 3.5pt;
  text-align: left;
  vertical-align: bottom;
  overflow-wrap: break-word;
  word-break: normal;
}
td {
  border-bottom: 0.6px solid #cbd5e1;
  padding: 2.5pt 3.5pt;
  text-align: left;
  vertical-align: top;
  overflow-wrap: break-word;
  word-break: normal;
}
tbody tr:nth-child(even) {
  background: #f8fafc;
}
tbody tr:last-child td {
  border-bottom: 1.2px solid #334155;
}
table code {
  font-size: 6.2pt;
  padding: 0.2pt 1.5pt;
}
figure.paper-figure {
  margin: 8pt 0 10pt 0;
  padding: 6pt 7pt;
  border: 1px solid #cbd5e1;
  border-radius: 4px;
  background: #ffffff;
  text-align: center;
  break-inside: avoid;
  page-break-inside: avoid;
}
figure.paper-figure img {
  max-width: 100%;
  max-height: 165mm;
  width: auto;
  height: auto;
  display: block;
  margin: 0 auto 5pt auto;
}
figure.paper-figure figcaption {
  font-family: "Liberation Sans", "DejaVu Sans", Arial, sans-serif;
  font-size: 7.7pt;
  line-height: 1.32;
  color: #1e293b;
  text-align: left;
  border-top: 0.6px solid #e2e8f0;
  padding-top: 4pt;
}
.math-inline {
  font-family: "Liberation Serif", "Times New Roman", Georgia, serif;
  font-style: italic;
  white-space: nowrap;
}
.math-display {
  font-family: "Liberation Serif", "Times New Roman", Georgia, serif;
  font-style: italic;
  text-align: center;
  margin: 6pt 0;
  padding: 4pt 0;
  font-size: 9.6pt;
  break-inside: avoid;
  page-break-inside: avoid;
}
.math-rm {
  font-style: normal;
}
.math-op {
  font-style: normal;
  font-size: 110%;
}
blockquote {
  margin: 6pt 0 10pt 0;
  padding: 7pt 10pt;
  border: 1px solid #cbd5e1;
  border-left: 3.5px solid #0072b2;
  border-radius: 4px;
  background: #f8fafc;
  font-size: 8.6pt;
  line-height: 1.38;
  color: #1e293b;
  break-inside: avoid;
  page-break-inside: avoid;
}
blockquote p {
  margin: 0 0 4pt 0;
}
blockquote ul {
  margin: 2pt 0 0 0;
  padding-left: 14pt;
}
"""


def _postprocess_body_html(body_html: str) -> str:
    """Pair figure blocks with their caption paragraphs and add clean identifier break hints."""
    # 1. Pair <figure class="paper-figure">...</figure> with following <p><strong>Figure N...</strong>...</p>
    body_html = re.sub(
        r'<figure class="paper-figure">\s*(<img[^>]+/>)\s*</figure>\s*<p>(<strong>Figure\s+\d+[\s\S]*?)</p>',
        r'<figure class="paper-figure">\n  \1\n  <figcaption>\2</figcaption>\n</figure>',
        body_html,
    )

    # 2. Protect <pre>...</pre> blocks before inserting zero-width break opportunities
    pre_blocks: list[str] = []

    def _save_pre(m: re.Match[str]) -> str:
        idx = len(pre_blocks)
        pre_blocks.append(m.group(0))
        return f"@@PREBLOCK{idx}@@"

    s = re.sub(r"<pre>[\s\S]*?</pre>", _save_pre, body_html)

    # 3. Insert zero-width space (&#8203;) after '_' and '/' inside text nodes of <td> and inline <code>
    def _add_zwsp_in_tags(html_fragment: str) -> str:
        parts = re.split(r"(<[^>]+>)", html_fragment)
        for i in range(0, len(parts), 2):
            parts[i] = re.sub(r"([_/])([A-Za-z0-9])", r"\1&#8203;\2", parts[i])
        return "".join(parts)

    s = re.sub(r"<td>([\s\S]*?)</td>", lambda m: f"<td>{_add_zwsp_in_tags(m.group(1))}</td>", s)
    s = re.sub(r"<code>([\s\S]*?)</code>", lambda m: f"<code>{_add_zwsp_in_tags(m.group(1))}</code>", s)

    for idx, block in enumerate(pre_blocks):
        s = s.replace(f"@@PREBLOCK{idx}@@", block)

    return s


def markdown_to_html_document(md_text: str, base_dir: Path, title: str) -> str:
    """Convert compiled Markdown into a complete, self-contained HTML5 document for PDF printing."""
    prepped = preprocess_markdown_for_html(md_text, base_dir)
    body_html = markdown.markdown(
        prepped,
        extensions=["tables", "fenced_code", "sane_lists"],
    )
    body_html = _postprocess_body_html(body_html)
    return (
        "<!DOCTYPE html>\n"
        '<html lang="en">\n'
        "<head>\n"
        '  <meta charset="utf-8" />\n'
        f"  <title>{html.escape(title)}</title>\n"
        f"  <style>{PRINT_CSS}</style>\n"
        "</head>\n"
        "<body>\n"
        f"{body_html}\n"
        "</body>\n"
        "</html>\n"
    )


def find_chrome_binary() -> str | None:
    """Locate headless Chrome/Chromium binary if installed."""
    for candidate in ("google-chrome", "chromium", "chromium-browser", "/usr/bin/google-chrome"):
        found = shutil.which(candidate)
        if found:
            return found
    return None


def render_markdown_to_pdf(md_text: str, pdf_path: Path, base_dir: Path, title: str) -> None:
    """Render Markdown text to a typeset PDF file via headless Chrome and optimize via Ghostscript."""
    html_doc = markdown_to_html_document(md_text, base_dir, title=title)
    chrome_bin = find_chrome_binary()
    if chrome_bin is None:
        if pdf_path.is_file() and pdf_path.stat().st_size > 0:
            return
        raise RuntimeError("No Chrome/Chromium binary found for PDF compilation and no pre-rendered PDF exists.")
    pdf_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="vardict_paper_pdf_") as tmpdir:
        tmp_html = Path(tmpdir) / "document.html"
        tmp_pdf = Path(tmpdir) / "output.pdf"
        opt_pdf = Path(tmpdir) / "optimized.pdf"
        tmp_html.write_text(html_doc, encoding="utf-8")

        cmd = [
            chrome_bin,
            "--headless",
            "--disable-gpu",
            "--no-sandbox",
            "--disable-dev-shm-usage",
            "--no-pdf-header-footer",
            f"--print-to-pdf={tmp_pdf}",
            tmp_html.as_uri(),
        ]
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if not tmp_pdf.is_file() or tmp_pdf.stat().st_size == 0:
            raise RuntimeError(f"Chrome failed to produce PDF at {tmp_pdf}")

        gs_bin = shutil.which("gs")
        if gs_bin:
            gs_cmd = [
                gs_bin,
                "-sDEVICE=pdfwrite",
                "-dCompatibilityLevel=1.5",
                "-dNOPAUSE",
                "-dQUIET",
                "-dBATCH",
                f"-sOutputFile={opt_pdf}",
                str(tmp_pdf),
            ]
            res = subprocess.run(gs_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if res.returncode == 0 and opt_pdf.is_file() and opt_pdf.stat().st_size > 0:
                shutil.copyfile(opt_pdf, pdf_path)
                return

        shutil.copyfile(tmp_pdf, pdf_path)
