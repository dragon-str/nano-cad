#!/usr/bin/env python3
"""Build a small offline HTML documentation site from the repository markdown.

Usage:
    python3 scripts/build_docs_site.py

The script reads every top-level ``*.md`` file and every ``docs/*.md`` file.
It writes one ``site/docs/<stem>.html`` page for each input and a
``site/docs/index.html`` list. Standard library only. No network access.

The converter is deliberately simple. It handles headings, paragraphs, code
fences, lists, tables, blockquotes, horizontal rules, and the common inline
spans. It is honest about what it does: it escapes all text first and then
applies a small set of substitutions.
"""

from __future__ import annotations

import html
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "site" / "docs"

# These files must exist. A missing core document is a hard error.
REQUIRED_TOP_LEVEL = (
    "AGENTS.md",
    "ARCHITECTURE.md",
    "DECISIONS.md",
    "PARAMETERS.md",
    "PLAN.md",
    "README.md",
    "TASKS.md",
)

PAGE_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1" />
    <title>{title} - Nano-CAD docs</title>
    <link rel="stylesheet" href="../style.css" />
  </head>
  <body class="doc-page">
    <header>
      <a href="index.html">Nano-CAD docs</a>
      <span><a href="../index.html">Three-scale viewer</a></span>
    </header>
    <main>
{body}
    </main>
    <footer>
      Generated from <code>{source}</code> by
      <code>scripts/build_docs_site.py</code>. The project reports simulated
      results. The viewer and its scene are schematic, not a built device. No
      medical claim.
    </footer>
  </body>
</html>
"""


def fail(message: str) -> "None":
    """Print an error and exit with a non-zero status."""
    print("build_docs_site: error: " + message, file=sys.stderr)
    raise SystemExit(1)


# --------------------------------------------------------------------------
# Inline conversion
# --------------------------------------------------------------------------


def inline(text: str) -> str:
    """Convert the inline markdown spans of one text fragment to HTML."""
    placeholders: list[str] = []

    def stash_code(match: "re.Match[str]") -> str:
        placeholders.append(match.group(1))
        return "\x00%d\x00" % (len(placeholders) - 1)

    text = html.escape(text, quote=False)
    text = re.sub(r"`([^`]+)`", stash_code, text)

    def link(match: "re.Match[str]") -> str:
        label, href = match.group(1), match.group(2)
        if href.endswith(".md"):
            href = href[: -len(".md")] + ".html"
        return '<a href="%s">%s</a>' % (href, label)

    text = re.sub(r"\[([^\]]+)\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)", link, text)
    text = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", text)
    text = re.sub(r"__([^_]+)__", r"<strong>\1</strong>", text)
    text = re.sub(r"(?<!\*)\*([^*]+)\*(?!\*)", r"<em>\1</em>", text)
    text = re.sub(r"(?<!_)_([^_]+)_(?!_)", r"<em>\1</em>", text)

    for index, code in enumerate(placeholders):
        text = text.replace("\x00%d\x00" % index, "<code>%s</code>" % html.escape(code))
    return text


# --------------------------------------------------------------------------
# Block conversion
# --------------------------------------------------------------------------


def is_hr(line: str) -> bool:
    return bool(re.match(r"^\s*(-{3,}|\*{3,}|_{3,})\s*$", line))


def is_heading(line: str) -> bool:
    return bool(re.match(r"^#{1,6}\s+\S", line))


def is_ul(line: str) -> bool:
    return bool(re.match(r"^\s*[-*+]\s+\S", line))


def is_ol(line: str) -> bool:
    return bool(re.match(r"^\s*\d+\.\s+\S", line))


def is_quote(line: str) -> bool:
    return line.lstrip().startswith(">")


def is_table_start(lines: "list[str]", index: int) -> bool:
    if index + 1 >= len(lines):
        return False
    if not lines[index].lstrip().startswith("|"):
        return False
    return bool(re.match(r"^\s*\|?\s*:?-{2,}", lines[index + 1]))


def is_block_start(lines: "list[str]", index: int) -> bool:
    if index >= len(lines):
        return True
    line = lines[index]
    if not line.strip():
        return True
    return (
        line.lstrip().startswith("```")
        or is_heading(line)
        or is_hr(line)
        or is_ul(line)
        or is_ol(line)
        or is_quote(line)
        or is_table_start(lines, index)
    )


def split_row(line: str) -> "list[str]":
    line = line.strip()
    if line.startswith("|"):
        line = line[1:]
    if line.endswith("|"):
        line = line[:-1]
    return [cell.strip() for cell in line.split("|")]


def parse_table(lines: "list[str]", index: int) -> "tuple[str, int]":
    header = split_row(lines[index])
    index += 2
    rows: list[list[str]] = []
    while index < len(lines) and lines[index].lstrip().startswith("|"):
        rows.append(split_row(lines[index]))
        index += 1
    parts = ["<table>", "<thead><tr>"]
    for cell in header:
        parts.append("<th>%s</th>" % inline(cell))
    parts.append("</tr></thead><tbody>")
    for row in rows:
        parts.append("<tr>")
        for cell in row:
            parts.append("<td>%s</td>" % inline(cell))
        parts.append("</tr>")
    parts.append("</tbody></table>")
    return "".join(parts), index


def parse_list(lines: "list[str]", index: int, ordered: bool) -> "tuple[str, int]":
    tag = "ol" if ordered else "ul"
    pattern = is_ol if ordered else is_ul
    parts = ["<%s>" % tag]
    while index < len(lines) and pattern(lines[index]):
        content = re.sub(r"^\s*(?:[-*+]|\d+\.)\s+", "", lines[index])
        index += 1
        continuation: list[str] = []
        while (
            index < len(lines)
            and lines[index].strip()
            and not is_block_start(lines, index)
            and (lines[index].startswith(" ") or lines[index].startswith("\t"))
        ):
            continuation.append(lines[index].strip())
            index += 1
        if continuation:
            content = content + " " + " ".join(continuation)
        parts.append("<li>%s</li>" % inline(content))
    parts.append("</%s>" % tag)
    return "".join(parts), index


def convert_markdown(text: str) -> str:
    lines = text.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    out: list[str] = []
    index = 0
    count = len(lines)

    while index < count:
        line = lines[index]

        if line.lstrip().startswith("```"):
            language = line.lstrip()[3:].strip()
            index += 1
            code: list[str] = []
            while index < count and not lines[index].lstrip().startswith("```"):
                code.append(lines[index])
                index += 1
            index += 1
            attribute = ' class="language-%s"' % html.escape(language) if language else ""
            out.append(
                "<pre><code%s>%s</code></pre>"
                % (attribute, html.escape("\n".join(code)))
            )
            continue

        if is_hr(line):
            out.append("<hr>")
            index += 1
            continue

        match = re.match(r"^(#{1,6})\s+(.*?)\s*#*\s*$", line)
        if match:
            level = len(match.group(1))
            out.append("<h%d>%s</h%d>" % (level, inline(match.group(2)), level))
            index += 1
            continue

        if is_table_start(lines, index):
            table, index = parse_table(lines, index)
            out.append(table)
            continue

        if is_quote(line):
            quoted: list[str] = []
            while index < count and is_quote(lines[index]):
                quoted.append(lines[index].lstrip()[1:].lstrip())
                index += 1
            out.append(
                "<blockquote>%s</blockquote>" % convert_markdown("\n".join(quoted))
            )
            continue

        if is_ul(line):
            items, index = parse_list(lines, index, ordered=False)
            out.append(items)
            continue

        if is_ol(line):
            items, index = parse_list(lines, index, ordered=True)
            out.append(items)
            continue

        if not line.strip():
            index += 1
            continue

        paragraph: list[str] = []
        while index < count and not is_block_start(lines, index):
            paragraph.append(lines[index].strip())
            index += 1
        if paragraph:
            out.append("<p>%s</p>" % inline(" ".join(paragraph)))

    return "\n".join(out)


def first_heading(text: str) -> "str | None":
    for line in text.splitlines():
        match = re.match(r"^#\s+(.*?)\s*#*\s*$", line)
        if match:
            return match.group(1)
    return None


# --------------------------------------------------------------------------
# Pages
# --------------------------------------------------------------------------


def render_page(title: str, body: str, source: str) -> str:
    return PAGE_TEMPLATE.format(
        title=html.escape(title),
        body=body,
        source=html.escape(source),
    )


def render_index(pages: "list[dict[str, str]]") -> str:
    repository = [page for page in pages if page["group"] == "root"]
    docs = [page for page in pages if page["group"] == "docs"]

    parts = [
        "<h1>Nano-CAD documentation</h1>",
        "<p>Generated from the repository markdown by "
        "<code>scripts/build_docs_site.py</code>. Every page is a static file. "
        "The project reports simulated results, not a built device.</p>",
    ]
    for heading, group in (("Repository", repository), ("docs/", docs)):
        parts.append('<div class="doc-index-group">')
        parts.append("<h2>%s</h2>" % heading)
        parts.append("<ul>")
        for page in group:
            parts.append(
                '<li><a href="%s.html">%s</a></li>'
                % (html.escape(page["stem"]), html.escape(page["title"]))
            )
        parts.append("</ul></div>")

    body = "\n".join(parts)
    return PAGE_TEMPLATE.format(title="Index", body=body, source="generated index")


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------


def collect_inputs() -> "list[dict[str, str]]":
    docs_dir = REPO_ROOT / "docs"
    if not docs_dir.is_dir():
        fail("missing input directory: docs/")

    missing = [
        name for name in REQUIRED_TOP_LEVEL if not (REPO_ROOT / name).is_file()
    ]
    if missing:
        fail("missing input file(s): " + ", ".join(missing))

    inputs: list[dict[str, str]] = []
    for path in sorted(REPO_ROOT.glob("*.md")):
        inputs.append({"path": str(path), "group": "root", "stem": path.stem})
    for path in sorted(docs_dir.glob("*.md")):
        inputs.append({"path": str(path), "group": "docs", "stem": path.stem})

    if not inputs:
        fail("no markdown inputs found")

    seen: dict[str, str] = {}
    for item in inputs:
        if item["stem"] in seen:
            fail(
                "output name collision for %r: %s and %s"
                % (item["stem"], seen[item["stem"]], item["path"])
            )
        seen[item["stem"]] = item["path"]
    return inputs


def main() -> "None":
    inputs = collect_inputs()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    pages: list[dict[str, str]] = []
    for item in inputs:
        source_path = Path(item["path"])
        text = source_path.read_text(encoding="utf-8")
        title = first_heading(text) or item["stem"]
        body = convert_markdown(text)
        relative = source_path.relative_to(REPO_ROOT)
        out_path = OUT_DIR / (item["stem"] + ".html")
        out_path.write_text(render_page(title, body, str(relative)), encoding="utf-8")
        pages.append(
            {
                "stem": item["stem"],
                "group": item["group"],
                "title": title,
                "path": str(relative),
            }
        )

    index_path = OUT_DIR / "index.html"
    index_path.write_text(render_index(pages), encoding="utf-8")

    for item in inputs:
        out_path = OUT_DIR / (item["stem"] + ".html")
        if not out_path.is_file() or out_path.stat().st_size == 0:
            fail("failed to write output: %s" % out_path)
    if not index_path.is_file() or index_path.stat().st_size == 0:
        fail("failed to write output: %s" % index_path)

    print(
        "build_docs_site: wrote %d page(s) plus index.html to %s"
        % (len(inputs), OUT_DIR)
    )
    for page in sorted(pages, key=lambda entry: entry["stem"]):
        print("  %s -> %s.html" % (page["path"], page["stem"]))


if __name__ == "__main__":
    main()
