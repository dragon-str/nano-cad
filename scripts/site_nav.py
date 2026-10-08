#!/usr/bin/env python3
"""The shared header of every Nano-CAD page, from one list (ADR-0135).

Each app page carries the header between two markers::

    <!-- site-header -->
    ...
    <!-- /site-header -->

`python3 scripts/site_nav.py` rewrites that block in every page.
`python3 scripts/site_nav.py --check` fails if any page differs from what this
file renders; `site/check.py` runs the check. `scripts/build_docs_site.py`
imports `render_header` for the generated docs pages.
"""

from __future__ import annotations

import html
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SITE = REPO_ROOT / "site"

START = "<!-- site-header -->"
END = "<!-- /site-header -->"

# (key, href relative to site/, label). Ordered from the whole device down to
# the atoms, then the documents.
NAV: tuple[tuple[str, str, str], ...] = (
    ("overview", "overview.html", "Overview"),
    ("device", "device.html", "Device explorer"),
    ("system", "system.html", "System model"),
    ("models", "models.html", "3D models"),
    ("viewer", "index.html", "Atomistic viewer"),
    ("docs", "docs/index.html", "Docs"),
)

# file -> (nav key, page title shown in the header, badge, id of the title span)
PAGES: dict[str, tuple[str, str, str, str | None]] = {
    "overview.html": ("overview", "What is here and how the views connect", "simulated", None),
    "device.html": ("device", "Device scale explorer", "simulated / schematic", None),
    "system.html": ("system", "Respirocyte system model (real size, 1 um hull)", "simulated / estimated", None),
    "models.html": ("models", "3D models at true relative scale", "simulated / schematic", None),
    "index.html": ("viewer", "Loading ...", "simulated / schematic", "scene-title"),
}


def render_header(current: str, title: str, badge: str, prefix: str = "", title_id: str | None = None,
                  indent: str = "    ") -> str:
    """Return the header block, markers included, for the page `current`."""
    links = []
    for key, href, label in NAV:
        mark = ' aria-current="page"' if key == current else ""
        links.append(f'{indent}    <a href="{prefix}{href}"{mark}>{html.escape(label)}</a>')
    tid = f' id="{title_id}"' if title_id else ""
    lines = [
        START,
        f'{indent}<header class="site-header">',
        f'{indent}  <a class="site-brand" href="{prefix}overview.html">Nano-CAD</a>',
        f'{indent}  <span class="site-page"{tid}>{html.escape(title)}</span>',
        f'{indent}  <span class="site-badge">{html.escape(badge)}</span>',
        f'{indent}  <nav class="site-nav" aria-label="Nano-CAD views">',
        *links,
        f"{indent}  </nav>",
        f"{indent}</header>",
        f"{indent}{END}",
    ]
    return "\n".join(lines)


def _expected(name: str, text: str) -> tuple[int, int, str]:
    start = text.find(START)
    end = text.find(END)
    if start < 0 or end < start:
        raise ValueError(f"{name}: no {START} ... {END} block")
    line_start = text.rfind("\n", 0, start) + 1
    indent = text[line_start:start]
    key, title, badge, title_id = PAGES[name]
    block = render_header(key, title, badge, title_id=title_id, indent=indent)
    return start, end + len(END), block


def check() -> list[str]:
    """Return one message for each page whose header is missing or stale."""
    problems = []
    for name in PAGES:
        path = SITE / name
        if not path.exists():
            problems.append(f"{name}: missing")
            continue
        text = path.read_text(encoding="utf-8")
        try:
            a, b, block = _expected(name, text)
        except ValueError as err:
            problems.append(str(err))
            continue
        if text[a:b] != block:
            problems.append(f"{name}: header differs from scripts/site_nav.py (run it to rewrite)")
        if 'href="nav.css"' not in text:
            problems.append(f"{name}: does not link nav.css")
    return problems


def write() -> None:
    for name in PAGES:
        path = SITE / name
        text = path.read_text(encoding="utf-8")
        a, b, block = _expected(name, text)
        if text[a:b] != block:
            path.write_text(text[:a] + block + text[b:], encoding="utf-8")
            print(f"site_nav: rewrote {name}")


def main(argv: list[str]) -> int:
    if "--check" in argv:
        problems = check()
        for p in problems:
            print("site_nav: " + p, file=sys.stderr)
        return 1 if problems else 0
    write()
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
