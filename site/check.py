#!/usr/bin/env python3
"""Check the generated static site against its acceptance criteria.

Usage:
    python3 site/check.py [site/scene.json]

Checks:
    - the scene JSON parses and has the default design (7 bodies, 6 joints,
      3 planets, gear ratio 3.5),
    - the viewer files exist and the page carries the honesty badge,
    - every markdown input has an HTML page under site/docs/.

This is a manual check. It is not part of ``just verify``. Standard library
only.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SITE_DIR = REPO_ROOT / "site"


def fail(message: str) -> "None":
    print("check: FAIL: " + message, file=sys.stderr)
    raise SystemExit(1)


def check_scene(path: Path) -> str:
    if not path.is_file():
        fail("missing scene file: %s" % path)
    data = json.loads(path.read_text(encoding="utf-8"))

    if data.get("schema") != "nanocad.scene":
        fail("scene schema is %r, expected 'nanocad.scene'" % data.get("schema"))

    device = data.get("device", {})
    design = data.get("design", {})
    atomistic = data.get("atomistic", {})
    coarse = data.get("coarse", {})

    checks = {
        "device.body_count": (device.get("body_count"), 7),
        "coarse.body_count": (coarse.get("body_count"), 7),
        "len(device.joints)": (len(device.get("joints", [])), 6),
        "design.planet_count": (design.get("planet_count"), 3),
    }
    for name, (actual, expected) in checks.items():
        if actual != expected:
            fail("%s is %r, expected %r" % (name, actual, expected))

    ratio = design.get("gear_ratio")
    if not isinstance(ratio, (int, float)) or abs(ratio - 3.5) > 1e-12:
        fail("gear_ratio is %r, expected 3.5" % ratio)

    atom_count = atomistic.get("atom_count")
    atoms = atomistic.get("atoms", [])
    if atom_count != len(atoms):
        fail("atom_count %r does not match len(atoms) %r" % (atom_count, len(atoms)))

    companion = path.with_name(path.stem + ".data.js")
    if not companion.is_file():
        fail("missing file:// companion: %s" % companion)
    if "window.NANOCAD_SCENE" not in companion.read_text(encoding="utf-8"):
        fail("%s does not assign window.NANOCAD_SCENE" % companion)

    return (
        "scene ok: bodies=%d joints=%d planets=%d ratio=%s atoms=%d"
        % (device["body_count"], len(device["joints"]), design["planet_count"], ratio, atom_count)
    )


def check_viewer() -> str:
    for name in ("index.html", "viewer.js", "style.css", "render_atoms.js"):
        if not (SITE_DIR / name).is_file():
            fail("missing viewer file: site/%s" % name)
    page = (SITE_DIR / "index.html").read_text(encoding="utf-8")
    if "simulated / schematic" not in page:
        fail("site/index.html does not carry the 'simulated / schematic' label")
    if "render_atoms.js" not in page:
        fail("site/index.html does not load render_atoms.js")
    return "viewer ok: index.html, viewer.js, style.css, render_atoms.js present"


def check_docs() -> str:
    docs = SITE_DIR / "docs"
    if not (docs / "index.html").is_file():
        fail("missing site/docs/index.html; run scripts/build_docs_site.py")

    inputs = sorted(REPO_ROOT.glob("*.md")) + sorted((REPO_ROOT / "docs").glob("*.md"))
    missing = [
        path.name
        for path in inputs
        if not (docs / (path.stem + ".html")).is_file()
    ]
    if missing:
        fail("missing HTML output for: " + ", ".join(missing))
    return "docs ok: %d markdown inputs have HTML pages" % len(inputs)


def main() -> "None":
    scene_path = Path(sys.argv[1]) if len(sys.argv) > 1 else SITE_DIR / "scene.json"
    print(check_scene(scene_path))
    print(check_viewer())
    print(check_docs())
    print("check: all checks passed")


if __name__ == "__main__":
    main()
