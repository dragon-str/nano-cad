#!/usr/bin/env python3
"""Render the nano-cad 90-second gearbox video from scratch.

This script draws every frame with Pillow, assembles the machine narration with
/usr/bin/say, and calls ffmpeg to write the mp4 and the srt. It invents no
number: every figure on screen comes from the files named in the script's
source column.

The video is schematic. It is not a photograph of a physical device. No
physical device exists.

Usage (normally through scripts/make_video.sh):
    python3 scripts/render_video.py \
        --out docs/media/nano-cad-gearbox.mp4 \
        --srt  docs/media/nano-cad-gearbox.srt \
        --build /tmp/ncad-video-build \
        --ffmpeg /opt/homebrew/bin/ffmpeg \
        --ffprobe /opt/homebrew/bin/ffprobe \
        --say /usr/bin/say
"""

from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
import sys
import textwrap

from PIL import Image, ImageDraw, ImageFont

W, H, FPS = 1920, 1080, 30
REPO_URL = "https://github.com/dragon-str/nano-cad"

# ---------------------------------------------------------------- palette
BG = (11, 15, 24)
PANEL = (20, 26, 39)
PANEL2 = (26, 34, 50)
LINE = (48, 60, 84)
FG = (232, 238, 248)
DIM = (128, 140, 163)
ACCENT = (86, 184, 255)
GREEN = (94, 220, 150)
AMBER = (255, 191, 82)
RED = (255, 116, 116)
CYAN = (120, 226, 232)
VIOLET = (178, 152, 255)

# ---------------------------------------------------------------- shots
# Start and end times follow docs/video-script.md exactly. Narration for shots
# 9 and 10 is shortened so that it fits the two short windows.
SHOTS = [
    {
        "n": 1, "start": 0.0, "end": 8.0,
        "title": "Cold open: the simulated gearbox",
        "scene": "gears",
        "narration": (
            "This is a planetary gearbox. It is simulated at the atomic scale "
            "and at the device scale. The sun gear turns. The planets turn. "
            "The carrier turns."
        ),
        "source": "TASKS.md M6-06; docs/reproduce.md",
    },
    {
        "n": 2, "start": 8.0, "end": 19.0,
        "title": "What this is: the five model levels",
        "scene": "levels",
        "narration": (
            "nano-cad is a multiscale design tool. It has five model levels. "
            "Only L1, the atomistic engine, and L2, the device layer, are "
            "implemented. L0, L3, and L4 are later work."
        ),
        "source": "ARCHITECTURE.md; docs/glossary.md",
    },
    {
        "n": 3, "start": 19.0, "end": 29.0,
        "title": "Clean checkout and one command",
        "scene": "reproduce",
        "narration": (
            "A clean checkout reproduces the whole demo with one command. "
            "The first run needs network access. It downloads the build and "
            "test tools."
        ),
        "source": "docs/reproduce.md; README.md",
    },
    {
        "n": 4, "start": 29.0, "end": 40.0,
        "title": "The Rust gates pass",
        "scene": "tests",
        "narration": (
            "Every force term has a finite-difference gradient test. Every "
            "integrator has an energy-conservation test. The tests are a "
            "merge gate, not a phase."
        ),
        "source": "TASKS.md M2-02, M2-03, M2-08, M2-10, M2-11",
    },
    {
        "n": 5, "start": 40.0, "end": 52.0,
        "title": "The agent generates the gear set",
        "scene": "generate",
        "narration": (
            "The agent calls typed tools. It never edits coordinates. It asks "
            "for a planetary set with three planets. The generator returns "
            "the atoms and the bonds."
        ),
        "source": "TASKS.md M4-06; docs/three-scale.md; mcp/README.md",
    },
    {
        "n": 6, "start": 52.0, "end": 62.0,
        "title": "The agent assembles the L2 device",
        "scene": "assemble",
        "narration": (
            "The agent assembles the L2 device. The device has seven bodies "
            "and six joints. The constraints are position-level projections."
        ),
        "source": "TASKS.md M7-04, M8-01; mcp/tests/test_agent_demo.py",
    },
    {
        "n": 7, "start": 62.0, "end": 72.0,
        "title": "The agent measures the gear ratio",
        "scene": "measure",
        "narration": (
            "The agent drives the sun and measures the sun-to-carrier ratio. "
            "The analytic ratio is three point five. The simulated ratio is "
            "three point five. The relative error is one point one four six "
            "times ten to the minus ten."
        ),
        "source": "TASKS.md M6-06; docs/reproduce.md",
    },
    {
        "n": 8, "start": 72.0, "end": 82.0,
        "title": "URDF and the three-scale scene",
        "scene": "urdf",
        "narration": (
            "The agent exports a URDF with seven links and six joints. The "
            "scene export carries three layers: atomistic, device, and "
            "coarse. A viewer animates across them."
        ),
        "source": "TASKS.md M8-01; docs/three-scale.md",
    },
    {
        "n": 9, "start": 82.0, "end": 88.0,
        "title": "Benchmark timings",
        "scene": "bench",
        "narration": (
            "Six hundred atoms. Serial: six hundred eighty microseconds. "
            "Four workers: three hundred seventy-one."
        ),
        "source": "benchmarks/results/2026-09-14-engine-timings.txt",
    },
    {
        "n": 10, "start": 88.0, "end": 90.0,
        "title": "Caveats and closing card",
        "scene": "closing",
        "narration": "All results are simulated.",
        "source": "PLAN.md; README.md; TASKS.md M2-11",
    },
]

TOTAL = SHOTS[-1]["end"]

# ---------------------------------------------------------------- fonts
FONT_DIR = "/System/Library/Fonts"
FONT_SUP = "/System/Library/Fonts/Supplemental"


def _load(path, size, index=0):
    try:
        return ImageFont.truetype(path, size, index=index)
    except Exception:
        return ImageFont.load_default()


FONTS = {}


def fonts():
    if FONTS:
        return FONTS
    FONTS["title"] = _load(f"{FONT_SUP}/Arial Bold.ttf", 60)
    FONTS["h1"] = _load(f"{FONT_SUP}/Arial Bold.ttf", 40)
    FONTS["h2"] = _load(f"{FONT_SUP}/Arial Bold.ttf", 28)
    FONTS["body"] = _load(f"{FONT_SUP}/Arial.ttf", 30)
    FONTS["body_b"] = _load(f"{FONT_SUP}/Arial Bold.ttf", 30)
    FONTS["small"] = _load(f"{FONT_SUP}/Arial.ttf", 22)
    FONTS["small_b"] = _load(f"{FONT_SUP}/Arial Bold.ttf", 22)
    FONTS["mono"] = _load(f"{FONT_DIR}/Menlo.ttc", 28, 0)
    FONTS["mono_b"] = _load(f"{FONT_DIR}/Menlo.ttc", 28, 1)
    FONTS["mono_s"] = _load(f"{FONT_DIR}/Menlo.ttc", 22, 0)
    FONTS["mono_l"] = _load(f"{FONT_DIR}/Menlo.ttc", 34, 0)
    FONTS["num"] = _load(f"{FONT_DIR}/Menlo.ttc", 62, 1)
    FONTS["caption"] = _load(f"{FONT_SUP}/Arial.ttf", 31)
    return FONTS


# ---------------------------------------------------------------- helpers
def txt(d, *a, **kw):
    """Accept txt(d, x, y, s, font, fill[, anchor]) or txt(d, (x, y), s, ...)."""
    if isinstance(a[0], (tuple, list)):
        xy, s, font, fill = a[0], a[1], a[2], a[3]
        anchor = a[4] if len(a) > 4 else kw.get("anchor", "la")
    else:
        xy, s, font, fill = (a[0], a[1]), a[2], a[3], a[4]
        anchor = a[5] if len(a) > 5 else kw.get("anchor", "la")
    d.text(tuple(xy), s, font=font, fill=fill, anchor=anchor)


def wrap(s, font, maxw, d):
    words = s.split()
    lines, cur = [], ""
    for w in words:
        trial = w if not cur else cur + " " + w
        if d.textlength(trial, font=font) <= maxw:
            cur = trial
        else:
            if cur:
                lines.append(cur)
            cur = w
    if cur:
        lines.append(cur)
    return lines


def panel(d, x0, y0, x1, y1, fill=PANEL, outline=LINE, r=16, width=2):
    d.rounded_rectangle([x0, y0, x1, y1], radius=r, fill=fill, outline=outline,
                        width=width)


def gear_points(cx, cy, r_tip, r_root, teeth, angle):
    pts = []
    step = 2.0 * math.pi / teeth
    for k in range(teeth):
        c = angle + k * step
        for frac, rad in ((0.06, r_root), (0.24, r_tip),
                          (0.76, r_tip), (0.94, r_root)):
            a = c + step * frac
            pts.append((cx + rad * math.cos(a), cy + rad * math.sin(a)))
    return pts


def draw_gear(d, cx, cy, r_tip, r_root, teeth, angle, fill, outline,
              bore=0.0, hub=0.0):
    d.polygon(gear_points(cx, cy, r_tip, r_root, teeth, angle),
              fill=fill, outline=outline)
    if bore > 0:
        d.ellipse([cx - bore, cy - bore, cx + bore, cy + bore], fill=BG)
    if hub > 0:
        d.ellipse([cx - hub, cy - hub, cx + hub, cy + hub],
                  fill=fill, outline=outline, width=2)


def arrow(d, x, y, dx, dy, fill, width=3, head=12):
    ex, ey = x + dx, y + dy
    d.line([x, y, ex, ey], fill=fill, width=width)
    ang = math.atan2(dy, dx)
    for off in (2.6, -2.6):
        hx = ex + head * math.cos(ang + off)
        hy = ey + head * math.sin(ang + off)
        d.line([ex, ey, hx, hy], fill=fill, width=width)


# ---------------------------------------------------------------- scene 1
def scene_gears(r, t, shot):
    d = r["d"]
    dur = shot["end"] - shot["start"]
    cx, cy = W * 0.5, H * 0.5 - 10
    scale = 300.0 / 1.5e-8
    r_sun, r_pl, r_ring, r_car = 6e-9 * scale, 4.5e-9 * scale, 1.5e-8 * scale, 1.05e-8 * scale
    a = t * 0.34

    d.ellipse([cx - r_car - 12, cy - r_car - 12, cx + r_car + 12, cy + r_car + 12],
              outline=DIM, width=2)
    txt(d, (cx, cy + r_car + 26), "carrier", fonts()["small"], DIM, "ma")

    for k in range(3):
        ang = a + 2 * math.pi * k / 3
        px, py = cx + r_car * math.cos(ang), cy + r_car * math.sin(ang)
        draw_gear(d, px, py, r_pl, r_pl * 0.80, 18, -ang * 3.5 + a,
                  (44, 70, 104), CYAN, bore=r_pl * 0.22, hub=r_pl * 0.30)
    draw_gear(d, cx, cy, r_sun, r_sun * 0.80, 24, -a * 3.5, (58, 84, 122),
              ACCENT, bore=r_sun * 0.22, hub=r_sun * 0.32)

    d.ellipse([cx - r_ring, cy - r_ring, cx + r_ring, cy + r_ring],
              outline=AMBER, width=20)
    d.ellipse([cx - r_ring + 20, cy - r_ring + 20, cx + r_ring - 20, cy + r_ring - 20],
              outline=AMBER, width=3)
    for k in range(60):
        ang = k * 2 * math.pi / 60
        x0 = cx + r_ring * math.cos(ang)
        y0 = cy + r_ring * math.sin(ang)
        x1 = cx + (r_ring - 24) * math.cos(ang)
        y1 = cy + (r_ring - 24) * math.sin(ang)
        d.line([x0, y0, x1, y1], fill=AMBER, width=3)

    fade = max(0.0, 1.0 - t / 3.0)
    if fade > 0.02:
        for k in range(60):
            ang = k * 2.399963
            rad = r_sun * math.sqrt((k + 0.5) / 60)
            col = tuple(int(BG[i] + (CYAN[i] - BG[i]) * fade) for i in range(3))
            x, y = cx + rad * math.cos(ang), cy + rad * math.sin(ang)
            d.ellipse([x - 4, y - 4, x + 4, y + 4], fill=col)
        txt(d, (cx, cy - r_ring - 46), "atomistic layer", fonts()["small"],
            tuple(int(BG[i] + (CYAN[i] - BG[i]) * fade) for i in range(3)), "ma")

    panel(d, 120, H - 210, 660, H - 120, PANEL, LINE)
    txt(d, 150, H - 196, "sun ratio", fonts()["small"], DIM)
    txt(d, 150, H - 170, "3.500000000", fonts()["mono_b"], GREEN)
    txt(d, 150, H - 132, "simulated  " + shot["source"], fonts()["small"], DIM)


# ---------------------------------------------------------------- scene 2
def scene_levels(r, t, shot):
    d = r["d"]
    levels = [
        ("L0", "quantum", "later", DIM),
        ("L1", "atomistic engine", "implemented", GREEN),
        ("L2", "device layer", "implemented", GREEN),
        ("L3", "continuum", "later", DIM),
        ("L4", "system", "later", DIM),
    ]
    x0, y0, bw, bh, gap = 220, 210, 1480, 118, 22
    lit = t / (shot["end"] - shot["start"])
    for i, (code, name, state, col) in enumerate(levels):
        y = y0 + i * (bh + gap)
        active = col is GREEN
        local = max(0.0, min(1.0, lit * 5 - i))
        c = tuple(int(DIM[j] + (col[j] - DIM[j]) * local) for j in range(3))
        fill = PANEL2 if active else PANEL
        panel(d, x0, y, x0 + bw, y + bh, fill, c, r=14, width=3)
        txt(d, x0 + 40, y + bh // 2, code, fonts()["title"], c, "lm")
        txt(d, x0 + 170, y + bh // 2 - 16, name, fonts()["h1"], FG, "lm")
        tag = state.upper()
        txt(d, x0 + bw - 40, y + bh // 2, tag, fonts()["small_b"],
            c if active else DIM, "rm")
    txt(d, x0, y0 - 46, "five model levels", fonts()["h2"], FG)
    txt(d, x0 + bw, y0 - 46, "simulated", fonts()["h2"], AMBER, "ra")


# ---------------------------------------------------------------- terminal
def term_window(r, x0, y0, x1, y1, title):
    d = r["d"]
    panel(d, x0, y0, x1, y1, (8, 11, 18), LINE, r=14, width=2)
    d.rounded_rectangle([x0, y0, x1, y0 + 44], radius=14, fill=PANEL2)
    d.rectangle([x0, y0 + 30, x1, y0 + 44], fill=PANEL2)
    for i, c in enumerate((RED, AMBER, GREEN)):
        d.ellipse([x0 + 22 + i * 30, y0 + 15, x0 + 38 + i * 30, y0 + 31], fill=c)
    txt(d, (x0 + (x1 - x0) / 2, y0 + 22), title, fonts()["small"], DIM, "mm")


def typed(s, t, t0, cps=28.0):
    if t < t0:
        return ""
    return s[: int((t - t0) * cps)]


def scene_reproduce(r, t, shot):
    d = r["d"]
    x0, y0, x1, y1 = 180, 190, 1180, 900
    term_window(r, x0, y0, x1, y1, "sh — nano-cad")
    lines = []
    lines.append(("$ " + typed(f"git clone {REPO_URL}", t, 0.3, 42), ACCENT,
                  t >= 0.3))
    lines.append(("$ " + typed("cd nano-cad", t, 1.9, 34), ACCENT, t >= 1.9))
    lines.append(("$ " + typed("just reproduce", t, 2.9, 30), ACCENT, t >= 2.9))
    steps = [
        "==> Step 1/7: check toolchain",
        "==> Step 2/7: create virtual environment",
        "==> Step 3/7: install maturin and pytest",
        "==> Step 4/7: build the extension",
        "==> Step 5/7: cargo test --workspace",
        "==> Step 6/7: pytest mcp/tests/test_agent_demo.py",
        "==> Step 7/7: key measured result",
    ]
    for i, s in enumerate(steps):
        st = 4.0 + i * 0.78
        lines.append((s, GREEN if t >= st else DIM, t >= st))
    yy = y0 + 78
    for s, c, on in lines:
        if on:
            txt(d, x0 + 30, yy, s, fonts()["mono_s"], c)
        yy += 40

    panel(d, 1240, 210, 1780, 420, PANEL, LINE)
    txt(d, 1272, 236, "7 steps", fonts()["num"], ACCENT)
    txt(d, 1272, 316, "the demo runs end to end", fonts()["body"], FG)
    txt(d, 1272, 356, "source: docs/reproduce.md", fonts()["small"], DIM)
    panel(d, 1240, 450, 1780, 610, PANEL, AMBER, width=3)
    txt(d, 1272, 476, "needs network", fonts()["h1"], AMBER)
    txt(d, 1272, 522, "on the first run", fonts()["h1"], AMBER)
    txt(d, 1272, 566, "downloads maturin and pytest", fonts()["small"], DIM)


# ---------------------------------------------------------------- scene 4
def scene_tests(r, t, shot):
    d = r["d"]
    x0, y0, x1, y1 = 150, 180, 1210, 900
    term_window(r, x0, y0, x1, y1, "cargo test --workspace")
    rows = [
        "running 440 tests",
        "test bond_stretch::the_analytic_gradient_matches_a_central_finite_difference ... ok",
        "test angle_bend::the_analytic_gradient_matches_a_central_finite_difference ... ok",
        "test torsion::the_analytic_gradient_matches_a_central_finite_difference ... ok",
        "test out_of_plane::the_analytic_gradient_matches_a_central_finite_difference ... ok",
        "test van_der_waals::the_analytic_gradient_matches_a_central_finite_difference ... ok",
        "test electrostatic::the_analytic_gradient_matches_a_central_finite_difference ... ok",
        "test system::the_nonbonded_pair_list_path_matches_the_all_pairs_reference ... ok",
        "test integrator::the_nve_total_energy_drift_stays_below_the_bound ... ok",
        "test result: ok. 440 passed; 0 failed",
    ]
    n = int(min(len(rows), max(0.0, (t - 0.6) * 1.15)))
    yy = y0 + 76
    for i, s in enumerate(rows[:n]):
        c = GREEN if ("finite_difference" in s or "ok." in s) else DIM
        if "result:" in s:
            c = GREEN
        txt(d, x0 + 24, yy, s[:96], fonts()["mono_s"], c)
        yy += 38
    txt(d, x0 + 24, y1 - 60,
        "440 Rust tests + 75 Python and MCP tests pass  (README.md)",
        fonts()["small"], DIM)

    panel(d, 1250, 180, 1800, 520, PANEL, LINE)
    txt(d, 1282, 206, "finite-difference gradient error", fonts()["small_b"], CYAN)
    txt(d, 1282, 240, "simulated", fonts()["small"], AMBER)
    rows2 = [
        ("bond stretch", "1.28e-17 N"),
        ("angle bend", "8.3e-18 N"),
        ("total (6 terms)", "5.7e-17 N"),
    ]
    yy = 286
    for name, val in rows2:
        txt(d, 1282, yy, name, fonts()["body"], FG)
        txt(d, 1790, yy, val, fonts()["mono_b"], GREEN, "ra")
        yy += 62
    txt(d, 1282, 466, "source: TASKS.md M2-02, M2-03, M2-08", fonts()["small"], DIM)

    panel(d, 1250, 545, 1800, 760, PANEL, LINE)
    txt(d, 1282, 571, "integrator", fonts()["small_b"], CYAN)
    txt(d, 1282, 606, "NVE drift", fonts()["body"], FG)
    txt(d, 1790, 606, "1.278e-5", fonts()["mono_b"], GREEN, "ra")
    txt(d, 1282, 646, "bound 1.0e-4", fonts()["small"], DIM)
    txt(d, 1282, 682, "source: TASKS.md M2-10", fonts()["small"], DIM)

    panel(d, 1250, 785, 1800, 900, PANEL, AMBER, width=3)
    txt(d, 1282, 806, "OpenMM cross-check: a coding check,", fonts()["small"], AMBER)
    txt(d, 1282, 836, "not a physics validation.", fonts()["small"], AMBER)
    txt(d, 1282, 866, "source: openmm-water-crosscheck.txt", fonts()["small"], DIM)


# ---------------------------------------------------------------- gear art
def planetary_art(r, cx, cy, scale, t, show_joints=False, show_ratio=None):
    d = r["d"]
    r_sun, r_pl, r_ring, r_car = 6e-9 * scale, 4.5e-9 * scale, 1.5e-8 * scale, 1.05e-8 * scale
    a = t * 0.34
    d.ellipse([cx - r_car - 10, cy - r_car - 10, cx + r_car + 10, cy + r_car + 10],
              outline=DIM, width=2)
    for k in range(3):
        ang = a + 2 * math.pi * k / 3
        px, py = cx + r_car * math.cos(ang), cy + r_car * math.sin(ang)
        draw_gear(d, px, py, r_pl, r_pl * 0.8, 18, -ang * 3.5, (44, 70, 104),
                  CYAN, bore=r_pl * 0.22, hub=r_pl * 0.30)
        if r_pl > 40:
            txt(d, (px, py + r_pl * 0.55), f"planet {k+1}", fonts()["small"], DIM, "ma")
    draw_gear(d, cx, cy, r_sun, r_sun * 0.8, 24, -a * 3.5, (58, 84, 122),
              ACCENT, bore=r_sun * 0.22, hub=r_sun * 0.32)
    if r_sun > 40:
        txt(d, (cx, cy + r_sun * 0.58), "sun 24", fonts()["small"], DIM, "ma")
    d.ellipse([cx - r_ring, cy - r_ring, cx + r_ring, cy + r_ring],
              outline=AMBER, width=14)
    for k in range(60):
        ang = k * 2 * math.pi / 60
        d.line([cx + r_ring * math.cos(ang), cy + r_ring * math.sin(ang),
                cx + (r_ring - 18) * math.cos(ang), cy + (r_ring - 18) * math.sin(ang)],
               fill=AMBER, width=3)
    if show_joints:
        for k in range(3):
            ang = a + 2 * math.pi * k / 3
            px, py = cx + r_car * math.cos(ang), cy + r_car * math.sin(ang)
            d.ellipse([px - 9, py - 9, px + 9, py + 9], outline=VIOLET, width=3)
        d.ellipse([cx - 9, cy - 9, cx + 9, cy + 9], outline=VIOLET, width=3)
        arrow(d, cx + r_ring + 20, cy, 70, 0, RED)
        arrow(d, cx + r_ring + 20, cy, 0, -70, RED)
        txt(d, (cx + r_ring + 30, cy - 96), "z axis", fonts()["small"], RED)
    if show_ratio is not None:
        txt(d, (cx, cy - r_ring - 40), show_ratio, fonts()["mono_l"], GREEN, "ma")


def json_box(r, x0, y0, x1, y1, header, body_lines, t, appear=0.0):
    d = r["d"]
    panel(d, x0, y0, x1, y1, (8, 11, 18), LINE)
    txt(d, x0 + 24, y0 + 20, header, fonts()["small_b"], ACCENT)
    yy = y0 + 66
    shown = int(max(0.0, (t - appear)) * 9)
    for i, s in enumerate(body_lines):
        if i > shown:
            break
        txt(d, x0 + 24, yy, s, fonts()["mono_s"], CYAN)
        yy += 30


# ---------------------------------------------------------------- scene 5
def scene_generate(r, t, shot):
    d = r["d"]
    body = [
        '{"jsonrpc":"2.0","id":2,',
        ' "method":"tools/call",',
        ' "params":{"name":"generate_gear",',
        '  "arguments":{"name":"planetary",',
        '   "specs":{"planet_count":3}}}}',
    ]
    json_box(r, 130, 210, 900, 560, "agent -> MCP  (JSON-RPC 2.0)", body, t, 0.2)
    txt(d, 154, 596, "handle: part-1", fonts()["mono_b"], VIOLET)
    txt(d, 154, 636, "source: mcp/tools.json; mcp/README.md", fonts()["small"], DIM)

    cx, cy, scale = 1440, 560, 150.0 / 1.5e-8
    grow = min(1.0, max(0.0, (t - 1.2) / 1.2))
    rr = r_ring = 1.5e-8 * scale * grow
    r_sun, r_pl = 6e-9 * scale * grow, 4.5e-9 * scale * grow
    r_car = 1.05e-8 * scale * grow
    if grow > 0.05:
        d.ellipse([cx - r_ring, cy - r_ring, cx + r_ring, cy + r_ring],
                  outline=AMBER, width=12)
        for k in range(3):
            ang = 2 * math.pi * k / 3
            px, py = cx + r_car * math.cos(ang), cy + r_car * math.sin(ang)
            draw_gear(d, px, py, r_pl, r_pl * 0.8, 18, 0, (44, 70, 104), CYAN)
        draw_gear(d, cx, cy, r_sun, r_sun * 0.8, 24, 0, (58, 84, 122), ACCENT)

    panel(d, 130, 830, 1200, 900, PANEL, LINE)
    for i, (lbl, val) in enumerate([("sun", "24"), ("planet", "18"),
                                    ("ring", "60"), ("planets", "3")]):
        x = 180 + i * 250
        txt(d, x, 852, lbl, fonts()["small"], DIM)
        txt(d, x, 872, val, fonts()["mono_b"], GREEN)
    txt(d, 1230, 862, "source: TASKS.md M4-06", fonts()["small"], DIM)
    txt(d, 1440, 210, "the generator returns atoms and bonds", fonts()["h2"], FG)


# ---------------------------------------------------------------- scene 6
def scene_assemble(r, t, shot):
    d = r["d"]
    body = [
        '{"jsonrpc":"2.0","id":3,',
        ' "method":"tools/call",',
        ' "params":{"name":"assemble_planetary",',
        '  "arguments":{"specs":{"planet_count":3}}}}',
    ]
    json_box(r, 130, 210, 900, 540, "agent -> MCP  (JSON-RPC 2.0)", body, t, 0.2)
    txt(d, 154, 576, "handle: assembly-1", fonts()["mono_b"], VIOLET)

    cx, cy, scale = 1450, 560, 150.0 / 1.5e-8
    grow = min(1.0, max(0.0, (t - 1.0) / 1.5))
    if grow > 0.05:
        planetary_art(r, cx, cy, scale * grow, t, show_joints=True)

    panel(d, 130, 640, 900, 900, PANEL, LINE)
    txt(d, 158, 664, "L2 device", fonts()["small_b"], CYAN)
    txt(d, 158, 700, "7 bodies", fonts()["h1"], GREEN)
    txt(d, 500, 700, "6 joints", fonts()["h1"], GREEN)
    txt(d, 158, 760, "constraints: position-level projection", fonts()["body"], FG)
    txt(d, 158, 806, "no velocity-level pass  (ADR-0014)", fonts()["small"], DIM)
    txt(d, 158, 846, "source: TASKS.md M7-04, M8-01", fonts()["small"], DIM)


# ---------------------------------------------------------------- scene 7
def scene_measure(r, t, shot):
    d = r["d"]
    body = [
        '{"jsonrpc":"2.0","id":4,',
        ' "method":"tools/call",',
        ' "params":{"name":"measure_gear_ratio",',
        '  "arguments":{"assembly_id":"assembly-1"}}}',
    ]
    json_box(r, 130, 210, 900, 540, "agent -> MCP  (JSON-RPC 2.0)", body, t, 0.2)

    cx, cy, scale = 1450, 540, 150.0 / 1.5e-8
    settled = t > 4.2
    planetary_art(r, cx, cy, scale, t, show_ratio=None)

    panel(d, 130, 640, 900, 900, PANEL, LINE)
    txt(d, 158, 662, "gear ratio  (sun -> carrier)", fonts()["small_b"], CYAN)
    txt(d, 158, 694, "analytic", fonts()["body"], DIM)
    txt(d, 560, 694, "3.500000000", fonts()["mono_b"], FG)
    txt(d, 158, 730, "measured", fonts()["body"], DIM)
    txt(d, 560, 730, "3.500000000" if settled else "measuring...",
        fonts()["mono_b"], GREEN if settled else AMBER)
    txt(d, 158, 778, "relative error 1.146e-10", fonts()["mono"], GREEN)
    txt(d, 158, 814, "tolerance 1.0e-3", fonts()["mono"], DIM)
    txt(d, 158, 856, "source: TASKS.md M6-06; docs/reproduce.md",
        fonts()["small"], DIM)
    # progress bar
    frac = min(1.0, t / 4.2)
    d.rounded_rectangle([920, 660, 1180, 682], radius=11, fill=PANEL2)
    if frac > 0:
        d.rounded_rectangle([920, 660, 920 + int(260 * frac), 682], radius=11,
                            fill=GREEN if settled else AMBER)


# ---------------------------------------------------------------- scene 8
def scene_urdf(r, t, shot):
    d = r["d"]
    body = [
        '{"jsonrpc":"2.0","id":5,',
        ' "method":"tools/call",',
        ' "params":{"name":"export_urdf",',
        '  "arguments":{"assembly_id":"assembly-1"}}}',
    ]
    json_box(r, 120, 200, 820, 500, "agent -> MCP  (JSON-RPC 2.0)", body, t, 0.2)
    panel(d, 120, 560, 820, 900, PANEL, LINE)
    txt(d, 148, 584, "URDF", fonts()["small_b"], CYAN)
    txt(d, 148, 620, "7 links", fonts()["h1"], GREEN)
    txt(d, 460, 620, "6 joints", fonts()["h1"], GREEN)
    txt(d, 148, 680, "source: TASKS.md M7-04, M8-01", fonts()["small"], DIM)
    txt(d, 148, 720, "round-trips bodies, joints, axes", fonts()["body"], FG)
    txt(d, 148, 766, "(parse_urdf  /  format_urdf)", fonts()["small"], DIM)
    txt(d, 148, 820, "simulated", fonts()["small_b"], AMBER)

    # URDF tree
    tree_x, tree_y = 960, 210
    nodes = [("ground", tree_y), ("sun", tree_y + 78), ("planet 1", tree_y + 156),
             ("planet 2", tree_y + 234), ("planet 3", tree_y + 312),
             ("ring", tree_y + 390), ("carrier", tree_y + 468)]
    txt(d, tree_x - 30, tree_y - 46, "URDF tree", fonts()["h2"], FG)
    d.line([tree_x + 40, tree_y + 10, tree_x + 40, nodes[-1][1] + 10],
           fill=LINE, width=3)
    for name, y in nodes:
        d.line([tree_x + 40, y + 10, tree_x + 90, y + 10], fill=LINE, width=3)
        box_fill = PANEL2
        panel(d, tree_x + 90, y - 22, tree_x + 360, y + 42, box_fill, LINE, r=10)
        txt(d, tree_x + 112, y + 10, name, fonts()["body"], FG, "lm")
    txt(d, tree_x - 30, 760, "7 links, 6 joints", fonts()["h1"], GREEN)

    # three-scale crossfade
    cf = (t - 3.0) / 6.0
    labels = [("atomistic", CYAN), ("device", ACCENT), ("coarse", AMBER)]
    idx = min(2, max(0, int(cf * 3)))
    lx = 1360
    for i, (lbl, c) in enumerate(labels):
        on = i == idx
        txt(d, lx + i * 150, 860, lbl, fonts()["small_b"], c if on else DIM)
    cx, cy, scale = 1580, 560, 90.0 / 1.5e-8
    if idx == 0:
        for k in range(90):
            a = k * 2.399963
            rad = 190 * math.sqrt((k + 0.5) / 90)
            x, y = cx + rad * math.cos(a), cy + rad * math.sin(a)
            d.ellipse([x - 5, y - 5, x + 5, y + 5], fill=CYAN)
        txt(d, cx, cy + 230, "static atom snapshot", fonts()["small"], DIM, "ma")
    elif idx == 1:
        planetary_art(r, cx, cy, scale, t)
    else:
        r_ring = 1.5e-8 * scale
        d.ellipse([cx - r_ring, cy - r_ring, cx + r_ring, cy + r_ring],
                  outline=AMBER, width=4)
        d.ellipse([cx - 6e-9 * scale, cy - 6e-9 * scale,
                   cx + 6e-9 * scale, cy + 6e-9 * scale], outline=AMBER, width=3)
        for k in range(3):
            a = 2 * math.pi * k / 3
            px, py = cx + 1.05e-8 * scale * math.cos(a), cy + 1.05e-8 * scale * math.sin(a)
            d.ellipse([px - 4.5e-9 * scale, py - 4.5e-9 * scale,
                       px + 4.5e-9 * scale, py + 4.5e-9 * scale],
                      outline=AMBER, width=3)
        txt(d, cx, cy + 230, "pitch circles + bounding cylinders", fonts()["small"], DIM, "ma")
    txt(d, 120, 940, "schematic three-scale cross-fade; source: docs/three-scale.md",
        fonts()["small"], DIM)


# ---------------------------------------------------------------- scene 9
def scene_bench(r, t, shot):
    d = r["d"]
    x0, y0, x1, y1 = 120, 200, 760, 640
    term_window(r, x0, y0, x1, y1, "cargo bench -p nanocad-engine")
    rows = [
        "$ cargo bench -p nanocad-engine",
        "    Finished bench profile [optimized]",
        "water_200_energy_and_gradient  time: [677.92 us 680.97 us 685.88 us]",
        "water_200_forces_serial        time: [678.45 us 679.83 us 681.14 us]",
        "water_200_forces_parallel_4    time: [369.07 us 370.97 us 373.05 us]",
    ]
    n = int(min(len(rows), max(0.0, (t - 0.3) * 1.6)))
    yy = y0 + 78
    for i, s in enumerate(rows[:n]):
        txt(d, x0 + 22, yy, s[:92], fonts()["mono_s"], GREEN if i else ACCENT)
        yy += 40
    txt(d, x0, y1 + 20, "./scripts/bench.sh  writes one record per run",
        fonts()["small"], DIM)

    # bar chart
    bx, by, bw, bh = 880, 250, 880, 470
    panel(d, bx - 30, by - 40, bx + bw + 30, by + bh + 70, PANEL, LINE)
    txt(d, bx, by - 30, "mean time  (lower is better)", fonts()["h2"], FG)
    txt(d, bx + bw, by - 30, "simulated", fonts()["small_b"], AMBER, "ra")
    bars = [("energy", 680.97, ACCENT), ("serial", 679.83, ACCENT),
            ("4-worker", 370.97, GREEN)]
    maxv = 720.0
    grow = min(1.0, max(0.0, (t - 1.4) / 1.4))
    for i, (name, val, c) in enumerate(bars):
        x = bx + 60 + i * 290
        hgt = int((val / maxv) * bh * grow)
        d.rectangle([x, by + bh - hgt, x + 170, by + bh], fill=c)
        d.line([x, by + bh, x + 170, by + bh], fill=FG, width=3)
        if grow > 0.6:
            txt(d, (x + 85, by + bh - hgt - 40), f"{val:.2f} us", fonts()["mono_b"], c, "ma")
        txt(d, (x + 85, by + bh + 14), name, fonts()["body"], FG, "ma")
    txt(d, bx, by + bh + 40, "Apple M4, 10 cores  -  cargo 1.98.1", fonts()["small"], FG)
    panel(d, bx, by + bh + 74, bx + 560, by + bh + 122, PANEL, AMBER, width=3, r=10)
    txt(d, bx + 20, by + bh + 98, "one host; your numbers will differ", fonts()["small_b"], AMBER, "lm")
    txt(d, bx + 600, by + bh + 98,
        "source: 2026-09-14-engine-timings.txt", fonts()["small"], DIM, "lm")


# ---------------------------------------------------------------- scene 10
def scene_closing(r, t, shot):
    d = r["d"]
    d.rectangle([0, 0, W, H], fill=BG)
    txt(d, (W / 2, 330), "nano-cad", fonts()["title"], FG, "ma")
    txt(d, (W / 2, 420), "planetary gearbox demo", fonts()["h1"], DIM, "ma")
    panel(d, W / 2 - 300, 500, W / 2 + 300, 600, PANEL, AMBER, width=3)
    txt(d, (W / 2, 550), "SIMULATED", fonts()["title"], AMBER, "mm")
    txt(d, (W / 2, 660), "All results are simulated. No physical device exists.",
        fonts()["h1"], FG, "ma")
    txt(d, (W / 2, 722), "A design hypothesis, not a validated device.",
        fonts()["body"], DIM, "ma")
    txt(d, (W / 2, 770), "OpenMM cross-check: a coding check, not a physics validation.",
        fonts()["body"], DIM, "ma")
    txt(d, (W / 2, 860), REPO_URL, fonts()["body_b"], ACCENT, "ma")
    txt(d, (W / 2, 910), "Apache-2.0", fonts()["body"], DIM, "ma")


SCENES = {
    "gears": scene_gears,
    "levels": scene_levels,
    "reproduce": scene_reproduce,
    "tests": scene_tests,
    "generate": scene_generate,
    "assemble": scene_assemble,
    "measure": scene_measure,
    "urdf": scene_urdf,
    "bench": scene_bench,
    "closing": scene_closing,
}


# ---------------------------------------------------------------- chrome
def draw_chrome(r, t, shot):
    d = r["d"]
    f = fonts()
    d.rectangle([0, 0, W, 70], fill=(8, 11, 18))
    d.line([0, 70, W, 70], fill=LINE, width=2)
    txt(d, 40, 35, "nano-cad", f["h2"], FG, "lm")
    txt(d, 180, 35, "|  gearbox demo", f["body"], DIM, "lm")
    # SIMULATED badge, always visible next to a result
    panel(d, W - 260, 16, W - 40, 56, (52, 38, 12), AMBER, r=10, width=2)
    txt(d, (W - 150, 36), "SIMULATED", f["small_b"], AMBER, "mm")
    if shot["scene"] not in ("reproduce", "tests", "bench", "closing"):
        txt(d, (W - 300, 36), "schematic", f["small"], CYAN, "rm")

    # shot title
    panel(d, 40, 96, 40 + int(d.textlength(shot["title"], font=f["h2"])) + 40, 148,
          PANEL, LINE, r=10)
    txt(d, 64, 122, shot["title"], f["h2"], FG, "lm")

    # caption bar
    cap = shot["narration"]
    lines = wrap(cap, f["caption"], W - 220, d)
    ch = 30 + len(lines) * 42
    y0 = H - 20 - ch
    d.rounded_rectangle([70, y0, W - 70, H - 20], radius=14,
                        fill=(6, 9, 15), outline=LINE, width=2)
    yy = y0 + 24
    for ln in lines:
        txt(d, (W / 2, yy + 16), ln, f["caption"], FG, "mm")
        yy += 42


def render_frame(r, frame):
    t = frame / FPS
    shot = SHOTS[0]
    for s in SHOTS:
        if s["start"] <= t < s["end"]:
            shot = s
            break
    local = t - shot["start"]
    img = r["img"]
    d = r["d"]
    d.rectangle([0, 0, W, H], fill=BG)
    r["d"] = d
    SCENES[shot["scene"]](r, local, shot)
    r["d"] = d
    draw_chrome(r, t, shot)
    return img


def hms(sec):
    m = int(sec // 60)
    s = sec - m * 60
    return f"{m:02d}:{s:06.3f}".replace(".", ",")


def write_srt(path):
    out = []
    for i, s in enumerate(SHOTS, 1):
        out.append(str(i))
        out.append(f"{hms(s['start'])} --> {hms(s['end'])}")
        lines = textwrap.wrap(s["narration"], width=58)
        out.extend(lines)
        out.append("")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(out))


# ---------------------------------------------------------------- audio
def run(cmd, **kw):
    return subprocess.run(cmd, check=True, **kw)


def build_audio(build, ffmpeg, ffprobe, say, use_say=True):
    os.makedirs(build, exist_ok=True)
    clips = []
    silent_shots = []
    for s in SHOTS:
        win = s["end"] - s["start"]
        raw = os.path.join(build, f"n{s['n']}.aiff")
        wav = os.path.join(build, f"n{s['n']}.wav")
        ok = False
        if use_say:
            try:
                run([say, "-v", "Samantha", "-r", "180", "-o", raw, s["narration"]])
                ok = os.path.exists(raw)
            except Exception as e:
                print(f"  say failed on shot {s['n']}: {e}", file=sys.stderr)
                ok = False
        if not ok:
            silent_shots.append(s["n"])
            run([ffmpeg, "-y", "-v", "error", "-f", "lavfi",
                 "-i", "anullsrc=r=44100:cl=mono", "-t", str(win), wav])
            clips.append(wav)
            continue
        dur = float(subprocess.run(
            [ffprobe, "-v", "error", "-show_entries", "format=duration",
             "-of", "default=noprint_wrappers=1:nokey=1", raw],
            capture_output=True, text=True, check=True).stdout.strip())
        margin = win - 0.25
        af = []
        if dur > margin:
            tempo = min(2.5, dur / margin)
            af.append(f"atempo={tempo:.4f}")
        af.append(f"apad=whole_dur={win}")
        run([ffmpeg, "-y", "-v", "error", "-i", raw,
             "-af", ",".join(af), "-t", f"{win}", "-ar", "44100", "-ac", "1",
             "-c:a", "pcm_s16le", wav])
        clips.append(wav)
    # concatenate exactly; the sum of windows is 90 s
    audio = os.path.join(build, "narration.wav")
    inputs = []
    for c in clips:
        inputs += ["-i", c]
    filt = "".join(f"[{i}:a]" for i in range(len(clips))) + \
        f"concat=n={len(clips)}:v=0:a=1[out]"
    run([ffmpeg, "-y", "-v", "error", *inputs, "-filter_complex", filt,
         "-map", "[out]", "-ar", "44100", "-ac", "1", "-c:a", "pcm_s16le", audio])
    return audio, silent_shots


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--srt", required=True)
    ap.add_argument("--build", required=True)
    ap.add_argument("--ffmpeg", default="/opt/homebrew/bin/ffmpeg")
    ap.add_argument("--ffprobe", default="/opt/homebrew/bin/ffprobe")
    ap.add_argument("--say", default="/usr/bin/say")
    ap.add_argument("--no-say", action="store_true")
    args = ap.parse_args()

    if not os.path.exists(args.ffmpeg):
        sys.exit(f"error: ffmpeg not found at {args.ffmpeg}")
    if not os.path.exists(args.ffprobe):
        sys.exit(f"error: ffprobe not found at {args.ffprobe}")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    os.makedirs(args.build, exist_ok=True)

    if os.path.isdir(args.build):
        for name in os.listdir(args.build):
            p = os.path.join(args.build, name)
            if os.path.isfile(p):
                os.remove(p)

    print("building narration audio...", file=sys.stderr)
    audio, silent = build_audio(args.build, args.ffmpeg, args.ffprobe, args.say,
                                use_say=not args.no_say)
    if silent:
        print(f"WARNING: no TTS for shots {silent}; those shots are silent.",
              file=sys.stderr)

    print("writing srt...", file=sys.stderr)
    write_srt(args.srt)

    img = Image.new("RGB", (W, H), BG)
    r = {"img": img, "d": ImageDraw.Draw(img)}
    fonts()

    nframes = int(round(TOTAL * FPS))
    cmd = [
        args.ffmpeg, "-y", "-v", "error",
        "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", f"{W}x{H}", "-r", str(FPS),
        "-i", "-",
        "-i", audio,
        "-map", "0:v:0", "-map", "1:a:0",
        "-c:v", "libx264", "-preset", "medium", "-crf", "23",
        "-pix_fmt", "yuv420p",
        "-c:a", "aac", "-b:a", "128k",
        "-t", f"{TOTAL}",
        "-movflags", "+faststart",
        args.out,
    ]
    print(f"rendering {nframes} frames and encoding...", file=sys.stderr)
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE)
    try:
        for i in range(nframes):
            frame = render_frame(r, i)
            proc.stdin.write(frame.tobytes())
            if i % 150 == 0:
                print(f"  frame {i}/{nframes}", file=sys.stderr)
    except BrokenPipeError:
        pass
    finally:
        if proc.stdin:
            proc.stdin.close()
    rc = proc.wait()
    if rc != 0:
        sys.exit(f"error: ffmpeg exited {rc}")
    print(f"done: {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
