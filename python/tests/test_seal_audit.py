"""The seal audit finds a gap in a wall and passes a closed wall."""

import json
import os
import subprocess
import sys

import pytest

pytest.importorskip("scipy")

SCRIPT = os.path.join(os.path.dirname(__file__), "..", "..", "scripts", "seal_audit.py")


def wall_scene(path, gap_a):
    """A carbon wall on x = 0 (body 0) with a square hole, and a carbon plug
    (body 1) in the hole that leaves `gap_a` of clearance on each side."""
    atoms = []
    spacing = 1.5
    hole = 6.0
    for iy in range(-16, 17):
        for iz in range(-16, 17):
            y, z = iy * spacing, iz * spacing
            for x in (-0.75, 0.75):
                if abs(y) < hole and abs(z) < hole:
                    if abs(y) < hole - gap_a and abs(z) < hole - gap_a:
                        atoms.append((x, y, z, 1))
                    continue
                atoms.append((x, y, z, 0))
    # A blood-side and a tank-side marker atom far out, so the box spans both.
    atoms.append((-12.0, 0.0, 0.0, 1))
    atoms.append((12.0, 0.0, 0.0, 1))
    scene = {
        "scene": {
            "atomistic": {
                "atoms": [
                    {"position_m": [x * 1e-10, y * 1e-10, z * 1e-10], "atomic_number": 6, "body": b}
                    for x, y, z, b in atoms
                ]
            }
        }
    }
    with open(path, "w") as handle:
        json.dump(scene, handle)


def run(path, probe):
    out = subprocess.run(
        [sys.executable, SCRIPT, str(path), "--probe", str(probe), "--keep-drain"],
        capture_output=True,
        text=True,
    )
    return out.stdout


def test_a_wide_gap_leaks_and_a_closed_wall_seals(tmp_path):
    wide = tmp_path / "wide.json"
    wall_scene(wide, gap_a=6.0)
    assert "LEAK" in run(wide, 1.4)
    closed = tmp_path / "closed.json"
    wall_scene(closed, gap_a=0.0)
    assert "SEALED" in run(closed, 1.4)
