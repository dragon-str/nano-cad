#!/usr/bin/env python3
"""Seal audit of the respirocyte hull-wall section (ADR-0117).

The hull wall lies on the plane x = 0: blood on the +x side, the tank on the
-x side. The rotor turns in a hole through the wall. This audit asks one
geometric question of a relaxed scene: can a rigid probe sphere travel from
the blood side to the tank side through the free space, with the rotor frozen
in its scene pose?

Method. Voxelize the box of the wall section. A voxel is blocked when its
centre is closer to an atom than the atom's van der Waals radius plus the
probe radius. Label the free voxels (6-connected). A component that touches
both the +x face and the -x face of the box is a leak. The box stops at the
slab's y and z extent, because the real hull continues past the section.
The drain bore is meant to join the two sides (ADR-0097, ADR-0098), so it is
masked by default.

For each leak the script reports where it crosses the wall plane, classed by
height: the disc interior (|z| < disc half-thickness), the face gap, or the
cam groove.

This is a static, rigid-sphere estimate. Thermal motion and soft contacts let
molecules pass narrower gaps than this; a path it finds is real, and a pass
here is necessary but not sufficient.

Usage:
    python scripts/seal_audit.py SCENE.json [--probe 1.73] [--grid 0.5]
        [--radius-scale 1.0] [--keep-drain] [--dry-run]

Run with an environment that has numpy and scipy
(for example ~/.cache/nano-cad/md/bin/python).
"""

import argparse
import json
import sys

import numpy as np
from scipy import ndimage

# Bondi van der Waals radii, in angstrom. Iron: an estimate.
VDW_A = {1: 1.20, 6: 1.70, 7: 1.55, 8: 1.52, 26: 2.00}
DISC_HALF_A = 12.0
CHAMBER_HALF_A = 19.0
SLAB_HALF_A = 12.0


def load(path):
    with open(path) as handle:
        data = json.load(handle)
    scene = data.get("scene", data)
    atoms = scene["atomistic"]["atoms"]
    xyz = np.array([a["position_m"] for a in atoms]) * 1e10
    z = np.array([a["atomic_number"] for a in atoms])
    body = np.array([a["body"] for a in atoms])
    return xyz, z, body


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("scene")
    parser.add_argument("--probe", type=float, default=1.73, help="probe radius, A (O2 1.73, water 1.40)")
    parser.add_argument("--grid", type=float, default=0.5, help="voxel size, A")
    parser.add_argument("--radius-scale", type=float, default=1.0, help="scale on every vdW radius")
    parser.add_argument("--keep-drain", action="store_true", help="do not mask the drain bore")
    parser.add_argument("--dry-run", action="store_true", help="print the grid plan and stop")
    parser.add_argument(
        "--mask",
        action="append",
        default=[],
        choices=["faces", "disc"],
        help="block a region to isolate a path: faces (|z| >= disc half inside the plates), "
        "disc (|z| < disc half inside the chamber)",
    )
    parser.add_argument(
        "--where",
        nargs=2,
        type=float,
        metavar=("R_MIN_A", "R_MAX_A"),
        help="print where each leak occupies the band R_MIN..R_MAX (|z| >= disc half), by azimuth",
    )
    args = parser.parse_args()

    xyz, znum, body = load(args.scene)
    wall = body == 0
    lo = xyz[wall].min(axis=0)
    hi = xyz[wall].max(axis=0)
    # x: the whole section plus a margin; y, z: the slab extent.
    margin = 6.0
    # Stay 3 A inside the slab's y and z extent, so no path runs along the
    # box edge outside the last lattice layer.
    lo = np.array([xyz[:, 0].min() - margin, lo[1] + 3.0, lo[2] + 3.0])
    hi = np.array([xyz[:, 0].max() + margin, hi[1] - 3.0, hi[2] - 3.0])
    shape = np.ceil((hi - lo) / args.grid).astype(int) + 1
    voxels = int(np.prod(shape))
    print(f"atoms {len(xyz)}, grid {args.grid} A, shape {tuple(shape)}, {voxels / 1e6:.1f} M voxels")
    print(f"memory plan: about {voxels * 6 / 1e9:.2f} GB (bool + int32 labels + mask)")
    if args.dry_run:
        return 0

    blocked = np.zeros(shape, dtype=bool)
    radii = np.array([VDW_A.get(int(n), 1.7) for n in znum]) * args.radius_scale + args.probe
    for r in np.unique(radii):
        sel = xyz[radii == r]
        idx = np.round((sel - lo) / args.grid).astype(int)
        reach = int(np.ceil(r / args.grid))
        span = np.arange(-reach, reach + 1)
        offsets = np.array(np.meshgrid(span, span, span, indexing="ij")).reshape(3, -1).T
        offsets = offsets[np.sum((offsets * args.grid) ** 2, axis=1) <= r * r]
        for off in offsets:
            p = idx + off
            ok = np.all((p >= 0) & (p < shape), axis=1)
            p = p[ok]
            blocked[p[:, 0], p[:, 1], p[:, 2]] = True

    if not args.keep_drain:
        oxy = xyz[wall & (znum == 8)]
        if len(oxy):
            dy = oxy[:, 1].mean()
            ys = lo[1] + np.arange(shape[1]) * args.grid
            band = np.abs(ys - dy) < 10.0
            blocked[:, band, :] = True
            print(f"drain bore masked: y = {dy / 10:.2f} nm +- 1 nm")

    xs = lo[0] + np.arange(shape[0]) * args.grid
    ys = lo[1] + np.arange(shape[1]) * args.grid
    zs = lo[2] + np.arange(shape[2]) * args.grid
    gx, gy, gz = np.meshgrid(xs, ys, zs, indexing="ij", sparse=True)
    rr = np.sqrt(gx**2 + gy**2)
    rotor_r = np.hypot(xyz[body == 1][:, 0], xyz[body == 1][:, 1]).max()
    for name in args.mask:
        if name == "faces":
            blocked |= (rr < rotor_r + 15.0) & (np.abs(gz) >= DISC_HALF_A)
        if name == "disc":
            blocked |= (rr < rotor_r + 5.0) & (np.abs(gz) < DISC_HALF_A)
        print(f"masked: {name}")
    free = ~blocked
    labels, count = ndimage.label(free)
    blood = set(np.unique(labels[-1, :, :])) - {0}
    tank = set(np.unique(labels[0, :, :])) - {0}
    leaks = sorted(blood & tank)
    print(f"probe {args.probe} A, radius scale {args.radius_scale}: {count} free components")
    if not leaks:
        print("SEALED: no free path from the blood side to the tank side")
        return 0
    plane = np.argmin(np.abs(xs))
    for leak in leaks:
        cut = labels[plane] == leak
        iy, iz = np.nonzero(cut)
        y, z = ys[iy], zs[iz]
        classes = {
            "disc interior (|z| < 12 A)": np.abs(z) < DISC_HALF_A,
            "face gap (12-19 A)": (np.abs(z) >= DISC_HALF_A) & (np.abs(z) < CHAMBER_HALF_A),
            "cam groove (|z| >= 19 A)": np.abs(z) >= CHAMBER_HALF_A,
        }
        print(f"LEAK component {leak}: {int((labels == leak).sum())} voxels, {len(y)} at the wall plane")
        for name, sel in classes.items():
            if sel.any():
                print(
                    f"  {name}: {int(sel.sum())} voxels, y {y[sel].min() / 10:.2f}..{y[sel].max() / 10:.2f} nm, "
                    f"z {z[sel].min() / 10:.2f}..{z[sel].max() / 10:.2f} nm"
                )
        if args.where:
            comp = labels == leak
            band = (rr >= args.where[0]) & (rr <= args.where[1]) & (np.abs(gz) >= DISC_HALF_A)
            ix, iy, iz = np.nonzero(comp & band)
            if len(ix):
                th = np.degrees(np.arctan2(ys[iy], xs[ix])) % 360
                hist, edges = np.histogram(th, bins=24, range=(0, 360))
                print(f"  in band r {args.where[0]}..{args.where[1]} A: {len(ix)} voxels, z {zs[iz].min():.1f}..{zs[iz].max():.1f} A")
                print("  by azimuth (15 deg bins): " + " ".join(f"{int(e)}:{h}" for e, h in zip(edges[:-1], hist) if h))
    return 1


if __name__ == "__main__":
    sys.exit(main())
