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
        "--energy",
        choices=["o2", "water"],
        help="instead of a hard probe, find the lowest Lennard-Jones energy barrier of a one-site "
        "probe (O2: sigma 3.467 A, eps 0.887 kJ/mol; water: SPC/E O) on any blood-to-tank path",
    )
    parser.add_argument(
        "--pair",
        nargs=6,
        type=float,
        metavar=("X1", "Y1", "Z1", "X2", "Y2", "Z2"),
        help="with --energy: the barrier between two points (A) instead of blood to tank",
    )
    parser.add_argument(
        "--box",
        nargs=6,
        type=float,
        metavar=("XMIN", "XMAX", "YMIN", "YMAX", "ZMIN", "ZMAX"),
        help="restrict the grid to this box (A); faster local checks",
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
    if args.box:
        lo = np.array(args.box[0::2], dtype=float)
        hi = np.array(args.box[1::2], dtype=float)
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
        # The drain bore: a clear line along x through the wall at z = 0,
        # beyond the rotor chamber (lined or H-capped).
        from scipy.spatial import cKDTree

        tree = cKDTree(xyz[wall])
        rotor_reach = np.hypot(xyz[body == 1][:, 0], xyz[body == 1][:, 1]).max()
        best = None
        for y in np.arange(rotor_reach + 10.0, xyz[wall][:, 1].max(), 0.25):
            for sign in (1.0, -1.0):
                line = np.column_stack([np.linspace(-10, 10, 41), np.full(41, sign * y), np.zeros(41)])
                clear = tree.query(line)[0].min()
                if best is None or clear > best[0]:
                    best = (clear, sign * y)
        if best and best[0] > 2.0:
            dy = best[1]
            ys = lo[1] + np.arange(shape[1]) * args.grid
            band = np.abs(ys - dy) < 10.0
            blocked[:, band, :] = True
            print(f"drain bore masked: y = {dy / 10:.2f} nm +- 1 nm (axis clearance {best[0]:.2f} A)")
        else:
            print("no drain bore found")

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
    if args.energy:
        return energy_barrier(args, xyz, znum, blocked, lo, shape)
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


# UFF sigma (A) and epsilon (kJ/mol) of the solid atoms.
UFF_LJ = {1: (2.571, 0.184), 6: (3.431, 0.439), 7: (3.261, 0.289), 8: (3.118, 0.251), 26: (2.594, 0.054)}
PROBES = {"o2": (3.467, 0.887), "water": (3.166, 0.650)}
KT_KJ = 0.0083145 * 310.0


def energy_barrier(args, xyz, znum, blocked, lo, shape):
    """Bisect the energy at which the region U < E first joins blood to tank.

    `blocked` here comes from a small hard probe, so it only removes voxels
    deep inside atoms; the energy decides the rest. Voxels farther than the
    cutoff from every atom have U = 0 (bulk).
    """
    from scipy.spatial import cKDTree

    sig_p, eps_p = PROBES[args.energy]
    sig = np.array([0.5 * (UFF_LJ.get(int(n), UFF_LJ[6])[0] + sig_p) for n in znum])
    eps = np.array([np.sqrt(UFF_LJ.get(int(n), UFF_LJ[6])[1] * eps_p) for n in znum])
    tree = cKDTree(xyz)
    grid = args.grid
    free_idx = np.argwhere(~blocked)
    cutoff = 8.0
    u = np.zeros(blocked.shape, dtype=np.float32)
    chunk = 100_000
    for start in range(0, len(free_idx), chunk):
        idx = free_idx[start:start + chunk]
        pts = lo + idx * grid
        d, j = tree.query(pts, k=96, distance_upper_bound=cutoff)
        ok = np.isfinite(d)
        jj = np.where(ok, j, 0)
        dd = np.where(ok, d, cutoff)
        sr6 = (sig[jj] / dd) ** 6
        e = np.where(ok, 4.0 * eps[jj] * (sr6 * sr6 - sr6), 0.0).sum(axis=1)
        u[idx[:, 0], idx[:, 1], idx[:, 2]] = e
    u[blocked] = np.inf

    def voxel(pt):
        # The lowest-energy voxel within 3 A of the point.
        c = np.round((np.array(pt) - lo) / grid).astype(int)
        r = int(round(3.0 / grid))
        sl = tuple(slice(max(0, c[i] - r), min(u.shape[i], c[i] + r + 1)) for i in range(3))
        sub = u[sl]
        k = np.unravel_index(np.argmin(sub), sub.shape)
        return tuple(sl[i].start + k[i] for i in range(3))

    def joined(e_max):
        labels, _ = ndimage.label(u < e_max)
        if args.pair:
            a_lab = labels[voxel(args.pair[:3])]
            b_lab = labels[voxel(args.pair[3:])]
            ok = a_lab != 0 and a_lab == b_lab
            return ok, labels, {a_lab} if ok else set()
        blood = set(np.unique(labels[-1, :, :])) - {0}
        tank = set(np.unique(labels[0, :, :])) - {0}
        return bool(blood & tank), labels, blood & tank

    if args.pair:
        for name, pt in (("from", args.pair[:3]), ("to", args.pair[3:])):
            v = voxel(pt)
            print(f"  {name} {pt} -> {np.round(lo + np.array(v) * grid, 1)}: U {u[v]:.2f} kJ/mol")
    lo_e, hi_e = -40.0, 200.0
    if not joined(hi_e)[0]:
        print(f"{args.energy}: no path below {hi_e} kJ/mol")
        return 0
    while hi_e - lo_e > 0.1:
        mid = 0.5 * (lo_e + hi_e)
        if joined(mid)[0]:
            hi_e = mid
        else:
            lo_e = mid
    _, labels, comps = joined(hi_e)
    # The saddle: voxels of the joining component within kT of the barrier.
    comp = np.isin(labels, list(comps))
    saddle = comp & (u > hi_e - KT_KJ)
    sx, sy, sz = np.nonzero(saddle)
    pts = lo + np.stack([sx, sy, sz], 1) * grid
    print(f"{args.energy}: barrier {hi_e:.1f} kJ/mol ({hi_e / KT_KJ:.1f} kT) above bulk")
    # Where the joining component crosses the wall plane, by radius from
    # the rotor axis (z).
    xs_ = lo[0] + np.arange(u.shape[0]) * grid
    band_x = np.abs(xs_) < 1.0
    if band_x.any() and not args.pair:
        sub = comp[band_x]
        bi, by, bz = np.nonzero(sub)
        yy, zz = lo[1] + by * grid, lo[2] + bz * grid
        rr_ = np.abs(yy)
        for name, sel in (("r < 33 A (rotor interior)", rr_ < 33), ("33-41 A (rim zone)", (rr_ >= 33) & (rr_ < 41)), ("r >= 41 A (wall)", rr_ >= 41)):
            if sel.any():
                print(f"  crossing x=0, {name}: {int(sel.sum())} voxels, |y| {rr_[sel].min():.1f}..{rr_[sel].max():.1f}, z {zz[sel].min():.1f}..{zz[sel].max():.1f}")
    plane = int(round((0.0 - lo[0]) / grid))
    cy, cz = (np.nonzero(comp[plane]) if 0 <= plane < u.shape[0] else (np.array([]), np.array([])))
    if len(cy) and not args.pair:
        py, pz = lo[1] + cy * grid, lo[2] + cz * grid
        near = tree.query(np.column_stack([np.zeros(len(py)), py, pz]))[0]
        print(f"  at the wall plane: {len(cy)} voxels, y {py.min():.1f}..{py.max():.1f} A, "
              f"z {pz.min():.1f}..{pz.max():.1f} A, nearest atom {near.min():.2f}..{near.max():.2f} A")
        order = np.argsort(near)[:5]
        for k in order:
            v = lo + np.array([plane, cy[k], cz[k]]) * grid
            dd, jj = tree.query(v, k=400, distance_upper_bound=cutoff)
            okk = np.isfinite(dd)
            sr = (sig[jj[okk]] / dd[okk]) ** 6
            direct = float((4.0 * eps[jj[okk]] * (sr * sr - sr)).sum())
            print(f"    voxel {np.round(v, 2)}: nearest atom {tree.query(v)[0]:.2f} A, "
                  f"U grid {u[plane, cy[k], cz[k]]:.2f}, U direct {direct:.2f}")
    if len(pts):
        th = np.degrees(np.arctan2(pts[:, 1], pts[:, 0])) % 360
        print(f"  saddle voxels {len(pts)}: x {pts[:, 0].min():.1f}..{pts[:, 0].max():.1f} A, "
              f"azimuth {th.min():.0f}..{th.max():.0f} deg, z {pts[:, 2].min():.1f}..{pts[:, 2].max():.1f} A")
    return 0


if __name__ == "__main__":
    sys.exit(main())
