"""Drain transit test (M28-10, M28-35): count water and oxygen that cross a lined
single-file channel, with and without a pressure difference.

Geometry from `channel_json`: a frozen diamond slab with one lined
through-channel. SPC/E water fills both sides (periodic box). The tank side
(below the slab) also holds dissolved oxygen at about 1 M, the solubility at
roughly 1000 atm. Oxygen is kept out of the periodic return path by two
flat-bottom walls in the bulk far from the slab, so an oxygen molecule above
the slab can only have come through the channel.

Pressure: a constant force f along +z on every water oxygen (Zhu, Tajkhorshid
and Schulten 2002). The pressure difference across the slab is
dP = N_water f / A. `--dp-mpa 0` is the idle tank.

Transit counting: each molecule gets a label every sample: below, in the
channel, or above. A transit is below -> channel -> above (or the reverse)
without leaving the channel region in between for the other side.

Oxygen model: two neutral Lennard-Jones sites 1.21 A apart (sigma 3.02 A,
epsilon 0.407 kJ/mol). Solid: UFF with model charges; Lorentz-Berthelot.
One process, below 1 GB. A simulated estimate.

    python scripts/drain_transit_md.py channel.json --dp-mpa 200 --production-ps 2000
"""

from __future__ import annotations

import argparse
import json
import math
import sys

import numpy as np
import openmm as mm
from openmm import unit
from scipy.spatial import cKDTree

KCAL = 4.184
LJ_FACTOR = 2.0 ** (1.0 / 6.0)
UFF = {"C": (3.851, 0.105), "H": (2.886, 0.044), "O": (3.500, 0.060)}
SPCE_SIGMA, SPCE_EPS, SPCE_Q_O, SPCE_Q_H = 0.3166, 0.6502, -0.8476, 0.4238
OH, HH = 0.1, 0.1633
WATER_VOLUME_A3 = 29.9
O2_SIGMA, O2_EPS, O2_BOND = 0.302, 0.407, 0.121
AVOGADRO = 6.02214076e23


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("geometry")
    p.add_argument("--dp-mpa", type=float, default=0.0)
    p.add_argument("--o2-molar", type=float, default=1.0)
    p.add_argument("--equilibrate-ps", type=float, default=50.0)
    p.add_argument("--production-ps", type=float, default=2000.0)
    p.add_argument("--sample-ps", type=float, default=0.5)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--platform", default="OpenCL")
    p.add_argument("--umbrella", action="store_true", help="PMF of one O2 along the channel axis")
    p.add_argument("--window-ps", type=float, default=200.0)
    p.add_argument("--window-step-nm", type=float, default=0.1)
    p.add_argument("--reverse", action="store_true", help="walk the windows from the far side back")
    p.add_argument("--save", default=None, help="write the window series to this .npz")
    p.add_argument("--combine", nargs="+", default=None, help="WHAM over saved .npz files; no MD")
    p.add_argument("--density-scale", type=float, default=None,
                   help="pack the water this much denser than bulk; default 1 + kappa dP / 2 "
                        "(kappa 4.6e-10 /Pa, SPC/E), so the low-pressure side stays above zero")
    p.add_argument("--dt-fs", type=float, default=2.0,
                   help="time step in fs; use 1 with --mobile-oh (a 2 fs step blew up once in 50-100 ps)")
    p.add_argument("--mobile-oh", action="store_true",
                   help="let each lining hydroxyl hydrogen turn about its C-O bond (heavy atoms stay frozen)")
    p.add_argument("--profile", action="store_true",
                   help="also print the mean water count near the axis (r < 0.25 nm) in 1 A slices along the bore")
    a = p.parse_args()

    g = json.load(open(a.geometry))
    if a.combine:
        centres, series = [], []
        for path in a.combine:
            d = np.load(path, allow_pickle=True)
            centres += list(d["centres"])
            series += list(d["series"])
        return report(a, centres, series, float(d["slab_lo"]), float(d["slab_hi"]), float(d["zc"]),
                      float(d["k_z"]), float(d["r_cyl"]), g)
    atoms = g["atoms"]
    solid = np.array([x[1:4] for x in atoms]) / 10.0
    elements = [x[0] for x in atoms]
    charges = np.array([x[4] for x in atoms])
    charged = np.abs(charges) > 1e-9
    if charged.any():
        charges[charged] -= charges.sum() / charged.sum()
    half_t = g["thickness"] / 20.0
    lo, hi = solid.min(axis=0), solid.max(axis=0)
    box_x, box_y = hi[0] - lo[0] + 0.25, hi[1] - lo[1] + 0.25
    reservoir = 2.0
    box_z = (hi[2] - lo[2]) + 2 * reservoir
    shift = np.array([-lo[0] + 0.125, -lo[1] + 0.125, reservoir - lo[2]])
    solid = solid + shift
    ax, ay = shift[0], shift[1]
    zc = 0.5 * (lo[2] + hi[2]) + shift[2]
    slab_lo, slab_hi = lo[2] + shift[2], hi[2] + shift[2]

    rng = np.random.default_rng(a.seed)
    tree = cKDTree(solid)
    # The constant force splits dP about the mean box pressure, so the low
    # side sits near -dP/2 and cavitates in an NVT box at bulk density
    # (found in M28-35: a vacuum gap opened above the slab at 500 MPa).
    # Packing the water denser raises the mean pressure by about dP/2.
    if a.density_scale is None:
        a.density_scale = 1.0 + 4.6e-10 * a.dp_mpa * 1e6 / 2.0
    spacing = (WATER_VOLUME_A3 / a.density_scale) ** (1 / 3) / 10.0
    # Round, not floor, the row count on each axis and stretch the rows to
    # fill the box. Flooring on a 1.8 nm box lost up to 14 percent of the
    # water per axis and left a vapour bubble at every pressure (M28-35).
    nx, ny, nz = (max(1, round(length / spacing)) for length in (box_x, box_y, box_z))
    sx, sy, sz = box_x / nx, box_y / ny, box_z / nz
    sites = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                o = np.array([(ix + 0.5) * sx, (iy + 0.5) * sy, (iz + 0.5) * sz])
                if tree.query(o)[0] >= 0.26:
                    sites.append(o)
    sites = np.array(sites)
    # Oxygen replaces water at sites on the tank side, in bulk, 0.4-1.6 nm
    # below the slab, at the requested molarity.
    tank = np.where((sites[:, 2] < slab_lo - 0.4) & (sites[:, 2] > slab_lo - 1.6))[0]
    tank_volume_l = box_x * box_y * 1.2 * 1e-24
    n_o2 = max(1, int(round(a.o2_molar * tank_volume_l * AVOGADRO)))
    if a.umbrella:
        n_o2 = 1
    o2_sites = set(rng.choice(tank, size=min(n_o2, len(tank)), replace=False).tolist())

    # Mobile hydroxyl hydrogens (M28-35). A frozen lining fixes each OH
    # hydrogen in one arbitrary direction, and the water transport hung on
    # that choice. Each lining H keeps a stiff O-H bond and a C-O-H angle and
    # turns freely about C-O. Mass 4 amu keeps the 2 fs step stable; it does
    # not change the equilibrium orientations.
    mobile = {}
    crowded = 0
    if a.mobile_oh:
        el = np.array(elements)
        heavy = np.where(el != "H")[0]
        heavy_tree = cKDTree(solid[heavy])
        lining = {}
        for hi in np.where(el == "H")[0]:
            d_h, k_h = heavy_tree.query(solid[hi])
            o_i = int(heavy[k_h])
            if el[o_i] != "O" or d_h > 0.11:
                continue
            carbons = [c for c in tree.query_ball_point(solid[o_i], 0.16) if el[c] == "C"]
            if carbons:
                c_i = int(min(carbons, key=lambda c: np.linalg.norm(solid[c] - solid[o_i])))
                lining[int(hi)] = (o_i, c_i)
        for hi, (o_i, c_i) in lining.items():
            # A hydrogen within 0.2 nm of a carbon or a cap hydrogen sits in a
            # crowded corner of the as-built lining and stays frozen. An
            # oxygen that close is a hydrogen bond; another hydroxyl hydrogen
            # turns away.
            others = [
                j for j in tree.query_ball_point(solid[hi], 0.2)
                if j not in (hi, o_i, c_i) and el[j] != "O" and j not in lining
            ]
            if others:
                crowded += 1
                continue
            mobile[hi] = (o_i, c_i)

    system = mm.System()
    system.setDefaultPeriodicBoxVectors(mm.Vec3(box_x, 0, 0), mm.Vec3(0, box_y, 0), mm.Vec3(0, 0, box_z))
    nb = mm.NonbondedForce()
    nb.setNonbondedMethod(mm.NonbondedForce.PME)
    nb.setCutoffDistance(0.9)
    pos = list(solid)
    for idx, (e, q) in enumerate(zip(elements, charges)):
        system.addParticle(4.0 if idx in mobile else 0.0)
        x1, d1 = UFF[e]
        nb.addParticle(float(q), x1 / LJ_FACTOR / 10.0, d1 * KCAL)
    water_o, o2_a = [], []
    for k, o in enumerate(sites):
        d = rng.normal(size=3)
        d /= np.linalg.norm(d)
        if k in o2_sites:
            i = system.addParticle(15.999)
            j = system.addParticle(15.999)
            nb.addParticle(0.0, O2_SIGMA, O2_EPS)
            nb.addParticle(0.0, O2_SIGMA, O2_EPS)
            nb.addException(i, j, 0.0, 0.1, 0.0)
            system.addConstraint(i, j, O2_BOND)
            pos += [o - 0.5 * O2_BOND * d, o + 0.5 * O2_BOND * d]
            o2_a.append(i)
            continue
        r = rng.normal(size=3)
        r -= np.dot(r, d) * d
        r /= np.linalg.norm(r)
        h = math.radians(109.47 / 2)
        i = system.addParticle(15.999)
        h1 = system.addParticle(1.008)
        h2 = system.addParticle(1.008)
        nb.addParticle(SPCE_Q_O, SPCE_SIGMA, SPCE_EPS)
        nb.addParticle(SPCE_Q_H, 0.1, 0.0)
        nb.addParticle(SPCE_Q_H, 0.1, 0.0)
        for u, v in ((i, h1), (i, h2), (h1, h2)):
            nb.addException(u, v, 0.0, 0.1, 0.0)
        system.addConstraint(i, h1, OH)
        system.addConstraint(i, h2, OH)
        system.addConstraint(h1, h2, HH)
        pos += [o, o + OH * (math.cos(h) * d + math.sin(h) * r), o + OH * (math.cos(h) * d - math.sin(h) * r)]
        water_o.append(i)
    for i, j in tree.query_pairs(0.9):
        i, j = int(i), int(j)
        # A mobile H feels the rest of the solid, except its own O and C.
        if i in mobile and j not in mobile[i]:
            continue
        if j in mobile and i not in mobile[j]:
            continue
        nb.addException(i, j, 0.0, 0.1, 0.0)
    system.addForce(nb)
    if mobile:
        bonds = mm.HarmonicBondForce()
        angles = mm.HarmonicAngleForce()
        for h_i, (o_i, c_i) in mobile.items():
            r0 = float(np.linalg.norm(solid[h_i] - solid[o_i]))
            v1, v2 = solid[c_i] - solid[o_i], solid[h_i] - solid[o_i]
            th0 = float(np.arccos(np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)))
            bonds.addBond(o_i, h_i, r0, 4.6e5)
            angles.addAngle(c_i, o_i, h_i, th0, 400.0)
        system.addForce(bonds)
        system.addForce(angles)

    area_m2 = box_x * box_y * 1e-18
    force_n = a.dp_mpa * 1e6 * area_m2 / len(water_o)
    force_kj_nm = force_n / 1.66054e-12
    if a.dp_mpa:
        push = mm.CustomExternalForce("-f*z")
        push.addGlobalParameter("f", force_kj_nm)
        for i in water_o:
            push.addParticle(i, [])
        system.addForce(push)
    # Oxygen walls in the bulk, far from the slab, block the periodic path.
    wall_lo, wall_hi = slab_lo - 1.8, slab_hi + 1.8
    wall = mm.CustomExternalForce("k*(step(zl-z)*(zl-z)^2+step(z-zh)*(z-zh)^2)")
    wall.addGlobalParameter("k", 5000.0)
    wall.addGlobalParameter("zl", float(wall_lo))
    wall.addGlobalParameter("zh", float(wall_hi))
    for i in o2_a:
        wall.addParticle(i, [])
        wall.addParticle(i + 1, [])
    system.addForce(wall)

    k_z = 2000.0  # kJ/(mol nm^2)
    r_cyl = 0.5  # nm, flat-bottom radial wall on the O2 centre
    if a.umbrella:
        um = mm.CustomCentroidBondForce(
            1, "0.5*kz*(z1-z0)^2 + kr*step(r-rc)*(r-rc)^2; r=sqrt((x1-cx)^2+(y1-cy)^2)"
        )
        for name, value in (("kz", k_z), ("z0", 0.0), ("kr", 5000.0), ("rc", r_cyl), ("cx", ax), ("cy", ay)):
            um.addGlobalParameter(name, float(value))
        um.addGroup([o2_a[0], o2_a[0] + 1])
        um.addBond([0], [])
        system.addForce(um)

    integ = mm.LangevinMiddleIntegrator(310 * unit.kelvin, 1 / unit.picosecond, a.dt_fs * 1e-3 * unit.picoseconds)
    integ.setRandomNumberSeed(a.seed)
    try:
        plat = mm.Platform.getPlatformByName(a.platform)
    except Exception:  # noqa: BLE001
        plat = mm.Platform.getPlatformByName("CPU")
    ctx = mm.Context(system, integ, plat)
    ctx.setPositions(np.array(pos) * unit.nanometer)
    mm.LocalEnergyMinimizer.minimize(ctx, 10.0, 1000)
    ctx.setVelocitiesToTemperature(310 * unit.kelvin, a.seed)
    integ.step(int(a.equilibrate_ps / (a.dt_fs * 1e-3)))

    if a.umbrella:
        return umbrella(a, ctx, integ, o2_a[0], slab_lo, slab_hi, zc, k_z, r_cyl, g)

    def labels(z, x, y):
        dx = (x - ax + box_x / 2) % box_x - box_x / 2
        dy = (y - ay + box_y / 2) % box_y - box_y / 2
        zz = z % box_z
        lab = np.zeros(len(z), dtype=int)  # 0 = elsewhere (bulk far away)
        near = np.hypot(dx, dy) < 0.4
        lab[(zz < slab_lo) & (zz > slab_lo - 1.0)] = -1
        lab[(zz > slab_hi) & (zz < slab_hi + 1.0)] = 1
        lab[near & (np.abs(zz - zc) < half_t)] = 2
        return lab

    def tracker(n):
        return {"last": np.zeros(n, dtype=int), "entered_from": np.zeros(n, dtype=int), "up": 0, "down": 0}

    def update(t, lab):
        entering = (lab == 2) & (t["last"] != 2)
        t["entered_from"][entering] = np.where(t["last"][entering] == 2, 0, t["last"][entering])
        leaving = (lab != 2) & (t["last"] == 2)
        for k in np.where(leaving & (lab == 1) & (t["entered_from"] == -1))[0]:
            t["up"] += 1
        for k in np.where(leaving & (lab == -1) & (t["entered_from"] == 1))[0]:
            t["down"] += 1
        t["last"] = np.where(lab == 0, t["last"], lab)

    wt, ot = tracker(len(water_o)), tracker(len(o2_a))
    o2_in_channel, water_in_channel, o2_above = [], [], []
    # Net crossings of the mid-plane near the axis: robust to a long bore.
    edges = np.arange(slab_lo - 0.5, slab_hi + 0.5001, 0.1)
    profile = np.zeros(len(edges) - 1)
    full_edges = np.linspace(0.0, box_z, 41)
    full = np.zeros(len(full_edges) - 1)
    mid_side = None
    mid_up = mid_down = 0
    for _ in range(int(a.production_ps / a.sample_ps)):
        integ.step(int(a.sample_ps / (a.dt_fs * 1e-3)))
        xyz = ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        w = xyz[water_o]
        o = 0.5 * (xyz[o2_a] + xyz[np.array(o2_a) + 1])
        lw, lo2 = labels(w[:, 2], w[:, 0], w[:, 1]), labels(o[:, 2], o[:, 0], o[:, 1])
        update(wt, lw)
        update(ot, lo2)
        dx = (w[:, 0] - ax + box_x / 2) % box_x - box_x / 2
        dy = (w[:, 1] - ay + box_y / 2) % box_y - box_y / 2
        zz = w[:, 2] % box_z
        axial = np.hypot(dx, dy) < 0.25
        profile += np.histogram(zz[axial], bins=edges)[0]
        full += np.histogram(zz, bins=full_edges)[0]
        side = np.sign(zz - zc)
        if mid_side is not None:
            near_mid = (np.abs(zz - zc) < 0.3) & (np.hypot(dx, dy) < 0.4)
            mid_up += int(np.sum(near_mid & (mid_side < 0) & (side > 0)))
            mid_down += int(np.sum(near_mid & (mid_side > 0) & (side < 0)))
        mid_side = side
        water_in_channel.append(int((lw == 2).sum()))
        o2_in_channel.append(int((lo2 == 2).sum()))
        o2_above.append(int((o[:, 2] > slab_hi).sum()))

    ns = a.production_ps / 1000.0
    print(json.dumps({
        "radius_A": g["radius"],
        "dp_mpa": a.dp_mpa,
        "density_scale": round(a.density_scale, 4),
        "mobile_oh": len(mobile),
        "crowded_oh_frozen": crowded,
        "force_kj_per_mol_nm": round(force_kj_nm, 4),
        "waters": len(water_o),
        "o2": len(o2_a),
        "production_ns": ns,
        "water_up": wt["up"],
        "water_down": wt["down"],
        "water_net_per_ns": round((wt["up"] - wt["down"]) / ns, 2),
        "o2_up": ot["up"],
        "o2_down": ot["down"],
        "o2_above_slab_final": o2_above[-1],
        "mean_water_in_channel": round(float(np.mean(water_in_channel)), 2),
        "fraction_o2_in_channel": round(float(np.mean(np.array(o2_in_channel) > 0)), 4),
        "mid_plane_up": mid_up,
        "mid_plane_down": mid_down,
        "platform": plat.getName(),
    }))
    if a.profile:
        n = int(a.production_ps / a.sample_ps)
        for lo_e, c in zip(edges[:-1], profile / n):
            print(f"z {10 * (lo_e - zc):+6.1f} A  axial waters {c:5.2f}")
        slab_vol = box_x * box_y * (full_edges[1] - full_edges[0])
        bulk = 1.0 / (WATER_VOLUME_A3 / 1000.0)
        for lo_e, c in zip(full_edges[:-1], full / n):
            print(f"box z {10 * (lo_e - zc):+6.1f} A  water density {c / slab_vol / bulk:4.2f} of bulk")
    return 0


def umbrella(a, ctx, integ, o2_first, slab_lo, slab_hi, zc, k_z, r_cyl, g) -> int:
    """Windows along z from the tank side to the far side; WHAM; leak rate."""
    kt = 0.0083145 * 310.0
    centres = np.arange(slab_lo - 0.7, slab_hi + 0.7 + 1e-9, a.window_step_nm)
    if a.reverse:
        centres = centres[::-1]
    series = []
    # Pull the O2 from its start to the first window.
    for z0 in centres:
        ctx.setParameter("z0", float(z0))
        integ.step(10000)  # 20 ps to settle in the window
        zs = []
        for _ in range(int(a.window_ps / 0.1)):
            integ.step(50)
            xyz = ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
            zs.append(0.5 * (xyz[o2_first][2] + xyz[o2_first + 1][2]))
        series.append(np.array(zs))
        print(json.dumps({"window": round(float(z0), 3), "mean_z": round(float(np.mean(zs)), 3),
                          "sd": round(float(np.std(zs)), 4)}), file=sys.stderr, flush=True)
    if a.save:
        np.savez(a.save, centres=centres, series=np.array(series, dtype=object), k_z=k_z,
                 slab_lo=slab_lo, slab_hi=slab_hi, zc=zc, r_cyl=r_cyl, thickness=g["thickness"],
                 radius=g["radius"])
    return report(a, list(centres), series, slab_lo, slab_hi, zc, k_z, r_cyl, g)


def report(a, centres, series, slab_lo, slab_hi, zc, k_z, r_cyl, g) -> int:
    kt = 0.0083145 * 310.0
    centres = np.array(centres)
    # WHAM on a grid.
    edges = np.arange(centres.min() - 0.3, centres.max() + 0.3, 0.02)
    mids = 0.5 * (edges[1:] + edges[:-1])
    hist = np.array([np.histogram(zs, edges)[0] for zs in series], float)
    n = hist.sum(axis=1)
    bias = np.array([0.5 * k_z * (mids - c) ** 2 for c in centres])
    f = np.zeros(len(centres))
    for _ in range(20000):
        denom = (n[:, None] * np.exp(-(bias - f[:, None]) / kt)).sum(axis=0)
        prob = hist.sum(axis=0) / np.where(denom > 0, denom, np.inf)
        f_new = -kt * np.log((prob[None, :] * np.exp(-bias / kt)).sum(axis=1))
        f_new -= f_new[0]
        if np.max(np.abs(f_new - f)) < 1e-7:
            f = f_new
            break
        f = f_new
    ok = prob > 0
    w = np.full(len(mids), np.nan)
    w[ok] = -kt * np.log(prob[ok])
    bulk = (mids < slab_lo - 0.3) | (mids > slab_hi + 0.3)
    w -= np.nanmean(w[bulk & ok])
    inside = (np.abs(mids - zc) < g["thickness"] / 20.0) & ok
    barrier = float(np.nanmax(w[ok & (mids > slab_lo - 0.3) & (mids < slab_hi + 0.3)]))
    # Solubility-diffusion estimate of the idle leak through one channel, from
    # c_bulk = 1 M on the tank side. D is an assumed 1e-9 m^2/s everywhere.
    d_m2_s = 1.0e-9
    span = ok & (mids > slab_lo - 0.3) & (mids < slab_hi + 0.3)
    resist_s_per_m = np.sum(np.exp(w[span] / kt)) * 0.02e-9 / d_m2_s
    area_m2 = math.pi * (r_cyl * 1e-9) ** 2
    c_per_m3 = 1.0e3 * AVOGADRO
    leak_per_s = c_per_m3 * area_m2 / resist_s_per_m
    print(json.dumps({
        "radius_A": g["radius"],
        "windows": len(centres),
        "window_ps": a.window_ps,
        "bulk_end_mismatch_kj": round(float(np.nanmean(w[(mids < slab_lo - 0.3) & ok]) - np.nanmean(w[(mids > slab_hi + 0.3) & ok])), 2),
        "barrier_kj": round(barrier, 1),
        "channel_mean_kj": round(float(np.nanmean(w[inside])), 1) if inside.any() else None,
        "empty_bins_in_span": int((~ok & (mids > slab_lo - 0.3) & (mids < slab_hi + 0.3)).sum()),
        "idle_leak_per_s_at_1M": float(f"{leak_per_s:.3e}"),
        "pmf": [[round(float(m), 3), None if not np.isfinite(x) else round(float(x), 2)] for m, x in zip(mids, w)],
    }))
    return 0


if __name__ == "__main__":
    sys.exit(main())
