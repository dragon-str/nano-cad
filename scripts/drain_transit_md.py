"""Drain transit test (M28-10): count water and oxygen that cross a lined
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
    a = p.parse_args()

    g = json.load(open(a.geometry))
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
    spacing = (WATER_VOLUME_A3 ** (1 / 3)) / 10.0
    sites = []
    for ix in range(int(box_x / spacing)):
        for iy in range(int(box_y / spacing)):
            for iz in range(int(box_z / spacing)):
                o = np.array([(ix + 0.5) * spacing, (iy + 0.5) * spacing, (iz + 0.5) * spacing])
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

    system = mm.System()
    system.setDefaultPeriodicBoxVectors(mm.Vec3(box_x, 0, 0), mm.Vec3(0, box_y, 0), mm.Vec3(0, 0, box_z))
    nb = mm.NonbondedForce()
    nb.setNonbondedMethod(mm.NonbondedForce.PME)
    nb.setCutoffDistance(0.9)
    pos = list(solid)
    for e, q in zip(elements, charges):
        system.addParticle(0.0)
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
        nb.addException(int(i), int(j), 0.0, 0.1, 0.0)
    system.addForce(nb)

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

    integ = mm.LangevinMiddleIntegrator(310 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
    integ.setRandomNumberSeed(a.seed)
    try:
        plat = mm.Platform.getPlatformByName(a.platform)
    except Exception:  # noqa: BLE001
        plat = mm.Platform.getPlatformByName("CPU")
    ctx = mm.Context(system, integ, plat)
    ctx.setPositions(np.array(pos) * unit.nanometer)
    mm.LocalEnergyMinimizer.minimize(ctx, 10.0, 1000)
    ctx.setVelocitiesToTemperature(310 * unit.kelvin, a.seed)
    integ.step(int(a.equilibrate_ps / 0.002))

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
    for _ in range(int(a.production_ps / a.sample_ps)):
        integ.step(int(a.sample_ps / 0.002))
        xyz = ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        w = xyz[water_o]
        o = 0.5 * (xyz[o2_a] + xyz[np.array(o2_a) + 1])
        lw, lo2 = labels(w[:, 2], w[:, 0], w[:, 1]), labels(o[:, 2], o[:, 0], o[:, 1])
        update(wt, lw)
        update(ot, lo2)
        water_in_channel.append(int((lw == 2).sum()))
        o2_in_channel.append(int((lo2 == 2).sum()))
        o2_above.append(int((o[:, 2] > slab_hi).sum()))

    ns = a.production_ps / 1000.0
    print(json.dumps({
        "radius_A": g["radius"],
        "dp_mpa": a.dp_mpa,
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
        "platform": plat.getName(),
    }))
    return 0


def umbrella(a, ctx, integ, o2_first, slab_lo, slab_hi, zc, k_z, r_cyl, g) -> int:
    """Windows along z from the tank side to the far side; WHAM; leak rate."""
    kt = 0.0083145 * 310.0
    centres = np.arange(slab_lo - 0.7, slab_hi + 0.7 + 1e-9, a.window_step_nm)
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
    # WHAM on a grid.
    edges = np.arange(centres[0] - 0.3, centres[-1] + 0.3, 0.02)
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
        "barrier_kj": round(barrier, 1),
        "channel_mean_kj": round(float(np.nanmean(w[inside])), 1) if inside.any() else None,
        "empty_bins_in_span": int((~ok & (mids > slab_lo - 0.3) & (mids < slab_hi + 0.3)).sum()),
        "idle_leak_per_s_at_1M": float(f"{leak_per_s:.3e}"),
        "pmf": [[round(float(m), 3), None if not np.isfinite(x) else round(float(x), 2)] for m, x in zip(mids, w)],
    }))
    return 0


if __name__ == "__main__":
    sys.exit(main())
