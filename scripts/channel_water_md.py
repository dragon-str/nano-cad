"""Drain-channel kill test, part 2 (M28-09): does a water-filled lined
channel exclude dissolved oxygen?

Water that drains out of the tank carries dissolved oxygen only in
proportion to the oxygen that dissolves in the water inside the channel. In
a single-file channel oxygen cannot overtake water, so the loss per drained
water is the oxygen/water ratio inside the channel. The test measures the
excess chemical potential of oxygen in the channel and in bulk water by
Widom test-particle insertion into MD snapshots. The difference is the
exclusion free energy: dG = mu_ex(channel) - mu_ex(bulk). The loss falls by
exp(-dG / RT) against plain water.

Geometry from `channel_json`: a frozen diamond slab with one through-channel
lined with a functional group (with its model charges). SPC/E water fills
both sides; the box is periodic, so the two reservoirs join. 310 K.

Oxygen model: two neutral Lennard-Jones sites, 1.21 A apart, TraPPE-like
(sigma 3.02 A, epsilon 0.407 kJ/mol). The bulk value checks the method: the
experimental hydration free energy of O2 is about +8.4 kJ/mol.

Memory: one process, below 1 GB.

    python scripts/channel_water_md.py channel.json --production-ps 1000
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
UFF = {"C": (3.851, 0.105), "H": (2.886, 0.044), "O": (3.500, 0.060), "N": (3.660, 0.069)}
SPCE_SIGMA = 0.3166
SPCE_EPS = 0.6502
SPCE_Q_O = -0.8476
SPCE_Q_H = 0.4238
OH = 0.1
HH = 0.1633
WATER_VOLUME_A3 = 29.9
O2_SIGMA = 0.302
O2_EPS = 0.407
O2_BOND = 0.121
KT = 0.0083145 * 310.0


def uff(element):
    x1, d1 = UFF[element]
    return x1 / LJ_FACTOR / 10.0, d1 * KCAL


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("geometry")
    p.add_argument("--equilibrate-ps", type=float, default=50.0)
    p.add_argument("--production-ps", type=float, default=1000.0)
    p.add_argument("--sample-ps", type=float, default=2.0)
    p.add_argument("--insertions", type=int, default=4000)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--platform", default="OpenCL")
    a = p.parse_args()

    g = json.load(open(a.geometry))
    atoms = g["atoms"]
    solid = np.array([x[1:4] for x in atoms]) / 10.0
    elements = [x[0] for x in atoms]
    charges = np.array([x[4] for x in atoms])
    charged = np.abs(charges) > 1e-9
    if charged.any():
        charges[charged] -= charges.sum() / charged.sum()  # neutral slab
    half_t = g["thickness"] / 20.0
    lo, hi = solid.min(axis=0), solid.max(axis=0)
    # The slab spans the box in x and y (periodic), with a small margin.
    box_x = hi[0] - lo[0] + 0.25
    box_y = hi[1] - lo[1] + 0.25
    reservoir = 1.6
    box_z = (hi[2] - lo[2]) + 2 * reservoir
    shift = np.array([-lo[0] + 0.125, -lo[1] + 0.125, reservoir - lo[2]])
    solid = solid + shift
    axis = np.array([shift[0], shift[1]])  # channel axis at x = y = 0
    zc = 0.5 * (lo[2] + hi[2]) + shift[2]

    tree = cKDTree(solid)
    spacing = (WATER_VOLUME_A3 ** (1 / 3)) / 10.0
    rng = np.random.default_rng(a.seed)
    waters = []
    for ix in range(int(box_x / spacing)):
        for iy in range(int(box_y / spacing)):
            for iz in range(int(box_z / spacing)):
                o = np.array([(ix + 0.5) * spacing, (iy + 0.5) * spacing, (iz + 0.5) * spacing])
                if tree.query(o)[0] < 0.26:
                    continue
                waters.append(o)
    pos = list(solid)
    for o in waters:
        q = rng.normal(size=3)
        q /= np.linalg.norm(q)
        r = rng.normal(size=3)
        r -= np.dot(r, q) * q
        r /= np.linalg.norm(r)
        h = math.radians(109.47 / 2)
        pos += [o, o + OH * (math.cos(h) * q + math.sin(h) * r), o + OH * (math.cos(h) * q - math.sin(h) * r)]
    pos = np.array(pos)

    system = mm.System()
    system.setDefaultPeriodicBoxVectors(mm.Vec3(box_x, 0, 0), mm.Vec3(0, box_y, 0), mm.Vec3(0, 0, box_z))
    nb = mm.NonbondedForce()
    nb.setNonbondedMethod(mm.NonbondedForce.PME)
    nb.setCutoffDistance(0.9)
    solid_sig_eps = []
    for e, q in zip(elements, charges):
        system.addParticle(0.0)
        s, ep = uff(e)
        solid_sig_eps.append((s, ep))
        nb.addParticle(float(q), s, ep)
    first = len(elements)
    oxygens = []
    for _ in waters:
        o = system.addParticle(15.999)
        h1 = system.addParticle(1.008)
        h2 = system.addParticle(1.008)
        nb.addParticle(SPCE_Q_O, SPCE_SIGMA, SPCE_EPS)
        nb.addParticle(SPCE_Q_H, 0.1, 0.0)
        nb.addParticle(SPCE_Q_H, 0.1, 0.0)
        for i, j in ((o, h1), (o, h2), (h1, h2)):
            nb.addException(i, j, 0.0, 0.1, 0.0)
        system.addConstraint(o, h1, OH)
        system.addConstraint(o, h2, OH)
        system.addConstraint(h1, h2, HH)
        oxygens.append(o)
    for i, j in tree.query_pairs(0.9):
        nb.addException(int(i), int(j), 0.0, 0.1, 0.0)
    system.addForce(nb)
    system.addForce(mm.MonteCarloAnisotropicBarostat(mm.Vec3(1, 1, 1) * unit.bar, 310 * unit.kelvin, False, False, True))

    integ = mm.LangevinMiddleIntegrator(310 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
    integ.setRandomNumberSeed(a.seed)
    try:
        plat = mm.Platform.getPlatformByName(a.platform)
    except Exception:  # noqa: BLE001
        plat = mm.Platform.getPlatformByName("CPU")
    ctx = mm.Context(system, integ, plat)
    ctx.setPositions(pos * unit.nanometer)
    mm.LocalEnergyMinimizer.minimize(ctx, 10.0, 1000)
    ctx.setVelocitiesToTemperature(310 * unit.kelvin, a.seed)
    integ.step(int(a.equilibrate_ps / 0.002))

    # Widom: O2 sites against water O (SPC/E LJ) and solid atoms (UFF), LB
    # mixing, 0.9 nm cutoff. Neutral O2, so no electrostatics.
    solid_sig = np.array([s for s, _ in solid_sig_eps])
    solid_eps = np.array([e for _, e in solid_sig_eps])

    def lj(r, sig, eps):
        x = (sig / r) ** 6
        return 4 * eps * (x * x - x)

    def insertion_energy(centres, dirs, ow, ow_tree, box):
        out = np.empty(len(centres))
        for k, (c, d) in enumerate(zip(centres, dirs)):
            total = 0.0
            for s in (-0.5, 0.5):
                site = c + s * O2_BOND * d
                site_w = site % box
                idx = ow_tree.query_ball_point(site_w, 0.9)
                if idx:
                    dv = ow[idx] - site_w
                    dv -= box * np.round(dv / box)
                    r = np.linalg.norm(dv, axis=1)
                    if r.min() < 0.2:
                        total = np.inf
                        break
                    total += lj(r, 0.5 * (O2_SIGMA + SPCE_SIGMA), math.sqrt(O2_EPS * SPCE_EPS)).sum()
                jdx = tree.query_ball_point(site, 0.9)
                if jdx:
                    r = np.linalg.norm(solid[jdx] - site, axis=1)
                    if r.min() < 0.2:
                        total = np.inf
                        break
                    total += lj(r, 0.5 * (O2_SIGMA + solid_sig[jdx]), np.sqrt(O2_EPS * solid_eps[jdx])).sum()
            out[k] = total
        return out

    chan_r = 0.12  # nm: the region an O2 centre can occupy in the channel
    bulk_sum = chan_sum = 0.0
    n_bulk = n_chan = 0
    water_in_channel = []
    samples = int(a.production_ps / a.sample_ps)
    for _ in range(samples):
        integ.step(int(a.sample_ps / 0.002))
        st = ctx.getState(getPositions=True, enforcePeriodicBox=False)
        box = np.array(st.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer).diagonal())
        allpos = st.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        ow = allpos[first::3].copy()
        ow %= box
        ow_tree = cKDTree(ow, boxsize=box)
        dx = (ow[:, 0] - axis[0] + box[0] / 2) % box[0] - box[0] / 2
        dy = (ow[:, 1] - axis[1] + box[1] / 2) % box[1] - box[1] / 2
        inside = (np.hypot(dx, dy) < 0.3) & (np.abs(ow[:, 2] - zc) < half_t)
        water_in_channel.append(int(inside.sum()))
        n = a.insertions
        dirs = rng.normal(size=(n, 3))
        dirs /= np.linalg.norm(dirs, axis=1)[:, None]
        # Bulk: random points in the reservoirs, far from the slab.
        zb = np.where(rng.random(n) < 0.5, zc - half_t - 0.5 - rng.random(n) * 0.5, zc + half_t + 0.5 + rng.random(n) * 0.5)
        cb = np.column_stack([rng.random(n) * box[0], rng.random(n) * box[1], zb])
        ub = insertion_energy(cb, dirs, ow, ow_tree, box)
        bulk_sum += np.exp(-ub / KT).sum()
        n_bulk += n
        # Channel: random points in a cylinder of radius chan_r about the axis.
        rr = chan_r * np.sqrt(rng.random(n))
        th = rng.random(n) * 2 * math.pi
        cc = np.column_stack([axis[0] + rr * np.cos(th), axis[1] + rr * np.sin(th), zc + (rng.random(n) - 0.5) * 2 * (half_t - 0.1)])
        uc = insertion_energy(cc, dirs, ow, ow_tree, box)
        chan_sum += np.exp(-uc / KT).sum()
        n_chan += n

    mu_bulk = -KT * math.log(bulk_sum / n_bulk) if bulk_sum > 0 else math.inf
    mu_chan = -KT * math.log(chan_sum / n_chan) if chan_sum > 0 else math.inf
    wc = np.array(water_in_channel)
    print(json.dumps({
        "radius_A": g["radius"],
        "production_ps": a.production_ps,
        "insertions_per_region": n_bulk,
        "mean_waters_in_channel": round(float(wc.mean()), 2),
        "fraction_channel_empty": round(float((wc == 0).mean()), 3),
        "mu_ex_bulk_kj": round(mu_bulk, 2),
        "mu_ex_channel_kj": round(mu_chan, 2) if math.isfinite(mu_chan) else None,
        "exclusion_dG_kj": round(mu_chan - mu_bulk, 2) if math.isfinite(mu_chan) else None,
        "channel_hits": int(chan_sum > 0),
        "platform": plat.getName(),
    }))
    return 0


if __name__ == "__main__":
    sys.exit(main())
