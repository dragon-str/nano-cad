"""M28-24: how fast does dissolved O2 enter the snug pocket?

The pocket geometry comes from `rod_block_json` (rod pulled back, the
loading state). The block and the rod are frozen. SPC/E water fills the space
above; mobile O2 molecules (two LJ sites, 1.21 A, sigma 3.02 A, epsilon
0.407 kJ/mol) are dissolved in it at `--o2-molar`. The iron keeps its Morse
well for water (the proxy, ADR-0086), because O2 must displace a coordinated
water. O2 has no Morse well to the iron: the run measures diffusive arrival
at the iron, not bond formation.

Counted events:
- entry: an O2 centre crosses from outside into the pocket (inside the bore
  radius, between the rod face and the mouth);
- dock: an O2 atom comes within `--dock-nm` of the iron, from a state where
  no O2 atom was within that distance.

k_entry = entries per second / [O2]. The physiological load rate is
k_dock * c with c = 130 uM (arterial pO2 100 mmHg).

One OpenMM process, below 0.5 GB. A simulated estimate.
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
UFF = {"C": (3.851, 0.105), "H": (2.886, 0.044), "Fe": (2.912, 0.013), "O": (3.500, 0.060)}
SPCE_SIGMA, SPCE_EPS, SPCE_Q_O, SPCE_Q_H = 0.3166, 0.6502, -0.8476, 0.4238
OH, HH = 0.1, 0.1633
WATER_VOLUME_A3 = 29.9
O2_SIGMA, O2_EPS, O2_BOND = 0.302, 0.407, 0.121
MORSE_D_KJ, MORSE_R0_NM, MORSE_A_PER_NM = 10.3, 0.270, 27.437
AVOGADRO = 6.02214076e23


def uff(element):
    x1, d1 = UFF[element]
    return x1 / LJ_FACTOR / 10.0, d1 * KCAL


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("geometry")
    p.add_argument("--o2-molar", type=float, default=0.5)
    p.add_argument("--dock-nm", type=float, default=0.40)
    p.add_argument("--equilibrate-ps", type=float, default=20.0)
    p.add_argument("--production-ps", type=float, default=2000.0)
    p.add_argument("--sample-ps", type=float, default=0.2)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--platform", default="OpenCL")
    a = p.parse_args()

    g = json.load(open(a.geometry))
    atoms = g["atoms"]
    iron = g["iron"]
    bore = g["bore_radius"] / 10.0
    solid = np.array([x[1:4] for x in atoms]) / 10.0
    elements = [x[0] for x in atoms]
    block = solid[: g["block_atoms"]]
    bmin, bmax = block.min(axis=0), block.max(axis=0)
    margin = 0.3
    box_x = bmax[0] - bmin[0] + margin
    box_y = bmax[1] - bmin[1] + margin
    water_top = bmax[2] + 2.5
    z_shift = -(bmin[2] - 1.0)
    box_z = water_top + z_shift + 1.5
    shift = np.array([-bmin[0] + margin / 2, -bmin[1] + margin / 2, z_shift])
    solid = solid + shift
    ax, ay = shift[0], shift[1]
    fe = solid[iron]
    lo = g["face_z"] / 10 + z_shift
    hi = g["mouth_z"] / 10 + z_shift
    floor_z = bmin[2] + z_shift

    rng = np.random.default_rng(a.seed)
    tree = cKDTree(solid)
    spacing = (WATER_VOLUME_A3 ** (1 / 3)) / 10.0
    sites = []
    for ix in range(int(box_x / spacing)):
        for iy in range(int(box_y / spacing)):
            for iz in range(int((water_top + z_shift - floor_z) / spacing)):
                o = np.array([(ix + 0.5) * spacing, (iy + 0.5) * spacing, floor_z + 0.05 + (iz + 0.5) * spacing])
                if tree.query(o)[0] >= 0.26:
                    sites.append(o)
    sites = np.array(sites)
    # O2 starts in bulk water at least 0.6 nm above the block top.
    bulk = np.where(sites[:, 2] > bmax[2] + z_shift + 0.6)[0]
    bulk_volume_l = box_x * box_y * (water_top - bmax[2] - 0.6) * 1e-24
    n_o2 = max(1, int(round(a.o2_molar * bulk_volume_l * AVOGADRO)))
    o2_set = set(rng.choice(bulk, size=min(n_o2, len(bulk)), replace=False).tolist())
    total_water_volume_l = len(sites) * WATER_VOLUME_A3 * 1e-27
    molar = n_o2 / AVOGADRO / total_water_volume_l

    system = mm.System()
    system.setDefaultPeriodicBoxVectors(mm.Vec3(box_x, 0, 0), mm.Vec3(0, box_y, 0), mm.Vec3(0, 0, box_z))
    nb = mm.NonbondedForce()
    nb.setNonbondedMethod(mm.NonbondedForce.PME)
    nb.setCutoffDistance(1.0)
    pos = list(solid)
    for e in elements:
        system.addParticle(0.0)
        s, ep = uff(e)
        nb.addParticle(0.0, s, ep)
    water_o, o2_a = [], []
    for k, o in enumerate(sites):
        d = rng.normal(size=3)
        d /= np.linalg.norm(d)
        if k in o2_set:
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
        nb.addException(iron, i, 0.0, 0.1, 0.0)
    for i, j in tree.query_pairs(1.0):
        nb.addException(int(i), int(j), 0.0, 0.1, 0.0)
    system.addForce(nb)
    morse = mm.CustomBondForce("D*((1-exp(-a*(r-r0)))^2-1)")
    morse.addGlobalParameter("D", MORSE_D_KJ)
    morse.addGlobalParameter("a", MORSE_A_PER_NM)
    morse.addGlobalParameter("r0", MORSE_R0_NM)
    morse.setUsesPeriodicBoundaryConditions(True)
    for i in water_o:
        morse.addBond(iron, i, [])
    system.addForce(morse)
    wall = mm.CustomExternalForce("k*(step(zl-z)*(zl-z)^2+step(z-zh)*(z-zh)^2)")
    wall.addGlobalParameter("k", 5000.0)
    wall.addGlobalParameter("zl", float(floor_z))
    wall.addGlobalParameter("zh", float(box_z - 0.3))
    for i in water_o:
        wall.addParticle(i, [])
    system.addForce(wall)
    # O2 is hydrophobic and collects at the free water surface at the top of
    # the box; seeds that did so gave almost no entries (2026-09-29). Keep
    # O2 at least 0.5 nm below that surface.
    o2_wall = mm.CustomExternalForce("k2*(step(zl2-z)*(zl2-z)^2+step(z-zh2)*(z-zh2)^2)")
    o2_wall.addGlobalParameter("k2", 5000.0)
    o2_wall.addGlobalParameter("zl2", float(floor_z))
    o2_wall.addGlobalParameter("zh2", float(water_top + z_shift - 0.5))
    for i in o2_a:
        o2_wall.addParticle(i, [])
        o2_wall.addParticle(i + 1, [])
    system.addForce(o2_wall)

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

    o2_a = np.array(o2_a)
    was_in = np.zeros(len(o2_a), bool)
    was_docked = False
    entries = docks = 0
    occupied = docked_samples = 0
    residence = []
    heights = []
    t_in = np.zeros(len(o2_a))
    samples = int(a.production_ps / a.sample_ps)
    for step in range(samples):
        integ.step(int(a.sample_ps / 0.002))
        xyz = ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        o1, o2p = xyz[o2_a], xyz[o2_a + 1]
        c = 0.5 * (o1 + o2p)
        dx = (c[:, 0] - ax + box_x / 2) % box_x - box_x / 2
        dy = (c[:, 1] - ay + box_y / 2) % box_y - box_y / 2
        inside = (np.hypot(dx, dy) < bore) & (c[:, 2] > lo) & (c[:, 2] < hi)
        new = inside & ~was_in
        entries += int(new.sum())
        t_in[new] = step
        for k in np.where(was_in & ~inside)[0]:
            residence.append((step - t_in[k]) * a.sample_ps)
        was_in = inside
        occupied += int(inside.any())
        heights.append(float(np.mean(c[:, 2]) - (bmax[2] + z_shift)))

        def fe_dist(p):
            v = p - fe
            v[:, 0] -= box_x * np.round(v[:, 0] / box_x)
            v[:, 1] -= box_y * np.round(v[:, 1] / box_y)
            return np.linalg.norm(v, axis=1)

        docked = bool((np.minimum(fe_dist(o1), fe_dist(o2p)) < a.dock_nm).any())
        if docked and not was_docked:
            docks += 1
        was_docked = docked
        docked_samples += int(docked)

    seconds = a.production_ps * 1e-12
    print(json.dumps({
        "o2": len(o2_a),
        "o2_molar": round(molar, 3),
        "production_ps": a.production_ps,
        "entries": entries,
        "docks": docks,
        "fraction_pocket_occupied": round(occupied / samples, 4),
        "fraction_docked": round(docked_samples / samples, 4),
        "mean_residence_ps": round(float(np.mean(residence)), 2) if residence else None,
        "k_entry_per_M_s": float(f"{entries / seconds / molar:.3e}"),
        "k_dock_per_M_s": float(f"{docks / seconds / molar:.3e}"),
        "load_rate_at_130uM_per_s": float(f"{docks / seconds / molar * 1.3e-4:.3e}"),
        "mean_o2_height_above_block_nm": round(float(np.mean(heights)), 3),
        "platform": plat.getName(),
    }))
    return 0


if __name__ == "__main__":
    sys.exit(main())
