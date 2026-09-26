"""Gas-space test of option A (M28-05): does water flood a diamond slot?

The collection chamber of the two-stage pump needs a gas space that the
liquid never floods. The gas space lives in a slot between hydrogen-capped
diamond plates. If water wets the capped diamond (contact angle below 90
degrees), capillary action fills the slot and the gas space is lost. If the
surface repels water, the slot stays dry.

Geometry from `slot_json`: two frozen plates with a gap between their facing
carbon planes. SPC/E water surrounds the plates, with a free liquid surface
above and vacuum over it, so the liquid sits near zero pressure. The slot
starts dry (`--start dry`, water removed) or filled (`--start wet`). The run
counts the water molecules inside the slot.

Model as in `pocket_water_md.py`: UFF Lennard-Jones on the solid,
Lorentz-Berthelot mixing, no charges on the solid. A simulated estimate.
Memory: one process, well under 1 GB. The GPU does the work.

    python scripts/slot_water_md.py slot.json --start dry --production-ps 1000
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
UFF = {"C": (3.851, 0.105), "H": (2.886, 0.044)}
SPCE_SIGMA_NM = 0.3166
SPCE_EPSILON_KJ = 0.6502
SPCE_Q_O = -0.8476
SPCE_Q_H = 0.4238
SPCE_OH_NM = 0.1
SPCE_HH_NM = 0.1633
WATER_VOLUME_A3 = 29.9


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("geometry")
    parser.add_argument("--start", choices=["wet", "dry"], default="dry")
    parser.add_argument("--equilibrate-ps", type=float, default=20.0)
    parser.add_argument("--production-ps", type=float, default=1000.0)
    parser.add_argument("--sample-ps", type=float, default=2.0)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--platform", default="OpenCL")
    args = parser.parse_args()

    geometry = json.load(open(args.geometry))
    gap = geometry["gap"] / 10.0
    atoms = geometry["atoms"]
    solid = np.array([a[1:4] for a in atoms]) / 10.0
    elements = [a[0] for a in atoms]
    lo, hi = solid.min(axis=0), solid.max(axis=0)

    margin_xy = 1.0
    water_above = 1.5
    vacuum = 3.0
    box_x = hi[0] - lo[0] + 2 * margin_xy
    box_y = hi[1] - lo[1] + 2 * margin_xy
    floor = 0.2
    shift = np.array([margin_xy - lo[0], margin_xy - lo[1], floor + water_above - lo[2]])
    solid_s = solid + shift
    top_liquid = hi[2] + shift[2] + water_above
    box_z = top_liquid + vacuum
    center = 0.5 * (lo + hi) + shift
    half_x = 0.5 * (hi[0] - lo[0]) - 0.15
    half_y = 0.5 * (hi[1] - lo[1]) - 0.15
    slot_half = 0.5 * gap - 0.1

    def in_slot(p):
        """p is an (n, 3) array; True where the point lies in the slot."""
        return (
            (np.abs(p[:, 0] - center[0]) < half_x)
            & (np.abs(p[:, 1] - center[1]) < half_y)
            & (np.abs(p[:, 2] - center[2]) < slot_half)
        )

    spacing = (WATER_VOLUME_A3 ** (1 / 3)) / 10.0
    tree = cKDTree(solid_s)
    waters = []
    for ix in range(int(box_x / spacing)):
        for iy in range(int(box_y / spacing)):
            for iz in range(int((top_liquid - floor) / spacing)):
                o = np.array([(ix + 0.5) * spacing, (iy + 0.5) * spacing, floor + (iz + 0.5) * spacing])
                if tree.query(o)[0] < 0.28:
                    continue
                if args.start == "dry" and in_slot(o[None, :])[0]:
                    continue
                waters.append(o)
    rng = np.random.default_rng(args.seed)
    positions = list(solid_s)
    for o in waters:
        q = rng.normal(size=3)
        q /= np.linalg.norm(q)
        r = rng.normal(size=3)
        r -= np.dot(r, q) * q
        r /= np.linalg.norm(r)
        half = math.radians(109.47 / 2)
        positions += [
            o,
            o + SPCE_OH_NM * (math.cos(half) * q + math.sin(half) * r),
            o + SPCE_OH_NM * (math.cos(half) * q - math.sin(half) * r),
        ]
    positions = np.array(positions)

    system = mm.System()
    system.setDefaultPeriodicBoxVectors(mm.Vec3(box_x, 0, 0), mm.Vec3(0, box_y, 0), mm.Vec3(0, 0, box_z))
    nb = mm.NonbondedForce()
    nb.setNonbondedMethod(mm.NonbondedForce.PME)
    nb.setCutoffDistance(1.0)
    for element in elements:
        system.addParticle(0.0)
        x1, d1 = UFF[element]
        nb.addParticle(0.0, x1 / LJ_FACTOR / 10.0, d1 * KCAL)
    first_water = len(elements)
    oxygens = []
    for _ in waters:
        o = system.addParticle(15.999)
        h1 = system.addParticle(1.008)
        h2 = system.addParticle(1.008)
        nb.addParticle(SPCE_Q_O, SPCE_SIGMA_NM, SPCE_EPSILON_KJ)
        nb.addParticle(SPCE_Q_H, 0.1, 0.0)
        nb.addParticle(SPCE_Q_H, 0.1, 0.0)
        for a, b in ((o, h1), (o, h2), (h1, h2)):
            nb.addException(a, b, 0.0, 0.1, 0.0)
        system.addConstraint(o, h1, SPCE_OH_NM)
        system.addConstraint(o, h2, SPCE_OH_NM)
        system.addConstraint(h1, h2, SPCE_HH_NM)
        oxygens.append(o)
    for i, j in cKDTree(solid_s).query_pairs(1.0):
        nb.addException(int(i), int(j), 0.0, 0.1, 0.0)
    system.addForce(nb)
    wall = mm.CustomExternalForce("k*step(zlo-z)*(zlo-z)^2")
    wall.addGlobalParameter("k", 5000.0)
    wall.addGlobalParameter("zlo", floor)
    for o in oxygens:
        wall.addParticle(o, [])
    system.addForce(wall)

    integrator = mm.LangevinMiddleIntegrator(310 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
    integrator.setRandomNumberSeed(args.seed)
    try:
        platform = mm.Platform.getPlatformByName(args.platform)
    except Exception:  # noqa: BLE001
        platform = mm.Platform.getPlatformByName("CPU")
    context = mm.Context(system, integrator, platform)
    context.setPositions(positions * unit.nanometer)
    mm.LocalEnergyMinimizer.minimize(context, 10.0, 500)
    context.setVelocitiesToTemperature(310 * unit.kelvin, args.seed)
    integrator.step(int(args.equilibrate_ps / 0.002))

    counts = []
    for _ in range(int(args.production_ps / args.sample_ps)):
        integrator.step(int(args.sample_ps / 0.002))
        pos = context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        ox = pos[first_water::3].copy()
        ox[:, 0] = (ox[:, 0] - center[0] + box_x / 2) % box_x - box_x / 2 + center[0]
        ox[:, 1] = (ox[:, 1] - center[1] + box_y / 2) % box_y - box_y / 2 + center[1]
        counts.append(int(in_slot(ox).sum()))
    counts = np.array(counts)
    slot_volume = (2 * half_x) * (2 * half_y) * (2 * slot_half)
    full = slot_volume / (WATER_VOLUME_A3 / 1000.0)
    tail = counts[len(counts) // 2 :]
    print(
        json.dumps(
            {
                "gap_nm": gap,
                "start": args.start,
                "waters": len(waters),
                "slot_full_count": round(full, 1),
                "first_50ps_mean": round(float(counts[:25].mean()), 1),
                "last_half_mean": round(float(tail.mean()), 1),
                "last_half_fill": round(float(tail.mean() / full), 3),
                "min": int(counts.min()),
                "max": int(counts.max()),
                "platform": platform.getName(),
            }
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
