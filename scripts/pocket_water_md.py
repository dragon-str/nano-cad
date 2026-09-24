"""Explicit-water test of the rod pocket (ADR-0090, test 2).

The question: does the loaded pocket stay dry in liquid water? The pocket is
the bore above the pulled-back rod, with one iron atom in the wall. The block
and the rod are frozen. SPC/E water fills the space above the block. The run
counts the water molecules inside the pocket and at the iron.

Input is the geometry JSON that `rod_block_json` writes:

    cargo run --release -p nanocad-meter --example rod_block_json -- pocket.json
    python scripts/pocket_water_md.py pocket.json --start wet --o2 none

`--start wet` fills the pocket with water at the start; `--start dry` leaves
it empty. The two starts bracket a slow filling or drying. `--o2 bound` puts
a frozen oxygen molecule on the iron, in the bound pose.

The model: UFF Lennard-Jones on carbon, hydrogen, iron and oxygen of O2;
Lorentz-Berthelot mixing (UFF uses geometric mixing, so this is an
approximation); no charges on the solid; a Morse term between the iron and
each water oxygen (10.3 kJ/mol at 2.7 A, the Fe(NH3)4 proxy of ADR-0086) that
replaces their Lennard-Jones pair. The result is a simulated estimate.

Needs OpenMM (MIT/LGPL), which is not a dependency of the package. Install it
in a separate environment.
"""

from __future__ import annotations

import argparse
import json
import math
import sys

import numpy as np
import openmm as mm
from openmm import unit

KCAL = 4.184  # kJ per kcal
LJ_FACTOR = 2.0 ** (1.0 / 6.0)
# UFF x1 (A) and D1 (kcal/mol), Rappe et al. 1992.
UFF = {"C": (3.851, 0.105), "H": (2.886, 0.044), "Fe": (2.912, 0.013), "O": (3.500, 0.060)}
# SPC/E water.
SPCE_SIGMA_NM = 0.3166
SPCE_EPSILON_KJ = 0.6502
SPCE_Q_O = -0.8476
SPCE_Q_H = 0.4238
SPCE_OH_NM = 0.1
SPCE_HH_NM = 0.1633
# Iron-water coordination (Morse), from the proxy.
MORSE_D_KJ = 10.3
MORSE_R0_NM = 0.270
MORSE_A_PER_NM = 27.437
WATER_VOLUME_A3 = 29.9  # one molecule of liquid water


def uff_sigma_epsilon(element: str) -> tuple[float, float]:
    x1, d1 = UFF[element]
    return x1 / LJ_FACTOR / 10.0, d1 * KCAL


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("geometry")
    parser.add_argument("--start", choices=["wet", "dry"], default="wet")
    parser.add_argument("--o2", choices=["none", "bound"], default="none")
    parser.add_argument("--equilibrate-ps", type=float, default=50.0)
    parser.add_argument("--production-ps", type=float, default=500.0)
    parser.add_argument("--sample-ps", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--platform", default="OpenCL")
    args = parser.parse_args()

    geometry = json.load(open(args.geometry))
    atoms = geometry["atoms"]
    iron = geometry["iron"]
    bore_a = geometry["bore_radius"]
    face_z = geometry["face_z"]
    mouth_z = geometry["mouth_z"]
    solid = np.array([a[1:] for a in atoms]) / 10.0  # nm
    elements = [a[0] for a in atoms]
    block_count = geometry["block_atoms"]
    block = solid[:block_count]
    block_min = block.min(axis=0)
    block_max = block.max(axis=0)

    # The iron's open side: away from its three bonded carbons.
    fe = solid[iron]
    neighbours = [
        p for i, p in enumerate(solid) if i != iron and elements[i] == "C"
        and np.linalg.norm(p - fe) < 0.17
    ]
    axis = -sum((p - fe) / np.linalg.norm(p - fe) for p in neighbours)
    axis /= np.linalg.norm(axis)

    extra = []  # frozen O2 atoms
    if args.o2 == "bound":
        donor = fe + 0.187 * axis
        # Bent end-on: Fe-O-O near 120 degrees, the far atom toward the mouth.
        up = np.array([0.0, 0.0, 1.0])
        side = up - np.dot(up, axis) * axis
        side = side / np.linalg.norm(side) if np.linalg.norm(side) > 1e-6 else np.array([1.0, 0, 0])
        far = donor + 0.121 * (0.5 * axis + math.sqrt(3) / 2 * side)
        extra = [donor, far]

    # The box: the block width plus a margin in x and y, periodic; a tall z
    # with water above the block and vacuum over it.
    margin = 0.3
    box_x = block_max[0] - block_min[0] + margin
    box_y = block_max[1] - block_min[1] + margin
    water_top = block_max[2] + 2.5
    box_z = (water_top - (block_min[2] - 1.0)) + 1.5
    z_shift = -(block_min[2] - 1.0)  # move so everything sits at z >= 0
    x_shift = -block_min[0] + margin / 2
    y_shift = -block_min[1] + margin / 2
    shift = np.array([x_shift, y_shift, z_shift])
    solid_s = solid + shift
    extra_s = [p + shift for p in extra]
    fe_s = fe + shift
    axis_center = np.array([x_shift, y_shift])  # bore axis is x = y = 0

    # Water on a cubic grid at liquid density, outside the solid.
    spacing = (WATER_VOLUME_A3 ** (1 / 3)) / 10.0
    heavy = np.vstack([solid_s] + ([np.array(extra_s)] if extra_s else []))
    rng = np.random.default_rng(args.seed)
    waters = []
    from scipy.spatial import cKDTree  # noqa: PLC0415

    tree = cKDTree(heavy)
    floor_z = z_shift + block_min[2]
    z0 = floor_z + 0.05
    for ix in range(int(box_x / spacing)):
        for iy in range(int(box_y / spacing)):
            for iz in range(int((water_top + z_shift - z0) / spacing)):
                o = np.array([(ix + 0.5) * spacing, (iy + 0.5) * spacing, z0 + (iz + 0.5) * spacing])
                if tree.query(o)[0] < 0.26:
                    continue
                radial = np.hypot(o[0] - axis_center[0], o[1] - axis_center[1])
                in_pocket = radial < bore_a / 10 and face_z / 10 + z_shift < o[2] < mouth_z / 10 + z_shift
                if args.start == "dry" and in_pocket:
                    continue
                # Only the region above the block top, the pocket, and the
                # gaps at the block sides.
                waters.append(o)
    # Random orientation of each water.
    positions = [p for p in solid_s] + list(extra_s)
    for o in waters:
        q = rng.normal(size=3)
        q /= np.linalg.norm(q)
        r = rng.normal(size=3)
        r -= np.dot(r, q) * q
        r /= np.linalg.norm(r)
        half = math.radians(109.47 / 2)
        h1 = o + SPCE_OH_NM * (math.cos(half) * q + math.sin(half) * r)
        h2 = o + SPCE_OH_NM * (math.cos(half) * q - math.sin(half) * r)
        positions += [o, h1, h2]
    positions = np.array(positions)

    system = mm.System()
    system.setDefaultPeriodicBoxVectors(
        mm.Vec3(box_x, 0, 0), mm.Vec3(0, box_y, 0), mm.Vec3(0, 0, box_z)
    )
    nb = mm.NonbondedForce()
    nb.setNonbondedMethod(mm.NonbondedForce.PME)
    nb.setCutoffDistance(1.0)
    nb.setUseDispersionCorrection(False)
    for element in elements:
        system.addParticle(0.0)  # frozen
        sigma, epsilon = uff_sigma_epsilon(element)
        nb.addParticle(0.0, sigma, epsilon)
    for _ in extra_s:
        system.addParticle(0.0)
        sigma, epsilon = uff_sigma_epsilon("O")
        nb.addParticle(0.0, sigma, epsilon)
    first_water = len(elements) + len(extra_s)
    water_oxygens = []
    for w in range(len(waters)):
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
        water_oxygens.append(o)
        # The Morse term replaces the Fe-Ow Lennard-Jones pair.
        nb.addException(iron, o, 0.0, 0.1, 0.0)
    # Frozen-frozen pairs do not move anything; exclude them to save time.
    frozen = len(elements) + len(extra_s)
    tree_all = cKDTree(positions[:frozen])
    for i, j in tree_all.query_pairs(1.0):
        nb.addException(int(i), int(j), 0.0, 0.1, 0.0)
    system.addForce(nb)

    morse = mm.CustomBondForce("D*((1-exp(-a*(r-r0)))^2-1)")
    morse.addGlobalParameter("D", MORSE_D_KJ)
    morse.addGlobalParameter("a", MORSE_A_PER_NM)
    morse.addGlobalParameter("r0", MORSE_R0_NM)
    morse.setUsesPeriodicBoundaryConditions(True)
    for o in water_oxygens:
        morse.addBond(iron, o, [])
    system.addForce(morse)

    # Keep water above the block floor and below the box top.
    wall = mm.CustomExternalForce("k*(step(zlo-z)*(zlo-z)^2 + step(z-zhi)*(z-zhi)^2)")
    wall.addGlobalParameter("k", 5000.0)
    wall.addGlobalParameter("zlo", float(floor_z))
    wall.addGlobalParameter("zhi", float(box_z - 0.3))
    for o in water_oxygens:
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

    steps_per_sample = int(args.sample_ps / 0.002)
    integrator.step(int(args.equilibrate_ps / 0.002))

    pocket_lo = face_z / 10 + z_shift
    pocket_hi = mouth_z / 10 + z_shift
    samples = []
    iron_contacts = []
    for _ in range(int(args.production_ps / args.sample_ps)):
        integrator.step(steps_per_sample)
        state = context.getState(getPositions=True)
        pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        ox = pos[first_water::3]
        # Wrap x and y into the box around the bore axis.
        dx = (ox[:, 0] - axis_center[0] + box_x / 2) % box_x - box_x / 2
        dy = (ox[:, 1] - axis_center[1] + box_y / 2) % box_y - box_y / 2
        radial = np.hypot(dx, dy)
        inside = (radial < bore_a / 10) & (ox[:, 2] > pocket_lo) & (ox[:, 2] < pocket_hi)
        samples.append(int(inside.sum()))
        d_fe = np.linalg.norm(ox - fe_s, axis=1)
        iron_contacts.append(int((d_fe < 0.32).sum()))

    samples = np.array(samples)
    contacts = np.array(iron_contacts)
    pocket_volume = math.pi * (bore_a / 10) ** 2 * (pocket_hi - pocket_lo)
    result = {
        "start": args.start,
        "o2": args.o2,
        "waters": len(waters),
        "production_ps": args.production_ps,
        "pocket_cylinder_nm3": round(pocket_volume, 4),
        "mean_waters_in_pocket": round(float(samples.mean()), 3),
        "max_waters_in_pocket": int(samples.max()),
        "fraction_empty": round(float((samples == 0).mean()), 3),
        "mean_iron_contacts": round(float(contacts.mean()), 3),
        "fraction_iron_wet": round(float((contacts > 0).mean()), 3),
        "platform": platform.getName(),
    }
    print(json.dumps(result))
    return 0


if __name__ == "__main__":
    sys.exit(main())
