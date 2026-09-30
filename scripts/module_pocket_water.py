"""M28-26: water in one pocket of the relaxed O2 module (ADR-0111).

The question: how much water does one as-built module pocket carry? The
pocket is the ring slot plus the bore to the rim, with the arm channel
behind it. The rotor, the pocket's heme site and its ligand arm are frozen.
SPC/E water fills the space around the piece. The run counts the water
oxygens in the pocket and in the arm gap.

Input is the scene JSON that `module_json relaxed <path>` writes:

    cargo run --release -p nanocad-jigs --example module_json -- relaxed module.json
    python scripts/module_pocket_water.py module.json --pocket 0

`--pocket 0` is a loaded pocket (ligand in, O2 bound); `--pocket 6` is a
released pocket (ligand out, empty). `--start wet` fills the pocket with
water at the start; `--start dry` leaves it empty. The two starts bracket a
slow filling or drying.

The model: UFF Lennard-Jones on the solid, Lorentz-Berthelot mixing, no
solid charges; a Morse term between the iron and each water oxygen (10.3
kJ/mol at 2.7 A, the Fe(NH3)4 proxy of ADR-0086) that replaces their
Lennard-Jones pair. One OpenMM process, below 0.5 GB. A simulated estimate.

Needs OpenMM, which is not a dependency of the package.
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
# UFF x1 (A) and D1 (kcal/mol), Rappe et al. 1992.
UFF = {"C": (3.851, 0.105), "H": (2.886, 0.044), "N": (3.660, 0.069), "O": (3.500, 0.060), "Fe": (2.912, 0.013)}
# SPC/E water.
SPCE_SIGMA_NM = 0.3166
SPCE_EPSILON_KJ = 0.6502
SPCE_Q_O = -0.8476
SPCE_Q_H = 0.4238
SPCE_OH_NM = 0.1
SPCE_HH_NM = 0.1633
MORSE_D_KJ = 10.3
MORSE_R0_NM = 0.270
MORSE_A_PER_NM = 27.437
WATER_VOLUME_A3 = 29.9
# Module geometry (heme_rotor defaults), nm. Local z is the pocket axis,
# +z toward the rim, the iron at z = 0.
DISC, SLOT_R, SLOT_H, BORE_R = 4.0, 0.90, 0.36, 0.56
ARM_R, ARM_CLR, ARM_L, LIG_OUT = 0.26, 0.55, 1.2, 0.40
IRON_R = DISC - SLOT_H - 0.30  # 3.34 nm, the iron radius
BORE_TOP = DISC - IRON_R  # the rim, local z = 0.66
CHAN_R = ARM_R + ARM_CLR  # 0.81
CHAN_BACK = -(LIG_OUT + ARM_L + 1.0)  # -2.6
POCKET_COUNT = 12


def uff_sigma_epsilon(element: str) -> tuple[float, float]:
    x1, d1 = UFF[element]
    return x1 / LJ_FACTOR / 10.0, d1 * KCAL


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("scene")
    parser.add_argument("--pocket", type=int, default=0)
    parser.add_argument("--cut-radius", type=float, default=1.6, help="nm")
    parser.add_argument("--z-back", type=float, default=-2.9, help="nm, past the arm channel")
    parser.add_argument("--z-front", type=float, default=1.0, help="nm, past the rim")
    parser.add_argument("--start", choices=["wet", "dry"], default="wet")
    parser.add_argument("--equilibrate-ps", type=float, default=50.0)
    parser.add_argument("--production-ps", type=float, default=500.0)
    parser.add_argument("--sample-ps", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--platform", default="OpenCL")
    args = parser.parse_args()

    scene = json.load(open(args.scene))
    atoms = scene["atomistic"]["atoms"]
    body_of = np.array([a["body"] for a in atoms])
    elements = [a["element"] for a in atoms]
    positions = np.array([a["position_m"] for a in atoms]) * 1e9  # nm

    k = args.pocket % POCKET_COUNT
    site_body, arm_body = 2 + k, 14 + k
    az = 2.0 * math.pi * k / POCKET_COUNT
    u = np.array([math.cos(az), math.sin(az), 0.0])
    t = np.array([-math.sin(az), math.cos(az), 0.0])
    w = np.array([0.0, 0.0, 1.0])
    fe = IRON_R * u

    rel = positions - fe
    local = np.column_stack([rel @ t, rel @ w, rel @ u])
    radial = np.hypot(local[:, 0], local[:, 1])
    keep = np.zeros(len(atoms), dtype=bool)
    keep |= (body_of == 1) & (radial < args.cut_radius) & (local[:, 2] > args.z_back) & (local[:, 2] < args.z_front)
    keep |= body_of == site_body
    keep |= body_of == arm_body
    solid = local[keep]
    solid_elements = [elements[i] for i in np.where(keep)[0]]
    if "Fe" not in solid_elements:
        print("no iron in the pocket; is this a module scene?", file=sys.stderr)
        return 1
    iron = solid_elements.index("Fe")

    # The box: the cut cylinder plus a margin, periodic in x and y.
    margin = 0.6
    box_xy = 2.0 * args.cut_radius + margin
    box_z = (args.z_front - args.z_back) + margin
    z_shift = -args.z_back + margin / 2
    solid_s = solid.copy()
    solid_s[:, 0] += box_xy / 2
    solid_s[:, 1] += box_xy / 2
    solid_s[:, 2] += z_shift

    spacing = (WATER_VOLUME_A3 ** (1.0 / 3.0)) / 10.0
    tree = cKDTree(solid_s)
    rng = np.random.default_rng(args.seed)
    waters = []
    n_xy = int(box_xy / spacing)
    n_z = int(box_z / spacing)
    for ix in range(n_xy):
        for iy in range(n_xy):
            for iz in range(n_z):
                o = np.array([(ix + 0.5) * spacing, (iy + 0.5) * spacing, (iz + 0.5) * spacing])
                if tree.query(o)[0] < 0.26:
                    continue
                lx, ly, lz = o[0] - box_xy / 2, o[1] - box_xy / 2, o[2] - z_shift
                r = math.hypot(lx, ly)
                in_pocket = (r < SLOT_R and -SLOT_H < lz < SLOT_H) or (r < BORE_R and SLOT_H < lz < BORE_TOP)
                in_arm = r < CHAN_R and CHAN_BACK < lz < -SLOT_H
                if args.start == "dry" and (in_pocket or in_arm):
                    continue
                waters.append(o)
    positions_all = [p for p in solid_s]
    for o in waters:
        q = rng.normal(size=3)
        q /= np.linalg.norm(q)
        r = rng.normal(size=3)
        r -= np.dot(r, q) * q
        r /= np.linalg.norm(r)
        half = math.radians(109.47 / 2)
        h1 = o + SPCE_OH_NM * (math.cos(half) * q + math.sin(half) * r)
        h2 = o + SPCE_OH_NM * (math.cos(half) * q - math.sin(half) * r)
        positions_all += [o, h1, h2]
    positions_all = np.array(positions_all)

    system = mm.System()
    system.setDefaultPeriodicBoxVectors(mm.Vec3(box_xy, 0, 0), mm.Vec3(0, box_xy, 0), mm.Vec3(0, 0, box_z))
    nb = mm.NonbondedForce()
    nb.setNonbondedMethod(mm.NonbondedForce.PME)
    nb.setCutoffDistance(1.0)
    nb.setUseDispersionCorrection(False)
    for element in solid_elements:
        system.addParticle(0.0)  # frozen
        sigma, epsilon = uff_sigma_epsilon(element)
        nb.addParticle(0.0, sigma, epsilon)
    first_water = len(solid_elements)
    water_oxygens = []
    for _ in range(len(waters)):
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
        nb.addException(iron, o, 0.0, 0.1, 0.0)  # replaced by the Morse term
    tree_frozen = cKDTree(positions_all[:first_water])
    for i, j in tree_frozen.query_pairs(1.0):
        nb.addException(int(i), int(j), 0.0, 0.1, 0.0)
    system.addForce(nb)

    morse = mm.CustomBondForce("D*((1-exp(-a*(r-r0)))^2-1)")
    morse.addPerBondParameter("D")
    morse.addGlobalParameter("a", MORSE_A_PER_NM)
    morse.addGlobalParameter("r0", MORSE_R0_NM)
    morse.setUsesPeriodicBoundaryConditions(True)
    for o in water_oxygens:
        morse.addBond(iron, o, [MORSE_D_KJ])
    system.addForce(morse)

    wall = mm.CustomExternalForce("k*(step(zlo-z)*(zlo-z)^2 + step(z-zhi)*(z-zhi)^2)")
    wall.addGlobalParameter("k", 5000.0)
    wall.addGlobalParameter("zlo", 0.0)
    wall.addGlobalParameter("zhi", float(box_z))
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
    context.setPositions(positions_all * unit.nanometer)
    mm.LocalEnergyMinimizer.minimize(context, 10.0, 500)
    context.setVelocitiesToTemperature(310 * unit.kelvin, args.seed)
    integrator.step(int(args.equilibrate_ps / 0.002))

    def oxygen_positions():
        state = context.getState(getPositions=True)
        pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        ox = pos[first_water::3]
        lx = ox[:, 0] - box_xy / 2
        ly = ox[:, 1] - box_xy / 2
        lz = ox[:, 2] - z_shift
        return np.hypot(lx, ly), lz

    steps_per_sample = int(args.sample_ps / 0.002)
    pocket_samples, arm_samples, iron_contacts = [], [], []
    for _ in range(int(args.production_ps / args.sample_ps)):
        integrator.step(steps_per_sample)
        r, lz = oxygen_positions()
        in_pocket = ((r < SLOT_R) & (lz > -SLOT_H) & (lz < SLOT_H)) | ((r < BORE_R) & (lz >= SLOT_H) & (lz < BORE_TOP))
        in_arm = (r < CHAN_R) & (lz > CHAN_BACK) & (lz < -SLOT_H)
        pocket_samples.append(int(in_pocket.sum()))
        arm_samples.append(int(in_arm.sum()))
        d_fe = np.hypot(r, lz)
        iron_contacts.append(int((d_fe < 0.32).sum()))

    pocket_samples = np.array(pocket_samples)
    arm_samples = np.array(arm_samples)
    result = {
        "pocket": k,
        "loaded": k < 6,
        "start": args.start,
        "waters": len(waters),
        "solid_atoms": first_water,
        "production_ps": args.production_ps,
        "mean_waters_in_pocket": round(float(pocket_samples.mean()), 3),
        "max_waters_in_pocket": int(pocket_samples.max()),
        "mean_waters_in_arm": round(float(arm_samples.mean()), 3),
        "fraction_pocket_empty": round(float((pocket_samples == 0).mean()), 3),
        "mean_iron_contacts": round(float(np.mean(iron_contacts)), 3),
        "platform": platform.getName(),
    }
    print(json.dumps(result))
    return 0


if __name__ == "__main__":
    sys.exit(main())
