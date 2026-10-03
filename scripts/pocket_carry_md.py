"""What does one module pocket carry back across the return seal? (ADR-0120)

A released pocket (ligand arm out, site empty) leaves the tank side through
the seal at 270 degrees and enters the blood. Whatever fills its dead volume
in the tank goes with it. This run puts one pocket of the relaxed hull-wall
scene in a fluid and counts the O2 and the water inside it.

Media:
  --medium gas    O2 gas, `--o2-per-nm3` molecules per nm^3 of free volume.
                  About 14 per nm^3 is 1000 atm at 310 K (an estimate: ideal
                  23.7 per nm^3, compressibility factor about 1.7).
  --medium water  SPC/E water with dissolved O2 at `--o2-molar` (about 1 M
                  at 1000 atm, a Henry's-law extrapolation).

The pocket is cut from the scene as in `module_pocket_water.py`: the rotor
lattice within `--cut-radius` of the pocket axis, the pocket's site and its
ligand arm, all frozen. The pocket frame comes from the iron position, so
any pocket index works. The pocket region is the ring slot plus the bore to
the rim; the arm channel is reported apart.

Model: UFF Lennard-Jones on the solid, no solid charges, Lorentz-Berthelot
mixing; O2 two neutral LJ sites (sigma 3.02 A, epsilon 0.407 kJ/mol); a
Morse term between the iron and each water oxygen (ADR-0086 proxy). No
Fe-O2 bond: a released site binds O2 weakly (ADR-0106). One process, below
0.5 GB. A simulated estimate.

    python scripts/pocket_carry_md.py relaxed.json --pocket 5 --medium gas
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
UFF = {"C": (3.851, 0.105), "H": (2.886, 0.044), "N": (3.660, 0.069), "O": (3.500, 0.060), "Fe": (2.912, 0.013)}
SPCE_SIGMA, SPCE_EPS, SPCE_Q_O, SPCE_Q_H = 0.3166, 0.6502, -0.8476, 0.4238
OH, HH = 0.1, 0.1633
O2_SIGMA, O2_EPS, O2_BOND = 0.302, 0.407, 0.121
MORSE_D_KJ, MORSE_R0_NM, MORSE_A_PER_NM = 10.3, 0.270, 27.437
WATER_VOLUME_A3 = 29.9
AVOGADRO = 6.02214076e23
# Pocket geometry (heme_rotor defaults), nm; local z along the pocket axis,
# +z toward the rim, the iron at z = 0.
DISC, SLOT_R, SLOT_H, BORE_R = 4.0, 0.90, 0.36, 0.56
BORE_TOP = 0.66
CHAN_R, CHAN_BACK = 0.81, -3.0


def uff(element: str) -> tuple[float, float]:
    x1, d1 = UFF[element]
    return x1 / LJ_FACTOR / 10.0, d1 * KCAL


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("scene")
    p.add_argument("--pocket", type=int, default=5)
    p.add_argument("--medium", choices=["gas", "water"], default="gas")
    p.add_argument("--o2-per-nm3", type=float, default=14.0)
    p.add_argument("--o2-molar", type=float, default=1.0)
    p.add_argument("--o2-spacing", type=float, default=0.34,
                   help="nm, the gas start lattice; 0.29 reaches a 1000 atm bulk after adsorption")
    p.add_argument("--cut-radius", type=float, default=1.6)
    p.add_argument("--bore-len", type=float, default=0.30, help="nm, the rotor's bore_length_m")
    p.add_argument("--add-bound-o2", action="store_true",
                   help="put a frozen O2 on the iron at the bound geometry (Fe-O 1.85 A on the pocket "
                        "axis, O-O 1.23 A at 60 degrees): the loaded state in the same pocket lattice")
    p.add_argument("--remove-bound-o2", action="store_true",
                   help="drop the site's bound O2 atoms: the released cavity in a loaded pocket's lattice")
    p.add_argument("--map", action="store_true", help="print where the O2 sits in the pocket (r, lz)")
    p.add_argument("--water-rotor", action="store_true",
                   help="the scene is the water-return rotor (no sites): take the pocket frame from the "
                        "rotor geometry, and give its hydroxyl groups the model charges")
    p.add_argument("--hydrophilic", action="store_true",
                   help="water-return pocket (ADR-0122): drop the heme site and turn every rotor cap "
                        "hydrogen that faces the pocket into a frozen hydroxyl with model charges")
    p.add_argument("--z-back", type=float, default=-3.1)
    p.add_argument("--z-front", type=float, default=1.4)
    p.add_argument("--equilibrate-ps", type=float, default=50.0)
    p.add_argument("--production-ps", type=float, default=500.0)
    p.add_argument("--sample-ps", type=float, default=1.0)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--platform", default="OpenCL")
    a = p.parse_args()

    data = json.load(open(a.scene))
    atoms = data.get("scene", data)["atomistic"]["atoms"]
    body = np.array([x["body"] for x in atoms])
    element = [x["element"] for x in atoms]
    pos = np.array([x["position_m"] for x in atoms]) * 1e9
    pockets = (len(set(body.tolist())) - 2) // 2
    if a.water_rotor:
        # The plunger rotor (ADR-0122): 6 pockets from 15 degrees; the frame
        # origin is on the pocket axis at the rim, and the plunger is body
        # 2 + k.
        k = a.pocket % 6
        site_body, arm_body = -1, 2 + k
        az0 = math.radians(15.0 + 60.0 * k)
        fe = DISC * np.array([math.cos(az0), math.sin(az0), 0.0])
        fe_n = float("nan")
    else:
        k = a.pocket % pockets
        site_body, arm_body = 2 + k, 2 + pockets + k
        iron_idx = [i for i in np.where(body == site_body)[0] if element[i] == "Fe"]
        if not iron_idx:
            print("no iron in that pocket", file=sys.stderr)
            return 1
        fe = pos[iron_idx[0]]
        arm_n = [i for i in np.where(body == arm_body)[0] if element[i] == "N"]
        fe_n = min(np.linalg.norm(pos[i] - fe) for i in arm_n)
    az = math.atan2(fe[1], fe[0])
    u = np.array([math.cos(az), math.sin(az), 0.0])
    t = np.array([-math.sin(az), math.cos(az), 0.0])
    w = np.array([0.0, 0.0, 1.0])
    rel = pos - fe
    local = np.column_stack([rel @ t, rel @ w, rel @ u])
    radial = np.hypot(local[:, 0], local[:, 1])
    keep = (body == 1) & (radial < a.cut_radius) & (local[:, 2] > a.z_back) & (local[:, 2] < a.z_front)
    keep |= (body == site_body) | (body == arm_body)
    if a.hydrophilic:
        keep &= body != site_body
    if a.remove_bound_o2:
        keep &= ~((body == site_body) & (np.array(element) == "O"))
    solid = local[keep]
    solid_el = [element[i] for i in np.where(keep)[0]]
    if a.add_bound_o2:
        # The local frame has the iron at the origin and +z toward the rim.
        o1 = np.array([0.0, 0.0, 0.185])
        o2b = o1 + 0.123 * np.array([math.sin(math.radians(60)), 0.0, math.cos(math.radians(60))])
        solid = np.vstack([solid, o1, o2b])
        solid_el = solid_el + ["O", "O"]
    solid_q = np.zeros(len(solid))
    iron = solid_el.index("Fe") if "Fe" in solid_el else None
    if a.water_rotor:
        # Model charges on the built hydroxyls (group_data.rs): host C
        # +0.0319 e, O -0.3996 e, H +0.2096 e.
        t0 = cKDTree(solid)
        n_oh = 0
        for i, e in enumerate(solid_el):
            if e != "O":
                continue
            nb_i = t0.query_ball_point(solid[i], 0.155)
            hs = [j for j in nb_i if solid_el[j] == "H" and np.linalg.norm(solid[j] - solid[i]) < 0.105]
            cs = [j for j in nb_i if solid_el[j] == "C"]
            if hs and cs:
                solid_q[i] += -0.3996
                solid_q[hs[0]] += 0.2096
                solid_q[cs[0]] += 0.0319
                n_oh += 1
        print(f"water rotor: {n_oh} hydroxyls charged", file=sys.stderr)
    if a.hydrophilic:
        # Hydroxyl charges of the parts library (group_data.rs, hydroxyl):
        # host C +0.0319 e, O -0.3996 e, H +0.2096 e.
        q_c, q_o, q_h = 0.0319, -0.3996, 0.2096
        tree0 = cKDTree(solid)
        extra_xyz, extra_el, extra_q = [], [], []
        n_oh = 0
        for i, e in enumerate(solid_el):
            if e != "H":
                continue
            x = solid[i]
            r_ax = math.hypot(x[0], x[1])
            if not (r_ax < SLOT_R + 0.15 and 0.0 <= x[2] < SLOT_H + a.bore_len + 0.1):
                continue
            nb_i = [j for j in tree0.query_ball_point(x, 0.12) if solid_el[j] == "C"]
            if not nb_i:
                continue
            c = solid[nb_i[0]]
            d = (x - c) / np.linalg.norm(x - c)
            o = c + 0.143 * d
            helper = np.array([0.0, 0.0, 1.0]) if abs(d[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
            e1 = np.cross(d, helper); e1 /= np.linalg.norm(e1)
            h = o + 0.096 * (math.cos(math.radians(71.5)) * d + math.sin(math.radians(71.5)) * e1)
            solid[i] = o
            solid_el[i] = "O"
            solid_q[i] = q_o
            solid_q[nb_i[0]] += q_c
            extra_xyz.append(h); extra_el.append("H"); extra_q.append(q_h)
            n_oh += 1
        if extra_xyz:
            solid = np.vstack([solid, np.array(extra_xyz)])
            solid_el += extra_el
            solid_q = np.concatenate([solid_q, np.array(extra_q)])
        print(f"hydrophilic pocket: {n_oh} hydroxyls", file=sys.stderr)

    margin = 0.6
    box_xy = 2.0 * a.cut_radius + margin
    box_z = (a.z_front - a.z_back) + margin
    shift = np.array([box_xy / 2, box_xy / 2, -a.z_back + margin / 2])
    solid_s = solid + shift
    tree = cKDTree(solid_s)
    rng = np.random.default_rng(a.seed)

    def lattice(spacing_nm: float, clearance_nm: float) -> np.ndarray:
        n = [max(1, round(L / spacing_nm)) for L in (box_xy, box_xy, box_z)]
        s = [L / m for L, m in zip((box_xy, box_xy, box_z), n)]
        g = np.stack(np.meshgrid(*[(np.arange(m) + 0.5) * d for m, d in zip(n, s)], indexing="ij"), -1).reshape(-1, 3)
        return g[tree.query(g)[0] >= clearance_nm]

    if a.water_rotor:
        # The rotor interior behind the plunger is sealed from the pocket
        # (ADR-0122); fluid placed in its cut-open back end only jams.
        base_lattice = lattice

        def lattice(spacing_nm: float, clearance_nm: float) -> np.ndarray:  # noqa: F811
            g = base_lattice(spacing_nm, clearance_nm)
            return g[(g[:, 2] - shift[2]) > -0.95]

    def random_dir() -> np.ndarray:
        d = rng.normal(size=3)
        return d / np.linalg.norm(d)

    o2_centres, water_sites = [], []
    if a.medium == "gas":
        free = lattice(0.1, 0.30)
        voxel_nm3 = box_xy * box_xy * box_z / (round(box_xy / 0.1) ** 2 * round(box_z / 0.1))
        free_volume = len(free) * voxel_nm3
        n_o2 = int(round(a.o2_per_nm3 * free_volume))
        sites = lattice(a.o2_spacing, 0.30)
        o2_centres = sites[rng.choice(len(sites), size=min(n_o2, len(sites)), replace=False)]
    else:
        sites = lattice((WATER_VOLUME_A3 ** (1 / 3)) / 10.0, 0.29)
        n_o2 = max(1, int(round(a.o2_molar * len(sites) / 55.5)))
        pick = set(rng.choice(len(sites), size=n_o2, replace=False).tolist())
        o2_centres = sites[sorted(pick)]
        water_sites = sites[[i for i in range(len(sites)) if i not in pick]]

    system = mm.System()
    system.setDefaultPeriodicBoxVectors(mm.Vec3(box_xy, 0, 0), mm.Vec3(0, box_xy, 0), mm.Vec3(0, 0, box_z))
    nb = mm.NonbondedForce()
    nb.setNonbondedMethod(mm.NonbondedForce.PME)
    nb.setCutoffDistance(1.0)
    coords = list(solid_s)
    for e, q in zip(solid_el, solid_q):
        system.addParticle(0.0)
        sg, ep = uff(e)
        nb.addParticle(float(q), sg, ep)
    o2_first = []
    for c in o2_centres:
        d = random_dir()
        i = system.addParticle(15.999)
        j = system.addParticle(15.999)
        nb.addParticle(0.0, O2_SIGMA, O2_EPS)
        nb.addParticle(0.0, O2_SIGMA, O2_EPS)
        nb.addException(i, j, 0.0, 0.1, 0.0)
        system.addConstraint(i, j, O2_BOND)
        coords += [c - 0.5 * O2_BOND * d, c + 0.5 * O2_BOND * d]
        o2_first.append(i)
    water_o = []
    half = math.radians(109.47 / 2)
    for o in water_sites:
        d = random_dir()
        r = random_dir()
        r -= np.dot(r, d) * d
        r /= np.linalg.norm(r)
        i = system.addParticle(15.999)
        h1 = system.addParticle(1.008)
        h2 = system.addParticle(1.008)
        nb.addParticle(SPCE_Q_O, SPCE_SIGMA, SPCE_EPS)
        nb.addParticle(SPCE_Q_H, 0.1, 0.0)
        nb.addParticle(SPCE_Q_H, 0.1, 0.0)
        for x, y in ((i, h1), (i, h2), (h1, h2)):
            nb.addException(x, y, 0.0, 0.1, 0.0)
        system.addConstraint(i, h1, OH)
        system.addConstraint(i, h2, OH)
        system.addConstraint(h1, h2, HH)
        coords += [o, o + OH * (math.cos(half) * d + math.sin(half) * r), o + OH * (math.cos(half) * d - math.sin(half) * r)]
        water_o.append(i)
        if iron is not None:
            nb.addException(iron, i, 0.0, 0.1, 0.0)
    for i, j in tree.query_pairs(1.0):
        nb.addException(int(i), int(j), 0.0, 0.1, 0.0)
    system.addForce(nb)
    if water_o and iron is not None:
        morse = mm.CustomBondForce("D*((1-exp(-a*(r-r0)))^2-1)")
        morse.addPerBondParameter("D")
        morse.addGlobalParameter("a", MORSE_A_PER_NM)
        morse.addGlobalParameter("r0", MORSE_R0_NM)
        morse.setUsesPeriodicBoundaryConditions(True)
        for o in water_o:
            morse.addBond(iron, o, [MORSE_D_KJ])
        system.addForce(morse)
    # Keep the fluid inside the box along z (no periodic short cut).
    wall = mm.CustomExternalForce("k*(step(zlo-z)*(zlo-z)^2 + step(z-zhi)*(z-zhi)^2)")
    wall.addGlobalParameter("k", 5000.0)
    wall.addGlobalParameter("zlo", 0.0)
    wall.addGlobalParameter("zhi", float(box_z))
    for i in o2_first:
        wall.addParticle(i, [])
        wall.addParticle(i + 1, [])
    for i in water_o:
        wall.addParticle(i, [])
    system.addForce(wall)
    integ = mm.LangevinMiddleIntegrator(310 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
    integ.setRandomNumberSeed(a.seed)
    try:
        plat = mm.Platform.getPlatformByName(a.platform)
    except Exception:  # noqa: BLE001
        plat = mm.Platform.getPlatformByName("CPU")
    ctx = mm.Context(system, integ, plat)
    ctx.setPositions(np.array(coords) * unit.nanometer)
    mm.LocalEnergyMinimizer.minimize(ctx, 10.0, 1000)
    ctx.setVelocitiesToTemperature(310 * unit.kelvin, a.seed)
    integ.step(int(a.equilibrate_ps / 0.002))

    def region(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        lx, ly, lz = x[:, 0] - shift[0], x[:, 1] - shift[1], x[:, 2] - shift[2]
        r = np.hypot(lx, ly)
        if a.water_rotor:
            # The bore in front of the plunger tip, up to the rim.
            pocket = (r < 0.66) & (lz > -0.80) & (lz < 0.0)
            arm = (r < 0.66) & (lz > CHAN_BACK) & (lz <= -0.80)
            return pocket, arm
        # The pocket that opens to the rim: above the ring plane and below
        # the rim. The slot below the ring belongs to the sealed interior
        # behind the arm collar (ADR-0121); it is counted with the arm side.
        pocket = (r < SLOT_R) & (lz >= 0.0) & (lz < SLOT_H + a.bore_len)
        arm = (r < SLOT_R) & (lz > CHAN_BACK) & (lz < 0.0)
        return pocket, arm

    o2_pocket, o2_arm, w_pocket, w_arm, bulk_o2 = [], [], [], [], []
    map_r_edges = np.arange(0.0, 1.01, 0.1)
    map_z_edges = np.arange(-0.4, 0.81, 0.1)
    o2_map = np.zeros((len(map_r_edges) - 1, len(map_z_edges) - 1))
    o2_idx = np.array(o2_first, dtype=int)
    for _ in range(int(a.production_ps / a.sample_ps)):
        integ.step(int(a.sample_ps / 0.002))
        x = ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        if len(o2_idx):
            c = 0.5 * (x[o2_idx] + x[o2_idx + 1])
            pk, am = region(c)
            o2_pocket.append(int(pk.sum()))
            if a.map:
                lx, ly, lzz = c[:, 0] - shift[0], c[:, 1] - shift[1], c[:, 2] - shift[2]
                o2_map += np.histogram2d(np.hypot(lx, ly), lzz, bins=[map_r_edges, map_z_edges])[0]
            o2_arm.append(int(am.sum()))
            # Bulk density far above the rim (local z 0.9-1.3 nm), for reference.
            lz = c[:, 2] - shift[2]
            bulk_o2.append(int(((lz > 0.9) & (lz < 1.3)).sum()))
        if water_o:
            pk, am = region(x[np.array(water_o)])
            w_pocket.append(int(pk.sum()))
            w_arm.append(int(am.sum()))

    def mean(v: list[int]) -> float:
        return round(float(np.mean(v)), 3) if v else 0.0

    slab = 0.4 * box_xy * box_xy
    print(json.dumps({
        "pocket": k,
        "pockets": 6 if a.water_rotor else pockets,
        "azimuth_deg": round(math.degrees(az) % 360, 1),
        "fe_n_A": None if a.water_rotor else round(fe_n * 10, 2),
        "bound_o2_added": a.add_bound_o2,
        "bound_o2_removed": a.remove_bound_o2,
        "hydrophilic": a.hydrophilic,
        "medium": a.medium,
        "o2": len(o2_first),
        "waters": len(water_o),
        "production_ps": a.production_ps,
        "o2_in_pocket": mean(o2_pocket),
        "o2_in_pocket_max": max(o2_pocket) if o2_pocket else 0,
        "o2_in_arm_channel": mean(o2_arm),
        "water_in_pocket": mean(w_pocket),
        "water_in_arm_channel": mean(w_arm),
        "bulk_o2_per_nm3_far": round(float(np.mean(bulk_o2)) / slab, 2) if bulk_o2 else 0.0,
        "platform": plat.getName(),
    }))
    if a.map:
        n = max(1, len(o2_pocket))
        print("O2 per sample; rows r (nm from the pocket axis), columns local z (nm above the iron)")
        print("r\\z   " + " ".join(f"{z:5.1f}" for z in map_z_edges[:-1]))
        for i, r0 in enumerate(map_r_edges[:-1]):
            print(f"{r0:4.1f}  " + " ".join(f"{v / n:5.2f}" for v in o2_map[i]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
