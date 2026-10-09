#!/usr/bin/env python3
"""Does the porphine mount hold, and what does it cost? An independent check in UFF.

The CAD builder places the tethers in a simple harmonic field (ADR-0137). This
script re-examines the built carrier site in a different force field, UFF as
implemented in RDKit, which has torsion and inversion (planarity) terms and
iron parameters. It does three things:

1. relax: the site (Fe-porphine, tethers, their hydrogens, the O2) relaxes
   inside a frozen shell of the diamond scaffold. Report the tether geometry and
   the closest non-bonded contacts.
2. strain: the site with its four anchor carbons kept as methyl caps is taken
   out of the scaffold. Its UFF energy at the relaxed in-pocket geometry, minus
   its energy after a free relaxation in vacuum, is the mount strain.
3. md (optional, --md-ps): Langevin dynamics of the site at 310 K in the frozen
   shell. Report the ring drift and tilt, and the iron's motion along the axis.

UFF is a generic force field; for the iron and the porphine it is a rough
model. These are simulated estimates.

    uv run --no-project --with rdkit --with numpy python scripts/mount_uff.py SCENE.json [--md-ps 20]

SCENE is a `carrier_site_json` export with a `bonds` list.
"""

from __future__ import annotations

import argparse
import collections
import json
import math
import sys
import time

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdForceFieldHelpers
from rdkit.Geometry import Point3D

SHELL_A = 7.0
KB_KCAL = 0.0019872043  # kcal/mol/K
KCAL_TO_KJ = 4.184
AMU = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "Fe": 55.845}


def load(path: str):
    d = json.load(open(path))
    el = [a[0] for a in d["atoms"]]
    x = np.array([a[1:] for a in d["atoms"]], dtype=float)
    nbr = collections.defaultdict(list)
    for u, v, _ in d["bonds"]:
        nbr[u].append(v)
        nbr[v].append(u)
    return d, el, x, nbr


def ring_atoms(el, nbr, site: set[int]) -> set[int]:
    """The porphine C and N: site atoms with three heavy neighbours and no
    hydrogen beyond the ring's own beta hydrogens, joined to the macrocycle."""
    ring = set()
    for i in site:
        if el[i] not in ("C", "N"):
            continue
        heavy = [k for k in nbr[i] if el[k] != "H"]
        if el[i] == "N" or (len(heavy) == 3 and sum(el[k] == "H" for k in nbr[i]) == 0):
            ring.add(i)
        elif len(heavy) == 2 and sum(el[k] == "H" for k in nbr[i]) == 1:
            ring.add(i)  # a beta carbon, C-H
    # Keep only atoms joined to another ring atom (drops sp3 tether carbons).
    return {i for i in ring if any(k in ring for k in nbr[i])}


def build_mol(el, x, nbr, keep: list[int], ring: set[int], caps: dict[int, str] | None = None):
    """An RDKit molecule of the atoms `keep`, typed for UFF by hand."""
    caps = caps or {}
    idx = {a: i for i, a in enumerate(keep)}
    m = Chem.RWMol()
    for a in keep:
        at = Chem.Atom(caps.get(a, el[a]))
        at.SetNoImplicit(True)
        m.AddAtom(at)
    seen = set()
    for a in keep:
        for b in nbr[a]:
            if b in idx and (min(a, b), max(a, b)) not in seen:
                seen.add((min(a, b), max(a, b)))
                arom = a in ring and b in ring
                oo = el[a] == "O" and el[b] == "O"
                bt = Chem.BondType.AROMATIC if arom else (Chem.BondType.DOUBLE if oo else Chem.BondType.SINGLE)
                m.AddBond(idx[a], idx[b], bt)
                if arom:
                    m.GetBondBetweenAtoms(idx[a], idx[b]).SetIsAromatic(True)
    for a in keep:
        at = m.GetAtomWithIdx(idx[a])
        e = caps.get(a, el[a])
        deg = at.GetDegree()
        if a in ring:
            at.SetIsAromatic(True)
            at.SetHybridization(Chem.HybridizationType.SP2)
        elif e == "Fe":
            at.SetHybridization(Chem.HybridizationType.SP3D2)
        elif e == "O" and any(el[k] == "O" for k in nbr[a]):
            at.SetHybridization(Chem.HybridizationType.SP2)
        elif e in ("C", "N", "O"):
            at.SetHybridization(Chem.HybridizationType.SP3 if deg >= 2 or e != "C" else Chem.HybridizationType.SP3)
    mol = m.GetMol()
    mol.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(mol)
    conf = Chem.Conformer(len(keep))
    for a in keep:
        p = x[a]
        conf.SetAtomPosition(idx[a], Point3D(*p))
    mol.AddConformer(conf, assignId=True)
    return mol, idx


def positions(mol) -> np.ndarray:
    return np.array(mol.GetConformer().GetPositions())


def ring_frame(p: np.ndarray, ring_idx: list[int]):
    c = p[ring_idx].mean(axis=0)
    _, _, vt = np.linalg.svd(p[ring_idx] - c)
    return c, vt[2]


def angle(a, b, c) -> float:
    u, v = a - b, c - b
    return math.degrees(math.acos(max(-1.0, min(1.0, float(u @ v / np.linalg.norm(u) / np.linalg.norm(v))))))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("scene")
    ap.add_argument("--md-ps", type=float, default=0.0)
    ap.add_argument("--dt-fs", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--hmr", action="store_true", help="hydrogen mass repartitioning (H 3 amu)")
    ap.add_argument("--cut-tethers", action="store_true", help="control: drop the four tether-wall bonds")
    ap.add_argument("--scan-shift-A", type=float, default=0.8)
    ap.add_argument("--scan-deg", type=float, default=16.0)
    ap.add_argument("--scan-points", type=int, default=4, help="points on each side of zero")
    ap.add_argument("--stiffness", action="store_true",
                    help="quasi-static hold test: rigid ring moves, relaxed mount, spring constants")
    a = ap.parse_args()

    d, el, x, nbr = load(a.scene)
    nb = d["block_atoms"]
    if a.cut_tethers:
        # The control: no covalent mount, the ring held by the slot walls alone.
        d["bonds"] = [b for b in d["bonds"] if (b[0] < nb) == (b[1] < nb)]
        nbr = collections.defaultdict(list)
        for u, v, _ in d["bonds"]:
            nbr[u].append(v)
            nbr[v].append(u)
    site = set(range(nb, len(el)))
    ring = ring_atoms(el, nbr, site)
    # The scaffold shell: atoms within SHELL_A of any site atom.
    from scipy.spatial import cKDTree

    tree = cKDTree(x[:nb])
    shell = sorted({j for i in site for j in tree.query_ball_point(x[i], SHELL_A)})
    keep = shell + sorted(site)
    mol, idx = build_mol(el, x, nbr, keep, ring)
    if not rdForceFieldHelpers.UFFHasAllMoleculeParams(mol):
        print("mount_uff: UFF lacks parameters for this molecule", file=sys.stderr)
        return 1
    ff = rdForceFieldHelpers.UFFGetMoleculeForceField(mol, ignoreInterfragInteractions=False)
    for j in shell:
        ff.AddFixedPoint(idx[j])
    ff.Initialize()
    e0 = ff.CalcEnergy()
    t0 = time.time()
    ff.Minimize(maxIts=20000)
    e1 = ff.CalcEnergy()
    p = np.array(ff.Positions()).reshape(-1, 3)
    conf = mol.GetConformer()
    for i in range(len(keep)):
        conf.SetAtomPosition(i, Point3D(*p[i]))
    out = {"scene": a.scene, "slot_radius_A": d.get("slot_radius"), "shell_atoms": len(shell),
           "site_atoms": len(site), "uff_kj_mol_built": round(e0 * KCAL_TO_KJ, 1),
           "uff_kj_mol_relaxed": round(e1 * KCAL_TO_KJ, 1), "relax_s": round(time.time() - t0, 1)}

    # Tether geometry after relaxation.
    anchors = [(u, v) for u, v, _ in d["bonds"] if (u < nb) != (v < nb)]
    tet = []
    for u, v in anchors:
        s, w = (u, v) if u >= nb else (v, u)
        path, seen = [w, s], {w, s}
        while path[-1] not in ring:
            nxt = [k for k in nbr[path[-1]] if k in site and el[k] != "H" and k not in seen]
            if not nxt:
                break
            path.append(nxt[0])
            seen.add(nxt[0])
        q = [p[idx[k]] for k in path]
        bonds = [round(float(np.linalg.norm(q[i + 1] - q[i])), 3) for i in range(len(q) - 1)]
        angs = [round(angle(q[i], q[i + 1], q[i + 2]), 1) for i in range(len(q) - 2)]
        tet.append({"bonds_A": bonds, "angles_deg": angs})
    out["tethers"] = tet
    # Closest site-scaffold non-bonded contacts (not 1-2 or 1-3).
    near12 = {(min(u, v), max(u, v)) for u, v, _ in d["bonds"]}
    near13 = set()
    for b in range(len(el)):
        for i in nbr[b]:
            for k in nbr[b]:
                if i < k:
                    near13.add((i, k))
    worst = []
    pk = {k: p[idx[k]] for k in keep}
    stree = cKDTree(np.array([pk[j] for j in shell]))
    for i in site:
        for jj in stree.query_ball_point(pk[i], 2.4):
            j = shell[jj]
            key = (min(i, j), max(i, j))
            if key in near12 or key in near13:
                continue
            worst.append((round(float(np.linalg.norm(pk[i] - pk[j])), 2), el[i], el[j]))
    worst.sort()
    out["contacts_below_2p4_A"] = len(worst)
    out["closest_contacts"] = worst[:6]

    # Mount strain: the site plus methyl caps on its four anchors, at the
    # in-pocket geometry, against the same molecule relaxed free.
    wall_atoms = [w for u, v in anchors for w in (u, v) if w < nb]
    if not anchors:
        wall_atoms = []
    cap_keep = sorted(site) + wall_atoms
    cap_h = []
    for w in wall_atoms:
        for k in nbr[w]:
            if k < nb and k not in cap_keep:
                cap_h.append((w, k))
    # Replace each anchor's other scaffold neighbours by hydrogens on the same bond line.
    xs = {k: pk[k] for k in cap_keep}
    extra = []
    for w, k in cap_h:
        dvec = pk[k] - pk[w] if k in pk else x[k] - x[w]
        extra.append(pk[w] + 1.09 * dvec / np.linalg.norm(dvec))
    sub_keep = list(cap_keep) + [-(i + 1) for i in range(len(extra))]
    el2 = dict((k, el[k]) for k in cap_keep)
    x2 = dict(xs)
    nbr2 = collections.defaultdict(list)
    for k in cap_keep:
        nbr2[k] = [j for j in nbr[k] if j in xs]
    for i, (w, _) in enumerate(cap_h):
        hid = -(i + 1)
        el2[hid] = "H"
        x2[hid] = extra[i]
        nbr2[hid] = [w]
        nbr2[w].append(hid)
    order = list(range(len(sub_keep)))
    remap = {k: n for n, k in zip(order, sub_keep)}
    el_l = [el2[k] for k in sub_keep]
    x_l = np.array([x2[k] for k in sub_keep])
    nbr_l = collections.defaultdict(list)
    for k in sub_keep:
        nbr_l[remap[k]] = [remap[j] for j in nbr2[k]]
    ring_l = {remap[k] for k in ring}
    mol2, _ = build_mol(el_l, x_l, nbr_l, order, ring_l)
    ff2 = rdForceFieldHelpers.UFFGetMoleculeForceField(mol2)
    e_pocket = ff2.CalcEnergy()
    best = None
    for trial in range(3):
        ff2b = rdForceFieldHelpers.UFFGetMoleculeForceField(mol2)
        if trial:
            pos = np.array(ff2b.Positions()).reshape(-1, 3)
            pos += np.random.default_rng(trial).normal(scale=0.05, size=pos.shape)
            c = mol2.GetConformer()
            for i in range(len(pos)):
                c.SetAtomPosition(i, Point3D(*pos[i]))
            ff2b = rdForceFieldHelpers.UFFGetMoleculeForceField(mol2)
        ff2b.Minimize(maxIts=50000)
        e = ff2b.CalcEnergy()
        best = e if best is None else min(best, e)
    out["mount_strain_kj_mol"] = round((e_pocket - best) * KCAL_TO_KJ, 1)

    # Quasi-static hold test: move the ring core rigidly, relax the rest, and
    # fit the energy. Thermal rms amplitude at 310 K is sqrt(kT / k).
    if a.stiffness:
        core = sorted(ring | {d["iron"]} | {k for k in site if el[k] == "O"}
                      | {k for r in ring for k in nbr[r] if el[k] == "H"})
        core_i = [idx[k] for k in core]
        c0, n0 = ring_frame(p, [idx[k] for k in ring])
        e_ax = n0
        e_in = np.cross(n0, [1.0, 0.0, 0.0])
        if np.linalg.norm(e_in) < 0.5:
            e_in = np.cross(n0, [0.0, 1.0, 0.0])
        e_in /= np.linalg.norm(e_in)

        def relaxed_energy(moved: np.ndarray) -> float:
            c = mol.GetConformer()
            for i in range(len(keep)):
                c.SetAtomPosition(i, Point3D(*moved[i]))
            f2 = rdForceFieldHelpers.UFFGetMoleculeForceField(mol, ignoreInterfragInteractions=False)
            for j in shell:
                f2.AddFixedPoint(idx[j])
            for i in core_i:
                f2.AddFixedPoint(i)
            f2.Initialize()
            f2.Minimize(maxIts=20000)
            return f2.CalcEnergy()

        def rotate(q: np.ndarray, axis: np.ndarray, deg: float) -> np.ndarray:
            th = math.radians(deg)
            k = axis / np.linalg.norm(axis)
            r = q - c0
            return c0 + r * math.cos(th) + np.cross(k, r) * math.sin(th) + np.outer(r @ k, k) * (1 - math.cos(th))

        kt = KB_KCAL * 310.0
        modes = {}
        shift = np.linspace(-a.scan_shift_A, a.scan_shift_A, 2 * a.scan_points + 1)
        turn = np.linspace(-a.scan_deg, a.scan_deg, 2 * a.scan_points + 1)
        for name, steps, move in (
            ("axial_shift_A", shift, lambda q, s: q + s * e_ax),
            ("in_plane_shift_A", shift, lambda q, s: q + s * e_in),
            ("tilt_deg", turn, lambda q, s: rotate(q, e_in, s)),
            ("spin_deg", turn, lambda q, s: rotate(q, e_ax, s)),
        ):
            es = []
            for st in steps:
                q = p.copy()
                q[core_i] = move(p[core_i], st)
                es.append(relaxed_energy(q))
            es = np.array(es) - min(es)
            # The rest position is the scan minimum; the spring constant is a
            # parabola through it and its two neighbours.
            i0 = int(np.argmin(es))
            j = [max(0, min(len(steps) - 3, i0 - 1)) + o for o in range(3)]
            k2 = np.polyfit(np.array(steps)[j], es[j], 2)[0] * 2.0  # kcal/mol per unit^2
            modes[name] = {"rest": round(float(steps[i0]), 3), "rest_at_scan_edge": i0 in (0, len(steps) - 1),
                           "k_kcal_per_unit2": round(float(k2), 3),
                           "rms_310K": round(math.sqrt(kt / k2), 3) if k2 > 0 else None,
                           "drop_to_rest_kj_mol": round(float(es[len(steps) // 2]) * KCAL_TO_KJ, 1),
                           "energies_kj_mol": [round(float(v) * KCAL_TO_KJ, 1) for v in es]}
        out["stiffness"] = modes
        c = mol.GetConformer()
        for i in range(len(keep)):
            c.SetAtomPosition(i, Point3D(*p[i]))

    # Optional Langevin dynamics of the site in the frozen shell.
    if a.md_ps > 0:
        rng = np.random.default_rng(a.seed)
        n = len(keep)
        free = np.array([k in site for k in keep])
        mass = np.array([AMU[el[k]] for k in keep])
        # Hydrogen mass repartitioning: each H takes 2 amu from its heavy
        # neighbour, so the total mass is kept and the fast C-H motion slows.
        # It changes how fast atoms move, not where the ring sits on average.
        if a.hmr:
            for i, k in enumerate(keep):
                if el[k] == "H" and k in site:
                    heavy = [j for j in nbr[k] if j in idx]
                    if heavy:
                        mass[i] += 2.0
                        mass[idx[heavy[0]]] -= 2.0
        T = 310.0
        dt = a.dt_fs * 1e-3  # ps
        gamma = 5.0  # 1/ps
        # Units: kcal/mol, A, amu, ps. Force kcal/mol/A -> acceleration A/ps^2 times 418.4.
        conv = 418.4
        pos = p.copy()
        sig = np.sqrt(KB_KCAL * T * conv / mass)[:, None]
        vel = rng.normal(size=(n, 3)) * sig
        vel[~free] = 0.0
        ring_i = [idx[k] for k in ring]
        fe_i = idx[d["iron"]]
        c0, n0 = ring_frame(pos, ring_i)
        fe0 = pos[fe_i].copy()
        drift, tilt, fe_ax = [], [], []

        nudges = [0]

        def force(q):
            # A UFF angle or inversion derivative is singular at an exactly
            # degenerate geometry and returns NaN; a 1e-5 A nudge of the free
            # atoms steps off it. The count is reported.
            for _ in range(5):
                g = np.array(ff.CalcGrad(q.ravel().tolist())).reshape(-1, 3)
                if np.isfinite(g).all():
                    break
                nudges[0] += 1
                q[free] += rng.normal(scale=1e-5, size=(int(free.sum()), 3))
            f = -g
            f[~free] = 0.0
            return f

        f = force(pos)
        steps = int(round(a.md_ps / dt))
        c1 = math.exp(-gamma * dt)
        c2 = np.sqrt((1 - c1 * c1) * KB_KCAL * T * conv / mass)[:, None]
        t_start = time.time()
        sample = max(1, int(round(0.05 / dt)))
        for step in range(steps):
            # BAOAB
            vel += 0.5 * dt * conv * f / mass[:, None]
            pos += 0.5 * dt * vel
            vel = c1 * vel + c2 * rng.normal(size=(n, 3))
            vel[~free] = 0.0
            pos += 0.5 * dt * vel
            f = force(pos)
            vel += 0.5 * dt * conv * f / mass[:, None]
            vel[~free] = 0.0
            if step % sample == 0:
                c, nn = ring_frame(pos, ring_i)
                drift.append(float(np.linalg.norm(c - c0)))
                tilt.append(math.degrees(math.acos(min(1.0, abs(float(nn @ n0))))))
                fe_ax.append(float((pos[fe_i] - fe0) @ n0))
        out["md"] = {
            "ps": a.md_ps, "seconds": round(time.time() - t_start), "samples": len(drift),
            "ring_drift_A_mean": round(float(np.mean(drift)), 3), "ring_drift_A_max": round(float(np.max(drift)), 3),
            "ring_tilt_deg_mean": round(float(np.mean(tilt)), 2), "ring_tilt_deg_max": round(float(np.max(tilt)), 2),
            "fe_axial_A_sd": round(float(np.std(fe_ax)), 3), "fe_axial_A_range": round(float(np.ptp(fe_ax)), 3),
            "nan_nudges": nudges[0], "dt_fs": a.dt_fs, "hmr": a.hmr,
        }
    print(json.dumps(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
