#!/usr/bin/env python3
"""Strain of the porphine-to-diamond tethers, and the best bridge length.

ADR-0113 bridges each porphine meso carbon to the nearest wall carbon with a
planar zigzag whose bond count is forced even. In the carrier site the span is
about 3.6 A, so four 1.54 A bonds fold to 72-79 degree angles. This script
asks the design question directly: for each meso carbon, which wall anchor and
how many bridge carbons (0 to 3) give the least strain?

For each candidate the two ends are fixed: the meso carbon with its two ring
neighbours, and a wall surface carbon that gives up one hydrogen. The bridge
carbons move to minimize a harmonic valence energy of MM3 scale, the same
constants as `crates/meter/src/valence_relax.rs`:

    bonds   E = 1/2 k (r - r0)^2, C-C 449 N/m, r0 1.544 A (1.51 A to the sp2 meso C)
    angles  E = 1/2 k (t - t0)^2, 0.67e-18 J/rad^2; t0 109.47 deg at sp3, 120 at the meso C
    plane   the meso substituent out of the ring plane, as an angle term toward 0
    steric  a soft wall, 1/2 k (3.0 A - d)^2 for each bridge carbon closer than
            3.0 A to a non-neighbour heavy atom, k 449 N/m

Bridge hydrogens are not modelled; the steric margin stands in for them. This
is a design estimate in a simple field, not quantum chemistry. Simulated.

    python scripts/tether_strain.py ~/.cache/nano-cad/md/carrier_site.json

The scene must carry a `bonds` list (`carrier_site_json` writes one).
"""

from __future__ import annotations

import argparse
import collections
import json
import math
import sys

import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree

AVOGADRO = 6.02214076e23
K_BOND = 449.0  # N/m
R0_CC = 1.544e-10
R0_SP2 = 1.51e-10
K_ANGLE = 0.67e-18  # J/rad^2
T_SP3 = math.acos(-1.0 / 3.0)
T_SP2 = math.radians(120.0)
STERIC_D = 3.0e-10
ANCHOR_RADIUS_A = 6.5


def j_to_kj_mol(e_j: float) -> float:
    return e_j * AVOGADRO / 1000.0


def angle(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    u, v = a - b, c - b
    cosv = float(u @ v / (np.linalg.norm(u) * np.linalg.norm(v)))
    return math.acos(max(-1.0, min(1.0, cosv)))


class Bridge:
    """One candidate: fixed meso and wall ends, n moving bridge carbons."""

    def __init__(self, x_m, meso, ring_nb, wall, wall_nb, n, obstacles):
        self.x = x_m
        self.meso, self.ring_nb, self.wall, self.wall_nb, self.n = meso, ring_nb, wall, wall_nb, n
        self.obstacles = obstacles  # positions, metres, of heavy atoms that are not chain neighbours
        c = x_m[ring_nb].mean(axis=0)
        u, _, vt = np.linalg.svd(np.vstack([x_m[ring_nb], x_m[meso]]) - np.vstack([x_m[ring_nb], x_m[meso]]).mean(0))
        self.ring_normal = vt[2]
        self.meso_dir = (x_m[meso] - c) / np.linalg.norm(x_m[meso] - c)

    def chain(self, flat: np.ndarray) -> list[np.ndarray]:
        pts = [self.x[self.meso]]
        pts += [flat[3 * i: 3 * i + 3] for i in range(self.n)]
        pts.append(self.x[self.wall])
        return pts

    def terms(self, flat: np.ndarray) -> dict:
        p = self.chain(flat)
        e = collections.defaultdict(float)
        for i in range(len(p) - 1):
            r0 = R0_SP2 if i == 0 else R0_CC
            e["bond"] += 0.5 * K_BOND * (np.linalg.norm(p[i + 1] - p[i]) - r0) ** 2
        # The meso carbon: two ring angles at 120 degrees, and in plane.
        for r in self.ring_nb:
            e["angle"] += 0.5 * K_ANGLE * (angle(self.x[r], p[0], p[1]) - T_SP2) ** 2
        d1 = (p[1] - p[0]) / np.linalg.norm(p[1] - p[0])
        e["plane"] += 0.5 * K_ANGLE * math.asin(max(-1.0, min(1.0, float(d1 @ self.ring_normal)))) ** 2
        for i in range(1, len(p) - 1):
            e["angle"] += 0.5 * K_ANGLE * (angle(p[i - 1], p[i], p[i + 1]) - T_SP3) ** 2
        for w in self.wall_nb:
            e["angle"] += 0.5 * K_ANGLE * (angle(self.x[w], p[-1], p[-2]) - T_SP3) ** 2
        if self.n and len(self.obstacles):
            for q in p[1:-1]:
                d = np.linalg.norm(self.obstacles - q, axis=1)
                close = d < STERIC_D
                e["steric"] += float(np.sum(0.5 * K_BOND * (STERIC_D - d[close]) ** 2))
        return e

    def energy(self, flat: np.ndarray) -> float:
        return sum(self.terms(flat).values())

    def solve(self, rng: np.random.Generator, starts: int = 12):
        if self.n == 0:
            flat = np.zeros(0)
            return flat, self.terms(flat)
        a, b = self.x[self.meso], self.x[self.wall]
        best = None
        for _ in range(starts):
            guess = []
            for i in range(1, self.n + 1):
                t = i / (self.n + 1)
                guess.append(a + t * (b - a) + rng.normal(scale=0.6e-10, size=3) + 0.4e-10 * self.meso_dir)
            x0 = np.concatenate(guess) * 1e10  # optimize in angstrom for conditioning
            res = minimize(lambda f: self.energy(f * 1e-10) * 1e20, x0, method="L-BFGS-B")
            if best is None or res.fun < best.fun:
                best = res
        flat = best.x * 1e-10
        return flat, self.terms(flat)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("scene")
    ap.add_argument("--max-bridge", type=int, default=3)
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()

    d = json.load(open(a.scene))
    if "bonds" not in d:
        print("tether_strain: the scene has no bonds list", file=sys.stderr)
        return 1
    el = [r[0] for r in d["atoms"]]
    x = np.array([r[1:] for r in d["atoms"]]) * 1e-10
    nb = d["block_atoms"]
    nbr = collections.defaultdict(list)
    for u, v, _ in d["bonds"]:
        nbr[u].append(v)
        nbr[v].append(u)
    heavy = [i for i, e in enumerate(el) if e != "H"]
    tree = cKDTree(x[heavy] * 1e10)

    # The meso carbons: site carbons bonded to a site sp3 carbon that leads to the wall.
    cross = [(u, v) for u, v, _ in d["bonds"] if (u < nb) != (v < nb)]
    tethers = []
    for u, v in cross:
        s, w = (u, v) if u >= nb else (v, u)
        path, seen = [s], {s, w}
        while True:
            nxt = [k for k in nbr[path[-1]] if k >= nb and el[k] == "C" and k not in seen]
            ring_like = [k for k in nbr[path[-1]] if el[k] != "H"]
            if len(ring_like) == 3 and el[path[-1]] == "C" and not any(el[k] == "H" for k in nbr[path[-1]]):
                break
            if not nxt:
                break
            path.append(nxt[0])
            seen.add(nxt[0])
        # The walk ran wall to meso; a Bridge orders its chain meso to wall.
        tethers.append({"wall": w, "chain": path[:-1][::-1], "meso": path[-1]})

    # Surface wall carbons carrying a hydrogen are anchor candidates.
    def wall_candidates(meso: int) -> list[int]:
        out = []
        for i in range(nb):
            if el[i] != "C":
                continue
            hs = [k for k in nbr[i] if el[k] == "H"]
            if hs and np.linalg.norm(x[i] - x[meso]) * 1e10 < ANCHOR_RADIUS_A:
                out.append(i)
        return out

    rng = np.random.default_rng(7)
    results = []
    for t_i, t in enumerate(tethers):
        meso = t["meso"]
        ring_nb = [k for k in nbr[meso] if k not in t["chain"] and el[k] != "H"]
        built = Bridge(x, meso, ring_nb, t["wall"], [k for k in nbr[t["wall"]] if k not in t["chain"]],
                       len(t["chain"]), np.zeros((0, 3)))
        built_e = j_to_kj_mol(sum(built.terms(np.concatenate([x[c] for c in t["chain"]])).values()))
        rows = []
        for w in wall_candidates(meso):
            w_heavy = [k for k in nbr[w] if el[k] != "H"]
            span = float(np.linalg.norm(x[w] - x[meso]) * 1e10)
            for n in range(0, a.max_bridge + 1):
                # Obstacles: heavy atoms within 6 A, minus the chain ends and their neighbours.
                near = [heavy[j] for j in tree.query_ball_point(x[meso] * 1e10, 7.0)]
                skip = {meso, w, *ring_nb, *w_heavy, *t["chain"]}
                obst = x[[k for k in near if k not in skip]]
                br = Bridge(x, meso, ring_nb, w, w_heavy, n, obst)
                flat, terms = br.solve(rng)
                rows.append({
                    "wall": w, "span_A": round(span, 2), "n_bridge": n,
                    "strain_kj_mol": round(j_to_kj_mol(sum(terms.values())), 1),
                    "terms_kj_mol": {k: round(j_to_kj_mol(v), 1) for k, v in terms.items()},
                })
        rows.sort(key=lambda r: r["strain_kj_mol"])
        best_by_n = {}
        for r in rows:
            best_by_n.setdefault(r["n_bridge"], r)
        results.append({"tether": t_i, "meso": meso, "built_wall": t["wall"], "built_n_bridge": len(t["chain"]),
                        "built_strain_kj_mol": round(built_e, 1), "best": rows[0],
                        "best_by_n": [best_by_n[n] for n in sorted(best_by_n)]})

    if a.json:
        print(json.dumps(results, indent=1))
        return 0
    for r in results:
        print(f"tether {r['tether']}: as built {r['built_n_bridge']} bridge C, strain {r['built_strain_kj_mol']} kJ/mol")
        for b in r["best_by_n"]:
            print(f"    best with {b['n_bridge']} bridge C: anchor {b['wall']} span {b['span_A']} A "
                  f"strain {b['strain_kj_mol']} kJ/mol {b['terms_kj_mol']}")
    return 0


if __name__ == "__main__" and "--ethynyl" not in sys.argv:
    raise SystemExit(main())


# --------------------------------------------------------------------------
# A rigid link: meso-C#C-C(wall), an ethynyl tether (ADR-0137).
# --------------------------------------------------------------------------

ETHYNYL_R0 = (1.43e-10, 1.20e-10, 1.46e-10)  # meso-Csp, C#C, Csp-Csp3
STERIC_H = 2.4e-10


class Ethynyl(Bridge):
    """Two sp carbons, rest angle 180 degrees, no hydrogens on the link."""

    def __init__(self, x_m, meso, ring_nb, wall, wall_nb, obstacles, hydrogens):
        super().__init__(x_m, meso, ring_nb, wall, wall_nb, 2, obstacles)
        self.hydrogens = hydrogens

    def terms(self, flat: np.ndarray) -> dict:
        p = self.chain(flat)
        e = collections.defaultdict(float)
        for i in range(3):
            e["bond"] += 0.5 * K_BOND * (np.linalg.norm(p[i + 1] - p[i]) - ETHYNYL_R0[i]) ** 2
        for r in self.ring_nb:
            e["angle"] += 0.5 * K_ANGLE * (angle(self.x[r], p[0], p[1]) - T_SP2) ** 2
        d1 = (p[1] - p[0]) / np.linalg.norm(p[1] - p[0])
        e["plane"] += 0.5 * K_ANGLE * math.asin(max(-1.0, min(1.0, float(d1 @ self.ring_normal)))) ** 2
        for i in (1, 2):
            e["angle"] += 0.5 * K_ANGLE * (angle(p[i - 1], p[i], p[i + 1]) - math.pi) ** 2
        for w in self.wall_nb:
            e["angle"] += 0.5 * K_ANGLE * (angle(self.x[w], p[3], p[2]) - T_SP3) ** 2
        for q in p[1:3]:
            if len(self.obstacles):
                d = np.linalg.norm(self.obstacles - q, axis=1)
                e["steric"] += float(np.sum(0.5 * K_BOND * (STERIC_D - d[d < STERIC_D]) ** 2))
            if len(self.hydrogens):
                d = np.linalg.norm(self.hydrogens - q, axis=1)
                e["steric"] += float(np.sum(0.5 * K_BOND * (STERIC_H - d[d < STERIC_H]) ** 2))
        return e

    def solve(self, rng: np.random.Generator, starts: int = 6):
        a, b = self.x[self.meso], self.x[self.wall]
        best = None
        for _ in range(starts):
            g = [a + t * (b - a) + rng.normal(scale=0.2e-10, size=3) for t in (0.35, 0.65)]
            x0 = np.concatenate(g) * 1e10
            res = minimize(lambda f: self.energy(f * 1e-10) * 1e20, x0, method="L-BFGS-B")
            if best is None or res.fun < best.fun:
                best = res
        flat = best.x * 1e-10
        return flat, self.terms(flat)


def ethynyl_screen(path: str) -> list[dict]:
    """For each meso carbon of a built site, the best ethynyl anchor in its scaffold."""
    d = json.load(open(path))
    el = [r[0] for r in d["atoms"]]
    x = np.array([r[1:] for r in d["atoms"]]) * 1e-10
    nb = d["block_atoms"]
    nbr = collections.defaultdict(list)
    for u, v, _ in d["bonds"]:
        nbr[u].append(v)
        nbr[v].append(u)
    site = set(range(nb, len(el)))
    # Ring atoms: site C/N with no sp3 character; tether carbons carry two H.
    tether = {i for i in site if el[i] == "C" and sum(el[k] == "H" for k in nbr[i]) == 2}
    tether |= {k for i in tether for k in nbr[i] if el[k] == "H"}
    core = site - tether
    meso = [i for i in core if el[i] == "C" and any(k in tether for k in nbr[i])]
    rng = np.random.default_rng(11)
    out = []
    used_h: set[int] = set()
    for m in meso:
        ring_nb = [k for k in nbr[m] if k in core and el[k] != "H"]
        fe = d["iron"]
        mdir = (x[m] - x[fe]) / np.linalg.norm(x[m] - x[fe])
        best = None
        for w in range(nb):
            if el[w] != "C":
                continue
            span = np.linalg.norm(x[w] - x[m]) * 1e10
            if not 3.4 <= span <= 4.8:
                continue
            if (x[w] - x[m]) @ mdir < 0.8 * np.linalg.norm(x[w] - x[m]):
                continue
            hs = [k for k in nbr[w] if el[k] == "H" and k not in used_h]
            toward = (x[m] - x[w]) / np.linalg.norm(x[m] - x[w])
            hs = [k for k in hs if (x[k] - x[w]) @ toward > 0]
            if not hs:
                continue
            h = max(hs, key=lambda k: float((x[k] - x[w]) @ toward))
            wall_nb = [k for k in nbr[w] if k != h]
            obst = [k for k in range(len(el)) if el[k] != "H" and k not in tether and k not in (m, w)
                    and k not in ring_nb and k not in wall_nb and np.linalg.norm(x[k] - x[m]) < 8e-10]
            hyd = [k for k in range(len(el)) if el[k] == "H" and k not in tether and k != h
                   and np.linalg.norm(x[k] - x[m]) < 8e-10]
            br = Ethynyl(x, m, ring_nb, w, wall_nb, x[obst], x[hyd])
            flat, terms = br.solve(rng)
            tot = sum(terms.values())
            if best is None or tot < best[0]:
                p = br.chain(flat)
                best = (tot, w, h, terms, [round(float(np.linalg.norm(p[i + 1] - p[i]) * 1e10), 3) for i in range(3)],
                        [round(math.degrees(angle(p[i - 1], p[i], p[i + 1])), 1) for i in (1, 2)]
                        + [round(math.degrees(angle(x[k], p[3], p[2])), 1) for k in wall_nb])
        if best is None:
            out.append({"meso": m, "anchor": None})
            continue
        used_h.add(best[2])
        out.append({"meso": m, "anchor": best[1], "span_A": round(float(np.linalg.norm(x[best[1]] - x[m]) * 1e10), 2),
                    "strain_kj_mol": round(j_to_kj_mol(best[3]["bond"] + best[3]["angle"] + best[3]["plane"]), 1),
                    "steric_kj_mol": round(j_to_kj_mol(best[3]["steric"]), 1),
                    "bonds_A": best[4], "angles_deg": best[5]})
    return out


if __name__ == "__main__" and "--ethynyl" in sys.argv:
    for f in [a for a in sys.argv[1:] if not a.startswith("--")]:
        rows = ethynyl_screen(f)
        print(json.dumps({"scene": f, "tethers": rows}))
