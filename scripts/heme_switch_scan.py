"""M28-14: the strain switch on a heme-like Fe-N4 site.

The site is an idealized Fe-porphine (D4h, Fe-N 2.00 A) with an axial NH3
below the iron at 2.10 A. A cam that pulls the axial ligand domes the iron h
Angstrom out of the ring plane, away from the O2 side; the ring stays planar
and rigid (the atoms that hold the ring are the rotor). O2 binds on top, end
on, bent at 60 degrees from the axis.

For each dome h: E_bind = min E(site + O2) - min E(site) - E(O2), and the
strain cost E(h) - E(0) of the empty and the loaded site.

Method: UKS with density fitting (def2-universal-jkfit), def2-SVP on Fe, N
and O, 6-31G on ring carbons, STO-3G on hydrogens, grid level 1. `--xc` chooses the functional. Rigid, idealized geometry.

Memory safety as in `fe_o2_switch_scan.py`: 1 GB PySCF budget per worker,
worker count capped at 40 percent of physical memory, a watchdog, one task
per worker process, `--dry-run`.
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

import argparse
import json
import multiprocessing as mp
import sys
import threading
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fe_o2_switch_scan import (  # noqa: E402
    MEMORY_FRACTION,
    WORKER_MAX_MEMORY_MB,
    WORKER_PLANNED_MB,
    physical_memory_mb,
    watchdog,
)

O2_BOND = 1.23
O2_TILT_DEG = 60.0
AXIAL_N = 2.10


def rot(p, k):
    c, s = np.cos(k * np.pi / 2), np.sin(k * np.pi / 2)
    return (c * p[0] - s * p[1], s * p[0] + c * p[1], p[2])


def porphine(core=2.00):
    """Idealized planar porphine dianion core, z = 0, scaled so that Fe-N is
    `core` Angstrom (a uniform squeeze or stretch of the whole ring)."""
    quarter = [
        ("N", (2.00, 0.00, 0.0)),
        ("C", (2.83, 1.10, 0.0)),
        ("C", (2.83, -1.10, 0.0)),
        ("C", (4.23, 0.675, 0.0)),
        ("C", (4.23, -0.675, 0.0)),
        ("H", (5.13, 1.27, 0.0)),
        ("H", (5.13, -1.27, 0.0)),
        ("C", (2.43, 2.43, 0.0)),
        ("H", (3.19, 3.19, 0.0)),
    ]
    f = core / 2.00
    quarter = [(e, (p[0] * f, p[1] * f, p[2])) for e, p in quarter]
    out = []
    for k in range(4):
        out += [(e, rot(p, k)) for e, p in quarter]
    return out


def nh3_below(fe_z, axial=AXIAL_N):
    n = np.array([0.0, 0.0, fe_z - axial])
    out = [("N", tuple(n))]
    for k in range(3):
        a = np.radians(120 * k)
        out.append(("H", tuple(n + 1.01 * np.array([0.94 * np.cos(a), 0.94 * np.sin(a), -0.34]))))
    return out


def site(h, axial=AXIAL_N, core=2.00):
    fe_z = -h
    return [("Fe", (0.0, 0.0, fe_z))] + porphine(core) + nh3_below(fe_z, axial), fe_z


def o2(fe_z, d):
    t = np.radians(O2_TILT_DEG)
    z = fe_z + d
    return [("O", (0.0, 0.0, z)), ("O", (O2_BOND * np.sin(t), 0.0, z + O2_BOND * np.cos(t)))]


BASIS = {"Fe": "def2-svp", "N": "def2-svp", "O": "def2-svp", "C": "6-31g", "H": "sto-3g"}


def scf(atoms, charge, spin2, xc):
    """Damped, level-shifted first pass, then a second-order (Newton) finish.

    The iron, its donors and the O2 carry def2-SVP; the ring carbons carry
    6-31G and the hydrogens STO-3G, to keep the basis near 300 functions.
    Peak resident memory measured at 1.7 GB per process (2026-09-27).
    """
    from pyscf import dft, gto

    mol = gto.M(atom=atoms, basis=BASIS, charge=charge, spin=spin2, verbose=0,
                max_memory=WORKER_MAX_MEMORY_MB)
    mf = dft.UKS(mol).density_fit(auxbasis="def2-universal-jkfit")
    mf.xc = xc
    mf.grids.level = 1
    mf.init_guess = "atom"
    mf.damp = 0.5
    mf.level_shift = 0.5
    mf.diis_start_cycle = 8
    mf.max_cycle = 40
    mf.conv_tol = 1e-5
    try:
        mf.kernel()
        nt = mf.newton()
        nt.conv_tol = 1e-6
        nt.max_cycle = 30
        e = nt.kernel(mf.mo_coeff, mf.mo_occ)
        return float(e), bool(nt.converged), float(nt.spin_square()[0])
    except Exception:  # noqa: BLE001
        return None, False, None


def work(task):
    t0 = time.time()
    e, ok, s2 = scf(task["atoms"], task["charge"], task["spin2"], task["xc"])
    out = {k: v for k, v in task.items() if k != "atoms"}
    out.update(E=e, conv=ok, s2=s2, seconds=round(time.time() - t0))
    return out


def tasks(xc, domes, distances):
    # The porphine core is a dianion; with Fe(II) and neutral NH3 and O2 the
    # site is neutral.
    out = [dict(kind="o2", xc=xc, h=0.0, d=0.0, charge=0, spin2=2,
                atoms=[("O", (0, 0, 0)), ("O", (0, 0, O2_BOND))])]
    for h in domes:
        atoms, fe_z = site(h)
        for spin2 in (0, 2, 4):
            out.append(dict(kind="site", xc=xc, h=h, d=0.0, charge=0, spin2=spin2, atoms=atoms))
        for d in distances:
            for spin2 in (0, 2):
                out.append(dict(kind="complex", xc=xc, h=h, d=d, charge=0, spin2=spin2,
                                atoms=atoms + o2(fe_z, d)))
    out.sort(key=lambda t: -len(t["atoms"]))
    return out


def pull_empty_tasks(xc, axials):
    """The empty flat site with the axial NH3 pulled away (M28-19)."""
    out = [dict(kind="nh3", xc=xc, h=0.0, d=0.0, axial=0.0, charge=0, spin2=0,
                atoms=nh3_below(0.0, 0.0))]
    for axial in axials:
        atoms, _ = site(0.0, axial)
        for spin2 in (0, 2, 4):
            out.append(dict(kind="pull_empty", xc=xc, h=0.0, d=0.0, axial=axial, charge=0,
                            spin2=spin2, atoms=atoms))
    out.sort(key=lambda t: -len(t["atoms"]))
    return out


def pull_tasks(xc, h, d, axials):
    """The loaded site with the axial NH3 pulled away (the cam lever test)."""
    out = [dict(kind="nh3", xc=xc, h=h, d=d, axial=0.0, charge=0, spin2=0,
                atoms=nh3_below(0.0, 0.0))]
    for axial in axials:
        atoms, fe_z = site(h, axial)
        for spin2 in (0, 2):
            out.append(dict(kind="pull", xc=xc, h=h, d=d, axial=axial, charge=0, spin2=spin2,
                            atoms=atoms + o2(fe_z, d)))
    out.sort(key=lambda t: -len(t["atoms"]))
    return out


def core_tasks(xc, cores, domes, d):
    """M28-18: ring size against dome, empty (S=1, 2) and loaded (S=0, 1)."""
    out = [dict(kind="o2", xc=xc, core=0.0, h=0.0, d=0.0, charge=0, spin2=2,
                atoms=[("O", (0, 0, 0)), ("O", (0, 0, O2_BOND))])]
    for core in cores:
        for h in domes:
            atoms, fe_z = site(h, AXIAL_N, core)
            for spin2 in (2, 4):
                out.append(dict(kind="site", xc=xc, core=core, h=h, d=0.0, charge=0, spin2=spin2,
                                atoms=atoms))
            for spin2 in (0, 2):
                out.append(dict(kind="complex", xc=xc, core=core, h=h, d=d, charge=0,
                                spin2=spin2, atoms=atoms + o2(fe_z, d)))
    out.sort(key=lambda t: -len(t["atoms"]))
    return out


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--workers", type=int, default=3)
    p.add_argument("--xc", default="pbe")
    p.add_argument("--domes", nargs="+", type=float, default=[0.0, 0.3, 0.5])
    p.add_argument("--distances", nargs="+", type=float, default=[1.80, 2.00])
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--only", type=int, default=None, help="run only the first N tasks (timing)")
    p.add_argument("--pull", nargs="+", type=float, default=None,
                   help="axial Fe-N distances for the lever test at dome 0 and O2 at 1.85 A")
    p.add_argument("--cores", nargs="+", type=float, default=None,
                   help="Fe-N core sizes for the ring-squeeze test (M28-18)")
    p.add_argument("--empty", action="store_true", help="with --pull: pull the empty site")
    a = p.parse_args()
    if a.cores:
        job = core_tasks(a.xc, a.cores, a.domes, 1.85)
    elif a.pull and a.empty:
        job = pull_empty_tasks(a.xc, a.pull)
    elif a.pull:
        job = pull_tasks(a.xc, 0.0, 1.85, a.pull)
    else:
        job = tasks(a.xc, a.domes, a.distances)
    if a.only:
        job = job[: a.only]
    total = physical_memory_mb()
    cap = MEMORY_FRACTION * total
    # A hybrid functional peaked at 2.0 GB for one process and 5.8 GB for
    # three (2026-09-28), close to the cap, so plan 2.8 GB for a hybrid.
    planned_mb = WORKER_PLANNED_MB if a.xc == "pbe" else 2800
    workers = max(1, min(a.workers, int(cap // planned_mb)))
    print(json.dumps({"plan": {"physical_memory_mb": round(total), "cap_mb": round(cap),
                               "workers": workers, "tasks": len(job), "xc": a.xc}}),
          file=sys.stderr, flush=True)
    if a.dry_run:
        return 0
    stop = threading.Event()
    threading.Thread(target=watchdog, args=(cap, stop), daemon=True).start()
    with mp.Pool(workers, maxtasksperchild=1) as pool:
        for r in pool.imap_unordered(work, job):
            print(json.dumps(r), flush=True)
    stop.set()
    return 0


if __name__ == "__main__":
    sys.exit(main())
