"""M29-10: does the distal guard lower the autoxidation of the carrier heme?

The carrier site (ADR-0129) has a distal C-O-H group on the pocket wall. Its
oxygen hydrogen-bonds the bound O2 (O...O 2.99 A, H...O 2.03 A). In
hemoglobin the distal histidine does the same, and the hydrogen bond is one
reason the bound O2 resists autoxidation.

Model. The idealized Fe-porphine of `heme_switch_scan.py` (D4h, Fe-N 2.00 A)
with an axial NH3 below the iron, O2 on top end on and bent 60 degrees, and a
water molecule as the hydrogen-bond donor in place of the guard. The water is
placed at the measured guard geometry: its oxygen 2.99 A from the terminal O2
oxygen, with one O-H bond pointing at that oxygen. The water is a proxy for
the C-O-H guard; it has the same O-H donor and fewer atoms.

For each state the scan reports the energy, the Mulliken charge on the O2
atoms, and the Mulliken spin on the O2 atoms. The spin is the superoxide
character of the bound O2: a lower spin on O2 is a less activated O2, and a
lower autoxidation rate.

Derived numbers:
    E_bind          = E(complex) - E(site) - E(O2)            (min over spins)
    E_bind_guard    = E(complex+water) - E(site) - E(O2) - E(water)
    H_bond          = E(complex+water) - E(complex) - E(water)

Method: UKS with density fitting (def2-universal-jkfit), def2-SVP on Fe, N
and O, 6-31G on ring carbons, STO-3G on hydrogens, grid level 1, PBE, gas
phase (the pocket is hydrophobic). Rigid, idealized geometry. Every energy is
an estimate for a proxy, not the diamond site. Needs pyscf (Apache-2.0).

Memory safety as in `fe_o2_switch_scan.py`: 1 GB PySCF budget per worker, the
worker count capped at 40 percent of physical memory, a watchdog, one task
per worker process, `--dry-run`.

    python scripts/heme_guard_scan.py > guard.jsonl
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
from heme_switch_scan import BASIS, O2_BOND, o2, site  # noqa: E402

# The measured guard geometry (ADR-0129): the guard oxygen 2.99 A from the
# terminal O2 oxygen. The unit direction comes from the CAD part.
GUARD_O_O = 2.99
GUARD_DIR = np.array([0.181, 0.598, 0.781])
GUARD_DIR /= np.linalg.norm(GUARD_DIR)
OH_BOND = 0.96
HOH_ANGLE_DEG = 104.5


def water(toward: np.ndarray) -> list:
    """A water molecule whose first O-H points at `toward`."""
    o = toward + GUARD_O_O * GUARD_DIR
    u = toward - o
    u = u / np.linalg.norm(u)
    perp = np.cross(u, [0.0, 0.0, 1.0])
    if np.linalg.norm(perp) < 1.0e-6:
        perp = np.cross(u, [1.0, 0.0, 0.0])
    perp = perp / np.linalg.norm(perp)
    ang = np.radians(HOH_ANGLE_DEG)
    h1 = o + OH_BOND * u
    h2 = o + OH_BOND * (np.cos(ang) * u + np.sin(ang) * perp)
    return [("O", tuple(o)), ("H", tuple(h1)), ("H", tuple(h2))]


def scf(atoms, charge, spin2, xc):
    """UKS with density fitting: a damped first pass, then a second-order
    finish. Returns the energy and the per-atom Mulliken charge and spin."""
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
        e = float(nt.kernel(mf.mo_coeff, mf.mo_occ))
        ok = bool(nt.converged)
        s2 = float(nt.spin_square()[0])
        dm = nt.make_rdm1()
        s = nt.get_ovlp()
        ps = dm @ s
        slices = mol.aoslice_by_atom()
        if ps.ndim == 2:
            pop = np.array([np.trace(ps[a:b, a:b]) for (_, _, a, b) in slices])
            spin = np.zeros_like(pop)
        else:
            pa, pb = ps[0], ps[1]
            pop = np.array([np.trace(pa[a:b, a:b]) + np.trace(pb[a:b, a:b]) for (_, _, a, b) in slices])
            spin = np.array([np.trace(pa[a:b, a:b]) - np.trace(pb[a:b, a:b]) for (_, _, a, b) in slices])
        charge_atom = mol.atom_charges() - pop
        return e, ok, s2, charge_atom, spin
    except Exception:  # noqa: BLE001
        return None, False, None, None, None


def work(task):
    t0 = time.time()
    e, ok, s2, charge, spin = scf(task["atoms"], task["charge"], task["spin2"], task["xc"])
    o2_idx = task["o2_idx"]
    out = {k: v for k, v in task.items() if k not in ("atoms", "o2_idx")}
    out.update(
        E=e,
        conv=ok,
        s2=s2,
        o2_charge=None if charge is None else round(float(charge[o2_idx].sum()), 3),
        o2_spin=None if spin is None else round(float(spin[o2_idx].sum()), 3),
        seconds=round(time.time() - t0),
    )
    return out


def tasks(xc, d):
    atoms_site, fe_z = site(0.0)
    oxy = o2(fe_z, d)
    complex_atoms = atoms_site + oxy
    o2_idx = list(range(len(atoms_site), len(complex_atoms)))
    terminal = np.array(oxy[1][1])
    w = water(terminal)
    guard_atoms = complex_atoms + w
    guard_o2_idx = list(range(len(atoms_site), len(atoms_site) + 2))
    out = [
        dict(kind="o2", xc=xc, charge=0, spin2=2, d=0.0,
             atoms=[("O", (0.0, 0.0, 0.0)), ("O", (0.0, 0.0, O2_BOND))], o2_idx=[0, 1]),
        dict(kind="water", xc=xc, charge=0, spin2=0, d=0.0, atoms=w, o2_idx=[]),
    ]
    for spin2 in (0, 2, 4):
        out.append(dict(kind="site", xc=xc, charge=0, spin2=spin2, d=0.0, atoms=atoms_site, o2_idx=[]))
    for spin2 in (0, 2):
        out.append(dict(kind="complex", xc=xc, charge=0, spin2=spin2, d=d, atoms=complex_atoms, o2_idx=o2_idx))
    for spin2 in (0, 2):
        out.append(dict(kind="complex_guard", xc=xc, charge=0, spin2=spin2, d=d, atoms=guard_atoms,
                        o2_idx=guard_o2_idx))
    site_water_atoms = atoms_site + w
    for spin2 in (0, 2):
        out.append(dict(kind="site_water", xc=xc, charge=0, spin2=spin2, d=0.0, atoms=site_water_atoms,
                        o2_idx=[]))
    out.sort(key=lambda t: -len(t["atoms"]))
    return out


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--workers", type=int, default=3)
    p.add_argument("--xc", default="pbe")
    p.add_argument("--d", type=float, default=1.85, help="Fe-O2 distance in Angstrom")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--only", type=int, default=None, help="run only the first N tasks (timing)")
    p.add_argument("--kinds", default=None, help="comma-separated kinds to run")
    a = p.parse_args()
    job = tasks(a.xc, a.d)
    if a.only:
        job = job[: a.only]
    if a.kinds:
        keep = set(a.kinds.split(","))
        job = [t for t in job if t["kind"] in keep]
    total = physical_memory_mb()
    cap = MEMORY_FRACTION * total
    workers = max(1, min(a.workers, int(cap // WORKER_PLANNED_MB)))
    print(json.dumps({"plan": {"physical_memory_mb": round(total), "cap_mb": round(cap),
                               "workers": workers, "tasks": len(job), "xc": a.xc, "d": a.d}}),
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
