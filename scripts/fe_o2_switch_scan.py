"""Option B kill test: can strain on the iron site switch off O2 binding?

A cam cannot pull the oxygen off the iron, but it can move the atoms that
hold the iron. The scan asks whether any such strain weakens the Fe-O2 bond by
the 45-50 kJ/mol that a chemical release switch needs (ADR-0092 context).

Model: the Fe(NH3)4 proxy of ADR-0084, UKS PBE/def2-SVP, grid level 1, no
dispersion, rigid fragments. Three strain modes:

* dome h: the iron moves h Angstrom out of the N4 plane, away from O2
  (the hemoglobin T-state motion);
* stretch r: the four Fe-N distances become r Angstrom;
* trans t: an NH3 is pushed onto the back of the iron to t Angstrom
  (trans influence). The unstrained reference is no trans ligand.

For each state: E_bind = min over d of E(frame+O2, S=1 or 2) - min over S of
E(frame) - E(O2). Every energy is an estimate for a proxy, not the diamond
site. Needs pyscf (Apache-2.0) in its own environment.

    python scripts/fe_o2_switch_scan.py > switch.jsonl

Memory safety. The run on 2026-09-25 started 10 workers with a 3 GB PySCF
budget each on a 16 GB machine, and the machine ran out of memory and
crashed. Now each worker has a 1 GB PySCF budget, the worker count is capped
so that the planned total stays under 40 percent of physical memory, the
second-order (Newton) fallback is gone because it needs more memory, and a
watchdog kills every worker when their summed resident memory passes the
cap. `--dry-run` prints the plan and starts nothing.
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

import argparse
import json
import multiprocessing as mp
import signal
import subprocess
import sys
import threading
import time

import numpy as np

HARTREE_KJ = 2625.5
# PySCF memory budget of one worker, in MB. PySCF treats it as a soft limit,
# so the planner assumes each worker can use twice this much.
WORKER_MAX_MEMORY_MB = 1000
WORKER_PLANNED_MB = 2 * WORKER_MAX_MEMORY_MB
# The fraction of physical memory that all workers together may use.
MEMORY_FRACTION = 0.40
WATCHDOG_PERIOD_S = 2.0
O2_BOND_A = 1.208
O2_TILT_DEG = 60.0


def nh3(center, toward, d_nh=1.01):
    n = np.array(center, float)
    t = np.array(toward, float)
    t /= np.linalg.norm(t)
    a = np.array([0.0, 0.0, 1.0]) - t * t[2]
    if np.linalg.norm(a) < 1e-6:
        a = np.array([1.0, 0.0, 0.0]) - t * t[0]
    a /= np.linalg.norm(a)
    b = np.cross(t, a)
    th = np.radians(112.0)
    out = [("N", tuple(n))]
    for k in range(3):
        ang = np.radians(120.0 * k)
        v = t * np.cos(th) + (np.cos(ang) * a + np.sin(ang) * b) * np.sin(th)
        out.append(("H", tuple(n + d_nh * v)))
    return out


def frame(mode: str, value: float):
    """The iron with its ligands. The iron sits at (0, 0, fe_z)."""
    r = value if mode == "stretch" else 2.15
    fe_z = -value if mode == "dome" else 0.0
    atoms = [("Fe", (0.0, 0.0, fe_z))]
    for d in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)):
        v = np.array(d, float)
        atoms += nh3(v * r, v)
    if mode == "trans":
        atoms += nh3((0.0, 0.0, -value), (0.0, 0.0, -1.0))
    return atoms, fe_z


def o2(fe_z: float, d: float):
    t = np.radians(O2_TILT_DEG)
    z = fe_z + d
    return [
        ("O", (0.0, 0.0, z)),
        ("O", (float(O2_BOND_A * np.sin(t)), 0.0, float(z + O2_BOND_A * np.cos(t)))),
    ]


def scf(atoms, spin2):
    from pyscf import dft, gto

    mol = gto.M(atom=atoms, basis="def2-svp", charge=0, spin=spin2, verbose=0, max_memory=WORKER_MAX_MEMORY_MB)
    mf = dft.UKS(mol)
    mf.xc = "pbe"
    mf.grids.level = 1
    mf.conv_tol = 1e-6
    mf.max_cycle = 100
    mf.level_shift = 0.2
    try:
        e = mf.kernel()
        if not mf.converged:
            mf.level_shift = 0.5
            mf.max_cycle = 200
            e = mf.kernel(mf.make_rdm1())
        return float(e), bool(mf.converged)
    except Exception:  # noqa: BLE001
        return None, False


def work(task):
    t0 = time.time()
    e, ok = scf(task["atoms"], task["spin2"])
    out = {k: v for k, v in task.items() if k != "atoms"}
    out.update(E=e, conv=ok, seconds=round(time.time() - t0))
    return out


STATES = [
    ("base", 0.0),
    ("dome", 0.3),
    ("dome", 0.6),
    ("stretch", 2.35),
    ("stretch", 2.55),
    ("trans", 2.10),
    ("trans", 1.95),
]
DISTANCES = (1.85, 2.05, 2.25)


def tasks():
    out = [dict(kind="o2", mode="-", value=0.0, d=0.0, spin2=2, atoms=o2(0.0, 5.0)[:2])]
    for mode, value in STATES:
        atoms, fe_z = frame(mode, value)
        for spin2 in (0, 2, 4):
            out.append(dict(kind="frame", mode=mode, value=value, d=0.0, spin2=spin2, atoms=atoms))
        for d in DISTANCES:
            out.append(dict(kind="complex", mode=mode, value=value, d=d, spin2=2, atoms=atoms + o2(fe_z, d)))
        out.append(dict(kind="complex", mode=mode, value=value, d=2.05, spin2=4, atoms=atoms + o2(fe_z, 2.05)))
    # Heavy jobs first, so the pool drains evenly.
    out.sort(key=lambda t: -len(t["atoms"]))
    return out


def physical_memory_mb() -> float:
    out = subprocess.run(["sysctl", "-n", "hw.memsize"], capture_output=True, text=True)
    if out.returncode == 0 and out.stdout.strip():
        return int(out.stdout.strip()) / 2**20
    return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 2**20


def child_rss_mb(parent: int) -> tuple[float, list[int]]:
    """The summed resident memory of all descendants of `parent`, in MB."""
    out = subprocess.run(["ps", "-A", "-o", "pid=,ppid=,rss="], capture_output=True, text=True)
    children: dict[int, list[tuple[int, int]]] = {}
    for line in out.stdout.splitlines():
        parts = line.split()
        if len(parts) == 3:
            pid, ppid, rss = (int(x) for x in parts)
            children.setdefault(ppid, []).append((pid, rss))
    total_kb, pids, stack = 0, [], [parent]
    while stack:
        for pid, rss in children.get(stack.pop(), []):
            total_kb += rss
            pids.append(pid)
            stack.append(pid)
    return total_kb / 1024, pids


def watchdog(cap_mb: float, stop: threading.Event) -> None:
    me = os.getpid()
    while not stop.wait(WATCHDOG_PERIOD_S):
        used_mb, pids = child_rss_mb(me)
        if used_mb > cap_mb:
            print(
                f"watchdog: workers use {used_mb:.0f} MB > cap {cap_mb:.0f} MB; stopping",
                file=sys.stderr,
                flush=True,
            )
            for pid in pids:
                try:
                    os.kill(pid, signal.SIGKILL)
                except OSError:
                    pass
            os._exit(3)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    total_mb = physical_memory_mb()
    cap_mb = MEMORY_FRACTION * total_mb
    workers = max(1, min(args.workers, int(cap_mb // WORKER_PLANNED_MB)))
    plan = {
        "physical_memory_mb": round(total_mb),
        "cap_mb": round(cap_mb),
        "requested_workers": args.workers,
        "workers": workers,
        "worker_max_memory_mb": WORKER_MAX_MEMORY_MB,
        "tasks": len(tasks()),
    }
    print(json.dumps({"plan": plan}), file=sys.stderr, flush=True)
    if args.dry_run:
        return 0
    stop = threading.Event()
    guard = threading.Thread(target=watchdog, args=(cap_mb, stop), daemon=True)
    guard.start()
    # One task per worker process, so memory returns to the system each time.
    with mp.Pool(workers, maxtasksperchild=1) as pool:
        for result in pool.imap_unordered(work, tasks()):
            print(json.dumps(result), flush=True)
    stop.set()
    return 0


if __name__ == "__main__":
    sys.exit(main())
