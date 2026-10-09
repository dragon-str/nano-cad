"""M29-03: longer, repeated same-pocket MD at the top stage (600-1000 atm).

`pump_balance.py` fits the free O2 a released pocket carries back as
n_c = s * rho. That fit used single 0.6 ns runs, so the top stage of the
cascade (600 -> 1000 atm) stayed uncertain. This driver runs
`pocket_carry_md.py` for every pocket of the O2 rotor, at 600 and 1000 atm,
with a longer production and several seeds, and writes one JSON line per run.
`pump_balance.py` then consumes the lines.

The O2 density comes from the barostat equation of state of ADR-0123. The
start lattice uses a 0.29 nm spacing so the far bulk reaches 1000 atm.

    python scripts/top_stage_sweep.py --scene ~/.cache/nano-cad/md/module.json \
        --out ~/.cache/nano-cad/md/top_stage.jsonl

Memory safety (AGENTS.md). Each `pocket_carry_md.py` process stays below
0.5 GB. The worker count is capped so the planned total stays under 40 percent
of physical memory, and a watchdog stops every child when the summed resident
memory passes the cap. `--dry-run` prints the plan and starts nothing.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EOS = [(1.0, 0.0237), (100.0, 3.12), (300.0, 6.82), (1000.0, 14.07)]
WORKER_MEMORY_MB = 500
MEMORY_FRACTION = 0.40
WATCHDOG_PERIOD_S = 5.0


def rho(p_atm: float) -> float:
    ps, rs = zip(*EOS)
    return float(np.interp(p_atm, ps, rs))


def physical_memory_mb() -> float:
    out = subprocess.run(["sysctl", "-n", "hw.memsize"], capture_output=True, text=True).stdout.strip()
    return int(out) / (1024.0 * 1024.0)


def child_rss_mb() -> float:
    out = subprocess.run(["ps", "-axo", "rss,command"], capture_output=True, text=True).stdout
    total = 0.0
    for line in out.splitlines():
        if "pocket_carry_md.py" in line:
            total += int(line.split()[0])
    return total / 1024.0


def watchdog(cap_mb: float, stop: threading.Event) -> None:
    while not stop.wait(WATCHDOG_PERIOD_S):
        if child_rss_mb() > cap_mb:
            print(f"watchdog: resident memory over {cap_mb:.0f} MB, stopping", file=sys.stderr, flush=True)
            # pkill is unreliable on this macOS host; kill by pid.
            out = subprocess.run(["ps", "-axo", "pid,command"], capture_output=True, text=True).stdout
            for line in out.splitlines():
                if "pocket_carry_md.py" in line:
                    try:
                        os.kill(int(line.split()[0]), 9)
                    except (ProcessLookupError, ValueError):
                        pass
            return


def run_one(job: dict) -> dict:
    # A NaN start is common at high density; retry with the next seed.
    for attempt in range(3):
        seed = job["seed"] + 100 * attempt
        cmd = [
            sys.executable, os.path.join(HERE, "pocket_carry_md.py"), job["scene"],
            "--pocket", str(job["pocket"]), "--medium", "gas",
            "--o2-per-nm3", f"{job['start_rho']:.3f}", "--o2-spacing", "0.29", "--z-front", str(job["z_front"]),
            "--production-ps", str(job["production_ps"]), "--seed", str(seed),
        ]
        # The module scene keeps a bound O2 in the cam's loaded pockets. A
        # released pocket must drop it; a loaded pocket gets exactly one, at
        # the standard geometry (remove, then add).
        cmd.append("--remove-bound-o2")
        if job["bound"]:
            cmd.append("--add-bound-o2")
        if job.get("prefill"):
            cmd += ["--prefill", str(job["prefill"])]
        else:
            cmd.append("--start-empty")
        if job.get("block_ps"):
            cmd += ["--block-ps", str(job["block_ps"])]
        t0 = time.time()
        proc = subprocess.run(cmd, capture_output=True, text=True)
        line = ""
        for row in proc.stdout.splitlines():
            if row.strip().startswith("{"):
                line = row
        if line:
            rec = json.loads(line)
            rec.update({"ok": True, "atm": job["atm"], "rho_target": job["rho"], "start_rho": job["start_rho"],
                        "z_front": job["z_front"], "seed": seed,
                        "bound": job["bound"], "protocol": "remove-bound-o2", "seconds": round(time.time() - t0)})
            return rec
    return {**job, "seed": seed, "ok": False, "seconds": round(time.time() - t0), "stderr": proc.stderr[-300:]}


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--scene", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--pressures", nargs="+", type=float, default=[600.0, 1000.0])
    p.add_argument("--pockets", nargs="+", type=int, default=list(range(8)))
    p.add_argument("--seeds", nargs="+", type=int, default=[1, 2])
    p.add_argument("--production-ps", type=float, default=2000.0)
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--bound", action="store_true", help="loaded pockets (a frozen O2 on the iron)")
    p.add_argument("--prefill", nargs="+", type=int, default=[0],
                   help="O2 started inside the pocket; 0 and a positive value bracket equilibrium")
    p.add_argument("--block-ps", type=float, default=0.0)
    p.add_argument("--z-front", type=float, default=1.4,
                   help="nm of gas above the rim; a tall box (4.4) keeps the reservoir at the target "
                        "density after the rotor surface adsorbs O2")
    p.add_argument("--start-density", nargs="+", type=float, default=None,
                   help="start O2 per nm3 for each pressure (calibrated so the reservoir reaches the EOS "
                        "density); default: the EOS density itself")
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args()

    starts = a.start_density or [rho(atm) for atm in a.pressures]
    if len(starts) != len(a.pressures):
        p.error("--start-density needs one value per pressure")
    jobs = [dict(scene=a.scene, pocket=k, atm=atm, rho=rho(atm), start_rho=st, z_front=a.z_front, seed=s,
                 production_ps=a.production_ps, bound=a.bound, prefill=f, block_ps=a.block_ps)
            for atm, st in zip(a.pressures, starts) for k in a.pockets for s in a.seeds for f in a.prefill]
    total = physical_memory_mb()
    cap = MEMORY_FRACTION * total
    workers = max(1, min(a.workers, int(cap // WORKER_MEMORY_MB)))
    plan = {"physical_memory_mb": round(total), "cap_mb": round(cap), "workers": workers,
            "runs": len(jobs), "production_ps": a.production_ps, "pressures": a.pressures}
    print(json.dumps({"plan": plan}), file=sys.stderr, flush=True)
    if a.dry_run:
        return 0

    stop = threading.Event()
    threading.Thread(target=watchdog, args=(cap, stop), daemon=True).start()
    done = 0
    with open(a.out, "a") as handle, ThreadPoolExecutor(max_workers=workers) as pool:
        for rec in pool.map(run_one, jobs):
            done += 1
            handle.write(json.dumps(rec) + "\n")
            handle.flush()
            print(f"[{done}/{len(jobs)}] pocket {rec.get('pocket')} {rec.get('atm')}atm seed {rec.get('seed')} "
                  f"prefill {rec.get('prefill')} blocks {rec.get('o2_in_pocket_blocks')} "
                  f"-> o2 {rec.get('o2_in_pocket')} far {rec.get('bulk_o2_per_nm3_far')} ({rec.get('seconds')}s)",
                  file=sys.stderr, flush=True)
    stop.set()
    return 0


if __name__ == "__main__":
    sys.exit(main())
