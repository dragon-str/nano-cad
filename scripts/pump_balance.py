"""O2 and water balance of the respirocyte rotor cascade (ADR-0123).

Inputs are the JSON lines that `pocket_carry_md.py` prints:

  --clearance FILES   O2-rotor returning pockets in O2 gas (medium gas): the
                      free O2 a released pocket carries back, against the far
                      bulk density.
  --loaded FILES      O2-rotor loaded pockets (O2 bound) in O2 gas: the free
                      O2 a loaded pocket carries forward on a gas-to-gas
                      stage. Without them the forward term is zero.
  --water-sump FILE   water-rotor pocket withdrawn in tank-sump water.
  --water-blood FILE  water-rotor pocket flush in blood water.
  --o2-blood FILE     O2-rotor loaded pocket in blood water (water carried in).

The clearance n_c is fitted as n_c = s * rho through the origin. One pocket
pass from a low side to a high side moves

    net = 1 (the bound O2) + n_l(rho_low) - n_c(rho_high)

and a liquid low side (blood) carries no free O2 (n_c(low) = 0). The O2
equation of state of this force field (barostat MD, 310 K) maps pressure to
density; it is interpolated linearly between the measured points.

    python scripts/pump_balance.py --clearance nc_*.txt ... --water-sump w_sump.txt
"""

from __future__ import annotations

import argparse
import json

import numpy as np

# Barostat MD of the two-site O2 model at 310 K (ADR-0123): atm -> O2/nm^3.
EOS = [(1.0, 0.0237), (100.0, 3.12), (300.0, 6.82), (1000.0, 14.07)]


def rho(p_atm: float) -> float:
    ps, rs = zip(*EOS)
    return float(np.interp(p_atm, ps, rs))


def first_json(path: str) -> dict:
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith("{"):
                return json.loads(line)
    raise SystemExit(f"no result in {path}")


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--clearance", nargs="+", required=True)
    p.add_argument("--loaded", nargs="+", default=[])
    p.add_argument("--water-sump")
    p.add_argument("--water-blood")
    p.add_argument("--o2-blood")
    p.add_argument("--stages", nargs="+", type=float, default=[100.0, 300.0, 600.0, 1000.0],
                   help="the high-side pressure of each stage, atm (the first stage loads from blood)")
    a = p.parse_args()

    def fit(paths: list[str], name: str) -> tuple[float, float]:
        pts = [first_json(path) for path in paths]
        x = np.array([d["bulk_o2_per_nm3_far"] for d in pts])
        y = np.array([d["o2_in_pocket"] for d in pts])
        slope = float((x @ y) / (x @ x))
        resid = y - slope * x
        err = float(np.sqrt((resid @ resid) / max(1, len(x) - 1) / (x @ x)))
        print(f"{name} fit: n = {slope:.4f} * rho (+- {err:.4f}), {len(x)} runs, "
              f"rms residual {np.sqrt(np.mean(resid ** 2)):.2f}")
        return slope, err

    s, s_err = fit(a.clearance, "clearance (released pocket)")
    sl, sl_err = fit(a.loaded, "forward (loaded pocket)") if a.loaded else (0.0, 0.0)
    single = 1.0 / s
    p_stall = float(np.interp(single, [r for _, r in EOS], [q for q, _ in EOS]))
    print(f"single stage from blood stalls at rho = {single:.2f} /nm^3, about {p_stall:.0f} atm")

    print("cascade (net O2 per loaded pocket pass):")
    low = 0.0
    nets = []
    for k, p_hi in enumerate(a.stages):
        n_low = 0.0 if k == 0 else sl * rho(low)
        n_high = s * rho(p_hi)
        net = 1.0 + n_low - n_high
        net_lo = 1.0 + (0.0 if k == 0 else (sl - sl_err) * rho(low)) - (s + s_err) * rho(p_hi)
        nets.append(net)
        print(f"  stage {k + 1}: {'blood' if k == 0 else f'{low:.0f} atm'} -> {p_hi:.0f} atm: "
              f"carried back {n_high:.2f}, forward {n_low:.2f}, net {net:.2f} (low bound {net_lo:.2f})")
        low = p_hi
    print(f"  pocket passes per O2 stored, summed over stages: {sum(1.0 / n for n in nets):.1f}")

    if a.o2_blood and a.water_sump:
        w_in = first_json(a.o2_blood)["water_in_pocket"]
        sump = first_json(a.water_sump)
        back = first_json(a.water_blood)["water_in_pocket"] if a.water_blood else 0.0
        w_out = sump["water_in_pocket"] - back
        o2_out = sump["o2_in_pocket"]
        net1 = nets[0]
        per_o2 = w_in / net1
        print("water (stage 1 delivers into the wet buffer):")
        print(f"  carried in {w_in:.2f} per loaded pocket pass, {per_o2:.2f} per O2 net")
        print(f"  water rotor: {sump['water_in_pocket']:.2f} out, {back:.2f} back, net {w_out:.2f} per pocket pass")
        print(f"  water-rotor passes per O2 net: {per_o2 / w_out:.3f}")
        print(f"  O2 in a sump pocket: {o2_out:.3f} (max {sump['o2_in_pocket_max']}) "
              f"-> O2 returned to blood per O2 net {o2_out * per_o2 / w_out:.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
