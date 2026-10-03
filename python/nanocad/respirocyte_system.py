"""Respirocyte system model (task M29-01): a real-size device built from
components whose parameters come from the atomistic work.

This is the L4 layer for the respirocyte. Nothing here is atomistic. Each
component carries parameters, and every parameter carries its provenance:
a value, a low-high range, a unit, a tag (``measured`` from a simulation or
CAD build in this repository, ``estimated`` from a calculation on such
results, ``assumed`` with no calculation behind it, ``reference`` from the
literature), and a source (an ADR or a citation). Results are simulated
estimates. Nothing here is built or validated.

The device
----------
A sphere of hull diameter ``hull.diameter_m`` (Freitas: 1 um). Inside it:

* an O2 tank (``tank``), the top of the cascade;
* ``n_stages - 1`` intermediate gas chambers between cascade stages;
* rotor banks. A bank is a count of identical rotor sections:

  - ``load``: the O2 rotor of ADR-0121 in the hull wall, blood to chamber 1;
  - ``cascade``: the same rotor on an internal bulkhead, chamber k-1 to k;
  - ``water``: the plunger rotor of ADR-0122, chamber 1 sump to blood;
  - ``discharge``: the O2 rotor with the cam reversed, tank to blood.

* a power supply with a peak output; rotor speeds are scaled down together
  when the demand exceeds it (the duty cycle).

Pocket balance (ADR-0123)
-------------------------
One pocket pass of an O2 rotor from a low side to a high side moves

    net = occupancy + s_fwd * rho_low - s_back * rho_high

O2 molecules, where rho is the free-gas number density on each side
(molecules per nm^3), ``s_back`` the free gas a released pocket carries
back per unit density and ``s_fwd`` the free gas a loaded pocket carries
forward. A liquid low side (blood) has rho_low = 0. Occupancy is the
fraction of pockets that bind O2 before leaving the low side: from blood it
is limited by the binding rate (ADR-0109); from gas it is taken as 1.
Density follows the O2 equation of state of the atomistic force field
(barostat MD, ADR-0123).

Usage
-----
    from nanocad import respirocyte_system as rs
    device = rs.default_device()
    print(rs.static_budget(device))
    course = rs.time_course(device)
    python -m nanocad.respirocyte_system --json report.json
"""

from __future__ import annotations

import argparse
import copy
import json
import math
from dataclasses import asdict, dataclass, field

AVOGADRO = 6.02214076e23
KJ_PER_MOL_TO_J = 1.0e3 / AVOGADRO
NM3 = 1.0e-27
ATM_PA = 101325.0

TAGS = ("measured", "estimated", "assumed", "reference")


@dataclass
class Param:
    """One model parameter with its provenance."""

    value: float
    low: float
    high: float
    unit: str
    tag: str
    source: str
    note: str = ""

    def __post_init__(self) -> None:
        if self.tag not in TAGS:
            raise ValueError(f"unknown tag {self.tag!r}")
        if not self.low <= self.value <= self.high:
            raise ValueError(f"value {self.value} outside [{self.low}, {self.high}]")


def P(value, low=None, high=None, unit="", tag="assumed", source="", note=""):  # noqa: N802
    return Param(float(value), float(value if low is None else low),
                 float(value if high is None else high), unit, tag, source, note)


# The O2 equation of state of the two-site O2 model at 310 K (ADR-0123):
# pressure (atm) -> number density (molecules per nm^3). Above 1000 atm the
# last segment is extended (an extrapolation).
EOS_ATM_RHO = ((0.0, 0.0), (1.0, 0.0237), (100.0, 3.12), (300.0, 6.82), (1000.0, 14.07))


def rho_of_p(p_atm: float) -> float:
    pts = EOS_ATM_RHO
    for (p0, r0), (p1, r1) in zip(pts, pts[1:]):
        if p_atm <= p1:
            return r0 + (r1 - r0) * (max(p_atm, 0.0) - p0) / (p1 - p0)
    (p0, r0), (p1, r1) = pts[-2], pts[-1]
    return r1 + (r1 - r0) * (p_atm - p1) / (p1 - p0)


def p_of_rho(rho: float) -> float:
    pts = EOS_ATM_RHO
    for (p0, r0), (p1, r1) in zip(pts, pts[1:]):
        if rho <= r1:
            return p0 + (p1 - p0) * (max(rho, 0.0) - r0) / (r1 - r0)
    (p0, r0), (p1, r1) = pts[-2], pts[-1]
    return p1 + (p1 - p0) * (rho - r1) / (r1 - r0)


def section_parameters() -> dict[str, Param]:
    """The rotor section (rotor in its hull-wall slab), from the CAD build."""
    return {
        "footprint_y_m": P(13.7e-9, unit="m", tag="measured", source="ADR-0121 relaxed scene",
                           note="slab extent along the wall, normal to the rotor axis"),
        "footprint_z_m": P(6.9e-9, unit="m", tag="measured", source="ADR-0121 relaxed scene",
                           note="slab extent along the rotor axis"),
        "depth_m": P(10.5e-9, unit="m", tag="measured", source="ADR-0121 relaxed scene",
                     note="end plates reach 5.2 nm into each side"),
        "atoms": P(80620, unit="", tag="measured", source="ADR-0121 relaxed scene"),
    }


@dataclass
class Bank:
    """A bank of identical rotor sections."""

    kind: str  # load, cascade, water, discharge
    count: int
    rev_per_s: float  # design speed, before the power cap
    stage: int = 1  # cascade stage this bank serves (load is stage 1)
    placement: str = "hull"  # hull or internal
    params: dict[str, Param] = field(default_factory=dict)


def o2_rotor_parameters() -> dict[str, Param]:
    return {
        "pockets": P(8, unit="", tag="measured", source="ADR-0121"),
        "binding_rate_per_s": P(2.0e3, 5.0e2, 5.0e3, "1/s", "estimated", "ADR-0109",
                                "k_on c per site in blood; sets loading occupancy"),
        "blood_arc_fraction": P(0.41, 0.35, 0.45, "", "measured", "ADR-0118, ADR-0121",
                                "fraction of a turn a pocket spends open to blood (seal arcs excluded)"),
        "s_back_nm3": P(0.109, 0.076, 0.142, "nm^3", "estimated", "ADR-0123, ADR-0124",
                        "free O2 carried back by a released pocket per unit gas density; fit of 16 runs, "
                        "range +-1 sd of the per-run ratio (pocket-to-pocket scatter)"),
        "s_fwd_nm3": P(0.081, 0.033, 0.141, "nm^3", "estimated", "ADR-0123, ADR-0124",
                       "free O2 carried forward by a loaded pocket per unit gas density; fit of 11 runs, "
                       "range +-1 sd"),
        "gas_occupancy": P(1.0, 0.9, 1.0, "", "assumed", "",
                           "largest fraction of pockets that bind from gas; not computed"),
        "gas_binding_per_s_per_rho": P(1.0e5, 1.0e4, 1.0e6, "nm^3/s", "assumed", "",
                                       "binding rate from gas per unit density: kinetic arrival "
                                       "(about 1e7 per s per molecule/nm^3) times a 1 percent capture"),
        "water_per_pass": P(0.41, 0.2, 1.4, "", "measured", "ADR-0121",
                            "waters a loaded pocket carries from blood into chamber 1 (1 ns MD)"),
        "energy_per_pass_kj_mol": P(50.0, 23.0, 63.0, "kJ/mol", "estimated", "ADR-0100, ADR-0106, ADR-0107",
                                    "cam work of one arm withdrawal; 23 is the reversible minimum"),
        "max_rev_per_s": P(1.0e4, 1.0e3, 1.0e5, "1/s", "assumed", "",
                           "mechanical limit; rim speed 0.25 m/s at 1e4 rev/s"),
    }


def water_rotor_parameters() -> dict[str, Param]:
    return {
        "pockets": P(6, unit="", tag="measured", source="ADR-0122"),
        "water_per_pass": P(5.9, 3.0, 7.4, "", "measured", "ADR-0122",
                            "net waters per pocket pass, sump to blood (1 ns MD)"),
        "o2_per_pass_in_gas_nm3": P(0.163, 0.1, 0.3, "nm^3", "measured", "ADR-0122",
                                    "O2 per pass per unit density if the plunger draws gas (dry sump)"),
        "energy_per_pass_kj_mol": P(5.0, 1.0, 20.0, "kJ/mol", "assumed", "",
                                    "plunger stroke; not computed"),
        "max_rev_per_s": P(1.0e4, 1.0e3, 1.0e5, "1/s", "assumed", ""),
    }


def discharge_rotor_parameters() -> dict[str, Param]:
    p = o2_rotor_parameters()
    p["water_per_pass"] = P(1.0, 0.41, 2.0, "", "assumed", "ADR-0111, ADR-0121",
                            "blood water a released pocket carries into the tank; not computed for the 4 A bore")
    return p


def default_device(n_stages: int = 4) -> dict:
    """A first proof-of-concept device. Bank counts are design choices
    (assumed), sized so the hull holds them."""
    hull = {
        "diameter_m": P(1.0e-6, unit="m", tag="reference", source="Freitas 1998, Respirocytes"),
        "wall_thickness_m": P(2.4e-9, unit="m", tag="measured", source="ADR-0118"),
        "hull_fill_max": P(0.5, 0.2, 0.8, "", "assumed", "",
                           "largest fraction of the hull area rotor sections may take"),
    }
    tank = {
        "volume_m3": P(1.19e-19, unit="m^3", tag="reference", source="Freitas 1998 (ADR-0100)"),
        "target_atm": P(1000.0, 500.0, 1000.0, "atm", "reference", "Freitas 1998"),
        "cushion_atm": P(10.0, 5.0, 20.0, "atm", "estimated", "ADR-0102",
                         "never drained below this, so the gas phase never has to nucleate"),
        "chamber_volume_m3": P(2.0e-21, 5.0e-22, 1.0e-20, "m^3", "assumed", "",
                               "each intermediate cascade chamber"),
    }
    power = {
        "peak_w": P(7.3e-12, 3.0e-13, 2.0e-11, "W", "reference", "Freitas 1998 (ADR-0100)",
                    "all 12 glucose engines; low is one engine"),
    }
    banks = [Bank("load", 10000, 300.0, 1, "hull", o2_rotor_parameters())]
    for k in range(2, n_stages + 1):
        banks.append(Bank("cascade", 6000, 1000.0, k, "internal", o2_rotor_parameters()))
    banks.append(Bank("water", 2000, 300.0, 1, "hull", water_rotor_parameters()))
    banks.append(Bank("discharge", 4000, 1000.0, n_stages, "hull", discharge_rotor_parameters()))
    return {"hull": hull, "tank": tank, "power": power, "section": section_parameters(),
            "n_stages": n_stages, "banks": banks}


def _v(p: Param, case: str, worse_is_high: bool) -> float:
    if case == "central":
        return p.value
    pessimistic = case == "pessimistic"
    return p.high if (pessimistic == worse_is_high) else p.low


def with_case(device: dict, case: str) -> dict:
    """A copy of the device with every uncertain parameter at the end of its
    range that hurts (``pessimistic``) or helps (``optimistic``)."""
    if case == "central":
        return device
    d = copy.deepcopy(device)
    hurts_high = {"s_back_nm3", "energy_per_pass_kj_mol", "water_per_pass", "o2_per_pass_in_gas_nm3"}
    hurts_low = {"binding_rate_per_s", "s_fwd_nm3", "gas_occupancy", "blood_arc_fraction", "peak_w"}
    for b in d["banks"]:
        for name, p in b.params.items():
            if b.kind == "water" and name == "water_per_pass":
                hi = False
            elif name in hurts_high:
                hi = True
            elif name in hurts_low:
                hi = False
            else:
                continue
            p.value = _v(p, case, hi)
    pw = d["power"]["peak_w"]
    pw.value = _v(pw, case, False)
    return d


def hull_area_m2(device: dict) -> float:
    return math.pi * device["hull"]["diameter_m"].value ** 2


def section_area_m2(device: dict) -> float:
    s = device["section"]
    return s["footprint_y_m"].value * s["footprint_z_m"].value


def load_occupancy(params: dict[str, Param], rev_per_s: float) -> float:
    dwell = params["blood_arc_fraction"].value / max(rev_per_s, 1e-12)
    return 1.0 - math.exp(-params["binding_rate_per_s"].value * dwell)


def gas_occupancy(params: dict[str, Param], rho: float, rev_per_s: float) -> float:
    dwell = params["blood_arc_fraction"].value / max(rev_per_s, 1e-12)
    rate = params["gas_binding_per_s_per_rho"].value * max(rho, 0.0)
    return params["gas_occupancy"].value * (1.0 - math.exp(-rate * dwell))


def net_per_pass(params: dict[str, Param], rho_low: float, rho_high: float, occupancy: float) -> float:
    return occupancy + params["s_fwd_nm3"].value * rho_low - params["s_back_nm3"].value * rho_high


def stall_atm(params: dict[str, Param], low_atm: float | None, occupancy: float = 1.0) -> float:
    """The high-side pressure at which a stage stops pumping (net = 0).
    ``low_atm`` None means a liquid low side (blood)."""
    rho_low = 0.0 if low_atm is None else rho_of_p(low_atm)
    rho_hi = (occupancy + params["s_fwd_nm3"].value * rho_low) / params["s_back_nm3"].value
    return p_of_rho(rho_hi)


def static_budget(device: dict) -> dict:
    """Areas, counts, capacities, speeds, energies and the stall pressures."""
    hull_a = hull_area_m2(device)
    sec_a = section_area_m2(device)
    hull_banks = [b for b in device["banks"] if b.placement == "hull"]
    used = sum(b.count for b in hull_banks) * sec_a
    tank_v = device["tank"]["volume_m3"].value
    target = device["tank"]["target_atm"].value
    capacity = rho_of_p(target) * tank_v / NM3
    cushion = rho_of_p(device["tank"]["cushion_atm"].value) * tank_v / NM3
    out = {
        "hull_area_m2": hull_area_m2(device),
        "section_area_m2": sec_a,
        "sections_that_fit": int(hull_a / sec_a),
        "hull_fraction_used": used / hull_a,
        "hull_fraction_allowed": device["hull"]["hull_fill_max"].value,
        "hull_fits": used / hull_a <= device["hull"]["hull_fill_max"].value,
        "tank_capacity_o2_at_target": capacity,
        "tank_cushion_o2": cushion,
        "deliverable_o2": capacity - cushion,
        "banks": [],
    }
    # Stall pressures: each stage with its low side at the stall of the one
    # below (the most a cascade can reach with that many stages).
    low = None
    stalls = []
    for k in range(1, device["n_stages"] + 1):
        b = next(x for x in device["banks"] if x.kind in ("load", "cascade") and x.stage == k)
        occ = load_occupancy(b.params, b.rev_per_s) if k == 1 else b.params["gas_occupancy"].value
        st = stall_atm(b.params, low, occ)
        stalls.append(st)
        low = st
    out["stage_stall_atm"] = stalls
    out["cascade_ceiling_atm"] = stalls[-1]
    # The limit of many close stages near the top: net = 1 - (s_back - s_fwd) rho.
    b1 = device["banks"][0].params
    ds = b1["s_back_nm3"].value - b1["s_fwd_nm3"].value
    out["infinite_cascade_net_at_target"] = 1.0 - ds * rho_of_p(target)
    for b in device["banks"]:
        e = b.params["energy_per_pass_kj_mol"].value * KJ_PER_MOL_TO_J
        passes = b.count * b.params["pockets"].value * b.rev_per_s
        out["banks"].append({
            "kind": b.kind, "stage": b.stage, "count": b.count, "placement": b.placement,
            "rev_per_s": b.rev_per_s, "pocket_passes_per_s": passes,
            "power_at_design_speed_w": passes * e,
        })
    out["power_demand_at_design_w"] = sum(x["power_at_design_speed_w"] for x in out["banks"])
    out["power_peak_w"] = device["power"]["peak_w"].value
    return out


def _rates(device: dict, n_o2: list[float], water: float, phase: str, scale: float):
    """O2 transfer rates between compartments, water rates and power.

    Compartments: 0 blood (infinite), 1..n_stages-1 chambers, n_stages tank.
    Returns (dN per compartment, dWater of chamber 1, power W, per-bank info).
    """
    ns = device["n_stages"]
    vols = [0.0] + [device["tank"]["chamber_volume_m3"].value] * (ns - 1) + [device["tank"]["volume_m3"].value]
    rho = [0.0] + [n_o2[k] * NM3 / vols[k] for k in range(1, ns + 1)]
    dn = [0.0] * (ns + 1)
    dw = 0.0
    power = 0.0
    info = []
    for b in device["banks"]:
        pr = b.params
        rev = min(b.rev_per_s, pr["max_rev_per_s"].value) * scale
        on = (phase == "fill" and b.kind in ("load", "cascade", "water")) or \
             (phase == "discharge" and b.kind in ("discharge", "water"))
        if not on or rev <= 0:
            info.append((b.kind, b.stage, 0.0, 0.0))
            continue
        passes = b.count * pr["pockets"].value * rev
        e = pr["energy_per_pass_kj_mol"].value * KJ_PER_MOL_TO_J
        if b.kind in ("load", "cascade"):
            lo, hi = b.stage - 1, b.stage
            occ = load_occupancy(pr, rev) if lo == 0 else gas_occupancy(pr, rho[lo], rev)
            net = net_per_pass(pr, rho[lo], rho[hi], occ)
            flow = passes * net
            dn[hi] += flow
            if lo > 0:
                dn[lo] -= flow
            if b.kind == "load":
                dw += passes * occ * pr["water_per_pass"].value
        elif b.kind == "discharge":
            # Loads in the tank gas, releases into blood; the free gas a
            # loaded pocket carries leaves too.
            occ = gas_occupancy(pr, rho[ns], rev)
            flow = passes * (occ + pr["s_fwd_nm3"].value * rho[ns])
            dn[ns] -= flow
            info.append((b.kind, b.stage, passes, flow))
            power += passes * e
            continue
        elif b.kind == "water":
            # Runs only as the sump holds water; a dry plunger draws gas.
            need = max(water, 0.0)
            full = passes * pr["water_per_pass"].value
            frac = min(1.0, need / max(full * 1e-3, 1e-30)) if need > 0 else 0.0
            passes *= frac
            dw -= passes * pr["water_per_pass"].value
            flow = 0.0
        power += passes * e
        info.append((b.kind, b.stage, passes, flow))
    return dn, dw, power, info


def time_course(device: dict, t_fill_max_s: float = 600.0, t_discharge_max_s: float = 120.0,
                samples: int = 400) -> dict:
    """Fill from the cushion until the tank reaches its target (or stalls, or
    the time runs out), then discharge to the cushion. State: O2 molecules in
    each chamber and the tank, water in the chamber-1 sump, and energy used.
    Stiff solver (scipy LSODA). Rotor speeds scale down together when the
    power demand exceeds the peak supply."""
    import numpy as np
    from scipy.integrate import solve_ivp

    ns = device["n_stages"]
    vols = [0.0] + [device["tank"]["chamber_volume_m3"].value] * (ns - 1) + [device["tank"]["volume_m3"].value]
    cushion_atm = device["tank"]["cushion_atm"].value
    target_atm = device["tank"]["target_atm"].value
    peak_w = device["power"]["peak_w"].value

    def p_atm(k, nn):
        return p_of_rho(nn[k] * NM3 / vols[k])

    def full_state(y):
        return [0.0] + [max(v, 0.0) for v in y[:ns]], max(y[ns], 0.0)

    def evaluate(y, phase):
        nn, w = full_state(y)
        _, _, demand, _ = _rates(device, nn, w, phase, 1.0)
        scale = min(1.0, peak_w / demand) if demand > 0 else 1.0
        dn, dw, pw, _ = _rates(device, nn, w, phase, scale)
        return dn, dw, pw, scale

    def rhs(phase):
        def f(_t, y):
            dn, dw, pw, _ = evaluate(y, phase)
            return dn[1:] + [dw, pw, -dn[ns] if phase == "discharge" else 0.0]
        return f

    y = [rho_of_p(cushion_atm) * vols[k] / NM3 for k in range(1, ns + 1)] + [0.0, 0.0, 0.0]
    log = {"t_s": [], "phase": [], "tank_atm": [], "chamber_atm": [], "water": [], "power_w": [],
           "duty": [], "tank_o2": []}
    events = {}
    t0 = 0.0
    water_peak = 0.0
    for phase, t_max in (("fill", t_fill_max_s), ("discharge", t_discharge_max_s)):
        if phase == "fill":
            def stop(_t, yy):
                return p_atm(ns, full_state(yy)[0]) - target_atm
        else:
            def stop(_t, yy):
                return p_atm(ns, full_state(yy)[0]) - cushion_atm
        stop.terminal = True
        t_eval = np.linspace(t0, t0 + t_max, max(2, samples // 2))
        sol = solve_ivp(rhs(phase), (t0, t0 + t_max), y, method="LSODA", events=stop,
                        t_eval=t_eval, rtol=1e-6, atol=1.0, max_step=t_max / 50)
        for i, t in enumerate(sol.t):
            yy = sol.y[:, i]
            nn, w = full_state(yy)
            _, _, pw, scale = evaluate(yy, phase)
            log["t_s"].append(float(t))
            log["phase"].append(phase)
            log["tank_atm"].append(p_atm(ns, nn))
            log["chamber_atm"].append([p_atm(k, nn) for k in range(1, ns)])
            log["water"].append(w)
            log["power_w"].append(pw)
            log["duty"].append(scale)
            log["tank_o2"].append(nn[ns])
            water_peak = max(water_peak, w)
        y_end = sol.y_events[0][0] if sol.status == 1 else sol.y[:, -1]
        t_end = sol.t_events[0][0] if sol.status == 1 else sol.t[-1]
        nn, _ = full_state(y_end)
        if phase == "fill":
            events["fill_end"] = "target" if sol.status == 1 else (
                "stalled" if abs(evaluate(y_end, phase)[0][ns]) < 1e-3 * nn[ns] / t_max else "time")
            events["fill_time_s"] = float(t_end - t0)
            events["tank_atm_after_fill"] = p_atm(ns, nn)
            events["tank_o2_after_fill"] = nn[ns]
            tank_series = [(t, pa) for t, pa, ph in zip(log["t_s"], log["tank_atm"], log["phase"]) if ph == "fill"]
            for mark in (100.0, 300.0, 500.0, 900.0, target_atm):
                hit = next((t for t, pa in tank_series if pa >= mark), None)
                if hit is not None:
                    events[f"t_tank_{int(mark)}_atm_s"] = float(hit - t0)
        else:
            events["discharge_end"] = "cushion" if sol.status == 1 else "time"
            events["discharge_time_s"] = float(t_end - t0)
            events["o2_delivered"] = float(y_end[ns + 2])
        y = list(y_end)
        y[ns + 2] = 0.0
        t0 = float(t_end)
    events["water_peak"] = water_peak
    # Liquid water takes 0.0299 nm^3 per molecule; the sump must hold it.
    events["water_peak_fraction_of_chamber"] = water_peak * 0.0299 * NM3 / vols[1]
    events["energy_per_delivered_o2_j"] = (
        float(y[ns + 1]) / events["o2_delivered"] if events.get("o2_delivered") else float("inf"))
    events["energy_j"] = float(y[ns + 1])
    events["mean_power_w"] = events["energy_j"] / max(t0, 1e-12)
    events["min_duty"] = min(log["duty"]) if log["duty"] else 1.0
    return {"events": events, "log": log}


def provenance_table(device: dict) -> list[dict]:
    """Every parameter of the device with its provenance."""
    rows = []
    for group in ("hull", "tank", "power", "section"):
        for name, p in device[group].items():
            rows.append({"component": group, "name": name, **asdict(p)})
    for b in device["banks"]:
        label = f"{b.kind}" + (f"_stage{b.stage}" if b.kind in ("load", "cascade") else "")
        for name, p in b.params.items():
            rows.append({"component": label, "name": name, **asdict(p)})
    return rows


def report(device: dict | None = None) -> dict:
    """Static budget and time course for the central, pessimistic and
    optimistic cases."""
    device = device or default_device()
    out = {"provenance": provenance_table(device), "cases": {}}
    for case in ("central", "pessimistic", "optimistic"):
        d = with_case(device, case)
        course = time_course(d)
        out["cases"][case] = {"budget": static_budget(d), "events": course["events"],
                              "log": course["log"] if case == "central" else None}
    return out


def _summary(rep: dict) -> str:
    lines = []
    for case, c in rep["cases"].items():
        b, e = c["budget"], c["events"]
        lines.append(f"[{case}]")
        lines.append(f"  hull: {b['sections_that_fit']} sections fit; {100 * b['hull_fraction_used']:.1f}% used "
                     f"(allowed {100 * b['hull_fraction_allowed']:.0f}%)")
        lines.append(f"  tank capacity at target {b['tank_capacity_o2_at_target']:.3g} O2, deliverable {b['deliverable_o2']:.3g}")
        lines.append("  stage stall (atm): " + ", ".join(f"{s:.0f}" for s in b["stage_stall_atm"]))
        lines.append(f"  power demand at design speed {b['power_demand_at_design_w']:.3g} W vs peak {b['power_peak_w']:.3g} W")
        lines.append(f"  fill: {e['fill_end']} after {e['fill_time_s']:.1f} s at {e['tank_atm_after_fill']:.0f} atm; "
                     + ", ".join(f"{k[7:-2]} at {v:.1f} s" for k, v in e.items() if k.startswith("t_tank_")))
        lines.append(f"  discharge: {e['discharge_end']} after {e['discharge_time_s']:.1f} s, delivered {e['o2_delivered']:.3g} O2")
        lines.append(f"  energy {e['energy_j']:.3g} J, mean power {e['mean_power_w']:.3g} W, min duty {e['min_duty']:.2f}, "
                     f"water peak {e['water_peak']:.3g} ({100 * e['water_peak_fraction_of_chamber']:.1f}% of chamber 1)")
        lines.append(f"  energy per delivered O2 {e['energy_per_delivered_o2_j']:.3g} J; "
                     f"infinite-cascade net per pass at target {b['infinite_cascade_net_at_target']:.2f}")
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description="Respirocyte system model (M29-01)")
    ap.add_argument("--stages", type=int, default=4)
    ap.add_argument("--json", help="write the full report here")
    a = ap.parse_args()
    rep = report(default_device(a.stages))
    print(_summary(rep))
    if a.json:
        with open(a.json, "w") as handle:
            json.dump(rep, handle, indent=1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
