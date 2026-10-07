"""Respirocyte top-down design (task M29-04): size a blood-resident device
from the circulation it must serve.

The system model (``respirocyte_system.py``) answers "given a device, what
does it do". This module answers the inverse question: "given the blood
circulation and the red-cell gas capacity, what device is best". It is the
L4 layer that Freitas did not do: he fixed 1 um and 1000 atm and asked
whether the device could be built. Here the diameter, the tank pressure, the
rotor count and the rotor speed are free, and the optimum is found.

Circulation model
-----------------
A device rides the blood. It loads O2 in the pulmonary capillaries and
releases it in the tissue capillaries. Each passage through a capillary bed
gives one loading window ``t_lung`` and one release window ``t_tissue``. The
rest of the circulation (arteries, veins) is storage. One full loop takes
``T_circ = blood_volume / cardiac_output``.

The red-cell target
-------------------
Normal blood carries about 201 mL O2 per litre (15 g/dL Hb, 1.34 mL/g). One
pass through the body releases the arterio-venous difference, about 50 mL per
litre at rest and 150 mL per litre at peak exercise. A respirocyte that
replaces the red cells must deliver that much O2 per litre of blood per loop,
at rest and at exercise.

Per-device model
----------------
A hull of diameter ``D`` has area ``pi D^2``. A measured rotor section takes
``13.7 x 6.9 nm`` of it. The sections that fit, split between the loading, the
discharge and the water bank, set the throughput. Each pocket pass costs the
cam work plus the compression work from the blood partial pressure to the tank
pressure. The power comes from glucose engines at the Freitas surface power
density; the loading peak (during ``t_lung``) and the discharge peak (during
``t_tissue``) must each fit under it. The tank must hold one cycle of O2 and
fit inside the rotor shell.

The optimum
-----------
Smaller devices need a smaller volume fraction (nanocrit), because power
scales with area and volume with ``D^3``. Lower pressure needs less
compression work but a bigger tank. The module sweeps diameter and pressure
and reports the design that meets the red-cell flux at the lowest nanocrit.

Results are simulated estimates. Nothing here is built or validated.
"""

from __future__ import annotations

import argparse
import math

from nanocad.respirocyte_system import (
    ATM_PA,
    AVOGADRO,
    EOS_ATM_RHO,
    KJ_PER_MOL_TO_J,
    NM3,
    P,
    Param,
    _implicit_step,
    _rates,
    load_occupancy,
    o2_rotor_parameters,
    p_of_rho,
    rho_of_p,
    section_parameters,
    stall_atm,
    with_case,
)

STD_MOLAR_VOLUME_M3 = 22.414e-3
MMHG_TO_ATM = 1.0 / 760.0


def physiology() -> dict[str, Param]:
    """Human circulation and blood-gas reference values."""
    return {
        "blood_volume_m3": P(5.0e-3, 4.5e-3, 5.5e-3, "m^3", "reference", "standard 70 kg adult; 5 L"),
        "cardiac_output_rest_m3_s": P(
            5.0e-3 / 60.0, 4.5e-3 / 60.0, 5.5e-3 / 60.0, "m^3/s", "reference", "5 L/min at rest"
        ),
        "cardiac_output_max_m3_s": P(
            25.0e-3 / 60.0, 20.0e-3 / 60.0, 35.0e-3 / 60.0, "m^3/s", "reference", "25 L/min at peak exercise"
        ),
        "t_lung_rest_s": P(0.75, 0.5, 1.0, "s", "reference", "pulmonary capillary transit at rest"),
        "t_lung_max_s": P(0.30, 0.25, 0.5, "s", "reference", "pulmonary capillary transit at peak exercise"),
        "t_tissue_rest_s": P(0.75, 0.5, 2.0, "s", "reference", "systemic capillary transit at rest"),
        "t_tissue_max_s": P(0.30, 0.25, 0.5, "s", "reference", "systemic capillary transit at peak exercise"),
        "arterial_po2_atm": P(
            100.0 * MMHG_TO_ATM,
            90.0 * MMHG_TO_ATM,
            104.0 * MMHG_TO_ATM,
            "atm",
            "reference",
            "PaO2 100 mmHg; the loading side",
        ),
        "venous_po2_atm": P(
            40.0 * MMHG_TO_ATM,
            35.0 * MMHG_TO_ATM,
            46.0 * MMHG_TO_ATM,
            "atm",
            "reference",
            "PvO2 40 mmHg; the release side",
        ),
        "hemoglobin_g_per_m3": P(150.0e3, 130.0e3, 180.0e3, "g/m^3", "reference", "15 g/dL"),
        "o2_per_hemoglobin_m3_per_g": P(1.34e-6, 1.30e-6, 1.39e-6, "m^3/g", "reference", "Hufner 1.34 mL O2 per g Hb"),
        "av_diff_rest_m3_per_m3": P(0.050, 0.040, 0.060, "m^3/m^3", "reference", "50 mL O2 per L blood at rest"),
        "av_diff_max_m3_per_m3": P(
            0.150, 0.130, 0.180, "m^3/m^3", "reference", "150 mL O2 per L blood at peak exercise"
        ),
        "glucose_energy_j_per_g": P(1.6e4, 1.5e4, 1.7e4, "J/g", "reference", "16 kJ/g; combustion"),
        "engine_efficiency": P(0.5, 0.3, 0.7, "", "assumed", "", "glucose engine chemical to mechanical"),
        "power_density_w_m2": P(
            7.3e-12 / (math.pi * (1.0e-6) ** 2),
            1.0e-12 / (math.pi * (1.0e-6) ** 2),
            2.0e-11 / (math.pi * (1.0e-6) ** 2),
            "W/m^2",
            "reference",
            "Freitas 12 glucose engines at 1 um (ADR-0100) over the hull area",
        ),
        "metabolic_budget_w": P(
            100.0, 20.0, 400.0, "W", "assumed", "", "chemical power the body can spare; basal metabolism is about 100 W"
        ),
    }


def rbc_reference(phys: dict[str, Param] | None = None) -> dict:
    """Red-cell O2 capacity and the per-litre per-loop delivery targets."""
    phys = phys or physiology()

    # O2 volumes are m^3 O2 per m^3 blood; a litre is 1e-3 m^3.
    def per_m3(vol_per_vol: float) -> float:
        return vol_per_vol / STD_MOLAR_VOLUME_M3 * AVOGADRO

    capacity_per_m3 = phys["hemoglobin_g_per_m3"].value * phys["o2_per_hemoglobin_m3_per_g"].value
    return {
        "capacity_o2_per_l_m3": capacity_per_m3 * 1.0e-3,
        "capacity_o2_per_m3_molecules": per_m3(capacity_per_m3),
        "capacity_o2_per_l_molecules": per_m3(capacity_per_m3) * 1.0e-3,
        "rest_deliver_per_m3_molecules": per_m3(phys["av_diff_rest_m3_per_m3"].value),
        "rest_deliver_per_l_molecules": per_m3(phys["av_diff_rest_m3_per_m3"].value) * 1.0e-3,
        "max_deliver_per_m3_molecules": per_m3(phys["av_diff_max_m3_per_m3"].value),
        "max_deliver_per_l_molecules": per_m3(phys["av_diff_max_m3_per_m3"].value) * 1.0e-3,
    }


def circulation(phys: dict[str, Param] | None = None) -> dict:
    """The transit model: loop time, capillary windows and the O2 the blood
    must carry per litre and per loop."""
    phys = phys or physiology()
    v = phys["blood_volume_m3"].value
    q_rest = phys["cardiac_output_rest_m3_s"].value
    q_max = phys["cardiac_output_max_m3_s"].value
    ref = rbc_reference(phys)
    return {
        "loop_time_rest_s": v / q_rest,
        "loop_time_max_s": v / q_max,
        "lung_window_rest_s": phys["t_lung_rest_s"].value,
        "lung_window_max_s": phys["t_lung_max_s"].value,
        "tissue_window_rest_s": phys["t_tissue_rest_s"].value,
        "tissue_window_max_s": phys["t_tissue_max_s"].value,
        "capacity_o2_per_l": ref["capacity_o2_per_l_m3"],
        "rest_deliver_per_l_molecules": ref["rest_deliver_per_l_molecules"],
        "max_deliver_per_l_molecules": ref["max_deliver_per_l_molecules"],
        "rest_flux_molecules_s": ref["rest_deliver_per_m3_molecules"] * v / (v / q_rest),
        "max_flux_molecules_s": ref["max_deliver_per_m3_molecules"] * v / (v / q_max),
    }


def _window(device: dict, y: list[float], phase: str, t_max_s: float, samples: int = 150) -> tuple[list[float], float, float]:
    """Integrate the compartment O2 for exactly ``t_max_s`` seconds in one
    phase (``fill`` or ``discharge``).

    State ``y`` is ``[n_1..n_ns, water, energy, delivered]``; it is mutated in
    place and returned. Backward Euler, fixed step, so the stiff small
    chambers stay stable without numpy. Returns ``(y, last_power_w,
    energy_j_added)``.
    """
    ns = device["n_stages"]
    peak_w = device["power"]["peak_w"].value

    def full_state(yy):
        return [0.0] + [max(v, 0.0) for v in yy[:ns]], max(yy[ns], 0.0)

    def evaluate(yy):
        nn, w = full_state(yy)
        _, _, demand, _ = _rates(device, nn, w, phase, 1.0)
        scale = min(1.0, peak_w / demand) if demand > 0 else 1.0
        return _rates(device, nn, w, phase, scale)

    def f(_t, yy):
        dn, dw, pw, _ = evaluate(yy)
        return dn[1:] + [dw, pw, (-dn[ns] if phase == "discharge" else 0.0)]

    e0 = y[ns + 1]
    dt = t_max_s / samples
    t = 0.0
    while t < t_max_s - 1e-12:
        step = min(dt, t_max_s - t)
        y = _implicit_step(f, t, y, step)
        t += step
    _, _, pw, _ = evaluate(y)
    return y, pw, y[ns + 1] - e0


def cycle_steady_state(
    device: dict, phys: dict[str, Param] | None = None, n_loops: int = 25, case: str = "central", samples: int = 150
) -> dict:
    """Replace the infinite blood reservoir with a circulation cycle.

    Each loop the device loads for ``t_lung`` and discharges for ``t_tissue``;
    the rest of the loop (arteries, veins) is storage and no O2 moves. The
    tank therefore sees short windows, not a steady reservoir, and only the
    O2 it loads in one lung transit can be delivered in one tissue transit.
    The loop repeats until the tank reaches a periodic steady state.

    Returns the per-loop history, the steady-state load and delivery, the
    tank-pressure swing, the metabolic power, and the nanocrit the device
    needs to meet the red-cell target at that delivery.
    """
    phys = phys or physiology()
    d = with_case(device, case)
    ns = d["n_stages"]
    vols = [0.0] + [d["tank"]["chamber_volume_m3"].value] * (ns - 1) + [d["tank"]["volume_m3"].value]
    cushion = d["tank"]["cushion_atm"].value
    y = [rho_of_p(cushion) * vols[k] / NM3 for k in range(1, ns + 1)] + [0.0, 0.0, 0.0]
    t_lung = phys["t_lung_rest_s"].value
    t_tissue = phys["t_tissue_rest_s"].value
    loop_time = phys["blood_volume_m3"].value / phys["cardiac_output_rest_m3_s"].value

    def tank_atm(yy):
        return p_of_rho(max(yy[ns - 1], 0.0) * NM3 / vols[ns])

    hist = []
    for i in range(n_loops):
        p_start = tank_atm(y)
        n0 = sum(y[:ns])
        y, pw_fill, e_fill = _window(d, y, "fill", t_lung, samples)
        n1 = sum(y[:ns])
        p_peak = tank_atm(y)
        y, pw_disch, e_disch = _window(d, y, "discharge", t_tissue, samples)
        p_end = tank_atm(y)
        delivered = y[ns + 2]
        y[ns + 2] = 0.0
        hist.append(
            {
                "loop": i + 1,
                "p_start_atm": p_start,
                "p_peak_atm": p_peak,
                "p_end_atm": p_end,
                "o2_loaded": n1 - n0,
                "o2_delivered": delivered,
                "fill_power_w": pw_fill,
                "discharge_power_w": pw_disch,
                "energy_j": e_fill + e_disch,
            }
        )
    last, prev = hist[-1], hist[-2] if len(hist) > 1 else hist[-1]
    delivered_ss = last["o2_delivered"]
    loaded_ss = last["o2_loaded"]
    conv = abs(delivered_ss - prev["o2_delivered"]) / delivered_ss if delivered_ss > 0 else float("inf")
    v_dev = math.pi * device["hull"]["diameter_m"].value**3 / 6.0
    ref = rbc_reference(phys)
    energy_ss = last["energy_j"]
    return {
        "n_loops": n_loops,
        "loop_time_s": loop_time,
        "t_lung_s": t_lung,
        "t_tissue_s": t_tissue,
        "history": hist,
        "steady_o2_loaded": loaded_ss,
        "steady_o2_delivered": delivered_ss,
        "steady_fill_fraction": loaded_ss / (delivered_ss if delivered_ss > 0 else float("nan")),
        "converged": conv < 0.01,
        "convergence": conv,
        "tank_swing_atm": [last["p_end_atm"], last["p_peak_atm"]],
        "target_atm": d["tank"]["target_atm"].value,
        "reaches_target": last["p_peak_atm"] >= 0.99 * d["tank"]["target_atm"].value,
        "energy_per_loop_j": energy_ss,
        "mean_mech_power_w": energy_ss / loop_time,
        "chem_power_w": energy_ss / loop_time / phys["engine_efficiency"].value,
        "metabolic_budget_w": phys["metabolic_budget_w"].value,
        "v_dev_m3": v_dev,
        "nanocrit_rest": ref["rest_deliver_per_m3_molecules"] * v_dev / delivered_ss if delivered_ss > 0 else float("inf"),
        "nanocrit_max": ref["max_deliver_per_m3_molecules"] * v_dev / delivered_ss if delivered_ss > 0 else float("inf"),
    }


def cycle_report(device: dict | None = None, phys: dict[str, Param] | None = None, n_loops: int = 25) -> str:
    """A text summary of the circulation-cycle result for one device."""
    device = device or _default_device_for_cycle()
    r = cycle_steady_state(device, phys, n_loops=n_loops)
    lines = [
        f"circulation cycle: {r['t_lung_s']:.2f} s lung, {r['t_tissue_s']:.2f} s tissue, "
        f"loop {r['loop_time_s']:.0f} s, {r['n_loops']} loops",
        f"  steady state: load {r['steady_o2_loaded']:.3g} O2, deliver {r['steady_o2_delivered']:.3g} O2 per loop "
        f"(converged {r['converged']}, drift {100 * r['convergence']:.2g}%)",
        f"  tank swings {r['tank_swing_atm'][0]:.1f} .. {r['tank_swing_atm'][1]:.1f} atm "
        f"of the {r['target_atm']:.0f} atm target; reaches target {r['reaches_target']}",
        f"  mechanical {r['mean_mech_power_w']:.3g} W, chemical {r['chem_power_w']:.3g} W "
        f"(budget {r['metabolic_budget_w']:.0f} W)",
        f"  nanocrit to meet the red-cell target: {100 * r['nanocrit_rest']:.3g}% rest, "
        f"{100 * r['nanocrit_max']:.3g}% peak exercise",
    ]
    return "\n".join(lines)


def _default_device_for_cycle() -> dict:
    from nanocad.respirocyte_system import default_device

    return default_device()


def compression_j_per_molecule(p_low_atm: float, p_high_atm: float) -> float:
    """Isothermal work to compress one O2 molecule from ``p_low`` to ``p_high``
    through the force-field equation of state (ADR-0123).

    The equation of state is piecewise linear in pressure, so ``1/rho`` is
    ``1/(a + b p)`` and each segment integrates in closed form. This avoids
    the coarse-step error of a fixed grid near the low-pressure knee.
    """
    if p_high_atm <= p_low_atm:
        return 0.0
    lo, hi = max(p_low_atm, 1.0e-9), p_high_atm
    total = 0.0
    for (p0, r0), (p1, r1) in zip(EOS_ATM_RHO, EOS_ATM_RHO[1:], strict=False):
        a0, a1 = max(lo, p0), min(hi, p1)
        if a1 <= a0:
            continue
        b = (r1 - r0) / (p1 - p0)
        a = r0 - b * p0  # rho(p)/nm^3 = a + b p
        if abs(b) < 1e-300:
            total += (a1 - a0) / (a * 1.0e27)
        else:
            total += math.log((a + b * a1) / (a + b * a0)) / (b * 1.0e27)
    # Extend the last segment above 1000 atm.
    if hi > EOS_ATM_RHO[-1][0]:
        (p0, r0), (p1, r1) = EOS_ATM_RHO[-2], EOS_ATM_RHO[-1]
        b = (r1 - r0) / (p1 - p0)
        a = r0 - b * p0
        total += math.log((a + b * hi) / (a + b * p1)) / (b * 1.0e27)
    return total * ATM_PA


def work_per_o2_j(
    phys: dict[str, Param], p_max_atm: float, params: dict[str, Param], cam_kj_mol: float | None = None
) -> dict:
    """Mechanical work per O2 delivered, per bank. ``cam_kj_mol`` overrides the
    cam work for a sensitivity test (23 kJ/mol is the reversible minimum)."""
    cam = (params["energy_per_pass_kj_mol"].value if cam_kj_mol is None else cam_kj_mol) * KJ_PER_MOL_TO_J
    load = cam + compression_j_per_molecule(phys["arterial_po2_atm"].value, p_max_atm)
    discharge = cam
    return {"cam_j": cam, "load_j": load, "discharge_j": discharge, "round_trip_j": load + discharge}


def design_point(
    diameter_m: float,
    p_max_atm: float,
    phys: dict[str, Param] | None = None,
    n_stages: int = 1,
    hull_fill: float = 0.6,
    f_load: float = 0.45,
    f_discharge: float = 0.45,
    f_water: float = 0.10,
    params: dict[str, Param] | None = None,
    cam_kj_mol: float | None = None,
) -> dict:
    """One candidate design. Returns per-device throughput, tank, power,
    nanocrit at rest and at exercise, and the feasibility flags.

    The power caps the rotor speed. Each bank uses the whole surface power at
    its own time (loading in the lung, discharge in the tissue), so the two do
    not compete for the peak.
    """
    phys = phys or physiology()
    params = params or o2_rotor_parameters()
    sec = section_parameters()
    a_sec = sec["footprint_y_m"].value * sec["footprint_z_m"].value
    depth = sec["depth_m"].value

    r = 0.5 * diameter_m
    a_hull = math.pi * diameter_m**2
    v_dev = math.pi * diameter_m**3 / 6.0
    p_peak = phys["power_density_w_m2"].value * a_hull

    n_sec = hull_fill * a_hull / a_sec
    n_load = f_load * n_sec
    n_disch = f_discharge * n_sec
    n_water = f_water * n_sec

    r_in = max(r - depth, 0.0)
    v_tank_avail = 4.0 / 3.0 * math.pi * r_in**3

    work = work_per_o2_j(phys, p_max_atm, params, cam_kj_mol)
    pockets = params["pockets"].value
    max_rev = params["max_rev_per_s"].value

    def bank_throughput(n_rotors, work_j, occ_at_rev):
        if n_rotors <= 0:
            return 0.0, 0.0, 0.0
        # A passive (downhill, zero-work) bank is not power limited; it turns
        # at the cap. A driven bank is capped by the surface power.
        rev = max_rev if work_j <= 0 else min(p_peak / (n_rotors * pockets * work_j), max_rev)
        occ = occ_at_rev(rev)
        return n_rotors * pockets * rev * occ, rev, occ

    load_occ = lambda rev: load_occupancy(params, rev)  # noqa: E731
    gas_occ = lambda rev: params["gas_occupancy"].value * (  # noqa: E731
        1.0
        - math.exp(
            -params["gas_binding_per_s_per_rho"].value
            * rho_of_p(p_max_atm)
            * params["blood_arc_fraction"].value
            / max(rev, 1e-12)
        )
    )

    load_rate, rev_load, occ_load = bank_throughput(n_load, work["load_j"], load_occ)
    disch_rate, rev_disch, occ_disch = bank_throughput(n_disch, work["discharge_j"], gas_occ)

    n_load_cycle = load_rate * phys["t_lung_rest_s"].value
    n_disch_cycle = disch_rate * phys["t_tissue_rest_s"].value
    n_tank = rho_of_p(p_max_atm) * v_tank_avail / NM3
    n_cycle = min(n_load_cycle, n_disch_cycle, n_tank)
    limiting = min((n_load_cycle, "load"), (n_disch_cycle, "discharge"), (n_tank, "tank"))[1]

    ref = rbc_reference(phys)

    def nanocrit(deliver_per_m3):
        return v_dev * deliver_per_m3 / n_cycle if n_cycle > 0 else float("inf")

    stall = stall_atm(params, None, occ_load)
    # Peak power the banks would like, and the metabolic average over a loop.
    t_rest = phys["blood_volume_m3"].value / phys["cardiac_output_rest_m3_s"].value
    t_max = phys["blood_volume_m3"].value / phys["cardiac_output_max_m3_s"].value
    peak_mech = max(n_load * pockets * rev_load * work["load_j"], n_disch * pockets * rev_disch * work["discharge_j"])
    rest_flux = ref["rest_deliver_per_m3_molecules"] * phys["blood_volume_m3"].value / t_rest
    max_flux = ref["max_deliver_per_m3_molecules"] * phys["blood_volume_m3"].value / t_max
    mech_rest = rest_flux * work["round_trip_j"] * n_stages
    mech_max = max_flux * work["round_trip_j"] * n_stages
    chem_rest = mech_rest / phys["engine_efficiency"].value
    chem_max = mech_max / phys["engine_efficiency"].value

    return {
        "diameter_m": diameter_m,
        "p_max_atm": p_max_atm,
        "n_stages": n_stages,
        "v_dev_m3": v_dev,
        "a_hull_m2": a_hull,
        "p_peak_w": p_peak,
        "n_sections": n_sec,
        "n_load": n_load,
        "n_discharge": n_disch,
        "n_water": n_water,
        "rev_load_per_s": rev_load,
        "rev_discharge_per_s": rev_disch,
        "occ_load": occ_load,
        "load_rate_o2_s": load_rate,
        "discharge_rate_o2_s": disch_rate,
        "n_cycle_o2": n_cycle,
        "limiting": limiting,
        "n_load_cycle": n_load_cycle,
        "n_discharge_cycle": n_disch_cycle,
        "n_tank": n_tank,
        "v_tank_needed_m3": n_cycle / rho_of_p(p_max_atm) * NM3 if rho_of_p(p_max_atm) > 0 else float("inf"),
        "v_tank_avail_m3": v_tank_avail,
        "tank_fits": n_cycle <= n_tank,
        "nanocrit_rest": nanocrit(ref["rest_deliver_per_m3_molecules"]),
        "nanocrit_max": nanocrit(ref["max_deliver_per_m3_molecules"]),
        "work_load_j": work["load_j"],
        "work_discharge_j": work["discharge_j"],
        "work_round_trip_j": work["round_trip_j"],
        "peak_mech_w": peak_mech,
        "mech_power_rest_w": mech_rest,
        "mech_power_max_w": mech_max,
        "chem_power_rest_w": chem_rest,
        "chem_power_max_w": chem_max,
        "metabolic_budget_w": phys["metabolic_budget_w"].value,
        "stall_atm": stall,
        "stall_ok": p_max_atm <= stall,
    }


def sweep(
    phys: dict[str, Param] | None = None,
    n_stages: int = 1,
    diameters_m: list[float] | None = None,
    pressures_atm: list[float] | None = None,
    hull_fill: float = 0.6,
) -> list[dict]:
    phys = phys or physiology()
    diameters_m = diameters_m or [0.1e-6, 0.2e-6, 0.3e-6, 0.5e-6, 0.75e-6, 1.0e-6, 1.5e-6, 2.0e-6]
    pressures_atm = pressures_atm or [50.0, 100.0, 200.0, 300.0, 500.0, 1000.0]
    return [
        design_point(d, p, phys, n_stages=n_stages, hull_fill=hull_fill) for d in diameters_m for p in pressures_atm
    ]


def optimize(
    phys: dict[str, Param] | None = None,
    n_stages: int = 1,
    budget_w: float | None = None,
    objective: str = "nanocrit_max",
    cam_kj_mol: float | None = None,
) -> dict:
    """The feasible design with the lowest nanocrit at peak exercise.

    Feasible means: the tank fits and the pressure is below the single-stage
    stall. If ``budget_w`` is given, the metabolic chemical power at peak
    exercise must also be inside it. The pressure is capped at the stall, so a
    cascade is never needed.
    """
    phys = phys or physiology()
    best = None
    for d in [0.05e-6 + 0.01e-6 * i for i in range(0, 196)]:  # 0.05 to 2.0 um
        for p in [25.0 + 5.0 * i for i in range(0, 90)]:  # 25 to 470 atm
            pt = design_point(d, p, phys, n_stages=n_stages, cam_kj_mol=cam_kj_mol)
            if not pt["tank_fits"] or not pt["stall_ok"]:
                continue
            if budget_w is not None and pt["chem_power_max_w"] > budget_w:
                continue
            if best is None or pt[objective] < best[objective]:
                best = pt
    return best or {}


def _fmt(x, unit="", scale=1.0, prec=3):
    if x is None:
        return "n/a"
    return f"{x / scale:.{prec}g}{unit}"


def summary(pt: dict) -> str:
    return "\n".join(
        [
            f"device {_fmt(pt['diameter_m'], ' um', 1e-6)} diameter, tank "
            f"{_fmt(pt['p_max_atm'], ' atm')}, {pt['n_stages']} stage(s); limited by {pt['limiting']}",
            f"  hull area {_fmt(pt['a_hull_m2'], ' um^2', 1e-12)} = {pt['n_sections']:.0f} rotor sections "
            f"(8 pockets each); load {pt['n_load']:.0f}, discharge {pt['n_discharge']:.0f}, water {pt['n_water']:.0f}",
            f"  speed: load {pt['rev_load_per_s']:.0f} rev/s (occ {pt['occ_load']:.2f}), "
            f"discharge {pt['rev_discharge_per_s']:.0f} rev/s",
            f"  per cycle {pt['n_cycle_o2']:.3g} O2 "
            f"(load {pt['n_load_cycle']:.3g}, discharge {pt['n_discharge_cycle']:.3g}, tank {pt['n_tank']:.3g})",
            f"  tank {_fmt(pt['v_tank_needed_m3'], ' nm^3', 1e-27)} needed of "
            f"{_fmt(pt['v_tank_avail_m3'], ' nm^3', 1e-27)} available; fits {pt['tank_fits']}",
            f"  work per O2: load {pt['work_load_j']:.2e} J, discharge {pt['work_discharge_j']:.2e} J",
            f"  nanocrit: rest {100 * pt['nanocrit_rest']:.3g}%, peak exercise {100 * pt['nanocrit_max']:.3g}%",
            f"  power: peak mechanical {pt['peak_mech_w']:.2e} W; metabolic "
            f"{pt['chem_power_rest_w']:.1f} W rest, {pt['chem_power_max_w']:.1f} W exercise "
            f"(budget {pt['metabolic_budget_w']:.0f} W)",
            f"  single-stage stall {pt['stall_atm']:.0f} atm; target ok {pt['stall_ok']}",
        ]
    )


def main() -> int:
    ap = argparse.ArgumentParser(description="Respirocyte top-down design (M29-04)")
    ap.add_argument("--stages", type=int, default=1)
    ap.add_argument("--budget", type=float, help="metabolic budget in W")
    ap.add_argument("--sweep", action="store_true", help="print the diameter-pressure table")
    ap.add_argument("--cycle", action="store_true", help="run the lung/tissue circulation cycle (M29-02)")
    ap.add_argument("--loops", type=int, default=25, help="loops for --cycle")
    a = ap.parse_args()
    if a.cycle:
        print(cycle_report(n_loops=a.loops))
        return 0
    if a.sweep:
        print("D(um)  P(atm)  sections  limited   n_cycle      nanocrit_rest  nanocrit_max  chem_max_W")
        for pt in sweep(n_stages=a.stages):
            print(
                f"{pt['diameter_m']/1e-6:5.2f}  {pt['p_max_atm']:6.0f}  "
                f"{pt['n_sections']:8.0f}  {pt['limiting']:9s} {pt['n_cycle_o2']:10.3g}  "
                f"{100*pt['nanocrit_rest']:8.3g}%  {100*pt['nanocrit_max']:8.3g}%  "
                f"{pt['chem_power_max_w']:9.3g}"
            )
        return 0
    best = optimize(n_stages=a.stages, budget_w=a.budget)
    if not best:
        print("no feasible design in the sweep")
        return 1
    print("lowest nanocrit at peak exercise:")
    print(summary(best))
    rev = optimize(n_stages=a.stages, budget_w=a.budget, cam_kj_mol=23.0)
    if rev:
        print("\nwith the reversible cam work (23 kJ/mol, the thermodynamic minimum):")
        print(summary(rev))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
