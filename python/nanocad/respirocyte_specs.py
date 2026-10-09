"""Device specifications, counts, and the site-lifetime model (task M29-06).

The design optimization (``respirocyte_design.py``) chooses an architecture
and its volume fraction. This module turns each design into a buildable
specification sheet: how many devices, how much material, how much O2 each
carries, what each costs in power, and how long the binding sites last before
they must be repaired or replaced.

The site lifetime is the least certain number in the whole program. A heme
loses O2 capacity when its Fe(II) oxidizes to Fe(III) (methemoglobin). A red
cell lasts 120 days because it carries a methemoglobin reductase that repairs
the iron. A synthetic carrier has no repair unless one is engineered in, so its
lifetime is set by the autoxidation rate alone. The reference rates below are
measured (each with its PubMed id); the rate of our synthetic pocket is not, so
``pocket_autoxidation`` estimates it with the two-channel mechanism of Brantley
et al. (1993), calibrated on hemoglobin in plasma, and scaled by two pocket
factors (ADR-0136). The ring damage by superoxide and peroxide (M28-22) is the
same risk seen from the atomistic side.

Results are simulated estimates. Nothing here is built or validated.
"""

from __future__ import annotations

import argparse
import math

from nanocad.respirocyte_design import (
    carrier_design,
    carrier_reference,
    physiology,
    pump_designs,
)
from nanocad.respirocyte_system import AVOGADRO

STD_MOLAR_VOLUME_M3 = 22.414e-3
SECONDS_PER_DAY = 86400.0

HOUR_S = 3600.0
LN2 = math.log(2.0)

# Measured autoxidation of heme O2 carriers, pH about 7, in air unless noted.
# (name, k per s, temperature C, PubMed id, note). Half-life is ln2 / k.
AUTOXIDATION_REFERENCES: tuple[tuple[str, float, float, str, str], ...] = (
    ("myoglobin H64A, no distal H-bond", 58.0 / HOUR_S, 37.0, "8463233",
     "the distal His removed; the fastest Brantley mutant"),
    ("albumin-heme, covalent, red-cell p50", LN2 / (5.0 * HOUR_S), 37.0, "15656571",
     "rHSA(FeP-Glu), p50 27 Torr"),
    ("myoglobin, native", 0.055 / HOUR_S, 37.0, "8463233", "sperm whale, Brantley et al. 1993"),
    ("hemoglobin A", 0.029 / HOUR_S, 37.0, "3422420", "HbO2 A, Hebbel et al. 1988; 0.025 /h in PMID 9753633"),
    ("albumin-heme FecycP(His), best synthetic", LN2 / (25.0 * HOUR_S), 37.0, "12009926",
     "longest-lived synthetic heme carrier found"),
    ("heme-cyclodextrin (hemoCD)", LN2 / (30.1 * HOUR_S), 25.0, "19178464",
     "measured at 25 C, so faster at 37 C"),
    ("myoglobin L29F, crowded pocket", 0.005 / HOUR_S, 37.0, "8463233",
     "a bulky residue keeps water from the iron"),
)
AUTOXIDATION_PER_S = {name: k for name, k, _, _, _ in AUTOXIDATION_REFERENCES}

# A practical daily dose needs the sites to last about a week. Not a measurement.
AUTOXIDATION_TARGET_PER_S = 1.0e-6

# Hemoglobin in plasma at 37 C, half-life against pO2 (Furlan and Bachofen
# 1976, PubMed 1248570). Autoxidation is fastest at low pO2.
HB_PLASMA_HALF_LIFE_H = ((126.0, 20.0), (57.0, 12.0), (23.0, 7.0))
HB_P50_MMHG, HB_HILL_N = 26.8, 2.7

# Where a circulating carrier spends its time: arterial against venous and
# capillary blood. Assumed from the usual blood-volume split.
BLOOD_PO2_TIME_WEIGHTS = ((95.0, 0.2), (40.0, 0.8))

# Pocket factors relative to hemoglobin (ADR-0136). (optimistic, central, pessimistic)
# f_ox scales the superoxide (HO2) dissociation from the oxy site. Brantley: the
# distal-His H-bond suppresses it about 1000-fold (H64A). Our guard OH makes a
# His-like H-bond (ADR-0131), but a hydroxyl is also a proton donor, and pocket
# polarity speeds protonation (V68T/S). Assumed.
POCKET_F_OX = (0.3, 1.0, 30.0)
# f_bi scales the bimolecular O2 + deoxy-site channel, which needs a water on
# the iron. Our pocket stays dry in MD (ADR-0130); Brantley's crowded L29F
# pocket is about 10 times slower. Estimated.
POCKET_F_BI = (0.01, 0.1, 1.0)

# Assumed cycles a diamond MEMS rotor survives. Not measured.
MEMS_CYCLES_TO_FAILURE = 1.0e11

# Red-cell reference values.
RBC_COUNT_PER_M3 = 5.0e12 * 1.0e3
RBC_LIFETIME_S = 120.0 * SECONDS_PER_DAY
RBC_HCT = 0.45
RBC_CAPACITY_M3_PER_M3 = 0.201


def rbc_spec() -> dict:
    """The red cell, the thing every device is compared with."""
    return {
        "name": "red blood cell",
        "count_per_l": 5.0e12,
        "count_body": 2.5e13,
        "volume_fraction": RBC_HCT,
        "capacity_o2_per_l": 0.201,
        "o2_per_l_cell": RBC_CAPACITY_M3_PER_M3 / RBC_HCT,
        "lifetime_s": RBC_LIFETIME_S,
        "lifetime_days": RBC_LIFETIME_S / SECONDS_PER_DAY,
        "replacement_per_day": 2.5e13 / (RBC_LIFETIME_S / SECONDS_PER_DAY),
        "power_w": 0.0,
    }


def autoxidation_lifetime_s(k_per_s: float) -> float:
    """Site half-life from an autoxidation rate."""
    return math.log(2.0) / k_per_s if k_per_s > 0 else float("inf")


def _hill(p_mmhg: float, p50_mmhg: float, hill_n: float) -> float:
    if p_mmhg <= 0.0:
        return 0.0
    x = (p_mmhg / p50_mmhg) ** hill_n
    return x / (1.0 + x)


def fit_two_channel() -> dict:
    """Fit k(p) = k_ox Y(p) + k_bi p (1 - Y(p)) to hemoglobin in plasma.

    The first channel is the dissociation of superoxide from the oxy site; the
    second is the bimolecular reaction of O2 with a deoxy site that holds a
    water (Brantley et al. 1993). Least squares over the three points. The
    fitted k_ox is checked against the independent HbA measurement in air.
    """
    rows = []
    for p, t_h in HB_PLASMA_HALF_LIFE_H:
        y = _hill(p, HB_P50_MMHG, HB_HILL_N)
        rows.append((y, p * (1.0 - y), LN2 / (t_h * HOUR_S)))
    s11 = sum(a * a for a, _, _ in rows)
    s12 = sum(a * b for a, b, _ in rows)
    s22 = sum(b * b for _, b, _ in rows)
    t1 = sum(a * k for a, _, k in rows)
    t2 = sum(b * k for _, b, k in rows)
    det = s11 * s22 - s12 * s12
    k_ox = (t1 * s22 - t2 * s12) / det
    k_bi = (s11 * t2 - s12 * t1) / det
    resid = [k - (k_ox * a + k_bi * b) for a, b, k in rows]
    return {
        "k_ox_per_s": k_ox,
        "k_bi_per_s_mmhg": k_bi,
        "max_rel_residual": max(abs(r) / k for r, (_, _, k) in zip(resid, rows)),
        "hba_air_check_ratio": k_ox / AUTOXIDATION_PER_S["hemoglobin A"],
    }


def autoxidation_rate(p_mmhg: float, p50_mmhg: float, hill_n: float, f_ox: float = 1.0,
                      f_bi: float = 1.0, fit: dict | None = None) -> float:
    """Autoxidation rate of a site at one pO2, from the two-channel model."""
    fit = fit or fit_two_channel()
    y = _hill(p_mmhg, p50_mmhg, hill_n)
    return f_ox * fit["k_ox_per_s"] * y + f_bi * fit["k_bi_per_s_mmhg"] * p_mmhg * (1.0 - y)


def pocket_autoxidation(p50_mmhg: float = HB_P50_MMHG, hill_n: float = 2.0) -> list[dict]:
    """The circulating carrier's rate, time-weighted over arterial and venous blood.

    Three cases from the pocket factor ranges, plus hemoglobin itself (both
    factors 1, the hemoglobin p50) as the calibration point.
    """
    fit = fit_two_channel()
    cases = [("hemoglobin-like pocket", HB_P50_MMHG, HB_HILL_N, 1.0, 1.0)]
    for label, i in (("optimistic", 0), ("central", 1), ("pessimistic", 2)):
        cases.append((label, p50_mmhg, hill_n, POCKET_F_OX[i], POCKET_F_BI[i]))
    out = []
    for label, p50, n, f_ox, f_bi in cases:
        k = sum(w * autoxidation_rate(p, p50, n, f_ox, f_bi, fit) for p, w in BLOOD_PO2_TIME_WEIGHTS)
        ox = sum(w * f_ox * fit["k_ox_per_s"] * _hill(p, p50, n) for p, w in BLOOD_PO2_TIME_WEIGHTS)
        out.append({
            "case": label, "p50_mmhg": p50, "hill_n": n, "f_ox": f_ox, "f_bi": f_bi,
            "k_per_s": k, "half_life_h": LN2 / k / HOUR_S, "superoxide_share": ox / k,
            "suppression_to_target": k / AUTOXIDATION_TARGET_PER_S,
        })
    return out


def replacement_per_day(mass_kg: float, lifetime_s: float) -> float:
    """Mass that must be replaced each day for an exponential loss."""
    if lifetime_s <= 0 or not math.isfinite(lifetime_s):
        return 0.0
    lost_fraction = 1.0 - 2.0 ** (-SECONDS_PER_DAY / lifetime_s)
    return mass_kg * lost_fraction


def carrier_spec(material: dict, phys=None) -> dict:
    """A carrier molecule: count, mass, capacity, and site lifetime."""
    phys = phys or physiology()
    d = carrier_design(material, phys)
    v_blood = phys["blood_volume_m3"].value
    kg_per_m3 = d["kg_per_m3"]
    mass_body = kg_per_m3 * v_blood
    mw = material["molar_mass_kg_per_mol"]
    n_molecules_per_m3 = kg_per_m3 / mw * AVOGADRO if mw > 0 else 0.0
    n_sites_per_m3 = n_molecules_per_m3 * material["o2_sites_per_molecule"]
    n = material["hill_n"] if material["p50_mmhg"] else 2.0
    k = next(r["k_per_s"] for r in pocket_autoxidation(material["p50_mmhg"] or HB_P50_MMHG, n)
             if r["case"] == "central")
    lifetime_s = autoxidation_lifetime_s(k)
    return {
        "architecture": "carrier",
        "name": material["name"],
        "provenance": material["provenance"],
        "volume_fraction": d["nanocrit"],
        "count_per_l": n_molecules_per_m3 * 1.0e-3,
        "count_body": n_molecules_per_m3 * v_blood,
        "sites_body": n_sites_per_m3 * v_blood,
        "mass_body_kg": mass_body,
        "mass_body_g": mass_body * 1.0e3,
        "capacity_o2_per_l": 0.201,
        "o2_per_l_device": 0.201 / d["nanocrit"],
        "power_body_w": 0.0,
        "power_per_device_w": 0.0,
        "site_lifetime_days": lifetime_s / SECONDS_PER_DAY,
        "replacement_per_day_g": replacement_per_day(mass_body, lifetime_s) * 1.0e3,
    }


def pump_spec(pt: dict, phys=None) -> dict:
    """A pump device: count, per-device power, tank, and site lifetime (the
    shorter of the rotor wear and the heme autoxidation)."""
    phys = phys or physiology()
    v_dev = pt["v_dev_m3"]
    n_per_m3 = pt["nanocrit_max"] / v_dev if v_dev > 0 else 0.0
    v_blood = phys["blood_volume_m3"].value
    n_body = n_per_m3 * v_blood
    mech_life = MEMS_CYCLES_TO_FAILURE / max(pt["rev_load_per_s"], 1e-9)
    heme_life = autoxidation_lifetime_s(
        next(r["k_per_s"] for r in pocket_autoxidation() if r["case"] == "central"))
    life = min(mech_life, heme_life)
    return {
        "architecture": "pump",
        "name": pt["name"],
        "provenance": pt["provenance"],
        "diameter_m": pt["diameter_m"],
        "volume_fraction": pt["nanocrit_max"],
        "count_per_l": n_per_m3 * 1.0e-3,
        "count_body": n_body,
        "mass_body_kg": float("nan"),
        "mass_body_g": float("nan"),
        "capacity_o2_per_l": 0.201,
        "o2_per_l_device": 0.201 / pt["nanocrit_max"],
        "power_body_w": pt["chem_power_max_w"],
        "power_rest_w": pt["chem_power_rest_w"],
        "power_per_device_w": pt["chem_power_max_w"] / n_body if n_body > 0 else 0.0,
        "tank_atm": pt["p_max_atm"],
        "tank_volume_m3": pt["v_tank_needed_m3"],
        "rev_per_s": pt["rev_load_per_s"],
        "site_lifetime_days": life / SECONDS_PER_DAY,
        "replacement_per_day_g": float("nan"),
    }


def lifetime_table(mass_kg: float) -> list[dict]:
    """Site lifetime and daily replacement for each measured reference, fastest first."""
    rows = []
    for name, k, temp_c, pmid, _ in sorted(AUTOXIDATION_REFERENCES, key=lambda r: -r[1]):
        life = autoxidation_lifetime_s(k)
        rows.append(
            {
                "case": name,
                "pmid": pmid,
                "temp_c": temp_c,
                "k_per_s": k,
                "lifetime_days": life / SECONDS_PER_DAY,
                "replacement_per_day_g": replacement_per_day(mass_kg, life) * 1.0e3,
            }
        )
    return rows


def report(phys=None) -> str:
    phys = phys or physiology()
    rbc = rbc_spec()
    lines = [
        "Device specifications and counts",
        f"reference red cell: {rbc['count_body']:.2g} cells, {100 * rbc['volume_fraction']:.0f}% by volume, "
        f"{1000 * rbc['capacity_o2_per_l']:.0f} mL O2/L, {rbc['lifetime_days']:.0f} d, no engine",
        "",
        "carriers:",
        f"{'material':<18} {'vol %':>6} {'molecules/L':>12} {'mass (g)':>9} {'O2/L dev':>9} "
        f"{'life (d)':>9} {'replace (g/d)':>13}",
    ]
    for m in carrier_reference():
        s = carrier_spec(m, phys)
        if not math.isfinite(s["count_per_l"]) or s["count_per_l"] == 0:
            continue
        lines.append(
            f"{s['name']:<18} {100 * s['volume_fraction']:6.2f} {s['count_per_l']:12.2e} "
            f"{s['mass_body_g']:9.2g} {1000 * s['o2_per_l_device']:8.0f}L {s['site_lifetime_days']:9.1f} "
            f"{s['replacement_per_day_g']:13.3g}"
        )
    lines += [
        "",
        "pumps (0.5 um):",
        f"{'mode':<20} {'vol %':>6} {'devices/L':>11} {'W body':>8} {'W/device':>10} "
        f"{'tank atm':>9} {'life (d)':>9}",
    ]
    for pt in pump_designs(phys):
        s = pump_spec(pt, phys)
        lines.append(
            f"{s['name']:<20} {100 * s['volume_fraction']:6.2f} {s['count_per_l']:11.2e} "
            f"{s['power_body_w']:8.1f} {s['power_per_device_w']:10.2e} {s['tank_atm']:9.0f} "
            f"{s['site_lifetime_days']:9.1f}"
        )
    dend = next(m for m in carrier_reference() if m["name"] == "dendrimer-heme")
    dend_mass = carrier_spec(dend, phys)["mass_body_kg"]
    lines += [
        "",
        f"site lifetime vs the autoxidation rate (dendrimer-heme, {dend_mass * 1e3:.0f} g body):",
        f"{'measured reference':<42} {'k (1/s)':>9} {'life':>10} {'replace':>12}",
    ]
    for r in lifetime_table(dend_mass):
        lines.append(
            f"{r['case']:<42} {r['k_per_s']:9.1e} {r['lifetime_days']:8.2f} d "
            f"{r['replacement_per_day_g']:10.3g} g/d  PMID {r['pmid']}"
        )
    fit = fit_two_channel()
    lines += [
        "",
        f"two-channel fit to hemoglobin in plasma: k_ox {fit['k_ox_per_s']:.2e} /s, "
        f"k_bi {fit['k_bi_per_s_mmhg']:.2e} /s/mmHg, worst residual {100 * fit['max_rel_residual']:.0f}%, "
        f"k_ox / HbA in air {fit['hba_air_check_ratio']:.2f}",
        f"{'our pocket (dendrimer-heme, p50 26.8)':<38} {'f_ox':>5} {'f_bi':>5} {'k (1/s)':>9} "
        f"{'half-life':>10} {'HO2 share':>9} {'x target':>9}",
    ]
    for r in pocket_autoxidation(26.8, 2.0):
        lines.append(
            f"{r['case']:<38} {r['f_ox']:5.2g} {r['f_bi']:5.2g} {r['k_per_s']:9.1e} "
            f"{r['half_life_h']:8.1f} h {100 * r['superoxide_share']:8.0f}% {r['suppression_to_target']:9.1f}"
        )
    lines.append("p50 lever (central pocket factors): a higher affinity spends less time deoxygenated")
    for p50 in (10.0, 26.8, 40.0):
        r = next(x for x in pocket_autoxidation(p50, 2.0) if x["case"] == "central")
        lines.append(f"  p50 {p50:5.1f} mmHg: half-life {r['half_life_h']:6.1f} h")
    lines.append(
        f"target: k below {AUTOXIDATION_TARGET_PER_S:.0e} /s (about 8 days). No measured heme carrier "
        "reaches it; the best is myoglobin L29F at 5.8 days."
    )
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description="Device specs and lifetimes (M29-06)")
    ap.parse_args()
    print(report())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
