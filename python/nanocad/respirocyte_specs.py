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
lifetime is set by the autoxidation rate alone. That rate is measured for free
heme and for hemoglobin, but not for a synthetic pocket, so this module treats
it as a parameter and reports the range. The ring damage by superoxide and
peroxide (M28-22) is the same risk seen from the atomistic side.

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

# Autoxidation of the heme iron, per second. Half-life is ln2 / k.
AUTOXIDATION_PER_S = {
    "free heme": 1.0e-3,  # minutes
    "hemoglobin, no reductase": 1.4e-5,  # about 14 h
    "myoglobin": 1.2e-7,  # about 67 days
    "engineered pocket": 1.0e-8,  # about 800 days
}

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
    lifetime_s = autoxidation_lifetime_s(AUTOXIDATION_PER_S["myoglobin"])
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
    heme_life = autoxidation_lifetime_s(AUTOXIDATION_PER_S["myoglobin"])
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
    """Site lifetime and daily replacement against the autoxidation rate."""
    rows = []
    for name, k in AUTOXIDATION_PER_S.items():
        life = autoxidation_lifetime_s(k)
        rows.append(
            {
                "case": name,
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
        f"{'case':<26} {'k (1/s)':>9} {'life':>10} {'replace':>12}",
    ]
    for r in lifetime_table(dend_mass):
        lines.append(
            f"{r['case']:<26} {r['k_per_s']:9.1e} {r['lifetime_days']:8.1f} d "
            f"{r['replacement_per_day_g']:10.3g} g/d"
        )
    lines.append(
        "the pocket must slow autoxidation below about 1e-6 /s (life under 8 d) " "for a practical daily dose."
    )
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description="Device specs and lifetimes (M29-06)")
    ap.parse_args()
    print(report())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
