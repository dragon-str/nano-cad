"""A self-sustaining repair strategy for the heme carrier site (M29-14).

The site degrades in two ways (ADR-0108, ADR-0136):

  * the iron oxidizes, Fe(II) -> Fe(III), and stops binding O2. This is
    reversible: a one-electron donor on the empty arc re-reduces it. The
    electrons are trivial and come from the device, so no material is used.
  * the superoxide from that same oxidation attacks the porphyrin ring. This
    is irreversible: a broken ring cannot be rebuilt without new carbon.

So a self-sustaining site cannot repair the ring. It must stop the ring from
breaking. Three permanent elements do that:

  1. a distal pocket that suppresses autoxidation (a low superoxide rate);
  2. a co-located catalyst that dismutates the superoxide before it reaches
     the ring, and then decomposes the peroxide;
  3. the iron-repair donor, so the site stays in the Fe(II) state.

None of the three is consumed. The only inputs are electrons (a charged
onboard store) and protons (buffered blood). The alternative, a consumable
antioxidant, would need new material each day and is rejected.

This module estimates the ring lifetime from the race between the catalyst,
escape, and the ring attack, and finds the design that outlives the rotor.
All rates here are estimates; the model shows the design, not a measurement.
"""

from __future__ import annotations

import math

from nanocad import respirocyte_specs as rs

ELEMENTARY_CHARGE_C = 1.602176634e-19
DAYS_PER_YEAR = 365.25
HOUR_S = 3600.0


def central_rate_per_s() -> float:
    """The ADR-0136 central pocket autoxidation rate (about 31 h half-life)."""
    fit = rs.fit_two_channel()
    return sum(
        w * rs.autoxidation_rate(p, rs.HB_P50_MMHG, 2.0, rs.POCKET_F_OX[1], rs.POCKET_F_BI[1], fit)
        for p, w in rs.BLOOD_PO2_TIME_WEIGHTS
    )


# The pocket suppresses the autoxidation rate by a factor. The two channels
# (ADR-0136) are suppressed together: a distal H-bond holds the bound O2 and
# slows the superoxide channel (Brantley, PMID 8463233, about 1000x across the
# His mutants), and a dry pocket keeps water off the iron and slows the water
# channel (ADR-0130). The central value is the ADR-0136 pocket; the design goal
# is a 100-fold suppression. (label, suppression).
POCKET_CASES: tuple[tuple[str, float], ...] = (
    ("central (as modelled)", 1.0),
    ("H-bonded, dry pocket", 10.0),
    ("optimized distal pocket", 100.0),
    ("native-mutant span (limit)", 1000.0),
)

# The superoxide is born on the iron. It has three fates, all assumed rates in
# events per second: it meets the catalyst, it escapes to the blood, or it
# attacks the meso carbon. Only the ring attack is irreversible.
SUPEROXIDE_RING_ATTACK_PER_S = 1.0e5
SUPEROXIDE_ESCAPE_PER_S = 1.0e4
CATALYST_CASES: tuple[tuple[str, float], ...] = (
    ("no catalyst", 0.0),
    ("Mn-porphyrin mimic", 1.0e4),
    ("enzyme-like catalyst", 1.0e6),
    ("engineered catalyst (limit)", 1.0e8),
)

# Iron repair (ADR-0108): a one-electron donor on the empty arc. Moser-Dutton
# tunnelling gives 0.1/s at 22 A and 25/s at 18 A. Use the conservative end.
IRON_REPAIR_PER_S = 0.1

# Rotor fatigue, for the comparison (ADR-0136). A diamond MEMS rotor at the
# design speed. The passive carrier has no rotor.
ROTOR_REV_PER_S = 1000.0


def superoxide_fates(k_dismut_per_s: float,
                     k_escape_per_s: float = SUPEROXIDE_ESCAPE_PER_S,
                     k_ring_per_s: float = SUPEROXIDE_RING_ATTACK_PER_S) -> dict:
    """Split the superoxide between the catalyst, escape, and the ring."""
    total = k_dismut_per_s + k_escape_per_s + k_ring_per_s
    return {
        "k_dismut_per_s": k_dismut_per_s,
        "k_escape_per_s": k_escape_per_s,
        "k_ring_per_s": k_ring_per_s,
        "f_dismutated": k_dismut_per_s / total,
        "f_escaped": k_escape_per_s / total,
        "f_ring_attack": k_ring_per_s / total,
    }


def ring_lifetime_s(k_autox_per_s: float, f_ring_attack: float) -> float:
    """Half-life of the ring against irreversible attack."""
    k = k_autox_per_s * f_ring_attack
    return math.log(2.0) / k if k > 0.0 else float("inf")


def iron_repair(k_autox_per_s: float, k_repair_per_s: float = IRON_REPAIR_PER_S) -> dict:
    """Steady-state fraction of oxidized iron with a donor station."""
    f = k_autox_per_s / (k_autox_per_s + k_repair_per_s)
    return {
        "k_repair_per_s": k_repair_per_s,
        "f_oxidized": f,
        "f_loaded": 1.0 - f,
        "margin": k_repair_per_s / k_autox_per_s if k_autox_per_s > 0 else float("inf"),
    }


def electron_budget(sites_per_device: float, k_autox_per_s: float) -> dict:
    """Electrons and current the repair station needs, per device."""
    rate = sites_per_device * k_autox_per_s  # one electron per autoxidation
    current = rate * ELEMENTARY_CHARGE_C
    return {
        "sites_per_device": sites_per_device,
        "electrons_per_s": rate,
        "current_a": current,
        "protons_per_s": rate,  # one H+ per superoxide dismutated
    }


def rotor_lifetime_s() -> float:
    """Cycles to failure over the design speed."""
    return rs.MEMS_CYCLES_TO_FAILURE / ROTOR_REV_PER_S


def site_design(suppression: float, k_dismut_per_s: float) -> dict:
    """Ring lifetime and the fate split for one design point."""
    k_autox = central_rate_per_s() / suppression
    fates = superoxide_fates(k_dismut_per_s)
    life = ring_lifetime_s(k_autox, fates["f_ring_attack"])
    rotor = rotor_lifetime_s()
    return {
        "suppression": suppression,
        "k_dismut_per_s": k_dismut_per_s,
        "k_autox_per_s": k_autox,
        "autoxidation_half_life_h": math.log(2.0) / k_autox / HOUR_S,
        "f_ring_attack": fates["f_ring_attack"],
        "f_dismutated": fates["f_dismutated"],
        "f_escaped": fates["f_escaped"],
        "ring_lifetime_s": life,
        "ring_lifetime_days": life / rs.SECONDS_PER_DAY,
        "ring_lifetime_years": life / (DAYS_PER_YEAR * rs.SECONDS_PER_DAY),
        "outlives_rotor": life > rotor,
    }


def design_table() -> list[dict]:
    """Every pocket case against every catalyst case."""
    out = []
    for pocket_label, suppression in POCKET_CASES:
        for cat_label, k_dismut in CATALYST_CASES:
            row = site_design(suppression, k_dismut)
            row["pocket"] = pocket_label
            row["catalyst"] = cat_label
            out.append(row)
    return out


def optimize() -> dict:
    """The design that makes the site outlive the rotor, with no new material."""
    rotor = rotor_lifetime_s()
    s_best = POCKET_CASES[-1][1]
    s_central = POCKET_CASES[0][1]
    k_cat_best = CATALYST_CASES[-1][1]
    return {
        "rotor_lifetime_days": rotor / rs.SECONDS_PER_DAY,
        "rotor_lifetime_years": rotor / (DAYS_PER_YEAR * rs.SECONDS_PER_DAY),
        "pocket_only": site_design(s_best, 0.0),
        "catalyst_only": site_design(s_central, k_cat_best),
        "both": site_design(s_best, k_cat_best),
        "target_met": site_design(s_best, k_cat_best)["outlives_rotor"],
    }


def carrier_repair_budget(material: dict) -> dict:
    """Electron budget for one carrier molecule (all its sites)."""
    sites = float(material.get("o2_sites_per_molecule", 1.0)) or 1.0
    return electron_budget(sites, central_rate_per_s() / POCKET_CASES[-1][1])


def report() -> str:
    """A short text table of the design study."""
    opt = optimize()
    lines = ["Self-sustaining site repair (M29-14)", ""]
    lines.append(
        f"rotor limit: {opt['rotor_lifetime_years']:.1f} years "
        f"({rs.MEMS_CYCLES_TO_FAILURE:.0e} cycles at {ROTOR_REV_PER_S:.0f} rev/s)"
    )
    lines.append(f"central autoxidation: {math.log(2.0) / central_rate_per_s() / HOUR_S:.1f} h half-life")
    lines.append("")
    lines.append(f"{'pocket':<26}{'catalyst':<24}{'ring life':>14}{'outlives rotor':>16}")
    for row in design_table():
        life = row["ring_lifetime_years"]
        life_s = f"{life:.1f} y" if life < 1000.0 else f"{life:.0e} y"
        lines.append(
            f"{row['pocket']:<26}{row['catalyst']:<24}{life_s:>14}"
            f"{('yes' if row['outlives_rotor'] else 'no'):>16}"
        )
    lines.append("")
    lines.append("no material is consumed; the inputs are electrons and protons.")
    return "\n".join(lines)


def main() -> int:
    print(report())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
