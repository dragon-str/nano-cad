"""A self-sustaining repair strategy for the heme carrier site (M29-14, M28-22).

The site degrades in two ways (ADR-0108, ADR-0136):

  * the iron oxidizes, Fe(II) -> Fe(III), and stops binding O2. This is
    reversible: a one-electron donor on the empty arc re-reduces it. The
    electrons are trivial and come from the device, so no material is used.
  * the porphyrin ring is damaged. This is irreversible: a broken ring cannot
    be rebuilt without new carbon.

So a self-sustaining site cannot repair the ring. It must stop the ring from
breaking. Three permanent elements do that:

  1. a distal pocket that suppresses autoxidation (a low superoxide rate);
  2. a co-located catalyst that dismutates the superoxide and then decomposes
     the peroxide before it reaches the ring;
  3. the iron-repair donor, so the site stays in the Fe(II) state.

None of the three is consumed. The only inputs are electrons (a charged
onboard store) and protons (buffered blood). The alternative, a consumable
antioxidant, would need new material each day and is rejected.

Mechanism of the ring damage (M28-22). The autoxidation superoxide does not
attack the ring directly. It converts almost completely into hydrogen peroxide
by dismutation (Wazawa 1992, PMID 1420896). The H2O2 is the oxidant: it reacts
with the Fe(II) heme to form a ferryl, which opens the ring (Yusa & Shikama
1987, PMID 3427037). The ring lifetime is set by the race for the H2O2 between
the heme and the co-located catalyst.

This module estimates the ring lifetime from that race, and finds the design
that outlives the rotor. All rates here are estimates or measured proxy values;
the model shows the design, not a measurement of this pocket.
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

# The ring damage is mediated by H2O2 (M28-22), so the rate is bounded by the
# reaction of the heme with H2O2. The measured rate of deoxy-myoglobin
# (Fe(II)) + H2O2 -> ferryl is 3.6e3 /M/s (Yusa & Shikama 1987, PMID 3427037).
# This is an upper bound on the ring-damage rate: the ferryl can return to the
# resting state instead of opening the ring, and the superoxide can escape
# before it dismutates, both of which lower the damage.
HEME_H2O2_PER_M_S = 3.6e3

# The catalyst's catalase-like rate (H2O2 -> O2 + H2O), taken as a
# pseudo-first-order rate with the co-located catalyst at one molar. A native
# catalase is one of the fastest enzymes (k_cat/K_M about 1e7 /M/s); a
# Mn-porphyrin or Mn-salen dual SOD/catalase mimic is slower, about 1e5 /M/s.
# (label, rate).
CATALYST_CASES: tuple[tuple[str, float], ...] = (
    ("no catalyst", 0.0),
    ("Mn-porphyrin mimic", 1.0e5),
    ("enzyme-like (catalase)", 1.0e7),
    ("engineered catalyst (limit)", 1.0e9),
)

# The literature basis for the ring damage (M28-22). Each entry is
# (name, rate_per_m_per_s, source).
RING_ATTACK_REFERENCES: tuple[tuple[str, float, str], ...] = (
    ("heme + H2O2 -> ferryl", HEME_H2O2_PER_M_S, "Yusa & Shikama 1987, PMID 3427037"),
    ("native catalase + H2O2", 1.0e7, "catalase k_cat/K_M, one of the fastest enzymes"),
    ("Mn-porphyrin catalase-like", 1.0e5, "Mn porphyrin/salen dual SOD/catalase mimic"),
    ("superoxide -> H2O2", float("inf"), "Wazawa 1992, PMID 1420896 (almost complete)"),
)

# Iron repair (ADR-0108): a one-electron donor on the empty arc. Moser-Dutton
# tunnelling gives 0.1/s at 22 A and 25/s at 18 A. Use the conservative end.
IRON_REPAIR_PER_S = 0.1

# The repair rate falls with the donor-iron distance (Moser-Dutton tunnelling).
# ADR-0108 fixes the slope with two points: 22 A repairs 0.1/s and 18 A repairs
# 25/s, so beta = ln(250)/4 A = 1.38 /A. IRON_REPAIR_PER_S keeps the far 22 A
# value as the model's worst case. The CAD donor is much closer (ADR-0142), so
# the modelled margin is a floor, not the design value.
MOSER_DUTTON_BETA_PER_M = 1.38e10
DONOR_REFERENCE_M = 22.0e-10

# The donor-iron distance in each CAD build (ADR-0139, ADR-0141).
DONOR_CAD_DISTANCES: tuple[tuple[str, float], ...] = (
    ("free carrier site (ADR-0139)", 3.61e-10),
    ("device oxygen rotor (ADR-0141)", 7.4e-10),
)

# Rotor fatigue, for the comparison (ADR-0136). A diamond MEMS rotor at the
# design speed. The passive carrier has no rotor.
ROTOR_REV_PER_S = 1000.0


def donor_repair_rate(
    distance_m: float, reference_rate_per_s: float = IRON_REPAIR_PER_S
) -> float:
    """Moser-Dutton repair rate at a donor-iron distance (per second)."""
    return reference_rate_per_s * math.exp(
        MOSER_DUTTON_BETA_PER_M * (DONOR_REFERENCE_M - distance_m)
    )


def donor_table(k_autox_per_s: float) -> list[dict]:
    """Repair rate and margin at each CAD donor distance, against one rate."""
    out = []
    for label, distance in DONOR_CAD_DISTANCES:
        rate = donor_repair_rate(distance)
        out.append(
            {
                "site": label,
                "donor_distance_m": distance,
                "k_repair_per_s": rate,
                "margin": rate / k_autox_per_s if k_autox_per_s > 0 else float("inf"),
            }
        )
    return out


def h2o2_fates(k_cat_per_s: float) -> dict:
    """Split the pocket H2O2 between the heme (ring damage) and the catalyst."""
    total = HEME_H2O2_PER_M_S + k_cat_per_s
    return {
        "k_heme_per_s": HEME_H2O2_PER_M_S,
        "k_cat_per_s": k_cat_per_s,
        "f_ring": HEME_H2O2_PER_M_S / total,
        "f_catalyzed": k_cat_per_s / total,
    }


def ring_lifetime_s(k_autox_per_s: float, f_ring: float) -> float:
    """Half-life of the ring against irreversible attack."""
    k = k_autox_per_s * f_ring
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


def ring_attack_rate_per_s(suppression: float, k_cat_per_s: float) -> float:
    """The irreversible ring-damage rate for one design point."""
    k_autox = central_rate_per_s() / suppression
    return k_autox * h2o2_fates(k_cat_per_s)["f_ring"]


def site_design(suppression: float, k_cat_per_s: float) -> dict:
    """Ring lifetime and the H2O2 split for one design point."""
    k_autox = central_rate_per_s() / suppression
    fates = h2o2_fates(k_cat_per_s)
    life = ring_lifetime_s(k_autox, fates["f_ring"])
    rotor = rotor_lifetime_s()
    return {
        "suppression": suppression,
        "k_cat_per_s": k_cat_per_s,
        "k_autox_per_s": k_autox,
        "autoxidation_half_life_h": math.log(2.0) / k_autox / HOUR_S,
        "f_ring": fates["f_ring"],
        "f_catalyzed": fates["f_catalyzed"],
        "ring_attack_per_s": k_autox * fates["f_ring"],
        "ring_lifetime_s": life,
        "ring_lifetime_days": life / rs.SECONDS_PER_DAY,
        "ring_lifetime_years": life / (DAYS_PER_YEAR * rs.SECONDS_PER_DAY),
        "outlives_rotor": life > rotor,
    }


def design_table() -> list[dict]:
    """Every pocket case against every catalyst case."""
    out = []
    for pocket_label, suppression in POCKET_CASES:
        for cat_label, k_cat in CATALYST_CASES:
            row = site_design(suppression, k_cat)
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


def ring_damage_scenarios(suppression: float) -> list[dict]:
    """The ring lifetime across the sourced catalyst range (M28-22).

    The heme + H2O2 rate is measured, so the only range is the catalyst. The
    no-catalyst case is the upper bound on the damage.
    """
    rotor = rotor_lifetime_s()
    out = []
    for label, k_cat in CATALYST_CASES:
        row = site_design(suppression, k_cat)
        out.append(
            {
                "scenario": label,
                "k_cat_per_s": k_cat,
                "f_ring": row["f_ring"],
                "ring_lifetime_years": row["ring_lifetime_years"],
                "outlives_rotor": row["ring_lifetime_years"] * DAYS_PER_YEAR * rs.SECONDS_PER_DAY > rotor,
            }
        )
    return out


def carrier_repair_budget(material: dict) -> dict:
    """Electron budget for one carrier molecule (all its sites)."""
    sites = float(material.get("o2_sites_per_molecule", 1.0)) or 1.0
    return electron_budget(sites, central_rate_per_s() / POCKET_CASES[-1][1])


def report() -> str:
    """A short text table of the design study."""
    opt = optimize()
    lines = ["Self-sustaining site repair (M29-14, ring damage M28-22)", ""]
    lines.append(
        f"rotor limit: {opt['rotor_lifetime_years']:.1f} years "
        f"({rs.MEMS_CYCLES_TO_FAILURE:.0e} cycles at {ROTOR_REV_PER_S:.0f} rev/s)"
    )
    lines.append(f"central autoxidation: {math.log(2.0) / central_rate_per_s() / HOUR_S:.1f} h half-life")
    lines.append(f"heme + H2O2 -> ferryl: {HEME_H2O2_PER_M_S:.1e} /M/s (measured, Yusa & Shikama 1987)")
    lines.append("")
    lines.append("iron repair rate against the donor-iron distance (Moser-Dutton, beta 1.38 /A):")
    lines.append(f"{'site':<32}{'d(Fe-donor)':>13}{'k_repair /s':>14}{'margin':>12}")
    for row in donor_table(central_rate_per_s()):
        lines.append(
            f"{row['site']:<32}{row['donor_distance_m'] * 1e10:>10.2f} A"
            f"{row['k_repair_per_s']:>14.1e}{row['margin']:>12.1e}"
        )
    lines.append(
        f"the model keeps the far 22 A reference ({IRON_REPAIR_PER_S:.1f} /s) as its worst case"
    )
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
    lines.append("sourced catalyst range (optimized pocket, 100x suppression):")
    lines.append(f"{'scenario':<26}{'k_cat /M/s':>12}{'f_ring':>12}{'ring life':>14}{'outlives':>10}")
    for sc in ring_damage_scenarios(POCKET_CASES[-2][1]):
        life = sc["ring_lifetime_years"]
        life_s = f"{life:.2f} y" if life < 1000.0 else f"{life:.0e} y"
        lines.append(
            f"{sc['scenario']:<26}{sc['k_cat_per_s']:>12.0e}{sc['f_ring']:>12.2e}"
            f"{life_s:>14}{('yes' if sc['outlives_rotor'] else 'no'):>10}"
        )
    lines.append("")
    lines.append("no material is consumed; the inputs are electrons and protons.")
    return "\n".join(lines)


def main() -> int:
    print(report())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
