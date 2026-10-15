"""Seal leak rates for the respirocyte (M28-36).

The seal audit (ADR-0118, ADR-0121, ADR-0122) gives the *barriers* of the
blood-to-tank paths: for the O2 rotor, O2 146 kJ/mol and water 57 kJ/mol; for
the water rotor, O2 200 kJ/mol and water 149 kJ/mol. A barrier is not a rate.
This module turns each barrier into a leak rate.

Model
-----
A molecule crosses a barrier of height dG when it strikes the leak area A with
enough energy. This is the Hertz-Knudsen permeation rate:

    J = (1/4) n <v> A exp(-dG / RT)                        [molecules / s]

where n is the number density on the high-pressure side and <v> the mean
thermal speed sqrt(8 R T / (pi M)). The exponential is the Boltzmann fraction
of the striking molecules that pass; the impingement rate (1/4) n <v> is the
number that strike unit area each second.

What this is and what it is not
-------------------------------
This is an estimate. The barrier is a static Lennard-Jones value (a one-site
probe on a frozen scene). It carries no entropy term; a free-energy profile
from a PMF would add one, and would be the better number. The leak area A is a
design input, not a measurement. The rate is exponential in dG, so a 10 kJ/mol
error in the barrier moves the rate by a factor of about 50. Treat a rate here
as an order of magnitude, not a value.

The bulk impingement rate also over-states the flux into a channel only a few
angstroms wide, so each rate here is an upper bound. Simulated estimates.
"""

from __future__ import annotations

import math

from nanocad.respirocyte_system import AVOGADRO, EOS_ATM_RHO

GAS_CONSTANT_KJ_PER_MOL_K = 0.0083145
T_K = 310.0
KT_KJ_PER_MOL = GAS_CONSTANT_KJ_PER_MOL_K * T_K

# O2 0.032 kg/mol, water 0.018 kg/mol.
MOLAR_MASS_KG_PER_MOL = {"O2": 0.032, "water": 0.018}

# Number density on the high-pressure side, per m^3.
O2_PER_NM3_1000ATM = dict(EOS_ATM_RHO)[1000.0]  # ADR-0123
O2_DENSITY_PER_M3 = O2_PER_NM3_1000ATM * 1.0e27
WATER_DENSITY_PER_M3 = 55.5 * 1.0e3 * AVOGADRO  # blood, ~55.5 mol/L
DENSITY_PER_M3 = {"O2": O2_DENSITY_PER_M3, "water": WATER_DENSITY_PER_M3}

# Blood-to-tank seal barriers (kJ/mol). ADR-0121 (O2 rotor), ADR-0122 (water).
SEAL_BARRIERS_KJ = (
    ("O2 rotor", "O2", 146.0),
    ("O2 rotor", "water", 57.0),
    ("water rotor", "O2", 200.0),
    ("water rotor", "water", 149.0),
)

# The unit leak area the table is reported against: 1 nm^2.
LEAK_AREA_M2 = 1.0e-18


def mean_speed_m_per_s(molecule: str, t_k: float = T_K) -> float:
    """The mean thermal speed, sqrt(8 R T / (pi M))."""
    return math.sqrt(8.0 * GAS_CONSTANT_KJ_PER_MOL_K * 1000.0 * t_k / (math.pi * MOLAR_MASS_KG_PER_MOL[molecule]))


def impingement_per_m2_s(molecule: str, t_k: float = T_K) -> float:
    """The Hertz-Knudsen impingement rate, (1/4) n <v>."""
    return 0.25 * DENSITY_PER_M3[molecule] * mean_speed_m_per_s(molecule, t_k)


def boltzmann_fraction(barrier_kj_per_mol: float, t_k: float = T_K) -> float:
    """The fraction of strikes that pass the barrier, exp(-dG / RT)."""
    return math.exp(-barrier_kj_per_mol / (GAS_CONSTANT_KJ_PER_MOL_K * t_k))


def leak_rate_per_s(
    barrier_kj_per_mol: float,
    molecule: str,
    area_m2: float = LEAK_AREA_M2,
    t_k: float = T_K,
) -> float:
    """The permeation rate through a barrier of the given height."""
    return impingement_per_m2_s(molecule, t_k) * area_m2 * boltzmann_fraction(barrier_kj_per_mol, t_k)


def leak_table(area_m2: float = LEAK_AREA_M2) -> list[dict]:
    """One row per seal and molecule."""
    rows = []
    for rotor, molecule, barrier in SEAL_BARRIERS_KJ:
        rate = leak_rate_per_s(barrier, molecule, area_m2)
        rows.append(
            {
                "rotor": rotor,
                "molecule": molecule,
                "barrier_kj_per_mol": barrier,
                "barrier_kt": barrier / KT_KJ_PER_MOL,
                "area_m2": area_m2,
                "leak_per_s": rate,
            }
        )
    return rows


def o2_leak_safety_factor() -> float:
    """How much slower the O2 leak is than the water leak on the O2 rotor.

    The O2 leak is the one that loses the product. The water leak is a load on
    the water-return rotor, not a loss (ADR-0122 returns about 5.9 waters per
    pocket pass).
    """
    o2 = leak_rate_per_s(146.0, "O2")
    water = leak_rate_per_s(57.0, "water")
    return water / o2


def report() -> str:
    """The M28-36 seal-leak estimate as text."""
    lines = [
        "Seal leak rates (M28-36)",
        "",
        f"temperature {T_K:.0f} K, RT {KT_KJ_PER_MOL:.2f} kJ/mol, leak area {LEAK_AREA_M2 * 1e18:.0f} nm^2",
        "Hertz-Knudsen permeation, (1/4) n <v> A exp(-dG/RT); static barrier, no entropy term",
        "",
        f"{'rotor':<14}{'molecule':<10}{'barrier':>10}{'kT':>7}{'leak /s':>12}",
    ]
    for row in leak_table():
        lines.append(
            f"{row['rotor']:<14}{row['molecule']:<10}"
            f"{row['barrier_kj_per_mol']:>8.0f} kJ{row['barrier_kt']:>7.0f}"
            f"{row['leak_per_s']:>12.1e}"
        )
    lines.append("")
    lines.append(f"O2 leak is {o2_leak_safety_factor():.0e}x slower than the water leak")
    lines.append("verdict: the O2 barrier (57 kT) holds, so the product does not leak. The")
    lines.append("water barrier (22 kT) is the soft one; the water-return rotor must match")
    lines.append("it (ADR-0122 returns about 6 waters per pocket pass). These are upper")
    lines.append("bounds; a PMF would add the entropy term, and the balance is a follow-up.")
    return "\n".join(lines)


def main() -> int:
    print(report())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
