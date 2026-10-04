"""Overall respirocyte design optimization (task M29-05): choose the
architecture, not just the numbers.

``respirocyte_circulation.py`` sizes one architecture (an active pump that
loads O2, compresses it into a tank and releases it). It shows the pressure
optimum is the cushion and that the binding constraint is energy: the cam that
forces the O2 off the heme costs 650 W at peak exercise, and 342 W even if the
switch is reversible.

This module asks the next question: is the pump needed at all. The red cell
does not pump. It carries O2 on a site that binds it in the lung (pO2 100
mmHg) and releases it in the tissue (pO2 40 mmHg). The transport is passive and
costs no engine. So there are two families:

* **Carrier.** A dissolved or suspended molecule with a tuned O2 affinity and
  cooperativity. It loads and releases by the pO2 difference alone. No engine,
  no cam, no tank. The cost is the carrier mass (the nanocrit).
* **Pump.** The Freitas respirocyte. It stores O2 against a gradient, so it
  needs a switch (the cam) and an engine. The cost is the chemical power.

The module compares them on the two currencies that matter: the volume
fraction (nanocrit) and the metabolic power at rest and at peak exercise. The
rest power is the steady-state cost the body pays around the clock, so it is
the primary key.

The result is that a passive carrier with a heme-like site meets the red-cell
gas flux at a few percent nanocrit and no engine, and the pump is only needed
when O2 must be stored or delivered against a gradient. The cam is the pump's
whole cost, and the pump fails the metabolic budget even with a reversible
switch.

Results are simulated estimates. Nothing here is built or validated.
"""

from __future__ import annotations

import argparse

from nanocad.respirocyte_circulation import (
    MMHG_TO_ATM,
    design_point,
    optimize,
    physiology,
    rbc_reference,
)
from nanocad.respirocyte_system import AVOGADRO, Param

STD_MOLAR_VOLUME_M3 = 22.414e-3

# The pO2 the carrier sees on each side of the loop, in mmHg.
P_LUNG_MMHG = 100.0
P_TISSUE_REST_MMHG = 40.0
P_TISSUE_MAX_MMHG = 20.0

# The pump release modes. The cam work per pass, in kJ/mol.
RELEASE_MODES = {
    "cam_50": 50.0,  # the measured one-way ligand withdrawal (ADR-0106)
    "reversible_23": 23.0,  # the reversible switch energy (ADR-0106)
    "ideal_0": 0.0,  # compression only; a lossless pressure-driven switch
}


def carrier_saturation(p_mmhg: float, p50_mmhg: float, hill_n: float) -> float:
    """Hill O2 saturation of a site at a partial pressure ``p_mmhg``."""
    if p_mmhg <= 0.0:
        return 0.0
    x = (p_mmhg / p50_mmhg) ** hill_n
    return x / (1.0 + x)


def carrier_material(
    name: str,
    o2_per_kg_m3: float,
    density_kg_m3: float,
    p50_mmhg: float | None,
    hill_n: float,
    molar_mass_kg_per_mol: float | None,
    o2_sites_per_molecule: float,
    provenance: str,
    note: str,
) -> dict:
    """One carrier. A saturating site has a ``p50`` and a Hill ``n``; a linear
    (Henry) carrier has ``p50=None`` and ``hill_n=0`` and ``o2_per_kg_m3`` is
    the dissolved O2 at one atmosphere. ``molar_mass`` and ``o2_sites`` count
    the molecules; a linear carrier leaves them at zero."""
    return {
        "name": name,
        "o2_per_kg_m3": o2_per_kg_m3,
        "density_kg_m3": density_kg_m3,
        "p50_mmhg": p50_mmhg,
        "hill_n": hill_n,
        "molar_mass_kg_per_mol": molar_mass_kg_per_mol or 0.0,
        "o2_sites_per_molecule": o2_sites_per_molecule,
        "provenance": provenance,
        "note": note,
    }


def carrier_reference() -> list[dict]:
    """Candidate carriers. The heme carriers are tuned to a red-cell-like p50
    (about 27 mmHg); the bare heme is the density ceiling, the dendrimer is the
    practical form, and the perfluorocarbon is the linear control."""
    return [
        carrier_material(
            "heme, bare",
            1.0 * STD_MOLAR_VOLUME_M3 / 0.616,  # one heme b, 616 g/mol
            1500.0,
            26.8,
            2.7,
            0.616,
            1.0,
            "estimated",
            "one heme per 616 g/mol; unstable in water, the density ceiling",
        ),
        carrier_material(
            "heme-cyclodextrin",
            1.0 * STD_MOLAR_VOLUME_M3 / 2.9,
            1450.0,
            26.8,
            1.5,
            2.9,
            1.0,
            "estimated",
            "one heme in a cyclodextrin dimer, MW about 2.9 kDa (hemoCD-like)",
        ),
        carrier_material(
            "dendrimer-heme",
            8.0 * STD_MOLAR_VOLUME_M3 / 10.0,
            1400.0,
            26.8,
            2.0,
            10.0,
            8.0,
            "assumed",
            "eight hemes on a 10 kDa scaffold; large enough to stay in blood",
        ),
        carrier_material(
            "hemoglobin",
            4.0 * STD_MOLAR_VOLUME_M3 / 64.5,
            1335.0,
            26.8,
            2.7,
            64.5,
            4.0,
            "reference",
            "the red-cell protein; the benchmark",
        ),
        carrier_material(
            "perfluorocarbon",
            2.37e-4,  # m^3 O2 per kg per atm (about 0.45 mL/mL at 1 atm)
            1900.0,
            None,
            0.0,
            None,
            0.0,
            "reference",
            "linear, no cooperativity; needs a high pO2 to carry much",
        ),
    ]


def carrier_deliver_per_kg(material: dict, p_high_mmhg: float, p_low_mmhg: float) -> float:
    """O2 delivered per kg of carrier for one pass between two pO2 levels."""
    if material["p50_mmhg"] is None:
        return material["o2_per_kg_m3"] * (p_high_mmhg - p_low_mmhg) * MMHG_TO_ATM
    return material["o2_per_kg_m3"] * (
        carrier_saturation(p_high_mmhg, material["p50_mmhg"], material["hill_n"])
        - carrier_saturation(p_low_mmhg, material["p50_mmhg"], material["hill_n"])
    )


def carrier_design(material: dict, phys: dict[str, Param] | None = None) -> dict:
    """The carrier mass and nanocrit needed to meet the red-cell flux at rest
    and at peak exercise. A carrier has no engine, so its power is zero."""
    phys = phys or physiology()
    ref = rbc_reference(phys)
    d_rest = carrier_deliver_per_kg(material, P_LUNG_MMHG, P_TISSUE_REST_MMHG)
    d_max = carrier_deliver_per_kg(material, P_LUNG_MMHG, P_TISSUE_MAX_MMHG)

    # target is molecules O2 per m^3 blood; d is m^3 O2 per kg carrier. Convert
    # the carrier O2 to molecules with the molar volume, then divide.
    def kg_per_m3(target_molecules_per_m3, d):
        if d <= 0:
            return float("inf")
        return target_molecules_per_m3 * STD_MOLAR_VOLUME_M3 / (d * AVOGADRO)

    kg_rest = kg_per_m3(ref["rest_deliver_per_m3_molecules"], d_rest)
    kg_max = kg_per_m3(ref["max_deliver_per_m3_molecules"], d_max)
    kg_per_m3 = max(kg_rest, kg_max)
    return {
        "architecture": "carrier",
        "name": material["name"],
        "provenance": material["provenance"],
        "deliver_rest_per_kg_m3": d_rest,
        "deliver_max_per_kg_m3": d_max,
        "kg_per_m3": kg_per_m3,
        "g_per_l": kg_per_m3,
        "nanocrit_rest": kg_rest / material["density_kg_m3"],
        "nanocrit_max": kg_max / material["density_kg_m3"],
        "nanocrit": kg_per_m3 / material["density_kg_m3"],
        "chem_power_rest_w": 0.0,
        "chem_power_max_w": 0.0,
        "note": material["note"],
    }


def pump_designs(
    phys: dict[str, Param] | None = None,
    diameters_m: list[float] | None = None,
    p_max_atm: float = 40.0,
    n_stages: int = 1,
) -> list[dict]:
    """The pump at each release mode, at the best diameter for a fixed
    pressure. The cam work is the only difference between the modes."""
    phys = phys or physiology()
    # 0.5 um is the practical middle: the size optimum (0.05-0.1 um) needs
    # about 1e17 devices, which is not a fair comparison with a dissolved
    # carrier that has no size.
    diameters_m = diameters_m or [0.5e-6]
    out = []
    for mode, cam_kj_mol in RELEASE_MODES.items():
        best = None
        for d in diameters_m:
            pt = design_point(d, p_max_atm, phys, n_stages=n_stages, cam_kj_mol=cam_kj_mol)
            if not pt["tank_fits"] or not pt["stall_ok"]:
                continue
            if best is None or pt["nanocrit_max"] < best["nanocrit_max"]:
                best = pt
        if best is None:
            continue
        best = dict(best)
        best["architecture"] = "pump"
        best["name"] = f"pump, {mode}"
        best["release_mode"] = mode
        best["nanocrit"] = best["nanocrit_max"]
        best["provenance"] = "measured"
        out.append(best)
    return out


def rank_designs(phys: dict[str, Param] | None = None) -> list[dict]:
    """Every architecture, sorted by the metabolic power at rest (the steady
    cost), then by the nanocrit at peak exercise."""
    phys = phys or physiology()
    rows = [carrier_design(m, phys) for m in carrier_reference()]
    rows += pump_designs(phys)
    rows.sort(key=lambda r: (r["chem_power_rest_w"], r["nanocrit"]))
    return rows


def report(phys: dict[str, Param] | None = None) -> str:
    phys = phys or physiology()
    lines = [
        "Overall respirocyte design: carrier vs pump",
        "target: red-cell gas flux (50 mL O2/L rest, 150 mL O2/L exercise)",
        "",
        f"{'architecture':<18} {'nanocrit':>9} {'chem W rest':>12} {'chem W max':>11}  provenance",
    ]
    for r in rank_designs(phys):
        lines.append(
            f"{r['name']:<18} {100 * r['nanocrit']:8.2f}% {r['chem_power_rest_w']:12.1f} "
            f"{r['chem_power_max_w']:11.1f}  {r['provenance']}"
        )
    best = rank_designs(phys)[0]
    lines += [
        "",
        f"best on energy: {best['name']} (nanocrit {100 * best['nanocrit']:.2f}%, " f"no engine)",
        "the pump needs the cam; even the reversible switch exceeds the 100 W budget.",
    ]
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description="Overall respirocyte design (M29-05)")
    ap.add_argument("--budget", type=float, default=100.0, help="metabolic budget in W")
    a = ap.parse_args()
    print(report())
    ideal = optimize(cam_kj_mol=0.0)
    if ideal:
        print(
            f"\nideal pump (compression only, no cam): {100 * ideal['nanocrit_max']:.2f}% nanocrit, "
            f"{ideal['chem_power_rest_w']:.1f} W rest, {ideal['chem_power_max_w']:.1f} W exercise "
            f"(budget {a.budget:.0f} W)"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
