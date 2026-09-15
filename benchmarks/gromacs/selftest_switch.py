#!/usr/bin/env python3
"""Self-test for the GROMACS switch correction.

GROMACS has no equivalent of the nanocad switching polynomial, so the cross-
check driver adds an analytic correction before it compares energies and forces.
This script checks that correction against a direct, independent evaluation of
the reference model in `crates/engine/src/reference.rs`.

It checks three things on the fixed water box:

1. An independent bonded plus switched non-bonded energy equals the potential
   energy that the Rust example prints.
2. `switch_difference` equals the plain-cutoff minus switched non-bonded energy.
3. The correction is zero when the switch is disabled.

Run it with:

    python3 benchmarks/gromacs/selftest_switch.py

It needs cargo, because it runs the Rust example. It does not need GROMACS.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import compare_water_box as drv  # noqa: E402


def minimum_image(delta_m: list[float], length_m: float) -> list[float]:
    """Reduces a pair displacement by the minimum-image convention."""
    out_m = list(delta_m)
    for axis in range(3):
        if length_m > 0.0:
            out_m[axis] -= length_m * round(out_m[axis] / length_m)
    return out_m


def switch_value(r_m: float, switch_on_m: float, cutoff_m: float) -> float:
    """Returns the nanocad switching value S(r)."""
    if r_m <= switch_on_m:
        return 1.0
    if r_m >= cutoff_m:
        return 0.0
    denom = cutoff_m * cutoff_m - switch_on_m * switch_on_m
    d = cutoff_m * cutoff_m - r_m * r_m
    e = cutoff_m * cutoff_m + 2.0 * r_m * r_m - 3.0 * switch_on_m * switch_on_m
    return d * d * e / (denom**3)


def bonded_energy_j(data: dict) -> float:
    """Returns the harmonic bond plus angle energy, computed from the JSON."""
    positions = data["positions_m"]
    bond_k = data["bond_k_n_per_m"]
    bond_r0 = data["bond_r0_m"]
    angle_k = data["angle_k_j_per_rad2"]
    angle_theta0 = data["angle_theta0_rad"]
    energy = 0.0
    for molecule in range(data["molecule_count"]):
        oxygen = 3 * molecule
        vectors = []
        for hydrogen in (oxygen + 1, oxygen + 2):
            delta = [
                positions[3 * hydrogen + axis] - positions[3 * oxygen + axis]
                for axis in range(3)
            ]
            vectors.append(delta)
            r_m = math.sqrt(sum(value * value for value in delta))
            energy += 0.5 * bond_k * (r_m - bond_r0) ** 2
        norm_one = math.sqrt(sum(value * value for value in vectors[0]))
        norm_two = math.sqrt(sum(value * value for value in vectors[1]))
        cosine = sum(a * b for a, b in zip(vectors[0], vectors[1])) / (
            norm_one * norm_two
        )
        cosine = max(-1.0, min(1.0, cosine))
        energy += 0.5 * angle_k * (math.acos(cosine) - angle_theta0) ** 2
    return energy


def nonbonded_energy_j(data: dict, switched: bool) -> float:
    """Returns the vdW plus Coulomb energy, computed from the JSON."""
    positions = data["positions_m"]
    atom_count = data["atom_count"]
    box = data["box_length_m"]
    cutoff = data["cutoff_m"]
    switch_on = data["switch_on_m"]
    ke = data["coulomb_constant_n_m2_per_c2"]
    vdw_a = data["vdw_a_j"]
    vdw_b = data["vdw_b_per_m"]
    vdw_c = data["vdw_c_j_m6"]
    charges = [
        data["oxygen_charge_c"] if atom % 3 == 0 else data["hydrogen_charge_c"]
        for atom in range(atom_count)
    ]
    excluded = {
        (3 * molecule + first, 3 * molecule + second)
        for molecule in range(data["molecule_count"])
        for first, second in ((0, 1), (0, 2), (1, 2))
    }
    energy = 0.0
    for i in range(atom_count):
        for j in range(i + 1, atom_count):
            if (i, j) in excluded:
                continue
            delta = minimum_image(
                [
                    positions[3 * i + axis] - positions[3 * j + axis]
                    for axis in range(3)
                ],
                box,
            )
            r_sq = sum(value * value for value in delta)
            if r_sq >= cutoff * cutoff:
                continue
            r_m = math.sqrt(r_sq)
            pair = ke * charges[i] * charges[j] / r_m
            if i % 3 == 0 and j % 3 == 0:
                pair += vdw_a * math.exp(-vdw_b * r_m) - vdw_c / r_m**6
            energy += pair * (switch_value(r_m, switch_on, cutoff) if switched else 1.0)
    return energy


def main() -> int:
    data = drv.run_rust_example()
    recorded = drv.read_recorded_energy_j()

    independent_j = bonded_energy_j(data) + nonbonded_energy_j(data, switched=True)
    relative_error = abs(independent_j - data["energy_j"]) / abs(data["energy_j"])
    print(f"independent bonded+nonbonded (J): {independent_j:.12e}")
    print(f"Rust example energy (J):          {data['energy_j']:.12e}")
    print(f"relative error:                   {relative_error:.3e}")
    if recorded is not None:
        recorded_error = abs(recorded - data["energy_j"]) / abs(recorded)
        print(f"recorded OpenMM result energy (J):{recorded:.12e}")
        print(f"recorded vs live relative error:  {recorded_error:.3e}")

    delta_energy_j, _ = drv.switch_difference(data)
    plain_j = nonbonded_energy_j(data, switched=False)
    switched_j = nonbonded_energy_j(data, switched=True)
    direct_delta_j = plain_j - switched_j
    delta_error_j = abs(delta_energy_j - direct_delta_j)
    print(f"switch_difference dE (J):         {delta_energy_j:.12e}")
    print(f"direct plain-minus-switched (J):  {direct_delta_j:.12e}")
    print(f"absolute difference (J):          {delta_error_j:.3e}")

    no_switch = dict(data)
    no_switch["switch_on_m"] = no_switch["cutoff_m"]
    no_switch_energy_j, no_switch_gradient = drv.switch_difference(no_switch)
    print(f"correction with switch disabled:  {no_switch_energy_j:.3e}")
    max_disabled = max(abs(value) for value in no_switch_gradient)
    print(f"max gradient with switch disabled:{max_disabled:.3e}")

    ok = (
        relative_error < 1.0e-12
        and delta_error_j < 1.0e-28
        and no_switch_energy_j == 0.0
        and all(value == 0.0 for value in no_switch_gradient)
    )
    if ok:
        print("PASS: the switch correction is consistent with the reference model.")
        return 0
    print("FAIL: the switch correction does not match the reference model.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
