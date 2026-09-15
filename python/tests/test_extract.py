"""Tests for the parameter extractors in the ``nanocad._core`` extension."""

from __future__ import annotations

import math

import pytest

import nanocad

R0_M = 1.5e-10
K_N_PER_M = 300.0
CARBON_MASS_KG = 1.99264687992e-26
AREA_M2 = 1.0e-20
BOND_COUNT = 20


def uniform_chain(
    atom_count: int = BOND_COUNT + 1,
) -> tuple[nanocad.System, list[float]]:
    """Build a straight harmonic chain along x."""
    system = nanocad.System(atom_count, [CARBON_MASS_KG] * atom_count)
    system.set_bond_stretch(
        [(i, i + 1, K_N_PER_M, R0_M) for i in range(atom_count - 1)]
    )
    positions: list[float] = []
    for atom in range(atom_count):
        positions += [atom * R0_M, 0.0, 0.0]
    return system, positions


def test_extract_stiffness_recovers_the_chain_modulus() -> None:
    system, positions = uniform_chain()
    result = nanocad.extract_stiffness(
        system,
        positions,
        reference_length_m=BOND_COUNT * R0_M,
        cross_section_area_m2=AREA_M2,
        derivative_step=1.0e-5,
    )
    expected_pa = K_N_PER_M * R0_M / AREA_M2
    assert result.elastic_modulus_pa.value_si == pytest.approx(expected_pa, rel=1.0e-6)
    assert result.r_squared > 0.999999
    assert result.samples == 11
    assert result.provenance.method == "fit"


def test_extract_stiffness_rejects_a_bad_sweep() -> None:
    system, positions = uniform_chain(3)
    with pytest.raises(ValueError):
        nanocad.extract_stiffness(system, positions, strain_max=-1.0)


def test_summarize_friction_profiles_a_sinusoid() -> None:
    samples = [math.sin(2.0 * math.pi * index / 500.0) for index in range(2000)]
    result = nanocad.summarize_friction(samples, 1.0e-10)
    assert result.samples == 2000
    assert result.friction_force_mean_n > 0.0
    assert result.friction_force_spread_n > 0.0
    assert result.friction_coefficient.value_si > 0.0
    assert result.provenance.method == "fit"


def test_extract_failure_stress_from_system_finds_the_rupture() -> None:
    system, positions = uniform_chain()
    result = nanocad.extract_failure_stress_from_system(
        system,
        positions,
        0.05,
        reference_length_m=BOND_COUNT * R0_M,
        cross_section_area_m2=AREA_M2,
        derivative_step=1.0e-5,
    )
    assert result.failure_strain == pytest.approx(0.05, abs=1.0e-9)
    assert 0 <= result.critical_bond < BOND_COUNT
    assert result.failure_stress_pa.value_si > 0.0
    assert result.provenance.method == "fit"


def test_extract_failure_stress_with_explicit_bonds() -> None:
    system, positions = uniform_chain()
    bonds = [(i, i + 1, R0_M, 0.05) for i in range(BOND_COUNT)]
    result = nanocad.extract_failure_stress(
        system,
        positions,
        bonds,
        reference_length_m=BOND_COUNT * R0_M,
        cross_section_area_m2=AREA_M2,
        derivative_step=1.0e-5,
    )
    assert result.failure_strain == pytest.approx(0.05, abs=1.0e-9)


def test_extract_friction_rejects_a_bad_slider() -> None:
    system, positions = uniform_chain()
    with pytest.raises(ValueError):
        nanocad.extract_friction(system, positions, slider_atom=999)


def test_thermal_extractors_are_reachable_and_validate() -> None:
    assert callable(nanocad.extract_specific_heat)
    assert callable(nanocad.extract_thermal_conductivity)
    assert callable(nanocad.extract_thermal)
    system, _positions = uniform_chain()
    with pytest.raises(ValueError):
        nanocad.extract_specific_heat(system, [0.0])
