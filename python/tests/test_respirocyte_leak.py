"""Seal leak rates (M28-36)."""

from __future__ import annotations

import math

import pytest

from nanocad import respirocyte_leak as rl


def test_rt_at_body_temperature_is_about_2_6_kj_per_mol():
    assert 2.5 < rl.KT_KJ_PER_MOL < 2.7


def test_the_mean_speed_of_oxygen_is_about_450_m_per_s():
    assert 400.0 < rl.mean_speed_m_per_s("O2") < 500.0


def test_the_impingement_rate_is_density_times_a_quarter_of_the_speed():
    speed = rl.mean_speed_m_per_s("water")
    assert rl.impingement_per_m2_s("water") == pytest.approx(0.25 * rl.WATER_DENSITY_PER_M3 * speed)


def test_the_barrier_heights_match_the_audit():
    barriers = {(rotor, mol): kj for rotor, mol, kj in rl.SEAL_BARRIERS_KJ}
    assert barriers[("O2 rotor", "O2")] == 146.0
    assert barriers[("O2 rotor", "water")] == 57.0
    assert barriers[("water rotor", "O2")] == 200.0
    assert barriers[("water rotor", "water")] == 149.0


def test_the_leak_table_covers_both_rotors_and_both_molecules():
    rows = rl.leak_table()
    assert len(rows) == 4
    assert {r["rotor"] for r in rows} == {"O2 rotor", "water rotor"}
    assert {r["molecule"] for r in rows} == {"O2", "water"}


def test_the_o2_leak_through_the_o2_rotor_seal_is_negligible():
    # ADR-0121 measures a 146 kJ/mol O2 barrier (57 kT).
    rate = rl.leak_rate_per_s(146.0, "O2")
    assert rate < 1.0e-9


def test_the_water_leak_is_much_faster_than_the_o2_leak():
    water = rl.leak_rate_per_s(57.0, "water")
    o2 = rl.leak_rate_per_s(146.0, "O2")
    assert water > 1.0e2
    assert rl.o2_leak_safety_factor() > 1.0e10


def test_a_ten_kj_barrier_moves_the_rate_by_more_than_an_order_of_magnitude():
    fast = rl.leak_rate_per_s(146.0, "O2")
    slow = rl.leak_rate_per_s(156.0, "O2")
    assert fast / slow > 10.0
    assert math.exp(10.0 / rl.KT_KJ_PER_MOL) == pytest.approx(fast / slow)


def test_the_report_names_both_rotors():
    text = rl.report()
    assert "O2 rotor" in text and "water rotor" in text
    assert "Hertz-Knudsen" in text
