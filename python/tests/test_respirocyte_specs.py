"""Device specs, counts and lifetimes (M29-06)."""

from __future__ import annotations

import math

import pytest

from nanocad import respirocyte_specs as rs
from nanocad.respirocyte_design import carrier_reference, pump_designs


def test_red_cell_reference():
    r = rs.rbc_spec()
    assert r["volume_fraction"] == pytest.approx(0.45)
    assert r["lifetime_days"] == pytest.approx(120.0)
    assert r["o2_per_l_cell"] == pytest.approx(0.447, abs=0.01)


def test_autoxidation_lifetime_is_ln2_over_k():
    assert rs.autoxidation_lifetime_s(1.0e-3) == pytest.approx(math.log(2) * 1.0e3)
    assert rs.autoxidation_lifetime_s(0.0) == float("inf")


def test_replacement_matches_the_loss_over_one_day():
    life = rs.SECONDS_PER_DAY  # half-life of one day
    assert rs.replacement_per_day(100.0, life) == pytest.approx(50.0)


def test_a_dendrimer_carrier_spec_is_physical():
    dend = next(m for m in carrier_reference() if m["name"] == "dendrimer-heme")
    s = rs.carrier_spec(dend)
    assert 0.0 < s["volume_fraction"] < 0.05
    assert s["count_per_l"] > 1.0e20
    assert s["mass_body_g"] == pytest.approx(73.0, rel=0.3)
    assert s["power_body_w"] == 0.0
    # A heme carrier is far denser in O2 than a red cell.
    assert s["o2_per_l_device"] > 10.0 * rs.rbc_spec()["o2_per_l_cell"]


def test_a_pump_spec_has_a_count_power_and_tank():
    pt = next(p for p in pump_designs() if p["release_mode"] == "cam_50")
    s = rs.pump_spec(pt)
    assert s["count_body"] > 1.0e14
    assert s["power_body_w"] > 100.0
    assert s["tank_atm"] == pytest.approx(40.0)
    assert 0.0 < s["power_per_device_w"] < 1.0e-12


def test_site_lifetime_falls_as_autoxidation_rises():
    rows = rs.lifetime_table(0.073)
    lives = [r["lifetime_days"] for r in rows]
    assert lives == sorted(lives)
    assert lives[0] < 1.0 and lives[-1] > 100.0


def test_engineered_pocket_beats_the_red_cell():
    engineered = rs.AUTOXIDATION_PER_S["engineered pocket"]
    assert rs.autoxidation_lifetime_s(engineered) > rs.RBC_LIFETIME_S
