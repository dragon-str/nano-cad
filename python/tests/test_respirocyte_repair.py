"""Self-sustaining site repair (M29-14)."""

from __future__ import annotations

import math

import pytest

from nanocad import respirocyte_repair as rr
from nanocad import respirocyte_specs as rs


def test_central_rate_matches_the_adr_0136_pocket():
    k = rr.central_rate_per_s()
    life_h = math.log(2.0) / k / 3600.0
    assert 25.0 < life_h < 35.0


def test_the_fates_sum_to_one_and_the_catalyst_wins():
    fates = rr.superoxide_fates(1.0e6)
    assert fates["f_dismutated"] + fates["f_escaped"] + fates["f_ring_attack"] == pytest.approx(1.0)
    assert fates["f_ring_attack"] < rr.superoxide_fates(0.0)["f_ring_attack"]


def test_ring_lifetime_is_ln2_over_the_damage_rate():
    assert rr.ring_lifetime_s(1.0e-3, 1.0) == pytest.approx(math.log(2.0) * 1.0e3)
    assert rr.ring_lifetime_s(1.0e-3, 0.0) == float("inf")


def test_iron_repair_keeps_the_site_loaded():
    r = rr.iron_repair(rr.central_rate_per_s())
    assert r["f_loaded"] > 0.999
    assert r["margin"] > 1000.0


def test_the_optimized_design_outlives_the_rotor_but_the_central_does_not():
    central = rr.site_design(rr.POCKET_CASES[0][1], rr.CATALYST_CASES[2][1])
    best = rr.site_design(rr.POCKET_CASES[2][1], rr.CATALYST_CASES[2][1])
    rotor_y = rr.rotor_lifetime_s() / (rr.DAYS_PER_YEAR * rs.SECONDS_PER_DAY)
    assert not central["outlives_rotor"]
    assert best["outlives_rotor"]
    assert best["ring_lifetime_years"] > rotor_y


def test_a_better_pocket_or_catalyst_extends_the_ring_lifetime():
    base = rr.site_design(1.0, 0.0)["ring_lifetime_s"]
    better_pocket = rr.site_design(100.0, 0.0)["ring_lifetime_s"]
    better_catalyst = rr.site_design(1.0, 1.0e8)["ring_lifetime_s"]
    assert better_pocket > base
    assert better_catalyst > base


def test_the_electron_budget_is_negligible():
    b = rr.electron_budget(8.0, rr.central_rate_per_s() / 100.0)
    assert b["current_a"] < 1.0e-20
    assert b["protons_per_s"] == pytest.approx(b["electrons_per_s"])


def test_the_optimize_summary_reports_a_target_met():
    opt = rr.optimize()
    assert opt["target_met"]
    assert opt["pocket_only"]["ring_lifetime_s"] < opt["both"]["ring_lifetime_s"]
