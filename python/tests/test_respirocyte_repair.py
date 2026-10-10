"""Self-sustaining site repair (M29-14) and the ring damage (M28-22)."""

from __future__ import annotations

import math

import pytest

from nanocad import respirocyte_repair as rr
from nanocad import respirocyte_specs as rs


def test_central_rate_matches_the_adr_0136_pocket():
    k = rr.central_rate_per_s()
    life_h = math.log(2.0) / k / 3600.0
    assert 25.0 < life_h < 35.0


def test_h2o2_fates_sum_to_one_and_the_catalyst_wins():
    fates = rr.h2o2_fates(1.0e7)
    assert fates["f_ring"] + fates["f_catalyzed"] == pytest.approx(1.0)
    assert fates["f_ring"] < rr.h2o2_fates(0.0)["f_ring"]
    assert rr.h2o2_fates(0.0)["f_ring"] == pytest.approx(1.0)


def test_ring_lifetime_is_ln2_over_the_damage_rate():
    assert rr.ring_lifetime_s(1.0e-3, 1.0) == pytest.approx(math.log(2.0) * 1.0e3)
    assert rr.ring_lifetime_s(1.0e-3, 0.0) == float("inf")


def test_iron_repair_keeps_the_site_loaded():
    r = rr.iron_repair(rr.central_rate_per_s())
    assert r["f_loaded"] > 0.999
    assert r["margin"] > 1000.0


def test_the_optimized_pocket_and_mimic_outlive_the_rotor_but_the_central_does_not():
    mimic = rr.CATALYST_CASES[1][1]
    central = rr.site_design(rr.POCKET_CASES[0][1], mimic)
    best = rr.site_design(rr.POCKET_CASES[2][1], mimic)
    rotor_y = rr.rotor_lifetime_s() / (rr.DAYS_PER_YEAR * rs.SECONDS_PER_DAY)
    assert not central["outlives_rotor"]
    assert best["outlives_rotor"]
    assert best["ring_lifetime_years"] > rotor_y


def test_a_better_pocket_or_catalyst_extends_the_ring_lifetime():
    base = rr.site_design(1.0, 0.0)["ring_lifetime_s"]
    better_pocket = rr.site_design(100.0, 0.0)["ring_lifetime_s"]
    better_catalyst = rr.site_design(1.0, 1.0e7)["ring_lifetime_s"]
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


def test_the_ring_attack_references_are_well_formed():
    names = [name for name, _, _ in rr.RING_ATTACK_REFERENCES]
    assert len(names) == len(set(names))
    for name, rate, source in rr.RING_ATTACK_REFERENCES:
        assert rate > 0.0, name
        assert source


def test_the_catalyst_range_lengthens_the_ring_lifetime():
    sc = rr.ring_damage_scenarios(rr.POCKET_CASES[-2][1])
    assert sc[0]["scenario"] == "no catalyst"
    life = [s["ring_lifetime_years"] for s in sc]
    assert life == sorted(life)
    assert sc[0]["outlives_rotor"] is False
    assert sc[-1]["outlives_rotor"] is True


def test_the_no_catalyst_scenario_matches_the_design_table():
    no_cat = rr.ring_damage_scenarios(rr.POCKET_CASES[-2][1])[0]
    row = rr.site_design(rr.POCKET_CASES[-2][1], 0.0)
    assert no_cat["ring_lifetime_years"] == pytest.approx(row["ring_lifetime_years"])
