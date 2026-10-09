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
    # The measured references span seconds (no distal H-bond) to about 6 days.
    assert lives[0] < 1.0 and 5.0 < lives[-1] < 7.0


def test_every_reference_rate_is_sourced_and_positive():
    for name, k, temp_c, pmid, _ in rs.AUTOXIDATION_REFERENCES:
        assert k > 0, name
        assert pmid.isdigit(), name
        assert 20.0 <= temp_c <= 40.0, name


def test_no_measured_heme_carrier_reaches_the_target_or_the_red_cell():
    slowest = min(rs.AUTOXIDATION_PER_S.values())
    assert slowest > rs.AUTOXIDATION_TARGET_PER_S
    assert rs.autoxidation_lifetime_s(slowest) < rs.RBC_LIFETIME_S


def test_native_myoglobin_matches_brantley():
    # 0.055 per hour at 37 C (PubMed 8463233): a half-life of about 12.6 h.
    k = rs.AUTOXIDATION_PER_S["myoglobin, native"]
    assert rs.autoxidation_lifetime_s(k) / 3600.0 == pytest.approx(12.6, rel=0.01)


def test_two_channel_fit_reproduces_hemoglobin():
    fit = rs.fit_two_channel()
    assert fit["k_ox_per_s"] > 0 and fit["k_bi_per_s_mmhg"] > 0
    assert fit["max_rel_residual"] < 0.10
    # The oxy channel alone must agree with HbA in air within a factor of 2.
    assert 0.5 < fit["hba_air_check_ratio"] < 2.0
    for p, t_h in rs.HB_PLASMA_HALF_LIFE_H:
        k = rs.autoxidation_rate(p, rs.HB_P50_MMHG, rs.HB_HILL_N, fit=fit)
        assert math.log(2) / k / 3600.0 == pytest.approx(t_h, rel=0.10)


def test_autoxidation_is_faster_in_venous_blood():
    art = rs.autoxidation_rate(95.0, rs.HB_P50_MMHG, rs.HB_HILL_N)
    ven = rs.autoxidation_rate(40.0, rs.HB_P50_MMHG, rs.HB_HILL_N)
    assert ven > art


def test_pocket_cases_are_ordered_and_bracket_the_central_value():
    cases = {r["case"]: r for r in rs.pocket_autoxidation(26.8, 2.0)}
    assert cases["optimistic"]["k_per_s"] < cases["central"]["k_per_s"] < cases["pessimistic"]["k_per_s"]
    assert cases["central"]["k_per_s"] < cases["hemoglobin-like pocket"]["k_per_s"]


def test_a_higher_affinity_site_lasts_longer():
    def central(p50):
        return next(r for r in rs.pocket_autoxidation(p50, 2.0) if r["case"] == "central")["half_life_h"]

    assert central(10.0) > central(26.8) > central(40.0)


def test_carrier_lifetime_comes_from_the_pocket_model():
    dend = next(m for m in rs.carrier_reference() if m["name"] == "dendrimer-heme")
    spec = rs.carrier_spec(dend)
    central = next(r for r in rs.pocket_autoxidation(26.8, 2.0) if r["case"] == "central")
    assert spec["site_lifetime_days"] * 24.0 == pytest.approx(central["half_life_h"], rel=1e-9)
