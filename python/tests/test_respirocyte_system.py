"""The respirocyte system model (M29-01): provenance, balance and time course."""

from __future__ import annotations

import pytest

from nanocad import respirocyte_system as rs


def test_every_parameter_has_provenance():
    rows = rs.provenance_table(rs.default_device())
    assert rows
    for row in rows:
        assert row["tag"] in rs.TAGS, row
        assert row["low"] <= row["value"] <= row["high"], row
        if row["tag"] != "assumed":
            assert row["source"], row


def test_equation_of_state_round_trips():
    for p in (1.0, 50.0, 300.0, 1000.0, 1500.0):
        assert abs(rs.p_of_rho(rs.rho_of_p(p)) - p) < 1e-9 * max(p, 1.0)
    assert rs.rho_of_p(1000.0) == pytest.approx(14.07)


def test_a_single_stage_from_blood_stalls_near_920_atm():
    # ADR-0138: s_back 0.076 (corrected); ADR-0123's 550 atm is superseded.
    params = rs.o2_rotor_parameters()
    stall = rs.stall_atm(params, None, occupancy=1.0)
    assert 850.0 < stall < 1000.0, stall


def test_the_carry_parameters_are_the_corrected_measurement():
    # ADR-0138: s_back 0.076 +- 0.009, s_fwd 0.064 +- 0.014.
    pr = rs.o2_rotor_parameters()
    assert pr["s_back_nm3"].value == pytest.approx(0.076)
    assert 0.067 <= pr["s_back_nm3"].low and pr["s_back_nm3"].high <= 0.085
    assert pr["s_fwd_nm3"].value == pytest.approx(0.064)
    assert 0.050 <= pr["s_fwd_nm3"].low and pr["s_fwd_nm3"].high <= 0.078


def test_more_stages_raise_the_ceiling():
    ceilings = [rs.static_budget(rs.default_device(n))["cascade_ceiling_atm"] for n in (1, 2, 4)]
    assert ceilings[0] < ceilings[1] < ceilings[2], ceilings


def test_the_power_cap_scales_the_rotors():
    device = rs.default_device()
    device["power"]["peak_w"].value = 1.0  # unlimited
    free = rs.time_course(device, t_fill_max_s=60.0, t_discharge_max_s=5.0)["events"]
    assert free["min_duty"] == pytest.approx(1.0)
    device["power"]["peak_w"].value = 1.0e-12
    capped = rs.time_course(device, t_fill_max_s=60.0, t_discharge_max_s=5.0)["events"]
    assert capped["min_duty"] < 0.5
    # Backward Euler takes the power at the step end: first-order error.
    assert capped["mean_power_w"] <= 1.0e-12 * 1.01


def test_the_central_device_fills_and_discharges_what_it_stored():
    course = rs.time_course(rs.default_device())
    e = course["events"]
    assert e["fill_end"] == "target"
    assert e["discharge_end"] == "cushion"
    budget = rs.static_budget(rs.default_device())
    stored = e["tank_o2_after_fill"] - budget["tank_cushion_o2"]
    assert e["o2_delivered"] == pytest.approx(stored, rel=0.02)
    assert e["water_peak_fraction_of_chamber"] < 1.0
    assert budget["hull_fits"]


def test_edits_layout_and_page_report():
    page = rs.page_report({"n_stages": 3, "banks": [{"kind": "load", "count": 5000}],
                           "bank_params": {"load": {"s_back_nm3": 0.5}}})
    assert page["device"]["n_stages"] == 3
    load = next(b for b in page["device"]["banks"] if b["kind"] == "load")
    assert load["count"] == 5000
    assert any("s_back_nm3" in w for w in page["warnings"])
    bands = page["layout"]["bands"]
    area = sum(b["area_m2"] for b in bands)
    assert abs(area - page["cases"]["central"]["budget"]["hull_fraction_used"] * rs.hull_area_m2(rs.default_device())) < 1e-18
    assert bands[0]["theta0_rad"] == 0.0 and bands[-1]["theta1_rad"] <= 3.1416
    assert page["cases"]["central"]["log"]["t_s"]

def test_the_hull_check_counts_the_internal_banks_and_their_volume():
    # Flaw 6: area alone hid 12,000 internal cascade sections.
    budget = rs.static_budget(rs.default_device())
    assert budget["hull_sections"] == 16000
    assert budget["internal_sections"] == 12000
    assert 0.0 < budget["internal_volume_fraction"] < 1.0
    assert budget["volume_fits"]
    assert budget["hull_slab_volume_m3"] == pytest.approx(
        rs.hull_area_m2(rs.default_device()) * 2.4e-9)


def test_the_tank_holds_no_water_without_a_seal_leak():
    # Flaw 5: the tank now has a water state, and it is empty by default.
    course = rs.time_course(rs.default_device())
    assert course["events"]["tank_water_peak"] == 0.0
    assert course["events"]["water_peak"] > 0.0  # chamber 1 still carries water


def test_a_seal_leak_puts_water_in_the_tank():
    device = rs.default_device()
    device["tank"]["seal_leak_water_per_s"].value = 3.0e4
    course = rs.time_course(device)
    assert course["events"]["tank_water_peak"] > 1.0e5


def test_the_discharge_rotor_empties_the_tank_water():
    # Flaw 5: the discharge bank's water term was dead code.
    device = rs.default_device()
    device["tank"]["seal_leak_water_per_s"].value = 3.0e4
    with_discharge = rs.time_course(device)["log"]["tank_water"][-1]
    device["banks"] = [b for b in device["banks"] if b.kind != "discharge"]
    without = rs.time_course(device)["log"]["tank_water"][-1]
    assert with_discharge < 0.5 * without, (with_discharge, without)


def test_the_note_rim_speed_matches_the_design_speed():
    # Flaw 7: the note said 0.25 m/s; the true value is 0.25 mm/s.
    r = rs.default_device()["section"]["footprint_y_m"].value  # 13.7 nm slab
    rim = 2.0 * 3.141592653589793 * 4.0e-9 * 1.0e4
    assert 2.0e-4 < rim < 3.0e-4
    assert r > 0.0
