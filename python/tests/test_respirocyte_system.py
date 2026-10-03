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


def test_a_single_stage_from_blood_stalls_near_550_atm():
    params = rs.o2_rotor_parameters()
    stall = rs.stall_atm(params, None, occupancy=1.0)
    assert 450.0 < stall < 650.0, stall


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
