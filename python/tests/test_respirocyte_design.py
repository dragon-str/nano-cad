"""The overall respirocyte design: carrier vs pump (M29-05)."""

from __future__ import annotations

import pytest

from nanocad import respirocyte_design as rd


def test_hill_saturation_is_bounded_and_monotone():
    assert rd.carrier_saturation(0.0, 26.8, 2.7) == 0.0
    s40 = rd.carrier_saturation(40.0, 26.8, 2.7)
    s100 = rd.carrier_saturation(100.0, 26.8, 2.7)
    assert 0.0 < s40 < s100 < 1.0
    assert rd.carrier_saturation(1e9, 26.8, 2.7) == pytest.approx(1.0, abs=1e-6)


def test_hemoglobin_carrier_matches_the_real_nanocrit():
    hb = next(m for m in rd.carrier_reference() if m["name"] == "hemoglobin")
    d = rd.carrier_design(hb)
    # 15 g/dL Hb in whole blood is about a 12% volume fraction of protein.
    assert d["nanocrit"] == pytest.approx(0.12, abs=0.03)


def test_a_denser_carrier_needs_less_volume():
    by_name = {m["name"]: m for m in rd.carrier_reference()}
    bare = rd.carrier_design(by_name["heme, bare"])["nanocrit"]
    dend = rd.carrier_design(by_name["dendrimer-heme"])["nanocrit"]
    hb = rd.carrier_design(by_name["hemoglobin"])["nanocrit"]
    assert bare < dend < hb


def test_perfluorocarbon_cannot_carry_the_flux_at_room_air():
    pfc = next(m for m in rd.carrier_reference() if m["name"] == "perfluorocarbon")
    assert rd.carrier_design(pfc)["nanocrit"] > 1.0


def test_carriers_use_no_engine():
    for m in rd.carrier_reference():
        d = rd.carrier_design(m)
        assert d["chem_power_rest_w"] == 0.0 and d["chem_power_max_w"] == 0.0


def test_pump_release_modes_are_ordered_by_energy():
    rows = {r["name"]: r for r in rd.pump_designs()}
    ideal = rows["pump, ideal_0"]
    rev = rows["pump, reversible_23"]
    cam = rows["pump, cam_50"]
    assert ideal["chem_power_max_w"] < rev["chem_power_max_w"] < cam["chem_power_max_w"]
    assert ideal["chem_power_rest_w"] < rev["chem_power_rest_w"] < cam["chem_power_rest_w"]


def test_the_ideal_pump_fits_the_budget_but_the_cam_pump_does_not():
    rows = {r["name"]: r for r in rd.pump_designs()}
    assert rows["pump, ideal_0"]["chem_power_max_w"] < 100.0
    assert rows["pump, cam_50"]["chem_power_max_w"] > 100.0


def test_the_best_architecture_is_a_carrier():
    best = rd.rank_designs()[0]
    assert best["architecture"] == "carrier"
    assert best["nanocrit"] < 0.05
