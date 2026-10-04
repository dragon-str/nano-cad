"""The respirocyte top-down design (M29-04): circulation, targets, optimum."""

from __future__ import annotations

import pytest

from nanocad import respirocyte_circulation as rc


def test_circulation_loop_and_windows():
    c = rc.circulation()
    assert c["loop_time_rest_s"] == pytest.approx(60.0, rel=0.05)
    assert c["loop_time_max_s"] == pytest.approx(12.0, rel=0.1)
    assert c["lung_window_rest_s"] == pytest.approx(0.75)
    assert c["lung_window_max_s"] < c["lung_window_rest_s"]


def test_red_cell_targets():
    r = rc.rbc_reference()
    assert r["capacity_o2_per_l_m3"] * 1.0e6 == pytest.approx(201.0, rel=0.02)  # 201 mL O2/L
    assert r["rest_deliver_per_l_molecules"] == pytest.approx(1.34e21, rel=0.02)
    assert r["max_deliver_per_l_molecules"] == pytest.approx(3.0 * 1.34e21, rel=0.02)


def test_compression_work_is_positive_and_rises_with_pressure():
    w50 = rc.compression_j_per_molecule(0.13, 50.0)
    w1000 = rc.compression_j_per_molecule(0.13, 1000.0)
    assert 0.0 < w50 < w1000
    # Below the ideal-gas value kT ln(P/p) = 3.8e-20 J, and above 1e-20 J.
    assert 1.0e-20 < w1000 < 3.83e-20


def test_a_small_device_fits_and_delivers():
    pt = rc.design_point(0.5e-6, 100.0)
    assert pt["tank_fits"]
    assert pt["n_sections"] > 100
    assert pt["nanocrit_max"] < 0.45
    assert pt["limiting"] in ("load", "discharge", "tank")


def test_the_optimum_beats_the_freitas_size_on_nanocrit():
    small = rc.design_point(0.2e-6, 100.0)
    freitas = rc.design_point(1.0e-6, 100.0)
    assert small["nanocrit_max"] < freitas["nanocrit_max"]


def test_optimize_returns_a_feasible_design():
    best = rc.optimize()
    assert best
    assert best["tank_fits"] and best["stall_ok"]
    assert best["nanocrit_max"] < 0.45
    assert best["p_max_atm"] < 550.0  # a single stage is enough


def test_the_reversible_cam_costs_less_energy():
    hi = rc.optimize(cam_kj_mol=50.0)
    lo = rc.optimize(cam_kj_mol=23.0)
    assert lo["chem_power_max_w"] < hi["chem_power_max_w"]


def test_load_work_exceeds_discharge_work():
    pt = rc.design_point(1.0e-6, 500.0)
    assert pt["work_load_j"] > pt["work_discharge_j"]
