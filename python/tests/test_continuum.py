"""Tests for the continuum adapter.

Every quantity is SI.  Lengths are in metres [m], areas in square metres
[m^2], pressure and modulus in pascal [Pa], force in newton [N], density in
[kg/m^3], viscosity in [Pa s], and velocity in [m/s].  Variable names carry a
unit suffix.
"""

from __future__ import annotations

import math

import pytest

np = pytest.importorskip("numpy")

from nanocad import continuum_adapter as ca  # noqa: E402

BAR_LENGTH_M = 1.0
BAR_AREA_M2 = 1.0e-4
STEEL_MODULUS_PA = 200.0e9
BAR_LOAD_N = 1000.0
STEEL_DENSITY_KG_M3 = 7850.0

BEAM_LENGTH_M = 1.0
BEAM_WIDTH_M = 0.05
BEAM_HEIGHT_M = 0.1
BEAM_INERTIA_M4 = BEAM_WIDTH_M * BEAM_HEIGHT_M**3 / 12.0
BEAM_OUTER_FIBER_M = BEAM_HEIGHT_M / 2.0
BEAM_TIP_LOAD_N = 100.0

FLOW_LENGTH_M = 0.1
WATER_VISCOSITY_PA_S = 1.0e-3
PRESSURE_DROP_PA = 100.0
TUBE_RADIUS_M = 1.0e-3
CHANNEL_HALF_HEIGHT_M = 5.0e-4
CHANNEL_WIDTH_M = 1.0e-2

HIGH_BINS = 131072


def test_solve_axial_bar_matches_the_closed_form() -> None:
    result = ca.solve_axial_bar(
        BAR_LENGTH_M,
        BAR_AREA_M2,
        STEEL_MODULUS_PA,
        BAR_LOAD_N,
        STEEL_DENSITY_KG_M3,
    )
    analytic_m = BAR_LOAD_N * BAR_LENGTH_M / (BAR_AREA_M2 * STEEL_MODULUS_PA)
    assert result["analytic_tip_displacement_m"] == pytest.approx(analytic_m, rel=1e-15)
    assert result["tip_displacement_m"] == pytest.approx(analytic_m, rel=1.0e-10)
    assert result["relative_error"] < 1.0e-10
    assert result["method"] == "axial_bar_linear_fem"
    assert result["validation"] == "cross-checked"


def test_axial_bar_mesh_is_exact_at_every_refinement() -> None:
    analytic_m = BAR_LOAD_N * BAR_LENGTH_M / (BAR_AREA_M2 * STEEL_MODULUS_PA)
    for elements in (1, 2, 5, 20, 50):
        result = ca.solve_axial_bar(
            BAR_LENGTH_M,
            BAR_AREA_M2,
            STEEL_MODULUS_PA,
            BAR_LOAD_N,
            STEEL_DENSITY_KG_M3,
            elements=elements,
        )
        assert result["elements"] == elements
        assert result["tip_displacement_m"] == pytest.approx(analytic_m, rel=1.0e-10)
        assert result["relative_error"] < 1.0e-10


def test_axial_bar_stress_and_strain_are_uniform() -> None:
    result = ca.solve_axial_bar(
        BAR_LENGTH_M,
        BAR_AREA_M2,
        STEEL_MODULUS_PA,
        BAR_LOAD_N,
        STEEL_DENSITY_KG_M3,
    )
    expected_stress_pa = BAR_LOAD_N / BAR_AREA_M2
    assert result["axial_stress_pa"] == pytest.approx(expected_stress_pa, rel=1.0e-12)
    assert result["axial_strain"] == pytest.approx(
        expected_stress_pa / STEEL_MODULUS_PA, rel=1.0e-12
    )
    profile_pa = np.asarray(result["stress_profile_pa"])
    assert np.allclose(profile_pa, expected_stress_pa, rtol=1.0e-12)


def test_axial_bar_reports_mass_and_first_frequency() -> None:
    result = ca.solve_axial_bar(
        BAR_LENGTH_M,
        BAR_AREA_M2,
        STEEL_MODULUS_PA,
        BAR_LOAD_N,
        STEEL_DENSITY_KG_M3,
    )
    expected_mass_kg = STEEL_DENSITY_KG_M3 * BAR_AREA_M2 * BAR_LENGTH_M
    expected_frequency_hz = math.sqrt(STEEL_MODULUS_PA / STEEL_DENSITY_KG_M3) / (
        4.0 * BAR_LENGTH_M
    )
    assert result["mass_kg"] == pytest.approx(expected_mass_kg, rel=1.0e-12)
    assert result["first_natural_frequency_hz"] == pytest.approx(
        expected_frequency_hz, rel=1.0e-12
    )


def test_solve_axial_bar_rejects_a_bad_geometry() -> None:
    with pytest.raises(ValueError):
        ca.solve_axial_bar(-1.0, BAR_AREA_M2, STEEL_MODULUS_PA, BAR_LOAD_N, 7850.0)
    with pytest.raises(ValueError):
        ca.solve_axial_bar(1.0, BAR_AREA_M2, STEEL_MODULUS_PA, BAR_LOAD_N, 7850.0, 0)


def test_solve_euler_bernoulli_matches_the_closed_form() -> None:
    result = ca.solve_euler_bernoulli(
        BEAM_LENGTH_M,
        STEEL_MODULUS_PA,
        BEAM_INERTIA_M4,
        BEAM_TIP_LOAD_N,
        BEAM_OUTER_FIBER_M,
    )
    analytic_m = (
        BEAM_TIP_LOAD_N * BEAM_LENGTH_M**3 / (3.0 * STEEL_MODULUS_PA * BEAM_INERTIA_M4)
    )
    assert result["analytic_tip_deflection_m"] == pytest.approx(analytic_m, rel=1e-15)
    assert result["tip_deflection_m"] == pytest.approx(analytic_m, rel=1.0e-10)
    assert result["relative_error"] < 1.0e-10
    assert result["method"] == "euler_bernoulli_hermite_fem"
    assert result["validation"] == "cross-checked"


def test_beam_mesh_keeps_the_exact_tip_deflection() -> None:
    analytic_m = (
        BEAM_TIP_LOAD_N * BEAM_LENGTH_M**3 / (3.0 * STEEL_MODULUS_PA * BEAM_INERTIA_M4)
    )
    for elements in (1, 2, 4, 10, 20):
        result = ca.solve_euler_bernoulli(
            BEAM_LENGTH_M,
            STEEL_MODULUS_PA,
            BEAM_INERTIA_M4,
            BEAM_TIP_LOAD_N,
            BEAM_OUTER_FIBER_M,
            elements=elements,
        )
        assert result["tip_deflection_m"] == pytest.approx(analytic_m, rel=1.0e-10)


def test_beam_root_moment_and_stress_match_the_closed_form() -> None:
    result = ca.solve_euler_bernoulli(
        BEAM_LENGTH_M,
        STEEL_MODULUS_PA,
        BEAM_INERTIA_M4,
        BEAM_TIP_LOAD_N,
        BEAM_OUTER_FIBER_M,
    )
    expected_moment_n_m = BEAM_TIP_LOAD_N * BEAM_LENGTH_M
    expected_stress_pa = expected_moment_n_m * BEAM_OUTER_FIBER_M / BEAM_INERTIA_M4
    assert result["root_moment_n_m"] == pytest.approx(expected_moment_n_m, rel=1.0e-9)
    assert result["max_bending_stress_pa"] == pytest.approx(
        expected_stress_pa, rel=1.0e-9
    )
    assert result["analytic_max_bending_stress_pa"] == pytest.approx(
        expected_stress_pa, rel=1.0e-12
    )


def test_poiseuille_tube_matches_hagen_poiseuille() -> None:
    result = ca.solve_poiseuille(
        FLOW_LENGTH_M,
        WATER_VISCOSITY_PA_S,
        PRESSURE_DROP_PA,
        radius_m=TUBE_RADIUS_M,
        bins=HIGH_BINS,
    )
    gradient_pa_per_m = PRESSURE_DROP_PA / FLOW_LENGTH_M
    analytic_q_m3_s = (
        math.pi * TUBE_RADIUS_M**4 * gradient_pa_per_m / (8.0 * WATER_VISCOSITY_PA_S)
    )
    assert result["geometry"] == "tube"
    assert result["analytic_flow_rate_m3_s"] == pytest.approx(
        analytic_q_m3_s, rel=1.0e-12
    )
    assert result["flow_rate_m3_s"] == pytest.approx(analytic_q_m3_s, rel=1.0e-10)
    assert result["relative_error"] < 1.0e-10
    assert result["max_velocity_m_s"] == pytest.approx(
        2.0 * result["mean_velocity_m_s"], rel=1.0e-9
    )
    assert result["validation"] == "cross-checked"


def test_poiseuille_channel_matches_plane_poiseuille() -> None:
    result = ca.solve_poiseuille(
        FLOW_LENGTH_M,
        WATER_VISCOSITY_PA_S,
        PRESSURE_DROP_PA,
        half_height_m=CHANNEL_HALF_HEIGHT_M,
        width_m=CHANNEL_WIDTH_M,
        bins=HIGH_BINS,
    )
    gradient_pa_per_m = PRESSURE_DROP_PA / FLOW_LENGTH_M
    analytic_q_m3_s = (
        2.0
        / 3.0
        * CHANNEL_WIDTH_M
        * CHANNEL_HALF_HEIGHT_M**3
        * gradient_pa_per_m
        / WATER_VISCOSITY_PA_S
    )
    assert result["geometry"] == "channel"
    assert result["analytic_flow_rate_m3_s"] == pytest.approx(
        analytic_q_m3_s, rel=1.0e-12
    )
    assert result["flow_rate_m3_s"] == pytest.approx(analytic_q_m3_s, rel=1.0e-10)
    assert result["relative_error"] < 1.0e-10
    assert result["max_velocity_m_s"] == pytest.approx(
        1.5 * result["mean_velocity_m_s"], rel=1.0e-9
    )


def test_poiseuille_mesh_converges_second_order() -> None:
    coarse = ca.solve_poiseuille(
        FLOW_LENGTH_M,
        WATER_VISCOSITY_PA_S,
        PRESSURE_DROP_PA,
        radius_m=TUBE_RADIUS_M,
        bins=64,
    )
    fine = ca.solve_poiseuille(
        FLOW_LENGTH_M,
        WATER_VISCOSITY_PA_S,
        PRESSURE_DROP_PA,
        radius_m=TUBE_RADIUS_M,
        bins=128,
    )
    assert fine["relative_error"] < coarse["relative_error"]
    ratio = coarse["relative_error"] / fine["relative_error"]
    assert 3.5 < ratio < 4.5


def test_poiseuille_hydraulic_resistance() -> None:
    result = ca.solve_poiseuille(
        FLOW_LENGTH_M,
        WATER_VISCOSITY_PA_S,
        PRESSURE_DROP_PA,
        radius_m=TUBE_RADIUS_M,
    )
    analytic_resistance_pa_s_m3 = (
        8.0 * WATER_VISCOSITY_PA_S * FLOW_LENGTH_M / (math.pi * TUBE_RADIUS_M**4)
    )
    assert result["hydraulic_resistance_pa_s_m3"] == pytest.approx(
        analytic_resistance_pa_s_m3, rel=1.0e-4
    )


def test_poiseuille_needs_a_geometry() -> None:
    with pytest.raises(ValueError):
        ca.solve_poiseuille(FLOW_LENGTH_M, WATER_VISCOSITY_PA_S, PRESSURE_DROP_PA)


def cubic_carbon_mapping() -> dict:
    positions_m: list[float] = []
    for ix in range(2):
        for iy in range(2):
            for iz in range(2):
                positions_m += [ix * 1.0e-10, iy * 1.0e-10, iz * 1.0e-10]
    return {"elements": ["C"] * 8, "positions_m": positions_m}


def test_voxelize_plain_atoms_reports_density_and_box() -> None:
    result = ca.voxelize_part(cubic_carbon_mapping())
    carbon_mass_kg = 12.011 * ca.AMU_TO_KG
    assert result["atom_count"] == 8
    assert result["mass_kg"] == pytest.approx(8.0 * carbon_mass_kg, rel=1.0e-12)
    assert result["bounding_box_size_m"] == pytest.approx([1.0e-10] * 3, rel=1e-12)
    assert result["bounding_box_volume_m3"] == pytest.approx(1.0e-30, rel=1.0e-12)
    assert result["center_of_mass_m"] == pytest.approx([5.0e-11] * 3, rel=1.0e-12)
    assert result["voxel_size_m"] == pytest.approx(5.0e-11, rel=1.0e-12)
    assert result["occupied_voxel_count"] == 8
    expected_density_kg_m3 = result["mass_kg"] / result["voxel_volume_m3"] / 8.0
    assert result["effective_density_kg_m3"] == pytest.approx(
        expected_density_kg_m3, rel=1.0e-12
    )
    assert result["method"] == "atomic_voxelization"
    assert result["validation"] == "unverified"


def test_voxelize_accepts_atom_dictionaries() -> None:
    mapping = {
        "atoms": [
            {"element": "C", "position_m": [0.0, 0.0, 0.0]},
            {"element": "H", "position_m": [1.0e-10, 0.0, 0.0]},
        ]
    }
    result = ca.voxelize_part(mapping, voxel_size_m=1.0e-10)
    assert result["atom_count"] == 2
    expected_mass_kg = (12.011 + 1.008) * ca.AMU_TO_KG
    assert result["mass_kg"] == pytest.approx(expected_mass_kg, rel=1.0e-12)


def test_voxelize_core_part_when_available() -> None:
    nanocad = pytest.importorskip("nanocad")
    if getattr(nanocad, "_core", None) is None:
        pytest.skip("nanocad._core is not built")
    topology = nanocad.Topology()
    topology.add_atom(6, [0.0, 0.0, 0.0], 0.0, "C")
    topology.add_atom(6, [1.0e-10, 0.0, 0.0], 0.0, "C")
    part = nanocad.Part("test", topology)
    result = ca.voxelize_part(part, voxel_size_m=1.0e-10)
    assert result["atom_count"] == 2
    assert result["bounding_box_size_m"][0] == pytest.approx(1.0e-10, rel=1.0e-12)


def test_voxelize_rejects_bad_input() -> None:
    with pytest.raises(ValueError):
        ca.voxelize_part({"elements": [], "positions_m": []})
    with pytest.raises(ValueError):
        ca.voxelize_atom_list(["Xx"], [0.0, 0.0, 0.0])
    with pytest.raises(TypeError):
        ca.voxelize_part(object())


PLATE_LENGTH_M = 1.0
PLATE_HEIGHT_M = 0.05
PLATE_THICKNESS_M = 1.0
PLATE_MODULUS_PA = 200.0e9
PLATE_POISSON = 0.3
PLATE_TIP_LOAD_N = 100.0

CHANNEL_2D_HEIGHT_M = 1.0e-3
CHANNEL_2D_LENGTH_M = 1.0
CHANNEL_2D_DROP_PA = 1.0


def test_solve_plane_stress_cantilever_matches_a_refined_reference() -> None:
    pytest.importorskip("scipy")
    result = ca.solve_plane_stress_cantilever(
        PLATE_LENGTH_M,
        PLATE_HEIGHT_M,
        PLATE_THICKNESS_M,
        PLATE_MODULUS_PA,
        PLATE_POISSON,
        PLATE_TIP_LOAD_N,
        elements_x=160,
    )
    assert result["method"] == "plane_stress_q4_fem_scipy"
    assert result["validation"] == "cross-checked"
    assert result["converged"] is True
    assert result["relative_error"] < 1.0e-2
    assert result["reference_tip_deflection_m"] > result["tip_deflection_m"]


def test_solve_plane_stress_cantilever_approaches_euler_bernoulli() -> None:
    pytest.importorskip("scipy")
    result = ca.solve_plane_stress_cantilever(
        PLATE_LENGTH_M,
        PLATE_HEIGHT_M,
        PLATE_THICKNESS_M,
        PLATE_MODULUS_PA,
        PLATE_POISSON,
        PLATE_TIP_LOAD_N,
        elements_x=160,
    )
    assert result["relative_error_vs_euler_bernoulli"] < 1.0e-2
    assert result["tip_deflection_m"] < result["analytic_tip_deflection_m"]


def test_solve_plane_stress_cantilever_converges_with_the_mesh() -> None:
    pytest.importorskip("scipy")
    coarse = ca.solve_plane_stress_cantilever(
        PLATE_LENGTH_M,
        PLATE_HEIGHT_M,
        PLATE_THICKNESS_M,
        PLATE_MODULUS_PA,
        PLATE_POISSON,
        PLATE_TIP_LOAD_N,
        elements_x=80,
        refinement=1,
    )
    fine = ca.solve_plane_stress_cantilever(
        PLATE_LENGTH_M,
        PLATE_HEIGHT_M,
        PLATE_THICKNESS_M,
        PLATE_MODULUS_PA,
        PLATE_POISSON,
        PLATE_TIP_LOAD_N,
        elements_x=160,
        refinement=1,
    )
    assert (
        fine["relative_error_vs_euler_bernoulli"]
        < coarse["relative_error_vs_euler_bernoulli"]
    )


def test_solve_plane_stress_cantilever_stress_matches_the_beam() -> None:
    pytest.importorskip("scipy")
    result = ca.solve_plane_stress_cantilever(
        PLATE_LENGTH_M,
        PLATE_HEIGHT_M,
        PLATE_THICKNESS_M,
        PLATE_MODULUS_PA,
        PLATE_POISSON,
        PLATE_TIP_LOAD_N,
        elements_x=320,
    )
    relative = (
        abs(result["max_axial_stress_pa"] - result["analytic_max_axial_stress_pa"])
        / result["analytic_max_axial_stress_pa"]
    )
    assert relative < 2.0e-2


def test_solve_plane_stress_cantilever_rejects_bad_input() -> None:
    with pytest.raises(ValueError):
        ca.solve_plane_stress_cantilever(
            -1.0, PLATE_HEIGHT_M, PLATE_THICKNESS_M, PLATE_MODULUS_PA, 0.3, 1.0
        )
    with pytest.raises(ValueError):
        ca.solve_plane_stress_cantilever(
            PLATE_LENGTH_M,
            PLATE_HEIGHT_M,
            PLATE_THICKNESS_M,
            PLATE_MODULUS_PA,
            0.6,
            1.0,
        )


def test_solve_poiseuille_2d_matches_the_parabola() -> None:
    pytest.importorskip("scipy")
    result = ca.solve_poiseuille_2d(
        CHANNEL_2D_LENGTH_M,
        CHANNEL_2D_HEIGHT_M,
        WATER_VISCOSITY_PA_S,
        CHANNEL_2D_DROP_PA,
        elements_x=16,
        elements_y=64,
    )
    assert result["method"] == "poiseuille_2d_stokes_finite_difference_scipy"
    assert result["validation"] == "cross-checked"
    assert result["converged"] is True
    assert result["relative_error"] < 1.0e-3
    assert result["max_velocity_m_s"] == pytest.approx(
        result["analytic_max_velocity_m_s"], rel=1.0e-3
    )
    assert result["mean_velocity_m_s"] == pytest.approx(
        result["analytic_mean_velocity_m_s"], rel=1.0e-3
    )
    profile_m_s = np.asarray(result["velocity_profile_m_s"])
    assert profile_m_s[0] == pytest.approx(profile_m_s[-1], rel=1.0e-12)


def test_solve_poiseuille_2d_converges_second_order() -> None:
    pytest.importorskip("scipy")
    coarse = ca.solve_poiseuille_2d(
        CHANNEL_2D_LENGTH_M,
        CHANNEL_2D_HEIGHT_M,
        WATER_VISCOSITY_PA_S,
        CHANNEL_2D_DROP_PA,
        elements_x=16,
        elements_y=32,
    )
    fine = ca.solve_poiseuille_2d(
        CHANNEL_2D_LENGTH_M,
        CHANNEL_2D_HEIGHT_M,
        WATER_VISCOSITY_PA_S,
        CHANNEL_2D_DROP_PA,
        elements_x=16,
        elements_y=64,
    )
    ratio = coarse["relative_error"] / fine["relative_error"]
    assert 3.5 < ratio < 4.5


def test_solve_poiseuille_2d_needs_a_mesh() -> None:
    with pytest.raises(ValueError):
        ca.solve_poiseuille_2d(1.0, 1.0e-3, 1.0e-3, 1.0, elements_y=1)


def test_continuum_falls_back_without_scipy(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(ca, "_SCIPY_AVAILABLE", False)
    plate = ca.solve_plane_stress_cantilever(
        PLATE_LENGTH_M,
        PLATE_HEIGHT_M,
        PLATE_THICKNESS_M,
        PLATE_MODULUS_PA,
        PLATE_POISSON,
        PLATE_TIP_LOAD_N,
    )
    assert plate["validation"] == "unverified"
    assert plate["method"] == "plane_stress_cantilever_closed_form_fallback"
    assert plate["tip_deflection_m"] == pytest.approx(
        plate["analytic_tip_deflection_m"], rel=1.0e-12
    )
    flow = ca.solve_poiseuille_2d(
        CHANNEL_2D_LENGTH_M,
        CHANNEL_2D_HEIGHT_M,
        WATER_VISCOSITY_PA_S,
        CHANNEL_2D_DROP_PA,
    )
    assert flow["validation"] == "unverified"
    assert flow["method"] == "poiseuille_2d_closed_form_fallback"
    assert flow["velocity_profile_m_s"] == pytest.approx(
        flow["analytic_velocity_profile_m_s"], rel=1.0e-12
    )
