"""Tests for the PME electrostatics adapter.

Every quantity is SI.  Positions and box lengths are in metres [m].  Energy is
in joules [J].  Charge is in units of the elementary charge.  Variable names
carry a unit suffix.

The tests need OpenMM.  They skip cleanly when OpenMM is absent.  The
comparison against the engine cutoff model also needs the compiled
``nanocad._core`` extension and skips when that is absent.
"""

from __future__ import annotations

import numpy as np
import pytest

openmm = pytest.importorskip("openmm")

import nanocad  # noqa: E402
from nanocad import pme_adapter as pa  # noqa: E402

NANOMETRE_M = 1.0e-9
ELEMENTARY_CHARGE_C = 1.602176634e-19
COULOMB_CONSTANT_N_M2_PER_C2 = 8.9875517923e9


def _has_core() -> bool:
    return getattr(nanocad, "_core", None) is not None


def test_openmm_backend_is_available() -> None:
    assert pa.available() is True


def test_pme_matches_direct_coulomb_in_a_large_box() -> None:
    separation_m = 0.5 * NANOMETRE_M
    positions_m = [
        4.0 * NANOMETRE_M,
        5.0 * NANOMETRE_M,
        5.0 * NANOMETRE_M,
        4.5 * NANOMETRE_M,
        5.0 * NANOMETRE_M,
        5.0 * NANOMETRE_M,
    ]
    result = pa.electrostatic_energy_pme(
        ["Na", "Cl"],
        positions_m,
        [1.0, -1.0],
        10.0 * NANOMETRE_M,
        2.0 * NANOMETRE_M,
    )
    direct_j = (
        COULOMB_CONSTANT_N_M2_PER_C2
        * (1.0 * -1.0)
        * ELEMENTARY_CHARGE_C**2
        / separation_m
    )
    relative = abs(result["energy_j"] - direct_j) / abs(direct_j)
    assert relative < 1.0e-2
    assert result["net_charge_e"] == pytest.approx(0.0)
    assert result["method"] == "openmm_pme"
    assert result["platform"] == "Reference"
    assert result["validation"] == "unverified"


def test_pme_is_deterministic_on_the_reference_platform() -> None:
    positions_m = [
        0.5 * NANOMETRE_M,
        0.5 * NANOMETRE_M,
        0.5 * NANOMETRE_M,
        1.2 * NANOMETRE_M,
        0.6 * NANOMETRE_M,
        0.5 * NANOMETRE_M,
    ]
    first = pa.electrostatic_energy_pme(
        ["Na", "Cl"], positions_m, [1.0, -1.0], 3.0 * NANOMETRE_M, 1.2 * NANOMETRE_M
    )
    second = pa.electrostatic_energy_pme(
        ["Na", "Cl"], positions_m, [1.0, -1.0], 3.0 * NANOMETRE_M, 1.2 * NANOMETRE_M
    )
    assert first["energy_j"] == second["energy_j"]
    assert first["pme_grid"] == second["pme_grid"]


def test_pme_forces_match_a_finite_difference() -> None:
    positions_m = np.array(
        [
            0.5 * NANOMETRE_M,
            0.5 * NANOMETRE_M,
            0.5 * NANOMETRE_M,
            1.2 * NANOMETRE_M,
            0.6 * NANOMETRE_M,
            0.5 * NANOMETRE_M,
        ]
    )
    result = pa.electrostatic_energy_pme(
        ["Na", "Cl"],
        positions_m,
        [1.0, -1.0],
        3.0 * NANOMETRE_M,
        1.2 * NANOMETRE_M,
        with_forces=True,
    )
    step_m = 1.0e-12
    plus = positions_m.copy()
    plus[0] += step_m
    minus = positions_m.copy()
    minus[0] -= step_m
    energy_plus = pa.electrostatic_energy_pme(
        ["Na", "Cl"], plus, [1.0, -1.0], 3.0 * NANOMETRE_M, 1.2 * NANOMETRE_M
    )["energy_j"]
    energy_minus = pa.electrostatic_energy_pme(
        ["Na", "Cl"], minus, [1.0, -1.0], 3.0 * NANOMETRE_M, 1.2 * NANOMETRE_M
    )["energy_j"]
    numerical_force_n = -(energy_plus - energy_minus) / (2.0 * step_m)
    assert result["forces_n"][0] == pytest.approx(numerical_force_n, rel=1.0e-3)


def test_cutoff_and_pme_report_the_long_range_difference() -> None:
    if not _has_core():
        pytest.skip("nanocad._core is not built")
    positions_m = [
        0.5 * NANOMETRE_M,
        0.5 * NANOMETRE_M,
        0.5 * NANOMETRE_M,
        1.2 * NANOMETRE_M,
        0.6 * NANOMETRE_M,
        0.5 * NANOMETRE_M,
    ]
    result = pa.compare_cutoff_and_pme(
        ["Na", "Cl"],
        positions_m,
        [1.0, -1.0],
        3.0 * NANOMETRE_M,
        1.2 * NANOMETRE_M,
    )
    assert result["validation"] == "cross-checked"
    assert result["cutoff"]["method"] == "engine_cutoff_electrostatic"
    assert result["pme"]["method"] == "openmm_pme"
    assert np.sign(result["difference_j"]) == np.sign(result["cutoff_energy_j"])
    assert result["difference_j"] != 0.0
    assert abs(result["relative_difference"]) < 0.1
    assert "long-range" in result["difference_explanation"]


def test_cutoff_and_pme_agree_as_the_box_grows() -> None:
    if not _has_core():
        pytest.skip("nanocad._core is not built")
    small = pa.compare_cutoff_and_pme(
        ["Na", "Cl"],
        [
            0.5 * NANOMETRE_M,
            0.5 * NANOMETRE_M,
            0.5 * NANOMETRE_M,
            1.2 * NANOMETRE_M,
            0.6 * NANOMETRE_M,
            0.5 * NANOMETRE_M,
        ],
        [1.0, -1.0],
        3.0 * NANOMETRE_M,
        1.2 * NANOMETRE_M,
    )
    large = pa.compare_cutoff_and_pme(
        ["Na", "Cl"],
        [
            4.0 * NANOMETRE_M,
            5.0 * NANOMETRE_M,
            5.0 * NANOMETRE_M,
            4.5 * NANOMETRE_M,
            5.0 * NANOMETRE_M,
            5.0 * NANOMETRE_M,
        ],
        [1.0, -1.0],
        10.0 * NANOMETRE_M,
        2.0 * NANOMETRE_M,
    )
    assert abs(large["relative_difference"]) < abs(small["relative_difference"])


def test_cutoff_energy_matches_the_direct_coulomb_pair() -> None:
    if not _has_core():
        pytest.skip("nanocad._core is not built")
    separation_m = 0.7071067811865476 * NANOMETRE_M
    result = pa.cutoff_electrostatic_energy_j(
        ["Na", "Cl"],
        [
            0.5 * NANOMETRE_M,
            0.5 * NANOMETRE_M,
            0.5 * NANOMETRE_M,
            1.2 * NANOMETRE_M,
            0.6 * NANOMETRE_M,
            0.5 * NANOMETRE_M,
        ],
        [1.0, -1.0],
        3.0 * NANOMETRE_M,
        1.2 * NANOMETRE_M,
    )
    direct_j = (
        COULOMB_CONSTANT_N_M2_PER_C2
        * (1.0 * -1.0)
        * ELEMENTARY_CHARGE_C**2
        / separation_m
    )
    assert result["energy_j"] == pytest.approx(direct_j, rel=1.0e-3)


def test_pme_rejects_a_cutoff_over_half_the_box() -> None:
    with pytest.raises(ValueError):
        pa.electrostatic_energy_pme(
            ["Na", "Cl"],
            [0.0, 0.0, 0.0, 1.0 * NANOMETRE_M, 0.0, 0.0],
            [1.0, -1.0],
            2.0 * NANOMETRE_M,
            1.5 * NANOMETRE_M,
        )


def test_pme_rejects_bad_input() -> None:
    with pytest.raises(ValueError):
        pa.electrostatic_energy_pme(
            ["Na"], [0.0, 0.0, 0.0], [1.0, -1.0], 1.0 * NANOMETRE_M, 0.3 * NANOMETRE_M
        )


def test_pme_backend_error_when_openmm_is_hidden(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(pa, "_OPENMM_AVAILABLE", False)
    with pytest.raises(pa.PmeBackendError):
        pa.electrostatic_energy_pme(
            ["Na", "Cl"],
            [0.0, 0.0, 0.0, 1.0 * NANOMETRE_M, 0.0, 0.0],
            [1.0, -1.0],
            3.0 * NANOMETRE_M,
            1.0 * NANOMETRE_M,
        )
