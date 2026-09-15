"""Tests for the quantum adapter.

The suite must pass with no backend installed. Backend-specific tests skip
through ``pytest.importorskip``. Tests that need only the compiled core skip
when ``nanocad._core`` is absent.
"""

from __future__ import annotations

import math

import pytest

from nanocad import quantum_adapter as qa

# H2 at 0.735 Angstrom, the STO-3G reference geometry.
H2_POSITIONS_M = [[0.0, 0.0, 0.0], [0.0, 0.0, 7.35e-11]]


def test_module_imports_without_a_backend() -> None:
    assert isinstance(qa.available_backends(), list)


def test_no_backend_raises_typed_error(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(qa, "available_backends", lambda: [])
    with pytest.raises(qa.NoQuantumBackendError):
        qa.compute_energy(["H", "H"], H2_POSITIONS_M)


def test_unknown_backend_raises_value_error() -> None:
    with pytest.raises(ValueError):
        qa.compute_energy(["H", "H"], H2_POSITIONS_M, backend="not-a-backend")


def test_bad_input_raises_value_error() -> None:
    with pytest.raises(ValueError):
        qa.compute_energy(["H", "H"], [[0.0, 0.0, 0.0]])
    with pytest.raises(ValueError):
        qa.compute_energy(["Xx"], [[0.0, 0.0, 0.0]])
    with pytest.raises(ValueError):
        qa.compute_energy(["H", "H"], H2_POSITIONS_M, multiplicity=4)
    with pytest.raises(ValueError):
        qa.compute_energy(["H", "H"], [[0.0, 0.0, float("nan")]] * 2)


@pytest.mark.parametrize("backend", ["pyscf", "xtb"])
def test_h2_energy_with_a_quantum_backend(backend: str) -> None:
    pytest.importorskip(backend)
    result = qa.compute_energy(["H", "H"], H2_POSITIONS_M, backend=backend)
    assert result["backend"] == backend
    assert result["is_quantum"] is True
    assert result["method"].startswith(backend)
    assert result["atom_count"] == 2
    assert math.isfinite(result["energy_j"])
    # A bound H2 has a negative energy of order a few times 10^-18 J.
    assert -1.0e-17 < result["energy_j"] < -1.0e-19


def test_h2_pyscf_matches_a_cited_reference() -> None:
    pytest.importorskip("pyscf")
    result = qa.compute_energy(["H", "H"], H2_POSITIONS_M, backend="pyscf")
    assert result["validation"].startswith("reference:")
    assert "Hehre" in result["validation"]
    assert result["level_of_theory"].startswith("RHF")


def test_ase_fallback_is_labelled_not_quantum() -> None:
    pytest.importorskip("ase")
    result = qa.compute_energy(["H", "H"], H2_POSITIONS_M, backend="ase")
    assert result["backend"] == "ase"
    assert result["is_quantum"] is False
    assert "not quantum" in result["method"].lower()
    assert result["validation"] == "unverified"
    assert math.isfinite(result["energy_j"])


def test_atomic_numbers_and_flat_positions_are_accepted() -> None:
    pytest.importorskip("ase")
    by_symbol = qa.compute_energy(["H", "H"], H2_POSITIONS_M, backend="ase")
    by_number = qa.compute_energy([1, 1], H2_POSITIONS_M, backend="ase")
    by_flat = qa.compute_energy(
        ["H", "H"], [0.0, 0.0, 0.0, 0.0, 0.0, 7.35e-11], backend="ase"
    )
    assert by_number["energy_j"] == pytest.approx(by_symbol["energy_j"])
    assert by_flat["energy_j"] == pytest.approx(by_symbol["energy_j"])


def test_topology_to_geometry_from_the_core() -> None:
    core = pytest.importorskip("nanocad._core")
    topology = core.Topology()
    topology.add_atom(1, [0.0, 0.0, 0.0], 0.0, "H")
    topology.add_atom(1, [0.0, 0.0, 7.35e-11], 0.0, "H")
    symbols, positions, charge = qa.topology_to_geometry(topology)
    assert symbols == ["H", "H"]
    assert positions[1] == pytest.approx([0.0, 0.0, 7.35e-11])
    assert charge == 0


def test_topology_to_geometry_unwraps_a_core_part() -> None:
    core = pytest.importorskip("nanocad._core")
    topology = core.Topology()
    topology.add_atom(1, [0.0, 0.0, 0.0], 0.0, "H")
    topology.add_atom(1, [0.0, 0.0, 7.35e-11], 0.0, "H")
    part = core.Part("h2", topology)
    symbols, positions, _charge = qa.topology_to_geometry(part)
    assert symbols == ["H", "H"]
    assert len(positions) == 2


def test_compute_energy_from_a_core_topology() -> None:
    core = pytest.importorskip("nanocad._core")
    if not qa.available_backends():
        pytest.skip("no backend installed")
    topology = core.Topology()
    topology.add_atom(1, [0.0, 0.0, 0.0], 0.0, "H")
    topology.add_atom(1, [0.0, 0.0, 7.35e-11], 0.0, "H")
    result = qa.compute_energy_from_topology(topology)
    assert math.isfinite(result["energy_j"])
    assert result["method"]
