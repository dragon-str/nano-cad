"""Tests for the RDKit adapter.

The suite must pass with no RDKit installed. Tests that call RDKit skip
through ``pytest.importorskip("rdkit")``. One test proves the typed error for
the absent case, and it runs everywhere.
"""

from __future__ import annotations

import sys

import pytest

from nanocad import rdkit_adapter as rd

WATER_M = rd.Molecule(
    elements=("O", "H", "H"),
    positions_m=(
        (0.0, 0.0, 0.0),
        (0.9572e-10, 0.0, 0.0),
        (-0.239987e-10, 0.927e-10, 0.0),
    ),
    bonds=(rd.Bond(0, 1, 1), rd.Bond(0, 2, 1)),
    name="water",
)


def test_module_imports_without_rdkit() -> None:
    assert isinstance(rd.rdkit_available(), bool)


def test_missing_rdkit_raises_typed_error(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "rdkit", None)
    assert rd.rdkit_available() is False
    with pytest.raises(rd.RdkitUnavailableError):
        rd.write_smiles(WATER_M)


def test_a_malformed_molecule_raises_typed_error() -> None:
    bad = rd.Molecule(elements=("C", "C"), bonds=(rd.Bond(0, 5, 1),))
    with pytest.raises(rd.MoleculeConversionError):
        rd.to_rdkit_mol(bad)
    with pytest.raises(rd.MoleculeConversionError):
        rd.to_rdkit_mol(rd.Molecule(elements=()))


def test_rdkit_mol_round_trip_keeps_the_graph() -> None:
    pytest.importorskip("rdkit")
    mol = rd.to_rdkit_mol(WATER_M)
    assert mol.GetNumAtoms() == 3
    assert mol.GetNumBonds() == 2
    back = rd.from_rdkit_mol(mol, name="water")
    assert back.elements == WATER_M.elements
    assert len(back.bonds) == len(WATER_M.bonds)
    assert back.name == "water"


def test_coordinates_round_trip_through_angstrom() -> None:
    pytest.importorskip("rdkit")
    mol = rd.to_rdkit_mol(WATER_M)
    back = rd.from_rdkit_mol(mol)
    for original, restored in zip(WATER_M.positions_m, back.positions_m, strict=True):
        for left, right in zip(original, restored, strict=True):
            assert left == pytest.approx(right, abs=1.0e-14)


def test_smiles_round_trip_is_canonical() -> None:
    pytest.importorskip("rdkit")
    molecule = rd.read_smiles("CCO")
    assert len(molecule.elements) == 9
    assert rd.write_smiles(molecule) == "CCO"


def test_read_smiles_rejects_bad_text() -> None:
    pytest.importorskip("rdkit")
    with pytest.raises(rd.StructureParseError):
        rd.read_smiles("C1CC")
    with pytest.raises(rd.StructureParseError):
        rd.read_smiles("")


def test_sdf_round_trip_keeps_the_topology() -> None:
    pytest.importorskip("rdkit")
    molecule = rd.read_smiles("CCO")
    text = rd.write_sdf(molecule)
    assert text.rstrip().endswith("$$$$")
    parsed = rd.read_sdf(text)
    assert len(parsed.elements) == len(molecule.elements)
    assert len(parsed.bonds) == len(molecule.bonds)


def test_smiles_positions_are_in_metres() -> None:
    pytest.importorskip("rdkit")
    molecule = rd.read_smiles("CCO")
    magnitude = max(
        abs(value) for position in molecule.positions_m for value in position
    )
    assert 0.0 < magnitude < 1.0e-9


def test_molecule_from_core_topology() -> None:
    core = pytest.importorskip("nanocad._core")
    pytest.importorskip("rdkit")
    topology = core.Topology()
    topology.add_atom(8, [0.0, 0.0, 0.0], 0.0, "O")
    topology.add_atom(1, [0.9572e-10, 0.0, 0.0], 0.0, "H")
    topology.add_bond(0, 1, 1, "single")
    molecule = rd.molecule_from_topology(topology)
    assert molecule.elements == ("O", "H")
    assert len(molecule.bonds) == 1
    assert molecule.bonds[0].order == 1
