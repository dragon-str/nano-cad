"""RDKit interop for small molecules.

The module gives NanoCAD one small, honest bridge to RDKit. It imports RDKit
lazily inside each function, so the package imports with no RDKit installed.
When RDKit is absent, every RDKit-touching call raises
:class:`RdkitUnavailableError`.

The bridge carries a :class:`Molecule`, a plain value: element symbols,
positions in SI metres, and bonds with integer order.

    to_rdkit_mol(molecule) -> rdkit.Chem.Mol
    from_rdkit_mol(mol) -> Molecule
    read_smiles(text) -> Molecule
    write_smiles(molecule) -> str
    read_sdf(text) -> Molecule
    write_sdf(molecule) -> str

RDKit works in Angstrom. The model works in metres. All conversions happen at
the boundary of this module.

Honest limits
-------------
* A SMILES string carries no coordinates. :func:`read_smiles` makes 2D
  depiction coordinates, which are not physical. It says so in the result.
* Aromatic bond order 1.5 becomes order 1. Aromaticity is not carried.
* The adapter reads the first molecule of an SDF file. It does not read a
  multi-molecule stream.
"""

from __future__ import annotations

import importlib.util
import math
from dataclasses import dataclass
from io import BytesIO
from typing import Any

__all__ = [
    "ANGSTROM_TO_METRE",
    "Bond",
    "Molecule",
    "MoleculeConversionError",
    "RdkitAdapterError",
    "RdkitUnavailableError",
    "StructureParseError",
    "from_rdkit_mol",
    "molecule_from_topology",
    "rdkit_available",
    "read_sdf",
    "read_smiles",
    "to_rdkit_mol",
    "write_sdf",
    "write_smiles",
]

ANGSTROM_TO_METRE = 1.0e-10
METRE_TO_ANGSTROM = 1.0e10


class RdkitAdapterError(RuntimeError):
    """Base class for every error raised by this module."""


class RdkitUnavailableError(RdkitAdapterError):
    """RDKit is not importable."""


class MoleculeConversionError(RdkitAdapterError):
    """A :class:`Molecule` value is malformed and cannot become an RDKit mol."""


class StructureParseError(RdkitAdapterError):
    """RDKit could not parse a SMILES or SDF string."""


@dataclass(frozen=True)
class Bond:
    """A bond between two atom indices with an integer order."""

    u: int
    v: int
    order: int = 1


@dataclass(frozen=True)
class Molecule:
    """A small molecule: element symbols, positions in metres, and bonds.

    ``positions_m`` may be empty when no coordinates exist. Otherwise its
    length equals the length of ``elements``.
    """

    elements: tuple[str, ...]
    positions_m: tuple[tuple[float, float, float], ...] = ()
    bonds: tuple[Bond, ...] = ()
    name: str = ""


def rdkit_available() -> bool:
    """Return True when RDKit imports. This is cheap and side-effect free."""
    try:
        if importlib.util.find_spec("rdkit") is None:
            return False
    except (ImportError, ValueError):
        return False
    try:
        _require_rdkit()
    except RdkitUnavailableError:
        return False
    return True


def _require_rdkit() -> tuple[Any, Any]:
    """Import RDKit lazily and return ``(Chem, Point3D)``.

    Raise :class:`RdkitUnavailableError` when RDKit is absent.
    """
    try:
        from rdkit import Chem
        from rdkit.Geometry import Point3D
    except Exception as exc:  # noqa: BLE001 - any import failure means absent
        raise RdkitUnavailableError(
            "RDKit is not importable; install it with 'pip install rdkit'"
        ) from exc
    return Chem, Point3D


def _bond_type(chem: Any, order: int) -> Any:
    table = {
        1: chem.BondType.SINGLE,
        2: chem.BondType.DOUBLE,
        3: chem.BondType.TRIPLE,
        4: getattr(chem.BondType, "QUADRUPLE", chem.BondType.SINGLE),
    }
    try:
        return table[order]
    except KeyError as exc:
        raise MoleculeConversionError(f"unsupported bond order {order}") from exc


def _as_molecule(value: Any) -> Molecule:
    """Accept a :class:`Molecule` or a duck-typed core topology."""
    if isinstance(value, Molecule):
        return value
    if hasattr(value, "elements") and hasattr(value, "positions_m"):
        return molecule_from_topology(value)
    raise MoleculeConversionError(
        "expected a Molecule or a topology with elements() and positions_m()"
    )


def molecule_from_topology(topology: Any) -> Molecule:
    """Build a :class:`Molecule` from a ``nanocad._core`` topology or part.

    The input must expose ``elements()`` (atomic numbers) and
    ``positions_m()`` (a flat list in metres). When it also exposes
    ``bond_count()``, ``bond_u()``, ``bond_v()`` and ``bond_order()``, the
    bonds are carried across too. RDKit provides the number-to-symbol table, so
    this function needs RDKit.
    """
    chem, _ = _require_rdkit()
    core = getattr(topology, "topology", topology)
    raw_elements = list(core.elements())
    table = chem.GetPeriodicTable()
    elements = tuple(table.GetElementSymbol(int(number)) for number in raw_elements)
    flat = [float(value) for value in core.positions_m()]
    positions: tuple[tuple[float, float, float], ...] = ()
    if flat:
        if len(flat) != 3 * len(elements):
            raise MoleculeConversionError(
                "the position buffer length is not three times the atom count"
            )
        positions = tuple(
            (flat[3 * i], flat[3 * i + 1], flat[3 * i + 2])
            for i in range(len(elements))
        )
    bonds: tuple[Bond, ...] = ()
    raw_count = getattr(core, "bond_count", 0)
    bond_count = int(raw_count() if callable(raw_count) else raw_count)
    if bond_count > 0:
        collected: list[Bond] = []
        for index in range(bond_count):
            raw = core.bond(index)
            if raw is None:
                continue
            u, v, order = raw[0], raw[1], raw[2]
            collected.append(Bond(int(u), int(v), int(order)))
        bonds = tuple(collected)
    return Molecule(elements=elements, positions_m=positions, bonds=bonds)


def _validate(molecule: Molecule) -> None:
    if not molecule.elements:
        raise MoleculeConversionError("the molecule has no atoms")
    if molecule.positions_m and len(molecule.positions_m) != len(molecule.elements):
        raise MoleculeConversionError(
            "the position count does not equal the atom count"
        )
    for bond in molecule.bonds:
        if not 0 <= bond.u < len(molecule.elements):
            raise MoleculeConversionError(f"bond endpoint {bond.u} is out of range")
        if not 0 <= bond.v < len(molecule.elements):
            raise MoleculeConversionError(f"bond endpoint {bond.v} is out of range")
    for position in molecule.positions_m:
        if not all(math.isfinite(value) for value in position):
            raise MoleculeConversionError("a position is not finite")


def to_rdkit_mol(molecule: Any) -> Any:
    """Convert a :class:`Molecule` to an RDKit mol with one conformer.

    The conformer holds positions in Angstrom when the molecule has positions.
    Raise :class:`RdkitUnavailableError` without RDKit and
    :class:`MoleculeConversionError` for a malformed molecule.
    """
    value = _as_molecule(molecule)
    _validate(value)
    chem, point3d = _require_rdkit()

    editable = chem.RWMol()
    for symbol in value.elements:
        try:
            editable.AddAtom(chem.Atom(symbol))
        except Exception as exc:  # noqa: BLE001 - RDKit rejects a bad symbol
            raise MoleculeConversionError(f"unknown element {symbol!r}") from exc
    for bond in value.bonds:
        try:
            editable.AddBond(bond.u, bond.v, _bond_type(chem, bond.order))
        except MoleculeConversionError:
            raise
        except Exception as exc:  # noqa: BLE001 - RDKit rejects a bad bond
            raise MoleculeConversionError(
                f"cannot add the bond {bond.u}-{bond.v}"
            ) from exc

    mol = editable.GetMol()
    try:
        chem.SanitizeMol(mol)
    except Exception as exc:  # noqa: BLE001 - RDKit rejects a bad valence
        raise MoleculeConversionError(f"RDKit sanitization failed: {exc}") from exc

    if value.positions_m:
        conformer = chem.Conformer(len(value.elements))
        for index, (x, y, z) in enumerate(value.positions_m):
            conformer.SetAtomPosition(
                index,
                point3d(
                    x * METRE_TO_ANGSTROM,
                    y * METRE_TO_ANGSTROM,
                    z * METRE_TO_ANGSTROM,
                ),
            )
        mol.AddConformer(conformer)
    return mol


def from_rdkit_mol(mol: Any, *, name: str = "") -> Molecule:
    """Convert an RDKit mol to a :class:`Molecule`.

    Positions come back in metres. A mol with no conformer gets 2D depiction
    coordinates, which are not physical.
    """
    chem, _ = _require_rdkit()
    copy = chem.Mol(mol)
    if copy.GetNumConformers() == 0:
        from rdkit.Chem import rdDepictor

        rdDepictor.Compute2DCoords(copy)
    conformer = copy.GetConformer()
    elements = tuple(atom.GetSymbol() for atom in copy.GetAtoms())
    positions = tuple(
        (
            conformer.GetAtomPosition(index).x * ANGSTROM_TO_METRE,
            conformer.GetAtomPosition(index).y * ANGSTROM_TO_METRE,
            conformer.GetAtomPosition(index).z * ANGSTROM_TO_METRE,
        )
        for index in range(copy.GetNumAtoms())
    )
    bonds = tuple(
        Bond(
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            int(round(bond.GetBondTypeAsDouble())),
        )
        for bond in copy.GetBonds()
    )
    return Molecule(elements=elements, positions_m=positions, bonds=bonds, name=name)


def read_smiles(text: str, *, add_hydrogens: bool = True) -> Molecule:
    """Parse a SMILES string into a :class:`Molecule`.

    When ``add_hydrogens`` is true, RDKit adds the implicit hydrogens as atoms,
    so the bond graph is complete. The positions are 2D depiction coordinates,
    not a physical geometry.
    """
    chem, _ = _require_rdkit()
    if not isinstance(text, str) or not text.strip():
        raise StructureParseError("the SMILES string is empty")
    mol = chem.MolFromSmiles(text)
    if mol is None:
        raise StructureParseError(f"RDKit cannot parse the SMILES {text!r}")
    if add_hydrogens:
        mol = chem.AddHs(mol)
    return from_rdkit_mol(mol)


def write_smiles(molecule: Any) -> str:
    """Write a :class:`Molecule` as a canonical SMILES string."""
    chem, _ = _require_rdkit()
    mol = to_rdkit_mol(molecule)
    try:
        without_hydrogens = chem.RemoveHs(mol)
        return chem.MolToSmiles(without_hydrogens)
    except Exception as exc:  # noqa: BLE001 - RDKit rejects what it cannot write
        raise MoleculeConversionError(f"cannot write SMILES: {exc}") from exc


def read_sdf(text: str) -> Molecule:
    """Read the first molecule of an SDF document."""
    chem, _ = _require_rdkit()
    if not isinstance(text, str) or not text.strip():
        raise StructureParseError("the SDF text is empty")
    supplier = chem.ForwardSDMolSupplier(BytesIO(text.encode("utf-8")), removeHs=False)
    for mol in supplier:
        if mol is not None:
            return from_rdkit_mol(mol)
    raise StructureParseError("the SDF text holds no readable molecule")


def write_sdf(molecule: Any) -> str:
    """Write a :class:`Molecule` as an SDF record."""
    chem, _ = _require_rdkit()
    mol = to_rdkit_mol(molecule)
    try:
        block = chem.MolToMolBlock(mol)
    except Exception as exc:  # noqa: BLE001 - RDKit rejects what it cannot write
        raise MoleculeConversionError(f"cannot write SDF: {exc}") from exc
    name = molecule.name if isinstance(molecule, Molecule) else ""
    if name:
        lines = block.splitlines()
        if lines:
            lines[0] = name
            block = "\n".join(lines) + "\n"
    return block + "$$$$\n"
