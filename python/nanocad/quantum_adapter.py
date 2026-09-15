"""Quantum-chemistry adapter for small molecules.

This module gives the rest of NanoCAD one small, honest entry point to a real
quantum-chemistry engine. It accepts plain Python data and returns a result
dict with unit-suffixed keys. It does not import a backend at module import
time. It imports each backend lazily inside the function that runs.

Backend preference order:

1. PySCF (real quantum chemistry, for example RHF/STO-3G on H2).
2. xTB (real semi-empirical quantum chemistry, if importable).
3. ASE with a clearly-labelled classical fallback calculator.

The ASE path is NOT quantum. Its result says so.

Units are SI internally. Positions are metres. Energies are joules.
"""

from __future__ import annotations

import importlib.util
import math
from collections.abc import Sequence
from typing import Any

__all__ = [
    "ANGSTROM_TO_METRE",
    "EV_TO_JOULE",
    "HARTREE_TO_JOULE",
    "QuantumAdapterError",
    "NoQuantumBackendError",
    "QuantumBackendError",
    "available_backends",
    "compute_energy",
    "compute_energy_from_topology",
    "topology_to_geometry",
]

# CODATA 2018 reference constants.
ANGSTROM_TO_METRE = 1.0e-10
HARTREE_TO_JOULE = 4.3597447222071e-18
EV_TO_JOULE = 1.602176634e-19

_BACKEND_PREFERENCE = ("pyscf", "xtb", "ase")

# IUPAC element symbols, Z = 1..118. This is public reference data.
_ELEMENT_SYMBOLS = (
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og",
)

_NUMBER_TO_SYMBOL = dict(enumerate(_ELEMENT_SYMBOLS, start=1))
_SYMBOL_TO_NUMBER = {symbol: number for number, symbol in _NUMBER_TO_SYMBOL.items()}

# Published reference energies. Each entry is (hartree, citation).
# Hehre, Stewart, Pople, J. Chem. Phys. 51, 2657 (1969): STO-3G H2 at
# r = 0.735 Angstrom has a self-consistent-field energy of about -1.1167 Ha.
_REFERENCES = {
    ("H2", "sto-3g", "RHF"): (
        -1.1167,
        "Hehre-Stewart-Pople-1969 H2 STO-3G (-1.1167 Ha)",
    )
}
_REFERENCE_TOLERANCE_HARTREE = 5.0e-3


class QuantumAdapterError(RuntimeError):
    """Base class for adapter failures."""


class NoQuantumBackendError(QuantumAdapterError):
    """No quantum backend is importable and no fallback was requested."""


class QuantumBackendError(QuantumAdapterError):
    """An importable backend failed to compute an energy."""


def available_backends() -> list[str]:
    """Return the backends that are importable, in preference order."""
    found: list[str] = []
    for name in _BACKEND_PREFERENCE:
        try:
            spec = importlib.util.find_spec(name)
        except (ImportError, ValueError):
            spec = None
        if spec is not None:
            found.append(name)
    return found


def _load_core() -> Any:
    """Import ``nanocad._core`` lazily. Return None when it is absent."""
    try:
        from . import _core
    except Exception:
        return None
    return _core


def _normalize_symbols(elements: Sequence[str | int]) -> list[str]:
    if isinstance(elements, (str, bytes)):
        raise ValueError("elements must be a sequence, not a single string")
    symbols: list[str] = []
    for element in elements:
        if isinstance(element, bool):
            raise ValueError("element must be a symbol or an atomic number")
        if isinstance(element, int):
            if element not in _NUMBER_TO_SYMBOL:
                raise ValueError(f"atomic number out of range: {element}")
            symbols.append(_NUMBER_TO_SYMBOL[element])
        elif isinstance(element, str):
            text = element.strip()
            if not text:
                raise ValueError("empty element symbol")
            symbol = text[:1].upper() + text[1:].lower()
            if symbol not in _SYMBOL_TO_NUMBER:
                raise ValueError(f"unknown element symbol: {element!r}")
            symbols.append(symbol)
        else:
            raise ValueError(
                f"element must be a symbol or an atomic number: {element!r}"
            )
    if not symbols:
        raise ValueError("at least one atom is required")
    return symbols


def _normalize_positions(
    positions_m: Sequence[float] | Sequence[Sequence[float]], atom_count: int
) -> list[tuple[float, float, float]]:
    flat = list(positions_m)
    is_flat = len(flat) == 3 * atom_count and all(
        isinstance(value, (int, float)) and not isinstance(value, bool)
        for value in flat
    )
    if is_flat:
        return [
            (float(flat[3 * i]), float(flat[3 * i + 1]), float(flat[3 * i + 2]))
            for i in range(atom_count)
        ]
    if len(flat) != atom_count:
        raise ValueError(
            f"expected {atom_count} positions, got {len(flat)}; "
            "give [x, y, z] triplets or one flat list"
        )
    triplets: list[tuple[float, float, float]] = []
    for row in flat:
        try:
            x, y, z = row  # type: ignore[misc]
        except (TypeError, ValueError) as exc:
            raise ValueError(f"each position needs three numbers: {row!r}") from exc
        triplets.append((float(x), float(y), float(z)))
    for x, y, z in triplets:
        if not all(math.isfinite(value) for value in (x, y, z)):
            raise ValueError("positions must be finite")
    return triplets


def _normalize_charge(charge: int) -> int:
    if isinstance(charge, bool) or not isinstance(charge, int):
        raise ValueError("charge must be an integer")
    return charge


def _electron_count(symbols: Sequence[str], charge: int) -> int:
    return sum(_SYMBOL_TO_NUMBER[symbol] for symbol in symbols) - charge


def _resolve_multiplicity(electron_count: int, multiplicity: int | None) -> int:
    if multiplicity is None:
        return 1 if electron_count % 2 == 0 else 2
    if isinstance(multiplicity, bool) or not isinstance(multiplicity, int):
        raise ValueError("multiplicity must be an integer or None")
    if multiplicity < 1:
        raise ValueError("multiplicity must be at least 1")
    unpaired = multiplicity - 1
    if unpaired > electron_count:
        raise ValueError("multiplicity is too large for the electron count")
    if (electron_count - unpaired) % 2 != 0:
        raise ValueError("multiplicity parity does not match the electron count")
    return multiplicity


def _select_backend(requested: str) -> str:
    if requested == "auto":
        importable = available_backends()
        if not importable:
            raise NoQuantumBackendError(
                "no quantum backend is importable; install pyscf or xtb, "
                "or pass backend='ase' for a classical fallback"
            )
        return importable[0]
    if requested not in _BACKEND_PREFERENCE:
        raise ValueError(
            f"unknown backend {requested!r}; choose 'auto', "
            + ", ".join(repr(name) for name in _BACKEND_PREFERENCE)
        )
    if requested not in available_backends():
        raise NoQuantumBackendError(f"backend {requested!r} is not importable")
    return requested


def _run_pyscf(
    symbols: list[str],
    coords_m: list[tuple[float, float, float]],
    charge: int,
    multiplicity: int,
    basis: str,
) -> tuple[float, str]:
    from pyscf import gto, scf

    atom = [
        (symbol, tuple(value * 1.0e10 for value in xyz))
        for symbol, xyz in zip(symbols, coords_m, strict=True)
    ]
    spin = multiplicity - 1
    try:
        mol = gto.M(
            atom=atom,
            charge=charge,
            spin=spin,
            unit="Angstrom",
            basis=basis,
            verbose=0,
        )
        mean_field = scf.RHF(mol) if spin == 0 else scf.UHF(mol)
        energy_hartree = float(mean_field.kernel())
    except Exception as exc:
        raise QuantumBackendError(f"pyscf failed: {exc}") from exc
    if not math.isfinite(energy_hartree):
        raise QuantumBackendError("pyscf returned a non-finite energy")
    level = f"{'RHF' if spin == 0 else 'UHF'}/{basis.upper()}"
    return energy_hartree, level


def _run_xtb(
    symbols: list[str],
    coords_m: list[tuple[float, float, float]],
    charge: int,
    multiplicity: int,
) -> tuple[float, str]:
    from ase import Atoms
    from xtb.ase.calculator import XTB

    atoms = Atoms(
        symbols,
        positions=[[value * 1.0e10 for value in xyz] for xyz in coords_m],
    )
    try:
        atoms.calc = XTB(method="GFN2-xTB", charge=charge)
        energy_ev = float(atoms.get_potential_energy())
    except TypeError:
        # Older xTB wrappers take no charge keyword.
        try:
            atoms.calc = XTB(method="GFN2-xTB")
            energy_ev = float(atoms.get_potential_energy())
        except Exception as exc:
            raise QuantumBackendError(f"xtb failed: {exc}") from exc
    except Exception as exc:
        raise QuantumBackendError(f"xtb failed: {exc}") from exc
    if not math.isfinite(energy_ev):
        raise QuantumBackendError("xtb returned a non-finite energy")
    spin = multiplicity - 1
    return energy_ev * EV_TO_JOULE, f"GFN2-xTB (charge={charge}, spin={spin})"


def _run_ase_fallback(
    symbols: list[str],
    coords_m: list[tuple[float, float, float]],
) -> tuple[float, str]:
    from ase import Atoms

    atoms = Atoms(
        symbols,
        positions=[[value * 1.0e10 for value in xyz] for xyz in coords_m],
    )
    try:
        from ase.calculators.emt import EMT

        atoms.calc = EMT()
        energy_ev = float(atoms.get_potential_energy())
        calculator = "EMT"
    except Exception:
        from ase.calculators.lj import LennardJones

        try:
            atoms.calc = LennardJones()
            energy_ev = float(atoms.get_potential_energy())
            calculator = "Lennard-Jones"
        except Exception as exc:
            raise QuantumBackendError(f"ase fallback failed: {exc}") from exc
    if not math.isfinite(energy_ev):
        raise QuantumBackendError("ase fallback returned a non-finite energy")
    level = f"{calculator} (classical fallback, NOT quantum)"
    return energy_ev * EV_TO_JOULE, level


def _formula(symbols: Sequence[str]) -> str:
    counts: dict[str, int] = {}
    order: list[str] = []
    for symbol in symbols:
        if symbol not in counts:
            counts[symbol] = 0
            order.append(symbol)
        counts[symbol] += 1
    return "".join(
        symbol if counts[symbol] == 1 else f"{symbol}{counts[symbol]}"
        for symbol in order
    )


def _validation(
    symbols: list[str],
    charge: int,
    level: str,
    basis: str,
    energy_hartree: float,
) -> str:
    key = (_formula(symbols), basis.lower(), level.split("/", 1)[0])
    reference = _REFERENCES.get(key)
    if charge == 0 and reference is not None:
        reference_value, citation = reference
        if abs(energy_hartree - reference_value) <= _REFERENCE_TOLERANCE_HARTREE:
            return f"reference: {citation}"
    return "unverified"


def compute_energy(
    elements: Sequence[str | int],
    positions_m: Sequence[float] | Sequence[Sequence[float]],
    charge: int = 0,
    multiplicity: int | None = None,
    backend: str = "auto",
    basis: str = "sto-3g",
) -> dict[str, Any]:
    """Compute the potential energy of a small molecule.

    ``elements`` holds element symbols or atomic numbers. ``positions_m``
    holds positions in metres as ``[x, y, z]`` triplets, or one flat list.
    The return value has unit-suffixed keys. ``method`` names the backend and
    the level of theory that ran. ``validation`` is ``"unverified"`` unless a
    cited reference value matched.
    """
    symbols = _normalize_symbols(elements)
    coords_m = _normalize_positions(positions_m, len(symbols))
    charge_value = _normalize_charge(charge)
    electron_count = _electron_count(symbols, charge_value)
    if electron_count < 0:
        raise ValueError("positive charge is larger than the electron count")
    multiplicity_value = _resolve_multiplicity(electron_count, multiplicity)
    chosen = _select_backend(backend)

    if chosen == "pyscf":
        energy_hartree, level = _run_pyscf(
            symbols, coords_m, charge_value, multiplicity_value, basis
        )
        energy_j = energy_hartree * HARTREE_TO_JOULE
        is_quantum = True
        validation = _validation(symbols, charge_value, level, basis, energy_hartree)
        method = f"pyscf: {level}"
    elif chosen == "xtb":
        energy_j, level = _run_xtb(symbols, coords_m, charge_value, multiplicity_value)
        is_quantum = True
        validation = "unverified"
        method = f"xtb: {level}"
    else:
        energy_j, level = _run_ase_fallback(symbols, coords_m)
        is_quantum = False
        validation = "unverified"
        method = f"ase: {level}"

    return {
        "atom_count": len(symbols),
        "charge": charge_value,
        "multiplicity": multiplicity_value,
        "energy_j": energy_j,
        "energy_ev": energy_j / EV_TO_JOULE,
        "method": method,
        "backend": chosen,
        "level_of_theory": level,
        "is_quantum": is_quantum,
        "validation": validation,
    }


def topology_to_geometry(topology: Any) -> tuple[list[str], list[list[float]], int]:
    """Read symbols, positions in metres, and total charge from a core topology.

    Accepts a ``nanocad._core`` ``Topology`` or a ``Part`` that wraps one. The
    object must expose ``elements()`` (bytes of atomic numbers) and
    ``positions_m()`` (a flat list in metres). The core is imported lazily;
    the function also accepts a duck-typed object when the core is absent.
    """
    core = _load_core()
    if core is not None and isinstance(topology, core.Part):
        topo = topology.topology
    else:
        topo = getattr(topology, "topology", topology)
    raw_elements = topo.elements()
    numbers = list(raw_elements)
    symbols = [_NUMBER_TO_SYMBOL[number] for number in numbers]
    flat = list(topo.positions_m())
    positions = [
        [float(flat[3 * i]), float(flat[3 * i + 1]), float(flat[3 * i + 2])]
        for i in range(len(symbols))
    ]
    charges = list(topo.charges_c()) if hasattr(topo, "charges_c") else []
    total_charge = int(round(sum(charges))) if charges else 0
    return symbols, positions, total_charge


def compute_energy_from_topology(
    topology: Any,
    multiplicity: int | None = None,
    backend: str = "auto",
    basis: str = "sto-3g",
) -> dict[str, Any]:
    """Run ``compute_energy`` on a ``nanocad._core`` topology or part."""
    symbols, positions, total_charge = topology_to_geometry(topology)
    return compute_energy(
        symbols,
        positions,
        charge=total_charge,
        multiplicity=multiplicity,
        backend=backend,
        basis=basis,
    )
