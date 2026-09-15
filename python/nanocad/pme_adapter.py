"""Particle-mesh Ewald (PME) electrostatics through OpenMM.

This module is the L1 engine backend for long-range electrostatics.
``ARCHITECTURE.md`` records that long-range electrostatics arrive later through
OpenMM.  This adapter adds that path.  It builds an OpenMM system with one
``NonbondedForce`` in PME mode and returns the potential energy, and
optionally the forces.

The adapter is independent of the compiled ``nanocad._core`` extension.  It
imports OpenMM lazily.  When OpenMM is absent, the module imports and every
entry point raises :class:`PmeBackendError`.

Units are SI at the boundary.  Positions are metres, length is metres, and
charge is the elementary charge.  OpenMM uses nanometres, the elementary
charge, and kilojoule per mole internally.  This module converts at the
boundary.  Every returned dictionary key carries a unit suffix.

Why the cutoff model differs
    The engine electrostatic term in ``crates/engine/src/electrostatic.rs``
    sums only atom pairs inside a cutoff.  It applies the minimum-image
    convention and a polynomial switching function.  That model omits every
    interaction beyond the cutoff, including the long-range tail and the
    periodic images beyond the first.  PME sums every periodic image and the
    reciprocal-space part.  The difference is the omitted long-range energy.
    The two results agree as the periodic box grows, because the nearest
    images move away.

Honest limits
    OpenMM is a third-party engine.  This adapter does not reimplement it.
    PME needs a neutral system and a periodic box.  A charged system without a
    counter charge gives a background term.  The comparison to the cutoff
    model measures a difference; it does not by itself prove either value
    correct.  The forces come from OpenMM and are reported as
    ``cross-checked`` only when a caller also checks them.
"""

from __future__ import annotations

import importlib.util
import math
from collections.abc import Sequence
from typing import Any

import numpy as np

from .continuum_adapter import _ATOMIC_MASS_AMU, _SYMBOL_BY_Z

__all__ = [
    "AVOGADRO_PER_MOL",
    "COULOMB_CONSTANT_N_M2_PER_C2",
    "ELEMENTARY_CHARGE_C",
    "KJ_PER_MOL_TO_J",
    "M_TO_NM",
    "PmeBackendError",
    "available",
    "compare_cutoff_and_pme",
    "cutoff_electrostatic_energy_j",
    "electrostatic_energy_pme",
]

# CODATA 2018 reference constants.
ELEMENTARY_CHARGE_C = 1.602176634e-19
AVOGADRO_PER_MOL = 6.02214076e23
COULOMB_CONSTANT_N_M2_PER_C2 = 8.9875517923e9

KJ_PER_MOL_TO_J = 1000.0 / AVOGADRO_PER_MOL
M_TO_NM = 1.0e9
KJ_PER_MOL_PER_NM_TO_N = KJ_PER_MOL_TO_J * M_TO_NM


class PmeBackendError(RuntimeError):
    """A PME backend is absent or failed to build a system."""


try:
    _OPENMM_AVAILABLE: bool = importlib.util.find_spec("openmm") is not None
except (ImportError, ValueError):  # pragma: no cover - unusual install state
    _OPENMM_AVAILABLE = False


def available() -> bool:
    """Return True when OpenMM is importable."""
    if not _OPENMM_AVAILABLE:
        return False
    try:
        import openmm  # noqa: F401
    except ImportError:
        return False
    return True


def _load_openmm() -> Any:
    """Import OpenMM lazily. Raise when it is absent."""
    if not _OPENMM_AVAILABLE:
        raise PmeBackendError(
            "OpenMM is not importable; install openmm to use the PME adapter"
        )
    try:
        import openmm
        import openmm.unit as unit
    except ImportError as error:  # pragma: no cover - guarded by find_spec
        raise PmeBackendError(f"OpenMM import failed: {error}") from error
    return openmm, unit


def _positive(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be a positive finite number")
    return value


def _normalize_positions_m(
    positions_m: Sequence[float] | Sequence[Sequence[float]], atom_count: int
) -> np.ndarray:
    array = np.asarray(positions_m, dtype=float)
    if array.size == 3 * atom_count and array.ndim == 1:
        array = array.reshape(-1, 3)
    if array.shape != (atom_count, 3):
        raise ValueError(
            f"positions must have shape ({atom_count}, 3) or a flat length "
            f"{3 * atom_count}"
        )
    if not np.all(np.isfinite(array)):
        raise ValueError("positions must be finite")
    return array


def _normalize_charges_e(charges_e: Sequence[float], atom_count: int) -> np.ndarray:
    charges = np.asarray(charges_e, dtype=float).reshape(-1)
    if charges.shape[0] != atom_count:
        raise ValueError("charge count must equal the atom count")
    if not np.all(np.isfinite(charges)):
        raise ValueError("charges must be finite")
    return charges


def _normalize_box_m(box_length_m: float | Sequence[float]) -> list[float]:
    if isinstance(box_length_m, (int, float)):
        lengths = [float(box_length_m)] * 3
    else:
        lengths = [float(value) for value in box_length_m]
        if len(lengths) != 3:
            raise ValueError("box_length_m must be a scalar or three lengths")
    lengths = [_positive(length, "box_length_m") for length in lengths]
    return lengths


def _mass_amu(element: Any) -> float:
    if isinstance(element, (int, np.integer)) and not isinstance(element, bool):
        symbol = _SYMBOL_BY_Z.get(int(element))
        if symbol is None:
            raise ValueError(f"unknown atomic number {element}")
    else:
        symbol = str(element)
    mass = _ATOMIC_MASS_AMU.get(symbol)
    if mass is None:
        raise ValueError(f"unknown element {element!r}")
    return float(mass)


def _to_float(value: Any, unit: Any) -> float:
    try:
        return float(value.value_in_unit(unit))
    except AttributeError:
        return float(value)


def electrostatic_energy_pme(
    elements: Sequence[str | int],
    positions_m: Sequence[float] | Sequence[Sequence[float]],
    charges_e: Sequence[float],
    box_length_m: float | Sequence[float],
    cutoff_m: float,
    periodic: bool = True,
    ewald_error_tolerance: float = 5.0e-4,
    with_forces: bool = False,
) -> dict:
    """Return the PME electrostatic energy of a periodic atom box.

    ``elements`` holds element symbols or atomic numbers.  ``positions_m``
    holds positions in metres, as ``[x, y, z]`` triplets or one flat list.
    ``charges_e`` holds partial charges in units of the elementary charge.
    ``box_length_m`` is a scalar or three orthorhombic box lengths in metres.
    ``cutoff_m`` is the real-space cutoff in metres.

    The system uses one OpenMM ``NonbondedForce`` in PME mode on the
    ``Reference`` platform, so the result is deterministic.  Lennard-Jones is
    switched off with a zero epsilon, so the energy is purely electrostatic.

    Returns a dictionary with ``energy_j`` and, when ``with_forces`` is true,
    a flat ``forces_n`` list.  ``method`` is ``openmm_pme``.  ``validation``
    is ``unverified``: the caller must compare the value.

    When OpenMM is absent this raises :class:`PmeBackendError`.
    """
    openmm, unit = _load_openmm()
    elements = list(elements)
    atom_count = len(elements)
    if atom_count == 0:
        raise ValueError("at least one atom is required")
    positions = _normalize_positions_m(positions_m, atom_count)
    charges = _normalize_charges_e(charges_e, atom_count)
    box = _normalize_box_m(box_length_m)
    cutoff_m = _positive(cutoff_m, "cutoff_m")
    ewald_error_tolerance = _positive(ewald_error_tolerance, "ewald_error_tolerance")
    if periodic:
        for length_m in box:
            if cutoff_m > 0.5 * length_m:
                raise ValueError("cutoff_m must not exceed half of any box length")

    try:
        system = openmm.System()
        for element in elements:
            system.addParticle(_mass_amu(element))
        if periodic:
            system.setDefaultPeriodicBoxVectors(
                openmm.Vec3(box[0] * M_TO_NM, 0.0, 0.0) * unit.nanometer,
                openmm.Vec3(0.0, box[1] * M_TO_NM, 0.0) * unit.nanometer,
                openmm.Vec3(0.0, 0.0, box[2] * M_TO_NM) * unit.nanometer,
            )
        force = openmm.NonbondedForce()
        if periodic:
            force.setNonbondedMethod(openmm.NonbondedForce.PME)
        else:
            force.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)
        force.setCutoffDistance(cutoff_m * M_TO_NM * unit.nanometer)
        force.setEwaldErrorTolerance(ewald_error_tolerance)
        for charge in charges:
            force.addParticle(float(charge), 0.3, 0.0)
        system.addForce(force)
        platform = openmm.Platform.getPlatformByName("Reference")
        context = openmm.Context(system, openmm.VerletIntegrator(0.001), platform)
        context.setPositions(
            [openmm.Vec3(*(row * M_TO_NM)) * unit.nanometer for row in positions]
        )
        state = context.getState(getEnergy=True, getForces=with_forces)
        energy_kj_mol = _to_float(state.getPotentialEnergy(), unit.kilojoule_per_mole)
    except PmeBackendError:
        raise
    except Exception as error:
        raise PmeBackendError(f"OpenMM PME failed: {error}") from error

    result: dict = {
        "atom_count": atom_count,
        "net_charge_e": float(charges.sum()),
        "box_length_m": box,
        "cutoff_m": cutoff_m,
        "periodic": bool(periodic),
        "ewald_error_tolerance": ewald_error_tolerance,
        "energy_kj_mol": energy_kj_mol,
        "energy_j": energy_kj_mol * KJ_PER_MOL_TO_J,
        "energy_per_atom_j": energy_kj_mol * KJ_PER_MOL_TO_J / atom_count,
        "method": "openmm_pme" if periodic else "openmm_direct_no_cutoff",
        "backend": "openmm",
        "platform": "Reference",
        "validation": "unverified",
    }
    if periodic:
        try:
            pme_parameters = force.getPMEParametersInContext(context)
            result["pme_alpha_per_nm"] = _to_float(
                pme_parameters[0], 1.0 / unit.nanometer
            )
            result["pme_grid"] = [
                int(pme_parameters[1]),
                int(pme_parameters[2]),
                int(pme_parameters[3]),
            ]
        except Exception:  # pragma: no cover - older OpenMM builds
            result["pme_alpha_per_nm"] = None
            result["pme_grid"] = None
    if with_forces:
        flat_forces_n: list[float] = []
        for vector in state.getForces():
            for axis in range(3):
                flat_forces_n.append(
                    _to_float(vector[axis], unit.kilojoule_per_mole / unit.nanometer)
                    * KJ_PER_MOL_PER_NM_TO_N
                )
        result["forces_n"] = flat_forces_n
    return result


def cutoff_electrostatic_energy_j(
    elements: Sequence[str | int],
    positions_m: Sequence[float] | Sequence[Sequence[float]],
    charges_e: Sequence[float],
    box_length_m: float | Sequence[float],
    cutoff_m: float,
    switch_on_fraction: float = 0.9,
) -> dict:
    """Return the engine cutoff Coulomb energy for a periodic atom box.

    This calls the compiled ``nanocad._core.System``.  It is the comparison
    baseline for :func:`electrostatic_energy_pme`.  The engine uses the
    minimum-image convention and a switching function, and it sums only pairs
    inside the cutoff.  ``switch_on_fraction`` sets the switch-on distance as
    a fraction of the cutoff.

    Raises :class:`PmeBackendError` when the compiled extension is absent.
    """
    elements = list(elements)
    atom_count = len(elements)
    if atom_count == 0:
        raise ValueError("at least one atom is required")
    positions = _normalize_positions_m(positions_m, atom_count)
    charges = _normalize_charges_e(charges_e, atom_count)
    box = _normalize_box_m(box_length_m)
    cutoff_m = _positive(cutoff_m, "cutoff_m")
    switch_on_fraction = float(switch_on_fraction)
    if not 0.0 <= switch_on_fraction < 1.0:
        raise ValueError("switch_on_fraction must be in [0, 1)")
    switch_on_m = switch_on_fraction * cutoff_m

    try:
        from . import _core
    except ImportError as error:  # pragma: no cover - extension absent
        raise PmeBackendError("nanocad._core is not importable") from error
    if _core is None:
        raise PmeBackendError("nanocad._core is not built")

    try:
        system = _core.System(atom_count)
        system.set_electrostatic(
            cutoff_m, switch_on_m, [float(q) * ELEMENTARY_CHARGE_C for q in charges]
        )
        system.set_periodic_box(box)
        energy_j = system.energy(positions.reshape(-1).tolist())
    except Exception as error:
        raise PmeBackendError(
            f"engine cutoff electrostatics failed: {error}"
        ) from error

    return {
        "atom_count": atom_count,
        "net_charge_e": float(charges.sum()),
        "box_length_m": box,
        "cutoff_m": cutoff_m,
        "switch_on_m": switch_on_m,
        "energy_j": float(energy_j),
        "method": "engine_cutoff_electrostatic",
        "backend": "nanocad._core",
        "validation": "unverified",
    }


def compare_cutoff_and_pme(
    elements: Sequence[str | int],
    positions_m: Sequence[float] | Sequence[Sequence[float]],
    charges_e: Sequence[float],
    box_length_m: float | Sequence[float],
    cutoff_m: float,
    switch_on_fraction: float = 0.9,
    ewald_error_tolerance: float = 5.0e-4,
) -> dict:
    """Compare the engine cutoff Coulomb energy with the OpenMM PME energy.

    Returns both energies, the signed difference, and the relative difference.
    The difference is the long-range contribution that the cutoff model
    omits, weighted by the switching function.  The value depends on the box
    size and the cutoff, so it is not a universal constant.
    """
    pme = electrostatic_energy_pme(
        elements,
        positions_m,
        charges_e,
        box_length_m,
        cutoff_m,
        periodic=True,
        ewald_error_tolerance=ewald_error_tolerance,
    )
    cutoff = cutoff_electrostatic_energy_j(
        elements,
        positions_m,
        charges_e,
        box_length_m,
        cutoff_m,
        switch_on_fraction=switch_on_fraction,
    )
    difference_j = pme["energy_j"] - cutoff["energy_j"]
    if cutoff["energy_j"] == 0.0:
        relative_difference = math.inf if difference_j != 0.0 else 0.0
    else:
        relative_difference = difference_j / abs(cutoff["energy_j"])
    return {
        "pme_energy_j": pme["energy_j"],
        "cutoff_energy_j": cutoff["energy_j"],
        "difference_j": difference_j,
        "relative_difference": relative_difference,
        "net_charge_e": pme["net_charge_e"],
        "box_length_m": pme["box_length_m"],
        "cutoff_m": pme["cutoff_m"],
        "pme": pme,
        "cutoff": cutoff,
        "method": "engine_cutoff_vs_openmm_pme",
        "validation": "cross-checked",
        "difference_explanation": (
            "The engine cutoff model sums only minimum-image pairs inside the "
            "cutoff and applies a switching function. PME sums every periodic "
            "image and the reciprocal-space part. The difference is the "
            "omitted long-range electrostatics; it shrinks as the box grows."
        ),
    }
