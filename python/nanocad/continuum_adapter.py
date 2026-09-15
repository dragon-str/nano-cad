"""Continuum adapter: lift atomistic or part data to continuum models.

Pure Python plus NumPy.  This module does not import ``nanocad._core`` at
import time, so a plain atom list works without the compiled extension.  The
core is imported lazily only when a caller passes a raw core object.

Units are SI.  Every returned dictionary key carries a unit suffix.

Structural models
    * Axial bar with linear finite elements.
    * Euler-Bernoulli cantilever beam with Hermite finite elements.

Flow models
    * Hagen-Poiseuille flow in a circular tube.
    * Plane Poiseuille flow in a straight channel.

The structural models compare the finite-element result with the exact closed
form.  The flow model compares a finite-difference solution with the exact
Hagen-Poiseuille or plane-Poiseuille result.

Every result carries ``method`` and ``validation``.  A result is
``cross-checked`` only when the code compares it with the exact analytic
solution.  Otherwise it is ``unverified``.
"""

from __future__ import annotations

import math

import numpy as np

AMU_TO_KG = 1.66053906660e-27

_ATOMIC_MASS_AMU: dict[str, float] = {
    "H": 1.008,
    "He": 4.0026,
    "Li": 6.94,
    "Be": 9.0122,
    "B": 10.81,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "F": 18.998,
    "Ne": 20.180,
    "Na": 22.990,
    "Mg": 24.305,
    "Al": 26.982,
    "Si": 28.085,
    "P": 30.974,
    "S": 32.06,
    "Cl": 35.45,
    "Ar": 39.948,
    "K": 39.098,
    "Ca": 40.078,
    "Ti": 47.867,
    "Cr": 51.996,
    "Mn": 54.938,
    "Fe": 55.845,
    "Co": 58.933,
    "Ni": 58.693,
    "Cu": 63.546,
    "Zn": 65.38,
    "Ge": 72.630,
    "As": 74.922,
    "Se": 78.971,
    "Br": 79.904,
    "Zr": 91.224,
    "Mo": 95.95,
    "Ag": 107.868,
    "Sn": 118.710,
    "I": 126.904,
    "Xe": 131.293,
    "W": 183.84,
    "Pt": 195.084,
    "Au": 196.967,
    "Hg": 200.592,
    "Pb": 207.2,
    "U": 238.029,
}

_SYMBOL_BY_Z: dict[int, str] = {
    1: "H",
    2: "He",
    3: "Li",
    4: "Be",
    5: "B",
    6: "C",
    7: "N",
    8: "O",
    9: "F",
    10: "Ne",
    11: "Na",
    12: "Mg",
    13: "Al",
    14: "Si",
    15: "P",
    16: "S",
    17: "Cl",
    18: "Ar",
    19: "K",
    20: "Ca",
    22: "Ti",
    24: "Cr",
    25: "Mn",
    26: "Fe",
    27: "Co",
    28: "Ni",
    29: "Cu",
    30: "Zn",
    32: "Ge",
    33: "As",
    34: "Se",
    35: "Br",
    40: "Zr",
    42: "Mo",
    47: "Ag",
    50: "Sn",
    53: "I",
    54: "Xe",
    74: "W",
    78: "Pt",
    79: "Au",
    80: "Hg",
    82: "Pb",
    92: "U",
}


def _require_positive(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be a positive finite number")
    return value


def _relative_error(computed: float, exact: float) -> float:
    if exact == 0.0:
        return abs(computed)
    return abs(computed - exact) / abs(exact)


def _assemble_bar_stiffness(
    elastic_modulus_pa: float, area_m2: float, length_m: float, elements: int
) -> np.ndarray:
    nodes = elements + 1
    element_length_m = length_m / elements
    stiffness_n_per_m = elastic_modulus_pa * area_m2 / element_length_m
    matrix = np.zeros((nodes, nodes))
    for index in range(elements):
        matrix[index, index] += stiffness_n_per_m
        matrix[index, index + 1] -= stiffness_n_per_m
        matrix[index + 1, index] -= stiffness_n_per_m
        matrix[index + 1, index + 1] += stiffness_n_per_m
    return matrix


def solve_axial_bar(
    length_m: float,
    area_m2: float,
    elastic_modulus_pa: float,
    axial_load_n: float,
    density_kg_m3: float,
    elements: int = 20,
) -> dict:
    """Solve a 1D axial bar with a fixed end and a tip axial load.

    Units: ``length_m`` [m], ``area_m2`` [m^2], ``elastic_modulus_pa`` [Pa],
    ``axial_load_n`` [N] (tension positive), ``density_kg_m3`` [kg/m^3].

    Returns axial displacement and stress in SI units.  The result is
    ``cross-checked`` against ``delta = P L / (A E)``.
    """
    length_m = _require_positive(length_m, "length_m")
    area_m2 = _require_positive(area_m2, "area_m2")
    elastic_modulus_pa = _require_positive(elastic_modulus_pa, "elastic_modulus_pa")
    density_kg_m3 = _require_positive(density_kg_m3, "density_kg_m3")
    if elements < 1:
        raise ValueError("elements must be at least 1")

    stiffness = _assemble_bar_stiffness(elastic_modulus_pa, area_m2, length_m, elements)
    free = np.arange(1, elements + 1)
    load_n = np.zeros(elements + 1)
    load_n[-1] = float(axial_load_n)
    displacement_m = np.zeros(elements + 1)
    displacement_m[free] = np.linalg.solve(stiffness[np.ix_(free, free)], load_n[free])

    element_length_m = length_m / elements
    strain = np.diff(displacement_m) / element_length_m
    stress_pa = elastic_modulus_pa * strain

    analytic_tip_m = axial_load_n * length_m / (area_m2 * elastic_modulus_pa)
    tip_m = float(displacement_m[-1])
    relative_error = _relative_error(tip_m, analytic_tip_m)

    mass_kg = density_kg_m3 * area_m2 * length_m
    wave_speed_m_s = math.sqrt(elastic_modulus_pa / density_kg_m3)
    frequency_hz = wave_speed_m_s / (4.0 * length_m)

    return {
        "length_m": length_m,
        "area_m2": area_m2,
        "elastic_modulus_pa": elastic_modulus_pa,
        "axial_load_n": float(axial_load_n),
        "density_kg_m3": density_kg_m3,
        "elements": elements,
        "tip_displacement_m": tip_m,
        "axial_stress_pa": float(stress_pa[0]),
        "axial_strain": float(strain[0]),
        "stress_profile_pa": stress_pa.tolist(),
        "stiffness_n_per_m": elastic_modulus_pa * area_m2 / length_m,
        "mass_kg": mass_kg,
        "first_natural_frequency_hz": frequency_hz,
        "analytic_tip_displacement_m": analytic_tip_m,
        "relative_error": relative_error,
        "converged": relative_error < 1.0e-10,
        "method": "axial_bar_linear_fem",
        "validation": "cross-checked",
    }


def solve_euler_bernoulli(
    length_m: float,
    elastic_modulus_pa: float,
    second_moment_of_area_m4: float,
    tip_load_n: float,
    outer_fiber_distance_m: float,
    density_kg_m3: float | None = None,
    area_m2: float | None = None,
    elements: int = 20,
) -> dict:
    """Solve a cantilever Euler-Bernoulli beam with a tip transverse load.

    Units: ``length_m`` [m], ``elastic_modulus_pa`` [Pa],
    ``second_moment_of_area_m4`` [m^4], ``tip_load_n`` [N],
    ``outer_fiber_distance_m`` [m] is the distance from the neutral axis to
    the extreme fiber (for a rectangle, half the height).

    Returns tip deflection and maximum bending stress.  The result is
    ``cross-checked`` against ``delta = P L^3 / (3 E I)`` and
    ``sigma = P L c / I``.
    """
    length_m = _require_positive(length_m, "length_m")
    elastic_modulus_pa = _require_positive(elastic_modulus_pa, "elastic_modulus_pa")
    inertia_m4 = _require_positive(second_moment_of_area_m4, "second_moment_of_area_m4")
    outer_fiber_m = _require_positive(outer_fiber_distance_m, "outer_fiber_distance_m")
    if elements < 1:
        raise ValueError("elements must be at least 1")

    nodes = elements + 1
    dofs = 2 * nodes
    element_length_m = length_m / elements
    factor = elastic_modulus_pa * inertia_m4 / element_length_m**3
    local = factor * np.array(
        [
            [12.0, 6.0 * element_length_m, -12.0, 6.0 * element_length_m],
            [
                6.0 * element_length_m,
                4.0 * element_length_m**2,
                -6.0 * element_length_m,
                2.0 * element_length_m**2,
            ],
            [-12.0, -6.0 * element_length_m, 12.0, -6.0 * element_length_m],
            [
                6.0 * element_length_m,
                2.0 * element_length_m**2,
                -6.0 * element_length_m,
                4.0 * element_length_m**2,
            ],
        ]
    )
    stiffness = np.zeros((dofs, dofs))
    for index in range(elements):
        map_index = [2 * index, 2 * index + 1, 2 * index + 2, 2 * index + 3]
        stiffness[np.ix_(map_index, map_index)] += local

    load_n = np.zeros(dofs)
    load_n[2 * (nodes - 1)] = float(tip_load_n)
    free = np.arange(2, dofs)
    displacement = np.zeros(dofs)
    displacement[free] = np.linalg.solve(stiffness[np.ix_(free, free)], load_n[free])
    reactions = stiffness @ displacement - load_n
    root_moment_n_m = abs(float(reactions[1]))

    analytic_tip_m = tip_load_n * length_m**3 / (3.0 * elastic_modulus_pa * inertia_m4)
    tip_m = float(displacement[2 * (nodes - 1)])
    analytic_stress_pa = abs(tip_load_n) * length_m * outer_fiber_m / inertia_m4
    stress_pa = root_moment_n_m * outer_fiber_m / inertia_m4
    relative_error = _relative_error(tip_m, analytic_tip_m)

    result: dict = {
        "length_m": length_m,
        "elastic_modulus_pa": elastic_modulus_pa,
        "second_moment_of_area_m4": inertia_m4,
        "tip_load_n": float(tip_load_n),
        "outer_fiber_distance_m": outer_fiber_m,
        "elements": elements,
        "tip_deflection_m": tip_m,
        "max_bending_stress_pa": float(stress_pa),
        "root_moment_n_m": root_moment_n_m,
        "analytic_tip_deflection_m": analytic_tip_m,
        "analytic_max_bending_stress_pa": analytic_stress_pa,
        "relative_error": relative_error,
        "converged": relative_error < 1.0e-10,
        "method": "euler_bernoulli_hermite_fem",
        "validation": "cross-checked",
    }
    if density_kg_m3 is not None and area_m2 is not None:
        density_kg_m3 = _require_positive(density_kg_m3, "density_kg_m3")
        area_m2 = _require_positive(area_m2, "area_m2")
        beta = 1.875104068711961
        frequency_hz = (
            beta**2
            / (2.0 * math.pi * length_m**2)
            * math.sqrt(elastic_modulus_pa * inertia_m4 / (density_kg_m3 * area_m2))
        )
        result["mass_kg"] = density_kg_m3 * area_m2 * length_m
        result["first_natural_frequency_hz"] = frequency_hz
    return result


def _solve_tridiagonal(
    lower: np.ndarray, diag: np.ndarray, upper: np.ndarray, rhs: np.ndarray
) -> np.ndarray:
    count = diag.shape[0]
    upper = upper.astype(float, copy=True)
    diag = diag.astype(float, copy=True)
    rhs = rhs.astype(float, copy=True)
    lower = lower.astype(float, copy=True)
    for index in range(1, count):
        weight = lower[index] / diag[index - 1]
        diag[index] -= weight * upper[index - 1]
        rhs[index] -= weight * rhs[index - 1]
    solution = np.empty(count)
    solution[-1] = rhs[-1] / diag[-1]
    for index in range(count - 2, -1, -1):
        solution[index] = (rhs[index] - upper[index] * solution[index + 1]) / diag[
            index
        ]
    return solution


def _poiseuille_profile(
    geometry: str,
    extent_m: float,
    gradient_pa_per_m: float,
    viscosity_pa_s: float,
    bins: int,
) -> tuple[np.ndarray, np.ndarray]:
    step_m = extent_m / bins
    coordinate_m = np.arange(bins) * step_m
    rhs = np.full(bins, -gradient_pa_per_m / viscosity_pa_s)
    lower = np.zeros(bins)
    diag = np.zeros(bins)
    upper = np.zeros(bins)

    inverse_step_sq = 1.0 / step_m**2
    if geometry == "tube":
        for index in range(bins):
            radius_m = coordinate_m[index]
            diag[index] = -2.0 * inverse_step_sq
            if index == 0:
                diag[index] = -4.0 * inverse_step_sq
                upper[index] = 4.0 * inverse_step_sq
            elif index < bins - 1:
                convection = 1.0 / (2.0 * step_m * radius_m)
                lower[index] = inverse_step_sq - convection
                upper[index] = inverse_step_sq + convection
            else:
                lower[index] = inverse_step_sq - 1.0 / (2.0 * step_m * radius_m)
    else:
        for index in range(bins):
            diag[index] = -2.0 * inverse_step_sq
            if index == 0:
                upper[index] = 2.0 * inverse_step_sq
            elif index < bins - 1:
                lower[index] = inverse_step_sq
                upper[index] = inverse_step_sq
            else:
                lower[index] = inverse_step_sq

    interior = _solve_tridiagonal(lower, diag, upper, rhs)
    coordinate_full = np.arange(bins + 1) * step_m
    profile = np.concatenate([interior, [0.0]])
    return coordinate_full, profile


def _poiseuille_analytic(
    geometry: str, extent_m: float, gradient_pa_per_m: float, viscosity_pa_s: float
) -> tuple[float, float, float]:
    if geometry == "tube":
        max_velocity_m_s = gradient_pa_per_m * extent_m**2 / (4.0 * viscosity_pa_s)
        mean_velocity_m_s = gradient_pa_per_m * extent_m**2 / (8.0 * viscosity_pa_s)
        flow_rate_m3_s = (
            math.pi * extent_m**4 * gradient_pa_per_m / (8.0 * viscosity_pa_s)
        )
    else:
        max_velocity_m_s = gradient_pa_per_m * extent_m**2 / (2.0 * viscosity_pa_s)
        mean_velocity_m_s = gradient_pa_per_m * extent_m**2 / (3.0 * viscosity_pa_s)
        flow_rate_m3_s = 2.0 / 3.0 * extent_m**3 * gradient_pa_per_m / viscosity_pa_s
    return flow_rate_m3_s, mean_velocity_m_s, max_velocity_m_s


def solve_poiseuille(
    length_m: float,
    viscosity_pa_s: float,
    pressure_drop_pa: float,
    radius_m: float | None = None,
    half_height_m: float | None = None,
    width_m: float | None = None,
    bins: int = 512,
) -> dict:
    """Solve pressure-driven viscous flow in a tube or a plane channel.

    Give ``radius_m`` for a circular tube, or ``half_height_m`` for a plane
    channel (``width_m`` defaults to 1 m).  Units: ``length_m`` [m],
    ``viscosity_pa_s`` [Pa s], ``pressure_drop_pa`` [Pa], radii [m].

    The profile is solved by finite differences on the cross-section.  The
    result is ``cross-checked`` against the exact Hagen-Poiseuille or
    plane-Poiseuille solution.
    """
    length_m = _require_positive(length_m, "length_m")
    viscosity_pa_s = _require_positive(viscosity_pa_s, "viscosity_pa_s")
    pressure_drop_pa = _require_positive(pressure_drop_pa, "pressure_drop_pa")
    if bins < 3:
        raise ValueError("bins must be at least 3")

    if radius_m is not None:
        geometry = "tube"
        extent_m = _require_positive(radius_m, "radius_m")
        area_m2 = math.pi * extent_m**2
        width_m = None
    elif half_height_m is not None:
        geometry = "channel"
        extent_m = _require_positive(half_height_m, "half_height_m")
        width_m = 1.0 if width_m is None else _require_positive(width_m, "width_m")
        area_m2 = 2.0 * extent_m * width_m
    else:
        raise ValueError("give radius_m for a tube or half_height_m for a channel")

    gradient_pa_per_m = pressure_drop_pa / length_m
    coordinate_m, profile_m_s = _poiseuille_profile(
        geometry, extent_m, gradient_pa_per_m, viscosity_pa_s, bins
    )

    if geometry == "tube":
        integrand = profile_m_s * coordinate_m
        flow_rate_m3_s = 2.0 * math.pi * float(np.trapezoid(integrand, coordinate_m))
    else:
        flow_rate_m3_s = width_m * 2.0 * float(np.trapezoid(profile_m_s, coordinate_m))
    mean_velocity_m_s = flow_rate_m3_s / area_m2
    max_velocity_m_s = float(profile_m_s[0])

    analytic_q, analytic_mean, analytic_max = _poiseuille_analytic(
        geometry, extent_m, gradient_pa_per_m, viscosity_pa_s
    )
    if geometry == "channel":
        analytic_q *= width_m

    resistance = pressure_drop_pa / flow_rate_m3_s
    relative_error = _relative_error(flow_rate_m3_s, analytic_q)
    result = {
        "geometry": geometry,
        "length_m": length_m,
        "viscosity_pa_s": viscosity_pa_s,
        "pressure_drop_pa": pressure_drop_pa,
        "extent_m": extent_m,
        "bins": bins,
        "flow_rate_m3_s": flow_rate_m3_s,
        "mean_velocity_m_s": mean_velocity_m_s,
        "max_velocity_m_s": max_velocity_m_s,
        "hydraulic_resistance_pa_s_m3": resistance,
        "analytic_flow_rate_m3_s": analytic_q,
        "analytic_mean_velocity_m_s": analytic_mean,
        "analytic_max_velocity_m_s": analytic_max,
        "relative_error": relative_error,
        "converged": relative_error < 1.0e-4,
        "method": "poiseuille_finite_difference",
        "validation": "cross-checked",
    }
    if geometry == "tube":
        result["radius_m"] = extent_m
    else:
        result["half_height_m"] = extent_m
        result["width_m"] = width_m
    return result


def _element_mass_kg(symbol_or_number: str | int) -> float:
    if isinstance(symbol_or_number, (int, np.integer)):
        symbol = _SYMBOL_BY_Z.get(int(symbol_or_number))
        if symbol is None:
            raise ValueError(f"unknown atomic number {symbol_or_number}")
    else:
        symbol = str(symbol_or_number)
    try:
        return _ATOMIC_MASS_AMU[symbol] * AMU_TO_KG
    except KeyError as error:
        raise ValueError(f"unknown element {symbol!r}") from error


def _atoms_from_plain_mapping(mapping: dict) -> tuple[list[str | int], np.ndarray]:
    if "atoms" in mapping:
        symbols: list[str | int] = []
        coordinates: list[list[float]] = []
        for atom in mapping["atoms"]:
            symbol = atom.get("element", atom.get("symbol"))
            position = atom.get("position_m", atom.get("position"))
            if symbol is None or position is None:
                raise ValueError("each atom needs an element and a position_m")
            symbols.append(symbol)
            coordinates.append([float(value) for value in position])
        positions = np.asarray(coordinates, dtype=float)
    elif "elements" in mapping and "positions_m" in mapping:
        symbols = list(mapping["elements"])
        flat = np.asarray(mapping["positions_m"], dtype=float).reshape(-1)
        positions = flat.reshape(-1, 3)
    else:
        raise ValueError("mapping needs 'atoms', or 'elements' and 'positions_m'")
    if positions.size == 0:
        raise ValueError("no atoms found")
    if positions.shape[1] != 3:
        raise ValueError("positions must have three components")
    if positions.shape[0] != len(symbols):
        raise ValueError("element and position counts differ")
    return symbols, positions


def _atoms_from_topology(topology: object) -> tuple[list[int], np.ndarray]:
    elements = topology.elements()
    positions_m = topology.positions_m()
    symbols = [int(value) for value in elements]
    positions = np.asarray(positions_m, dtype=float).reshape(-1, 3)
    if positions.shape[0] != len(symbols):
        raise ValueError("topology element and position counts differ")
    return symbols, positions


def _coerce_atoms(
    document_or_part: object,
) -> tuple[list[str | int], np.ndarray]:
    if isinstance(document_or_part, dict):
        return _atoms_from_plain_mapping(document_or_part)
    topology = getattr(document_or_part, "topology", None)
    if topology is not None:
        return _atoms_from_topology(topology)
    try:
        from . import _core
    except ImportError:
        _core = None
    if _core is not None and isinstance(document_or_part, _core.Part):
        return _atoms_from_topology(document_or_part.topology)
    raise TypeError("expected a plain atom mapping or a part with a topology")


def voxelize_part(document_or_part: object, voxel_size_m: float | None = None) -> dict:
    """Map an atom list to an effective continuum density and a bounding box.

    ``document_or_part`` is a plain mapping with ``atoms`` (each atom has an
    ``element`` and a ``position_m``), or a mapping with ``elements`` and a
    flat ``positions_m``, or a part object with a ``topology``.  The core
    extension is imported lazily and is not required.

    When ``voxel_size_m`` is ``None`` a size of half the smallest bounding-box
    extent is used.  The effective density is the atom mass divided by the
    volume of the occupied voxels.  This result is ``unverified`` because no
    exact continuum value exists for a discrete atom list.
    """
    symbols, positions = _coerce_atoms(document_or_part)
    masses_kg = np.array([_element_mass_kg(symbol) for symbol in symbols])
    total_mass_kg = float(masses_kg.sum())

    box_min_m = positions.min(axis=0)
    box_max_m = positions.max(axis=0)
    box_size_m = box_max_m - box_min_m
    box_volume_m3 = float(np.prod(box_size_m))
    center_of_mass_m = (positions * masses_kg[:, None]).sum(axis=0) / total_mass_kg

    if voxel_size_m is None:
        positive = box_size_m[box_size_m > 0.0]
        voxel_size_m = float(positive.min() / 2.0) if positive.size else 1.0
    else:
        voxel_size_m = _require_positive(voxel_size_m, "voxel_size_m")

    grid = np.floor((positions - box_min_m) / voxel_size_m).astype(int)
    occupied = int(np.unique(grid, axis=0).shape[0])
    voxel_volume_m3 = voxel_size_m**3
    effective_density_kg_m3 = total_mass_kg / (occupied * voxel_volume_m3)
    box_density_kg_m3 = (
        total_mass_kg / box_volume_m3 if box_volume_m3 > 0.0 else math.inf
    )

    return {
        "atom_count": int(positions.shape[0]),
        "mass_kg": total_mass_kg,
        "center_of_mass_m": center_of_mass_m.tolist(),
        "bounding_box_min_m": box_min_m.tolist(),
        "bounding_box_max_m": box_max_m.tolist(),
        "bounding_box_size_m": box_size_m.tolist(),
        "bounding_box_volume_m3": box_volume_m3,
        "voxel_size_m": voxel_size_m,
        "occupied_voxel_count": occupied,
        "voxel_volume_m3": voxel_volume_m3,
        "effective_density_kg_m3": effective_density_kg_m3,
        "bounding_box_density_kg_m3": box_density_kg_m3,
        "method": "atomic_voxelization",
        "validation": "unverified",
    }


def voxelize_atom_list(
    elements: list[str | int],
    positions_m: list[float],
    voxel_size_m: float | None = None,
) -> dict:
    """Convenience wrapper: voxelize a flat element and position list."""
    return voxelize_part(
        {"elements": elements, "positions_m": positions_m},
        voxel_size_m=voxel_size_m,
    )
