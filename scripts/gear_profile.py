#!/usr/bin/env python3
"""The exact full-depth involute gear profile, in Python.

This module is a line-for-line port of the 2D geometry in
`crates/parts/src/gear_profile.rs` and `reflect_to_internal` in
`crates/parts/src/planetary.rs`. The video renderer and the atom-geometry
check import it, so the schematic outline and the atomic layer use one
formula. When the Rust profile changes, change this file in the same commit.

All lengths are SI metres. The pressure angle is 20 degrees.
"""

from __future__ import annotations

import json
import math

PRESSURE_ANGLE_RAD = math.radians(20.0)

# Full-depth addendum and dedendum coefficients `h_a*` and `h_f*`.
# Source: J. E. Shigley, Mechanical Engineering Design.
ADDENDUM_COEFF = 1.0
DEDENDUM_COEFF = 1.25


def involute_function(alpha_rad: float) -> float:
    return math.tan(alpha_rad) - alpha_rad


def pitch_radius_m(module_m: float, teeth: int) -> float:
    return module_m * teeth / 2.0


def base_radius_m(module_m: float, teeth: int,
                  pressure_angle_rad: float = PRESSURE_ANGLE_RAD) -> float:
    return pitch_radius_m(module_m, teeth) * math.cos(pressure_angle_rad)


def addendum_m(module_m: float) -> float:
    return ADDENDUM_COEFF * module_m


def dedendum_m(module_m: float) -> float:
    return DEDENDUM_COEFF * module_m


def outer_radius_m(module_m: float, teeth: int) -> float:
    return pitch_radius_m(module_m, teeth) + addendum_m(module_m)


def root_radius_m(module_m: float, teeth: int) -> float:
    return pitch_radius_m(module_m, teeth) - dedendum_m(module_m)


def flank_angle_rad(radius_m: float, sign: float, module_m: float, teeth: int,
                    pressure_angle_rad: float = PRESSURE_ANGLE_RAD) -> float:
    base_m = base_radius_m(module_m, teeth, pressure_angle_rad)
    if radius_m < base_m:
        raise ValueError(f"radius {radius_m} m is inside the base circle")
    cosine = max(-1.0, min(1.0, base_m / radius_m))
    alpha_r = math.acos(cosine)
    half_tooth_rad = math.pi / (2.0 * teeth)
    return sign * (half_tooth_rad + involute_function(pressure_angle_rad)
                   - involute_function(alpha_r))


def outline_points(module_m: float, teeth: int, flank_samples: int = 6,
                   arc_samples: int = 3,
                   pressure_angle_rad: float = PRESSURE_ANGLE_RAD,
                   ) -> list[tuple[float, float]]:
    """The closed counter-clockwise outline of one external gear, in metres."""
    flank_samples = max(1, flank_samples)
    arc_samples = max(1, arc_samples)
    teeth_f = float(teeth)
    tooth_pitch_rad = 2.0 * math.pi / teeth_f
    half_pitch_rad = math.pi / teeth_f
    outer_m = outer_radius_m(module_m, teeth)
    root_m = root_radius_m(module_m, teeth)
    base_m = base_radius_m(module_m, teeth, pressure_angle_rad)
    if root_m <= 0.0:
        raise ValueError(f"root radius {root_m} m must be positive")

    inv_pitch = involute_function(pressure_angle_rad)
    tip_ratio = max(-1.0, min(1.0, base_m / outer_m))
    tip_half_angle = half_pitch_rad + inv_pitch - involute_function(
        math.acos(tip_ratio))

    flank_start_m = max(root_m, base_m)
    flank_start_angle = flank_angle_rad(flank_start_m, -1.0, module_m, teeth,
                                        pressure_angle_rad)

    local: list[tuple[float, float]] = []
    if root_m < base_m:
        for i in range(arc_samples + 1):
            t = i / arc_samples
            angle = -half_pitch_rad + (half_pitch_rad + flank_start_angle) * t
            local.append((root_m * math.cos(angle), root_m * math.sin(angle)))
    local.append((flank_start_m * math.cos(flank_start_angle),
                  flank_start_m * math.sin(flank_start_angle)))
    for i in range(1, flank_samples + 1):
        t = i / flank_samples
        radius_m = flank_start_m + (outer_m - flank_start_m) * t
        angle = flank_angle_rad(radius_m, -1.0, module_m, teeth,
                                pressure_angle_rad)
        local.append((radius_m * math.cos(angle), radius_m * math.sin(angle)))
    for i in range(arc_samples + 1):
        t = i / arc_samples
        angle = -tip_half_angle + 2.0 * tip_half_angle * t
        local.append((outer_m * math.cos(angle), outer_m * math.sin(angle)))
    for i in range(flank_samples - 1, -1, -1):
        t = i / flank_samples
        radius_m = flank_start_m + (outer_m - flank_start_m) * t
        angle = flank_angle_rad(radius_m, 1.0, module_m, teeth,
                                pressure_angle_rad)
        local.append((radius_m * math.cos(angle), radius_m * math.sin(angle)))
    if root_m < base_m:
        for i in range(arc_samples + 1):
            t = i / arc_samples
            angle = -flank_start_angle + (half_pitch_rad + flank_start_angle) * t
            local.append((root_m * math.cos(angle), root_m * math.sin(angle)))

    points: list[tuple[float, float]] = []
    for tooth in range(teeth):
        rotation = tooth * tooth_pitch_rad
        cos_rot, sin_rot = math.cos(rotation), math.sin(rotation)
        for x, y in local:
            points.append((x * cos_rot - y * sin_rot, x * sin_rot + y * cos_rot))

    deduped: list[tuple[float, float]] = []
    for point in points:
        if deduped and abs(deduped[-1][0] - point[0]) < 1e-15 \
                and abs(deduped[-1][1] - point[1]) < 1e-15:
            continue
        deduped.append(point)
    return deduped


def reflect_to_internal(points_m: list[tuple[float, float]],
                        pitch_radius_m_value: float) -> list[tuple[float, float]]:
    """Map an external outline to an internal one by reflection at the pitch circle."""
    result = []
    for x, y in points_m:
        radius_m = math.hypot(x, y)
        if radius_m <= 0.0:
            result.append((x, y))
            continue
        scale = (2.0 * pitch_radius_m_value - radius_m) / radius_m
        result.append((x * scale, y * scale))
    return result


def design_from_scene(path: str) -> dict:
    with open(path, encoding="utf-8") as handle:
        document = json.load(handle)
    return document["design"]


def profile_radii(module_m: float, teeth: int) -> dict:
    """The pitch, base, outer (tip), and root radii of an external gear."""
    return {
        "pitch_m": pitch_radius_m(module_m, teeth),
        "base_m": base_radius_m(module_m, teeth),
        "outer_m": outer_radius_m(module_m, teeth),
        "root_m": root_radius_m(module_m, teeth),
    }
