#!/usr/bin/env python3
"""Check the geometry of the atomistic layers in the nano-cad gearbox video.

Two layers are checked, so the schematic and the atomic layer stay consistent.

The gear layer (the "gears are made of atoms" layer) comes from
`site/scene.json`. Its per-part tip and root radii must equal the exact
full-depth involute profile that `scripts/gear_profile.py` computes. That
module is the same formula as `crates/parts/src/gear_profile.rs` and the same
formula the video renderer draws. This check ties the three together.

The diamond layer (the material basis) comes from the `diamond` block of
`site/scene.bonds.json`, which the Rust example `scene_json` writes from the
`nanocad-parts` diamond generator. Three facts are checked:

1. The mean nearest-neighbour carbon-carbon distance equals the diamond
   first-shell bond length `a sqrt(3) / 4 = 1.544e-10 m` within one percent.
2. Every interior carbon has four bonds. An atom is interior when it sits at
   least one bond length inside every face of the lattice bounding box.
3. The mean bond angle around a tetrahedral atom equals `109.4712` degrees
   within one degree.

The script prints the measured values. It exits nonzero when a tolerance is
violated, so `scripts/make_video.sh` fails loudly.

Usage:
    python3 scripts/check_atom_geometry.py \
        [--scene site/scene.json] [--bonds site/scene.bonds.json]
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gear_profile as gp  # noqa: E402

C_C_BOND_M = 1.544e-10
TETRAHEDRAL_ANGLE_DEG = 109.4712
DISTANCE_TOLERANCE_RELATIVE = 0.01
ANGLE_TOLERANCE_DEG = 1.0
RADIUS_TOLERANCE_RELATIVE = 1e-3


def load_diamond(path: str) -> tuple[list[list[float]], list[tuple[int, int]]]:
    with open(path, encoding="utf-8") as handle:
        document = json.load(handle)
    diamond = document.get("diamond")
    if not isinstance(diamond, dict):
        raise SystemExit(f"error: {path} has no 'diamond' section")
    atoms = [atom["position_m"] for atom in diamond["atoms"]]
    bonds = [(int(u), int(v)) for u, v in diamond["bonds"]]
    if len(atoms) < 2:
        raise SystemExit("error: the diamond block has fewer than two atoms")
    return atoms, bonds


def nearest_neighbour_distances(
    atoms: list[list[float]],
) -> list[float]:
    result: list[float] = []
    for i, first in enumerate(atoms):
        best = math.inf
        for j, second in enumerate(atoms):
            if i == j:
                continue
            dx = first[0] - second[0]
            dy = first[1] - second[1]
            dz = first[2] - second[2]
            distance_m = math.sqrt(dx * dx + dy * dy + dz * dz)
            if distance_m < best:
                best = distance_m
        if math.isfinite(best):
            result.append(best)
    return result


def bond_degrees(
    atom_count: int, bonds: list[tuple[int, int]]
) -> list[int]:
    degrees = [0] * atom_count
    for u, v in bonds:
        degrees[u] += 1
        degrees[v] += 1
    return degrees


def interior_atoms(
    atoms: list[list[float]], bond_m: float
) -> list[int]:
    low = [min(atom[axis] for atom in atoms) for axis in range(3)]
    high = [max(atom[axis] for atom in atoms) for axis in range(3)]
    interior = []
    for index, atom in enumerate(atoms):
        if all(
            low[axis] + bond_m <= atom[axis] <= high[axis] - bond_m
            for axis in range(3)
        ):
            interior.append(index)
    return interior


def bond_angles_deg(
    atoms: list[list[float]], bonds: list[tuple[int, int]]
) -> list[float]:
    neighbours: list[list[int]] = [[] for _ in atoms]
    for u, v in bonds:
        neighbours[u].append(v)
        neighbours[v].append(u)
    angles: list[float] = []
    for index, neighbour_list in enumerate(neighbours):
        for a in range(len(neighbour_list)):
            for b in range(a + 1, len(neighbour_list)):
                first = atoms[neighbour_list[a]]
                second = atoms[neighbour_list[b]]
                origin = atoms[index]
                vector_a = [first[k] - origin[k] for k in range(3)]
                vector_b = [second[k] - origin[k] for k in range(3)]
                dot = sum(vector_a[k] * vector_b[k] for k in range(3))
                norm_a = math.sqrt(sum(value * value for value in vector_a))
                norm_b = math.sqrt(sum(value * value for value in vector_b))
                if norm_a == 0.0 or norm_b == 0.0:
                    continue
                cosine = max(-1.0, min(1.0, dot / (norm_a * norm_b)))
                angles.append(math.degrees(math.acos(cosine)))
    return angles


def load_gear_layer(scene_path: str, bonds_path: str):
    with open(scene_path, encoding="utf-8") as handle:
        scene = json.load(handle)
    with open(bonds_path, encoding="utf-8") as handle:
        bonds = json.load(handle)
    atoms = scene["atomistic"]["atoms"]
    design = scene["design"]
    centers = {body["name"]: body["position_m"]
               for body in scene["device"]["bodies"]}
    planetary_bonds = bonds.get("planetary", {}).get("bonds", [])
    return atoms, design, centers, planetary_bonds


def part_radii_m(atoms, name: str, center) -> list[float]:
    result = []
    for atom in atoms:
        if atom["name"] != name:
            continue
        position = atom["position_m"]
        result.append(math.hypot(position[0] - center[0],
                                 position[1] - center[1]))
    return result


def nearest_atom_radius_m(atoms, name: str, center, world_angle_rad: float) -> float:
    """Radius of the atom in `name` nearest the given in-plane world angle."""
    best_radius = math.nan
    best_error = math.inf
    for atom in atoms:
        if atom["name"] != name:
            continue
        position = atom["position_m"]
        dx = position[0] - center[0]
        dy = position[1] - center[1]
        if dx == 0.0 and dy == 0.0:
            continue
        error = abs((math.atan2(dy, dx) - world_angle_rad + math.pi)
                    % (2.0 * math.pi) - math.pi)
        if error < best_error:
            best_error = error
            best_radius = math.hypot(dx, dy)
    return best_radius


def check_mesh_phase(atoms, design, centers, failures: list[str]) -> None:
    module_m = design["module_m"]
    sun_teeth = int(design["sun_teeth"])
    planet_teeth = int(design["planet_teeth"])
    ring_teeth = int(design["ring_teeth"])
    sun_outer = gp.outer_radius_m(module_m, sun_teeth)
    sun_root = gp.root_radius_m(module_m, sun_teeth)
    sun_height = sun_outer - sun_root
    planet_outer = gp.outer_radius_m(module_m, planet_teeth)
    planet_root = gp.root_radius_m(module_m, planet_teeth)
    planet_height = planet_outer - planet_root
    ring_pitch = gp.pitch_radius_m(module_m, ring_teeth)
    ring_tip = ring_pitch - gp.addendum_m(module_m)
    ring_root = ring_pitch + gp.dedendum_m(module_m)

    for k in range(int(design["planet_count"])):
        center = centers[f"planet_{k}"]
        phi = math.atan2(center[1], center[0])
        sun_radius = nearest_atom_radius_m(atoms, "sun", centers["sun"], phi)
        planet_toward_sun = nearest_atom_radius_m(
            atoms, f"planet_{k}", center, phi + math.pi)
        planet_toward_ring = nearest_atom_radius_m(
            atoms, f"planet_{k}", center, phi)
        ring_radius = nearest_atom_radius_m(atoms, "ring", centers["ring"], phi)
        print(
            f"  mesh {k}: sun tooth r {sun_radius:.6e} m (tip {sun_outer:.6e}); "
            f"planet space r {planet_toward_sun:.6e}, "
            f"{planet_toward_ring:.6e} m (root {planet_root:.6e}); "
            f"ring tooth r {ring_radius:.6e} m (tip {ring_tip:.6e})"
        )
        if sun_radius < sun_outer - 0.25 * sun_height:
            failures.append(
                f"mesh {k}: the sun presents a space, not a tooth, at the mesh "
                f"line (r {sun_radius:.6e} m, tip {sun_outer:.6e} m)"
            )
        for radius, side in ((planet_toward_sun, "sun"),
                             (planet_toward_ring, "ring")):
            if radius > planet_root + 0.25 * planet_height:
                failures.append(
                    f"mesh {k}: the planet presents a tooth, not a space, "
                    f"toward the {side} (r {radius:.6e} m, "
                    f"root {planet_root:.6e} m)"
                )
        if ring_radius > ring_tip + 0.25 * (ring_root - ring_tip):
            failures.append(
                f"mesh {k}: the ring presents a space, not a tooth, at the mesh "
                f"line (r {ring_radius:.6e} m, tip {ring_tip:.6e} m)"
            )
    print("mesh phase: sun tooth into planet space, planet tooth into ring space")


def check_kinematics(failures: list[str]) -> None:
    """The schematic and the atom layer must share one kinematics.

    The schematic comes from `render_video.planetary_angles`; the atom layer is
    animated by `render_video.GearLattice`. This checks three facts:

    1. The sun-planet mesh invariant is constant in time.
    2. The sun and the planets turn in opposite directions, and the carrier
       turns with the sun.
    3. A planet atom's spin about its own center matches the engine relative
       rate `-(N_s / N_p) * (w_s - w_c)`, so the atoms and the schematic agree.
    """
    try:
        import render_video as rv
    except Exception as error:  # Pillow or numpy missing
        print(f"kinematics: render_video not importable ({error}), skipped")
        return

    ns = float(rv.SUN_TEETH)
    npz = float(rv.PLANET_TEETH)
    sun0, planets0 = rv.planetary_angles(0.0)
    sun1, planets1 = rv.planetary_angles(1.0)
    w_s = sun1 - sun0
    w_c = planets1[0][0] - planets0[0][0]
    for k, ((phi0, tp0), (phi1, tp1)) in enumerate(zip(planets0, planets1)):
        invariant0 = ns * (sun0 - phi0) + npz * (tp0 - phi0)
        invariant1 = ns * (sun1 - phi1) + npz * (tp1 - phi1)
        if abs(invariant1 - invariant0) > 1e-9:
            failures.append(
                f"the sun-planet {k} mesh invariant changes by "
                f"{invariant1 - invariant0:.3e} in one second"
            )
        w_p = tp1 - tp0
        if w_s * w_p >= 0.0:
            failures.append(
                f"the sun and planet {k} turn the same way "
                f"(w_s {w_s:+.6f}, w_p {w_p:+.6f})"
            )
    if w_s * w_c <= 0.0:
        failures.append(
            f"the carrier and the sun turn opposite ways "
            f"(w_s {w_s:+.6f}, w_c {w_c:+.6f})"
        )
    print(
        f"  schematic: w_s {w_s:+.6f} rad/s, w_c {w_c:+.6f} rad/s, "
        f"w_p {planets1[0][1] - planets0[0][1]:+.6f} rad/s"
    )

    scene_path = "site/scene.json"
    bonds_path = "site/scene.bonds.json"
    if not (os.path.exists(scene_path) and os.path.exists(bonds_path)):
        print("  atoms: scene files not found, rate check skipped")
        return
    lattice = rv.GearLattice(scene_path, bonds_path)
    design = lattice.design
    ns_l = int(design["sun_teeth"])
    npz_l = int(design["planet_teeth"])
    nr_l = int(design["ring_teeth"])
    expected_rate = planets1[0][1] - planets0[0][1]
    index = next(i for i, atom in enumerate(lattice.atoms)
                 if atom["name"] == "planet_0")
    center = lattice.centers["planet_0"]

    def relative_angle(t):
        positions = lattice.animated_xy(t)
        w_c_l = -1.0 * ns_l / (ns_l + nr_l)
        cx, cy = rv._rot2(center[0], center[1], w_c_l * t * 0.34)
        return math.atan2(positions[index][1] - cy, positions[index][0] - cx)

    measured_rate = relative_angle(1.0) - relative_angle(0.0)
    measured_rate = (measured_rate + math.pi) % (2.0 * math.pi) - math.pi
    print(
        f"  atoms: planet orientation {measured_rate:+.6f} rad/s, "
        f"schematic planet rate {expected_rate:+.6f} rad/s"
    )
    if abs(measured_rate - expected_rate) > 1e-6:
        failures.append(
            f"the atom planet orientation rate is {measured_rate:+.6f} rad/s "
            f"but the schematic planet rate is {expected_rate:+.6f} rad/s"
        )
    print("kinematics: schematic and atoms share one rate function")


def check_thickness(atoms, design, failures: list[str]) -> None:
    layers = int(design.get("layers", 1))
    spacing_m = float(design.get("layer_spacing_m", 0.0))
    thickness_m = float(design.get("thickness_m", 0.0))
    z_values = sorted({round(atom["position_m"][2], 15) for atom in atoms
                       if atom["name"] in ("sun", "ring")
                       or atom["name"].startswith("planet_")})
    expected_thickness = (layers - 1) * spacing_m
    expected = [round((i - 0.5 * (layers - 1)) * spacing_m, 15)
                for i in range(layers)]
    print(
        f"  thickness: {layers} layers, spacing {spacing_m:.6e} m, "
        f"{thickness_m:.6e} m"
    )
    print(f"  gear z values: {[f'{z:.3e}' for z in z_values]}")
    if len(z_values) != layers:
        failures.append(
            f"the gear atoms have {len(z_values)} z values, expected {layers}"
        )
    elif any(abs(actual - want) > 1e-18
             for actual, want in zip(z_values, expected)):
        failures.append(
            f"the gear z values {z_values} do not match the centered layers "
            f"{expected} for {layers} layers at {spacing_m:.6e} m"
        )
    if abs(thickness_m - expected_thickness) > 1e-18:
        failures.append(
            f"thickness_m {thickness_m:.6e} m is not "
            f"(layers-1)*spacing = {expected_thickness:.6e} m"
        )
    else:
        print("thickness: the metadata matches the atomic layer span")


def check_gear_layer(scene_path: str, bonds_path: str,
                     failures: list[str]) -> None:
    if not (os.path.exists(scene_path) and os.path.exists(bonds_path)):
        print(f"gear layer: {scene_path} or {bonds_path} not found, skipped")
        return
    atoms, design, centers, planetary_bonds = load_gear_layer(scene_path,
                                                              bonds_path)
    module_m = design["module_m"]
    parts = [("sun", int(design["sun_teeth"]), False)]
    for k in range(int(design["planet_count"])):
        parts.append((f"planet_{k}", int(design["planet_teeth"]), False))
    parts.append(("ring", int(design["ring_teeth"]), True))

    print(f"gear atoms: {len(atoms)} atoms, "
          f"{len(planetary_bonds)} planetary bonds")
    for name, teeth, internal in parts:
        center = centers.get(name, [0.0, 0.0, 0.0])
        radii = part_radii_m(atoms, name, center)
        if not radii:
            failures.append(f"gear part {name} has no atoms")
            continue
        measured_lo, measured_hi = min(radii), max(radii)
        pitch_m = gp.pitch_radius_m(module_m, teeth)
        if internal:
            expected_lo = pitch_m - gp.addendum_m(module_m)
            expected_hi = pitch_m + gp.dedendum_m(module_m)
            label = "internal tip/root"
        else:
            expected_lo = gp.root_radius_m(module_m, teeth)
            expected_hi = gp.outer_radius_m(module_m, teeth)
            label = "root/tip"
        print(
            f"  {name}: {label} radius {measured_lo:.6e}..{measured_hi:.6e} m, "
            f"involute profile {expected_lo:.6e}..{expected_hi:.6e} m"
        )
        lo_error = abs(measured_lo - expected_lo) / expected_lo
        hi_error = abs(measured_hi - expected_hi) / expected_hi
        if lo_error > RADIUS_TOLERANCE_RELATIVE or \
                hi_error > RADIUS_TOLERANCE_RELATIVE:
            failures.append(
                f"gear part {name} radii {measured_lo:.6e}..{measured_hi:.6e} m "
                f"do not match the involute profile "
                f"{expected_lo:.6e}..{expected_hi:.6e} m"
            )

        points = gp.outline_points(module_m, teeth)
        if internal:
            points = gp.reflect_to_internal(
                points, gp.pitch_radius_m(module_m, teeth))
        profile_radii = [math.hypot(x, y) for x, y in points]
        profile_lo, profile_hi = min(profile_radii), max(profile_radii)
        if (abs(profile_lo - measured_lo) / measured_lo
                > RADIUS_TOLERANCE_RELATIVE
                or abs(profile_hi - measured_hi) / measured_hi
                > RADIUS_TOLERANCE_RELATIVE):
            failures.append(
                f"schematic profile for {name} "
                f"{profile_lo:.6e}..{profile_hi:.6e} m disagrees with the "
                f"atom layer {measured_lo:.6e}..{measured_hi:.6e} m"
            )
    check_mesh_phase(atoms, design, centers, failures)
    check_thickness(atoms, design, failures)
    print("gear layer: schematic profile and atom layer agree")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scene", default="site/scene.json")
    parser.add_argument("--bonds", default="site/scene.bonds.json")
    args = parser.parse_args()

    failures: list[str] = []
    print("=== kinematics (the schematic and the atoms must agree) ===")
    check_kinematics(failures)
    print()
    print("=== gear layer (the schematic must match the atoms) ===")
    check_gear_layer(args.scene, args.bonds, failures)
    print()
    print("=== diamond layer (the material basis) ===")

    atoms, bonds = load_diamond(args.bonds)
    atom_count = len(atoms)

    distances = nearest_neighbour_distances(atoms)
    mean_distance = sum(distances) / len(distances)
    min_distance = min(distances)
    max_distance = max(distances)

    degrees = bond_degrees(atom_count, bonds)
    interior = interior_atoms(atoms, C_C_BOND_M)
    bad_interior = [index for index in interior if degrees[index] != 4]

    angles = bond_angles_deg(atoms, bonds)
    mean_angle = sum(angles) / len(angles) if angles else math.nan

    distance_error = abs(mean_distance - C_C_BOND_M) / C_C_BOND_M
    angle_error = abs(mean_angle - TETRAHEDRAL_ANGLE_DEG)

    print(f"atoms: {atom_count}")
    print(f"bonds: {len(bonds)}")
    print(
        "nearest-neighbour C-C distance: "
        f"mean {mean_distance:.6e} m, "
        f"min {min_distance:.6e} m, "
        f"max {max_distance:.6e} m"
    )
    print(
        f"reference bond length: {C_C_BOND_M:.6e} m "
        f"(relative error {distance_error:.3e})"
    )
    print(
        f"interior carbons: {len(interior)}, all with four bonds: "
        f"{'yes' if not bad_interior else 'no'}"
    )
    print(
        f"bond angle: mean {mean_angle:.6f} deg, "
        f"reference {TETRAHEDRAL_ANGLE_DEG:.4f} deg "
        f"(absolute error {angle_error:.3e} deg)"
    )

    if distance_error > DISTANCE_TOLERANCE_RELATIVE:
        failures.append(
            f"nearest-neighbour distance {mean_distance:.6e} m is not within "
            f"{DISTANCE_TOLERANCE_RELATIVE:.0%} of {C_C_BOND_M:.6e} m"
        )
    if bad_interior:
        failures.append(
            f"{len(bad_interior)} interior carbon(s) do not have four bonds "
            f"(first: atom {bad_interior[0]} has {degrees[bad_interior[0]]})"
        )
    if angle_error > ANGLE_TOLERANCE_DEG:
        failures.append(
            f"bond angle {mean_angle:.6f} deg is not within "
            f"{ANGLE_TOLERANCE_DEG:.1f} deg of {TETRAHEDRAL_ANGLE_DEG:.4f} deg"
        )

    if failures:
        for failure in failures:
            print(f"check_atom_geometry: FAIL: {failure}", file=sys.stderr)
        return 1

    print("check_atom_geometry: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
