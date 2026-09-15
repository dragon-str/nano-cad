#!/usr/bin/env python3
"""Check the geometry of the atomistic layer in the nano-cad gearbox video.

The video's atomistic layer must be real diamondoid carbon. This script reads
the companion bond document that the Rust example `scene_json` writes, and it
checks three facts against the diamond lattice block:

1. The mean nearest-neighbour carbon-carbon distance equals the diamond
   first-shell bond length `a sqrt(3) / 4 = 1.544e-10 m` within one percent.
2. Every interior carbon has four bonds. An atom is interior when it sits at
   least one bond length inside every face of the lattice bounding box.
3. The mean bond angle around a tetrahedral atom equals `109.4712` degrees
   within one degree.

The script prints the measured mean, minimum, and maximum distance, and the
measured mean angle. It exits nonzero when a tolerance is violated, so
`scripts/make_video.sh` fails loudly.

Usage:
    python3 scripts/check_atom_geometry.py [--bonds site/scene.bonds.json]
"""

from __future__ import annotations

import argparse
import json
import math
import sys

C_C_BOND_M = 1.544e-10
TETRAHEDRAL_ANGLE_DEG = 109.4712
DISTANCE_TOLERANCE_RELATIVE = 0.01
ANGLE_TOLERANCE_DEG = 1.0


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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bonds", default="site/scene.bonds.json")
    args = parser.parse_args()

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

    failures: list[str] = []
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
