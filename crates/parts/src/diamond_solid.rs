//! Solid hydrogen-capped diamond filling for a gear profile.
//!
//! A gear is not a wire outline. It is a solid block of diamond cubic
//! carbon, cut to the involute profile. This module fills the profile with
//! the real diamond lattice, bonds every carbon to its first-shell
//! neighbours, and caps each surface carbon with hydrogen so that every
//! carbon keeps a valence of four.
//!
//! The geometry follows the diamond cubic structure. The conventional
//! lattice constant is `a = 3.567e-10 m`. The first-shell bond length is
//! `a * sqrt(3) / 4 = 1.544e-10 m`, and the tetrahedral bond angle is
//! `acos(-1/3) = 109.4712 deg`.

use std::collections::HashMap;

use nanocad_model::{Atom, Bond, BondType, Element, Topology};

use crate::error::PartError;

/// The diamond conventional cubic lattice constant at 300 K.
///
/// Source: CRC Handbook of Chemistry and Physics; N. W. Ashcroft and
/// N. D. Mermin, Solid State Physics (1976), chapter 4.
pub(crate) const DIAMOND_LATTICE_CONSTANT_M: f64 = 3.567e-10;

/// The first-shell carbon-carbon bond length in diamond.
pub(crate) fn diamond_bond_length_m() -> f64 {
    DIAMOND_LATTICE_CONSTANT_M * (3.0_f64).sqrt() / 4.0
}

/// The carbon-hydrogen bond length of a hydrogen-capped surface.
///
/// Source: the C-H bond length in methane, 1.09e-10 m.
const C_H_BOND_M: f64 = 1.09e-10;

/// The spacing between adjacent diamond (001) atomic planes, `a / 4`.
pub(crate) fn diamond_plane_spacing_m() -> f64 {
    DIAMOND_LATTICE_CONSTANT_M / 4.0
}

/// The eight basis atoms of the diamond cubic cell, in fractional
/// coordinates. Two interpenetrating FCC lattices, offset by (1/4,1/4,1/4).
pub(crate) const BASIS: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0],
    [0.0, 0.5, 0.5],
    [0.5, 0.0, 0.5],
    [0.5, 0.5, 0.0],
    [0.25, 0.25, 0.25],
    [0.25, 0.75, 0.75],
    [0.75, 0.25, 0.75],
    [0.75, 0.75, 0.25],
];

/// The four tetrahedral neighbour directions of a sublattice-A carbon.
/// A sublattice-B carbon uses the negative of each direction.
const TETRAHEDRAL_A: [[f64; 3]; 4] = [
    [1.0, 1.0, 1.0],
    [1.0, -1.0, -1.0],
    [-1.0, 1.0, -1.0],
    [-1.0, -1.0, 1.0],
];

/// Returns the diamond carbon positions for one gear.
///
/// The outline is a closed polygon in the local `xy` plane, in metres. For an
/// external gear the solid is the interior of the outline. For an internal
/// ring the solid is the annulus between the outline and `outer_radius_m`.
///
/// The gear has `layers` atomic (001) planes. Consecutive planes are spaced
/// `a / 4`, and the slab is centred on `z = 0`. Each returned atom sits on a
/// real diamond site.
pub(crate) fn fill_profile(
    outline_m: &[[f64; 2]],
    internal: bool,
    outer_radius_m: Option<f64>,
    layers: usize,
) -> Vec<[f64; 3]> {
    if layers == 0 || outline_m.len() < 3 {
        return Vec::new();
    }
    let a = DIAMOND_LATTICE_CONSTANT_M;
    let quarter_m = a / 4.0;

    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for point in outline_m {
        min_x = min_x.min(point[0]);
        min_y = min_y.min(point[1]);
        max_x = max_x.max(point[0]);
        max_y = max_y.max(point[1]);
    }
    if let Some(radius_m) = outer_radius_m {
        min_x = min_x.min(-radius_m);
        min_y = min_y.min(-radius_m);
        max_x = max_x.max(radius_m);
        max_y = max_y.max(radius_m);
    }

    let min_m = [min_x, min_y, 0.0];
    let max_m = [max_x, max_y, (layers - 1) as f64 * quarter_m];
    crate::lattice_fill::lattice_atoms(min_m, max_m, &mut |point_m| {
        let radius_m = (point_m[0] * point_m[0] + point_m[1] * point_m[1]).sqrt();
        if internal {
            let within_rim = outer_radius_m.is_none_or(|rim| radius_m <= rim);
            within_rim && !point_in_polygon(point_m[0], point_m[1], outline_m)
        } else {
            point_in_polygon(point_m[0], point_m[1], outline_m)
        }
    })
}

/// One free tetrahedral direction of a surface carbon.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct FreeDirection {
    /// The surface carbon.
    pub host: u32,
    /// The slot of the direction, from 0 to 3.
    pub slot: usize,
    /// The unit direction from the carbon outwards.
    pub direction_m: [f64; 3],
}

/// A group at one free direction, keyed by host and slot.
///
/// A free direction outside the plan takes one hydrogen.
pub(crate) type CapPlan = HashMap<(u32, usize), crate::group::FunctionalGroup>;

/// Bonds every first-shell carbon pair in `atoms` and caps the surface.
///
/// The bonds use the diamond first-shell length. A carbon with fewer than
/// four bonds receives one hydrogen for each missing tetrahedral direction,
/// so every carbon reaches a valence of four.
pub(crate) fn bond_and_cap(topology: &mut Topology, atoms: &[u32]) -> Result<Vec<u32>, PartError> {
    bond_atoms(topology, atoms)?;
    let free = free_directions(topology, atoms);
    cap_free(topology, &free, &CapPlan::new())
}

/// Adds one bond for every first-shell carbon pair in `atoms`.
///
/// The pass adds no cap. Call [`free_directions`] and then [`cap_free`] to
/// fill the surface, or use [`bond_and_cap`] to do both.
pub(crate) fn bond_atoms(topology: &mut Topology, atoms: &[u32]) -> Result<(), PartError> {
    let bond_m = diamond_bond_length_m();
    let tolerance_m = bond_m * 1.0e-3;
    let cell_m = bond_m;

    let mut grid: HashMap<(i64, i64, i64), Vec<u32>> = HashMap::new();
    for &index in atoms {
        if let Some(position_m) = topology.position_m(index as usize) {
            grid.entry(cell_of(position_m, cell_m))
                .or_default()
                .push(index);
        }
    }

    for &u in atoms {
        let position_u_m = match topology.position_m(u as usize) {
            Some(position) => position,
            None => continue,
        };
        let (cx, cy, cz) = cell_of(position_u_m, cell_m);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let Some(neighbours) = grid.get(&(cx + dx, cy + dy, cz + dz)) else {
                        continue;
                    };
                    for &v in neighbours {
                        if v <= u {
                            continue;
                        }
                        let Some(position_v_m) = topology.position_m(v as usize) else {
                            continue;
                        };
                        if (distance_m(position_u_m, position_v_m) - bond_m).abs() <= tolerance_m {
                            topology.add_bond(Bond::new(u, v, 1, BondType::Single))?;
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

/// Returns every free tetrahedral direction of the atoms in `atoms`.
///
/// A direction is free when the carbon has fewer than four bonds inside
/// `atoms`, and when no atom of `atoms` already sits one bond along that
/// direction.
pub(crate) fn free_directions(topology: &Topology, atoms: &[u32]) -> Vec<FreeDirection> {
    let bond_m = diamond_bond_length_m();
    let tolerance_m = bond_m * 1.0e-3;
    let cell_m = bond_m;
    let members: std::collections::HashSet<u32> = atoms.iter().copied().collect();

    let mut grid: HashMap<(i64, i64, i64), Vec<u32>> = HashMap::new();
    for &index in atoms {
        if let Some(position_m) = topology.position_m(index as usize) {
            grid.entry(cell_of(position_m, cell_m))
                .or_default()
                .push(index);
        }
    }

    // One pass over the bonds gives the degree of every member. A per-atom
    // scan would be quadratic, and a gear has more than a hundred thousand
    // atoms.
    let mut degree: HashMap<u32, u32> = HashMap::new();
    for bond in topology.bonds() {
        if members.contains(&bond.u) && members.contains(&bond.v) {
            *degree.entry(bond.u).or_insert(0) += 1;
            *degree.entry(bond.v).or_insert(0) += 1;
        }
    }

    let mut free = Vec::new();
    for &index in atoms {
        let Some(position_m) = topology.position_m(index as usize) else {
            continue;
        };
        if *degree.get(&index).unwrap_or(&0) >= 4 {
            continue;
        }
        for (slot, direction) in tetrahedral_directions(position_m, DIAMOND_LATTICE_CONSTANT_M)
            .iter()
            .enumerate()
        {
            // An existing carbon neighbour sits at the C-C bond length along
            // this direction, so probe there, not at the C-H length.
            let neighbour_m = [
                position_m[0] + bond_m * direction[0],
                position_m[1] + bond_m * direction[1],
                position_m[2] + bond_m * direction[2],
            ];
            // The probe point can sit exactly on a cell boundary, where
            // floating error moves it into a neighbouring cell. Scan the
            // surrounding cells, as the bonding pass does.
            let (px, py, pz) = cell_of(neighbour_m, cell_m);
            let mut occupied = false;
            'cells: for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let Some(neighbours) = grid.get(&(px + dx, py + dy, pz + dz)) else {
                            continue;
                        };
                        for &other in neighbours {
                            if other == index {
                                continue;
                            }
                            if topology.position_m(other as usize).is_some_and(|position| {
                                distance_m(position, neighbour_m) <= tolerance_m
                            }) {
                                occupied = true;
                                break 'cells;
                            }
                        }
                    }
                }
            }
            if occupied {
                continue;
            }
            free.push(FreeDirection {
                host: index,
                slot,
                direction_m: *direction,
            });
        }
    }
    free
}

/// Fills every free direction with a group from `plan`, or one hydrogen.
///
/// The returned indices are the atoms that fill the directions, in the order
/// of `free`.
pub(crate) fn cap_free(
    topology: &mut Topology,
    free: &[FreeDirection],
    plan: &CapPlan,
) -> Result<Vec<u32>, PartError> {
    let mut added = Vec::new();
    for direction in free {
        let Some(host_position_m) = topology.position_m(direction.host as usize) else {
            continue;
        };
        if let Some(group) = plan.get(&(direction.host, direction.slot)) {
            let group_atoms = crate::group::add_group(
                topology,
                direction.host,
                host_position_m,
                direction.direction_m,
                *group,
            )?;
            added.extend(group_atoms);
            continue;
        }
        let candidate_m = [
            host_position_m[0] + C_H_BOND_M * direction.direction_m[0],
            host_position_m[1] + C_H_BOND_M * direction.direction_m[1],
            host_position_m[2] + C_H_BOND_M * direction.direction_m[2],
        ];
        let hydrogen = topology.add_atom(Atom::new(Element::HYDROGEN, candidate_m, 0.0, "H"));
        topology.add_bond(Bond::new(direction.host, hydrogen, 1, BondType::Single))?;
        added.push(hydrogen);
    }
    Ok(added)
}

/// The four tetrahedral directions of the carbon at `position_m`.
fn tetrahedral_directions(position_m: [f64; 3], a: f64) -> [[f64; 3]; 4] {
    let inv = 1.0 / 3.0_f64.sqrt();
    let fx = position_m[0] / a;
    let fy = position_m[1] / a;
    let fz = position_m[2] / a;
    let sum = fx + fy + fz;
    let residual = sum - sum.floor();
    let is_sublattice_a = residual < 0.125 || residual > 0.875;
    let mut directions = [[0.0; 3]; 4];
    for (slot, direction) in directions.iter_mut().zip(TETRAHEDRAL_A) {
        let sign = if is_sublattice_a { 1.0 } else { -1.0 };
        *slot = [
            sign * direction[0] * inv,
            sign * direction[1] * inv,
            sign * direction[2] * inv,
        ];
    }
    directions
}

fn point_in_polygon(x: f64, y: f64, polygon: &[[f64; 2]]) -> bool {
    let mut inside = false;
    let count = polygon.len();
    let mut j = count - 1;
    for i in 0..count {
        let (xi, yi) = (polygon[i][0], polygon[i][1]);
        let (xj, yj) = (polygon[j][0], polygon[j][1]);
        let crosses = (yi > y) != (yj > y) && x < (xj - xi) * (y - yi) / (yj - yi) + xi;
        if crosses {
            inside = !inside;
        }
        j = i;
    }
    inside
}

fn cell_of(position_m: [f64; 3], cell_m: f64) -> (i64, i64, i64) {
    (
        (position_m[0] / cell_m).floor() as i64,
        (position_m[1] / cell_m).floor() as i64,
        (position_m[2] / cell_m).floor() as i64,
    )
}

fn distance_m(a_m: [f64; 3], b_m: [f64; 3]) -> f64 {
    let dx = a_m[0] - b_m[0];
    let dy = a_m[1] - b_m[1];
    let dz = a_m[2] - b_m[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}
