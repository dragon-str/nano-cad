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
const BASIS: [[f64; 3]; 8] = [
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

    let margin_m = a;
    let ix0 = ((min_x - margin_m) / a).floor() as i64;
    let ix1 = ((max_x + margin_m) / a).ceil() as i64;
    let iy0 = ((min_y - margin_m) / a).floor() as i64;
    let iy1 = ((max_y + margin_m) / a).ceil() as i64;
    let iz1 = ((layers as i64) + 3) / 4;

    let mut positions_m = Vec::new();
    for iz in 0..=iz1 {
        for ix in ix0..=ix1 {
            for iy in iy0..=iy1 {
                for basis in BASIS {
                    let plane = (iz as f64 + basis[2]) * 4.0;
                    let plane_index = plane.round() as i64;
                    if plane_index < 0 || plane_index >= layers as i64 {
                        continue;
                    }
                    let x_m = (ix as f64 + basis[0]) * a;
                    let y_m = (iy as f64 + basis[1]) * a;
                    let z_m = plane_index as f64 * quarter_m;
                    let radius_m = (x_m * x_m + y_m * y_m).sqrt();
                    let inside = if internal {
                        let within_rim = outer_radius_m.is_none_or(|rim| radius_m <= rim);
                        within_rim && !point_in_polygon(x_m, y_m, outline_m)
                    } else {
                        point_in_polygon(x_m, y_m, outline_m)
                    };
                    if inside {
                        positions_m.push([x_m, y_m, z_m]);
                    }
                }
            }
        }
    }
    positions_m
}

/// Bonds every first-shell carbon pair in `atoms` and caps the surface.
///
/// The bonds use the diamond first-shell length. A carbon with fewer than
/// four bonds receives one hydrogen for each missing tetrahedral direction,
/// so every carbon reaches a valence of four.
pub(crate) fn bond_and_cap(topology: &mut Topology, atoms: &[u32]) -> Result<Vec<u32>, PartError> {
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

    let mut degree: HashMap<u32, u32> = HashMap::new();
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
                            *degree.entry(u).or_insert(0) += 1;
                            *degree.entry(v).or_insert(0) += 1;
                        }
                    }
                }
            }
        }
    }

    let mut capped: Vec<u32> = Vec::new();
    for &index in atoms {
        if *degree.get(&index).unwrap_or(&0) >= 4 {
            continue;
        }
        let position_m = match topology.position_m(index as usize) {
            Some(position) => position,
            None => continue,
        };
        for direction in tetrahedral_directions(position_m, DIAMOND_LATTICE_CONSTANT_M) {
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
            let candidate_m = [
                position_m[0] + C_H_BOND_M * direction[0],
                position_m[1] + C_H_BOND_M * direction[1],
                position_m[2] + C_H_BOND_M * direction[2],
            ];
            let hydrogen = topology.add_atom(Atom::new(Element::HYDROGEN, candidate_m, 0.0, "H"));
            topology.add_bond(Bond::new(index, hydrogen, 1, BondType::Single))?;
            capped.push(hydrogen);
        }
    }
    Ok(capped)
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
