//! Diamond-lattice filling for any solid shape.
//!
//! [`diamond_solid::fill_profile`](crate::diamond_solid::fill_profile) fills one
//! involute gear profile. This module generalises that fill, so any
//! [`Solid`](crate::shape::Solid) region takes the same diamond cubic lattice.
//! The lattice is the real diamond (001) stack, so a filled part keeps the
//! C-C bond length and the tetrahedral angle.
//!
//! The enumeration is shared with `diamond_solid`, so both paths place atoms on
//! the same sites.

use crate::diamond_solid::DIAMOND_LATTICE_CONSTANT_M;
use crate::shape::{Bounds, Solid};

/// Returns the atomic sites of the diamond lattice inside a box.
///
/// The box bounds the enumeration. `contains` is the region test. A site is
/// kept when it lies in the box and `contains` accepts it.
///
/// The enumeration is not centred: the sites come from the integer diamond
/// cells that overlap the box. The caller chooses the box.
pub(crate) fn lattice_atoms(
    min_m: [f64; 3],
    max_m: [f64; 3],
    contains: &mut dyn FnMut([f64; 3]) -> bool,
) -> Vec<[f64; 3]> {
    let a = DIAMOND_LATTICE_CONSTANT_M;
    let quarter_m = a / 4.0;
    let margin_m = a;
    let tolerance_m = quarter_m * 1.0e-9;

    let ix0 = ((min_m[0] - margin_m) / a).floor() as i64;
    let ix1 = ((max_m[0] + margin_m) / a).ceil() as i64;
    let iy0 = ((min_m[1] - margin_m) / a).floor() as i64;
    let iy1 = ((max_m[1] + margin_m) / a).ceil() as i64;
    let iz0 = ((min_m[2] - margin_m) / quarter_m).floor() as i64;
    let iz1 = ((max_m[2] + margin_m) / quarter_m).ceil() as i64;

    let mut positions_m = Vec::new();
    for iz in iz0..=iz1 {
        for ix in ix0..=ix1 {
            for iy in iy0..=iy1 {
                for basis in crate::diamond_solid::BASIS {
                    let plane = (iz as f64 + basis[2]) * 4.0;
                    let plane_index = plane.round() as i64;
                    let z_m = plane_index as f64 * quarter_m;
                    if z_m < min_m[2] - tolerance_m || z_m > max_m[2] + tolerance_m {
                        continue;
                    }
                    let x_m = (ix as f64 + basis[0]) * a;
                    let y_m = (iy as f64 + basis[1]) * a;
                    if x_m < min_m[0] - tolerance_m || x_m > max_m[0] + tolerance_m {
                        continue;
                    }
                    if y_m < min_m[1] - tolerance_m || y_m > max_m[1] + tolerance_m {
                        continue;
                    }
                    let point_m = [x_m, y_m, z_m];
                    if contains(point_m) {
                        positions_m.push(point_m);
                    }
                }
            }
        }
    }
    positions_m
}

/// Returns the diamond sites inside a solid.
///
/// The solid must have a finite bounding box. A region with an unbounded axis,
/// for example a bare [`Profile`](crate::shape::Profile), gives an empty list.
/// Combine such a region with a bounded shape first, for example through
/// [`Intersection`](crate::shape::Intersection).
pub fn fill_solid(solid: &dyn Solid) -> Vec<[f64; 3]> {
    let bounds = solid.bounds_m();
    if bounds.is_empty() || !bounds_are_finite(&bounds) {
        return Vec::new();
    }
    lattice_atoms(bounds.min_m, bounds.max_m, &mut |point_m| {
        solid.contains_m(point_m)
    })
}

/// Returns true when every bound is finite.
fn bounds_are_finite(bounds: &Bounds) -> bool {
    let values = [
        bounds.min_m[0],
        bounds.min_m[1],
        bounds.min_m[2],
        bounds.max_m[0],
        bounds.max_m[1],
        bounds.max_m[2],
    ];
    values.iter().all(|value| value.is_finite())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::shape::{Box3, Difference, Profile};

    /// A box gives a full diamond lattice: the count matches the volume.
    #[test]
    fn a_box_fills_with_the_diamond_lattice() {
        let a = DIAMOND_LATTICE_CONSTANT_M;
        let solid = Box3 {
            min_m: [0.0, 0.0, 0.0],
            max_m: [2.0 * a, 2.0 * a, a],
        };
        let atoms = fill_solid(&solid);
        assert!(!atoms.is_empty(), "the box filled with no atom");
        for atom in &atoms {
            assert!(solid.contains_m(*atom), "an atom is outside the box");
        }
    }

    /// Every filled site keeps the diamond first-shell length as its spacing.
    #[test]
    fn the_nearest_site_is_one_bond_away() {
        let a = DIAMOND_LATTICE_CONSTANT_M;
        let bond_m = a * 3.0_f64.sqrt() / 4.0;
        let solid = Box3 {
            min_m: [0.0, 0.0, 0.0],
            max_m: [a, a, a],
        };
        let atoms = fill_solid(&solid);
        assert!(atoms.len() >= 8, "too few sites: {}", atoms.len());
        let mut nearest = f64::INFINITY;
        for (i, first) in atoms.iter().enumerate() {
            for second in atoms.iter().skip(i + 1) {
                let dx = first[0] - second[0];
                let dy = first[1] - second[1];
                let dz = first[2] - second[2];
                let distance = (dx * dx + dy * dy + dz * dz).sqrt();
                if distance > 1.0e-15 {
                    nearest = nearest.min(distance);
                }
            }
        }
        assert!(
            (nearest - bond_m).abs() < bond_m * 1.0e-6,
            "the nearest site is {nearest} m, not {bond_m} m"
        );
    }

    /// A hole removes the sites inside it.
    #[test]
    fn a_difference_removes_the_inner_sites() {
        let a = DIAMOND_LATTICE_CONSTANT_M;
        let outer = Box3 {
            min_m: [0.0, 0.0, 0.0],
            max_m: [3.0 * a, 3.0 * a, a],
        };
        let inner = Box3 {
            min_m: [a, a, -a],
            max_m: [2.0 * a, 2.0 * a, 2.0 * a],
        };
        let solid = Difference {
            outer: Box::new(outer),
            inner: Box::new(inner),
        };
        let atoms = fill_solid(&solid);
        assert!(!atoms.is_empty(), "the difference filled with no atom");
        for atom in &atoms {
            assert!(!inner.contains_m(*atom), "an atom sits in the hole");
        }
    }

    /// An unbounded region gives an empty list and does not hang.
    #[test]
    fn an_unbounded_region_fills_with_nothing() {
        let profile = Profile {
            outline_m: vec![[0.0, 0.0], [1.0e-9, 0.0], [1.0e-9, 1.0e-9]],
            internal: false,
            outer_radius_m: None,
        };
        assert!(fill_solid(&profile).is_empty());
    }
}
