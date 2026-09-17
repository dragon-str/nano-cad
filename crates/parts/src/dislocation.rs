//! Elastic angular displacement map of a wedge disclination.
//!
//! A wedge disclination in a 2D hexagonal lattice removes or inserts a
//! 60-degree wedge. Removing one 60-degree wedge turns the two shared 6-rings
//! into one 5-ring and one 7-ring: the classic 5-7 defect. This module applies
//! the elastic angular displacement map. All quantities are SI: metres and
//! radians.

use std::f64::consts::TAU;

/// A wedge disclination about the z axis through `origin_m`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WedgeDisclination {
    /// The point the map turns about, in metres.
    pub origin_m: [f64; 3],
    /// The removed wedge angle, in radians. Use `five_seven_wedge_rad`.
    pub wedge_rad: f64,
    /// Atoms nearer than this to the axis are the core and are not mapped.
    pub core_radius_m: f64,
}

impl Default for WedgeDisclination {
    fn default() -> Self {
        Self {
            origin_m: [0.0; 3],
            wedge_rad: five_seven_wedge_rad(),
            core_radius_m: 0.0,
        }
    }
}

/// The removed wedge angle that gives a 5-7 pair, in radians. It is TAU/6.
pub fn five_seven_wedge_rad() -> f64 {
    TAU / 6.0
}

/// Returns the mapped position of one point, in metres.
///
/// The map preserves the radius in the xy plane and the z coordinate. It
/// stretches the angle by `1 / (1 - strength)`, where
/// `strength = wedge_rad / TAU`.
///
/// The map is singular when `strength >= 1.0` or
/// `strength <= -1.0 + 1e-12`. The function returns the point unchanged in
/// that case, so it never panics and never returns a NaN.
///
/// A point at the origin has no defined angle. A point nearer to the axis than
/// `core_radius_m` is core material. Both return the point unchanged. The test
/// is `r < core_radius_m`, so a point at exactly `core_radius_m` is mapped.
pub fn displace_point_m(position_m: [f64; 3], disclination: &WedgeDisclination) -> [f64; 3] {
    let dx_m = position_m[0] - disclination.origin_m[0];
    let dy_m = position_m[1] - disclination.origin_m[1];
    let r_m = dx_m.hypot(dy_m);

    if r_m == 0.0 || r_m < disclination.core_radius_m {
        return position_m;
    }

    let strength = disclination.wedge_rad / TAU;
    if strength >= 1.0 || strength <= -1.0 + 1e-12 {
        return position_m;
    }

    let theta_rad = dy_m.atan2(dx_m);
    let theta_new_rad = theta_rad / (1.0 - strength);

    [
        disclination.origin_m[0] + r_m * theta_new_rad.cos(),
        disclination.origin_m[1] + r_m * theta_new_rad.sin(),
        position_m[2],
    ]
}

/// Maps every position in place.
pub fn displace(positions_m: &mut [[f64; 3]], disclination: &WedgeDisclination) {
    for position_m in positions_m.iter_mut() {
        *position_m = displace_point_m(*position_m, disclination);
    }
}

/// Marks the atoms inside the core radius. The caller removes them.
///
/// The value at each index is `true` when the distance from `origin_m` in the
/// xy plane is less than `core_radius_m`.
pub fn core_atom_mask(positions_m: &[[f64; 3]], disclination: &WedgeDisclination) -> Vec<bool> {
    positions_m
        .iter()
        .map(|position_m| {
            let dx_m = position_m[0] - disclination.origin_m[0];
            let dy_m = position_m[1] - disclination.origin_m[1];
            dx_m.hypot(dy_m) < disclination.core_radius_m
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_wedge_is_sixty_degrees() {
        assert_eq!(five_seven_wedge_rad(), TAU / 6.0);
    }

    #[test]
    fn a_point_on_the_axis_does_not_move() {
        let disclination = WedgeDisclination {
            origin_m: [0.0; 3],
            wedge_rad: five_seven_wedge_rad(),
            core_radius_m: 1e-10,
        };
        let position_m = [0.0, 0.0, 0.0];
        assert_eq!(displace_point_m(position_m, &disclination), position_m);
    }

    #[test]
    fn a_point_inside_the_core_does_not_move() {
        let disclination = WedgeDisclination {
            origin_m: [0.0; 3],
            wedge_rad: five_seven_wedge_rad(),
            core_radius_m: 1e-10,
        };
        let position_m = [5e-11, 0.0, 0.0];
        assert_eq!(displace_point_m(position_m, &disclination), position_m);
    }

    #[test]
    fn a_mapped_point_keeps_its_radius_and_z() {
        let disclination = WedgeDisclination {
            origin_m: [0.0; 3],
            wedge_rad: five_seven_wedge_rad(),
            core_radius_m: 1e-10,
        };
        let mapped_m = displace_point_m([1e-9, 0.0, 2e-9], &disclination);
        let r_m = mapped_m[0].hypot(mapped_m[1]);
        assert!((r_m - 1e-9).abs() < 1e-12);
        assert!((mapped_m[2] - 2e-9).abs() < 1e-12);
    }

    #[test]
    fn the_angle_stretches_by_one_over_one_minus_strength() {
        let disclination = WedgeDisclination {
            origin_m: [0.0; 3],
            wedge_rad: five_seven_wedge_rad(),
            core_radius_m: 1e-10,
        };
        let position_m = [1e-9 * 1.0_f64.cos(), 1e-9 * 1.0_f64.sin(), 0.0];
        let mapped_m = displace_point_m(position_m, &disclination);
        let theta_new_rad = mapped_m[1].atan2(mapped_m[0]);
        let expected_rad = 1.0 / (1.0 - 1.0 / 6.0);
        assert!((theta_new_rad - expected_rad).abs() < 1e-9);
    }

    #[test]
    fn a_singular_wedge_is_refused() {
        let disclination = WedgeDisclination {
            origin_m: [0.0; 3],
            wedge_rad: TAU,
            core_radius_m: 1e-10,
        };
        let position_m = [3e-9, 4e-9, 5e-9];
        assert_eq!(displace_point_m(position_m, &disclination), position_m);
    }

    #[test]
    fn the_core_mask_marks_only_the_core() {
        let disclination = WedgeDisclination {
            origin_m: [0.0; 3],
            wedge_rad: five_seven_wedge_rad(),
            core_radius_m: 1e-10,
        };
        let positions_m = [[5e-11, 0.0, 0.0], [2e-10, 0.0, 0.0]];
        assert_eq!(
            core_atom_mask(&positions_m, &disclination),
            vec![true, false]
        );
    }
}
