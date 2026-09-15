use nanocad_model::{Bond, BondType, Topology};

use crate::error::PartError;

/// The relative tolerance for a first-shell bond. Distances within this
/// fraction of the expected bond length become bonds.
pub(crate) const BOND_TOLERANCE_RELATIVE: f64 = 1.0e-3;

/// Returns the Euclidean distance between two points, in metres.
pub(crate) fn distance_m(a_m: [f64; 3], b_m: [f64; 3]) -> f64 {
    let dx = a_m[0] - b_m[0];
    let dy = a_m[1] - b_m[1];
    let dz = a_m[2] - b_m[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Returns the Euclidean distance between two planar points, in metres.
pub(crate) fn distance_2d_m(a_m: [f64; 2], b_m: [f64; 2]) -> f64 {
    let dx = a_m[0] - b_m[0];
    let dy = a_m[1] - b_m[1];
    (dx * dx + dy * dy).sqrt()
}

/// Adds one single bond for every atom pair at the first-shell distance.
///
/// The search visits each unordered pair once, so a bond cannot duplicate.
pub(crate) fn first_shell_bonds(topology: &mut Topology, expected_m: f64) -> Result<(), PartError> {
    let low_m = expected_m * (1.0 - BOND_TOLERANCE_RELATIVE);
    let high_m = expected_m * (1.0 + BOND_TOLERANCE_RELATIVE);
    let atom_count = topology.atom_count();
    for u in 0..atom_count {
        for v in (u + 1)..atom_count {
            let position_u_m = match topology.position_m(u) {
                Some(position) => position,
                None => continue,
            };
            let position_v_m = match topology.position_m(v) {
                Some(position) => position,
                None => continue,
            };
            let length_m = distance_m(position_u_m, position_v_m);
            if length_m >= low_m && length_m <= high_m {
                topology.add_bond(Bond::new(u as u32, v as u32, 1, BondType::Single))?;
            }
        }
    }
    Ok(())
}

/// Resamples a closed polyline at a uniform arc-length spacing.
///
/// The input is treated as a closed loop. The first point is emitted, then one
/// point for every `spacing_m` of arc length. The loop start is not repeated.
pub(crate) fn resample_closed_polyline(points_m: &[[f64; 2]], spacing_m: f64) -> Vec<[f64; 2]> {
    let mut sampled = Vec::new();
    if points_m.len() < 2 || spacing_m <= 0.0 {
        return sampled;
    }
    if let Some(first) = points_m.first().copied() {
        sampled.push(first);
    }
    let mut carried_m = 0.0;
    for window in points_m.windows(2) {
        let a_m = window[0];
        let b_m = window[1];
        let segment_m = distance_2d_m(a_m, b_m);
        if segment_m <= 0.0 {
            continue;
        }
        let mut cursor_m = 0.0;
        let mut remaining_m = segment_m;
        while carried_m + remaining_m >= spacing_m {
            let step_m = spacing_m - carried_m;
            cursor_m += step_m;
            remaining_m -= step_m;
            let t = cursor_m / segment_m;
            sampled.push([
                a_m[0] + (b_m[0] - a_m[0]) * t,
                a_m[1] + (b_m[1] - a_m[1]) * t,
            ]);
            carried_m = 0.0;
        }
        carried_m += remaining_m;
    }
    sampled
}
