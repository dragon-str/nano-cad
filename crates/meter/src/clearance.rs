//! The clearance metric: the least distance between two bodies as they move.
//!
//! Two meshing bodies turn about different centres. The clearance is the least
//! atom-to-atom distance over one motion period. A value below the non-bonded
//! diamond spacing means the bodies would interpenetrate.

use std::collections::HashMap;
use std::ops::Range;

use nanocad_parts::planetary::{PlanetaryDesign, PlanetarySet};

use crate::{Fidelity, MetricValue};

type Cell = (i64, i64, i64);

/// The rigid motion of one body in the plane.
///
/// The body orbits about `orbit_center_m` and spins about its own centre. The
/// spin centre is the body centre at zero time. The z axis is the rotation
/// axis, so the motion is planar.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BodyMotion {
    /// The body centre at zero time, in metres.
    pub spin_center_m: [f64; 3],
    /// The spin rate about the body centre, in radians per second.
    pub spin_rad_per_s: f64,
    /// The orbit centre, in metres.
    pub orbit_center_m: [f64; 3],
    /// The orbit rate about the orbit centre, in radians per second.
    pub orbit_rad_per_s: f64,
}

impl BodyMotion {
    /// A body that does not move.
    pub const FIXED: BodyMotion = BodyMotion {
        spin_center_m: [0.0, 0.0, 0.0],
        spin_rad_per_s: 0.0,
        orbit_center_m: [0.0, 0.0, 0.0],
        orbit_rad_per_s: 0.0,
    };

    /// Returns the world position of a body point at the given time.
    pub fn pose(&self, point_m: [f64; 3], time_s: f64) -> [f64; 3] {
        let orbit = self.orbit_rad_per_s * time_s;
        let (sin_orbit, cos_orbit) = orbit.sin_cos();
        let arm_x = self.spin_center_m[0] - self.orbit_center_m[0];
        let arm_y = self.spin_center_m[1] - self.orbit_center_m[1];
        let center_x = self.orbit_center_m[0] + arm_x * cos_orbit - arm_y * sin_orbit;
        let center_y = self.orbit_center_m[1] + arm_x * sin_orbit + arm_y * cos_orbit;
        let spin = self.spin_rad_per_s * time_s;
        let (sin_spin, cos_spin) = spin.sin_cos();
        let dx = point_m[0] - self.spin_center_m[0];
        let dy = point_m[1] - self.spin_center_m[1];
        [
            center_x + dx * cos_spin - dy * sin_spin,
            center_y + dx * sin_spin + dy * cos_spin,
            point_m[2],
        ]
    }
}

/// A set of atoms grouped into moving bodies.
///
/// Each group is a contiguous range into `positions_m`. The motions and the
/// groups are parallel vectors: motion `k` moves group `k`.
#[derive(Clone, Debug, PartialEq)]
pub struct MovingAtoms {
    /// The atom positions at zero time, in metres.
    pub positions_m: Vec<[f64; 3]>,
    /// The body groups, each with a name and a contiguous atom range.
    pub groups: Vec<(String, Range<usize>)>,
    /// The motion of each group, parallel to `groups`.
    pub motions: Vec<BodyMotion>,
    /// The period of the motion in seconds. The sweep covers one period.
    pub period_s: f64,
}

impl MovingAtoms {
    /// Builds the moving atoms of a planetary set with the ring fixed.
    ///
    /// The sun spins about the origin. Each planet orbits about the origin and
    /// spins about its own centre. The ring does not move.
    pub fn planetary(set: &PlanetarySet, design: &PlanetaryDesign, sun_rad_per_s: f64) -> Self {
        let topology = &set.part.topology;
        let mut positions_m = Vec::with_capacity(topology.atom_count());
        let mut groups = Vec::new();

        let mut push_group = |name: String, range: Range<usize>, positions: &mut Vec<[f64; 3]>| {
            let start = positions.len();
            for index in range {
                positions.push(topology.position_m(index).unwrap_or([0.0; 3]));
            }
            groups.push((name, start..positions.len()));
        };

        push_group("sun".to_string(), set.sun_atoms.clone(), &mut positions_m);
        for (index, range) in set.planet_atoms.iter().enumerate() {
            push_group(format!("planet_{index}"), range.clone(), &mut positions_m);
        }
        push_group("ring".to_string(), set.ring_atoms.clone(), &mut positions_m);

        let sun_teeth = design.sun_teeth() as f64;
        let ring_teeth = design.ring_teeth() as f64;
        let planet_teeth = design.planet_teeth() as f64;
        let carrier_rate = sun_rad_per_s * sun_teeth / (sun_teeth + ring_teeth);
        // The planet spins on the sun in the carrier frame. Its rate relative
        // to the carrier is below; its absolute rate adds the carrier rate.
        let planet_relative_rate = -(sun_teeth / planet_teeth) * (sun_rad_per_s - carrier_rate);

        let origin = [0.0, 0.0, 0.0];
        let mut motions = Vec::with_capacity(groups.len());
        motions.push(BodyMotion {
            spin_center_m: origin,
            spin_rad_per_s: sun_rad_per_s,
            orbit_center_m: origin,
            orbit_rad_per_s: 0.0,
        });
        for index in 0..set.planet_atoms.len() {
            let center = design.planet_center_m(index);
            motions.push(BodyMotion {
                spin_center_m: [center[0], center[1], 0.0],
                spin_rad_per_s: carrier_rate + planet_relative_rate,
                orbit_center_m: origin,
                orbit_rad_per_s: carrier_rate,
            });
        }
        motions.push(BodyMotion::FIXED);

        let period_s = if carrier_rate.abs() > 0.0 {
            std::f64::consts::TAU / carrier_rate.abs()
        } else {
            std::f64::consts::TAU / planet_relative_rate.abs()
        };

        Self {
            positions_m,
            groups,
            motions,
            period_s,
        }
    }

    /// Returns every atom position at the given time.
    pub fn transformed(&self, time_s: f64) -> Vec<[f64; 3]> {
        let mut out = self.positions_m.clone();
        for (group, (_, range)) in self.groups.iter().enumerate() {
            let motion = self.motions[group];
            for index in range.clone() {
                out[index] = motion.pose(self.positions_m[index], time_s);
            }
        }
        out
    }
}

/// The settings of a clearance measurement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ClearanceTarget {
    /// The least allowed atom-to-atom distance in metres.
    pub minimum_m: f64,
    /// The number of time samples over one period.
    pub samples: usize,
    /// The search radius in metres. A pair farther apart is not reported.
    pub search_m: f64,
}

impl Default for ClearanceTarget {
    fn default() -> Self {
        Self {
            minimum_m: 2.52e-10,
            samples: 48,
            search_m: 6.0e-10,
        }
    }
}

/// The result of a clearance measurement.
#[derive(Clone, Debug, PartialEq)]
pub struct ClearanceReport {
    /// The least distance found in metres. It equals `search_m` when nothing
    /// was found inside the search radius.
    pub minimum_m: f64,
    /// True when the sweep found a pair inside the search radius.
    pub found: bool,
    /// The group name of the first body at the least distance.
    pub body_a: String,
    /// The group name of the second body at the least distance.
    pub body_b: String,
    /// The time of the least distance, in seconds.
    pub time_s: f64,
    /// The number of time samples.
    pub samples: usize,
    /// The search radius in metres.
    pub search_m: f64,
}

impl ClearanceReport {
    /// Returns true when the clearance meets the target.
    pub fn passes(&self, target: &ClearanceTarget) -> bool {
        self.found && self.minimum_m >= target.minimum_m
    }

    /// Converts the report to a metric value.
    pub fn to_metric_value(&self, target: &ClearanceTarget) -> MetricValue {
        let verdict = if self.passes(target) {
            "passes"
        } else {
            "fails"
        };
        let note = if self.found {
            format!(
                "least distance between {} and {} at {:.4e} s over {} samples; {verdict} the {:.2e} m target",
                self.body_a, self.body_b, self.time_s, self.samples, target.minimum_m
            )
        } else {
            format!(
                "no pair within the {:.2e} m search radius over {} samples; {verdict} the {:.2e} m target",
                self.search_m, self.samples, target.minimum_m
            )
        };
        MetricValue {
            name: "clearance".to_string(),
            value: self.minimum_m,
            unit: "m".to_string(),
            fidelity: Fidelity::Geometric,
            note,
        }
    }
}

/// The clearance metric.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Clearance {
    /// The settings of the measurement.
    pub target: ClearanceTarget,
}

impl Clearance {
    /// A clearance metric with the given settings.
    pub fn new(target: ClearanceTarget) -> Self {
        Self { target }
    }

    /// Measures the least cross-body distance over one motion period.
    pub fn measure(&self, moving: &MovingAtoms) -> ClearanceReport {
        let samples = self.target.samples.max(1);
        let cell = self.target.search_m;
        let mut minimum_m = f64::INFINITY;
        let mut found = false;
        let mut body_a = String::new();
        let mut body_b = String::new();
        let mut time_s = 0.0;

        for sample in 0..samples {
            let time = moving.period_s * sample as f64 / samples as f64;
            let positions = moving.transformed(time);
            for (index, (name_a, range_a)) in moving.groups.iter().enumerate() {
                for (name_b, range_b) in moving.groups.iter().skip(index + 1) {
                    let a = &positions[range_a.clone()];
                    let b = &positions[range_b.clone()];
                    if let Some(distance) = min_pair_distance(a, b, cell) {
                        if distance < minimum_m {
                            minimum_m = distance;
                            found = true;
                            body_a = name_a.clone();
                            body_b = name_b.clone();
                            time_s = time;
                        }
                    }
                }
            }
        }

        ClearanceReport {
            minimum_m: if found {
                minimum_m
            } else {
                self.target.search_m
            },
            found,
            body_a,
            body_b,
            time_s,
            samples,
            search_m: self.target.search_m,
        }
    }
}

fn min_pair_distance(a: &[[f64; 3]], b: &[[f64; 3]], cell: f64) -> Option<f64> {
    let mut grid: HashMap<Cell, Vec<u32>> = HashMap::new();
    for (index, point) in a.iter().enumerate() {
        grid.entry(cell_of(*point, cell))
            .or_default()
            .push(index as u32);
    }
    let mut best = f64::INFINITY;
    for point in b {
        let (cx, cy, cz) = cell_of(*point, cell);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    if let Some(indices) = grid.get(&(cx + dx, cy + dy, cz + dz)) {
                        for index in indices {
                            let distance = distance(a[*index as usize], *point);
                            if distance < best {
                                best = distance;
                            }
                        }
                    }
                }
            }
        }
    }
    best.is_finite().then_some(best)
}

fn cell_of(point: [f64; 3], cell: f64) -> Cell {
    (
        (point[0] / cell).floor() as i64,
        (point[1] / cell).floor() as i64,
        (point[2] / cell).floor() as i64,
    )
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_body_spins_about_its_centre() {
        let motion = BodyMotion {
            spin_center_m: [0.0, 0.0, 0.0],
            spin_rad_per_s: std::f64::consts::PI,
            orbit_center_m: [0.0, 0.0, 0.0],
            orbit_rad_per_s: 0.0,
        };
        let moved = motion.pose([1.0, 0.0, 0.0], 1.0);
        assert!((moved[0] + 1.0).abs() < 1e-9);
        assert!(moved[1].abs() < 1e-9);
    }

    #[test]
    fn a_body_orbits_about_the_orbit_centre() {
        let motion = BodyMotion {
            spin_center_m: [1.0, 0.0, 0.0],
            spin_rad_per_s: 0.0,
            orbit_center_m: [0.0, 0.0, 0.0],
            orbit_rad_per_s: std::f64::consts::PI,
        };
        let moved = motion.pose([1.0, 0.0, 0.0], 1.0);
        assert!((moved[0] + 1.0).abs() < 1e-9);
        assert!(moved[1].abs() < 1e-9);
    }

    #[test]
    fn the_pair_search_finds_the_closest_pair() {
        let a = [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let b = [[5.0, 0.2, 0.0]];
        let distance = min_pair_distance(&a, &b, 1.0).unwrap();
        assert!((distance - 0.2).abs() < 1e-9);
    }

    #[test]
    fn the_pair_search_reports_nothing_beyond_the_search_radius() {
        let a = [[0.0, 0.0, 0.0]];
        let b = [[10.0, 0.0, 0.0]];
        assert!(min_pair_distance(&a, &b, 1.0).is_none());
    }

    #[test]
    fn a_static_pair_reports_its_distance() {
        let moving = MovingAtoms {
            positions_m: vec![[0.0, 0.0, 0.0], [3.0e-10, 0.0, 0.0]],
            groups: vec![("a".to_string(), 0..1), ("b".to_string(), 1..2)],
            motions: vec![BodyMotion::FIXED, BodyMotion::FIXED],
            period_s: 1.0,
        };
        let report = Clearance::default().measure(&moving);
        assert!(report.found);
        assert!((report.minimum_m - 3.0e-10).abs() < 1e-15);
        assert_eq!(report.body_a, "a");
        assert_eq!(report.body_b, "b");
    }

    #[test]
    fn the_sweep_reduces_the_distance() {
        let moving = MovingAtoms {
            positions_m: vec![[0.0, 0.0, 0.0], [2.0e-10, 0.0, 0.0]],
            groups: vec![("a".to_string(), 0..1), ("b".to_string(), 1..2)],
            motions: vec![
                BodyMotion {
                    spin_center_m: [1.0e-10, 0.0, 0.0],
                    spin_rad_per_s: std::f64::consts::PI,
                    orbit_center_m: [0.0, 0.0, 0.0],
                    orbit_rad_per_s: 0.0,
                },
                BodyMotion::FIXED,
            ],
            period_s: 1.0,
        };
        let report = Clearance::default().measure(&moving);
        assert!(report.found);
        assert!(report.minimum_m < 1.0e-10);
    }
}
