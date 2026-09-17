//! The slip-barrier metric: the energy barrier to sliding at a gear mesh.
//!
//! Two meshing gears slide over each other. If the potential energy changes
//! little as one tooth passes, the surfaces glide and the friction is low. If
//! the energy swings by much more than the thermal energy, the surfaces jam
//! and the motion is stick-slip.
//!
//! The metric sweeps one sun tooth pitch with the real gear kinematics. At
//! every step it sums a shifted Lennard-Jones interaction between the sun and
//! the first planet. The barrier is the difference between the largest and the
//! smallest energy on the sweep.
//!
//! The model is rigid: the atoms follow the gear motion, and the lattice does
//! not relax. A relaxed barrier is lower, so this value is an upper bound. The
//! interaction uses one carbon-like well for every atom, so the value compares
//! designs, and it is not an absolute friction coefficient.

use std::collections::HashMap;

use crate::clearance::MovingAtoms;
use crate::{Fidelity, MetricValue};

/// The Boltzmann constant in joules per kelvin.
pub const BOLTZMANN_J_PER_K: f64 = 1.380_649e-23;

/// A carbon-like Lennard-Jones well depth in joules.
pub(crate) const WELL_DEPTH_J: f64 = 5.98e-22;

/// A carbon-like Lennard-Jones zero crossing in metres.
pub(crate) const ZERO_CROSSING_M: f64 = 3.4e-10;

pub(crate) type Cell = (i64, i64, i64);

/// The settings of a slip-barrier measurement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlipBarrierTarget {
    /// The number of time samples over one sun tooth pitch.
    ///
    /// The samples must resolve the atom spacing, or the barrier aliases.
    /// One sun tooth pitch spans about fifty atoms, so 240 samples is the
    /// default count.
    pub steps: usize,
    /// The interaction cutoff in metres. A farther pair does not interact.
    pub cutoff_m: f64,
    /// The temperature of the thermal energy, in kelvin.
    pub temperature_k: f64,
}

impl Default for SlipBarrierTarget {
    fn default() -> Self {
        Self {
            steps: 240,
            cutoff_m: 1.2e-9,
            temperature_k: 300.0,
        }
    }
}

/// The result of a slip-barrier measurement.
#[derive(Clone, Debug, PartialEq)]
pub struct SlipBarrierReport {
    /// The energy barrier in joules, the largest energy minus the smallest.
    pub barrier_j: f64,
    /// The smallest energy on the sweep, in joules.
    pub minimum_j: f64,
    /// The largest energy on the sweep, in joules.
    pub maximum_j: f64,
    /// The step of the largest energy.
    pub step: usize,
    /// The number of time samples.
    pub steps: usize,
    /// The number of interacting pairs at the step of the largest energy.
    pub pairs: usize,
}

impl SlipBarrierReport {
    /// Returns the barrier in units of the thermal energy `k T`.
    pub fn barrier_over_kt(&self, temperature_k: f64) -> f64 {
        self.barrier_j / (BOLTZMANN_J_PER_K * temperature_k)
    }

    /// Converts the report to a metric value.
    pub fn to_metric_value(&self, target: &SlipBarrierTarget) -> MetricValue {
        let ratio = self.barrier_over_kt(target.temperature_k);
        let verdict = if ratio < 5.0 {
            "thermal motion smooths the slip"
        } else if ratio < 20.0 {
            "the slip is marginal"
        } else {
            "the surfaces jam"
        };
        MetricValue {
            name: "slip_barrier".to_string(),
            value: self.barrier_j,
            unit: "J".to_string(),
            fidelity: Fidelity::QuasiStatic,
            note: format!(
                "sun to planet_0 over one sun tooth pitch in {} steps; the barrier is {:.1} kT at {:.0} K, so {verdict}; {} pairs at the peak; the lattice is rigid, so this is an upper bound",
                self.steps, ratio, target.temperature_k, self.pairs
            ),
        }
    }
}

/// The slip-barrier metric.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SlipBarrier {
    /// The settings of the measurement.
    pub target: SlipBarrierTarget,
}

impl SlipBarrier {
    /// A slip-barrier metric with the given settings.
    pub fn new(target: SlipBarrierTarget) -> Self {
        Self { target }
    }

    /// Measures the barrier between the sun and the first planet of a set.
    ///
    /// The motion of `moving` must list the sun first and then the planets.
    /// `sun_teeth` is the tooth count of the sun.
    /// The sweep covers one sun tooth pitch.
    pub fn measure(&self, moving: &MovingAtoms, sun_teeth: usize) -> SlipBarrierReport {
        let steps = self.target.steps.max(1);
        let sun = &moving.groups[0].1;
        let planet = &moving.groups[1].1;
        let sun_rate = moving.motions[0].spin_rad_per_s.abs();
        if sun_rate == 0.0 {
            return SlipBarrierReport {
                barrier_j: 0.0,
                minimum_j: 0.0,
                maximum_j: 0.0,
                step: 0,
                steps,
                pairs: 0,
            };
        }
        let pitch_rad = std::f64::consts::TAU / sun_teeth.max(1) as f64;
        let span_s = pitch_rad / sun_rate;

        let mut minimum_j = f64::INFINITY;
        let mut maximum_j = f64::NEG_INFINITY;
        let mut peak_step = 0;
        let mut peak_pairs = 0;
        for step in 0..=steps {
            let time_s = span_s * step as f64 / steps as f64;
            let positions = moving.transformed(time_s);
            let (energy_j, pairs) = self.pair_energy(&positions, sun, planet);
            if energy_j < minimum_j {
                minimum_j = energy_j;
            }
            if energy_j > maximum_j {
                maximum_j = energy_j;
                peak_step = step;
                peak_pairs = pairs;
            }
        }
        SlipBarrierReport {
            barrier_j: maximum_j - minimum_j,
            minimum_j,
            maximum_j,
            step: peak_step,
            steps,
            pairs: peak_pairs,
        }
    }

    /// Sums the shifted Lennard-Jones energy between two body groups.
    ///
    /// The potential is shifted so that it is zero at the cutoff. The sum
    /// returns the energy in joules and the number of pairs inside the cutoff.
    fn pair_energy(
        &self,
        positions_m: &[[f64; 3]],
        group_a: &std::ops::Range<usize>,
        group_b: &std::ops::Range<usize>,
    ) -> (f64, usize) {
        let cutoff_m = self.target.cutoff_m;
        let cutoff_squared = cutoff_m * cutoff_m;
        let shift = lennard_jones(cutoff_m);

        let mut grid: HashMap<Cell, Vec<u32>> = HashMap::new();
        for index in group_b.clone() {
            grid.entry(cell_of(positions_m[index], cutoff_m))
                .or_default()
                .push(index as u32);
        }

        let mut energy_j = 0.0;
        let mut pairs = 0usize;
        for index in group_a.clone() {
            let point = positions_m[index];
            let (cx, cy, cz) = cell_of(point, cutoff_m);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let Some(bucket) = grid.get(&(cx + dx, cy + dy, cz + dz)) else {
                            continue;
                        };
                        for &other in bucket {
                            let other_point = positions_m[other as usize];
                            let rx = point[0] - other_point[0];
                            let ry = point[1] - other_point[1];
                            let rz = point[2] - other_point[2];
                            let squared = rx * rx + ry * ry + rz * rz;
                            if squared >= cutoff_squared {
                                continue;
                            }
                            energy_j += lennard_jones(squared.sqrt()) - shift;
                            pairs += 1;
                        }
                    }
                }
            }
        }
        (energy_j, pairs)
    }
}

/// Returns the Lennard-Jones energy of one pair, in joules.
pub(crate) fn lennard_jones(distance_m: f64) -> f64 {
    let ratio = ZERO_CROSSING_M / distance_m;
    let ratio6 = ratio * ratio * ratio * ratio * ratio * ratio;
    4.0 * WELL_DEPTH_J * (ratio6 * ratio6 - ratio6)
}

/// Returns the grid cell that holds a point, for the given cell size.
pub(crate) fn cell_of(point_m: [f64; 3], cell_m: f64) -> Cell {
    (
        (point_m[0] / cell_m).floor() as i64,
        (point_m[1] / cell_m).floor() as i64,
        (point_m[2] / cell_m).floor() as i64,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::clearance::BodyMotion;

    /// Two atoms at the well minimum have the well depth, as a negative energy.
    #[test]
    fn the_well_depth_is_the_minimum() {
        let distance = 2.0_f64.powf(1.0 / 6.0) * ZERO_CROSSING_M;
        let energy = lennard_jones(distance);
        assert!(
            (energy + WELL_DEPTH_J).abs() < WELL_DEPTH_J * 1e-9,
            "the minimum is {energy} J, not {} J",
            -WELL_DEPTH_J
        );
    }

    /// Two bodies farther apart than the cutoff do not interact.
    #[test]
    fn a_far_pair_has_no_barrier() {
        let moving = two_atom_pair(5.0e-9);
        let report = SlipBarrier::default().measure(&moving, 12);
        assert_eq!(report.pairs, 0);
        assert_eq!(report.barrier_j, 0.0);
    }

    /// The barrier is positive when two close bodies turn past each other.
    #[test]
    fn a_turning_pair_has_a_barrier() {
        let moving = two_atom_pair(4.0e-10);
        let report = SlipBarrier::default().measure(&moving, 12);
        assert!(report.pairs > 0, "no pairs inside the cutoff");
        assert!(
            report.barrier_j > 0.0,
            "the barrier is {} J",
            report.barrier_j
        );
    }

    /// Builds a moving set with one sun atom and one planet atom.
    fn two_atom_pair(separation_m: f64) -> MovingAtoms {
        let arm_m = 0.5e-9;
        let positions_m = vec![[arm_m, 0.0, 0.0], [arm_m + separation_m, 0.0, 0.0]];
        let groups = vec![("sun".to_string(), 0..1), ("planet_0".to_string(), 1..2)];
        let motions = vec![
            BodyMotion {
                spin_center_m: [0.0, 0.0, 0.0],
                spin_rad_per_s: -0.34,
                orbit_center_m: [0.0, 0.0, 0.0],
                orbit_rad_per_s: 0.0,
            },
            BodyMotion::FIXED,
        ];
        MovingAtoms {
            positions_m,
            groups,
            motions,
            period_s: 1.0,
        }
    }
}
