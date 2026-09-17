//! The relaxed slip barrier: the same sweep, with the contact atoms free.
//!
//! The rigid slip barrier holds every atom on its gear path, so it is an upper
//! bound. This metric lets the atoms that face the mesh move. Each free atom
//! feels the interaction of the other body and a harmonic tether to its own
//! rigid position. The tether is the elastic response of the diamond lattice.
//!
//! The metric sweeps one sun tooth pitch and relaxes the contact atoms at every
//! step. The barrier is the difference between the largest and the smallest
//! interaction energy. A stiff lattice relaxes little, so the relaxed barrier
//! stays close to the rigid upper bound. The difference between the two values
//! measures how much the lattice compliance matters.

use std::collections::HashMap;

use crate::clearance::MovingAtoms;
use crate::slip::{cell_of, lennard_jones, Cell, BOLTZMANN_J_PER_K, WELL_DEPTH_J, ZERO_CROSSING_M};
use crate::{Fidelity, MetricValue};

/// The largest displacement of a free atom in one relaxation step, in metres.
const MAX_RELAX_STEP_M: f64 = 4.0e-12;

/// The settings of a relaxed slip-barrier measurement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RelaxedSlipBarrierTarget {
    /// The number of time samples over one sun tooth pitch.
    ///
    /// The samples must resolve the atom spacing, or the barrier aliases.
    /// One sun tooth pitch spans about fifty atoms, so 120 samples is the
    /// minimum useful count.
    pub steps: usize,
    /// The interaction cutoff in metres. A farther pair does not interact.
    pub cutoff_m: f64,
    /// The radius that selects the atoms that face the mesh, in metres.
    pub contact_m: f64,
    /// The tether stiffness of one free atom, in newtons per metre.
    pub tether_n_per_m: f64,
    /// The number of relaxation steps at each time sample.
    pub relax_iterations: usize,
    /// The temperature of the thermal energy, in kelvin.
    pub temperature_k: f64,
}

impl Default for RelaxedSlipBarrierTarget {
    fn default() -> Self {
        Self {
            steps: 120,
            cutoff_m: 1.2e-9,
            contact_m: 1.5e-9,
            tether_n_per_m: 300.0,
            relax_iterations: 40,
            temperature_k: 300.0,
        }
    }
}

/// The result of a relaxed slip-barrier measurement.
#[derive(Clone, Debug, PartialEq)]
pub struct RelaxedSlipBarrierReport {
    /// The interaction barrier in joules.
    pub barrier_j: f64,
    /// The smallest interaction energy on the sweep, in joules.
    pub minimum_j: f64,
    /// The largest interaction energy on the sweep, in joules.
    pub maximum_j: f64,
    /// The step of the largest interaction energy.
    pub step: usize,
    /// The number of time samples.
    pub steps: usize,
    /// The number of free atoms at the mesh.
    pub free_atoms: usize,
    /// The number of interacting pairs at the step of the largest energy.
    pub pairs: usize,
    /// The largest displacement of a free atom from its rigid path, in metres.
    pub relaxation_m: f64,
    /// The barrier of the rigid path on the same samples, in joules.
    pub rigid_barrier_j: f64,
}

impl RelaxedSlipBarrierReport {
    /// Returns the barrier in units of the thermal energy `k T`.
    pub fn barrier_over_kt(&self, temperature_k: f64) -> f64 {
        self.barrier_j / (BOLTZMANN_J_PER_K * temperature_k)
    }

    /// Converts the report to a metric value.
    pub fn to_metric_value(&self, target: &RelaxedSlipBarrierTarget) -> MetricValue {
        let ratio = self.barrier_over_kt(target.temperature_k);
        let verdict = if ratio < 5.0 {
            "thermal motion smooths the slip"
        } else if ratio < 20.0 {
            "the slip is marginal"
        } else {
            "the surfaces jam"
        };
        let change = if self.rigid_barrier_j > 0.0 {
            100.0 * (self.barrier_j - self.rigid_barrier_j) / self.rigid_barrier_j
        } else {
            0.0
        };
        MetricValue {
            name: "relaxed_slip_barrier".to_string(),
            value: self.barrier_j,
            unit: "J".to_string(),
            fidelity: Fidelity::QuasiStatic,
            note: format!(
                "sun to planet_0 over one sun tooth pitch in {} steps with {} free atoms on a {:.0} N/m tether; the barrier is {:.1} kT at {:.0} K, so {verdict}; relaxation changes the rigid bound by {:+.2} percent; {} pairs at the peak",
                self.steps,
                self.free_atoms,
                target.tether_n_per_m,
                ratio,
                target.temperature_k,
                change,
                self.pairs
            ),
        }
    }
}

/// The relaxed slip-barrier metric.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RelaxedSlipBarrier {
    /// The settings of the measurement.
    pub target: RelaxedSlipBarrierTarget,
}

impl RelaxedSlipBarrier {
    /// A relaxed slip-barrier metric with the given settings.
    pub fn new(target: RelaxedSlipBarrierTarget) -> Self {
        Self { target }
    }

    /// Measures the barrier between the sun and the first planet of a set.
    ///
    /// The motion of `moving` must list the sun first and then the planets.
    /// The sweep covers one sun tooth pitch.
    pub fn measure(&self, moving: &MovingAtoms, sun_teeth: usize) -> RelaxedSlipBarrierReport {
        let steps = self.target.steps.max(1);
        let sun_rate = moving.motions[0].spin_rad_per_s.abs();
        if sun_rate == 0.0 {
            return RelaxedSlipBarrierReport {
                barrier_j: 0.0,
                minimum_j: 0.0,
                maximum_j: 0.0,
                step: 0,
                steps,
                free_atoms: 0,
                pairs: 0,
                relaxation_m: 0.0,
                rigid_barrier_j: 0.0,
            };
        }
        let sun_range = moving.groups[0].1.clone();
        let planet_range = moving.groups[1].1.clone();
        let reference = moving.transformed(0.0);
        let (sun_contact, planet_contact) =
            self.contact_atoms(&reference, &sun_range, &planet_range);

        let pitch_rad = std::f64::consts::TAU / sun_teeth.max(1) as f64;
        let span_s = pitch_rad / sun_rate;

        let mut minimum_j = f64::INFINITY;
        let mut maximum_j = f64::NEG_INFINITY;
        let mut rigid_minimum_j = f64::INFINITY;
        let mut rigid_maximum_j = f64::NEG_INFINITY;
        let mut peak_step = 0;
        let mut peak_pairs = 0;
        let mut max_relaxation_m = 0.0;
        for step in 0..=steps {
            let time_s = span_s * step as f64 / steps as f64;
            let positions = moving.transformed(time_s);
            let sun_points: Vec<[f64; 3]> = sun_contact.iter().map(|&i| positions[i]).collect();
            let anchors: Vec<[f64; 3]> = planet_contact.iter().map(|&i| positions[i]).collect();
            let relaxed = self.relax(&sun_points, &anchors);
            let full_sun: Vec<[f64; 3]> = sun_range.clone().map(|i| positions[i]).collect();
            let full_planet: Vec<[f64; 3]> = planet_range.clone().map(|i| positions[i]).collect();
            let mut planet_relaxed = full_planet.clone();
            for (slot, &atom) in planet_contact.iter().enumerate() {
                planet_relaxed[atom - planet_range.start] = relaxed[slot];
            }
            for (atom, point) in planet_contact.iter().zip(&relaxed) {
                let dx = point[0] - positions[*atom][0];
                let dy = point[1] - positions[*atom][1];
                let dz = point[2] - positions[*atom][2];
                let distance = (dx * dx + dy * dy + dz * dz).sqrt();
                if distance > max_relaxation_m {
                    max_relaxation_m = distance;
                }
            }
            let (interaction_j, pairs) = self.interaction_energy(&full_sun, &planet_relaxed);
            let tether_j = self.tether_energy(&anchors, &relaxed);
            let energy_j = interaction_j + tether_j;
            let (rigid_j, _) = self.interaction_energy(&full_sun, &full_planet);
            if energy_j < minimum_j {
                minimum_j = energy_j;
            }
            if energy_j > maximum_j {
                maximum_j = energy_j;
                peak_step = step;
                peak_pairs = pairs;
            }
            if rigid_j < rigid_minimum_j {
                rigid_minimum_j = rigid_j;
            }
            if rigid_j > rigid_maximum_j {
                rigid_maximum_j = rigid_j;
            }
        }
        RelaxedSlipBarrierReport {
            barrier_j: maximum_j - minimum_j,
            minimum_j,
            maximum_j,
            step: peak_step,
            steps,
            free_atoms: planet_contact.len(),
            pairs: peak_pairs,
            relaxation_m: max_relaxation_m,
            rigid_barrier_j: rigid_maximum_j - rigid_minimum_j,
        }
    }

    /// Returns the tether energy of a displaced point list, in joules.
    fn tether_energy(&self, anchors: &[[f64; 3]], free: &[[f64; 3]]) -> f64 {
        let tether = self.target.tether_n_per_m;
        let mut energy_j = 0.0;
        for (anchor, point) in anchors.iter().zip(free) {
            let apart = subtract(*point, *anchor);
            energy_j += 0.5 * tether * dot(apart, apart);
        }
        energy_j
    }

    /// Selects the atoms of each body that lie near the other body.
    fn contact_atoms(
        &self,
        positions_m: &[[f64; 3]],
        sun_range: &std::ops::Range<usize>,
        planet_range: &std::ops::Range<usize>,
    ) -> (Vec<usize>, Vec<usize>) {
        let radius = self.target.contact_m;
        let sun = nearest_within(positions_m, sun_range, planet_range, radius);
        let planet = nearest_within(positions_m, planet_range, sun_range, radius);
        (sun, planet)
    }

    /// Relaxes the free atoms against the fixed points and the tether.
    fn relax(&self, fixed: &[[f64; 3]], anchors: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let mut free = anchors.to_vec();
        let cutoff_m = self.target.cutoff_m;
        let cutoff_squared = cutoff_m * cutoff_m;
        let tether = self.target.tether_n_per_m;
        let grid = grid_of(fixed, cutoff_m);
        for _ in 0..self.target.relax_iterations {
            let mut forces = vec![[0.0; 3]; free.len()];
            for (index, point) in free.iter().enumerate() {
                let (cx, cy, cz) = cell_of(*point, cutoff_m);
                for dx in -1..=1 {
                    for dy in -1..=1 {
                        for dz in -1..=1 {
                            let Some(bucket) = grid.get(&(cx + dx, cy + dy, cz + dz)) else {
                                continue;
                            };
                            for &other in bucket {
                                let apart = subtract(fixed[other as usize], *point);
                                let squared = dot(apart, apart);
                                if squared >= cutoff_squared || squared == 0.0 {
                                    continue;
                                }
                                let distance = squared.sqrt();
                                let force = lennard_jones_force(distance);
                                forces[index] = add(forces[index], scale(apart, force / distance));
                            }
                        }
                    }
                }
            }
            for (index, point) in free.iter_mut().enumerate() {
                let pull = subtract(anchors[index], *point);
                let mut move_vector = add(pull, scale(forces[index], 1.0 / tether));
                let length = norm(move_vector);
                if length > MAX_RELAX_STEP_M {
                    move_vector = scale(move_vector, MAX_RELAX_STEP_M / length);
                }
                *point = add(*point, move_vector);
            }
        }
        free
    }

    /// Sums the shifted interaction energy between two point lists.
    fn interaction_energy(&self, fixed: &[[f64; 3]], free: &[[f64; 3]]) -> (f64, usize) {
        let cutoff_m = self.target.cutoff_m;
        let cutoff_squared = cutoff_m * cutoff_m;
        let shift = lennard_jones(cutoff_m);
        let grid = grid_of(fixed, cutoff_m);
        let mut energy_j = 0.0;
        let mut pairs = 0usize;
        for point in free {
            let (cx, cy, cz) = cell_of(*point, cutoff_m);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let Some(bucket) = grid.get(&(cx + dx, cy + dy, cz + dz)) else {
                            continue;
                        };
                        for &other in bucket {
                            let apart = subtract(fixed[other as usize], *point);
                            let squared = dot(apart, apart);
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

/// Returns the atoms of `probe` with a neighbour in `target` inside `radius`.
fn nearest_within(
    positions_m: &[[f64; 3]],
    probe: &std::ops::Range<usize>,
    target: &std::ops::Range<usize>,
    radius: f64,
) -> Vec<usize> {
    let radius_squared = radius * radius;
    let grid = grid_of(&positions_m[target.clone()], radius);
    let mut selected = Vec::new();
    for index in probe.clone() {
        let point = positions_m[index];
        let (cx, cy, cz) = cell_of(point, radius);
        let mut near = false;
        'cells: for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let Some(bucket) = grid.get(&(cx + dx, cy + dy, cz + dz)) else {
                        continue;
                    };
                    for &other in bucket {
                        let apart = subtract(positions_m[target.start + other as usize], point);
                        if dot(apart, apart) < radius_squared {
                            near = true;
                            break 'cells;
                        }
                    }
                }
            }
        }
        if near {
            selected.push(index);
        }
    }
    selected
}

/// Builds a cell grid over a point list.
fn grid_of(points: &[[f64; 3]], cell_m: f64) -> HashMap<Cell, Vec<u32>> {
    let mut grid: HashMap<Cell, Vec<u32>> = HashMap::new();
    for (index, point) in points.iter().enumerate() {
        grid.entry(cell_of(*point, cell_m))
            .or_default()
            .push(index as u32);
    }
    grid
}

/// Returns the Lennard-Jones force of one pair, in newtons.
fn lennard_jones_force(distance_m: f64) -> f64 {
    let sigma_over_r = ZERO_CROSSING_M / distance_m;
    let ratio6 = sigma_over_r.powi(6);
    4.0 * WELL_DEPTH_J * (12.0 * ratio6 * ratio6 - 6.0 * ratio6) / distance_m
}

fn subtract(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale(a: [f64; 3], factor: f64) -> [f64; 3] {
    [a[0] * factor, a[1] * factor, a[2] * factor]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::clearance::BodyMotion;

    /// A pair farther apart than the contact radius has no free atom.
    #[test]
    fn a_far_pair_has_no_free_atom() {
        let moving = two_atom_pair(5.0e-9);
        let report = RelaxedSlipBarrier::default().measure(&moving, 12);
        assert_eq!(report.free_atoms, 0);
        assert_eq!(report.barrier_j, 0.0);
    }

    /// A close pair relaxes and shows a barrier.
    #[test]
    fn a_turning_pair_relaxes() {
        let moving = two_atom_pair(4.0e-10);
        let report = RelaxedSlipBarrier::default().measure(&moving, 12);
        assert!(report.free_atoms > 0, "no free atom at the mesh");
        assert!(report.pairs > 0, "no interacting pair");
        assert!(
            report.barrier_j > 0.0,
            "the barrier is {} J",
            report.barrier_j
        );
        assert!(report.relaxation_m > 0.0, "the atoms did not relax");
        assert!(report.relaxation_m <= MAX_RELAX_STEP_M * 2.0);
    }

    /// A stiff tether keeps the relaxed barrier near the rigid barrier.
    #[test]
    fn a_stiff_tether_keeps_the_barrier() {
        let moving = two_atom_pair(4.5e-10);
        let report = RelaxedSlipBarrier::default().measure(&moving, 12);
        assert!(report.rigid_barrier_j > 0.0);
        let change = (report.barrier_j - report.rigid_barrier_j).abs() / report.rigid_barrier_j;
        assert!(change < 0.5, "the barrier moved by {change}");
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
