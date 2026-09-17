//! The loaded-contact metric: friction and wear under a stated normal load.
//!
//! The slip barrier reports a barrier with no normal load. This metric presses
//! the sun and the first planet together with a stated normal load, then slides
//! them through one sun tooth pitch. It reports the tangential friction force
//! and a wear count.
//!
//! At every time sample the metric finds the closest sun-planet pair and sets
//! the contact normal from that pair. It displaces the planet group along the
//! normal until the derivative of the interaction energy equals the target
//! load. The pressed energy along the sweep gives the lateral force.
//!
//! The bodies are rigid in this model. The atoms follow the gear motion, and
//! the lattice does not relax. The interaction uses one carbon-like well for
//! every atom, so the value compares designs. The result is a screening
//! estimate, not a measurement.

use std::collections::{HashMap, HashSet};

use crate::clearance::MovingAtoms;
use crate::slip::{cell_of, lennard_jones, Cell};
use crate::{Fidelity, MetricValue};

/// The finite-difference step of the normal force, in metres.
const FORCE_STEP_M: f64 = 1.0e-12;

/// The lower end of the press-offset bracket, in metres.
const PRESS_MIN_M: f64 = -1.0e-10;

/// The upper end of the press-offset bracket, in metres.
const PRESS_MAX_M: f64 = 3.0e-10;

/// The number of bisection iterations of the press offset.
const PRESS_ITERATIONS: usize = 40;

/// The settings of a loaded-contact measurement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LoadedContactTarget {
    /// The number of time samples over one sun tooth pitch.
    pub steps: usize,
    /// The interaction cutoff in metres. A farther pair does not interact.
    pub cutoff_m: f64,
    /// The normal load that presses the two bodies together, in newtons.
    pub load_n: f64,
    /// The wear distance threshold in metres. A pair that comes closer wears.
    pub wear_m: f64,
    /// The temperature of the thermal energy, in kelvin.
    pub temperature_k: f64,
}

impl Default for LoadedContactTarget {
    fn default() -> Self {
        Self {
            steps: 60,
            cutoff_m: 1.2e-9,
            load_n: 1.0e-11,
            wear_m: 3.0e-10,
            temperature_k: 300.0,
        }
    }
}

/// The result of a loaded-contact measurement.
#[derive(Clone, Debug, PartialEq)]
pub struct LoadedContactReport {
    /// The mean absolute lateral force over the sweep, in newtons.
    pub friction_force_n: f64,
    /// The largest absolute lateral force over the sweep, in newtons.
    pub peak_force_n: f64,
    /// The mean friction force over the normal load. Zero when the load is zero.
    pub friction_coefficient: f64,
    /// The number of unique sun-planet pairs that came within `wear_m`.
    pub wear_pairs: usize,
    /// The number of sun-planet pairs inside the cutoff at the last sample.
    pub contact_pairs: usize,
    /// The least sun-planet distance after the press, in metres.
    pub separation_m: f64,
    /// The number of time samples over one sun tooth pitch.
    pub steps: usize,
    /// The normal load in newtons.
    pub load_n: f64,
    /// The wear distance threshold in metres.
    pub wear_m: f64,
}

impl LoadedContactReport {
    /// Converts the report to a metric value.
    pub fn to_metric_value(&self) -> MetricValue {
        MetricValue {
            name: "contact_friction".to_string(),
            value: self.friction_force_n,
            unit: "N".to_string(),
            fidelity: Fidelity::QuasiStatic,
            note: format!(
                "sun to planet_0 over one sun tooth pitch in {} steps at a {:.2e} N load; the mean friction force is {:.3e} N ({:.4} of the load); {} pairs came within {:.1e} m during the slide; the bodies are rigid, so this is a screening estimate",
                self.steps,
                self.load_n,
                self.friction_force_n,
                self.friction_coefficient,
                self.wear_pairs,
                self.wear_m
            ),
        }
    }
}

/// The loaded-contact metric.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct LoadedContact {
    /// The settings of the measurement.
    pub target: LoadedContactTarget,
}

impl LoadedContact {
    /// A loaded-contact metric with the given settings.
    pub fn with_target(target: LoadedContactTarget) -> Self {
        Self { target }
    }

    /// Measures friction and wear between the sun and the first planet.
    ///
    /// The motion of `moving` must list the sun first and then the planets.
    /// `sun_teeth` is the tooth count of the sun. The sweep covers one sun
    /// tooth pitch.
    pub fn measure(&self, moving: &MovingAtoms, sun_teeth: usize) -> LoadedContactReport {
        let steps = self.target.steps.max(1);
        let Some(sun_group) = moving.groups.first() else {
            return self.empty_report(steps);
        };
        let Some(planet_group) = moving.groups.get(1) else {
            return self.empty_report(steps);
        };
        let Some(sun_motion) = moving.motions.first() else {
            return self.empty_report(steps);
        };
        let sun_rate = sun_motion.spin_rad_per_s.abs();
        if sun_rate == 0.0 {
            return self.empty_report(steps);
        }
        let sun_range = sun_group.1.clone();
        let planet_range = planet_group.1.clone();
        let spin_center_m = sun_motion.spin_center_m;

        let span_s = (std::f64::consts::TAU / sun_teeth.max(1) as f64) / sun_rate;
        let cutoff_m = self.target.cutoff_m;
        let load_n = self.target.load_n;
        let wear_m = self.target.wear_m;

        let mut samples: Vec<(f64, f64)> = Vec::new();
        let mut wear: HashSet<(usize, usize)> = HashSet::new();
        let mut scratch: HashSet<(usize, usize)> = HashSet::new();
        let mut separation_m = f64::INFINITY;
        let mut contact_pairs = 0usize;

        for step in 0..=steps {
            let time_s = span_s * step as f64 / steps as f64;
            let positions = moving.transformed(time_s);
            let sun = &positions[sun_range.clone()];
            let planet = &positions[planet_range.clone()];

            if step == steps {
                contact_pairs = 0;
            }

            let Some((sun_local, planet_local, distance_m)) =
                closest_contact(sun, planet, cutoff_m)
            else {
                continue;
            };
            if distance_m == 0.0 {
                continue;
            }
            let normal_m = unit(subtract(sun[sun_local], planet[planet_local]));

            let scan = Interaction {
                grid: grid_of(sun, cutoff_m),
                sun,
                planet,
                sun_start: sun_range.start,
                planet_start: planet_range.start,
                cutoff_m,
                wear_m,
            };
            scratch.clear();
            let mut low = PRESS_MIN_M;
            let mut high = PRESS_MAX_M;
            let f_low = scan.normal_force(normal_m, low, &mut scratch) - load_n;
            let f_high = scan.normal_force(normal_m, high, &mut scratch) - load_n;
            let offset_m = if f_low * f_high <= 0.0 {
                for _ in 0..PRESS_ITERATIONS {
                    let mid = 0.5 * (low + high);
                    let f_mid = scan.normal_force(normal_m, mid, &mut scratch) - load_n;
                    if f_mid < 0.0 {
                        low = mid;
                    } else {
                        high = mid;
                    }
                }
                0.5 * (low + high)
            } else {
                0.0
            };
            scratch.clear();

            let (energy_j, pairs) = scan.energy_and_pairs(normal_m, offset_m, &mut wear);

            let pressed: Vec<[f64; 3]> = planet
                .iter()
                .map(|point| add(*point, scale(normal_m, offset_m)))
                .collect();
            if let Some((_, _, pressed_m)) = closest_contact(sun, &pressed, cutoff_m) {
                if pressed_m < separation_m {
                    separation_m = pressed_m;
                }
            }

            let contact_radius_m = distance(sun[sun_local], spin_center_m);
            let angle_rad = sun_motion.spin_rad_per_s * time_s;
            samples.push((energy_j, contact_radius_m * angle_rad.abs()));
            if step == steps {
                contact_pairs = pairs;
            }
        }

        let mut sum_abs = 0.0;
        let mut count = 0usize;
        let mut peak = 0.0;
        for window in samples.windows(2) {
            let (energy_0, arc_0) = window[0];
            let (energy_1, arc_1) = window[1];
            let ds = arc_1 - arc_0;
            if ds == 0.0 {
                continue;
            }
            let force = -(energy_1 - energy_0) / ds;
            let magnitude = force.abs();
            sum_abs += magnitude;
            count += 1;
            if magnitude > peak {
                peak = magnitude;
            }
        }
        let friction_force_n = if count > 0 {
            sum_abs / count as f64
        } else {
            0.0
        };
        let friction_coefficient = if load_n > 0.0 {
            friction_force_n / load_n
        } else {
            0.0
        };

        LoadedContactReport {
            friction_force_n,
            peak_force_n: peak,
            friction_coefficient,
            wear_pairs: wear.len(),
            contact_pairs,
            separation_m: if separation_m.is_finite() {
                separation_m
            } else {
                0.0
            },
            steps,
            load_n,
            wear_m,
        }
    }

    /// Returns a report with no contact.
    fn empty_report(&self, steps: usize) -> LoadedContactReport {
        LoadedContactReport {
            friction_force_n: 0.0,
            peak_force_n: 0.0,
            friction_coefficient: 0.0,
            wear_pairs: 0,
            contact_pairs: 0,
            separation_m: 0.0,
            steps,
            load_n: self.target.load_n,
            wear_m: self.target.wear_m,
        }
    }
}

/// One sun-to-planet interaction scan at a fixed time sample.
struct Interaction<'a> {
    /// The cell grid over the sun atoms.
    grid: HashMap<Cell, Vec<u32>>,
    /// The sun atom positions, in metres.
    sun: &'a [[f64; 3]],
    /// The planet atom positions, in metres.
    planet: &'a [[f64; 3]],
    /// The global index of the first sun atom.
    sun_start: usize,
    /// The global index of the first planet atom.
    planet_start: usize,
    /// The interaction cutoff in metres.
    cutoff_m: f64,
    /// The wear distance threshold in metres.
    wear_m: f64,
}

impl Interaction<'_> {
    /// Sums the shifted interaction energy with the planet displaced.
    ///
    /// The planet group moves by `offset_m * normal_m`. The function returns
    /// the energy in joules and the number of pairs inside the cutoff. It
    /// records every pair closer than `wear_m` in the wear set.
    fn energy_and_pairs(
        &self,
        normal_m: [f64; 3],
        offset_m: f64,
        wear: &mut HashSet<(usize, usize)>,
    ) -> (f64, usize) {
        let cutoff_squared = self.cutoff_m * self.cutoff_m;
        let shift = lennard_jones(self.cutoff_m);
        let mut energy_j = 0.0;
        let mut pairs = 0usize;
        for (planet_index, point) in self.planet.iter().enumerate() {
            let moved = add(*point, scale(normal_m, offset_m));
            let (cx, cy, cz) = cell_of(moved, self.cutoff_m);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let Some(bucket) = self.grid.get(&(cx + dx, cy + dy, cz + dz)) else {
                            continue;
                        };
                        for &sun_index in bucket {
                            let apart = subtract(self.sun[sun_index as usize], moved);
                            let squared = dot(apart, apart);
                            if squared >= cutoff_squared {
                                continue;
                            }
                            let distance_m = squared.sqrt();
                            energy_j += lennard_jones(distance_m) - shift;
                            pairs += 1;
                            if distance_m < self.wear_m {
                                wear.insert((
                                    self.sun_start + sun_index as usize,
                                    self.planet_start + planet_index,
                                ));
                            }
                        }
                    }
                }
            }
        }
        (energy_j, pairs)
    }

    /// Returns the normal derivative of the energy, in newtons.
    fn normal_force(
        &self,
        normal_m: [f64; 3],
        offset_m: f64,
        wear: &mut HashSet<(usize, usize)>,
    ) -> f64 {
        let plus = self
            .energy_and_pairs(normal_m, offset_m + FORCE_STEP_M, wear)
            .0;
        let minus = self
            .energy_and_pairs(normal_m, offset_m - FORCE_STEP_M, wear)
            .0;
        (plus - minus) / (2.0 * FORCE_STEP_M)
    }
}

/// Returns the closest sun-planet pair and its distance, in metres.
///
/// The search builds a grid over the planet atoms, so it does not scan every
/// pair.
fn closest_contact(
    sun: &[[f64; 3]],
    planet: &[[f64; 3]],
    cutoff_m: f64,
) -> Option<(usize, usize, f64)> {
    let mut grid: HashMap<Cell, Vec<u32>> = HashMap::new();
    for (index, point) in planet.iter().enumerate() {
        grid.entry(cell_of(*point, cutoff_m))
            .or_default()
            .push(index as u32);
    }
    let mut best = f64::INFINITY;
    let mut pair = None;
    for (sun_index, sun_point) in sun.iter().enumerate() {
        let (cx, cy, cz) = cell_of(*sun_point, cutoff_m);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let Some(bucket) = grid.get(&(cx + dx, cy + dy, cz + dz)) else {
                        continue;
                    };
                    for &planet_index in bucket {
                        let distance_m = distance(*sun_point, planet[planet_index as usize]);
                        if distance_m < best {
                            best = distance_m;
                            pair = Some((sun_index, planet_index as usize, distance_m));
                        }
                    }
                }
            }
        }
    }
    pair
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

/// Returns the distance between two points, in metres.
fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    dot(subtract(a, b), subtract(a, b)).sqrt()
}

/// Returns the unit vector of a vector. A zero vector maps to zero.
fn unit(a: [f64; 3]) -> [f64; 3] {
    let length = dot(a, a).sqrt();
    if length == 0.0 {
        [0.0; 3]
    } else {
        scale(a, 1.0 / length)
    }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::clearance::BodyMotion;

    /// A pair beyond the cutoff has no contact, no friction, and no wear.
    #[test]
    fn a_far_pair_has_no_friction() {
        let report = LoadedContact::default().measure(&two_atom_pair(5.0e-9), 12);
        assert_eq!(report.contact_pairs, 0);
        assert_eq!(report.friction_force_n, 0.0);
        assert_eq!(report.wear_pairs, 0);
    }

    /// A pair inside the cutoff is pressed and shows a lateral force.
    #[test]
    fn a_pressed_pair_has_friction() {
        let report = LoadedContact::default().measure(&two_atom_pair(3.2e-10), 12);
        assert!(report.contact_pairs > 0, "no contact pair");
        assert!(
            report.peak_force_n > 0.0,
            "the peak force is {}",
            report.peak_force_n
        );
    }

    /// A heavier load gives a smaller or equal pressed separation.
    #[test]
    fn a_heavier_load_presses_the_bodies_together() {
        let light = LoadedContact {
            target: LoadedContactTarget {
                load_n: 1.0e-12,
                ..Default::default()
            },
        }
        .measure(&two_atom_pair(3.2e-10), 12);
        let heavy = LoadedContact {
            target: LoadedContactTarget {
                load_n: 1.0e-10,
                ..Default::default()
            },
        }
        .measure(&two_atom_pair(3.2e-10), 12);
        assert!(
            heavy.separation_m <= light.separation_m,
            "heavy {} m, light {} m",
            heavy.separation_m,
            light.separation_m
        );
    }

    /// A touching surface pair counts as wear.
    #[test]
    fn the_wear_count_is_positive_when_the_surfaces_touch() {
        let contact = LoadedContact {
            target: LoadedContactTarget {
                wear_m: 4.0e-10,
                ..Default::default()
            },
        };
        let report = contact.measure(&two_atom_pair(3.2e-10), 12);
        assert!(report.wear_pairs > 0, "no wear pair");
    }

    /// The coefficient is the friction force over the load.
    #[test]
    fn the_coefficient_divides_by_the_load() {
        let report = LoadedContact::default().measure(&two_atom_pair(3.2e-10), 12);
        let expected = report.friction_force_n / report.load_n;
        assert!((report.friction_coefficient - expected).abs() <= 1.0e-9 * expected.abs());
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
