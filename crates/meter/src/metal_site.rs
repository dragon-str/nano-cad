//! The binding search at a metal site (M27-04, was M25-02).
//!
//! A tooth-gap pocket of the gear pump carries one iron atom at its root. A
//! guest binds when one of its donor atoms (oxygen) coordinates the iron. The
//! search places the donor on the open side of the iron, turns the guest over
//! a sphere of directions and a spin, and keeps the lowest energy.
//!
//! The energy of a pose is the classical wall term of
//! [`crate::binding`] (UFF Lennard-Jones and point-charge Coulomb) over every
//! nearby atom, plus a Morse coordination term between the donor and the iron:
//! `U(r) = D (1 - exp(-a (r - r0)))^2 - D`, which is `-D` at `r0` and zero at
//! long range. The donor-iron pair has no Lennard-Jones term: the Morse term
//! replaces it. The well depths come from the quantum-chemistry proxy of
//! ADR-0084 and ADR-0086 (Fe(NH3)4, PBE/def2-svp): oxygen 80 kJ/mol at
//! 1.90e-10 m, water 10.3 kJ/mol at 2.7e-10 m, carbon dioxide none. The proxy
//! is not the diamond root site, so every number is an estimate.
//!
//! The search is rigid: the wall, the guest and the iron do not move, and
//! there is no solvent. The energies are enthalpies; the competition model of
//! [`crate::competition`] adds the entropy.

use nanocad_model::{pair_params, Element, Part};
use nanocad_parts::Guest;

use crate::competition::AVOGADRO_PER_MOL;
use nanocad_engine::COULOMB_CONSTANT_N_M2_PER_C2;

/// The Morse width of the iron-oxygen well, in reciprocal metres.
///
/// It is `sqrt(k / (2 D))` for the oxygen well of the metal-ligand example
/// (k = 200 N/m, D = 80 kJ/mol). Every guest uses the same width, so the
/// wells differ only in depth and position.
pub const IRON_OXYGEN_MORSE_A_PER_M: f64 = 2.743_7e10;

/// The coordination of one guest at an iron site.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Coordination {
    /// The donor element of the guest.
    pub donor: Element,
    /// The well depth `D`, in joules. Positive binds.
    pub well_depth_j: f64,
    /// The iron-donor distance at the well bottom, in metres.
    pub r0_m: f64,
    /// The Morse width, in reciprocal metres.
    pub a_per_m: f64,
}

impl Coordination {
    /// The coordination from a molar well depth in kJ/mol.
    pub fn from_molar(donor: Element, well_depth_kj_per_mol: f64, r0_m: f64) -> Self {
        Self {
            donor,
            well_depth_j: well_depth_kj_per_mol * 1.0e3 / AVOGADRO_PER_MOL,
            r0_m,
            a_per_m: IRON_OXYGEN_MORSE_A_PER_M,
        }
    }

    /// The coordination energy at a distance, in joules. Negative binds.
    pub fn energy_j(&self, r_m: f64) -> f64 {
        let x = 1.0 - (-self.a_per_m * (r_m - self.r0_m)).exp();
        self.well_depth_j * (x * x - 1.0)
    }
}

/// The calibrated coordination of a guest at an iron site, when it has one.
///
/// Oxygen and water come from the quantum-chemistry proxy (ADR-0084,
/// ADR-0086). Carbon dioxide does not coordinate the proxy (+0.6 kJ/mol at
/// 3.4e-10 m), so it has none.
pub fn iron_coordination(guest_id: &str) -> Option<Coordination> {
    match guest_id {
        "o2" => Some(Coordination::from_molar(Element::OXYGEN, 80.0, 1.90e-10)),
        "water" => Some(Coordination::from_molar(Element::OXYGEN, 10.3, 2.70e-10)),
        _ => None,
    }
}

/// The settings of a metal-site search.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MetalSiteTarget {
    /// The pair cutoff, in metres.
    pub cutoff_m: f64,
    /// The iron-donor distances to try, from this value.
    pub distance_min_m: f64,
    /// The iron-donor distances to try, up to this value.
    pub distance_max_m: f64,
    /// The number of iron-donor distances.
    pub distances: usize,
    /// The number of guest directions on the sphere.
    pub directions: usize,
    /// The number of spins about each direction.
    pub spins: usize,
    /// A pair closer than this fraction of sigma blocks the pose.
    pub overlap_ratio: f64,
}

impl Default for MetalSiteTarget {
    fn default() -> Self {
        Self {
            cutoff_m: 1.0e-9,
            distance_min_m: 1.7e-10,
            distance_max_m: 3.6e-10,
            distances: 12,
            directions: 48,
            spins: 8,
            overlap_ratio: 0.85,
        }
    }
}

/// The best pose of one guest at a metal site.
#[derive(Clone, Debug, PartialEq)]
pub struct MetalSiteReport {
    /// The guest id.
    pub guest: String,
    /// The lowest pose energy, in joules. Negative binds.
    pub minimum_energy_j: f64,
    /// The coordination part of that energy, in joules.
    pub coordination_energy_j: f64,
    /// The wall part of that energy, in joules.
    pub wall_energy_j: f64,
    /// The iron-donor distance of that pose, in metres.
    pub donor_distance_m: f64,
    /// The number of poses tried.
    pub poses: usize,
    /// The number of poses blocked by an overlap.
    pub blocked: usize,
}

impl MetalSiteReport {
    /// The lowest pose energy in kJ/mol.
    pub fn minimum_energy_kj_per_mol(&self) -> f64 {
        self.minimum_energy_j * AVOGADRO_PER_MOL / 1.0e3
    }
}

/// Searches the poses of `guest` at the iron atom `iron` of `part`.
///
/// The site axis runs from the atom bonded to the iron through the iron. The
/// donor is the first guest atom of the donor element; with no coordination
/// it is the first oxygen, or atom 0. Every part atom within the cutoff of
/// the iron is the wall.
pub fn metal_site_binding(
    part: &Part,
    iron: usize,
    guest: &Guest,
    coordination: Option<&Coordination>,
    target: &MetalSiteTarget,
) -> Option<MetalSiteReport> {
    let iron_m = part.topology.position_m(iron)?;
    let host = part.topology.bonds().find_map(|bond| {
        if bond.u as usize == iron {
            Some(bond.v as usize)
        } else if bond.v as usize == iron {
            Some(bond.u as usize)
        } else {
            None
        }
    })?;
    let host_m = part.topology.position_m(host)?;
    let axis = unit(sub(iron_m, host_m))?;

    let reach_m = target.cutoff_m + target.distance_max_m + guest.extent_m();
    let mut wall = Vec::new();
    for index in 0..part.atom_count() {
        let (Some(p), Some(element)) = (
            part.topology.position_m(index),
            part.topology.element(index),
        ) else {
            continue;
        };
        if norm(sub(p, iron_m)) <= reach_m {
            wall.push((
                p,
                part.topology.charge_c(index).unwrap_or(0.0),
                element,
                index == iron,
            ));
        }
    }

    let donor_element = coordination.map_or(Element::OXYGEN, |c| c.donor);
    let donor = guest
        .atoms
        .iter()
        .position(|atom| atom.element == donor_element)
        .unwrap_or(0);
    let donor_m = guest.atoms[donor].position_m;
    // The guest body direction: from the donor to the far atoms.
    let mut body = [0.0; 3];
    for (index, atom) in guest.atoms.iter().enumerate() {
        if index != donor {
            body = add(body, sub(atom.position_m, donor_m));
        }
    }
    let body = unit(body).unwrap_or([0.0, 0.0, 1.0]);

    let mut report = MetalSiteReport {
        guest: guest.id.to_string(),
        minimum_energy_j: f64::INFINITY,
        coordination_energy_j: 0.0,
        wall_energy_j: 0.0,
        donor_distance_m: 0.0,
        poses: 0,
        blocked: 0,
    };
    let cutoff_sq = target.cutoff_m * target.cutoff_m;
    let directions = sphere_points(target.directions.max(1));
    for step in 0..target.distances.max(1) {
        let r_m = if target.distances <= 1 {
            target.distance_min_m
        } else {
            target.distance_min_m
                + (target.distance_max_m - target.distance_min_m) * step as f64
                    / (target.distances - 1) as f64
        };
        let donor_at = add(iron_m, scale(axis, r_m));
        for direction in &directions {
            let align = rotation_between(body, *direction);
            for spin in 0..target.spins.max(1) {
                let angle = std::f64::consts::TAU * spin as f64 / target.spins.max(1) as f64;
                let turn = axis_angle(*direction, angle);
                report.poses += 1;
                let mut wall_j = 0.0;
                let mut blocked = false;
                let mut coordination_j = 0.0;
                'atoms: for (index, atom) in guest.atoms.iter().enumerate() {
                    let local = sub(atom.position_m, donor_m);
                    let p = add(donor_at, apply(&turn, apply(&align, local)));
                    for &(wall_m, wall_charge_c, wall_element, is_iron) in &wall {
                        let d = sub(wall_m, p);
                        let r_sq = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                        if r_sq >= cutoff_sq || r_sq <= 0.0 {
                            continue;
                        }
                        if is_iron && index == donor {
                            if let Some(c) = coordination {
                                coordination_j = c.energy_j(r_sq.sqrt());
                                continue;
                            }
                        }
                        let Some((epsilon_j, sigma_m)) = pair_params(atom.element, wall_element)
                        else {
                            continue;
                        };
                        let floor_m = target.overlap_ratio * sigma_m;
                        if r_sq < floor_m * floor_m {
                            blocked = true;
                            break 'atoms;
                        }
                        let ratio_sq = sigma_m * sigma_m / r_sq;
                        let inverse_6 = ratio_sq * ratio_sq * ratio_sq;
                        wall_j += 4.0 * epsilon_j * (inverse_6 * inverse_6 - inverse_6);
                        wall_j += COULOMB_CONSTANT_N_M2_PER_C2 * atom.charge_c * wall_charge_c
                            / r_sq.sqrt();
                    }
                }
                if blocked {
                    report.blocked += 1;
                    continue;
                }
                let total_j = wall_j + coordination_j;
                if total_j < report.minimum_energy_j {
                    report.minimum_energy_j = total_j;
                    report.coordination_energy_j = coordination_j;
                    report.wall_energy_j = wall_j;
                    report.donor_distance_m = r_m;
                }
            }
        }
    }
    Some(report)
}

type Matrix = [[f64; 3]; 3];

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn unit(a: [f64; 3]) -> Option<[f64; 3]> {
    let n = norm(a);
    (n > 0.0).then(|| scale(a, 1.0 / n))
}

fn apply(m: &Matrix, v: [f64; 3]) -> [f64; 3] {
    [dot(m[0], v), dot(m[1], v), dot(m[2], v)]
}

/// The rotation by `angle` about the unit `axis` (Rodrigues).
fn axis_angle(axis: [f64; 3], angle: f64) -> Matrix {
    let (s, c) = angle.sin_cos();
    let t = 1.0 - c;
    let [x, y, z] = axis;
    [
        [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
        [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
        [t * x * z - s * y, t * y * z + s * x, t * z * z + c],
    ]
}

/// The rotation that turns the unit `from` onto the unit `to`.
fn rotation_between(from: [f64; 3], to: [f64; 3]) -> Matrix {
    let c = dot(from, to).clamp(-1.0, 1.0);
    let axis = cross(from, to);
    match unit(axis) {
        Some(axis) => axis_angle(axis, c.acos()),
        None if c > 0.0 => axis_angle([0.0, 0.0, 1.0], 0.0),
        None => {
            // Opposite vectors: turn by PI about any perpendicular axis.
            let helper = if from[0].abs() < 0.9 {
                [1.0, 0.0, 0.0]
            } else {
                [0.0, 1.0, 0.0]
            };
            axis_angle(
                unit(cross(from, helper)).unwrap_or([0.0, 0.0, 1.0]),
                std::f64::consts::PI,
            )
        }
    }
}

/// Near-uniform unit vectors on a sphere (a Fibonacci lattice).
fn sphere_points(count: usize) -> Vec<[f64; 3]> {
    let golden = std::f64::consts::PI * (3.0 - 5.0_f64.sqrt());
    (0..count)
        .map(|i| {
            let z = 1.0 - 2.0 * (i as f64 + 0.5) / count as f64;
            let r = (1.0 - z * z).sqrt();
            let phi = golden * i as f64;
            [r * phi.cos(), r * phi.sin(), z]
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_morse_well_has_its_depth_at_r0_and_vanishes_far_away() {
        let c = iron_coordination("o2").expect("oxygen coordinates");
        assert!((c.energy_j(c.r0_m) + c.well_depth_j).abs() < 1.0e-30);
        assert!(c.energy_j(2.0e-9).abs() < 1.0e-3 * c.well_depth_j);
        assert!((c.well_depth_j * AVOGADRO_PER_MOL / 1.0e3 - 80.0).abs() < 1.0e-9);
        assert!(iron_coordination("co2").is_none());
    }

    #[test]
    fn the_width_matches_the_example_force_constant() {
        let d_j = 80.0e3 / AVOGADRO_PER_MOL;
        let a = (200.0 / (2.0 * d_j)).sqrt();
        assert!((a - IRON_OXYGEN_MORSE_A_PER_M).abs() / a < 1.0e-3, "{a}");
    }

    /// On the pump root iron, oxygen binds through the coordination well and
    /// water binds far more weakly. A coarse search keeps the test quick.
    #[test]
    fn the_pump_root_iron_prefers_oxygen_to_water() {
        use nanocad_parts::{guest, GearPumpGenerator, ParameterSet};
        let set = GearPumpGenerator
            .build(&ParameterSet::new())
            .expect("the pump builds");
        let iron = set.metal_sites[3];
        let target = MetalSiteTarget {
            distances: 6,
            directions: 16,
            spins: 4,
            ..MetalSiteTarget::default()
        };
        let energy = |id: &str| {
            let coordination = iron_coordination(id);
            metal_site_binding(
                &set.part,
                iron,
                guest(id).expect("guest"),
                coordination.as_ref(),
                &target,
            )
            .expect("search")
            .minimum_energy_kj_per_mol()
        };
        let oxygen = energy("o2");
        let water = energy("water");
        assert!(oxygen < -60.0, "{oxygen}");
        assert!(water - oxygen > 50.0, "{oxygen} {water}");
    }

    #[test]
    fn a_rotation_turns_one_vector_onto_another() {
        for (from, to) in [
            ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
            ([0.0, 0.0, 1.0], [0.0, 0.0, -1.0]),
            ([0.0, 0.6, 0.8], [0.0, 0.6, 0.8]),
        ] {
            let m = rotation_between(from, to);
            let got = apply(&m, from);
            assert!(norm(sub(got, to)) < 1.0e-12, "{got:?} {to:?}");
        }
    }
}
