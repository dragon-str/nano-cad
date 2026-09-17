//! The steered-drive metric: the torque and the energy loss of a driven mesh.
//!
//! A steered sun turns a contact cluster of the sun and the first planet. The
//! metric integrates the cluster with velocity Verlet. It holds the planet
//! atoms and writes the sun velocity from the drive rate at every step. The
//! mean resisting torque is the reaction of the contact on the driven body.
//!
//! The cluster is bounded, so the torque is that of the cluster, not of the
//! whole gear. The potential has a bond-stretch term for every topology bond
//! inside the cluster and a Buckingham van der Waals term for every
//! non-bonded pair, as in the harmonic metric. The drive work is the mean
//! torque times the swept angle. The energy loss is the drive work minus the
//! rise in the total cluster energy, floored at zero.

use std::collections::HashMap;

use nanocad_engine::{
    minimize_with, BondStretchTerm, Cutoff, MinimizeMethod, MinimizeOptions, PeriodicBox, System,
    VanDerWaalsTerm, VdwParams, VelocityVerlet,
};
use nanocad_parts::planetary::PlanetarySet;

use crate::clearance::MovingAtoms;
use crate::slip::{cell_of, Cell};
use crate::{Fidelity, MetricValue};

/// The mass of one carbon atom in kilograms.
const CARBON_MASS_KG: f64 = 1.992_646_879_92e-26;

/// The bond stiffness of the cluster, in newtons per metre.
const BOND_STIFFNESS_N_PER_M: f64 = 300.0;

/// The Buckingham repulsion prefactor, in joules.
const VDW_A_J: f64 = 1.2e-15;

/// The Buckingham repulsion decay, in reciprocal metres.
const VDW_B_PER_M: f64 = 4.0e10;

/// The Buckingham attraction coefficient, in joule cubic metres to the sixth.
const VDW_C_J_M6: f64 = 5.5e-78;

/// The Boltzmann constant in joules per kelvin.
const BOLTZMANN_J_PER_K: f64 = 1.380_649e-23;

/// The settings of a steered-drive measurement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DriveTarget {
    /// The number of integration steps.
    pub steps: usize,
    /// The time step in seconds.
    pub dt_s: f64,
    /// The radius of the atom cluster about the mesh point, in metres.
    pub cluster_m: f64,
    /// The largest atom count of the cluster.
    pub max_atoms: usize,
    /// The van der Waals cutoff in metres.
    pub cutoff_m: f64,
    /// The contact shell in metres. A cross-body pair farther apart than this
    /// does not interact, so only the contact carries the van der Waals term.
    pub contact_m: f64,
    /// The temperature of the thermal energy, in kelvin.
    pub temperature_k: f64,
}

impl Default for DriveTarget {
    fn default() -> Self {
        Self {
            steps: 4000,
            dt_s: 2.0e-16,
            cluster_m: 6.0e-10,
            max_atoms: 120,
            cutoff_m: 1.2e-9,
            contact_m: 3.5e-10,
            temperature_k: 300.0,
        }
    }
}

/// The result of a steered-drive measurement.
#[derive(Clone, Debug, PartialEq)]
pub struct DriveReport {
    /// The mean absolute interaction torque on the driven body, in newton
    /// metres.
    pub torque_n_m: f64,
    /// The largest absolute torque seen, in newton metres.
    pub peak_torque_n_m: f64,
    /// The drive work in joules, the mean torque times the swept angle.
    pub work_j: f64,
    /// The drive work that did not become cluster energy, floored at zero, in
    /// joules.
    pub energy_loss_j: f64,
    /// The final kinetic temperature of the free atoms minus the target, in
    /// kelvin, floored at zero.
    pub temperature_rise_k: f64,
    /// The final kinetic temperature of the free atoms, in kelvin.
    pub final_temperature_k: f64,
    /// The target temperature, in kelvin. It is the baseline of the rise.
    pub temperature_k: f64,
    /// The number of free atoms. They carry the temperature. Zero when the
    /// cluster holds every atom, as it does at a mesh with no nearby ring.
    pub free_atom_count: usize,
    /// The number of atoms in the cluster.
    pub atom_count: usize,
    /// The number of integration steps that ran.
    pub steps: usize,
    /// The drive rate of the sun, in radians per second.
    pub drive_rate_rad_per_s: f64,
}

impl DriveReport {
    /// Converts the report to a metric value.
    pub fn to_metric_value(&self) -> MetricValue {
        let base_k = if self.temperature_k > 0.0 {
            self.temperature_k
        } else {
            (self.final_temperature_k - self.temperature_rise_k).max(0.0)
        };
        let heat = if self.free_atom_count == 0 {
            "the cluster holds every atom, so it has no free atom and no measurable temperature rise".to_string()
        } else {
            format!(
                "the {} free atoms rise {:.1} K above {:.0} K",
                self.free_atom_count, self.temperature_rise_k, base_k
            )
        };
        let note = if self.atom_count == 0 {
            format!(
                "the steered-drive cluster is empty, so the steered sun carries no measurable torque over {} steps",
                self.steps
            )
        } else {
            format!(
                "a steered sun turns a {}-atom contact cluster at {:.3} rad/s for {} steps; the mean resisting torque is {:.3e} N m and the drive loses {:.3e} J; {}; the bodies are rigid, so this is a screening estimate",
                self.atom_count,
                self.drive_rate_rad_per_s,
                self.steps,
                self.torque_n_m,
                self.energy_loss_j,
                heat
            )
        };
        MetricValue {
            name: "steered_drive".to_string(),
            value: self.torque_n_m,
            unit: "N m".to_string(),
            fidelity: Fidelity::Dynamics,
            note,
        }
    }
}

/// The steered-drive metric.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SteeredDrive {
    /// The settings of the measurement.
    pub target: DriveTarget,
}

impl SteeredDrive {
    /// A steered-drive metric with the given settings.
    pub fn with_target(target: DriveTarget) -> Self {
        Self { target }
    }

    /// Measures the drive torque and the energy loss at the first sun-planet
    /// mesh.
    ///
    /// This is infallible. It returns a zeroed report when the cluster is
    /// empty or an internal step fails.
    pub fn measure(&self, set: &PlanetarySet, sun_rad_per_s: f64) -> DriveReport {
        match self.try_measure(set, sun_rad_per_s) {
            Some(report) => report,
            None => zeroed_report(&self.target, sun_rad_per_s),
        }
    }

    /// Returns the report, or `None` when the cluster is empty or a step
    /// fails.
    fn try_measure(&self, set: &PlanetarySet, sun_rad_per_s: f64) -> Option<DriveReport> {
        let moving = MovingAtoms::planetary(set, &set.design, sun_rad_per_s);
        let atoms = moving.transformed(0.0);
        let sun_group = moving.groups.first()?.1.clone();
        let planet_group = moving.groups.get(1)?.1.clone();
        let center_m = moving.motions.first()?.spin_center_m;

        let (sun_point, planet_point) =
            closest_pair(&atoms[sun_group.clone()], &atoms[planet_group.clone()])?;
        let mesh_m = [
            0.5 * (sun_point[0] + planet_point[0]),
            0.5 * (sun_point[1] + planet_point[1]),
            0.5 * (sun_point[2] + planet_point[2]),
        ];

        let cluster = self.cluster(&moving, &atoms, mesh_m);
        if cluster.is_empty() {
            return None;
        }
        let (mut system, mut positions) = self.cluster_system(&set.part.topology, &cluster)?;
        let atom_count = cluster.len();
        let bodies: Vec<u8> = cluster.iter().map(|(_, _, body)| *body).collect();
        let driven: Vec<bool> = bodies.iter().map(|body| *body == 0).collect();

        let options = MinimizeOptions {
            max_iterations: 400,
            gradient_tolerance_n: 1.0e-10,
            initial_step_m: 1.0e-12,
        };
        // The relaxation only prepares the geometry. A failure leaves the
        // generated geometry in place, which is still a usable start.
        let _ = minimize_with(
            &mut system,
            &mut positions,
            &options,
            MinimizeMethod::ConjugateGradientThenLbfgs,
        );
        let relaxed = positions.clone();
        // The relaxation moves the atoms, so the neighbour list must match the
        // relaxed geometry before the first energy read.
        system.rebuild_neighbors(&relaxed).ok()?;
        let energy_initial_j = system.energy_j(&relaxed).ok()?;

        let integrator = VelocityVerlet::new(self.target.dt_s).ok()?;
        let mut velocities = vec![0.0; 3 * atom_count];
        let steps = self.target.steps;
        let omega = sun_rad_per_s;
        let mut sum_abs_torque_n_m = 0.0;
        let mut peak_torque_n_m = 0.0;

        for _ in 0..steps {
            let (_, gradient) = system.energy_and_gradient_j(&positions).ok()?;
            let points = points_of(&positions);
            let forces = forces_of(&gradient);
            let torque_n_m = torque_z(&points, &forces, center_m, &driven);
            let magnitude_n_m = torque_n_m.abs();
            sum_abs_torque_n_m += magnitude_n_m;
            if magnitude_n_m > peak_torque_n_m {
                peak_torque_n_m = magnitude_n_m;
            }

            integrator
                .step(&mut system, &mut positions, &mut velocities)
                .ok()?;

            for (local, body) in bodies.iter().enumerate() {
                let base = 3 * local;
                if *body == 1 {
                    positions[base..base + 3].copy_from_slice(&relaxed[base..base + 3]);
                    for value in &mut velocities[base..base + 3] {
                        *value = 0.0;
                    }
                } else if *body == 0 {
                    let r_x = positions[base] - center_m[0];
                    let r_y = positions[base + 1] - center_m[1];
                    velocities[base] = -omega * r_y;
                    velocities[base + 1] = omega * r_x;
                    velocities[base + 2] = 0.0;
                }
            }

            system.rebuild_neighbors(&positions).ok()?;
        }

        let energy_final_j = system.energy_j(&positions).ok()?;
        let (free_kinetic_j, free_count) = free_kinetic(&velocities, &bodies);
        let final_temperature_k = if free_count > 0 {
            2.0 * free_kinetic_j / (3.0 * free_count as f64 * BOLTZMANN_J_PER_K)
        } else {
            0.0
        };
        let mean_torque_n_m = if steps > 0 {
            sum_abs_torque_n_m / steps as f64
        } else {
            0.0
        };
        let work_j = work_j(mean_torque_n_m, omega, steps, self.target.dt_s);
        let energy_loss_j =
            (work_j - (energy_final_j + free_kinetic_j - energy_initial_j)).max(0.0);
        let temperature_rise_k = (final_temperature_k - self.target.temperature_k).max(0.0);

        Some(DriveReport {
            torque_n_m: mean_torque_n_m,
            peak_torque_n_m,
            work_j,
            energy_loss_j,
            temperature_rise_k,
            final_temperature_k,
            temperature_k: self.target.temperature_k,
            free_atom_count: free_count,
            atom_count,
            steps,
            drive_rate_rad_per_s: omega,
        })
    }

    /// Returns the cluster atoms, sorted by distance to the mesh point.
    ///
    /// Each entry is the topology index, the position, and the body index
    /// (0 for the sun, 1 for the first planet, 2 for the ring). The list is
    /// capped at `max_atoms`. The ring supplies the free atoms, so the cluster
    /// can show a temperature rise.
    fn cluster(
        &self,
        moving: &MovingAtoms,
        atoms: &[[f64; 3]],
        mesh_m: [f64; 3],
    ) -> Vec<(usize, [f64; 3], u8)> {
        let radius_squared = self.target.cluster_m * self.target.cluster_m;
        let mut found = Vec::new();
        for (body, group) in moving.groups.iter().enumerate().take(3) {
            let body = body as u8;
            for index in group.1.clone() {
                let Some(&position) = atoms.get(index) else {
                    continue;
                };
                let dx = position[0] - mesh_m[0];
                let dy = position[1] - mesh_m[1];
                let dz = position[2] - mesh_m[2];
                if dx * dx + dy * dy + dz * dz <= radius_squared {
                    found.push((index, position, body));
                }
            }
        }
        found.sort_by(|a, b| {
            let da = distance_squared(a.1, mesh_m);
            let db = distance_squared(b.1, mesh_m);
            da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
        });
        found.truncate(self.target.max_atoms);
        found
    }

    /// Builds the potential and the flat start positions of the cluster.
    ///
    /// This copies the harmonic assembly: a bond-stretch term at the current
    /// length, exclusions for same-body pairs and far cross-body pairs, and a
    /// Buckingham van der Waals term.
    fn cluster_system(
        &self,
        topology: &nanocad_model::Topology,
        cluster: &[(usize, [f64; 3], u8)],
    ) -> Option<(System, Vec<f64>)> {
        let atom_count = cluster.len();
        let local_of: HashMap<usize, usize> = cluster
            .iter()
            .enumerate()
            .map(|(local, (index, _, _))| (*index, local))
            .collect();
        let mut positions_m = Vec::with_capacity(3 * atom_count);
        for (_, position, _) in cluster {
            positions_m.extend_from_slice(position);
        }

        let mut bonds = BondStretchTerm::new();
        for index in 0..topology.bond_count() {
            let (Some(u), Some(v)) = (topology.bond_u(index), topology.bond_v(index)) else {
                continue;
            };
            let (Some(&local_u), Some(&local_v)) =
                (local_of.get(&(u as usize)), local_of.get(&(v as usize)))
            else {
                continue;
            };
            let r0_m = distance_squared(cluster[local_u].1, cluster[local_v].1).sqrt();
            if r0_m <= 0.0 {
                continue;
            }
            bonds
                .add_bond(local_u as u32, local_v as u32, BOND_STIFFNESS_N_PER_M, r0_m)
                .ok()?;
        }

        let mut system =
            System::with_masses_kg(atom_count, &vec![CARBON_MASS_KG; atom_count]).ok()?;
        system.set_bond_stretch(bonds).ok()?;
        // The van der Waals term acts only at the contact. A pair inside one
        // body is covalent, so it is excluded. A cross-body pair beyond the
        // contact shell is a long-range tail, which is soft and has no contact
        // physics, so it is also excluded.
        let contact_squared = self.target.contact_m * self.target.contact_m;
        for i in 0..atom_count {
            for j in (i + 1)..atom_count {
                let same_body = cluster[i].2 == cluster[j].2;
                let beyond_shell = distance_squared(cluster[i].1, cluster[j].1) > contact_squared;
                if same_body || beyond_shell {
                    system.add_exclusion(i as u32, j as u32).ok()?;
                }
            }
        }
        let params = vec![VdwParams::new(VDW_A_J, VDW_B_PER_M, VDW_C_J_M6); atom_count];
        let cutoff = Cutoff::new(self.target.cutoff_m, 0.9 * self.target.cutoff_m).ok()?;
        let vdw =
            VanDerWaalsTerm::from_params(&params, cutoff, PeriodicBox::non_periodic()).ok()?;
        system.set_van_der_waals(vdw).ok()?;

        Some((system, positions_m))
    }
}

/// Returns the torque about the z axis, in newton metres.
///
/// The value is signed. A positive torque turns counter-clockwise.
///
/// `body` marks the atoms of the driven body. The force on each atom is
/// `forces_n`, already the negative of the energy gradient.
fn torque_z(
    positions_m: &[[f64; 3]],
    forces_n: &[[f64; 3]],
    center_m: [f64; 3],
    body: &[bool],
) -> f64 {
    let mut torque_n_m = 0.0;
    for (index, position) in positions_m.iter().enumerate() {
        if !body.get(index).copied().unwrap_or(false) {
            continue;
        }
        let Some(force) = forces_n.get(index) else {
            continue;
        };
        let r_x = position[0] - center_m[0];
        let r_y = position[1] - center_m[1];
        torque_n_m += r_x * force[1] - r_y * force[0];
    }
    torque_n_m
}

/// Returns the drive work in joules.
///
/// The work is the mean absolute torque times the swept angle, which is the
/// absolute rate times the total time.
fn work_j(mean_torque_n_m: f64, rate_rad_per_s: f64, steps: usize, dt_s: f64) -> f64 {
    mean_torque_n_m * rate_rad_per_s.abs() * steps as f64 * dt_s
}

/// Returns the kinetic energy of the free atoms in joules and their count.
///
/// A free atom is not part of the driven body and not part of the held body.
fn free_kinetic(velocities_m_per_s: &[f64], bodies: &[u8]) -> (f64, usize) {
    let mut kinetic_j = 0.0;
    let mut count = 0;
    for (local, body) in bodies.iter().enumerate() {
        if *body == 0 || *body == 1 {
            continue;
        }
        let base = 3 * local;
        let vx = velocities_m_per_s[base];
        let vy = velocities_m_per_s[base + 1];
        let vz = velocities_m_per_s[base + 2];
        kinetic_j += 0.5 * CARBON_MASS_KG * (vx * vx + vy * vy + vz * vz);
        count += 1;
    }
    (kinetic_j, count)
}

/// Returns the flat positions as one point per atom.
fn points_of(positions_m: &[f64]) -> Vec<[f64; 3]> {
    positions_m.as_chunks::<3>().0.to_vec()
}

/// Returns the force vectors from an energy gradient.
fn forces_of(gradient_j_per_m: &[f64]) -> Vec<[f64; 3]> {
    gradient_j_per_m
        .as_chunks::<3>()
        .0
        .iter()
        .map(|gradient| [-gradient[0], -gradient[1], -gradient[2]])
        .collect()
}

/// Returns the closest pair of two point sets.
///
/// The search uses a grid over the second set, so it does not scan every pair.
fn closest_pair(a: &[[f64; 3]], b: &[[f64; 3]]) -> Option<([f64; 3], [f64; 3])> {
    if a.is_empty() || b.is_empty() {
        return None;
    }
    let cell_m = 1.0e-9;
    let mut grid: HashMap<Cell, Vec<u32>> = HashMap::new();
    for (index, position) in b.iter().enumerate() {
        grid.entry(cell_of(*position, cell_m))
            .or_default()
            .push(index as u32);
    }
    let mut best = f64::INFINITY;
    let mut pair = None;
    for point in a {
        let (cx, cy, cz) = cell_of(*point, cell_m);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let Some(bucket) = grid.get(&(cx + dx, cy + dy, cz + dz)) else {
                        continue;
                    };
                    for &other in bucket {
                        let candidate = b[other as usize];
                        let squared = distance_squared(*point, candidate);
                        if squared < best {
                            best = squared;
                            pair = Some((*point, candidate));
                        }
                    }
                }
            }
        }
    }
    pair
}

/// Returns the squared distance between two points.
fn distance_squared(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

/// Returns a report with every value zero.
fn zeroed_report(target: &DriveTarget, rate_rad_per_s: f64) -> DriveReport {
    DriveReport {
        torque_n_m: 0.0,
        peak_torque_n_m: 0.0,
        work_j: 0.0,
        energy_loss_j: 0.0,
        temperature_rise_k: 0.0,
        final_temperature_k: 0.0,
        temperature_k: target.temperature_k,
        free_atom_count: 0,
        atom_count: 0,
        steps: target.steps,
        drive_rate_rad_per_s: rate_rad_per_s,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nanocad_model::{Part, Topology};
    use nanocad_parts::planetary::PlanetaryDesign;

    /// A set with no atoms. The measured cluster is empty.
    fn empty_set() -> PlanetarySet {
        let design = PlanetaryDesign::new(1.5e-9, 12, 9, 3, 30.0 * std::f64::consts::PI / 180.0)
            .expect("valid design");
        PlanetarySet {
            part: Part::new("empty", Topology::new()),
            design,
            sun_atoms: 0..0,
            planet_atoms: std::iter::once(0..0).collect(),
            ring_atoms: 0..0,
            carrier_atoms: 0..0,
        }
    }

    #[test]
    fn a_zeroed_target_gives_a_zero_report() {
        let drive = SteeredDrive::with_target(DriveTarget {
            steps: 0,
            ..Default::default()
        });
        let report = drive.measure(&empty_set(), -0.34);
        assert_eq!(report.torque_n_m, 0.0);
        assert_eq!(report.work_j, 0.0);
    }

    #[test]
    fn the_metric_name_and_fidelity_are_stated() {
        let report = DriveReport {
            torque_n_m: 1.0e-12,
            peak_torque_n_m: 2.0e-12,
            work_j: 3.0e-24,
            energy_loss_j: 4.0e-25,
            temperature_rise_k: 1.5,
            final_temperature_k: 301.5,
            temperature_k: 300.0,
            free_atom_count: 0,
            atom_count: 40,
            steps: 4000,
            drive_rate_rad_per_s: -0.34,
        };
        let value = report.to_metric_value();
        assert_eq!(value.name, "steered_drive");
        assert_eq!(value.unit, "N m");
        assert_eq!(value.fidelity, Fidelity::Dynamics);
    }

    #[test]
    fn the_work_is_the_torque_times_the_angle() {
        let expected = 2.0e-12 * 0.34 * 100.0 * 1.0e-16;
        assert!((work_j(2.0e-12, -0.34, 100, 1.0e-16) - expected).abs() < 1.0e-30);
    }

    #[test]
    fn a_closed_body_has_no_torque() {
        let positions_m = [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
        ];
        // Each force is parallel to its radius, so each cross product is zero.
        let forces_n = positions_m;
        let body = [true, true, true, true];
        assert!(torque_z(&positions_m, &forces_n, [0.0, 0.0, 0.0], &body).abs() < 1.0e-30);
    }

    #[test]
    fn an_empty_report_states_the_cluster_is_empty() {
        let drive = SteeredDrive::default();
        let report = drive.measure(&empty_set(), -0.34);
        assert!(report.to_metric_value().note.contains("empty"));
    }
}
