//! The harmonic metric: the normal modes of the atoms at a gear mesh.
//!
//! The Hessian of the potential gives the normal modes. A soft mode at the
//! mesh means the surfaces can slide with little force. A stiff mode means the
//! teeth hold each other. The metric builds a small cluster around the closest
//! sun and planet atoms, expands the potential to second order, and reports
//! the frequency of the softest internal mode.
//!
//! The cluster is bounded, so the mode is that of the cluster, not of the
//! whole gear. The potential has a bond-stretch term for every topology bond
//! inside the cluster and a Buckingham van der Waals term for every
//! non-bonded pair. The bond rest length is the current length, so the
//! expansion is about the generated geometry.
//!
//! [`HarmonicMesh::measure_part`] measures a whole part instead of a cut-out
//! cluster. A whole part has no dangling crystal boundary, so its torsion and
//! out-of-plane terms are physical and the potential includes them.

use std::collections::HashMap;

use nanocad_engine::{
    hessian_finite_difference, mass_weighted, minimize_with, spectrum, symmetric_eigenvalues,
    BondStretchTerm, Cutoff, MinimizeMethod, MinimizeOptions, OutOfPlaneTerm, PeriodicBox,
    Spectrum, System, TorsionTerm, VanDerWaalsTerm, VdwParams,
};
use nanocad_parts::planetary::PlanetarySet;

use crate::bonded::{add_impropers, add_torsions};
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

/// The separation below which a non-bonded pair counts as covalent, in metres.
///
/// The diamond second and third shells sit below this distance. A whole part
/// is one covalent network, so its van der Waals term must not act on those
/// pairs. Without this exclusion the minimizer lands two atoms on the same
/// point and fails with `CoincidentNonbondedAtoms`.
const COVALENT_EXCLUSION_M: f64 = 3.1e-10;

/// The speed of light in centimetres per second, for the wavenumber.
const LIGHT_CM_PER_S: f64 = 2.997_924_58e10;

/// The settings of a harmonic measurement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicTarget {
    /// The radius of the atom cluster about the mesh point, in metres.
    pub cluster_m: f64,
    /// The largest atom count of the cluster.
    pub max_atoms: usize,
    /// The van der Waals cutoff in metres.
    pub cutoff_m: f64,
    /// The contact shell in metres. A cross-body pair farther apart than this
    /// does not interact, so only the contact carries the van der Waals term.
    pub contact_m: f64,
    /// The finite-difference step of the Hessian, in metres.
    pub difference_step_m: f64,
}

impl Default for HarmonicTarget {
    fn default() -> Self {
        Self {
            cluster_m: 8.0e-10,
            max_atoms: 180,
            cutoff_m: 1.2e-9,
            contact_m: 3.5e-10,
            difference_step_m: 1.0e-13,
        }
    }
}

/// The result of a harmonic measurement.
#[derive(Clone, Debug, PartialEq)]
pub struct HarmonicReport {
    /// The frequency of the softest internal mode, in hertz.
    pub lowest_hz: f64,
    /// The mass-weighted eigenvalue of the softest internal mode.
    pub lowest_eigenvalue_per_s2: f64,
    /// The number of unstable modes below the tolerance.
    pub unstable_count: usize,
    /// The number of atoms in the cluster.
    pub atom_count: usize,
    /// The number of topology bonds in the cluster.
    pub bond_count: usize,
    /// The contact shell of the van der Waals term, in metres.
    pub contact_m: f64,
}

impl HarmonicReport {
    /// Returns the wavenumber of the softest internal mode, in reciprocal
    /// centimetres.
    pub fn lowest_wavenumber_per_cm(&self) -> f64 {
        self.lowest_hz / LIGHT_CM_PER_S
    }

    /// Converts the report to a metric value.
    pub fn to_metric_value(&self) -> MetricValue {
        let wavenumber = self.lowest_wavenumber_per_cm();
        let verdict = if self.atom_count == 0 {
            "the cluster is empty or the part is too large for the atom cap, so the mesh has no measured mode"
        } else if self.lowest_hz < 1.0e11 {
            "the mesh has a soft mode, so the teeth can slide"
        } else if self.lowest_hz < 1.0e12 {
            "the mesh has a mode of moderate stiffness"
        } else {
            "the mesh is stiff"
        };
        MetricValue {
            name: "mesh_mode".to_string(),
            value: self.lowest_hz,
            unit: "Hz".to_string(),
            fidelity: Fidelity::Harmonic,
            note: format!(
                "the softest internal mode of a {}-atom cluster at the sun-planet mesh is {:.3e} Hz ({:.1} per cm) over {} bonds; {verdict}; {} unstable modes; the cluster is relaxed to a local minimum first; the model is bond stretch plus a Buckingham van der Waals term inside a {:.1e} m contact shell",
                self.atom_count,
                self.lowest_hz,
                wavenumber,
                self.bond_count,
                self.unstable_count,
                self.contact_m
            ),
        }
    }
}

/// The harmonic mesh metric.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HarmonicMesh {
    /// The settings of the measurement.
    pub target: HarmonicTarget,
}

impl HarmonicMesh {
    /// A harmonic metric with the given settings.
    pub fn new(target: HarmonicTarget) -> Self {
        Self { target }
    }

    /// Measures the softest internal mode at the first sun-planet mesh.
    ///
    /// Returns an empty report when the cluster is too small or the potential
    /// does not build.
    pub fn measure(&self, set: &PlanetarySet) -> HarmonicReport {
        let topology = &set.part.topology;
        let Some(planet_range) = set.planet_atoms.first() else {
            return empty_report();
        };
        let sun_range = set.sun_atoms.clone();
        let sun_positions: Vec<[f64; 3]> = sun_range
            .clone()
            .filter_map(|index| topology.position_m(index))
            .collect();
        let planet_positions: Vec<[f64; 3]> = planet_range
            .clone()
            .filter_map(|index| topology.position_m(index))
            .collect();

        let Some((sun_point, planet_point)) = closest_pair(&sun_positions, &planet_positions)
        else {
            return empty_report();
        };
        let midpoint = [
            0.5 * (sun_point[0] + planet_point[0]),
            0.5 * (sun_point[1] + planet_point[1]),
            0.5 * (sun_point[2] + planet_point[2]),
        ];

        let cluster = self.cluster(&sun_range, planet_range, topology, midpoint);
        if cluster.len() < 3 {
            return empty_report();
        }
        let Some(report) = self.spectrum_of_cluster(topology, &cluster) else {
            return empty_report();
        };
        report
    }

    /// Measures the softest internal mode of a whole part.
    ///
    /// Every topology atom enters the potential in topology index order, with
    /// the carbon-like mass of the cluster path. The potential has bond
    /// stretch, torsion, out-of-plane and Buckingham van der Waals terms. The
    /// system is relaxed to a local minimum first, then the Hessian is
    /// finite-differenced at the relaxed geometry.
    ///
    /// Returns an empty report when the part is empty or larger than
    /// `target.max_atoms`, or when the potential does not build. When the part
    /// is too large, the caller must generate a smaller part.
    pub fn measure_part(&self, part: &nanocad_model::Part) -> HarmonicReport {
        let topology = &part.topology;
        let atom_count = topology.atom_count();
        if atom_count == 0 || atom_count > self.target.max_atoms {
            return empty_report();
        }
        let Some((found, _torsion_count, _improper_count)) = self.spectrum_of_part(topology) else {
            return empty_report();
        };

        // The six rigid modes are zero. The softest real mode is the first
        // frequency above the floor.
        let floor_hz = 1.0e9;
        let index = found
            .frequencies_hz
            .iter()
            .position(|frequency| *frequency > floor_hz)
            .unwrap_or(found.frequencies_hz.len() - 1);
        HarmonicReport {
            lowest_hz: found.frequencies_hz[index],
            lowest_eigenvalue_per_s2: found.eigenvalues[index],
            unstable_count: found.unstable_count,
            atom_count,
            bond_count: topology.bond_count(),
            contact_m: self.target.contact_m,
        }
    }

    /// Returns the cluster atoms, sorted by distance to the mesh point.
    ///
    /// Each entry is the topology index and the position. The list is capped
    /// at `max_atoms`.
    fn cluster(
        &self,
        sun_range: &std::ops::Range<usize>,
        planet_range: &std::ops::Range<usize>,
        topology: &nanocad_model::Topology,
        midpoint: [f64; 3],
    ) -> Vec<(usize, [f64; 3], u8)> {
        let radius_squared = self.target.cluster_m * self.target.cluster_m;
        let mut found = Vec::new();
        for (body, range) in [(0u8, sun_range.clone()), (1u8, planet_range.clone())] {
            for index in range {
                let Some(position) = topology.position_m(index) else {
                    continue;
                };
                let dx = position[0] - midpoint[0];
                let dy = position[1] - midpoint[1];
                let dz = position[2] - midpoint[2];
                if dx * dx + dy * dy + dz * dz <= radius_squared {
                    found.push((index, position, body));
                }
            }
        }
        found.sort_by(|a, b| {
            let da = distance_squared(a.1, midpoint);
            let db = distance_squared(b.1, midpoint);
            da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
        });
        found.truncate(self.target.max_atoms);
        found
    }

    /// Builds the potential of the cluster and returns its spectrum.
    fn spectrum_of_cluster(
        &self,
        topology: &nanocad_model::Topology,
        cluster: &[(usize, [f64; 3], u8)],
    ) -> Option<HarmonicReport> {
        let atom_count = cluster.len();
        let local_of: HashMap<usize, usize> = cluster
            .iter()
            .enumerate()
            .map(|(local, (topology_index, _, _))| (*topology_index, local))
            .collect();
        let mut positions_m = Vec::with_capacity(3 * atom_count);
        for (_, position, _) in cluster {
            positions_m.extend_from_slice(position);
        }

        let mut bonds = BondStretchTerm::new();
        let mut bond_count = 0;
        for index in 0..topology.bond_count() {
            let (Some(u), Some(v)) = (topology.bond_u(index), topology.bond_v(index)) else {
                continue;
            };
            let (Some(&local_u), Some(&local_v)) =
                (local_of.get(&(u as usize)), local_of.get(&(v as usize)))
            else {
                continue;
            };
            let a = cluster[local_u].1;
            let b = cluster[local_v].1;
            let r0_m = distance_squared(a, b).sqrt();
            if r0_m <= 0.0 {
                continue;
            }
            bonds
                .add_bond(local_u as u32, local_v as u32, BOND_STIFFNESS_N_PER_M, r0_m)
                .ok()?;
            bond_count += 1;
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

        let mut relaxed_m = positions_m.clone();
        let options = MinimizeOptions {
            max_iterations: 400,
            gradient_tolerance_n: 1.0e-10,
            initial_step_m: 1.0e-12,
        };
        minimize_with(
            &mut system,
            &mut relaxed_m,
            &options,
            MinimizeMethod::ConjugateGradientThenLbfgs,
        )
        .ok()?;

        let hessian =
            hessian_finite_difference(&mut system, &relaxed_m, self.target.difference_step_m)
                .ok()?;
        let masses = vec![CARBON_MASS_KG; atom_count];
        let weighted = mass_weighted(&hessian, &masses).ok()?;
        let eigenvalues = symmetric_eigenvalues(&weighted, 3 * atom_count).ok()?;
        let found = spectrum(eigenvalues);

        // The six rigid modes are zero, and an atom with no bond and no
        // contact adds three more zero modes. The softest real mode is the
        // first above the floor.
        let floor_hz = 1.0e9;
        let index = found
            .frequencies_hz
            .iter()
            .position(|frequency| *frequency > floor_hz)
            .unwrap_or(found.frequencies_hz.len() - 1);
        Some(HarmonicReport {
            lowest_hz: found.frequencies_hz[index],
            lowest_eigenvalue_per_s2: found.eigenvalues[index],
            unstable_count: found.unstable_count,
            atom_count,
            bond_count,
            contact_m: self.target.contact_m,
        })
    }

    /// Builds the potential of a whole part and returns its spectrum with the
    /// torsion and out-of-plane counts.
    ///
    /// The van der Waals term excludes every bonded pair and every pair inside
    /// the covalent shell, because a whole part has no separate bodies. The
    /// torsion and out-of-plane terms are added here: a whole part has no
    /// under-coordinated crystal boundary, so those terms are physical. A
    /// cut-out cluster cannot use them because its boundary carbons have too
    /// few neighbours.
    fn spectrum_of_part(
        &self,
        topology: &nanocad_model::Topology,
    ) -> Option<(Spectrum, usize, usize)> {
        let atom_count = topology.atom_count();
        let mut positions_m = Vec::with_capacity(atom_count);
        for index in 0..atom_count {
            let (Some(position_m), Some(_element)) =
                (topology.position_m(index), topology.element(index))
            else {
                return None;
            };
            positions_m.push(position_m);
        }
        let flat: Vec<f64> = positions_m.iter().flat_map(|p| p.iter().copied()).collect();

        let mut bonds = BondStretchTerm::new();
        let mut neighbours: Vec<Vec<u32>> = vec![Vec::new(); atom_count];
        for index in 0..topology.bond_count() {
            let (Some(u), Some(v)) = (topology.bond_u(index), topology.bond_v(index)) else {
                continue;
            };
            let (u, v) = (u as usize, v as usize);
            if u >= atom_count || v >= atom_count {
                continue;
            }
            let r0_m = distance_squared(positions_m[u], positions_m[v]).sqrt();
            if r0_m <= 0.0 {
                continue;
            }
            bonds
                .add_bond(u as u32, v as u32, BOND_STIFFNESS_N_PER_M, r0_m)
                .ok()?;
            neighbours[u].push(v as u32);
            neighbours[v].push(u as u32);
        }

        let mut system =
            System::with_masses_kg(atom_count, &vec![CARBON_MASS_KG; atom_count]).ok()?;
        system.set_bond_stretch(bonds).ok()?;

        // A whole part has no bodies, so the van der Waals term excludes every
        // bonded pair and every covalent-shell pair by hand. The remaining
        // pairs interact only inside the contact shell.
        let contact_squared = self.target.contact_m * self.target.contact_m;
        let covalent_squared = COVALENT_EXCLUSION_M * COVALENT_EXCLUSION_M;
        let mut bonded_pair = vec![false; atom_count * atom_count];
        for (u, list) in neighbours.iter().enumerate() {
            for &v in list {
                bonded_pair[u * atom_count + v as usize] = true;
            }
        }
        for i in 0..atom_count {
            for j in (i + 1)..atom_count {
                let distance = distance_squared(positions_m[i], positions_m[j]);
                let beyond_shell = distance > contact_squared;
                let covalent = distance <= covalent_squared;
                if bonded_pair[i * atom_count + j] || beyond_shell || covalent {
                    system.add_exclusion(i as u32, j as u32).ok()?;
                }
            }
        }
        let params = vec![VdwParams::new(VDW_A_J, VDW_B_PER_M, VDW_C_J_M6); atom_count];
        let cutoff = Cutoff::new(self.target.cutoff_m, 0.9 * self.target.cutoff_m).ok()?;
        let vdw =
            VanDerWaalsTerm::from_params(&params, cutoff, PeriodicBox::non_periodic()).ok()?;
        system.set_van_der_waals(vdw).ok()?;

        let mut torsions = TorsionTerm::new();
        let torsion_count = add_torsions(&mut torsions, &neighbours).ok()?;
        if torsion_count > 0 {
            system.set_torsion(torsions).ok()?;
        }
        let mut impropers = OutOfPlaneTerm::new();
        let improper_count = add_impropers(&mut impropers, &neighbours, &positions_m).ok()?;
        if improper_count > 0 {
            system.set_out_of_plane(impropers).ok()?;
        }

        let mut relaxed_m = flat.clone();
        let options = MinimizeOptions {
            max_iterations: 400,
            gradient_tolerance_n: 1.0e-10,
            initial_step_m: 1.0e-12,
        };
        minimize_with(
            &mut system,
            &mut relaxed_m,
            &options,
            MinimizeMethod::ConjugateGradientThenLbfgs,
        )
        .ok()?;

        // The minimizer moves the atoms, so the vdW pair list must match the
        // relaxed geometry before the Hessian reads it.
        system.rebuild_neighbors(&relaxed_m).ok()?;

        let hessian =
            hessian_finite_difference(&mut system, &relaxed_m, self.target.difference_step_m)
                .ok()?;
        let masses = vec![CARBON_MASS_KG; atom_count];
        let weighted = mass_weighted(&hessian, &masses).ok()?;
        let eigenvalues = symmetric_eigenvalues(&weighted, 3 * atom_count).ok()?;
        Some((spectrum(eigenvalues), torsion_count, improper_count))
    }
}

/// Returns the closest sun-planet atom pair.
///
/// The search uses a grid over the planet atoms, so it does not scan every
/// pair.
fn closest_pair(
    sun_positions: &[[f64; 3]],
    planet_positions: &[[f64; 3]],
) -> Option<([f64; 3], [f64; 3])> {
    if sun_positions.is_empty() || planet_positions.is_empty() {
        return None;
    }
    let cell_m = 1.0e-9;
    let mut grid: HashMap<Cell, Vec<u32>> = HashMap::new();
    for (index, position) in planet_positions.iter().enumerate() {
        grid.entry(cell_of(*position, cell_m))
            .or_default()
            .push(index as u32);
    }
    let mut best = f64::INFINITY;
    let mut pair = None;
    for sun in sun_positions {
        let (cx, cy, cz) = cell_of(*sun, cell_m);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let Some(bucket) = grid.get(&(cx + dx, cy + dy, cz + dz)) else {
                        continue;
                    };
                    for &other in bucket {
                        let planet = planet_positions[other as usize];
                        let squared = distance_squared(*sun, planet);
                        if squared < best {
                            best = squared;
                            pair = Some((*sun, planet));
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

/// Returns an empty report.
fn empty_report() -> HarmonicReport {
    HarmonicReport {
        lowest_hz: 0.0,
        lowest_eigenvalue_per_s2: 0.0,
        unstable_count: 0,
        atom_count: 0,
        bond_count: 0,
        contact_m: 0.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nanocad_parts::{DiamondGenerator, ParameterSet, PartGenerator, PlateGenerator};

    /// Builds a small uncapped diamond block.
    ///
    /// A tiny part keeps the test fast: the diagonalisation is cubic in the
    /// coordinate count and the test profile is unoptimised.
    fn small_diamond(cells: (f64, f64, f64)) -> nanocad_model::Part {
        DiamondGenerator
            .generate(
                &ParameterSet::new()
                    .with("cells_x", cells.0)
                    .with("cells_y", cells.1)
                    .with("cells_z", cells.2),
            )
            .expect("valid diamond")
    }

    /// A small whole part has a spectrum, and the metric finds a soft one.
    ///
    /// Diamond `2x1x1` gives 16 atoms. Two of them carry no bond, so the
    /// potential has more zero modes than the six rigid ones. The measured
    /// count is two.
    #[test]
    fn a_small_part_has_a_spectrum() {
        let part = small_diamond((2.0, 1.0, 1.0));
        let target = HarmonicTarget {
            max_atoms: 64,
            ..HarmonicTarget::default()
        };
        let report = HarmonicMesh::new(target).measure_part(&part);
        assert_eq!(report.atom_count, 16);
        assert!(
            report.lowest_hz > 0.0,
            "the softest mode is not positive: {} Hz",
            report.lowest_hz
        );
        assert_eq!(report.unstable_count, 2);
    }

    /// A part above the atom cap is refused with a "too large" note.
    #[test]
    fn a_part_above_the_cap_is_refused() {
        let part = PlateGenerator
            .generate(
                &ParameterSet::new()
                    .with("size_x_m", 1.5e-9)
                    .with("size_y_m", 1.5e-9)
                    .with("thickness_m", 7.0e-10),
            )
            .expect("valid plate");
        let target = HarmonicTarget {
            max_atoms: 4,
            ..HarmonicTarget::default()
        };
        let report = HarmonicMesh::new(target).measure_part(&part);
        assert_eq!(report.atom_count, 0);
        assert!(
            report.to_metric_value().note.contains("too large"),
            "the note does not say the part is too large"
        );
    }

    /// The whole-part report carries the harmonic `mesh_mode` identity.
    #[test]
    fn the_whole_part_mode_is_at_the_harmonic_fidelity() {
        let part = small_diamond((2.0, 1.0, 1.0));
        let target = HarmonicTarget {
            max_atoms: 64,
            ..HarmonicTarget::default()
        };
        let value = HarmonicMesh::new(target)
            .measure_part(&part)
            .to_metric_value();
        assert_eq!(value.fidelity, Fidelity::Harmonic);
        assert_eq!(value.name, "mesh_mode");
    }

    /// A larger whole part returns a non-empty report.
    #[test]
    fn a_bigger_part_does_not_crash() {
        let part = small_diamond((2.0, 2.0, 1.0));
        let target = HarmonicTarget {
            max_atoms: 64,
            ..HarmonicTarget::default()
        };
        let report = HarmonicMesh::new(target).measure_part(&part);
        assert!(report.atom_count > 0, "the part returned an empty report");
        assert!(report.lowest_hz > 0.0);
    }

    /// A pair of atoms has a known mode, and the metric finds a soft one.
    #[test]
    fn a_small_cluster_has_a_spectrum() {
        let mut system = System::with_masses_kg(2, &[CARBON_MASS_KG; 2]).unwrap();
        let mut bonds = BondStretchTerm::new();
        bonds
            .add_bond(0, 1, BOND_STIFFNESS_N_PER_M, 1.5e-10)
            .unwrap();
        system.set_bond_stretch(bonds).unwrap();
        let positions_m = vec![0.0, 0.0, 0.0, 1.5e-10, 0.0, 0.0];
        let hessian = hessian_finite_difference(&mut system, &positions_m, 1.0e-13).unwrap();
        let weighted = mass_weighted(&hessian, &[CARBON_MASS_KG; 2]).unwrap();
        let found = spectrum(symmetric_eigenvalues(&weighted, 6).unwrap());
        let expected =
            (2.0 * BOND_STIFFNESS_N_PER_M / CARBON_MASS_KG).sqrt() / std::f64::consts::TAU;
        assert!((found.frequencies_hz[5] - expected).abs() / expected < 1.0e-6);
    }

    /// An empty report carries no mode.
    #[test]
    fn an_empty_report_is_zero() {
        let report = empty_report();
        assert_eq!(report.atom_count, 0);
        assert_eq!(report.lowest_hz, 0.0);
        assert!(report.to_metric_value().note.contains("empty"));
    }

    /// The wavenumber is the frequency over the speed of light.
    #[test]
    fn the_wavenumber_divides_by_the_speed_of_light() {
        let report = HarmonicReport {
            lowest_hz: 9.0e12,
            lowest_eigenvalue_per_s2: 0.0,
            unstable_count: 0,
            atom_count: 10,
            bond_count: 12,
            contact_m: 3.5e-10,
        };
        assert!((report.lowest_wavenumber_per_cm() - 300.207_5).abs() < 0.01);
    }
}
