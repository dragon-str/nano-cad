//! The binding energy of one guest in one pocket.
//!
//! The pocket is held rigid and the guest is rigid. The energy is the sum of
//! the Lennard-Jones and Coulomb interaction over every guest and wall atom
//! pair inside a cutoff. The search moves the guest along the pocket axis and
//! turns it about that axis, and it reports the lowest energy that it finds.
//!
//! This is a model of a model. The guest charges, the wall charges and the
//! nonbonded parameters each come from a stated source, and none of them is
//! validated against experiment. The search is coarse: it does not relax the
//! wall, it does not relax the guest, and it does not model a solvent. Two
//! values are comparable only when they come from the same target. See
//! ADR-0063, ADR-0064 and ADR-0066.

use nanocad_model::{pair_params, Element, Part};
use nanocad_parts::Guest;

use crate::{Fidelity, MetricValue};

/// The Coulomb constant, in newton square metres per square coulomb.
const COULOMB_CONSTANT_N_M2_PER_C2: f64 = 8.987_551_792_3e9;

/// The Boltzmann constant, in joules per kelvin.
const BOLTZMANN_J_PER_K: f64 = 1.380_649e-23;

/// The search settings of one binding measurement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BindingTarget {
    /// The pair cutoff, in metres.
    pub cutoff_m: f64,
    /// The number of positions along the pocket axis.
    pub displacements: usize,
    /// The travel along the pocket axis, in metres.
    pub span_m: f64,
    /// The number of turns about the pocket axis.
    pub azimuths: usize,
    /// The temperature of the reported thermal energy, in kelvin.
    pub temperature_k: f64,
    /// The closest allowed pair distance, as a fraction of sigma.
    ///
    /// A position where a guest atom and a wall atom approach closer than this
    /// is not a physical bound state, so the search rejects it. Without this
    /// guard the rigid search finds interpenetrating positions and reports a
    /// meaningless energy.
    pub overlap_ratio: f64,
}

impl Default for BindingTarget {
    fn default() -> Self {
        Self {
            cutoff_m: 1.0e-9,
            displacements: 15,
            span_m: 2.0e-10,
            azimuths: 12,
            temperature_k: 300.0,
            overlap_ratio: 0.85,
        }
    }
}

/// The result of one binding measurement.
#[derive(Clone, Debug, PartialEq)]
pub struct BindingReport {
    /// The guest id.
    pub guest: String,
    /// The lowest interaction energy that the search found, in joules.
    pub minimum_energy_j: f64,
    /// The interaction energy with the guest at the pocket centre, in joules.
    pub centre_energy_j: f64,
    /// The displacement along the pocket axis at the minimum, in metres.
    pub displacement_m: f64,
    /// The turn about the pocket axis at the minimum, in radians.
    pub azimuth_rad: f64,
    /// The number of wall atoms inside the cutoff at the minimum.
    pub wall_atom_count: usize,
    /// The number of guest atoms.
    pub guest_atom_count: usize,
    /// The number of search positions that the overlap guard rejected.
    pub blocked_positions: usize,
    /// The temperature used for the thermal ratio, in kelvin.
    pub temperature_k: f64,
}

impl BindingReport {
    /// The thermal energy at the reported temperature, in joules.
    pub fn thermal_energy_j(&self) -> f64 {
        BOLTZMANN_J_PER_K * self.temperature_k
    }

    /// The binding energy in units of the thermal energy.
    ///
    /// A positive value means that the guest is held. A value below about one
    /// means that thermal motion releases it.
    pub fn over_kt(&self) -> f64 {
        -self.minimum_energy_j / self.thermal_energy_j()
    }

    /// The value as a metric, at quasi-static fidelity.
    pub fn to_metric_value(&self) -> MetricValue {
        MetricValue {
            name: format!("binding_{}", self.guest),
            value: self.minimum_energy_j,
            unit: "J".to_string(),
            fidelity: Fidelity::QuasiStatic,
            note: format!(
                "{} in a rigid pocket, {:.2} kT at {:.0} K",
                self.guest,
                self.over_kt(),
                self.temperature_k
            ),
        }
    }
}

/// Measure the binding energy of `guest` in `pocket`.
///
/// `centre_m` is the point that the guest centroid should occupy, in metres.
/// The search turns the guest about the `z` axis and moves it along `z`, and
/// it reports the lowest energy that it finds.
pub fn binding_energy_j(
    pocket: &Part,
    centre_m: [f64; 3],
    guest: &Guest,
    target: &BindingTarget,
) -> BindingReport {
    let centroid_m = guest.centroid_m();
    let wall = pocket_wall(pocket, centre_m, guest.extent_m() + target.cutoff_m);

    let mut report = BindingReport {
        guest: guest.id.to_string(),
        minimum_energy_j: f64::INFINITY,
        centre_energy_j: f64::INFINITY,
        displacement_m: 0.0,
        azimuth_rad: 0.0,
        wall_atom_count: wall.len(),
        guest_atom_count: guest.atom_count(),
        blocked_positions: 0,
        temperature_k: target.temperature_k,
    };

    for step in 0..target.displacements.max(1) {
        let offset_m = if target.displacements <= 1 {
            0.0
        } else {
            -0.5 * target.span_m + target.span_m * step as f64 / (target.displacements - 1) as f64
        };
        for turn in 0..target.azimuths.max(1) {
            let azimuth_rad = std::f64::consts::TAU * turn as f64 / target.azimuths.max(1) as f64;
            let energy_j = interaction_energy_j(
                &wall,
                guest,
                centroid_m,
                centre_m,
                offset_m,
                azimuth_rad,
                target.cutoff_m,
                target.overlap_ratio,
            );
            if !energy_j.is_finite() {
                report.blocked_positions += 1;
                continue;
            }
            if energy_j < report.minimum_energy_j {
                report.minimum_energy_j = energy_j;
                report.displacement_m = offset_m;
                report.azimuth_rad = azimuth_rad;
            }
            if step == target.displacements / 2 && turn == 0 {
                report.centre_energy_j = energy_j;
            }
        }
    }
    report
}

/// Compare two binding energies and report the stronger one.
///
/// The ratio is the energy of the stronger binder over the energy of the
/// weaker one, so a ratio above one means that the pocket prefers the first
/// guest. Both energies must come from the same target.
pub fn selectivity_ratio(strong_j: f64, weak_j: f64) -> f64 {
    if weak_j == 0.0 {
        return 0.0;
    }
    strong_j / weak_j
}

#[allow(clippy::too_many_arguments)]
/// Returns the pocket atoms within `reach_m` of `centre_m`.
///
/// Each entry holds the position in metres, the charge in coulombs, and the
/// element. This is the wall that faces a guest at `centre_m`.
pub(crate) fn pocket_wall(
    pocket: &Part,
    centre_m: [f64; 3],
    reach_m: f64,
) -> Vec<([f64; 3], f64, Element)> {
    let mut wall: Vec<([f64; 3], f64, Element)> = Vec::new();
    for index in 0..pocket.atom_count() {
        let (Some(position_m), Some(element)) = (
            pocket.topology.position_m(index),
            pocket.topology.element(index),
        ) else {
            continue;
        };
        let dx = position_m[0] - centre_m[0];
        let dy = position_m[1] - centre_m[1];
        let dz = position_m[2] - centre_m[2];
        if dx * dx + dy * dy + dz * dz > reach_m * reach_m {
            continue;
        }
        wall.push((
            position_m,
            pocket.topology.charge_c(index).unwrap_or(0.0),
            element,
        ));
    }
    wall
}

// Eight plain arguments stay cheaper than a wrapper struct for the pair kernel.
#[allow(clippy::too_many_arguments)]
pub(crate) fn interaction_energy_j(
    wall: &[([f64; 3], f64, Element)],
    guest: &Guest,
    centroid_m: [f64; 3],
    centre_m: [f64; 3],
    offset_m: f64,
    azimuth_rad: f64,
    cutoff_m: f64,
    overlap_ratio: f64,
) -> f64 {
    let cos = azimuth_rad.cos();
    let sin = azimuth_rad.sin();
    let anchor_m = [centre_m[0], centre_m[1], centre_m[2] + offset_m];
    let cutoff_sq = cutoff_m * cutoff_m;
    let mut energy_j = 0.0;

    for atom in guest.atoms {
        let rx = atom.position_m[0] - centroid_m[0];
        let ry = atom.position_m[1] - centroid_m[1];
        let rz = atom.position_m[2] - centroid_m[2];
        let position_m = [
            anchor_m[0] + cos * rx - sin * ry,
            anchor_m[1] + sin * rx + cos * ry,
            anchor_m[2] + rz,
        ];
        for &(wall_m, wall_charge_c, wall_element) in wall {
            let dx = wall_m[0] - position_m[0];
            let dy = wall_m[1] - position_m[1];
            let dz = wall_m[2] - position_m[2];
            let r_sq = dx * dx + dy * dy + dz * dz;
            if r_sq >= cutoff_sq || r_sq <= 0.0 {
                continue;
            }
            let Some((epsilon_j, sigma_m)) = pair_params(atom.element, wall_element) else {
                continue;
            };
            let floor_m = overlap_ratio * sigma_m;
            if r_sq < floor_m * floor_m {
                return f64::INFINITY;
            }
            energy_j += pair_energy_j(epsilon_j, sigma_m, r_sq);
            let r_m = r_sq.sqrt();
            energy_j += COULOMB_CONSTANT_N_M2_PER_C2 * atom.charge_c * wall_charge_c / r_m;
        }
    }
    energy_j
}

fn pair_energy_j(epsilon_j: f64, sigma_m: f64, r_sq_m2: f64) -> f64 {
    let ratio_sq = sigma_m * sigma_m / r_sq_m2;
    let inverse_6 = ratio_sq * ratio_sq * ratio_sq;
    let inverse_12 = inverse_6 * inverse_6;
    4.0 * epsilon_j * (inverse_12 - inverse_6)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nanocad_parts::{
        BindingPocketGenerator, FunctionalGroup, ParameterSet, PartGenerator, FUNCTIONAL_GROUPS,
    };

    /// The well radius of the test pocket, in metres.
    ///
    /// One contact distance past ethanol, so a bound guest can sit in the
    /// middle of the well without touching the wall.
    const TEST_RADIUS_M: f64 = 7.0e-10;

    fn pocket(group: FunctionalGroup) -> Part {
        let mut parameters = ParameterSet::new();
        parameters.set("pocket_radius_m", TEST_RADIUS_M);
        parameters.set("pocket_depth_m", 1.2e-9);
        parameters.set("wall_m", 0.8e-9);
        parameters.set("thickness_m", 2.0e-9);
        parameters.set("wall_group", group.index() as f64);
        BindingPocketGenerator
            .generate(&parameters)
            .expect("the test pocket builds")
    }

    fn measured(group: FunctionalGroup, id: &str) -> BindingReport {
        let pocket = pocket(group);
        let centre_m = BindingPocketGenerator.cavity_center_m();
        let guest = nanocad_parts::guest(id).expect("the test guest exists");
        binding_energy_j(&pocket, centre_m, guest, &BindingTarget::default())
    }

    #[test]
    fn a_guest_in_the_well_has_a_bound_state() {
        let report = measured(FunctionalGroup::Hydroxyl, "methanol");
        assert!(
            report.minimum_energy_j < 0.0,
            "the guest must sit in a well, not on a barrier: {}",
            report.minimum_energy_j
        );
        assert!(report.over_kt() > 5.0, "a bound guest is worth many kT");
    }

    #[test]
    fn a_guest_that_cannot_fit_is_rejected() {
        let mut parameters = ParameterSet::new();
        parameters.set("pocket_radius_m", 2.3e-10);
        parameters.set("pocket_depth_m", 1.2e-9);
        parameters.set("wall_m", 0.8e-9);
        parameters.set("thickness_m", 2.0e-9);
        parameters.set("wall_group", 0.0);
        let pocket = BindingPocketGenerator
            .generate(&parameters)
            .expect("the small pocket builds");
        let centre_m = BindingPocketGenerator.cavity_center_m();
        let guest = nanocad_parts::guest("ethanol").expect("ethanol exists");
        let report = binding_energy_j(&pocket, centre_m, guest, &BindingTarget::default());
        let positions = BindingTarget::default().displacements * BindingTarget::default().azimuths;
        assert_eq!(
            report.blocked_positions, positions,
            "a well that is too small leaves no free position"
        );
        assert!(report.minimum_energy_j.is_infinite());
    }

    #[test]
    fn a_wider_well_binds_the_same_guest_more_weakly() {
        let narrow = measured(FunctionalGroup::Hydroxyl, "methanol").over_kt();
        let mut parameters = ParameterSet::new();
        parameters.set("pocket_radius_m", 1.2e-9);
        parameters.set("pocket_depth_m", 1.2e-9);
        parameters.set("wall_m", 0.8e-9);
        parameters.set("thickness_m", 2.0e-9);
        parameters.set("wall_group", 0.0);
        let wide = BindingPocketGenerator
            .generate(&parameters)
            .expect("the wide pocket builds");
        let centre_m = BindingPocketGenerator.cavity_center_m();
        let guest = nanocad_parts::guest("methanol").expect("methanol exists");
        let report = binding_energy_j(&wide, centre_m, guest, &BindingTarget::default());
        assert!(
            report.over_kt() < narrow,
            "a wider well must hold the guest less tightly: {} against {narrow}",
            report.over_kt()
        );
    }

    #[test]
    fn the_wall_charge_changes_the_energy() {
        let charged = measured(FunctionalGroup::Hydroxyl, "methanol").minimum_energy_j;
        let mut neutral = pocket(FunctionalGroup::Hydroxyl);
        for index in 0..neutral.topology.atom_count() {
            let _ = neutral.topology.set_charge_c(index, 0.0);
        }
        let centre_m = BindingPocketGenerator.cavity_center_m();
        let guest = nanocad_parts::guest("methanol").expect("methanol exists");
        let bare = binding_energy_j(&neutral, centre_m, guest, &BindingTarget::default());
        assert!(
            charged < bare.minimum_energy_j,
            "the stated charges must add attraction: {charged} against {}",
            bare.minimum_energy_j
        );
    }

    #[test]
    fn a_polar_wall_prefers_the_hydrogen_bond_donor() {
        let ethanol = measured(FunctionalGroup::Hydroxyl, "ethanol").minimum_energy_j;
        let ether = measured(FunctionalGroup::Hydroxyl, "dimethyl_ether").minimum_energy_j;
        assert!(
            ethanol < ether,
            "a hydroxyl wall must hold ethanol more strongly: {ethanol} against {ether}"
        );
        assert!(selectivity_ratio(ethanol, ether) > 1.0);
    }

    #[test]
    fn a_fluorinated_wall_prefers_the_ether() {
        let ethanol = measured(FunctionalGroup::Fluoro, "ethanol").minimum_energy_j;
        let ether = measured(FunctionalGroup::Fluoro, "dimethyl_ether").minimum_energy_j;
        assert!(
            ether < ethanol,
            "a fluoro wall must hold the ether more strongly: {ether} against {ethanol}"
        );
        assert!(selectivity_ratio(ethanol, ether) < 1.0);
    }

    #[test]
    fn the_wall_group_changes_the_selectivity() {
        let ratios: Vec<f64> = FUNCTIONAL_GROUPS
            .iter()
            .map(|group| {
                selectivity_ratio(
                    measured(*group, "ethanol").minimum_energy_j,
                    measured(*group, "dimethyl_ether").minimum_energy_j,
                )
            })
            .collect();
        let smallest = ratios.iter().fold(f64::INFINITY, |a, b| a.min(*b));
        let largest = ratios.iter().fold(f64::NEG_INFINITY, |a, b| a.max(*b));
        assert!(
            largest / smallest > 2.0,
            "the wall chemistry must change the choice: {smallest} against {largest}"
        );
    }

    #[test]
    fn the_selectivity_ratio_is_the_quotient() {
        assert!((selectivity_ratio(-8.0e-21, -4.0e-21) - 2.0).abs() < 1.0e-12);
        assert_eq!(selectivity_ratio(-8.0e-21, 0.0), 0.0);
    }

    #[test]
    fn the_report_states_its_search() {
        let target = BindingTarget::default();
        let report = measured(FunctionalGroup::Hydroxyl, "ethanol");
        assert_eq!(report.guest, "ethanol");
        assert!(report.wall_atom_count > 0);
        assert_eq!(report.guest_atom_count, 9);
        assert_eq!(report.temperature_k, target.temperature_k);
        assert!(report.displacement_m.abs() <= 0.5 * target.span_m);
        assert!((0.0..std::f64::consts::TAU).contains(&report.azimuth_rad));
        assert!(report.blocked_positions <= target.displacements * target.azimuths);
    }

    #[test]
    fn the_metric_value_names_the_guest() {
        let report = measured(FunctionalGroup::Thiol, "ethanol");
        let value = report.to_metric_value();
        assert!(value.name.contains("binding"));
        assert!(value.note.contains("ethanol"));
        assert_eq!(value.fidelity, Fidelity::QuasiStatic);
    }
}
