//! The ejection work.
//!
//! The sorting rotor in Freitas, *Nanomedicine* Volume I, Section 3.4.2, holds
//! its bound molecule until the pocket reaches the inner chamber. There the
//! molecule is forced out: the reference says the bound molecules "are
//! forcibly ejected by rods thrust outward by the cam surface".
//!
//! This module measures the work the rod must supply. The guest sits at the
//! bound position that [`crate::binding`] found. The rod then pushes the guest
//! along the pocket axis, out of the well. The energy rises while the guest
//! passes the wall, and it falls again when the guest leaves. The work is the
//! height of that barrier above the bound energy.
//!
//! The reference states a sortation cost of about 10 to 40 zJ for each
//! molecule, at 310 K. That figure is the reversible work at a stated
//! concentration ratio, so it is a lower bound and not the raw binding energy.
//! [`EjectionReport::work_j`] is the barrier of this rigid model, and
//! [`EjectionReport::over_kt`] states it in thermal units.
//!
//! The pocket and the guest are both rigid, the wall carries stated model
//! charges, and there is no solvent. Every value is a model value. See
//! ADR-0063, ADR-0064 and ADR-0066. All lengths are SI metres.

use nanocad_model::Part;
use nanocad_parts::Guest;

use crate::binding::{binding_energy_j, interaction_energy_j, pocket_wall, BindingTarget};
use crate::{Fidelity, MetricValue};

const BOLTZMANN_J_PER_K: f64 = 1.380_649e-23;

/// The search along the ejection path.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EjectionTarget {
    /// How far the rod pushes the guest, in metres.
    pub travel_m: f64,
    /// How many samples the path takes.
    pub steps: usize,
}

impl Default for EjectionTarget {
    fn default() -> Self {
        Self {
            travel_m: 1.5e-9,
            steps: 48,
        }
    }
}

/// The measured ejection of one guest from one pocket.
#[derive(Clone, Debug, PartialEq)]
pub struct EjectionReport {
    /// The guest id.
    pub guest: String,
    /// The energy at the bound position, in joules.
    pub bound_energy_j: f64,
    /// The offset of the bound position along the pocket axis, in metres.
    pub bound_displacement_m: f64,
    /// The turn of the bound position about the pocket axis, in radians.
    pub bound_azimuth_rad: f64,
    /// The highest energy on the outward path, in joules.
    pub peak_energy_j: f64,
    /// Where the highest energy sits, in metres from the cavity centre.
    pub peak_travel_m: f64,
    /// The energy at the end of the path, in joules.
    pub released_energy_j: f64,
    /// The work the rod must supply, in joules. It is the barrier height.
    pub work_j: f64,
    /// How far the rod pushes the guest, in metres.
    pub travel_m: f64,
    /// How many samples the path took.
    pub steps: usize,
    /// How many samples the guest could not reach, because the wall blocks it.
    pub blocked_positions: usize,
    /// The temperature of the stated thermal energy, in kelvin.
    pub temperature_k: f64,
}

impl EjectionReport {
    /// Returns the thermal energy at the stated temperature, in joules.
    pub fn thermal_energy_j(&self) -> f64 {
        BOLTZMANN_J_PER_K * self.temperature_k
    }

    /// Returns the ejection work in units of the thermal energy.
    pub fn over_kt(&self) -> f64 {
        let thermal_j = self.thermal_energy_j();
        if thermal_j <= 0.0 {
            return 0.0;
        }
        self.work_j / thermal_j
    }

    /// Returns the metric form of this report.
    pub fn to_metric_value(&self) -> MetricValue {
        MetricValue {
            name: format!("ejection_{}", self.guest),
            value: self.work_j,
            unit: "J".to_string(),
            fidelity: Fidelity::QuasiStatic,
            note: format!(
                "{} pushed out of a rigid pocket, {:.2} kT at {:.0} K",
                self.guest,
                self.over_kt(),
                self.temperature_k
            ),
        }
    }
}

/// Measures the work to push `guest` out of `pocket`.
///
/// The bound position comes from [`binding_energy_j`] with the same `binding`
/// target. The path then continues along the pocket axis from that position
/// for `target.travel_m`, at the bound turn. A sample that the wall blocks is
/// counted and left out of the barrier.
pub fn ejection_work_j(
    pocket: &Part,
    centre_m: [f64; 3],
    guest: &Guest,
    binding: &BindingTarget,
    target: &EjectionTarget,
) -> EjectionReport {
    let bound = binding_energy_j(pocket, centre_m, guest, binding);
    let centroid_m = guest.centroid_m();
    let reach_m = guest.extent_m() + binding.cutoff_m + target.travel_m.max(0.0);
    let wall = pocket_wall(pocket, centre_m, reach_m);

    let steps = target.steps.max(2);
    let mut peak_energy_j = f64::NEG_INFINITY;
    let mut peak_travel_m = 0.0;
    let mut released_energy_j = f64::INFINITY;
    let mut blocked_positions = 0;

    for step in 0..steps {
        let along_m = target.travel_m * step as f64 / (steps - 1) as f64;
        let offset_m = bound.displacement_m + along_m;
        let energy_j = interaction_energy_j(
            &wall,
            guest,
            centroid_m,
            centre_m,
            offset_m,
            bound.azimuth_rad,
            binding.cutoff_m,
            binding.overlap_ratio,
        );
        if !energy_j.is_finite() {
            blocked_positions += 1;
            continue;
        }
        if energy_j > peak_energy_j {
            peak_energy_j = energy_j;
            peak_travel_m = along_m;
        }
        released_energy_j = energy_j;
    }

    let work_j = if peak_energy_j.is_finite() && bound.minimum_energy_j.is_finite() {
        peak_energy_j - bound.minimum_energy_j
    } else {
        f64::INFINITY
    };

    EjectionReport {
        guest: guest.id.to_string(),
        bound_energy_j: bound.minimum_energy_j,
        bound_displacement_m: bound.displacement_m,
        bound_azimuth_rad: bound.azimuth_rad,
        peak_energy_j,
        peak_travel_m,
        released_energy_j,
        work_j,
        travel_m: target.travel_m,
        steps,
        blocked_positions,
        temperature_k: binding.temperature_k,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nanocad_parts::{
        guest, BindingPocketGenerator, FunctionalGroup, ParameterSet, PartGenerator,
    };

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
            .expect("the pocket builds")
    }

    fn centre_m() -> [f64; 3] {
        BindingPocketGenerator.cavity_center_m()
    }

    fn measured(group: FunctionalGroup, id: &str) -> EjectionReport {
        let part = pocket(group);
        let guest = guest(id).expect("a guest");
        ejection_work_j(
            &part,
            centre_m(),
            guest,
            &BindingTarget::default(),
            &EjectionTarget::default(),
        )
    }

    #[test]
    fn the_rod_must_supply_work_to_free_the_guest() {
        let report = measured(FunctionalGroup::Hydroxyl, "ethanol");
        assert!(report.work_j.is_finite());
        assert!(report.work_j > 0.0, "the barrier is positive");
        assert!(report.over_kt() > 1.0, "the barrier is above the noise");
    }

    #[test]
    fn the_barrier_rises_above_the_bound_state() {
        let report = measured(FunctionalGroup::Hydroxyl, "ethanol");
        assert!(report.peak_energy_j >= report.bound_energy_j);
        assert!(report.peak_travel_m > 0.0);
    }

    #[test]
    fn the_guest_is_free_at_the_end_of_the_path() {
        let report = measured(FunctionalGroup::Hydroxyl, "ethanol");
        assert!(report.released_energy_j > report.bound_energy_j);
    }

    #[test]
    fn a_deeper_well_costs_more_work() {
        let shallow = measured(FunctionalGroup::Methyl, "methanol");
        let deep = measured(FunctionalGroup::Thiol, "methanol");
        assert!(deep.work_j > shallow.work_j);
    }

    #[test]
    fn the_report_states_its_search() {
        let report = measured(FunctionalGroup::Amino, "dimethyl_ether");
        assert_eq!(report.guest, "dimethyl_ether");
        assert_eq!(report.steps, EjectionTarget::default().steps);
        assert_eq!(report.travel_m, EjectionTarget::default().travel_m);
        assert_eq!(report.temperature_k, 300.0);
    }

    #[test]
    fn the_metric_value_names_the_guest() {
        let report = measured(FunctionalGroup::Hydroxyl, "ethanol");
        let metric = report.to_metric_value();
        assert_eq!(metric.name, "ejection_ethanol");
        assert_eq!(metric.unit, "J");
        assert!(metric.note.contains("ethanol"));
    }

    #[test]
    fn a_zero_travel_measures_only_the_bound_state() {
        let part = pocket(FunctionalGroup::Hydroxyl);
        let guest = guest("ethanol").expect("a guest");
        let report = ejection_work_j(
            &part,
            centre_m(),
            guest,
            &BindingTarget::default(),
            &EjectionTarget {
                travel_m: 0.0,
                steps: 2,
            },
        );
        assert!(report.work_j.abs() < 1.0e-20, "no travel means no work");
    }
}
