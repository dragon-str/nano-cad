//! Fidelity-tiered performance metrics for generated parts.
//!
//! A metric measures one property of a moving assembly. It reports the value
//! together with the fidelity at which it was computed. Two values of
//! different fidelity are not comparable, so every comparison must read the
//! fidelity field.
#![forbid(unsafe_code)]

pub mod binding;
pub mod bonded;
pub mod clearance;
pub mod drive;
pub mod geometry;
pub mod harmonic;
pub mod loaded;
pub mod relaxed;
pub mod slip;

pub use binding::{binding_energy_j, selectivity_ratio, BindingReport, BindingTarget};
pub use clearance::{BodyMotion, Clearance, ClearanceReport, ClearanceTarget, MovingAtoms};
pub use drive::{DriveReport, DriveTarget, SteeredDrive};
pub use geometry::{atom_count, contact_ratio};
pub use harmonic::{HarmonicMesh, HarmonicReport, HarmonicTarget};
pub use loaded::{LoadedContact, LoadedContactReport, LoadedContactTarget};
pub use relaxed::{RelaxedSlipBarrier, RelaxedSlipBarrierReport, RelaxedSlipBarrierTarget};
pub use slip::{SlipBarrier, SlipBarrierReport, SlipBarrierTarget};

use nanocad_parts::planetary::PlanetarySet;
use nanocad_parts::PartError;

/// The fidelity at which a metric was computed.
///
/// The order is the cost order. A cheaper tier screens, and a dearer tier
/// verifies. Never compare two values of a different tier.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Fidelity {
    /// Pure geometry. No forces and no time integration.
    Geometric,
    /// Relaxed configurations along a path. No time integration.
    QuasiStatic,
    /// A Hessian and its normal modes.
    Harmonic,
    /// Time integration with a potential.
    Dynamics,
}

/// One measured value, with its unit, its fidelity, and a short note.
#[derive(Clone, Debug, PartialEq)]
pub struct MetricValue {
    /// The metric name, for example `clearance`.
    pub name: String,
    /// The measured value in the units of the metric.
    pub value: f64,
    /// The unit symbol, for example `m`.
    pub unit: String,
    /// The fidelity at which the value was computed.
    pub fidelity: Fidelity,
    /// A short note on the method or the limit of the value.
    pub note: String,
}

/// A named collection of metric values for one design.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Score {
    /// The metric values, in insertion order.
    pub values: Vec<MetricValue>,
}

impl Score {
    /// Appends one metric value.
    pub fn push(&mut self, value: MetricValue) {
        self.values.push(value);
    }

    /// Returns the first value with the given name.
    pub fn get(&self, name: &str) -> Option<&MetricValue> {
        self.values.iter().find(|value| value.name == name)
    }
}

/// Scores a planetary set with every metric that is available now.
///
/// The sun rate sets the motion. The ring stays fixed.
pub fn score_planetary(set: &PlanetarySet, sun_rad_per_s: f64) -> Result<Score, PartError> {
    let mut score = Score::default();
    score.push(geometry::atom_count(set));
    score.push(geometry::contact_ratio(&set.design)?);
    let moving = MovingAtoms::planetary(set, &set.design, sun_rad_per_s);
    let clearance = Clearance::default();
    score.push(
        clearance
            .measure(&moving)
            .to_metric_value(&clearance.target),
    );
    let slip = SlipBarrier::default();
    let rigid = slip.measure(&moving, set.design.sun_teeth());
    score.push(rigid.to_metric_value(&slip.target));
    let relaxed = RelaxedSlipBarrier::default();
    score.push(
        relaxed
            .measure(&moving, set.design.sun_teeth())
            .to_metric_value(&relaxed.target),
    );
    score.push(HarmonicMesh::default().measure(set).to_metric_value());
    let contact = LoadedContact::default();
    score.push(
        contact
            .measure(&moving, set.design.sun_teeth())
            .to_metric_value(),
    );
    Ok(score)
}
