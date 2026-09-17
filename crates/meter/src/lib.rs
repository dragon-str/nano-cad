//! Fidelity-tiered performance metrics for generated parts.
//!
//! A metric measures one property of a moving assembly. It reports the value
//! together with the fidelity at which it was computed. Two values of
//! different fidelity are not comparable, so every comparison must read the
//! fidelity field.
#![forbid(unsafe_code)]

pub mod clearance;

pub use clearance::{BodyMotion, Clearance, ClearanceReport, ClearanceTarget, MovingAtoms};

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
