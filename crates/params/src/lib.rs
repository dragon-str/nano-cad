//! Multiscale parameter store, provenance, verification.
//!
//! This crate implements the central data contract of `PARAMETERS.md`. A
//! [`PartRecord`] carries the measured or computed properties of one part. Each
//! property carries a unit and an uncertainty. The record carries a method, a
//! validation status, and a [`Provenance`] block.
//!
//! Values are stored in SI base units. The unit label is only a presentation
//! detail. This follows ADR-0003: SI internally, and conversion at the boundary.
//!
//! Unknown uncertainty is [`None`], never zero. This follows Rule 3 of
//! `PARAMETERS.md`.
#![forbid(unsafe_code)]

mod check;
mod error;
mod failure;
mod friction;
mod library;
mod provenance;
mod quantity;
mod record;
mod stats;
mod stiffness;
mod strain;
mod thermal;

#[cfg(test)]
mod test_support;

pub use check::{check, is_consistent, Problem};
pub use error::ParamError;
pub use failure::{extract_failure_stress, BondLimit, FailureConfig, FailureResult};
pub use friction::{extract_friction, summarize_friction, FrictionConfig, FrictionResult};
pub use library::{LibraryEntry, ParameterLibrary, LIBRARY_SCHEMA, LIBRARY_VERSION};
pub use provenance::{Method, Provenance, Validation};
pub use quantity::Quantity;
pub use record::{PartRecord, SCHEMA, VERSION};
pub use stiffness::{extract_stiffness, StiffnessConfig, StiffnessResult};
pub use thermal::{
    extract_specific_heat, extract_thermal, extract_thermal_conductivity, ConductivityResult,
    SpecificHeatResult, ThermalConfig, ThermalResult,
};

/// Returns the crate version string.
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_is_not_empty() {
        assert!(!version().is_empty());
    }

    #[test]
    fn public_constants_are_stable() {
        assert_eq!(SCHEMA, "nanocad.param.part");
        assert_eq!(VERSION, 1);
    }
}
