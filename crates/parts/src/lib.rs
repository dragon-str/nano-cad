//! Parametric part generators: gear, nanotube, lattice.
//!
//! A generator implements [`PartGenerator`]. It declares its parameter specs
//! and builds a [`nanocad_model::Part`] from a validated [`ParameterSet`]. The
//! crate holds SI positions and unit-suffixed names, following ADR-0003.
#![forbid(unsafe_code)]

pub mod error;
pub mod gear_profile;
pub mod generator;
pub mod geometry;
pub mod lattice;
pub mod nanotube;
pub mod parameter;
pub mod planetary;
pub mod registry;
pub mod respirocyte;
pub mod spur_gear;
pub mod validation;

pub use error::PartError;
pub use gear_profile::{GearProfile, GearProfileGenerator};
pub use generator::{PartGenerator, PartSchema};
pub use lattice::{DiamondGenerator, GraphiteGenerator};
pub use nanotube::NanotubeGenerator;
pub use parameter::{ParameterSet, ParameterSpec};
pub use planetary::{
    planetary_constraint_holds, PlanetaryDesign, PlanetaryGenerator, PlanetarySet,
};
pub use registry::{gear_generator, gear_generators, generate_gear};
pub use respirocyte::{
    respirocyte_pump, respirocyte_rotor, respirocyte_tank, AnnularGapFlow, RespirocytePump,
    RespirocytePumpGenerator, RespirocyteRotor, RespirocyteRotorGenerator, RespirocyteTank,
    RespirocyteTankGenerator,
};
pub use spur_gear::SpurGearGenerator;
pub use validation::{
    element_valence, validate_part, ValidationReport, Violation, DEFAULT_STRAIN_TOLERANCE_RELATIVE,
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
}
