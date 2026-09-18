//! Parametric part generators: gear, nanotube, lattice.
//!
//! A generator implements [`PartGenerator`]. It declares its parameter specs
//! and builds a [`nanocad_model::Part`] from a validated [`ParameterSet`]. The
//! crate holds SI positions and unit-suffixed names, following ADR-0003.
#![forbid(unsafe_code)]

pub mod axle;
pub mod bearing;
pub mod block;
pub mod clutch;
pub mod diamond_solid;
pub mod dislocation;
pub mod error;
pub mod gear_profile;
pub mod generator;
pub mod geometry;
mod group;
mod group_data;
mod guest;
mod guest_data;
pub mod housing;
pub mod lattice;
pub mod lattice_fill;
pub mod nanotube;
pub mod parameter;
pub mod placement;
pub mod planetary;
mod pocket;
pub mod port;
pub mod registry;
pub mod respirocyte;
pub mod rod;
pub mod rotor;
pub mod shape;
pub mod spur_gear;
pub mod validation;

pub use axle::{HexAxleGenerator, PlainShaftGenerator};
pub use bearing::{BushingGenerator, RadialBearingGenerator};
pub use block::{BeamGenerator, BracketGenerator, PlateGenerator};
pub use clutch::{ClutchPlateGenerator, RatchetGenerator};
pub use dislocation::{
    core_atom_mask, displace, displace_point_m, five_seven_wedge_rad, WedgeDisclination,
};
pub use error::PartError;
pub use gear_profile::{GearProfile, GearProfileGenerator};
pub use generator::{PartGenerator, PartSchema};
pub use group::{FunctionalGroup, FUNCTIONAL_GROUPS};
pub use guest::{
    build_guest, build_guest_part, guest, guests, Guest, GuestAtom, GuestBond, AVOGADRO_PER_MOL,
};
pub use housing::RotorHousingGenerator;
pub use lattice::{DiamondGenerator, GraphiteGenerator};
pub use lattice_fill::fill_solid;
pub use nanotube::NanotubeGenerator;
pub use parameter::{ParameterSet, ParameterSpec};
pub use placement::{assembly_document, place, PlacedPart};
pub use planetary::{
    planetary_constraint_holds, PlanetaryDesign, PlanetaryGenerator, PlanetarySet,
};
pub use pocket::{
    pocket_radius_for, wall_contact_distance_m, wall_well_depth_j, BindingPocketGenerator,
};
pub use port::{connect, mate_offset_m, Dof, Port, PortFrame};
pub use registry::{
    gear_generator, gear_generators, generate, generate_gear, generator, library,
    library_categories, LibraryEntry,
};
pub use respirocyte::{
    respirocyte_pump, respirocyte_rotor, respirocyte_tank, AnnularGapFlow, RespirocytePump,
    RespirocytePumpGenerator, RespirocyteRotor, RespirocyteRotorGenerator, RespirocyteTank,
    RespirocyteTankGenerator,
};
pub use rod::EjectionRodGenerator;
pub use rotor::SortingRotorGenerator;
pub use shape::{
    Bounds, Box3, Cylinder, Difference, HexPrism, Intersection, Placed, Profile, Solid, Union,
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
