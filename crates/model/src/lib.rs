//! Document model: atoms, bonds, topology, parts.
//!
//! The model stores atoms and bonds in struct-of-arrays form and addresses
//! them by integer index. Positions are in SI metres and charges are in SI
//! coulombs. `Document` is the versioned, serializable design graph.
#![forbid(unsafe_code)]

mod atom;
mod bond;
mod chemistry;
mod document;
mod element;
mod encoding;
mod error;
mod nonbonded;
mod part;
mod selection;
mod topology;

pub use atom::Atom;
pub use bond::{Bond, BondType};
pub use chemistry::{
    atomic_mass_kg, bond_length_m, chemistry, covalent_radius_m, valence, ElementChemistry,
    ATOMIC_MASS_UNIT_KG,
};
pub use document::{Document, SCHEMA_VERSION};
pub use element::Element;
pub use error::ModelError;
pub use nonbonded::{
    nonbonded, pair_params, vdw_distance_m, well_depth_j, NonbondedParams, KILOCALORIE_PER_MOL_J,
    LENNARD_JONES_LENGTH_FACTOR,
};
pub use part::Part;
pub use selection::{parse_selection, Selection, SelectionResult};
pub use topology::Topology;

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
