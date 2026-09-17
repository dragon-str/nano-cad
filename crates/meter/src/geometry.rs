//! Cheap geometric metrics.
//!
//! These metrics read the generator output only. They need no forces and no
//! time integration, so they cost almost nothing and they carry the
//! [`Fidelity::Geometric`] label.

use std::f64::consts::PI;

use nanocad_parts::planetary::{PlanetaryDesign, PlanetarySet};
use nanocad_parts::{GearProfile, PartError};

use crate::{Fidelity, MetricValue};

/// Counts the atoms of a generated part.
pub fn atom_count(set: &PlanetarySet) -> MetricValue {
    let atoms = set.part.atom_count();
    MetricValue {
        name: "atom_count".to_string(),
        value: atoms as f64,
        unit: "atoms".to_string(),
        fidelity: Fidelity::Geometric,
        note: format!("{atoms} atoms in the part"),
    }
}

/// The contact ratio of the sun and planet mesh.
///
/// The contact ratio is the mean number of tooth pairs in contact. A value
/// below one means the drive can skip a tooth. A value near two means a smooth
/// drive with more sliding.
pub fn contact_ratio(design: &PlanetaryDesign) -> Result<MetricValue, PartError> {
    let module_m = design.module_m();
    let pressure_angle_rad = design.pressure_angle_rad();
    let sun = GearProfile::new(module_m, design.sun_teeth(), pressure_angle_rad)?;
    let planet = GearProfile::new(module_m, design.planet_teeth(), pressure_angle_rad)?;
    let centre_m = sun.pitch_radius_m() + planet.pitch_radius_m();
    let path_of_contact_m = (sun.outer_radius_m().powi(2) - sun.base_radius_m().powi(2)).sqrt()
        + (planet.outer_radius_m().powi(2) - planet.base_radius_m().powi(2)).sqrt()
        - centre_m * pressure_angle_rad.sin();
    let base_pitch_m = PI * module_m * pressure_angle_rad.cos();
    let ratio = path_of_contact_m / base_pitch_m;
    Ok(MetricValue {
        name: "contact_ratio".to_string(),
        value: ratio,
        unit: "".to_string(),
        fidelity: Fidelity::Geometric,
        note: format!(
            "sun and planet mesh; {:.2} tooth pairs in contact on average",
            ratio
        ),
    })
}
