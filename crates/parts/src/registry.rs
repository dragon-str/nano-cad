//! Agent-facing registry of the gear generators.
//!
//! An agent names a gear and supplies a parameter map. The registry dispatches
//! to the matching generator, so the agent never edits coordinates.
//!
//! Note: the Python exposure of this registry is pending M7-01. This module is
//! the Rust surface only.

use nanocad_model::Part;

use crate::error::PartError;
use crate::gear_profile::GearProfileGenerator;
use crate::generator::PartGenerator;
use crate::parameter::ParameterSet;
use crate::planetary::PlanetaryGenerator;
use crate::spur_gear::SpurGearGenerator;

static SPUR_GEAR: SpurGearGenerator = SpurGearGenerator;
static GEAR_PROFILE: GearProfileGenerator = GearProfileGenerator;
static PLANETARY: PlanetaryGenerator = PlanetaryGenerator;

static GEAR_GENERATORS: &[&(dyn PartGenerator + Send + Sync)] =
    &[&SPUR_GEAR, &GEAR_PROFILE, &PLANETARY];

/// Returns the registered gear generators, in a stable order.
pub fn gear_generators() -> &'static [&'static (dyn PartGenerator + Send + Sync)] {
    GEAR_GENERATORS
}

/// Finds a gear generator by its name, or `None`.
pub fn gear_generator(name: &str) -> Option<&'static (dyn PartGenerator + Send + Sync)> {
    GEAR_GENERATORS
        .iter()
        .copied()
        .find(|generator| generator.id() == name)
}

/// Generates a gear from a name and a parameter map.
///
/// The name is the generator id: `"spur_gear"`, `"gear_profile"`, or
/// `"planetary"`. An unknown name is [`PartError::UnknownGenerator`]. The
/// parameter map is passed to the generator, which fills the defaults, rejects
/// unknown names, and checks every range.
pub fn generate_gear(name: &str, specs: &ParameterSet) -> Result<Part, PartError> {
    let generator =
        gear_generator(name).ok_or_else(|| PartError::UnknownGenerator(name.to_owned()))?;
    generator.generate(specs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_registered_generator_has_a_unique_name() {
        let mut names: Vec<&str> = gear_generators().iter().map(|g| g.id()).collect();
        names.sort_unstable();
        let unique = names.len();
        names.dedup();
        assert_eq!(names.len(), unique);
        assert_eq!(unique, 3);
    }

    #[test]
    fn a_spur_gear_dispatches_by_name() {
        let part = generate_gear("spur_gear", &ParameterSet::new()).expect("generate");
        assert!(part.name.starts_with("spur-gear-"));
        assert_eq!(part.material, "diamondoid");
    }

    #[test]
    fn a_gear_profile_dispatches_by_name() {
        let part = generate_gear("gear_profile", &ParameterSet::new()).expect("generate");
        assert!(part.name.starts_with("gear-profile-"));
        assert_eq!(part.material, "geometry");
    }

    #[test]
    fn a_planetary_set_dispatches_by_name() {
        let part = generate_gear("planetary", &ParameterSet::new().with("planet_count", 3.0))
            .expect("generate");
        assert!(part.name.starts_with("planetary-"));
        assert_eq!(part.material, "diamondoid");
    }

    #[test]
    fn an_unknown_name_is_a_typed_error() {
        assert_eq!(
            generate_gear("worm_gear", &ParameterSet::new()),
            Err(PartError::UnknownGenerator("worm_gear".to_owned()))
        );
        assert!(gear_generator("worm_gear").is_none());
    }

    #[test]
    fn a_bad_parameter_reaches_the_dispatched_generator() {
        let bad = ParameterSet::new().with("tooth_count", 5.0);
        assert!(matches!(
            generate_gear("spur_gear", &bad),
            Err(PartError::ToothCountBelowUndercut { .. })
        ));
    }
}
