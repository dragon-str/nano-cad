//! Agent-facing registry of the part generators.
//!
//! An agent names a part and supplies a parameter map. The registry dispatches
//! to the matching generator, so the agent never edits coordinates.
//!
//! Note: the Python exposure of this registry is pending M7-01. This module is
//! the Rust surface only.

use std::sync::OnceLock;

use nanocad_model::Part;

use crate::axle::{HexAxleGenerator, PlainShaftGenerator};
use crate::bearing::{BushingGenerator, RadialBearingGenerator};
use crate::block::{BeamGenerator, BracketGenerator, PlateGenerator};
use crate::clutch::{ClutchPlateGenerator, RatchetGenerator};
use crate::error::PartError;
use crate::gear_profile::GearProfileGenerator;
use crate::generator::PartGenerator;
use crate::housing::RotorHousingGenerator;
use crate::lattice::{DiamondGenerator, GraphiteGenerator};
use crate::nanotube::NanotubeGenerator;
use crate::parameter::ParameterSet;
use crate::planetary::PlanetaryGenerator;
use crate::respirocyte::{
    RespirocytePumpGenerator, RespirocyteRotorGenerator, RespirocyteTankGenerator,
};
use crate::rotor::SortingRotorGenerator;
use crate::spur_gear::SpurGearGenerator;

static SPUR_GEAR: SpurGearGenerator = SpurGearGenerator;
static GEAR_PROFILE: GearProfileGenerator = GearProfileGenerator;
static PLANETARY: PlanetaryGenerator = PlanetaryGenerator;
static DIAMOND: DiamondGenerator = DiamondGenerator;
static GRAPHITE: GraphiteGenerator = GraphiteGenerator;
static NANOTUBE: NanotubeGenerator = NanotubeGenerator;
static HEX_AXLE: HexAxleGenerator = HexAxleGenerator;
static PLAIN_SHAFT: PlainShaftGenerator = PlainShaftGenerator;
static PLATE: PlateGenerator = PlateGenerator;
static BEAM: BeamGenerator = BeamGenerator;
static BRACKET: BracketGenerator = BracketGenerator;
static RADIAL_BEARING: RadialBearingGenerator = RadialBearingGenerator;
static BUSHING: BushingGenerator = BushingGenerator;
static CLUTCH_PLATE: ClutchPlateGenerator = ClutchPlateGenerator;
static RATCHET: RatchetGenerator = RatchetGenerator;
static RESPIROCYTE_ROTOR: RespirocyteRotorGenerator = RespirocyteRotorGenerator;
static RESPIROCYTE_PUMP: RespirocytePumpGenerator = RespirocytePumpGenerator;
static RESPIROCYTE_TANK: RespirocyteTankGenerator = RespirocyteTankGenerator;
static SORTING_ROTOR: SortingRotorGenerator = SortingRotorGenerator;
static ROTOR_HOUSING: RotorHousingGenerator = RotorHousingGenerator;

/// A generator plus its library category.
struct RegisteredGenerator {
    generator: &'static (dyn PartGenerator + Send + Sync),
    category: &'static str,
}

static GENERATORS: &[RegisteredGenerator] = &[
    RegisteredGenerator {
        generator: &SPUR_GEAR,
        category: "gears",
    },
    RegisteredGenerator {
        generator: &GEAR_PROFILE,
        category: "gears",
    },
    RegisteredGenerator {
        generator: &PLANETARY,
        category: "gears",
    },
    RegisteredGenerator {
        generator: &DIAMOND,
        category: "lattice",
    },
    RegisteredGenerator {
        generator: &GRAPHITE,
        category: "lattice",
    },
    RegisteredGenerator {
        generator: &NANOTUBE,
        category: "lattice",
    },
    RegisteredGenerator {
        generator: &HEX_AXLE,
        category: "structure",
    },
    RegisteredGenerator {
        generator: &PLAIN_SHAFT,
        category: "structure",
    },
    RegisteredGenerator {
        generator: &PLATE,
        category: "structure",
    },
    RegisteredGenerator {
        generator: &BEAM,
        category: "structure",
    },
    RegisteredGenerator {
        generator: &BRACKET,
        category: "structure",
    },
    RegisteredGenerator {
        generator: &RADIAL_BEARING,
        category: "device",
    },
    RegisteredGenerator {
        generator: &BUSHING,
        category: "device",
    },
    RegisteredGenerator {
        generator: &CLUTCH_PLATE,
        category: "device",
    },
    RegisteredGenerator {
        generator: &RATCHET,
        category: "device",
    },
    RegisteredGenerator {
        generator: &RESPIROCYTE_ROTOR,
        category: "device",
    },
    RegisteredGenerator {
        generator: &RESPIROCYTE_PUMP,
        category: "device",
    },
    RegisteredGenerator {
        generator: &RESPIROCYTE_TANK,
        category: "device",
    },
    RegisteredGenerator {
        generator: &SORTING_ROTOR,
        category: "device",
    },
    RegisteredGenerator {
        generator: &ROTOR_HOUSING,
        category: "device",
    },
];

static GEAR_GENERATORS: &[&(dyn PartGenerator + Send + Sync)] =
    &[&SPUR_GEAR, &GEAR_PROFILE, &PLANETARY];

static LIBRARY: OnceLock<Vec<LibraryEntry>> = OnceLock::new();

/// One entry in the general part library.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LibraryEntry {
    /// The stable generator id.
    pub id: &'static str,
    /// The human-readable generator name.
    pub name: &'static str,
    /// The library category, for example `"gears"`.
    pub category: &'static str,
}

/// Returns every library entry, in the registry order.
///
/// The table is built once from [`GENERATORS`], so the entry ids and names come
/// from the generators themselves and cannot drift.
pub fn library() -> &'static [LibraryEntry] {
    LIBRARY
        .get_or_init(|| {
            GENERATORS
                .iter()
                .map(|registered| LibraryEntry {
                    id: registered.generator.id(),
                    name: registered.generator.name(),
                    category: registered.category,
                })
                .collect()
        })
        .as_slice()
}

/// Finds a generator by its id, or `None`.
pub fn generator(name: &str) -> Option<&'static (dyn PartGenerator + Send + Sync)> {
    GENERATORS
        .iter()
        .map(|registered| registered.generator)
        .find(|generator| generator.id() == name)
}

/// Generates a part from a generator id and a parameter map.
///
/// An unknown id is [`PartError::UnknownGenerator`]. The parameter map is
/// passed to the generator, which fills the defaults, rejects unknown names,
/// and checks every range.
pub fn generate(name: &str, specs: &ParameterSet) -> Result<Part, PartError> {
    let generator = generator(name).ok_or_else(|| PartError::UnknownGenerator(name.to_owned()))?;
    generator.generate(specs)
}

/// Returns the library categories, sorted and unique.
pub fn library_categories() -> Vec<&'static str> {
    let mut categories: Vec<&'static str> = library().iter().map(|entry| entry.category).collect();
    categories.sort_unstable();
    categories.dedup();
    categories
}

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

    #[test]
    fn the_library_lists_every_generator_with_a_unique_id() {
        let mut ids: Vec<&str> = library().iter().map(|entry| entry.id).collect();
        assert_eq!(library().len(), GENERATORS.len());
        let unique = ids.len();
        ids.sort_unstable();
        ids.dedup();
        assert_eq!(ids.len(), unique, "library ids are not unique");
        assert_eq!(unique, GENERATORS.len());
    }

    #[test]
    fn every_library_entry_names_a_generator() {
        for entry in library() {
            let found = generator(entry.id).expect("entry id resolves");
            assert_eq!(found.id(), entry.id);
            assert_eq!(found.name(), entry.name);
            assert!(
                library_categories().contains(&entry.category),
                "category {} is missing from library_categories",
                entry.category
            );
        }
    }

    #[test]
    fn an_unknown_id_is_a_typed_error() {
        assert_eq!(
            generate("worm_gear", &ParameterSet::new()),
            Err(PartError::UnknownGenerator("worm_gear".to_owned()))
        );
        assert!(generator("worm_gear").is_none());
    }

    #[test]
    fn the_categories_are_sorted_and_unique() {
        let categories = library_categories();
        let mut sorted = categories.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(categories, sorted);
        assert!(categories.contains(&"gears"));
        assert!(categories.contains(&"structure"));
        assert!(categories.contains(&"lattice"));
        assert!(categories.contains(&"device"));
    }
}
