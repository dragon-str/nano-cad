//! Builds a sorting rotor and its housing, mates them, and prints the facts.
//!
//! The rotor follows Freitas, *Nanomedicine* Volume I, Section 3.4.2: a disk
//! with a row of binding pockets on its rim. The example reports the simulated
//! geometry, the mass and the kinematics at a stated turn rate. It also prints
//! the figures from the reference next to the simulated ones, so the reader can
//! see the agreement and the difference.
//!
//! The example does not model molecular selectivity and it does not model a
//! solvent. Those need a chemistry force field.

use std::f64::consts::TAU;

use nanocad_model::Element;
use nanocad_parts::{place, PartGenerator, RotorHousingGenerator, SortingRotorGenerator};

/// The atomic mass unit in kilograms.
const ATOMIC_MASS_UNIT_KG: f64 = 1.660_539_066_60e-27;

/// The standard atomic weight of carbon.
const CARBON_ATOMIC_WEIGHT: f64 = 12.011;

/// The standard atomic weight of hydrogen.
const HYDROGEN_ATOMIC_WEIGHT: f64 = 1.008;

/// The reference turn rate from Freitas, in revolutions per second.
const REFERENCE_REVOLUTIONS_PER_S: f64 = 86_000.0;

/// The reference rim speed from Freitas, in metres per second.
const REFERENCE_RIM_SPEED_M_PER_S: f64 = 2.7e-3;

/// The reference atom count from Freitas.
const REFERENCE_ATOMS: f64 = 1.0e5;

/// The reference mass from Freitas, in kilograms.
const REFERENCE_MASS_KG: f64 = 2.0e-21;

fn mass_kg(element: Element) -> f64 {
    if element == Element::CARBON {
        CARBON_ATOMIC_WEIGHT * ATOMIC_MASS_UNIT_KG
    } else {
        HYDROGEN_ATOMIC_WEIGHT * ATOMIC_MASS_UNIT_KG
    }
}

fn total_mass_kg(part: &nanocad_model::Part) -> f64 {
    (0..part.atom_count())
        .filter_map(|index| part.topology.element(index))
        .map(mass_kg)
        .sum()
}

fn main() {
    let rotor = match SortingRotorGenerator.generate_with_defaults() {
        Ok(part) => part,
        Err(error) => {
            eprintln!("the rotor failed: {error}");
            std::process::exit(1);
        }
    };
    let housing = match RotorHousingGenerator.generate_with_defaults() {
        Ok(part) => part,
        Err(error) => {
            eprintln!("the housing failed: {error}");
            std::process::exit(1);
        }
    };

    let frame = place(
        rotor,
        &SortingRotorGenerator.axis_port(),
        &RotorHousingGenerator.chamber_port(),
        0.0,
    );
    let placed = frame.transformed_part();
    let rotor_radius_m = 7.0e-9;
    let pocket_count = 12.0_f64;

    let rate_rad_per_s = TAU * REFERENCE_REVOLUTIONS_PER_S;
    let revolutions_per_s = rate_rad_per_s / TAU;
    let rim_speed_m_per_s = rate_rad_per_s * rotor_radius_m;
    let pocket_rate_per_s = revolutions_per_s * pocket_count;
    let pocket_cycle_s = 1.0 / revolutions_per_s;

    let rotor_atoms = placed.atom_count();
    let housing_atoms = housing.atom_count();
    let rotor_mass_kg = total_mass_kg(&placed);
    let housing_mass_kg = total_mass_kg(&housing);

    println!(
        "{{\"rotor_atoms\":{rotor_atoms},\"rotor_bonds\":{},\"rotor_radius_m\":{rotor_radius_m:e},\
         \"pocket_count\":{pocket_count},\"rotor_mass_kg\":{rotor_mass_kg:e},\
         \"housing_atoms\":{housing_atoms},\"housing_mass_kg\":{housing_mass_kg:e},\
         \"total_atoms\":{},\"total_mass_kg\":{:e},\
         \"rate_rad_per_s\":{rate_rad_per_s:e},\"revolutions_per_s\":{revolutions_per_s:e},\
         \"rim_speed_m_per_s\":{rim_speed_m_per_s:e},\"pocket_cycle_s\":{pocket_cycle_s:e},\
         \"pocket_rate_per_s\":{pocket_rate_per_s:e},\
         \"reference_revolutions_per_s\":{REFERENCE_REVOLUTIONS_PER_S:e},\
         \"reference_rim_speed_m_per_s\":{REFERENCE_RIM_SPEED_M_PER_S:e},\
         \"reference_atoms\":{REFERENCE_ATOMS:e},\"reference_mass_kg\":{REFERENCE_MASS_KG:e}}}",
        placed.bond_count(),
        rotor_atoms + housing_atoms,
        rotor_mass_kg + housing_mass_kg,
    );
}
