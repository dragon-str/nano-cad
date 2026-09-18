//! Builds a sorting rotor, its housing, its cam plate and its drive shaft, and
//! prints the facts.
//!
//! The rotor follows Freitas, *Nanomedicine* Volume I, Section 3.4.2: a disk
//! with a row of binding pockets on its rim, and rods that the cam surface
//! thrusts outward to eject a bound guest. The example reports the simulated
//! geometry, the mass and the kinematics at a stated turn rate. It also prints
//! the figures from the reference next to the simulated ones, so the reader can
//! see the agreement and the difference.
//!
//! The example does not model molecular selectivity and it does not model a
//! solvent. Those need a chemistry force field.
//!
//! With one argument, the example also writes the rotor scene to that path, in
//! the schema of the gearbox scene, so a viewer shows the rotor.

use std::f64::consts::TAU;

use nanocad_model::Element;
use nanocad_parts::{
    place, CamPlateGenerator, DriveShaftGenerator, EjectionRodGenerator, FollowerPinGenerator,
    ParameterSet, PartGenerator, RotorHousingGenerator, SortingRotorGenerator,
};

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

/// The number of ejection rods, one for each pocket.
const ROD_COUNT: f64 = 12.0;

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

fn resolved_m(generator: &impl PartGenerator, name: &str) -> f64 {
    generator
        .resolve(&ParameterSet::new())
        .ok()
        .and_then(|values| values.get(name))
        .unwrap_or(0.0)
}

fn main() {
    let scene_path = std::env::args().nth(1);
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
    let cam = match CamPlateGenerator.generate_with_defaults() {
        Ok(part) => part,
        Err(error) => {
            eprintln!("the cam plate failed: {error}");
            std::process::exit(1);
        }
    };
    let shaft = match DriveShaftGenerator.generate_with_defaults() {
        Ok(part) => part,
        Err(error) => {
            eprintln!("the drive shaft failed: {error}");
            std::process::exit(1);
        }
    };
    let pin = match FollowerPinGenerator.generate_with_defaults() {
        Ok(part) => part,
        Err(error) => {
            eprintln!("the follower pin failed: {error}");
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

    let pocket_orbit_m = resolved_m(&SortingRotorGenerator, "pocket_orbit_m");
    let pocket_radius_m = resolved_m(&SortingRotorGenerator, "pocket_radius_m");
    let pocket_inner_m = pocket_orbit_m - pocket_radius_m;
    let tip_length_m = resolved_m(&EjectionRodGenerator, "tip_length_m");
    let cam_retract_m = CamPlateGenerator.groove_base_radius_m();
    let stroke_m = CamPlateGenerator.groove_rise_m();
    let ramp_rad = CamPlateGenerator.groove_ramp_half_angle_rad();
    let rod_length_m = pocket_inner_m - cam_retract_m;
    let rod = match EjectionRodGenerator
        .generate(&ParameterSet::new().with("shaft_length_m", rod_length_m - tip_length_m))
    {
        Ok(rod) => rod,
        Err(error) => {
            eprintln!("the rod failed: {error}");
            std::process::exit(1);
        }
    };

    let rotor_atoms = placed.atom_count();
    let housing_atoms = housing.atom_count();
    let cam_atoms = cam.atom_count();
    let rod_atoms = rod.atom_count();
    let pin_atoms = pin.atom_count();
    let shaft_atoms = shaft.atom_count();
    let rotor_mass_kg = total_mass_kg(&placed);
    let housing_mass_kg = total_mass_kg(&housing);
    let cam_mass_kg = total_mass_kg(&cam);
    let rod_mass_kg = total_mass_kg(&rod);
    let pin_mass_kg = total_mass_kg(&pin);
    let shaft_mass_kg = total_mass_kg(&shaft);
    let total_atoms = rotor_atoms
        + housing_atoms
        + cam_atoms
        + shaft_atoms
        + (ROD_COUNT as usize) * (rod_atoms + pin_atoms);
    let total_mass_kg = rotor_mass_kg
        + housing_mass_kg
        + cam_mass_kg
        + shaft_mass_kg
        + ROD_COUNT * (rod_mass_kg + pin_mass_kg);

    if let Some(path) = scene_path {
        let scene = match nanocad_jigs::build_rotor_scene() {
            Ok(scene) => scene,
            Err(error) => {
                eprintln!("the rotor scene failed: {error}");
                std::process::exit(1);
            }
        };
        if let Err(error) =
            nanocad_jigs::write_scene_json(&scene, std::path::Path::new(&path), false)
        {
            eprintln!("the rotor scene was not written: {error}");
            std::process::exit(1);
        }
    }

    println!(
        "{{\"rotor_atoms\":{rotor_atoms},\"rotor_bonds\":{},\"rotor_radius_m\":{rotor_radius_m:e},\
         \"pocket_count\":{pocket_count},\"rotor_mass_kg\":{rotor_mass_kg:e},\
         \"housing_atoms\":{housing_atoms},\"housing_mass_kg\":{housing_mass_kg:e},\
         \"cam_atoms\":{cam_atoms},\"cam_mass_kg\":{cam_mass_kg:e},\
         \"cam_plate_atoms\":{cam_atoms},\"cam_plate_mass_kg\":{cam_mass_kg:e},\
         \"shaft_atoms\":{shaft_atoms},\"shaft_mass_kg\":{shaft_mass_kg:e},\
         \"pin_atoms\":{pin_atoms},\"pin_mass_kg\":{pin_mass_kg:e},\
         \"rod_count\":{ROD_COUNT},\"rod_atoms\":{rod_atoms},\"rod_mass_kg\":{rod_mass_kg:e},\
         \"rod_length_m\":{rod_length_m:e},\"rod_stroke_m\":{stroke_m:e},\"rod_ramp_rad\":{ramp_rad:e},\
         \"total_atoms\":{total_atoms},\"total_mass_kg\":{total_mass_kg:e},\
         \"rate_rad_per_s\":{rate_rad_per_s:e},\"revolutions_per_s\":{revolutions_per_s:e},\
         \"rim_speed_m_per_s\":{rim_speed_m_per_s:e},\"pocket_cycle_s\":{pocket_cycle_s:e},\
         \"pocket_rate_per_s\":{pocket_rate_per_s:e},\
         \"reference_revolutions_per_s\":{REFERENCE_REVOLUTIONS_PER_S:e},\
         \"reference_rim_speed_m_per_s\":{REFERENCE_RIM_SPEED_M_PER_S:e},\
         \"reference_atoms\":{REFERENCE_ATOMS:e},\"reference_mass_kg\":{REFERENCE_MASS_KG:e}}}",
        placed.bond_count(),
    );
}
