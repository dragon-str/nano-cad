//! Assemble a planetary gearbox into the L2 device layer.
//!
//! The assembly builds the sun, planet, ring, and carrier geometry with the
//! `nanocad-parts` generator. It converts each body to a
//! [`nanocad_params::PartRecord`] and stores the records in a
//! [`nanocad_params::ParameterLibrary`]. It then builds the
//! [`RigidBodySystem`] from those records. This is the M6-05 path: L1 geometry
//! to an L2 PartRecord to an L2 device.
//!
//! # Device topology
//!
//! The bodies are `ground`, `sun`, one `planet` per planet, `ring`, and
//! `carrier`. The ring is fixed to the ground. The revolute joints are:
//!
//! - sun to ground, about the z axis at the origin,
//! - carrier to ground, about the z axis at the origin,
//! - ring to ground, about the z axis at the origin,
//! - each planet to the carrier, about the z axis at its own centre.
//!
//! The gear mesh uses [`GearConstraint`] terms with signed coefficients. With
//! `theta` a revolute joint angle, the external sun-planet mesh and the
//! internal planet-ring mesh are
//!
//! ```text
//! r_sun * theta_sun - r_sun * theta_carrier + r_planet * theta_planet = 0
//! r_ring * theta_ring - (r_carrier + r_planet) * theta_carrier
//!     - r_planet * theta_planet = 0
//! ```
//!
//! With the ring fixed, these give `omega_carrier = r_sun * omega_sun /
//! (2 r_carrier)` and the sun-to-carrier ratio `(N_sun + N_ring) / N_sun`.

use nanocad_model::Topology;
use nanocad_params::{Method, ParameterLibrary, PartRecord, Provenance, Quantity, Validation};
use nanocad_parts::{ParameterSet, PartError, PlanetaryDesign, PlanetaryGenerator, PlanetarySet};
use nanocad_units::Unit;
use thiserror::Error;

use crate::device::{DeviceError, GearConstraint, Quat, RevoluteJoint, RigidBody, RigidBodySystem};
use crate::rotor::CARBON_ATOM_MASS_KG;

/// The library version used for the assembled part records.
pub const ASSEMBLY_RECORD_VERSION: u32 = 1;

/// The diamond C-C bond length, in metres. Sets the carrier body sampling.
const CARBON_BOND_LENGTH_M: f64 = 1.544e-10;

/// The part id of the sun record.
pub const SUN_PART_ID: &str = "planetary.sun";
/// The part id prefix of the planet records, with the planet index appended.
pub const PLANET_PART_ID_PREFIX: &str = "planetary.planet.";
/// The part id of the ring record.
pub const RING_PART_ID: &str = "planetary.ring";
/// The part id of the carrier record.
pub const CARRIER_PART_ID: &str = "planetary.carrier";

/// Errors from gearbox assembly.
#[derive(Debug, Error, PartialEq)]
pub enum AssemblyError {
    #[error(transparent)]
    Part(#[from] PartError),
    #[error(transparent)]
    Device(#[from] DeviceError),
    #[error("parameter library error: {0}")]
    Param(String),
    #[error("missing part record {part_id:?} in the library")]
    MissingRecord { part_id: String },
    #[error("part record {part_id:?} is invalid: {reason}")]
    InvalidRecord { part_id: String, reason: String },
    #[error("part record {part_id:?} has a non-finite {field}")]
    NonFiniteRecord {
        part_id: String,
        field: &'static str,
    },
}

/// The role of one body in the assembled gearbox.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BodyRole {
    /// The fixed world body.
    Ground,
    /// The central sun gear.
    Sun,
    /// One planet gear, with its index.
    Planet(usize),
    /// The fixed internal ring gear.
    Ring,
    /// The output carrier.
    Carrier,
}

/// An assembled planetary gearbox.
#[derive(Clone, Debug)]
pub struct PlanetaryAssembly {
    /// The device.
    pub system: RigidBodySystem,
    /// The parameter library that holds one [`PartRecord`] per body.
    pub library: ParameterLibrary,
    /// The derived geometry.
    pub design: PlanetaryDesign,
    /// The role of each body, in body-index order.
    pub roles: Vec<BodyRole>,
    /// The ground body index.
    pub ground_body: usize,
    /// The sun body index.
    pub sun_body: usize,
    /// The sun joint index.
    pub sun_joint: usize,
    /// The carrier body index.
    pub carrier_body: usize,
    /// The carrier joint index.
    pub carrier_joint: usize,
    /// The ring body index.
    pub ring_body: usize,
    /// The ring joint index.
    pub ring_joint: usize,
    /// The planet body indices, in planet order.
    pub planet_bodies: Vec<usize>,
    /// The planet joint indices, in planet order.
    pub planet_joints: Vec<usize>,
    /// The gear constraint indices, in construction order.
    pub gear_constraints: Vec<usize>,
}

impl PlanetaryAssembly {
    /// Returns the analytic sun-to-carrier gear ratio with the ring fixed.
    pub fn analytic_gear_ratio(&self) -> f64 {
        self.design.gear_ratio()
    }

    /// Returns the sun joint angle in radians.
    pub fn sun_angle_rad(&self) -> Result<f64, DeviceError> {
        self.system.revolute_joints()[self.sun_joint].angle_rad(self.system.bodies())
    }

    /// Returns the carrier joint angle in radians.
    pub fn carrier_angle_rad(&self) -> Result<f64, DeviceError> {
        self.system.revolute_joints()[self.carrier_joint].angle_rad(self.system.bodies())
    }

    /// Returns the ring joint angle in radians.
    pub fn ring_angle_rad(&self) -> Result<f64, DeviceError> {
        self.system.revolute_joints()[self.ring_joint].angle_rad(self.system.bodies())
    }

    /// Returns one planet joint angle in radians.
    pub fn planet_angle_rad(&self, planet: usize) -> Result<f64, DeviceError> {
        let joint =
            self.planet_joints
                .get(planet)
                .copied()
                .ok_or(DeviceError::JointIndexOutOfBounds {
                    joint: planet,
                    joint_count: self.planet_joints.len(),
                })?;
        self.system.revolute_joints()[joint].angle_rad(self.system.bodies())
    }

    /// Sets every body to a rotation that satisfies the gear constraints.
    ///
    /// The sun turns at `sun_angular_velocity_rad_per_s`. The carrier and
    /// planets follow the analytic ratio, so the device starts on the
    /// constraint manifold and the ratio test has no start-up transient.
    pub fn set_consistent_rotation(
        &mut self,
        sun_angular_velocity_rad_per_s: f64,
    ) -> Result<(), AssemblyError> {
        if !sun_angular_velocity_rad_per_s.is_finite() {
            return Err(AssemblyError::NonFiniteRecord {
                part_id: "rotation".to_owned(),
                field: "sun angular velocity",
            });
        }
        let radius_sun_m = self.design.sun_pitch_radius_m();
        let radius_planet_m = self.design.planet_pitch_radius_m();
        let radius_carrier_m = self.design.carrier_radius_m();
        let carrier_angular_velocity_rad_per_s =
            radius_sun_m * sun_angular_velocity_rad_per_s / (2.0 * radius_carrier_m);
        let planet_relative_angular_velocity_rad_per_s = -(radius_carrier_m + radius_planet_m)
            / radius_planet_m
            * carrier_angular_velocity_rad_per_s;
        let planet_angular_velocity_rad_per_s =
            carrier_angular_velocity_rad_per_s + planet_relative_angular_velocity_rad_per_s;

        set_body_omega(
            &mut self.system,
            self.sun_body,
            sun_angular_velocity_rad_per_s,
        )?;
        set_body_omega(
            &mut self.system,
            self.carrier_body,
            carrier_angular_velocity_rad_per_s,
        )?;
        for body in &self.planet_bodies {
            set_body_omega(&mut self.system, *body, planet_angular_velocity_rad_per_s)?;
        }
        Ok(())
    }
}

fn set_body_omega(
    system: &mut RigidBodySystem,
    body: usize,
    angular_velocity_rad_per_s: f64,
) -> Result<(), AssemblyError> {
    let body_count = system.body_count();
    let body = system
        .body_mut(body)
        .ok_or(DeviceError::BodyIndexOutOfBounds { body, body_count })?;
    body.set_angular_velocity_rad_per_s([0.0, 0.0, angular_velocity_rad_per_s])?;
    Ok(())
}

fn set_body_fixed(system: &mut RigidBodySystem, body: usize) -> Result<(), AssemblyError> {
    let body_count = system.body_count();
    let body = system
        .body_mut(body)
        .ok_or(DeviceError::BodyIndexOutOfBounds { body, body_count })?;
    body.set_fixed(true);
    Ok(())
}

/// The mass and inertia of one geometric body, about its centre of mass.
#[derive(Clone, Copy, Debug, PartialEq)]
struct BodyMassProperties {
    mass_kg: f64,
    center_of_mass_m: [f64; 3],
    inertia_diagonal_kg_m2: [f64; 3],
}

fn body_mass_properties(
    topology: &Topology,
    atom_range: std::ops::Range<usize>,
) -> Result<BodyMassProperties, AssemblyError> {
    let mut mass_kg = 0.0;
    let mut center_of_mass_m = [0.0; 3];
    let mut atom_count = 0usize;
    for atom in atom_range.clone() {
        let position_m = topology
            .position_m(atom)
            .ok_or_else(|| AssemblyError::InvalidRecord {
                part_id: "geometry".to_owned(),
                reason: format!("atom {atom} is outside the topology"),
            })?;
        mass_kg += CARBON_ATOM_MASS_KG;
        for axis in 0..3 {
            center_of_mass_m[axis] += CARBON_ATOM_MASS_KG * position_m[axis];
        }
        atom_count += 1;
    }
    if atom_count == 0 {
        return Err(AssemblyError::InvalidRecord {
            part_id: "geometry".to_owned(),
            reason: "the atom range is empty".to_owned(),
        });
    }
    for value in center_of_mass_m.iter_mut() {
        *value /= mass_kg;
    }
    let mut inertia_diagonal_kg_m2 = [0.0; 3];
    for atom in atom_range {
        let position_m = topology
            .position_m(atom)
            .ok_or_else(|| AssemblyError::InvalidRecord {
                part_id: "geometry".to_owned(),
                reason: format!("atom {atom} is outside the topology"),
            })?;
        let dx = position_m[0] - center_of_mass_m[0];
        let dy = position_m[1] - center_of_mass_m[1];
        let dz = position_m[2] - center_of_mass_m[2];
        inertia_diagonal_kg_m2[0] += CARBON_ATOM_MASS_KG * (dy * dy + dz * dz);
        inertia_diagonal_kg_m2[1] += CARBON_ATOM_MASS_KG * (dx * dx + dz * dz);
        inertia_diagonal_kg_m2[2] += CARBON_ATOM_MASS_KG * (dx * dx + dy * dy);
    }
    Ok(BodyMassProperties {
        mass_kg,
        center_of_mass_m,
        inertia_diagonal_kg_m2,
    })
}

fn provenance() -> Provenance {
    Provenance {
        source: "nanocad-jigs assembly".to_owned(),
        method: Method::Estimate,
        code_version: crate::version().to_owned(),
        force_field: None,
        timestamp: "2026-09-14T00:00:00Z".to_owned(),
        uncertainty: None,
        validation: Validation::Unverified,
        notes: "mass and inertia from the generated carbon geometry".to_owned(),
    }
}

fn record_for_body(
    part_id: &str,
    topology: &Topology,
    atom_range: std::ops::Range<usize>,
    material: &str,
) -> Result<PartRecord, AssemblyError> {
    let atoms = u64::try_from(atom_range.len()).unwrap_or(u64::MAX);
    let properties = body_mass_properties(topology, atom_range)?;
    record_from_properties(part_id, atoms, material, properties)
}

fn record_from_properties(
    part_id: &str,
    atoms: u64,
    material: &str,
    properties: BodyMassProperties,
) -> Result<PartRecord, AssemblyError> {
    let mass_kg = Quantity::from_si(properties.mass_kg, Unit::Kilogram);
    let inertia_kg_m2 = [
        Quantity::derived(properties.inertia_diagonal_kg_m2[0], "kg*m^2")
            .map_err(|error| AssemblyError::Param(error.to_string()))?,
        Quantity::derived(properties.inertia_diagonal_kg_m2[1], "kg*m^2")
            .map_err(|error| AssemblyError::Param(error.to_string()))?,
        Quantity::derived(properties.inertia_diagonal_kg_m2[2], "kg*m^2")
            .map_err(|error| AssemblyError::Param(error.to_string()))?,
    ];
    let geometry_ref = format!("nanocad-jigs:{part_id}:{}", mass_kg.value_si().to_bits());
    Ok(PartRecord {
        schema: nanocad_params::SCHEMA.to_owned(),
        version: nanocad_params::VERSION,
        part_id: part_id.to_owned(),
        geometry_ref,
        material: material.to_owned(),
        atoms,
        mass_kg,
        inertia_kg_m2,
        elastic_modulus_pa: Quantity::derived(1.05e12, "Pa")
            .map_err(|error| AssemblyError::Param(error.to_string()))?,
        shear_modulus_pa: Quantity::derived(4.4e11, "Pa")
            .map_err(|error| AssemblyError::Param(error.to_string()))?,
        poisson_ratio: Quantity::derived(0.2, "1")
            .map_err(|error| AssemblyError::Param(error.to_string()))?,
        failure_stress_pa: Quantity::derived(3.0e10, "Pa")
            .map_err(|error| AssemblyError::Param(error.to_string()))?,
        friction_coefficient: Quantity::derived(0.05, "1")
            .map_err(|error| AssemblyError::Param(error.to_string()))?,
        thermal_conductivity_w_m_k: Quantity::derived(1000.0, "W/(m*K)")
            .map_err(|error| AssemblyError::Param(error.to_string()))?,
        specific_heat_j_kg_k: Quantity::derived(500.0, "J/(kg*K)")
            .map_err(|error| AssemblyError::Param(error.to_string()))?,
        method: Method::Estimate,
        validation: Validation::Unverified,
        provenance: provenance(),
    })
}

/// Builds one [`PartRecord`] per body from a generated planetary set.
///
/// The mass is the carbon atom count times the carbon-12 mass. The inertia is
/// about the body centre of mass. The order is sun, one planet per planet,
/// ring, then carrier.
pub fn planetary_records(set: &PlanetarySet) -> Result<Vec<(BodyRole, PartRecord)>, AssemblyError> {
    let topology = &set.part.topology;
    let material = &set.part.material;
    let mut records = Vec::new();
    records.push((
        BodyRole::Sun,
        record_for_body(SUN_PART_ID, topology, set.sun_atoms.clone(), material)?,
    ));
    for (planet, atoms) in set.planet_atoms.iter().enumerate() {
        records.push((
            BodyRole::Planet(planet),
            record_for_body(
                &format!("{PLANET_PART_ID_PREFIX}{planet}"),
                topology,
                atoms.clone(),
                material,
            )?,
        ));
    }
    records.push((
        BodyRole::Ring,
        record_for_body(RING_PART_ID, topology, set.ring_atoms.clone(), material)?,
    ));
    records.push((BodyRole::Carrier, carrier_record(set, material)?));
    Ok(records)
}

/// The carrier carries no atoms, so its mass and inertia are analytic.
///
/// The carrier is a diamondoid ring at the carrier radius, with one carbon
/// per C-C bond length around the circumference. The z inertia is `m R^2`
/// for a thin ring; the x and y inertia are half that.
fn carrier_record(set: &PlanetarySet, material: &str) -> Result<PartRecord, AssemblyError> {
    let radius_m = set.design.carrier_radius_m();
    if !radius_m.is_finite() || radius_m <= 0.0 {
        return Err(AssemblyError::InvalidRecord {
            part_id: CARRIER_PART_ID.to_owned(),
            reason: "the carrier radius must be positive".to_owned(),
        });
    }
    let circumference_m = 2.0 * std::f64::consts::PI * radius_m;
    let count = (circumference_m / CARBON_BOND_LENGTH_M).round().max(3.0) as u64;
    let mass_kg = count as f64 * CARBON_ATOM_MASS_KG;
    let properties = BodyMassProperties {
        mass_kg,
        center_of_mass_m: [0.0; 3],
        inertia_diagonal_kg_m2: [
            0.5 * mass_kg * radius_m * radius_m,
            0.5 * mass_kg * radius_m * radius_m,
            mass_kg * radius_m * radius_m,
        ],
    };
    record_from_properties(CARRIER_PART_ID, count, material, properties)
}

fn record_mass_properties(record: &PartRecord) -> Result<BodyMassProperties, AssemblyError> {
    let mass_kg = record.mass_kg.value_si();
    if !mass_kg.is_finite() || mass_kg <= 0.0 {
        return Err(AssemblyError::NonFiniteRecord {
            part_id: record.part_id.clone(),
            field: "mass",
        });
    }
    let mut inertia_diagonal_kg_m2 = [0.0; 3];
    for (axis, inertia) in record.inertia_kg_m2.iter().enumerate() {
        let value = inertia.value_si();
        if !value.is_finite() || value <= 0.0 {
            return Err(AssemblyError::NonFiniteRecord {
                part_id: record.part_id.clone(),
                field: "inertia",
            });
        }
        inertia_diagonal_kg_m2[axis] = value;
    }
    Ok(BodyMassProperties {
        mass_kg,
        center_of_mass_m: [0.0; 3],
        inertia_diagonal_kg_m2,
    })
}

fn rigid_body_from_record(
    record: &PartRecord,
    position_m: [f64; 3],
) -> Result<RigidBody, AssemblyError> {
    let properties = record_mass_properties(record)?;
    Ok(RigidBody::new(
        properties.mass_kg,
        [
            [properties.inertia_diagonal_kg_m2[0], 0.0, 0.0],
            [0.0, properties.inertia_diagonal_kg_m2[1], 0.0],
            [0.0, 0.0, properties.inertia_diagonal_kg_m2[2]],
        ],
        properties.center_of_mass_m,
        position_m,
        Quat::IDENTITY,
    )?)
}

/// Assembles the device from a parameter library and a design.
///
/// The library must hold the records that [`planetary_records`] produced. The
/// design places the bodies. The function reads only mass and inertia from the
/// records; the geometry comes from the design.
pub fn assemble_from_records(
    library: &ParameterLibrary,
    design: PlanetaryDesign,
) -> Result<PlanetaryAssembly, AssemblyError> {
    let sun = require_record(library, SUN_PART_ID)?;
    let ring = require_record(library, RING_PART_ID)?;
    let carrier = require_record(library, CARRIER_PART_ID)?;

    let mut system = RigidBodySystem::new();
    let mut roles = Vec::new();

    let ground =
        RigidBody::with_diagonal_inertia(1.0, [1.0, 1.0, 1.0], [0.0, 0.0, 0.0], Quat::IDENTITY)?;
    let ground_body = system.add_body(ground);
    roles.push(BodyRole::Ground);
    set_body_fixed(&mut system, ground_body)?;

    let sun_body = system.add_body(rigid_body_from_record(sun, [0.0, 0.0, 0.0])?);
    roles.push(BodyRole::Sun);

    let mut planet_bodies = Vec::with_capacity(design.planet_count());
    for planet in 0..design.planet_count() {
        let part_id = format!("{PLANET_PART_ID_PREFIX}{planet}");
        let record = require_record(library, &part_id)?;
        let center_m = design.planet_center_m(planet);
        let body = system.add_body(rigid_body_from_record(
            record,
            [center_m[0], center_m[1], 0.0],
        )?);
        roles.push(BodyRole::Planet(planet));
        planet_bodies.push(body);
    }

    let ring_body = system.add_body(rigid_body_from_record(ring, [0.0, 0.0, 0.0])?);
    roles.push(BodyRole::Ring);
    set_body_fixed(&mut system, ring_body)?;

    let carrier_body = system.add_body(rigid_body_from_record(carrier, [0.0, 0.0, 0.0])?);
    roles.push(BodyRole::Carrier);

    let axis = [0.0, 0.0, 1.0];
    let sun_joint = system.add_revolute_joint(RevoluteJoint::from_world(
        system.bodies(),
        ground_body,
        sun_body,
        [0.0, 0.0, 0.0],
        axis,
    )?);
    let carrier_joint = system.add_revolute_joint(RevoluteJoint::from_world(
        system.bodies(),
        ground_body,
        carrier_body,
        [0.0, 0.0, 0.0],
        axis,
    )?);
    let ring_joint = system.add_revolute_joint(RevoluteJoint::from_world(
        system.bodies(),
        ground_body,
        ring_body,
        [0.0, 0.0, 0.0],
        axis,
    )?);
    let mut planet_joints = Vec::with_capacity(design.planet_count());
    for (planet, body) in planet_bodies.iter().enumerate() {
        let center_m = design.planet_center_m(planet);
        let joint = system.add_revolute_joint(RevoluteJoint::from_world(
            system.bodies(),
            carrier_body,
            *body,
            [center_m[0], center_m[1], 0.0],
            axis,
        )?);
        planet_joints.push(joint);
    }

    let radius_sun_m = design.sun_pitch_radius_m();
    let radius_planet_m = design.planet_pitch_radius_m();
    let radius_ring_m = design.ring_pitch_radius_m();
    let radius_carrier_m = design.carrier_radius_m();

    let mut gear_constraints = Vec::new();
    for joint in &planet_joints {
        let external = GearConstraint::new([
            (sun_joint, radius_sun_m),
            (carrier_joint, -radius_sun_m),
            (*joint, radius_planet_m),
        ])?;
        gear_constraints.push(system.add_gear_constraint(external)?);
        let internal = GearConstraint::new([
            (ring_joint, radius_ring_m),
            (carrier_joint, -(radius_carrier_m + radius_planet_m)),
            (*joint, -radius_planet_m),
        ])?;
        gear_constraints.push(system.add_gear_constraint(internal)?);
    }

    Ok(PlanetaryAssembly {
        system,
        library: library.clone(),
        design,
        roles,
        ground_body,
        sun_body,
        sun_joint,
        carrier_body,
        carrier_joint,
        ring_body,
        ring_joint,
        planet_bodies,
        planet_joints,
        gear_constraints,
    })
}

fn require_record<'a>(
    library: &'a ParameterLibrary,
    part_id: &str,
) -> Result<&'a PartRecord, AssemblyError> {
    library
        .get(part_id)
        .ok_or_else(|| AssemblyError::MissingRecord {
            part_id: part_id.to_owned(),
        })
}

/// Generates a planetary set, writes its records, and assembles the device.
///
/// This is the full M6-05 path. The returned assembly holds the device and the
/// library of part records that produced it.
pub fn assemble_planetary(parameters: &ParameterSet) -> Result<PlanetaryAssembly, AssemblyError> {
    let set = PlanetaryGenerator.build(parameters)?;
    let records = planetary_records(&set)?;
    let mut library = ParameterLibrary::new();
    for (_, record) in &records {
        library
            .insert(record.clone(), ASSEMBLY_RECORD_VERSION)
            .map_err(|error| AssemblyError::Param(error.to_string()))?;
    }
    let mut assembly = assemble_from_records(&library, set.design)?;
    assembly.library = library;
    Ok(assembly)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_parameters() -> ParameterSet {
        ParameterSet::new().with("planet_count", 3.0)
    }

    #[test]
    fn every_body_has_a_record_and_the_assembly_builds() {
        let assembly = assemble_planetary(&default_parameters()).expect("assemble");
        assert_eq!(assembly.library.len(), 6);
        assert_eq!(assembly.system.body_count(), 7);
        assert_eq!(assembly.system.revolute_joints().len(), 6);
        assert_eq!(assembly.system.gear_constraints().len(), 6);
        assert_eq!(assembly.system.gear_couplings().len(), 0);
        assert_eq!(assembly.roles.len(), 7);
        assert!(assembly
            .system
            .body(assembly.ring_body)
            .expect("ring")
            .is_fixed());
        assert!(assembly
            .system
            .body(assembly.ground_body)
            .expect("ground")
            .is_fixed());
        assert!(!assembly
            .system
            .body(assembly.carrier_body)
            .expect("carrier")
            .is_fixed());
    }

    #[test]
    fn the_assembly_starts_satisfied() {
        let assembly = assemble_planetary(&default_parameters()).expect("assemble");
        assert!(
            assembly
                .system
                .max_revolute_anchor_error_m()
                .expect("valid")
                < 1.0e-15
        );
        assert!(
            assembly
                .system
                .max_revolute_axis_error_rad()
                .expect("valid")
                < 1.0e-12
        );
        assert!(
            assembly
                .system
                .max_gear_constraint_error_m()
                .expect("valid")
                < 1.0e-24
        );
    }

    #[test]
    fn the_records_come_from_the_generated_geometry() {
        let set = PlanetaryGenerator
            .build(&default_parameters())
            .expect("generate");
        let records = planetary_records(&set).expect("records");
        assert_eq!(records.len(), 6);
        let sun_mass_kg = records[0].1.mass_kg.value_si();
        let expected_kg = set.sun_atoms.len() as f64 * CARBON_ATOM_MASS_KG;
        assert!((sun_mass_kg - expected_kg).abs() / expected_kg < 1.0e-12);
        for (_, record) in &records {
            assert!(record.check_header().is_ok());
            assert!(record.mass_kg.value_si() > 0.0);
        }
    }

    #[test]
    fn a_missing_record_is_a_typed_error() {
        let library = ParameterLibrary::new();
        let design = PlanetaryGenerator
            .build(&default_parameters())
            .expect("generate")
            .design;
        assert!(matches!(
            assemble_from_records(&library, design),
            Err(AssemblyError::MissingRecord { .. })
        ));
    }

    #[test]
    fn the_gear_ratio_matches_the_analytic_ratio_within_tolerance() {
        let mut assembly = assemble_planetary(&default_parameters()).expect("assemble");
        let analytic = assembly.analytic_gear_ratio();
        assert!((analytic - 3.5).abs() < 1.0e-12);

        let sun_angular_velocity_rad_per_s = 0.5;
        assembly
            .set_consistent_rotation(sun_angular_velocity_rad_per_s)
            .expect("consistent rotation");

        let sun_start_rad = assembly.sun_angle_rad().expect("valid");
        let carrier_start_rad = assembly.carrier_angle_rad().expect("valid");

        let dt_s = 1.0e-3;
        let steps = 2_000;
        let mut max_constraint_error_m = 0.0_f64;
        let mut max_anchor_error_m = 0.0_f64;
        let mut max_axis_error_rad = 0.0_f64;
        for _ in 0..steps {
            assembly.system.step(dt_s, 64).expect("stable step");
            max_constraint_error_m = max_constraint_error_m.max(
                assembly
                    .system
                    .max_gear_constraint_error_m()
                    .expect("valid"),
            );
            max_anchor_error_m = max_anchor_error_m.max(
                assembly
                    .system
                    .max_revolute_anchor_error_m()
                    .expect("valid"),
            );
            max_axis_error_rad = max_axis_error_rad.max(
                assembly
                    .system
                    .max_revolute_axis_error_rad()
                    .expect("valid"),
            );
        }

        let sun_delta_rad = assembly.sun_angle_rad().expect("valid") - sun_start_rad;
        let carrier_delta_rad = assembly.carrier_angle_rad().expect("valid") - carrier_start_rad;
        let measured_ratio = sun_delta_rad / carrier_delta_rad;
        let relative_error = (measured_ratio - analytic).abs() / analytic;
        let tolerance = 1.0e-3;
        println!(
            "gearbox: analytic_ratio {analytic:.9} measured_ratio {measured_ratio:.9} relative_error {relative_error:.3e} tolerance {tolerance:.1e} max_constraint_error {max_constraint_error_m:.3e} m max_anchor_error {max_anchor_error_m:.3e} m max_axis_error {max_axis_error_rad:.3e} rad"
        );
        assert!(
            relative_error < tolerance,
            "measured ratio {measured_ratio} vs analytic {analytic}"
        );
        assert!(
            max_constraint_error_m < 1.0e-9,
            "maximum gear constraint error {max_constraint_error_m} m"
        );
        assert_eq!(assembly.ring_angle_rad().expect("valid"), 0.0);
    }
}
