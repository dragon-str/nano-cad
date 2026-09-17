//! Three-scale scene export for the gearbox visualization.
//!
//! The export has three layers. A viewer can animate one layer and cross-fade
//! to the next.
//!
//! - [`AtomisticLayer`] is the L1 view: every atom with its element and its
//!   position.
//! - [`DeviceLayer`] is the L2 view: every rigid body with its pose and its
//!   velocities, plus the joints, the gear couplings, and the gear
//!   constraints.
//! - [`CoarseLayer`] is the handoff view: a gear pitch circle for each gear
//!   and a bounding cylinder for each body.
//!
//! # Coordinate system and units
//!
//! The frame is right-handed and Cartesian. The z axis is the gear axis, so
//! the gears lie in the x-y plane. Length is SI metres. Angle is SI radians.
//! An orientation is a unit quaternion `[w, x, y, z]`. It maps a body-frame
//! vector to the world frame.
//!
//! The atomistic positions are the world-frame positions at the zero
//! configuration, exactly as the `nanocad-parts` generator emitted them. The
//! device poses are the current poses of the [`RigidBodySystem`]. The two
//! agree at the zero configuration. To animate an atomistic body, a viewer
//! applies the device transform of the body named by [`SceneAtom::body`].
//!
//! # Determinism
//!
//! The export is deterministic. The body order, the joint order, and the atom
//! order are the construction order. The code does not iterate a map.
//!
//! [`RigidBodySystem`]: crate::device::RigidBodySystem

use std::ops::Range;
use std::path::Path;

use serde::{Deserialize, Serialize};
use thiserror::Error;

use nanocad_parts::PlanetarySet;

use crate::assembly::{
    BodyRole, PlanetaryAssembly, CARRIER_PART_ID, PLANET_PART_ID_PREFIX, RING_PART_ID, SUN_PART_ID,
};
use crate::device::DeviceError;

/// The schema name of a scene document.
pub const SCENE_SCHEMA: &str = "nanocad.scene";
/// The schema version of a scene document.
pub const SCENE_VERSION: u32 = 1;

/// Errors from scene building and from scene serialization.
#[derive(Debug, Error, PartialEq)]
pub enum SceneError {
    #[error(transparent)]
    Device(#[from] DeviceError),
    #[error("the assembly has {roles} roles but {bodies} bodies")]
    BodyCountMismatch { roles: usize, bodies: usize },
    #[error("the assembly design does not match the generated planetary set")]
    DesignMismatch,
    #[error("body {body} is missing from the device system")]
    MissingBody { body: usize },
    #[error("planet {planet} is missing from the generated set")]
    MissingPlanet { planet: usize },
    #[error("the atom range {start}..{end} for role {role:?} exceeds the atom count {atom_count}")]
    AtomRangeOutOfBounds {
        role: String,
        start: usize,
        end: usize,
        atom_count: usize,
    },
    #[error("atom {index} is missing from the topology")]
    MissingAtom { index: usize },
    #[error("atom {index} belongs to more than one body")]
    DuplicateAtom { index: usize },
    #[error("atom {index} belongs to no body")]
    UnassignedAtom { index: usize },
    #[error("atom {index} has a non-finite position {position_m:?} m")]
    NonFiniteAtom { index: usize, position_m: [f64; 3] },
    #[error("body {body} has a non-finite position {position_m:?} m")]
    NonFiniteBody { body: usize, position_m: [f64; 3] },
    #[error("JSON error: {0}")]
    Json(String),
    #[error("I/O error: {0}")]
    Io(String),
}

/// The derived geometry of the exported mechanism.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct SceneDesign {
    /// The gear module in metres.
    pub module_m: f64,
    /// The sun tooth count.
    pub sun_teeth: usize,
    /// One planet tooth count.
    pub planet_teeth: usize,
    /// The ring tooth count.
    pub ring_teeth: usize,
    /// The number of planets.
    pub planet_count: usize,
    /// The sun-to-carrier gear ratio with the ring fixed.
    pub gear_ratio: f64,
    /// The sun pitch radius in metres.
    pub sun_pitch_radius_m: f64,
    /// One planet pitch radius in metres.
    pub planet_pitch_radius_m: f64,
    /// The ring pitch radius in metres.
    pub ring_pitch_radius_m: f64,
    /// The carrier pin-circle radius in metres.
    pub carrier_radius_m: f64,
    /// The number of atomic layers in the axial gear thickness.
    pub layers: usize,
    /// The separation between adjacent atomic layers in metres.
    pub layer_spacing_m: f64,
    /// The gear thickness in metres, equal to `(layers - 1) * layer_spacing_m`.
    pub thickness_m: f64,
}

/// One atom in the L1 layer.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SceneAtom {
    /// The chemical symbol, or `Z<n>` when the symbol is unknown.
    pub element: String,
    /// The atomic number.
    pub atomic_number: u8,
    /// The world position at the zero configuration, in metres.
    pub position_m: [f64; 3],
    /// The index of the device body that owns this atom.
    pub body: usize,
    /// The unique body name, for example `planet_0`.
    pub name: String,
    /// The part id of the owning body, or an empty string for the ground.
    pub part_id: String,
}

/// The L1 atomistic layer.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct AtomisticLayer {
    /// The number of atoms. It equals `atoms.len()`.
    pub atom_count: usize,
    /// The atoms in topology index order.
    pub atoms: Vec<SceneAtom>,
}

/// One rigid body in the L2 layer.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SceneBody {
    /// The body index, which is the index in the device system.
    pub index: usize,
    /// The unique body name, for example `planet_0`.
    pub name: String,
    /// The body role, for example `planet`.
    pub role: String,
    /// The part id, or `null` for the ground.
    pub part_id: Option<String>,
    /// The world position of the centre of mass, in metres.
    pub position_m: [f64; 3],
    /// The orientation as a unit quaternion `[w, x, y, z]`.
    pub orientation_wxyz: [f64; 4],
    /// The mass in kilograms.
    pub mass_kg: f64,
    /// The principal moments of inertia in kilograms metre squared.
    pub inertia_diagonal_kg_m2: [f64; 3],
    /// True when the body does not move.
    pub fixed: bool,
    /// The world linear velocity in metres per second.
    pub linear_velocity_m_per_s: [f64; 3],
    /// The world angular velocity in radians per second.
    pub angular_velocity_rad_per_s: [f64; 3],
}

/// One joint in the L2 layer.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SceneJoint {
    /// The unique joint name.
    pub name: String,
    /// The joint kind, `revolute` or `prismatic`.
    pub kind: String,
    /// The first body index.
    pub body_a: usize,
    /// The second body index.
    pub body_b: usize,
    /// The world anchor point in metres. It is `null` for a prismatic joint.
    pub anchor_world_m: Option<[f64; 3]>,
    /// The common axis in the world frame.
    pub axis_world: [f64; 3],
    /// The unit axis in the first body frame.
    pub axis_a_body: [f64; 3],
    /// The unit axis in the second body frame. It is `null` for a prismatic
    /// joint, which does not expose the value.
    pub axis_b_body: Option<[f64; 3]>,
}

/// One gear coupling in the L2 layer.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SceneGearCoupling {
    /// The first revolute joint index.
    pub joint_a: usize,
    /// The second revolute joint index.
    pub joint_b: usize,
    /// The first pitch radius in metres.
    pub radius_a_m: f64,
    /// The second pitch radius in metres.
    pub radius_b_m: f64,
}

/// One term of a gear constraint.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct SceneGearTerm {
    /// The revolute joint index.
    pub joint: usize,
    /// The signed coefficient in metres.
    pub coefficient_m: f64,
}

/// One gear constraint in the L2 layer.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SceneGearConstraint {
    /// The signed joint terms.
    pub terms: Vec<SceneGearTerm>,
}

/// The L2 device layer.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DeviceLayer {
    /// The number of bodies. It equals `bodies.len()`.
    pub body_count: usize,
    /// The bodies in device index order.
    pub bodies: Vec<SceneBody>,
    /// The joints, revolute then prismatic.
    pub joints: Vec<SceneJoint>,
    /// The gear couplings.
    pub gear_couplings: Vec<SceneGearCoupling>,
    /// The gear constraints.
    pub gear_constraints: Vec<SceneGearConstraint>,
}

/// A circle in the coarse layer.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct CoarseCircle {
    /// The circle centre in metres.
    pub center_m: [f64; 3],
    /// The circle normal.
    pub axis: [f64; 3],
    /// The circle radius in metres.
    pub radius_m: f64,
}

/// A bounding cylinder in the coarse layer.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct CoarseCylinder {
    /// The cylinder centre in metres.
    pub center_m: [f64; 3],
    /// The cylinder axis.
    pub axis: [f64; 3],
    /// The cylinder radius in metres.
    pub radius_m: f64,
    /// The cylinder half length in metres.
    pub half_length_m: f64,
}

/// One body in the coarse layer.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CoarseBody {
    /// The device body index.
    pub index: usize,
    /// The unique body name.
    pub name: String,
    /// The body role.
    pub role: String,
    /// The gear pitch circle, or `null` when the body is not a gear.
    pub pitch_circle: Option<CoarseCircle>,
    /// The bounding cylinder of the body atoms, or `null` for the ground.
    pub bounding_cylinder: Option<CoarseCylinder>,
}

/// The coarse layer. This is the handoff representation.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CoarseLayer {
    /// The number of bodies. It equals `bodies.len()`.
    pub body_count: usize,
    /// The bodies in device index order.
    pub bodies: Vec<CoarseBody>,
}

/// A three-scale scene.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Scene {
    /// The schema name, [`SCENE_SCHEMA`].
    pub schema: String,
    /// The schema version, [`SCENE_VERSION`].
    pub version: u32,
    /// The length unit, always `m`.
    pub length_unit: String,
    /// The angle unit, always `rad`.
    pub angle_unit: String,
    /// A written description of the coordinate frame.
    pub frame: String,
    /// The derived mechanism geometry.
    pub design: SceneDesign,
    /// The L1 atomistic layer.
    pub atomistic: AtomisticLayer,
    /// The L2 device layer.
    pub device: DeviceLayer,
    /// The coarse handoff layer.
    pub coarse: CoarseLayer,
}

/// Builds the three-scale scene from an assembly and its generated parts.
///
/// The `set` must be the [`PlanetarySet`] that produced the records in
/// `assembly.library`. The function reads atom positions and atom ranges from
/// `set`, and it reads poses, joints, and couplings from `assembly`.
pub fn build_scene(assembly: &PlanetaryAssembly, set: &PlanetarySet) -> Result<Scene, SceneError> {
    let body_count = assembly.system.body_count();
    if assembly.roles.len() != body_count {
        return Err(SceneError::BodyCountMismatch {
            roles: assembly.roles.len(),
            bodies: body_count,
        });
    }
    if assembly.design != set.design {
        return Err(SceneError::DesignMismatch);
    }

    let topology = &set.part.topology;
    let atom_count = topology.atom_count();

    let layers = set
        .part
        .metadata
        .get("layers")
        .and_then(|value| value.parse().ok())
        .unwrap_or(1);
    let layer_spacing_m = set
        .part
        .metadata
        .get("layer_spacing_m")
        .and_then(|value| value.parse().ok())
        .unwrap_or(0.0);
    let thickness_m = set
        .part
        .metadata
        .get("thickness_m")
        .and_then(|value| value.parse().ok())
        .unwrap_or(0.0);

    let atomistic = build_atomistic(assembly, set, atom_count)?;
    let device = build_device(assembly, body_count)?;
    let coarse = build_coarse(assembly, set, atom_count)?;

    Ok(Scene {
        schema: SCENE_SCHEMA.to_owned(),
        version: SCENE_VERSION,
        length_unit: "m".to_owned(),
        angle_unit: "rad".to_owned(),
        frame: "right-handed, z is the gear axis, x-y is the gear plane, \
                metres and radians, quaternion [w, x, y, z]"
            .to_owned(),
        design: SceneDesign {
            module_m: set.design.module_m(),
            sun_teeth: set.design.sun_teeth(),
            planet_teeth: set.design.planet_teeth(),
            ring_teeth: set.design.ring_teeth(),
            planet_count: set.design.planet_count(),
            gear_ratio: set.design.gear_ratio(),
            sun_pitch_radius_m: set.design.sun_pitch_radius_m(),
            planet_pitch_radius_m: set.design.planet_pitch_radius_m(),
            ring_pitch_radius_m: set.design.ring_pitch_radius_m(),
            carrier_radius_m: set.design.carrier_radius_m(),
            layers,
            layer_spacing_m,
            thickness_m,
        },
        atomistic,
        device,
        coarse,
    })
}

fn build_atomistic(
    assembly: &PlanetaryAssembly,
    set: &PlanetarySet,
    atom_count: usize,
) -> Result<AtomisticLayer, SceneError> {
    let topology = &set.part.topology;
    let mut covered = vec![false; atom_count];
    let mut atoms = Vec::with_capacity(atom_count);
    for (body_index, role) in assembly.roles.iter().enumerate() {
        let Some(range) = atom_range(role, set)? else {
            continue;
        };
        if range.end > atom_count {
            return Err(SceneError::AtomRangeOutOfBounds {
                role: role_label(*role).to_owned(),
                start: range.start,
                end: range.end,
                atom_count,
            });
        }
        let part_id = role_part_id(*role).unwrap_or_default();
        let name = role_name(*role);
        for index in range {
            let slot = covered
                .get_mut(index)
                .ok_or(SceneError::MissingAtom { index })?;
            if *slot {
                return Err(SceneError::DuplicateAtom { index });
            }
            *slot = true;
            let atom = topology
                .atom(index)
                .ok_or(SceneError::MissingAtom { index })?;
            if atom.position_m.iter().any(|value| !value.is_finite()) {
                return Err(SceneError::NonFiniteAtom {
                    index,
                    position_m: atom.position_m,
                });
            }
            let element = atom
                .element
                .symbol()
                .map(str::to_owned)
                .unwrap_or_else(|| format!("Z{}", atom.element.atomic_number()));
            atoms.push(SceneAtom {
                element,
                atomic_number: atom.element.atomic_number(),
                position_m: atom.position_m,
                body: body_index,
                name: name.clone(),
                part_id: part_id.clone(),
            });
        }
    }
    if let Some(index) = covered.iter().position(|assigned| !assigned) {
        return Err(SceneError::UnassignedAtom { index });
    }
    Ok(AtomisticLayer { atom_count, atoms })
}

fn build_device(
    assembly: &PlanetaryAssembly,
    body_count: usize,
) -> Result<DeviceLayer, SceneError> {
    let mut bodies = Vec::with_capacity(body_count);
    for (index, role) in assembly.roles.iter().enumerate() {
        let body = assembly
            .system
            .body(index)
            .ok_or(SceneError::MissingBody { body: index })?;
        let position_m = body.position_m();
        if position_m.iter().any(|value| !value.is_finite()) {
            return Err(SceneError::NonFiniteBody {
                body: index,
                position_m,
            });
        }
        let orientation = body.orientation();
        let inertia = body.inertia_kg_m2();
        bodies.push(SceneBody {
            index,
            name: role_name(*role),
            role: role_label(*role).to_owned(),
            part_id: role_part_id(*role),
            position_m,
            orientation_wxyz: [orientation.w, orientation.x, orientation.y, orientation.z],
            mass_kg: body.mass_kg(),
            inertia_diagonal_kg_m2: [inertia[0][0], inertia[1][1], inertia[2][2]],
            fixed: body.is_fixed(),
            linear_velocity_m_per_s: body.linear_velocity_m_per_s(),
            angular_velocity_rad_per_s: body.angular_velocity_rad_per_s(),
        });
    }

    let all_bodies = assembly.system.bodies();
    let joint_names = assembly_joint_names(assembly);

    let mut joints = Vec::with_capacity(
        assembly.system.revolute_joints().len() + assembly.system.prismatic_joints().len(),
    );
    for (index, joint) in assembly.system.revolute_joints().iter().enumerate() {
        let anchor_world_m = joint.anchor_a_world_m(all_bodies)?;
        let axis_world = joint.axis_a_world(all_bodies)?;
        joints.push(SceneJoint {
            name: joint_names
                .get(index)
                .and_then(Clone::clone)
                .unwrap_or_else(|| format!("revolute_{index}")),
            kind: "revolute".to_owned(),
            body_a: joint.body_a(),
            body_b: joint.body_b(),
            anchor_world_m: Some(anchor_world_m),
            axis_world,
            axis_a_body: joint.axis_a_body(),
            axis_b_body: Some(joint.axis_b_body()),
        });
    }
    for (index, joint) in assembly.system.prismatic_joints().iter().enumerate() {
        joints.push(SceneJoint {
            name: format!("prismatic_{index}"),
            kind: "prismatic".to_owned(),
            body_a: joint.body_a(),
            body_b: joint.body_b(),
            anchor_world_m: None,
            axis_world: joint.axis_a_world(all_bodies)?,
            axis_a_body: joint.axis_a_body(),
            axis_b_body: None,
        });
    }

    let gear_couplings = assembly
        .system
        .gear_couplings()
        .iter()
        .map(|coupling| SceneGearCoupling {
            joint_a: coupling.joint_a(),
            joint_b: coupling.joint_b(),
            radius_a_m: coupling.radius_a_m(),
            radius_b_m: coupling.radius_b_m(),
        })
        .collect();

    let gear_constraints = assembly
        .system
        .gear_constraints()
        .iter()
        .map(|constraint| SceneGearConstraint {
            terms: constraint
                .terms()
                .iter()
                .map(|term| SceneGearTerm {
                    joint: term.joint(),
                    coefficient_m: term.coefficient_m(),
                })
                .collect(),
        })
        .collect();

    Ok(DeviceLayer {
        body_count,
        bodies,
        joints,
        gear_couplings,
        gear_constraints,
    })
}

fn build_coarse(
    assembly: &PlanetaryAssembly,
    set: &PlanetarySet,
    atom_count: usize,
) -> Result<CoarseLayer, SceneError> {
    let mut bodies = Vec::with_capacity(assembly.roles.len());
    for (index, role) in assembly.roles.iter().enumerate() {
        let body = assembly
            .system
            .body(index)
            .ok_or(SceneError::MissingBody { body: index })?;
        let center_m = body.position_m();
        let pitch_circle = gear_pitch_circle(*role, &set.design, center_m);
        let bounding_cylinder = match atom_range(role, set)? {
            Some(range) if range.end <= atom_count => {
                bounding_cylinder(&set.part.topology, range, center_m)
            }
            Some(range) => {
                return Err(SceneError::AtomRangeOutOfBounds {
                    role: role_label(*role).to_owned(),
                    start: range.start,
                    end: range.end,
                    atom_count,
                })
            }
            None => None,
        };
        bodies.push(CoarseBody {
            index,
            name: role_name(*role),
            role: role_label(*role).to_owned(),
            pitch_circle,
            bounding_cylinder,
        });
    }
    Ok(CoarseLayer {
        body_count: assembly.roles.len(),
        bodies,
    })
}

fn bounding_cylinder(
    topology: &nanocad_model::Topology,
    range: Range<usize>,
    center_m: [f64; 3],
) -> Option<CoarseCylinder> {
    let mut min_z_m = f64::INFINITY;
    let mut max_z_m = f64::NEG_INFINITY;
    let mut max_radius_m: f64 = 0.0;
    let mut any = false;
    for index in range {
        let Some(position_m) = topology.position_m(index) else {
            continue;
        };
        any = true;
        min_z_m = min_z_m.min(position_m[2]);
        max_z_m = max_z_m.max(position_m[2]);
        let dx = position_m[0] - center_m[0];
        let dy = position_m[1] - center_m[1];
        max_radius_m = max_radius_m.max((dx * dx + dy * dy).sqrt());
    }
    if !any {
        return None;
    }
    Some(CoarseCylinder {
        center_m: [center_m[0], center_m[1], 0.5 * (min_z_m + max_z_m)],
        axis: [0.0, 0.0, 1.0],
        radius_m: max_radius_m,
        half_length_m: 0.5 * (max_z_m - min_z_m),
    })
}

fn gear_pitch_circle(
    role: BodyRole,
    design: &nanocad_parts::PlanetaryDesign,
    center_m: [f64; 3],
) -> Option<CoarseCircle> {
    let radius_m = match role {
        BodyRole::Sun => design.sun_pitch_radius_m(),
        BodyRole::Planet(_) => design.planet_pitch_radius_m(),
        BodyRole::Ring => design.ring_pitch_radius_m(),
        BodyRole::Ground | BodyRole::Carrier => return None,
    };
    Some(CoarseCircle {
        center_m,
        axis: [0.0, 0.0, 1.0],
        radius_m,
    })
}

fn atom_range(role: &BodyRole, set: &PlanetarySet) -> Result<Option<Range<usize>>, SceneError> {
    Ok(match role {
        BodyRole::Ground => None,
        BodyRole::Sun => Some(set.sun_atoms.clone()),
        BodyRole::Planet(planet) => Some(
            set.planet_atoms
                .get(*planet)
                .cloned()
                .ok_or(SceneError::MissingPlanet { planet: *planet })?,
        ),
        BodyRole::Ring => Some(set.ring_atoms.clone()),
        BodyRole::Carrier => Some(set.carrier_atoms.clone()),
    })
}

fn assembly_joint_names(assembly: &PlanetaryAssembly) -> Vec<Option<String>> {
    let mut names = vec![None; assembly.system.revolute_joints().len()];
    for role in &assembly.roles {
        let joint = match role {
            BodyRole::Ground => continue,
            BodyRole::Sun => assembly.sun_joint,
            BodyRole::Carrier => assembly.carrier_joint,
            BodyRole::Ring => assembly.ring_joint,
            BodyRole::Planet(planet) => match assembly.planet_joints.get(*planet) {
                Some(joint) => *joint,
                None => continue,
            },
        };
        if let Some(slot) = names.get_mut(joint) {
            *slot = Some(format!("{}_joint", role_name(*role)));
        }
    }
    names
}

fn role_name(role: BodyRole) -> String {
    match role {
        BodyRole::Ground => "ground".to_owned(),
        BodyRole::Sun => "sun".to_owned(),
        BodyRole::Planet(planet) => format!("planet_{planet}"),
        BodyRole::Ring => "ring".to_owned(),
        BodyRole::Carrier => "carrier".to_owned(),
    }
}

fn role_label(role: BodyRole) -> &'static str {
    match role {
        BodyRole::Ground => "ground",
        BodyRole::Sun => "sun",
        BodyRole::Planet(_) => "planet",
        BodyRole::Ring => "ring",
        BodyRole::Carrier => "carrier",
    }
}

fn role_part_id(role: BodyRole) -> Option<String> {
    match role {
        BodyRole::Ground => None,
        BodyRole::Sun => Some(SUN_PART_ID.to_owned()),
        BodyRole::Planet(planet) => Some(format!("{PLANET_PART_ID_PREFIX}{planet}")),
        BodyRole::Ring => Some(RING_PART_ID.to_owned()),
        BodyRole::Carrier => Some(CARRIER_PART_ID.to_owned()),
    }
}

/// Serializes a scene to compact JSON.
pub fn scene_to_json(scene: &Scene) -> Result<String, SceneError> {
    serde_json::to_string(scene).map_err(|error| SceneError::Json(error.to_string()))
}

/// Serializes a scene to indented JSON.
pub fn scene_to_json_pretty(scene: &Scene) -> Result<String, SceneError> {
    serde_json::to_string_pretty(scene).map_err(|error| SceneError::Json(error.to_string()))
}

/// Parses a scene from JSON.
pub fn scene_from_json(json: &str) -> Result<Scene, SceneError> {
    serde_json::from_str(json).map_err(|error| SceneError::Json(error.to_string()))
}

/// Writes a scene to a JSON file. `pretty` selects indented output.
pub fn write_scene_json(
    scene: &Scene,
    path: impl AsRef<Path>,
    pretty: bool,
) -> Result<(), SceneError> {
    let text = if pretty {
        scene_to_json_pretty(scene)?
    } else {
        scene_to_json(scene)?
    };
    std::fs::write(path, text).map_err(|error| SceneError::Io(error.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::assemble_planetary;
    use nanocad_parts::{ParameterSet, PlanetaryGenerator};

    fn default_parameters() -> ParameterSet {
        ParameterSet::new().with("planet_count", 3.0)
    }

    fn fixture() -> (PlanetaryAssembly, PlanetarySet) {
        let parameters = default_parameters();
        let assembly = assemble_planetary(&parameters).expect("assemble");
        let set = PlanetaryGenerator.build(&parameters).expect("generate");
        (assembly, set)
    }

    fn scene() -> Scene {
        let (assembly, set) = fixture();
        build_scene(&assembly, &set).expect("scene")
    }

    #[test]
    fn layer_counts_match_the_default_design() {
        let (assembly, set) = fixture();
        let scene = build_scene(&assembly, &set).expect("scene");
        assert_eq!(scene.schema, SCENE_SCHEMA);
        assert_eq!(scene.version, SCENE_VERSION);
        assert_eq!(scene.atomistic.atom_count, set.part.atom_count());
        assert_eq!(scene.atomistic.atoms.len(), set.part.atom_count());
        assert_eq!(scene.device.body_count, assembly.system.body_count());
        assert_eq!(scene.device.body_count, 7);
        assert_eq!(scene.device.bodies.len(), 7);
        assert_eq!(scene.device.joints.len(), 6);
        assert_eq!(scene.device.gear_constraints.len(), 6);
        assert_eq!(scene.device.gear_couplings.len(), 0);
        assert_eq!(scene.coarse.body_count, 7);
        assert_eq!(scene.coarse.bodies.len(), 7);
        assert_eq!(scene.design.planet_count, 3);
        assert_eq!(scene.design.ring_teeth, 30);
        assert_eq!(scene.design.gear_ratio, 3.5);
    }

    #[test]
    fn every_position_is_finite() {
        let scene = scene();
        for atom in &scene.atomistic.atoms {
            assert!(
                atom.position_m.iter().all(|value| value.is_finite()),
                "atom {} is not finite",
                atom.name
            );
        }
        for body in &scene.device.bodies {
            assert!(
                body.position_m.iter().all(|value| value.is_finite()),
                "body {} is not finite",
                body.name
            );
            assert!(body.orientation_wxyz.iter().all(|value| value.is_finite()));
            assert!(body.mass_kg.is_finite() && body.mass_kg > 0.0);
        }
        for body in &scene.coarse.bodies {
            if let Some(cylinder) = &body.bounding_cylinder {
                assert!(cylinder.center_m.iter().all(|value| value.is_finite()));
                assert!(cylinder.radius_m.is_finite() && cylinder.radius_m >= 0.0);
                assert!(cylinder.half_length_m.is_finite() && cylinder.half_length_m >= 0.0);
            }
        }
    }

    #[test]
    fn the_atomistic_count_matches_the_parts() {
        let (assembly, set) = fixture();
        let scene = build_scene(&assembly, &set).expect("scene");
        let covered = set.sun_atoms.len()
            + set
                .planet_atoms
                .iter()
                .map(std::ops::Range::len)
                .sum::<usize>()
            + set.ring_atoms.len()
            + set.carrier_atoms.len();
        assert_eq!(covered, set.part.atom_count());
        assert_eq!(scene.atomistic.atom_count, covered);
        assert_eq!(scene.atomistic.atoms.len(), set.part.atom_count());
        for atom in &scene.atomistic.atoms {
            assert!(atom.body < scene.device.body_count);
            assert!(!atom.name.is_empty());
        }
    }

    #[test]
    fn the_atomistic_positions_match_the_topology() {
        let (assembly, set) = fixture();
        let scene = build_scene(&assembly, &set).expect("scene");
        for (index, atom) in scene.atomistic.atoms.iter().enumerate() {
            let topology = set.part.topology.atom(index).expect("atom");
            assert_eq!(atom.position_m, topology.position_m);
            assert_eq!(atom.atomic_number, topology.element.atomic_number());
        }
    }

    #[test]
    fn the_scene_round_trips_through_json() {
        let scene = scene();
        let compact = scene_to_json(&scene).expect("compact");
        assert_eq!(scene_from_json(&compact).expect("parse"), scene);
        let pretty = scene_to_json_pretty(&scene).expect("pretty");
        assert_eq!(scene_from_json(&pretty).expect("parse"), scene);
        assert!(pretty.contains('\n'));
        assert!(!compact.contains('\n'));
    }

    #[test]
    fn the_export_is_deterministic() {
        let (assembly, set) = fixture();
        let first = build_scene(&assembly, &set).expect("scene");
        let second = build_scene(&assembly, &set).expect("scene");
        assert_eq!(first, second);
        assert_eq!(
            scene_to_json(&first).expect("json"),
            scene_to_json(&second).expect("json")
        );
    }

    #[test]
    fn each_gear_has_a_pitch_circle_and_each_body_an_index() {
        let scene = scene();
        let teeth = scene.design.sun_teeth;
        assert!(teeth > 0);
        for body in &scene.coarse.bodies {
            assert!(body.index < scene.device.body_count);
            match body.role.as_str() {
                "sun" | "planet" | "ring" => {
                    let circle = body.pitch_circle.as_ref().expect("pitch circle");
                    assert!(circle.radius_m > 0.0);
                    assert_eq!(circle.axis, [0.0, 0.0, 1.0]);
                }
                _ => assert!(body.pitch_circle.is_none()),
            }
        }
        let ground = scene
            .coarse
            .bodies
            .iter()
            .find(|body| body.role == "ground")
            .expect("ground");
        assert!(ground.bounding_cylinder.is_none());
    }

    #[test]
    fn the_scene_writes_to_a_file() {
        let scene = scene();
        let path = std::env::temp_dir().join(format!(
            "nanocad-scene-{}-{:?}.json",
            std::process::id(),
            std::thread::current().id()
        ));
        write_scene_json(&scene, &path, true).expect("write");
        let text = std::fs::read_to_string(&path).expect("read");
        assert_eq!(scene_from_json(&text).expect("parse"), scene);
        std::fs::remove_file(&path).expect("remove");
    }
}
