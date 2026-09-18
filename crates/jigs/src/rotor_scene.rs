//! A three-scale scene for the sorting rotor.
//!
//! The scene uses the same schema as the planetary gearbox, so one viewer
//! shows both. The two mechanisms differ, so the scene states the difference
//! in its data instead of in a new schema.
//!
//! - The housing is body 0 and it is fixed. The rotor is body 1. The drive
//!   shaft turns with the rotor, because a key locks the shaft to the rotor's
//!   keyway, so the shaft atoms belong to body 1 and need no joint of their own.
//! - Bodies 2 through 13 are the twelve ejection rods, one for each pocket.
//!   Each rod is captive in a radial ejection bore in the rotor, so its
//!   prismatic joint names the rotor as its first body and the rod slides along
//!   that pocket's radius.
//! - Body 14 is the cam plate. It is fixed to the housing and shares the rotor
//!   axis. The plate lies below the rotor, and its face laps under the housing
//!   ring, so the housing holds the plate from below and around its rim. One
//!   groove is cut into the plate face. The groove dwells at a base radius for
//!   most of the turn, then rises over a short ramp at the outlet, so a rod
//!   stays retracted until its pocket reaches the outlet. The two walls of the
//!   groove push the rod out and pull it back, with no spring.
//! - Each rod carries a follower pin. The pin hangs from the rod into the
//!   groove, so the groove walls drive the rod.
//! - The design block is zero. Its fields describe a gear set, and a rotor has
//!   no gear set. A viewer reads a zero module as "not a gearbox".
//! - The rotor turns about the world `z` axis through the origin, because the
//!   housing chamber port and the rotor axis port share that axis.
//!
//! The device layer carries the joints that a viewer needs to animate the
//! device: the rotor turns about `z`, and each rod slides along its pocket
//! radius. It carries no gear coupling and no gear constraint.
//!
//! # Scope limit
//!
//! The scene does not include a binding site, a guest molecule, a solvent or
//! an ion, so it does not show molecular selectivity. The profiled groove gives
//! one out stroke and one return stroke for each rotor turn, and it holds the
//! stroke to a short ramp, so the rod pushes the guest out only as the pocket
//! meets the housing outlet. It is a radial key and not a tuned motion law, so
//! the rod acceleration is not designed. The pin slides in the groove; a roller
//! follower would roll. The scene omits the contact force between the pin and
//! the groove, the torque on the drive shaft, and the machine that turns the
//! shaft. The work of one ejection stays a measured barrier in
//! [`nanocad_meter`], not a rod force.

use nanocad_model::Part;
use nanocad_parts::{
    place, CamPlateGenerator, DriveShaftGenerator, EjectionRodGenerator, FollowerPinGenerator,
    ParameterSet, PartError, PartGenerator, RotorHousingGenerator, SortingRotorGenerator,
};

use crate::scene::{
    AtomisticLayer, CoarseBody, CoarseCylinder, CoarseLayer, DeviceLayer, Scene, SceneAtom,
    SceneBody, SceneError, SceneJoint, SCENE_SCHEMA, SCENE_VERSION,
};

/// The device body index of the fixed housing.
pub const HOUSING_BODY: usize = 0;
/// The device body index of the rotating rotor.
pub const ROTOR_BODY: usize = 1;
/// The device body index of the first ejection rod.
pub const ROD_BODY_FIRST: usize = 2;
/// The number of ejection rods, one for each pocket.
pub const ROD_COUNT: usize = 12;
/// The device body index of the fixed cam plate.
pub const CAM_BODY: usize = ROD_BODY_FIRST + ROD_COUNT;
/// The body role of the fixed housing.
pub const HOUSING_ROLE: &str = "housing";
/// The body role of the rotating rotor.
pub const ROTOR_ROLE: &str = "rotor";
/// The body role of an ejection rod.
pub const ROD_ROLE: &str = "ejection_rod";
/// The body role of the fixed cam plate.
pub const CAM_ROLE: &str = "cam_plate";
/// The part id of the drive shaft, whose atoms ride on the rotor body.
pub const SHAFT_PART_ID: &str = "drive_shaft";
/// The part id of a follower pin, whose atoms ride on a rod body.
pub const PIN_PART_ID: &str = "follower_pin";

/// A written description of the rotor scene frame.
pub const ROTOR_FRAME: &str =
    "right-handed Cartesian, z is the rotor axis, the rotor turns about z through the origin";

fn scene_atoms(part: &Part, body: usize, name: &str, part_id: &str) -> Vec<SceneAtom> {
    scene_atoms_mapped(part, body, name, part_id, |point_m| point_m)
}

/// Collects the atoms of one part, with a map from the part frame to the scene.
fn scene_atoms_mapped(
    part: &Part,
    body: usize,
    name: &str,
    part_id: &str,
    map: impl Fn([f64; 3]) -> [f64; 3],
) -> Vec<SceneAtom> {
    (0..part.atom_count())
        .filter_map(|index| {
            let element = part.topology.element(index)?;
            let position_m = part.topology.position_m(index)?;
            Some(SceneAtom {
                element: element
                    .symbol()
                    .map(str::to_string)
                    .unwrap_or_else(|| format!("Z{}", element.atomic_number())),
                atomic_number: element.atomic_number(),
                position_m: map(position_m),
                body,
                name: name.to_string(),
                part_id: part_id.to_string(),
            })
        })
        .collect()
}

/// The largest radius and the half height of the atoms of one part.
fn extent_m(part: &Part) -> Option<(f64, f64)> {
    let mut radius_m: f64 = 0.0;
    let mut half_height_m: f64 = 0.0;
    for index in 0..part.atom_count() {
        let position_m = part.topology.position_m(index)?;
        let radius = (position_m[0] * position_m[0] + position_m[1] * position_m[1]).sqrt();
        radius_m = radius_m.max(radius);
        half_height_m = half_height_m.max(position_m[2].abs());
    }
    Some((radius_m, half_height_m))
}

/// The default value of one parameter of a generator.
fn default_m(generator: &impl PartGenerator, name: &str) -> f64 {
    generator
        .resolve(&ParameterSet::new())
        .ok()
        .and_then(|values| values.get(name))
        .unwrap_or(0.0)
}

/// A scene body on the rotor axis, with no motion of its own.
fn scene_body(index: usize, role: &str, part_id: &str, fixed: bool) -> SceneBody {
    SceneBody {
        index,
        name: role.to_string(),
        role: role.to_string(),
        part_id: Some(part_id.to_string()),
        position_m: [0.0, 0.0, 0.0],
        orientation_wxyz: [1.0, 0.0, 0.0, 0.0],
        mass_kg: 0.0,
        inertia_diagonal_kg_m2: [0.0, 0.0, 0.0],
        fixed,
        linear_velocity_m_per_s: [0.0, 0.0, 0.0],
        angular_velocity_rad_per_s: [0.0, 0.0, 0.0],
    }
}

/// Builds the sorting rotor scene.
///
/// The rotor, the housing, the cam plate and the drive shaft use their default
/// parameters. The rotor mates to the housing through the rotor axis port and
/// the housing chamber port. The cam plate laps under the housing ring, so its
/// top face is flush with the housing lower face. Each ejection rod is built
/// with the shaft length that reaches from the cam groove out to the inner wall
/// of its pocket. Each follower pin hangs from the inner end of its rod into the
/// groove. The drive shaft shares the rotor body, because the key locks the two.
///
/// # Errors
///
/// Returns [`SceneError::Part`] when a generator rejects its defaults, when the
/// rod is too short to reach its pocket, when the cam groove or the cam bore
/// does not clear the rotor bore, or when a pin does not reach the groove.
/// Returns [`SceneError::MissingAtom`] when an atom has no element or position.
pub fn build_rotor_scene() -> Result<Scene, SceneError> {
    let rotor = SortingRotorGenerator.generate_with_defaults()?;
    let housing = RotorHousingGenerator.generate_with_defaults()?;
    let cam = CamPlateGenerator.generate_with_defaults()?;
    let shaft = DriveShaftGenerator.generate_with_defaults()?;
    let pin = FollowerPinGenerator.generate_with_defaults()?;
    let rotor = place(
        rotor,
        &SortingRotorGenerator.axis_port(),
        &RotorHousingGenerator.chamber_port(),
        0.0,
    )
    .transformed_part();

    let pocket_count = default_m(&SortingRotorGenerator, "pocket_count").round() as usize;
    let pocket_orbit_m = default_m(&SortingRotorGenerator, "pocket_orbit_m");
    let pocket_radius_m = default_m(&SortingRotorGenerator, "pocket_radius_m");
    let pocket_inner_m = pocket_orbit_m - pocket_radius_m;
    let rotor_bore_m = default_m(&SortingRotorGenerator, "bore_radius_m");
    let rod_radius_m = default_m(&EjectionRodGenerator, "shaft_radius_m");
    let tip_length_m = default_m(&EjectionRodGenerator, "tip_length_m");
    let groove_base_radius_m = CamPlateGenerator.groove_base_radius_m();
    let groove_rise_m = CamPlateGenerator.groove_rise_m();
    let groove_depth_m = CamPlateGenerator.groove_depth_m();
    let groove_width_m = CamPlateGenerator.groove_width_m();
    let cam_retract_m = groove_base_radius_m;
    if cam_retract_m <= rotor_bore_m + rod_radius_m {
        return Err(SceneError::Part(format!(
            "the cam groove retracts to {cam_retract_m} m, which reaches the rotor bore \
             {rotor_bore_m} m plus the rod radius {rod_radius_m} m"
        )));
    }
    let disc_radius_m = default_m(&SortingRotorGenerator, "disc_radius_m");
    let rod_length_m = pocket_inner_m - cam_retract_m;
    let rod_tip_m = pocket_inner_m + groove_rise_m;
    if rod_tip_m > disc_radius_m + 1.0e-15 {
        return Err(SceneError::Part(format!(
            "the rod reaches {rod_tip_m} m at full stroke, past the rotor rim {disc_radius_m} m, \
             so it would cross the housing"
        )));
    }
    let shaft_length_m = rod_length_m - tip_length_m;
    if shaft_length_m <= 0.0 {
        return Err(SceneError::Part(
            "the cam groove reaches the pocket inner wall, so no rod fits".to_string(),
        ));
    }
    let rod = EjectionRodGenerator
        .generate(&ParameterSet::new().with("shaft_length_m", shaft_length_m))?;
    let half_shaft_m = 0.5 * shaft_length_m;
    let radial_offset_m = cam_retract_m + half_shaft_m;

    let shaft_radius_m = default_m(&DriveShaftGenerator, "shaft_radius_m");
    if shaft_radius_m >= rotor_bore_m || shaft_radius_m >= CamPlateGenerator.inner_radius_m() {
        return Err(SceneError::Part(format!(
            "the drive shaft radius {shaft_radius_m} m does not clear the rotor bore \
             {rotor_bore_m} m and the cam plate bore"
        )));
    }

    let housing_half_thickness_m = default_m(&RotorHousingGenerator, "thickness_m") / 2.0;
    let rotor_half_thickness_m = default_m(&SortingRotorGenerator, "thickness_m") / 2.0;
    let cam_top_m = -housing_half_thickness_m;
    let cam_offset_m = cam_top_m - 0.5 * CamPlateGenerator.thickness_m();
    let cam_bottom_m = cam_top_m - groove_depth_m;
    let pin_length_m = FollowerPinGenerator.length_m();
    let pin_bottom_m = -rotor_half_thickness_m - pin_length_m;
    if pin_bottom_m >= cam_top_m || pin_bottom_m <= cam_bottom_m {
        return Err(SceneError::Part(format!(
            "the follower pin ends at {pin_bottom_m} m, outside the groove between \
             {cam_bottom_m} m and {cam_top_m} m"
        )));
    }
    let pin_head_radius_m = FollowerPinGenerator.head_radius_m();
    if 2.0 * pin_head_radius_m >= groove_width_m {
        return Err(SceneError::Part(format!(
            "the follower head diameter {} m does not fit the groove width {groove_width_m} m",
            2.0 * pin_head_radius_m
        )));
    }
    let housing_atoms = scene_atoms(&housing, HOUSING_BODY, HOUSING_ROLE, "rotor_housing");
    let mut rotor_atoms = scene_atoms(&rotor, ROTOR_BODY, ROTOR_ROLE, "sorting_rotor");
    rotor_atoms.extend(scene_atoms(&shaft, ROTOR_BODY, ROTOR_ROLE, SHAFT_PART_ID));
    let cam_atoms = scene_atoms_mapped(&cam, CAM_BODY, CAM_ROLE, "cam_plate", |point_m| {
        [point_m[0], point_m[1], point_m[2] + cam_offset_m]
    });

    let mut atoms = housing_atoms;
    atoms.extend(rotor_atoms.iter().cloned());

    let mut rod_bodies = Vec::with_capacity(ROD_COUNT);
    let mut rod_coarse = Vec::with_capacity(ROD_COUNT);
    let mut rod_joints = Vec::with_capacity(ROD_COUNT);
    for index in 0..ROD_COUNT {
        let angle_rad = 2.0 * std::f64::consts::PI * index as f64 / pocket_count.max(1) as f64;
        let (sin_rad, cos_rad) = angle_rad.sin_cos();
        let radial_m = [cos_rad, sin_rad, 0.0];
        let pocket_center_m =
            SortingRotorGenerator::pocket_center_m(index, pocket_count, pocket_orbit_m);
        let body = ROD_BODY_FIRST + index;
        let rod_atoms = scene_atoms_mapped(&rod, body, ROD_ROLE, "ejection_rod", |point_m| {
            [
                point_m[2] * cos_rad - point_m[1] * sin_rad + radial_offset_m * cos_rad,
                point_m[2] * sin_rad + point_m[1] * cos_rad + radial_offset_m * sin_rad,
                -point_m[0],
            ]
        });
        let pin_atoms = scene_atoms_mapped(&pin, body, PIN_PART_ID, PIN_PART_ID, |point_m| {
            [
                cam_retract_m * cos_rad + point_m[0],
                cam_retract_m * sin_rad + point_m[1],
                -rotor_half_thickness_m + point_m[2],
            ]
        });
        let mut rod_atom_set = rod_atoms;
        rod_atom_set.extend(pin_atoms);
        if let Some(cylinder) = local_cylinder(&rod_atom_set, radial_m) {
            rod_coarse.push(CoarseBody {
                index: body,
                name: ROD_ROLE.to_string(),
                role: ROD_ROLE.to_string(),
                pitch_circle: None,
                bounding_cylinder: Some(cylinder),
            });
        }
        atoms.extend(rod_atom_set);
        rod_bodies.push(scene_body(body, ROD_ROLE, "ejection_rod", false));
        rod_joints.push(SceneJoint {
            name: format!("ejection_rod_{index}"),
            kind: "prismatic".to_string(),
            body_a: ROTOR_BODY,
            body_b: body,
            anchor_world_m: Some(pocket_center_m),
            axis_world: radial_m,
            axis_a_body: radial_m,
            axis_b_body: Some(radial_m),
        });
    }
    atoms.extend(cam_atoms.iter().cloned());
    if atoms.is_empty() {
        return Err(SceneError::MissingAtom { index: 0 });
    }

    let axis = [0.0, 0.0, 1.0];
    let mut bodies = vec![
        scene_body(HOUSING_BODY, HOUSING_ROLE, "rotor_housing", true),
        scene_body(ROTOR_BODY, ROTOR_ROLE, "sorting_rotor", false),
    ];
    bodies.extend(rod_bodies);
    bodies.push(scene_body(CAM_BODY, CAM_ROLE, "cam_plate", true));

    let mut joints = Vec::with_capacity(1 + ROD_COUNT);
    joints.push(SceneJoint {
        name: "rotor_axis".to_string(),
        kind: "revolute".to_string(),
        body_a: HOUSING_BODY,
        body_b: ROTOR_BODY,
        anchor_world_m: Some([0.0, 0.0, 0.0]),
        axis_world: axis,
        axis_a_body: axis,
        axis_b_body: Some(axis),
    });
    joints.extend(rod_joints);

    let mut coarse = vec![
        CoarseBody {
            index: HOUSING_BODY,
            name: HOUSING_ROLE.to_string(),
            role: HOUSING_ROLE.to_string(),
            pitch_circle: None,
            bounding_cylinder: bounding_cylinder(&housing, axis, housing_half_thickness_m),
        },
        CoarseBody {
            index: ROTOR_BODY,
            name: ROTOR_ROLE.to_string(),
            role: ROTOR_ROLE.to_string(),
            pitch_circle: None,
            bounding_cylinder: local_cylinder(&rotor_atoms, axis),
        },
    ];
    coarse.extend(rod_coarse);
    coarse.push(CoarseBody {
        index: CAM_BODY,
        name: CAM_ROLE.to_string(),
        role: CAM_ROLE.to_string(),
        pitch_circle: None,
        bounding_cylinder: local_cylinder(&cam_atoms, axis),
    });

    let body_count = bodies.len();
    let coarse_count = coarse.len();
    Ok(Scene {
        schema: SCENE_SCHEMA.to_string(),
        version: SCENE_VERSION,
        length_unit: "m".to_string(),
        angle_unit: "rad".to_string(),
        frame: ROTOR_FRAME.to_string(),
        design: zero_design(),
        atomistic: AtomisticLayer {
            atom_count: atoms.len(),
            atoms,
        },
        device: DeviceLayer {
            body_count,
            bodies,
            joints,
            gear_couplings: Vec::new(),
            gear_constraints: Vec::new(),
        },
        coarse: CoarseLayer {
            body_count: coarse_count,
            bodies: coarse,
        },
    })
}

fn zero_design() -> crate::scene::SceneDesign {
    crate::scene::SceneDesign {
        module_m: 0.0,
        sun_teeth: 0,
        planet_teeth: 0,
        ring_teeth: 0,
        planet_count: 0,
        gear_ratio: 0.0,
        sun_pitch_radius_m: 0.0,
        planet_pitch_radius_m: 0.0,
        ring_pitch_radius_m: 0.0,
        carrier_radius_m: 0.0,
        layers: 0,
        layer_spacing_m: 0.0,
        thickness_m: 0.0,
    }
}

/// The rotor scene reports a generator error.
impl From<PartError> for SceneError {
    fn from(error: PartError) -> Self {
        SceneError::Part(error.to_string())
    }
}

/// The bounding cylinder of the atoms of one body.
fn bounding_cylinder(part: &Part, axis: [f64; 3], half_length_m: f64) -> Option<CoarseCylinder> {
    let (radius_m, atom_half_height_m) = extent_m(part)?;
    Some(CoarseCylinder {
        center_m: [0.0, 0.0, 0.0],
        axis,
        radius_m,
        half_length_m: half_length_m.max(atom_half_height_m),
    })
}

/// The tight cylinder around a set of atoms, along a stated axis.
fn local_cylinder(atoms: &[SceneAtom], axis: [f64; 3]) -> Option<CoarseCylinder> {
    let mut along_min_m = f64::INFINITY;
    let mut along_max_m = f64::NEG_INFINITY;
    let mut radius_m: f64 = 0.0;
    for atom in atoms {
        let point_m = atom.position_m;
        let along_m = point_m[0] * axis[0] + point_m[1] * axis[1] + point_m[2] * axis[2];
        along_min_m = along_min_m.min(along_m);
        along_max_m = along_max_m.max(along_m);
        let square_m2 = point_m[0] * point_m[0] + point_m[1] * point_m[1] + point_m[2] * point_m[2];
        radius_m = radius_m.max((square_m2 - along_m * along_m).max(0.0).sqrt());
    }
    if !along_min_m.is_finite() {
        return None;
    }
    let center_along_m = 0.5 * (along_min_m + along_max_m);
    Some(CoarseCylinder {
        center_m: [
            center_along_m * axis[0],
            center_along_m * axis[1],
            center_along_m * axis[2],
        ],
        axis,
        radius_m,
        half_length_m: 0.5 * (along_max_m - along_min_m),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::OnceLock;

    fn scene() -> &'static Scene {
        static SCENE: OnceLock<Scene> = OnceLock::new();
        SCENE.get_or_init(|| build_rotor_scene().expect("the rotor scene builds"))
    }

    #[test]
    fn the_rotor_scene_uses_the_shared_schema() {
        assert_eq!(scene().schema, SCENE_SCHEMA);
        assert_eq!(scene().version, SCENE_VERSION);
        assert_eq!(scene().length_unit, "m");
        assert_eq!(scene().angle_unit, "rad");
    }

    #[test]
    fn the_design_block_marks_a_non_gearbox() {
        assert_eq!(scene().design.module_m, 0.0);
        assert_eq!(scene().design.layers, 0);
    }

    #[test]
    fn the_rotor_scene_has_fifteen_bodies_and_thirteen_joints() {
        assert_eq!(scene().device.body_count, 15);
        assert_eq!(scene().coarse.body_count, 15);
        assert_eq!(scene().device.joints.len(), 13);
        let joint = &scene().device.joints[0];
        assert_eq!(joint.kind, "revolute");
        assert_eq!(joint.body_a, HOUSING_BODY);
        assert_eq!(joint.body_b, ROTOR_BODY);
        assert_eq!(joint.axis_world, [0.0, 0.0, 1.0]);
        for (index, ejection) in scene().device.joints[1..].iter().enumerate() {
            assert_eq!(ejection.kind, "prismatic");
            assert_eq!(ejection.body_a, ROTOR_BODY, "the rotor holds rod {index}");
            assert_eq!(ejection.body_b, ROD_BODY_FIRST + index);
            assert_eq!(ejection.axis_world[2], 0.0, "the rod axis is radial");
            let radius = (ejection.axis_world[0].powi(2) + ejection.axis_world[1].powi(2)).sqrt();
            assert!((radius - 1.0).abs() < 1.0e-9, "the axis is radial");
            assert!(ejection.anchor_world_m.is_some());
        }
        assert!(scene().device.gear_couplings.is_empty());
        assert!(scene().device.gear_constraints.is_empty());
    }

    #[test]
    fn every_rod_is_free_and_slides_on_the_rotor() {
        let scene = scene();
        for index in 0..ROD_COUNT {
            let body = ROD_BODY_FIRST + index;
            let rod = &scene.device.bodies[body];
            assert!(!rod.fixed);
            assert_eq!(rod.role, ROD_ROLE);
            assert_eq!(rod.part_id.as_deref(), Some("ejection_rod"));
            assert_eq!(rod.index, body, "the body index field matches its place");
        }
        assert_eq!(scene.device.bodies[CAM_BODY].index, CAM_BODY);
        assert_eq!(scene.device.bodies[HOUSING_BODY].index, HOUSING_BODY);
        assert_eq!(scene.device.bodies[ROTOR_BODY].index, ROTOR_BODY);
    }

    #[test]
    fn the_cam_is_fixed_to_the_housing() {
        let scene = scene();
        let cam = &scene.device.bodies[CAM_BODY];
        assert!(cam.fixed, "the cam plate does not turn");
        assert_eq!(cam.role, CAM_ROLE);
        assert_eq!(cam.part_id.as_deref(), Some("cam_plate"));
        assert!(
            scene
                .device
                .joints
                .iter()
                .all(|joint| joint.body_a != CAM_BODY && joint.body_b != CAM_BODY),
            "the cam plate has no joint of its own"
        );
    }

    #[test]
    fn the_housing_is_fixed_and_the_rotor_is_free() {
        let bodies = &scene().device.bodies;
        assert!(bodies[HOUSING_BODY].fixed);
        assert!(!bodies[ROTOR_BODY].fixed);
        assert_eq!(bodies[ROTOR_BODY].role, ROTOR_ROLE);
        assert_eq!(bodies[HOUSING_BODY].role, HOUSING_ROLE);
        assert_eq!(bodies[ROTOR_BODY].orientation_wxyz, [1.0, 0.0, 0.0, 0.0]);
    }

    #[test]
    fn the_drive_shaft_turns_with_the_rotor() {
        let scene = scene();
        let shaft_atoms = scene
            .atomistic
            .atoms
            .iter()
            .filter(|atom| atom.part_id == SHAFT_PART_ID)
            .count();
        assert!(shaft_atoms > 0, "the scene holds no drive shaft atom");
        assert!(
            scene
                .atomistic
                .atoms
                .iter()
                .filter(|atom| atom.part_id == SHAFT_PART_ID)
                .all(|atom| atom.body == ROTOR_BODY),
            "a drive shaft atom is not on the rotor body"
        );
        assert!(
            scene
                .device
                .joints
                .iter()
                .filter(|joint| joint.body_b == ROTOR_BODY && joint.kind == "revolute")
                .count()
                == 1,
            "the rotor body carries more than one rotor-axis joint"
        );
    }

    #[test]
    fn every_rod_carries_a_follower_pin() {
        let scene = scene();
        for index in 0..ROD_COUNT {
            let body = ROD_BODY_FIRST + index;
            let pins = scene
                .atomistic
                .atoms
                .iter()
                .filter(|atom| atom.body == body && atom.part_id == PIN_PART_ID)
                .count();
            assert!(pins > 0, "rod {index} carries no follower pin");
        }
    }

    #[test]
    fn every_atom_belongs_to_one_of_the_fifteen_bodies() {
        let scene = scene();
        assert_eq!(scene.atomistic.atom_count, scene.atomistic.atoms.len());
        let mut counts = [0_usize; 15];
        for atom in &scene.atomistic.atoms {
            assert!(atom.body < 15, "atom on unexpected body {}", atom.body);
            counts[atom.body] += 1;
            assert!(
                matches!(atom.atomic_number, 1 | 6 | 7 | 8 | 9 | 16 | 17),
                "a part atom is a lattice carbon, a cap, or a wall group atom"
            );
        }
        assert_eq!(counts[HOUSING_BODY], 39848);
        assert_eq!(counts[ROTOR_BODY], 54844);
        assert_eq!(counts[CAM_BODY], 55173);
        assert_eq!(counts[ROD_BODY_FIRST], 666);
        for index in 0..ROD_COUNT {
            assert_eq!(
                counts[ROD_BODY_FIRST + index],
                counts[ROD_BODY_FIRST],
                "every rod has the same atoms"
            );
            assert!(counts[ROD_BODY_FIRST] > 0);
        }
        assert_eq!(counts.iter().sum::<usize>(), scene.atomistic.atom_count);
    }

    #[test]
    fn the_coarse_bodies_wrap_the_atoms_of_their_body() {
        let scene = scene();
        for coarse in &scene.coarse.bodies {
            let cylinder = coarse.bounding_cylinder.expect("a cylinder");
            for atom in scene
                .atomistic
                .atoms
                .iter()
                .filter(|atom| atom.body == coarse.index)
            {
                let point_m = atom.position_m;
                let relative_m = [
                    point_m[0] - cylinder.center_m[0],
                    point_m[1] - cylinder.center_m[1],
                    point_m[2] - cylinder.center_m[2],
                ];
                let along_relative_m = relative_m[0] * cylinder.axis[0]
                    + relative_m[1] * cylinder.axis[1]
                    + relative_m[2] * cylinder.axis[2];
                let square_m2 = relative_m[0] * relative_m[0]
                    + relative_m[1] * relative_m[1]
                    + relative_m[2] * relative_m[2];
                let across_m = (square_m2 - along_relative_m * along_relative_m)
                    .max(0.0)
                    .sqrt();
                assert!(
                    across_m <= cylinder.radius_m + 1.0e-18,
                    "an atom of body {} sits at {across_m:e} m across, past the radius {} m",
                    coarse.index,
                    cylinder.radius_m
                );
                assert!(
                    along_relative_m.abs() <= cylinder.half_length_m + 1.0e-9,
                    "an atom of body {} sits past the cylinder ends",
                    coarse.index
                );
            }
        }
        let rotor = scene.coarse.bodies[ROTOR_BODY]
            .bounding_cylinder
            .expect("a cylinder");
        let housing = scene.coarse.bodies[HOUSING_BODY]
            .bounding_cylinder
            .expect("a cylinder");
        assert!(
            housing.radius_m > rotor.radius_m,
            "the housing is not wider than the rotor"
        );
    }

    #[test]
    fn the_rotor_scene_round_trips_through_json() {
        let text = crate::scene::scene_to_json(scene()).expect("the scene serializes");
        let parsed = crate::scene::scene_from_json(&text).expect("the scene parses");
        assert_eq!(&parsed, scene());
    }
}
