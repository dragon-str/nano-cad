//! The sorting rotor scene.
//!
//! This module builds the whole rotor from the part generators in
//! [`nanocad_parts`]. Body 0 is the housing, which is fixed. Body 1 is the
//! rotor. The drive shaft is keyed into the rotor bore, so the shaft atoms
//! ride body 1 and no joint is needed between them. Bodies 2 to 13 are the
//! twelve ejection rods. Each rod sits in a radial ejection bore in the rotor,
//! so it is captive in the rotor and its prismatic joint names the rotor as
//! `body_a`. Body 14 is the cam hub, which is fixed and shares the rotor axis.
//!
//! The cam hub lies below the rotor. Its outer surface dwells at the base
//! radius and rises by a short raised-cosine ramp at the outlet azimuth, so a
//! rod is pushed outward only as its pocket meets the housing outlet. The hub
//! is one-sided: the push comes from the hub, and the return comes from a leaf
//! spring on the rotor that holds the follower pin against the hub. The hub
//! carries a flange below the hub that laps under the housing ring, so the
//! housing holds the hub and the hub does not float.
//!
//! The scene uses the schema of the gearbox. The design block is zero, and the
//! viewer reads a zero module as "not a gearbox". It then turns the rotor about
//! `z` from the joint instead of from the gear kinematics.
//!
//! Scope limits. The scene holds no binding site, guest, solvent or ion. The
//! hub surface is a radial key, not a tuned motion law. The follower pin
//! slides on the hub; a roller follower would roll. The scene models no contact
//! force between the pin and the hub, no spring force, no torque on the drive
//! shaft and no drive machine. The ejection work is measured in
//! [`nanocad_meter`].

use std::f64::consts::PI;

use nanocad_model::Part;
use nanocad_parts::{
    place, CamHubGenerator, DriveShaftGenerator, EjectionRodGenerator, FollowerPinGenerator,
    LeafSpringGenerator, ParameterSet, PartError, PartGenerator, RotorHousingGenerator,
    SortingRotorGenerator,
};

use crate::scene::{
    AtomisticLayer, CoarseBody, CoarseCylinder, CoarseLayer, DeviceLayer, Scene, SceneAtom,
    SceneBody, SceneError, SceneJoint, SCENE_SCHEMA, SCENE_VERSION,
};

/// The housing body index.
pub const HOUSING_BODY: usize = 0;
/// The rotor body index.
pub const ROTOR_BODY: usize = 1;
/// The body index of the first ejection rod.
pub const ROD_BODY_FIRST: usize = 2;
/// The number of ejection rods, which equals the number of pockets.
pub const ROD_COUNT: usize = 12;
/// The body index of the cam hub.
pub const CAM_BODY: usize = ROD_BODY_FIRST + ROD_COUNT;
/// The role of the housing body.
pub const HOUSING_ROLE: &str = "housing";
/// The role of the rotor body.
pub const ROTOR_ROLE: &str = "rotor";
/// The role of an ejection rod body.
pub const ROD_ROLE: &str = "ejection_rod";
/// The role of the cam hub body.
pub const CAM_ROLE: &str = "cam_hub";
/// The part id of the drive shaft.
pub const SHAFT_PART_ID: &str = "drive_shaft";
/// The part id of a follower pin.
pub const PIN_PART_ID: &str = "follower_pin";
/// The part id of a leaf spring.
pub const SPRING_PART_ID: &str = "leaf_spring";

/// The frame of the rotor scene.
pub const ROTOR_FRAME: &str = "rotor";

/// The thickness of the housing along z, in metres.
///
/// The chamber is taller than the rotor, so the cam hub fits below the rotor
/// with a running clearance on each face.
pub const HOUSING_THICKNESS_M: f64 = 4.6e-9;

/// The running clearance between two facing surfaces, in metres.
///
/// A C-H bond reaches about 0.11 nm past the nominal carbon surface, so a
/// clearance below about 0.38 nm lets the hydrogen caps overlap. This value
/// keeps every interface clear of the clash test.
pub const INTERFACE_CLEARANCE_M: f64 = 4.0e-10;

/// The running clearance between the ejection rod and the rotor eject bore, in
/// metres.
///
/// The bore must stay below the pocket radius, so it cannot use the full
/// interface clearance. This value still leaves an atom gap above the clash
/// threshold.
pub const EJECTION_BORE_CLEARANCE_M: f64 = 2.9e-10;

/// The free length of one leaf spring, in metres.
pub const SPRING_LENGTH_M: f64 = 6.0e-10;

/// The width of one leaf spring, in metres.
pub const SPRING_WIDTH_M: f64 = 4.0e-10;

/// The thickness of one leaf spring, in metres.
pub const SPRING_THICKNESS_M: f64 = 2.0e-10;

/// Returns the scene atoms of a part, in its own frame.
fn scene_atoms(part: &Part, body: usize, name: &str, part_id: &str) -> Vec<SceneAtom> {
    scene_atoms_mapped(part, body, name, part_id, |point_m| point_m)
}

/// Returns the scene atoms of a part after a point map.
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
            let atomic_number = element.atomic_number();
            let symbol = element.symbol();
            Some(SceneAtom {
                element: match symbol {
                    Some(text) if !text.is_empty() => text.to_string(),
                    _ => format!("Z{atomic_number}"),
                },
                atomic_number,
                position_m: map(position_m),
                body,
                name: name.to_string(),
                part_id: part_id.to_string(),
            })
        })
        .collect()
}

/// Returns the maximum radius and the maximum half-height of a part, in metres.
fn extent_m(part: &Part) -> Option<(f64, f64)> {
    let mut radius_m = 0.0_f64;
    let mut half_height_m = 0.0_f64;
    for index in 0..part.atom_count() {
        let position_m = part.topology.position_m(index)?;
        radius_m = radius_m.max(position_m[0].hypot(position_m[1]));
        half_height_m = half_height_m.max(position_m[2].abs());
    }
    Some((radius_m, half_height_m))
}

/// Returns a default parameter of a generator, or zero when it is absent.
fn default_m(generator: &impl PartGenerator, name: &str) -> f64 {
    generator
        .resolve(&ParameterSet::new())
        .ok()
        .and_then(|resolved| resolved.get(name))
        .unwrap_or(0.0)
}

/// Returns a scene body with no initial motion.
fn scene_body(index: usize, role: &str, part_id: &str, fixed: bool) -> SceneBody {
    SceneBody {
        index,
        name: part_id.to_string(),
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
/// The scene holds the housing, the rotor with its keyed drive shaft, the
/// twelve ejection rods with their follower pins and leaf springs, and the
/// fixed cam hub with its flange.
pub fn build_rotor_scene() -> Result<Scene, SceneError> {
    let rod_radius_m = default_m(&EjectionRodGenerator, "shaft_radius_m");
    let rotor = SortingRotorGenerator.generate(&ParameterSet::new().with(
        "ejection_bore_radius_m",
        rod_radius_m + EJECTION_BORE_CLEARANCE_M,
    ))?;
    let housing = RotorHousingGenerator
        .generate(&ParameterSet::new().with("thickness_m", HOUSING_THICKNESS_M))?;
    let shaft = DriveShaftGenerator.generate_with_defaults()?;
    let pin = FollowerPinGenerator.generate_with_defaults()?;
    let spring = LeafSpringGenerator.generate(
        &ParameterSet::new()
            .with("leaf_length_m", SPRING_LENGTH_M)
            .with("leaf_width_m", SPRING_WIDTH_M)
            .with("leaf_thickness_m", SPRING_THICKNESS_M),
    )?;
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
    let tip_length_m = default_m(&EjectionRodGenerator, "tip_length_m");
    let disc_radius_m = default_m(&SortingRotorGenerator, "disc_radius_m");
    let rotor_half_thickness_m = default_m(&SortingRotorGenerator, "thickness_m") / 2.0;
    let housing_half_thickness_m = HOUSING_THICKNESS_M / 2.0;

    let hub_base_m = CamHubGenerator.base_radius_m();
    let hub_rise_m = CamHubGenerator.rise_m();
    let pin_head_radius_m = FollowerPinGenerator.head_radius_m();
    let pin_length_m = FollowerPinGenerator.length_m();
    let cam_retract_m = hub_base_m + pin_head_radius_m + INTERFACE_CLEARANCE_M;

    let hub_top_m = -rotor_half_thickness_m - INTERFACE_CLEARANCE_M;
    let hub_bottom_m = -housing_half_thickness_m + INTERFACE_CLEARANCE_M;
    let hub_thickness_m = hub_top_m - hub_bottom_m;
    if hub_thickness_m <= 0.0 {
        return Err(SceneError::Part(format!(
            "the housing half thickness {housing_half_thickness_m} m leaves no room \
             for the cam hub below the rotor"
        )));
    }
    let cam =
        CamHubGenerator.generate(&ParameterSet::new().with("thickness_m", hub_thickness_m))?;
    let hub_offset_m = hub_top_m - 0.5 * hub_thickness_m;

    if cam_retract_m <= rotor_bore_m + rod_radius_m {
        return Err(SceneError::Part(format!(
            "the cam hub retracts to {cam_retract_m} m, which reaches the rotor bore \
             {rotor_bore_m} m plus the rod radius {rod_radius_m} m"
        )));
    }
    let rod_length_m = pocket_inner_m - cam_retract_m - INTERFACE_CLEARANCE_M;
    let rod_tip_m = pocket_inner_m - INTERFACE_CLEARANCE_M + hub_rise_m;
    if rod_tip_m > disc_radius_m + 1.0e-15 {
        return Err(SceneError::Part(format!(
            "the rod reaches {rod_tip_m} m at full stroke, past the rotor rim \
             {disc_radius_m} m, so it would cross the housing"
        )));
    }
    let shaft_length_m = rod_length_m - tip_length_m;
    if shaft_length_m <= 0.0 {
        return Err(SceneError::Part(
            "the cam hub reaches the pocket inner wall, so no rod fits".to_string(),
        ));
    }
    let rod = EjectionRodGenerator
        .generate(&ParameterSet::new().with("shaft_length_m", shaft_length_m))?;
    let half_shaft_m = 0.5 * shaft_length_m;

    let shaft_radius_m = default_m(&DriveShaftGenerator, "shaft_radius_m");
    if shaft_radius_m >= rotor_bore_m || shaft_radius_m >= CamHubGenerator.bore_radius_m() {
        return Err(SceneError::Part(format!(
            "the drive shaft radius {shaft_radius_m} m does not clear the rotor bore \
             {rotor_bore_m} m and the cam hub bore"
        )));
    }

    let pin_top_m = -rotor_half_thickness_m - INTERFACE_CLEARANCE_M;
    let pin_bottom_m = pin_top_m - pin_length_m;
    if pin_top_m < hub_bottom_m || pin_bottom_m > hub_top_m {
        return Err(SceneError::Part(format!(
            "the follower pin spans z {pin_bottom_m}..{pin_top_m} m and does not overlap \
             the cam hub {hub_bottom_m}..{hub_top_m} m, so the pin cannot follow the lobe"
        )));
    }

    let housing_atoms = scene_atoms(&housing, HOUSING_BODY, HOUSING_ROLE, "rotor_housing");
    let mut rotor_atoms = scene_atoms(&rotor, ROTOR_BODY, ROTOR_ROLE, "sorting_rotor");
    rotor_atoms.extend(scene_atoms(&shaft, ROTOR_BODY, ROTOR_ROLE, SHAFT_PART_ID));
    let cam_atoms = scene_atoms_mapped(&cam, CAM_BODY, CAM_ROLE, "cam_hub", |point_m| {
        [point_m[0], point_m[1], point_m[2] + hub_offset_m]
    });

    let mut atoms = housing_atoms;
    atoms.extend(rotor_atoms.iter().cloned());

    let mut rod_bodies = Vec::with_capacity(ROD_COUNT);
    let mut rod_joints = Vec::with_capacity(ROD_COUNT);
    let mut rod_coarse = Vec::with_capacity(ROD_COUNT);
    for index in 0..ROD_COUNT {
        let angle_rad = 2.0 * PI * index as f64 / pocket_count.max(1) as f64;
        let (sin_rad, cos_rad) = angle_rad.sin_cos();
        let radial_m = [cos_rad, sin_rad, 0.0];
        let pocket_center_m =
            SortingRotorGenerator::pocket_center_m(index, pocket_count, pocket_orbit_m);
        let body = ROD_BODY_FIRST + index;
        let pin_center_m = CamHubGenerator.pin_center_radius_m(angle_rad, pin_head_radius_m)
            + INTERFACE_CLEARANCE_M;
        let shaft_center_m = pin_center_m + half_shaft_m;
        let spring_center_m =
            pin_center_m + pin_head_radius_m + INTERFACE_CLEARANCE_M + 0.5 * SPRING_WIDTH_M;

        let rod_atoms = scene_atoms_mapped(&rod, body, ROD_ROLE, "ejection_rod", |point_m| {
            [
                point_m[2] * cos_rad - point_m[1] * sin_rad + shaft_center_m * cos_rad,
                point_m[2] * sin_rad + point_m[1] * cos_rad + shaft_center_m * sin_rad,
                -point_m[0],
            ]
        });
        let pin_atoms = scene_atoms_mapped(&pin, body, ROD_ROLE, PIN_PART_ID, |point_m| {
            [
                pin_center_m * cos_rad + point_m[0],
                pin_center_m * sin_rad + point_m[1],
                pin_top_m + point_m[2],
            ]
        });
        let spring_atoms = scene_atoms_mapped(&spring, body, ROD_ROLE, SPRING_PART_ID, |point_m| {
            [
                spring_center_m * cos_rad - point_m[1] * cos_rad - point_m[0] * sin_rad,
                spring_center_m * sin_rad - point_m[1] * sin_rad + point_m[0] * cos_rad,
                -rotor_half_thickness_m - 0.5 * SPRING_LENGTH_M + point_m[2],
            ]
        });

        let mut rod_atom_set = rod_atoms;
        rod_atom_set.extend(pin_atoms);
        rod_atom_set.extend(spring_atoms);
        rod_coarse.push(CoarseBody {
            index: body,
            name: format!("ejection_rod_{index}"),
            role: ROD_ROLE.to_string(),
            pitch_circle: None,
            bounding_cylinder: local_cylinder(&rod_atom_set, radial_m),
        });
        atoms.extend(rod_atom_set.iter().cloned());
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
    bodies.push(scene_body(CAM_BODY, CAM_ROLE, "cam_hub", true));

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
            name: "rotor_housing".to_string(),
            role: HOUSING_ROLE.to_string(),
            pitch_circle: None,
            bounding_cylinder: bounding_cylinder(&housing, axis, housing_half_thickness_m),
        },
        CoarseBody {
            index: ROTOR_BODY,
            name: "sorting_rotor".to_string(),
            role: ROTOR_ROLE.to_string(),
            pitch_circle: None,
            bounding_cylinder: local_cylinder(&rotor_atoms, axis),
        },
    ];
    coarse.extend(rod_coarse);
    coarse.push(CoarseBody {
        index: CAM_BODY,
        name: "cam_hub".to_string(),
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

/// Returns the design block of a non-gearbox scene.
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

impl From<PartError> for SceneError {
    fn from(error: PartError) -> Self {
        SceneError::Part(error.to_string())
    }
}

/// Returns a coarse cylinder that wraps a part along an axis.
fn bounding_cylinder(part: &Part, axis: [f64; 3], half_length_m: f64) -> Option<CoarseCylinder> {
    let (radius_m, atom_half_height_m) = extent_m(part)?;
    Some(CoarseCylinder {
        center_m: [0.0, 0.0, 0.0],
        axis,
        radius_m,
        half_length_m: half_length_m.max(atom_half_height_m),
    })
}

/// Returns a coarse cylinder that wraps a set of atoms along an axis.
fn local_cylinder(atoms: &[SceneAtom], axis: [f64; 3]) -> Option<CoarseCylinder> {
    if atoms.is_empty() {
        return None;
    }
    let mut radius_m = 0.0_f64;
    let mut min_along_m = f64::INFINITY;
    let mut max_along_m = f64::NEG_INFINITY;
    for atom in atoms {
        let along_m = atom.position_m[0] * axis[0]
            + atom.position_m[1] * axis[1]
            + atom.position_m[2] * axis[2];
        let perpendicular_m = [
            atom.position_m[0] - along_m * axis[0],
            atom.position_m[1] - along_m * axis[1],
            atom.position_m[2] - along_m * axis[2],
        ];
        radius_m = radius_m.max(
            perpendicular_m[0]
                .hypot(perpendicular_m[1])
                .hypot(perpendicular_m[2]),
        );
        min_along_m = min_along_m.min(along_m);
        max_along_m = max_along_m.max(along_m);
    }
    let center_along_m = 0.5 * (min_along_m + max_along_m);
    Some(CoarseCylinder {
        center_m: [
            center_along_m * axis[0],
            center_along_m * axis[1],
            center_along_m * axis[2],
        ],
        axis,
        radius_m,
        half_length_m: 0.5 * (max_along_m - min_along_m),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::OnceLock;

    use crate::scene::{scene_from_json, scene_to_json};

    fn scene() -> &'static Scene {
        static SCENE: OnceLock<Scene> = OnceLock::new();
        SCENE.get_or_init(|| build_rotor_scene().expect("the rotor scene builds"))
    }

    fn body_atoms(scene: &Scene, body: usize) -> usize {
        scene
            .atomistic
            .atoms
            .iter()
            .filter(|atom| atom.body == body)
            .count()
    }

    fn part_atoms(scene: &Scene, body: usize, part_id: &str) -> usize {
        scene
            .atomistic
            .atoms
            .iter()
            .filter(|atom| atom.body == body && atom.part_id == part_id)
            .count()
    }

    #[test]
    fn the_rotor_scene_uses_the_shared_schema() {
        let scene = scene();
        assert_eq!(scene.schema, SCENE_SCHEMA);
        assert_eq!(scene.version, SCENE_VERSION);
        assert_eq!(scene.length_unit, "m");
        assert_eq!(scene.angle_unit, "rad");
        assert_eq!(scene.frame, "rotor");
    }

    #[test]
    fn the_design_block_marks_a_non_gearbox() {
        let scene = scene();
        assert_eq!(scene.design.module_m, 0.0);
        assert_eq!(scene.design.planet_count, 0);
    }

    #[test]
    fn the_rotor_scene_has_fifteen_bodies_and_thirteen_joints() {
        let scene = scene();
        assert_eq!(scene.device.body_count, 15);
        assert_eq!(scene.device.bodies.len(), 15);
        assert_eq!(scene.device.joints.len(), 13);
        assert_eq!(scene.coarse.body_count, 15);
        let axis_joint = &scene.device.joints[0];
        assert_eq!(axis_joint.kind, "revolute");
        assert_eq!(axis_joint.body_a, HOUSING_BODY);
        assert_eq!(axis_joint.body_b, ROTOR_BODY);
        assert_eq!(axis_joint.axis_world, [0.0, 0.0, 1.0]);
        for (index, joint) in scene.device.joints[1..].iter().enumerate() {
            assert_eq!(joint.kind, "prismatic");
            assert_eq!(joint.body_a, ROTOR_BODY);
            assert_eq!(joint.body_b, ROD_BODY_FIRST + index);
            assert!(joint.anchor_world_m.is_some());
        }
    }

    #[test]
    fn every_rod_is_free_and_slides_on_the_rotor() {
        let scene = scene();
        for index in 0..ROD_COUNT {
            let body = &scene.device.bodies[ROD_BODY_FIRST + index];
            assert!(!body.fixed);
            assert_eq!(body.role, ROD_ROLE);
        }
    }

    #[test]
    fn the_cam_hub_is_fixed_to_the_housing() {
        let scene = scene();
        let body = &scene.device.bodies[CAM_BODY];
        assert!(body.fixed);
        assert_eq!(body.role, CAM_ROLE);
        assert_eq!(body.part_id.as_deref(), Some("cam_hub"));
        assert!(scene
            .device
            .joints
            .iter()
            .all(|joint| joint.body_a != CAM_BODY && joint.body_b != CAM_BODY));
    }

    #[test]
    fn the_housing_is_fixed_and_the_rotor_is_free() {
        let scene = scene();
        assert!(scene.device.bodies[HOUSING_BODY].fixed);
        assert!(!scene.device.bodies[ROTOR_BODY].fixed);
        assert_eq!(scene.device.bodies[HOUSING_BODY].role, HOUSING_ROLE);
        assert_eq!(scene.device.bodies[ROTOR_BODY].role, ROTOR_ROLE);
    }

    #[test]
    fn the_drive_shaft_turns_with_the_rotor() {
        let scene = scene();
        let shaft_atoms = part_atoms(scene, ROTOR_BODY, SHAFT_PART_ID);
        assert!(shaft_atoms > 0);
        let shaft_joints = scene
            .device
            .joints
            .iter()
            .filter(|joint| joint.body_b == ROTOR_BODY)
            .count();
        assert_eq!(shaft_joints, 1);
    }

    #[test]
    fn every_rod_carries_a_follower_pin() {
        let scene = scene();
        for index in 0..ROD_COUNT {
            let body = ROD_BODY_FIRST + index;
            assert!(
                part_atoms(scene, body, PIN_PART_ID) > 0,
                "rod {index} has no pin"
            );
        }
    }

    #[test]
    fn every_rod_carries_a_leaf_spring() {
        let scene = scene();
        for index in 0..ROD_COUNT {
            let body = ROD_BODY_FIRST + index;
            assert!(
                part_atoms(scene, body, SPRING_PART_ID) > 0,
                "rod {index} has no spring"
            );
        }
    }

    #[test]
    fn every_atom_belongs_to_one_of_the_fifteen_bodies() {
        let scene = scene();
        assert_eq!(scene.atomistic.atom_count, scene.atomistic.atoms.len());
        for atom in &scene.atomistic.atoms {
            assert!(atom.body < 15);
            assert!(atom.atomic_number > 0);
        }
        assert_eq!(scene.atomistic.atom_count, 127_036);
        assert_eq!(body_atoms(scene, HOUSING_BODY), 71_624);
        assert_eq!(body_atoms(scene, ROTOR_BODY), 42_732);
        assert_eq!(body_atoms(scene, CAM_BODY), 1_184);
        for index in 0..ROD_COUNT {
            assert_eq!(body_atoms(scene, ROD_BODY_FIRST + index), 958);
            assert_eq!(
                part_atoms(scene, ROD_BODY_FIRST + index, SPRING_PART_ID),
                36
            );
        }
    }

    #[test]
    fn the_coarse_bodies_wrap_the_atoms_of_their_body() {
        let scene = scene();
        for body in &scene.coarse.bodies {
            let cylinder = body.bounding_cylinder.as_ref().expect("a coarse cylinder");
            let axis = cylinder.axis;
            let half_length_m = cylinder.half_length_m + 1.0e-15;
            for atom in scene
                .atomistic
                .atoms
                .iter()
                .filter(|a| a.body == body.index)
            {
                let along_m = atom.position_m[0] * axis[0]
                    + atom.position_m[1] * axis[1]
                    + atom.position_m[2] * axis[2]
                    - (cylinder.center_m[0] * axis[0]
                        + cylinder.center_m[1] * axis[1]
                        + cylinder.center_m[2] * axis[2]);
                let perpendicular_m = [
                    atom.position_m[0] - cylinder.center_m[0] - along_m * axis[0],
                    atom.position_m[1] - cylinder.center_m[1] - along_m * axis[1],
                    atom.position_m[2] - cylinder.center_m[2] - along_m * axis[2],
                ];
                let radius_m = perpendicular_m[0]
                    .hypot(perpendicular_m[1])
                    .hypot(perpendicular_m[2]);
                assert!(
                    radius_m <= cylinder.radius_m + 1.0e-15,
                    "an atom of body {} sits outside its cylinder",
                    body.index
                );
                assert!(along_m.abs() <= half_length_m);
            }
        }
    }

    #[test]
    fn the_rotor_scene_round_trips_through_json() {
        let scene = scene();
        let json = scene_to_json(scene).expect("the scene writes");
        let restored = scene_from_json(&json).expect("the scene reads");
        assert_eq!(restored.atomistic.atom_count, scene.atomistic.atom_count);
        assert_eq!(restored.device.body_count, scene.device.body_count);
        assert_eq!(restored.device.joints.len(), scene.device.joints.len());
        assert_eq!(restored.coarse.body_count, scene.coarse.body_count);
    }

    #[test]
    fn the_rod_and_the_cam_hub_relax_without_a_clash() {
        let scene = scene();
        let atoms: Vec<nanocad_meter::RelaxAtom> = scene
            .atomistic
            .atoms
            .iter()
            .filter(|atom| atom.body == ROD_BODY_FIRST || atom.body == CAM_BODY)
            .map(|atom| nanocad_meter::RelaxAtom {
                position_m: atom.position_m,
                element: nanocad_model::Element::from_atomic_number(atom.atomic_number)
                    .unwrap_or(nanocad_model::Element::CARBON),
                body: atom.body,
            })
            .collect();
        let report =
            nanocad_meter::relax_subassembly(&atoms, &nanocad_meter::RelaxTarget::default())
                .expect("the sub-assembly relaxes");
        assert!(report.converged);
        assert_eq!(report.initial_clash_count, 0);
        assert_eq!(report.clash_count, 0);
        assert!(report.max_bond_strain < 1.0e-3);
        assert!(report.max_displacement_m < 1.0e-11);
    }
}
