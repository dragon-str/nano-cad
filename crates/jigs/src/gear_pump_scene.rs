//! The selective gear pump scene (M27-02, ADR-0086).
//!
//! Body 0 is the figure-eight housing, which is fixed. It carries the two
//! posts, so it is one lattice. Body 1 is gear A, the driver, at `-x`. Body 2
//! is gear B, the follower, at `+x`. Each gear has a revolute joint to the
//! housing at its post. One gear coupling ties the two joints: the gears are
//! equal, so they turn at the same rate in opposite senses.
//!
//! The design block is zero, so the viewer does not read the scene as the
//! planetary gearbox. It turns each gear about its own joint anchor at the
//! rate that the gear coupling sets.
//!
//! Scope limits. The scene holds no binding site, guest or solvent. It models
//! no drive torque and no contact force.

use nanocad_model::Part;
use nanocad_parts::{GearPumpGenerator, GearPumpSet, ParameterSet};
use std::ops::Range;

use crate::rotor_scene::zero_design;
use crate::scene::{
    AtomisticLayer, CoarseBody, CoarseCircle, CoarseCylinder, CoarseLayer, DeviceLayer, Scene,
    SceneAtom, SceneBody, SceneError, SceneGearCoupling, SceneJoint, SCENE_SCHEMA, SCENE_VERSION,
};

/// The id of the machine scene of the gear pump.
pub const GEAR_PUMP_SCENE_ID: &str = "gear_pump";

/// The frame of the gear-pump scene.
pub const GEAR_PUMP_FRAME: &str = "gear_pump";

/// The body index of the pump housing.
pub const PUMP_HOUSING_BODY: usize = 0;

/// The body index of gear A, the driver.
pub const PUMP_GEAR_A_BODY: usize = 1;

/// The body index of gear B, the follower.
pub const PUMP_GEAR_B_BODY: usize = 2;

/// The role of the pump housing.
pub const PUMP_HOUSING_ROLE: &str = "pump_housing";

/// The role of a pump gear.
pub const PUMP_GEAR_ROLE: &str = "pump_gear";

const AXIS: [f64; 3] = [0.0, 0.0, 1.0];

/// Returns the scene atoms of one atom range of a part.
fn range_atoms(part: &Part, range: &Range<usize>, body: usize, part_id: &str) -> Vec<SceneAtom> {
    range
        .clone()
        .filter_map(|index| {
            let element = part.topology.element(index)?;
            let position_m = part.topology.position_m(index)?;
            let atomic_number = element.atomic_number();
            Some(SceneAtom {
                element: element
                    .symbol()
                    .filter(|symbol| !symbol.is_empty())
                    .map(str::to_string)
                    .unwrap_or_else(|| format!("Z{atomic_number}")),
                atomic_number,
                position_m,
                body,
                name: part_id.to_string(),
                part_id: part_id.to_string(),
            })
        })
        .collect()
}

fn body(index: usize, role: &str, name: &str, center_m: [f64; 2], fixed: bool) -> SceneBody {
    SceneBody {
        index,
        name: name.to_string(),
        role: role.to_string(),
        part_id: Some(name.to_string()),
        position_m: [center_m[0], center_m[1], 0.0],
        orientation_wxyz: [1.0, 0.0, 0.0, 0.0],
        mass_kg: 0.0,
        inertia_diagonal_kg_m2: [0.0, 0.0, 0.0],
        fixed,
        linear_velocity_m_per_s: [0.0, 0.0, 0.0],
        angular_velocity_rad_per_s: [0.0, 0.0, 0.0],
    }
}

fn revolute(name: &str, gear_body: usize, center_m: [f64; 2]) -> SceneJoint {
    SceneJoint {
        name: name.to_string(),
        kind: "revolute".to_string(),
        body_a: PUMP_HOUSING_BODY,
        body_b: gear_body,
        anchor_world_m: Some([center_m[0], center_m[1], 0.0]),
        axis_world: AXIS,
        axis_a_body: AXIS,
        axis_b_body: Some(AXIS),
    }
}

/// Builds the gear-pump scene from the default pump parameters.
pub fn build_gear_pump_scene() -> Result<Scene, SceneError> {
    let set = GearPumpGenerator.build(&ParameterSet::new())?;
    Ok(gear_pump_scene(&set))
}

/// Builds the gear-pump scene from a generated pump.
pub fn gear_pump_scene(set: &GearPumpSet) -> Scene {
    let design = &set.design;
    let a_center_m = design.gear_a_center_m();
    let b_center_m = design.gear_b_center_m();

    let mut atoms = range_atoms(
        &set.part,
        &set.housing_atoms,
        PUMP_HOUSING_BODY,
        "pump_housing",
    );
    atoms.extend(range_atoms(
        &set.part,
        &set.gear_a_atoms,
        PUMP_GEAR_A_BODY,
        "pump_gear_a",
    ));
    atoms.extend(range_atoms(
        &set.part,
        &set.gear_b_atoms,
        PUMP_GEAR_B_BODY,
        "pump_gear_b",
    ));

    let bodies = vec![
        body(
            PUMP_HOUSING_BODY,
            PUMP_HOUSING_ROLE,
            "pump_housing",
            [0.0, 0.0],
            true,
        ),
        body(
            PUMP_GEAR_A_BODY,
            PUMP_GEAR_ROLE,
            "pump_gear_a",
            a_center_m,
            false,
        ),
        body(
            PUMP_GEAR_B_BODY,
            PUMP_GEAR_ROLE,
            "pump_gear_b",
            b_center_m,
            false,
        ),
    ];
    let joints = vec![
        revolute("gear_a_axis", PUMP_GEAR_A_BODY, a_center_m),
        revolute("gear_b_axis", PUMP_GEAR_B_BODY, b_center_m),
    ];
    let gear_couplings = vec![SceneGearCoupling {
        joint_a: 0,
        joint_b: 1,
        radius_a_m: design.pitch_radius_m(),
        radius_b_m: design.pitch_radius_m(),
    }];

    let gear_half_m = design.gear_half_thickness_m();
    let gear_coarse = |index: usize, name: &str, center_m: [f64; 2]| CoarseBody {
        index,
        name: name.to_string(),
        role: PUMP_GEAR_ROLE.to_string(),
        pitch_circle: Some(CoarseCircle {
            center_m: [center_m[0], center_m[1], 0.0],
            axis: AXIS,
            radius_m: design.pitch_radius_m(),
        }),
        bounding_cylinder: Some(CoarseCylinder {
            center_m: [center_m[0], center_m[1], 0.0],
            axis: AXIS,
            radius_m: design.outer_radius_m(),
            half_length_m: gear_half_m,
        }),
    };
    let coarse = vec![
        CoarseBody {
            index: PUMP_HOUSING_BODY,
            name: "pump_housing".to_string(),
            role: PUMP_HOUSING_ROLE.to_string(),
            pitch_circle: None,
            bounding_cylinder: Some(CoarseCylinder {
                center_m: [0.0, 0.0, 0.0],
                axis: AXIS,
                radius_m: b_center_m[0] + design.housing_radius_m(),
                half_length_m: design.housing_half_height_m(),
            }),
        },
        gear_coarse(PUMP_GEAR_A_BODY, "pump_gear_a", a_center_m),
        gear_coarse(PUMP_GEAR_B_BODY, "pump_gear_b", b_center_m),
    ];

    Scene {
        schema: SCENE_SCHEMA.to_string(),
        version: SCENE_VERSION,
        length_unit: "m".to_string(),
        angle_unit: "rad".to_string(),
        frame: GEAR_PUMP_FRAME.to_string(),
        design: zero_design(),
        atomistic: AtomisticLayer {
            atom_count: atoms.len(),
            atoms,
        },
        device: DeviceLayer {
            body_count: bodies.len(),
            bodies,
            joints,
            gear_couplings,
            gear_constraints: Vec::new(),
        },
        coarse: CoarseLayer {
            body_count: coarse.len(),
            bodies: coarse,
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sense::sense_report;

    #[test]
    fn the_pump_has_three_bodies_two_joints_and_one_coupling() {
        let scene = build_gear_pump_scene().expect("the pump scene builds");
        assert_eq!(scene.device.body_count, 3);
        assert_eq!(scene.device.joints.len(), 2);
        assert_eq!(scene.device.gear_couplings.len(), 1);
        assert!(scene.device.bodies[PUMP_HOUSING_BODY].fixed);
        assert!(!scene.device.bodies[PUMP_GEAR_A_BODY].fixed);
        assert_eq!(scene.design.module_m, 0.0);
        for body in 0..3 {
            assert!(scene.atomistic.atoms.iter().any(|atom| atom.body == body));
        }
    }

    #[test]
    fn the_pump_bodies_are_whole_held_and_apart() {
        let scene = build_gear_pump_scene().expect("the pump scene builds");
        let report = sense_report(&scene);
        assert!(report.floating_bodies.is_empty(), "{report:?}");
        assert!(report.disjoint_bodies.is_empty(), "{report:?}");
        assert_eq!(report.cross_body_bonds, 0);
        // Only same-body pairs may overlap. Those are the known hydrogen cap
        // clash at a shared missing lattice site (ADR-0075), not a gear fault.
        assert!(
            report
                .overlap_bodies
                .iter()
                .all(|pair| pair.body_a == pair.body_b),
            "{report:?}"
        );
    }

    #[test]
    fn the_pump_scene_survives_json() {
        let scene = build_gear_pump_scene().expect("the pump scene builds");
        let json = crate::scene::scene_to_json(&scene).expect("json");
        let back = crate::scene::scene_from_json(&json).expect("parse");
        assert_eq!(back.device.gear_couplings, scene.device.gear_couplings);
    }
}
