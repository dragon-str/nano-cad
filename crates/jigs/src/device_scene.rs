//! The assembled respirocyte device (ADR-0115).
//!
//! Body 0 is the module housing (fixed). Body 1 is the oxygen-loading rotor,
//! with a revolute joint about z. Bodies 2..2+n are the heme sites and bodies
//! 2+n..2+2n are the ligand arms, as in the module scene. Body 2+2n is the
//! gas-tank shell, which encloses the module. Bodies 2+2n+1..2+2n+3 are the
//! pump cylinder and piston, parked above the module. The pump seal ring is
//! a sub-bond-thickness lattice cut, so it is not one connected body; it is
//! left out of the scene (see ADR-0115).
//!
//! Every subsystem relaxes on its own, because the parts do not touch. The
//! layout is schematic: the tank and the pump are drawn at the scale of the
//! module, not at the scale of a real respirocyte hull.

use nanocad_meter::ValenceRelaxOptions;
use nanocad_model::Part;
use nanocad_parts::{
    respirocyte_pump, respirocyte_tank, ParameterSet, RespirocytePump, RespirocyteTank,
};

use crate::module_scene::{module_parts, relax_module_parts, relax_part, ModuleParts};
use crate::rotor_scene::zero_design;
use crate::scene::{
    AtomisticLayer, CoarseBody, CoarseCylinder, CoarseLayer, DeviceLayer, Scene, SceneAtom,
    SceneBody, SceneError, SceneJoint, SCENE_SCHEMA, SCENE_VERSION,
};

/// The scene id of the assembled device.
pub const DEVICE_SCENE_ID: &str = "respirocyte_device";

/// The radial clearance between the module and the tank inner wall, in metres.
pub const DEVICE_TANK_CLEARANCE_M: f64 = 0.6e-9;

/// The tank shell wall thickness in the device, in metres.
pub const DEVICE_TANK_WALL_M: f64 = 0.6e-9;

/// The world position of the pump cylinder centre, in metres.
pub const DEVICE_PUMP_OFFSET_M: [f64; 3] = [0.0, 0.0, 3.6e-9];

const Z: [f64; 3] = [0.0, 0.0, 1.0];

/// The generated bodies of the device.
pub struct DeviceParts {
    /// The module (rotor, sites, arms and housing).
    pub module: ModuleParts,
    /// The gas tank.
    pub tank: RespirocyteTank,
    /// The pump cylinder, piston and seal.
    pub pump: RespirocytePump,
}

/// Generates the device parts, unrelaxed, with the pump placed above the module.
pub fn device_parts() -> Result<DeviceParts, SceneError> {
    let module = module_parts()?;
    let module_radius_m = module_radius_m(&module);
    let tank = respirocyte_tank(
        &ParameterSet::new()
            .with(
                "tank_inner_radius_m",
                module_radius_m + DEVICE_TANK_CLEARANCE_M,
            )
            .with("tank_wall_m", DEVICE_TANK_WALL_M),
    )?;
    let mut pump = respirocyte_pump(&ParameterSet::new())?;
    translate_part(&mut pump.cylinder, DEVICE_PUMP_OFFSET_M)?;
    translate_part(&mut pump.piston, DEVICE_PUMP_OFFSET_M)?;
    translate_part(&mut pump.seal, DEVICE_PUMP_OFFSET_M)?;
    Ok(DeviceParts { module, tank, pump })
}

fn module_radius_m(module: &ModuleParts) -> f64 {
    let mut radius_m: f64 = 0.0;
    let parts = std::iter::once(&module.housing)
        .chain(std::iter::once(&module.set.rotor))
        .chain(module.set.sites.iter())
        .chain(module.set.arms.iter());
    for part in parts {
        for i in 0..part.atom_count() {
            if let Some(p) = part.topology.position_m(i) {
                radius_m = radius_m.max((p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt());
            }
        }
    }
    radius_m
}

fn translate_part(part: &mut Part, offset_m: [f64; 3]) -> Result<(), SceneError> {
    for i in 0..part.atom_count() {
        if let Some(p) = part.topology.position_m(i) {
            part.topology
                .set_position_m(
                    i,
                    [p[0] + offset_m[0], p[1] + offset_m[1], p[2] + offset_m[2]],
                )
                .map_err(|e| SceneError::Part(e.to_string()))?;
        }
    }
    Ok(())
}

fn atoms_of(part: &Part, body: usize, name: &str) -> Vec<SceneAtom> {
    (0..part.atom_count())
        .filter_map(|i| {
            let e = part.topology.element(i)?;
            let n = e.atomic_number();
            Some(SceneAtom {
                element: e
                    .symbol()
                    .map(str::to_string)
                    .unwrap_or_else(|| format!("Z{n}")),
                atomic_number: n,
                position_m: part.topology.position_m(i)?,
                body,
                name: name.to_string(),
                part_id: name.to_string(),
            })
        })
        .collect()
}

fn body(index: usize, role: &str, name: &str, fixed: bool) -> SceneBody {
    SceneBody {
        index,
        name: name.to_string(),
        role: role.to_string(),
        part_id: Some(name.to_string()),
        position_m: [0.0; 3],
        orientation_wxyz: [1.0, 0.0, 0.0, 0.0],
        mass_kg: 0.0,
        inertia_diagonal_kg_m2: [0.0; 3],
        fixed,
        linear_velocity_m_per_s: [0.0; 3],
        angular_velocity_rad_per_s: [0.0; 3],
    }
}

fn joint(
    name: &str,
    kind: &str,
    a: usize,
    b: usize,
    anchor: [f64; 3],
    axis: [f64; 3],
) -> SceneJoint {
    SceneJoint {
        name: name.to_string(),
        kind: kind.to_string(),
        body_a: a,
        body_b: b,
        anchor_world_m: Some(anchor),
        axis_world: axis,
        axis_a_body: axis,
        axis_b_body: Some(axis),
    }
}

/// Builds the assembled device scene from its relaxed parts.
pub fn device_scene(parts: &DeviceParts) -> Scene {
    let set = &parts.module.set;
    let n = set.sites.len();
    let tank_body = 2 + 2 * n;
    let pump_cylinder_body = tank_body + 1;
    let pump_piston_body = tank_body + 2;

    let mut atoms = atoms_of(&parts.module.housing, 0, "module_housing");
    atoms.extend(atoms_of(&set.rotor, 1, "heme_rotor"));
    let mut bodies = vec![
        body(0, "housing", "module_housing", true),
        body(1, "rotor", "heme_rotor", false),
    ];
    let mut joints = vec![joint("rotor_axis", "revolute", 0, 1, [0.0; 3], Z)];
    for k in 0..n {
        let name = format!("heme_site_{k}");
        atoms.extend(atoms_of(&set.sites[k], 2 + k, &name));
        bodies.push(body(2 + k, "heme_site", &name, false));
        let (s, c) = set.azimuth_rad[k].sin_cos();
        let fe = [set.iron_radius_m * c, set.iron_radius_m * s, 0.0];
        joints.push(joint(
            &format!("site_{k}_attachment"),
            "fixed",
            1,
            2 + k,
            fe,
            Z,
        ));
    }
    for k in 0..n {
        let name = format!("ligand_arm_{k}");
        let index = 2 + n + k;
        atoms.extend(atoms_of(&set.arms[k], index, &name));
        bodies.push(body(index, "ligand_arm", &name, false));
        let (s, c) = set.azimuth_rad[k].sin_cos();
        let fe = [set.iron_radius_m * c, set.iron_radius_m * s, 0.0];
        joints.push(joint(
            &format!("arm_{k}_slide"),
            "prismatic",
            1,
            index,
            fe,
            [c, s, 0.0],
        ));
    }
    atoms.extend(atoms_of(&parts.tank.shell, tank_body, "gas_tank"));
    bodies.push(body(tank_body, "tank", "gas_tank", true));
    joints.push(joint("tank_mount", "fixed", 0, tank_body, [0.0; 3], Z));

    atoms.extend(atoms_of(
        &parts.pump.cylinder,
        pump_cylinder_body,
        "pump_cylinder",
    ));
    bodies.push(body(
        pump_cylinder_body,
        "pump_cylinder",
        "pump_cylinder",
        true,
    ));
    joints.push(joint(
        "pump_mount",
        "fixed",
        0,
        pump_cylinder_body,
        DEVICE_PUMP_OFFSET_M,
        Z,
    ));
    atoms.extend(atoms_of(
        &parts.pump.piston,
        pump_piston_body,
        "pump_piston",
    ));
    bodies.push(body(pump_piston_body, "pump_piston", "pump_piston", false));
    joints.push(joint(
        "piston_stroke",
        "prismatic",
        pump_cylinder_body,
        pump_piston_body,
        DEVICE_PUMP_OFFSET_M,
        Z,
    ));
    let coarse = vec![
        CoarseBody {
            index: 1,
            name: "heme_rotor".to_string(),
            role: "rotor".to_string(),
            pitch_circle: None,
            bounding_cylinder: Some(CoarseCylinder {
                center_m: [0.0; 3],
                axis: Z,
                radius_m: set.disc_radius_m,
                half_length_m: set.half_thickness_m,
            }),
        },
        CoarseBody {
            index: tank_body,
            name: "gas_tank".to_string(),
            role: "tank".to_string(),
            pitch_circle: None,
            bounding_cylinder: Some(CoarseCylinder {
                center_m: [0.0; 3],
                axis: Z,
                radius_m: parts.tank.outer_radius_m,
                half_length_m: parts.tank.outer_radius_m,
            }),
        },
        CoarseBody {
            index: pump_cylinder_body,
            name: "pump_cylinder".to_string(),
            role: "pump_cylinder".to_string(),
            pitch_circle: None,
            bounding_cylinder: Some(CoarseCylinder {
                center_m: DEVICE_PUMP_OFFSET_M,
                axis: Z,
                radius_m: parts.pump.outer_radius_m,
                half_length_m: 0.5 * parts.pump.cylinder_length_m,
            }),
        },
    ];

    Scene {
        schema: SCENE_SCHEMA.to_string(),
        version: SCENE_VERSION,
        length_unit: "m".to_string(),
        angle_unit: "rad".to_string(),
        frame: DEVICE_SCENE_ID.to_string(),
        design: zero_design(),
        atomistic: AtomisticLayer {
            atom_count: atoms.len(),
            atoms,
        },
        device: DeviceLayer {
            body_count: bodies.len(),
            bodies,
            joints,
            gear_couplings: Vec::new(),
            gear_constraints: Vec::new(),
        },
        coarse: CoarseLayer {
            body_count: coarse.len(),
            bodies: coarse,
        },
    }
}

/// Relaxes every subsystem and returns the scene and the reports.
pub fn build_relaxed_device_scene(
    options: &ValenceRelaxOptions,
) -> Result<(Scene, Vec<String>), SceneError> {
    let mut parts = device_parts()?;
    let mut reports = Vec::new();
    for (name, _) in relax_module_parts(&mut parts.module, options)? {
        reports.push(name);
    }
    reports.push("gas_tank".to_string());
    relax_part(&mut parts.tank.shell, options)?;
    for part in [&mut parts.pump.cylinder, &mut parts.pump.piston] {
        reports.push(part.name.clone());
        relax_part(part, options)?;
    }
    Ok((device_scene(&parts), reports))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sense::sense_report;

    #[test]
    fn the_device_has_its_bodies_and_joints() {
        let scene = device_scene(&device_parts().expect("the device builds"));
        assert_eq!(scene.device.body_count, 29);
        assert_eq!(scene.device.joints.len(), 28);
        let report = sense_report(&scene);
        assert!(report.floating_bodies.is_empty(), "{report:?}");
    }

    #[test]
    #[ignore = "builds and relaxes about 110,000 atoms"]
    fn the_relaxed_device_passes_the_audit() {
        let options = ValenceRelaxOptions {
            max_iterations: 3000,
            ..ValenceRelaxOptions::default()
        };
        let (scene, reports) = build_relaxed_device_scene(&options).expect("the device relaxes");
        assert_eq!(reports.len(), 5, "{reports:?}");
        let report = sense_report(&scene);
        assert_eq!(report.overlap_count, 0, "{report:?}");
        assert!(report.disjoint_bodies.is_empty(), "{report:?}");
        assert!(report.floating_bodies.is_empty(), "{report:?}");
        // Only the module has intended cross-body contacts: twelve tethers and
        // six Fe-N bonds. The tank and the pump touch nothing.
        let contacts = crate::module_scene::body_contacts(&scene, 2.2e-10);
        assert_eq!(contacts.len(), 18, "{contacts:?}");
        assert!(contacts.iter().all(|(a, b, _)| *a < 26 && *b < 26));
    }
}
