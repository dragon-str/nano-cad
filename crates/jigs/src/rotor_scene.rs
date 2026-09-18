//! A three-scale scene for the sorting rotor.
//!
//! The scene uses the same schema as the planetary gearbox, so one viewer
//! shows both. The two mechanisms differ, so the scene states the difference
//! in its data instead of in a new schema.
//!
//! - The housing is body 0 and it is fixed. The rotor is body 1.
//! - The design block is zero. Its fields describe a gear set, and a rotor has
//!   no gear set. A viewer reads a zero module as "not a gearbox".
//! - The rotor turns about the world `z` axis through the origin, because the
//!   housing chamber port and the rotor axis port share that axis.
//!
//! The device layer carries the joint that a viewer needs to animate the
//! rotor. It carries no gear coupling and no gear constraint.
//!
//! # Scope limit
//!
//! The scene does not include a binding site, a target molecule, a solvent or
//! an ion. It therefore does not show molecular selectivity.

use nanocad_model::Part;
use nanocad_parts::{
    place, ParameterSet, PartError, PartGenerator, RotorHousingGenerator, SortingRotorGenerator,
};

use crate::scene::{
    AtomisticLayer, CoarseBody, CoarseCylinder, CoarseLayer, DeviceLayer, Scene, SceneAtom,
    SceneBody, SceneError, SceneJoint, SCENE_SCHEMA, SCENE_VERSION,
};

/// The device body index of the fixed housing.
pub const HOUSING_BODY: usize = 0;
/// The device body index of the rotating rotor.
pub const ROTOR_BODY: usize = 1;
/// The body role of the fixed housing.
pub const HOUSING_ROLE: &str = "housing";
/// The body role of the rotating rotor.
pub const ROTOR_ROLE: &str = "rotor";

/// A written description of the rotor scene frame.
pub const ROTOR_FRAME: &str =
    "right-handed Cartesian, z is the rotor axis, the rotor turns about z through the origin";

fn scene_atoms(part: &Part, body: usize, name: &str, part_id: &str) -> Vec<SceneAtom> {
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
                position_m,
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

/// Builds the sorting rotor scene.
///
/// The rotor and the housing use their default parameters. The rotor mates to
/// the housing through the rotor axis port and the housing chamber port.
///
/// # Errors
///
/// Returns [`SceneError::Part`] when a generator rejects its defaults, and
/// [`SceneError::MissingAtom`] when an atom has no element or position.
pub fn build_rotor_scene() -> Result<Scene, SceneError> {
    let rotor = SortingRotorGenerator.generate_with_defaults()?;
    let housing = RotorHousingGenerator.generate_with_defaults()?;
    let rotor = place(
        rotor,
        &SortingRotorGenerator.axis_port(),
        &RotorHousingGenerator.chamber_port(),
        0.0,
    )
    .transformed_part();

    let rotor_radius_m = default_m(&SortingRotorGenerator, "disc_radius_m");
    let housing_radius_m = default_m(&RotorHousingGenerator, "chamber_radius_m")
        + default_m(&RotorHousingGenerator, "wall_m");
    let rotor_half_thickness_m = default_m(&SortingRotorGenerator, "thickness_m") / 2.0;
    let housing_half_thickness_m = default_m(&RotorHousingGenerator, "thickness_m") / 2.0;
    let _ = (
        rotor_radius_m,
        housing_radius_m,
        rotor_half_thickness_m,
        housing_half_thickness_m,
    );

    let mut atoms = scene_atoms(&housing, HOUSING_BODY, HOUSING_ROLE, "rotor_housing");
    atoms.extend(scene_atoms(&rotor, ROTOR_BODY, ROTOR_ROLE, "sorting_rotor"));
    if atoms.is_empty() {
        return Err(SceneError::MissingAtom { index: 0 });
    }

    let axis = [0.0, 0.0, 1.0];
    let identity = [1.0, 0.0, 0.0, 0.0];
    let ground = SceneBody {
        index: HOUSING_BODY,
        name: HOUSING_ROLE.to_string(),
        role: HOUSING_ROLE.to_string(),
        part_id: Some("rotor_housing".to_string()),
        position_m: [0.0, 0.0, 0.0],
        orientation_wxyz: identity,
        mass_kg: 0.0,
        inertia_diagonal_kg_m2: [0.0, 0.0, 0.0],
        fixed: true,
        linear_velocity_m_per_s: [0.0, 0.0, 0.0],
        angular_velocity_rad_per_s: [0.0, 0.0, 0.0],
    };
    let body = SceneBody {
        index: ROTOR_BODY,
        name: ROTOR_ROLE.to_string(),
        role: ROTOR_ROLE.to_string(),
        part_id: Some("sorting_rotor".to_string()),
        position_m: [0.0, 0.0, 0.0],
        orientation_wxyz: identity,
        mass_kg: 0.0,
        inertia_diagonal_kg_m2: [0.0, 0.0, 0.0],
        fixed: false,
        linear_velocity_m_per_s: [0.0, 0.0, 0.0],
        angular_velocity_rad_per_s: [0.0, 0.0, 0.0],
    };

    let joint = SceneJoint {
        name: "rotor_axis".to_string(),
        kind: "revolute".to_string(),
        body_a: HOUSING_BODY,
        body_b: ROTOR_BODY,
        anchor_world_m: Some([0.0, 0.0, 0.0]),
        axis_world: axis,
        axis_a_body: axis,
        axis_b_body: Some(axis),
    };

    let coarse = vec![
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
            bounding_cylinder: bounding_cylinder(&rotor, axis, rotor_half_thickness_m),
        },
    ];

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
            body_count: 2,
            bodies: vec![ground, body],
            joints: vec![joint],
            gear_couplings: Vec::new(),
            gear_constraints: Vec::new(),
        },
        coarse: CoarseLayer {
            body_count: coarse.len(),
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
    fn the_rotor_scene_has_two_bodies_and_one_revolute_joint() {
        assert_eq!(scene().device.body_count, 2);
        assert_eq!(scene().coarse.body_count, 2);
        assert_eq!(scene().device.joints.len(), 1);
        let joint = &scene().device.joints[0];
        assert_eq!(joint.kind, "revolute");
        assert_eq!(joint.body_a, HOUSING_BODY);
        assert_eq!(joint.body_b, ROTOR_BODY);
        assert_eq!(joint.axis_world, [0.0, 0.0, 1.0]);
        assert!(scene().device.gear_couplings.is_empty());
        assert!(scene().device.gear_constraints.is_empty());
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
    fn every_atom_belongs_to_one_of_the_two_bodies() {
        let scene = scene();
        assert_eq!(scene.atomistic.atom_count, scene.atomistic.atoms.len());
        let mut housing = 0;
        let mut rotor = 0;
        for atom in &scene.atomistic.atoms {
            match atom.body {
                HOUSING_BODY => housing += 1,
                ROTOR_BODY => rotor += 1,
                other => panic!("atom on unexpected body {other}"),
            }
            assert!(
                atom.atomic_number == 6 || atom.atomic_number == 1,
                "the default parts are carbon with hydrogen caps"
            );
        }
        assert_eq!(rotor, 53964);
        assert_eq!(housing, 39848);
    }

    #[test]
    fn the_coarse_bodies_wrap_the_atoms_of_their_body() {
        let scene = scene();
        let radius_m = |index: usize, atoms: &[SceneAtom]| {
            atoms
                .iter()
                .filter(|atom| atom.body == index)
                .map(|atom| (atom.position_m[0].powi(2) + atom.position_m[1].powi(2)).sqrt())
                .fold(0.0_f64, f64::max)
        };
        for coarse in &scene.coarse.bodies {
            let cylinder = coarse.bounding_cylinder.expect("a cylinder");
            let atoms = radius_m(coarse.index, &scene.atomistic.atoms);
            assert!(
                cylinder.radius_m >= atoms - 1.0e-18,
                "the cylinder radius {} m is below the atom radius {atoms} m",
                cylinder.radius_m
            );
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
