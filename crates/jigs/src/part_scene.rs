//! A one-body scene for a single library part.
//!
//! The part viewer shows one generated part at a time. The scene holds one
//! fixed body, its atoms, and a coarse cylinder that wraps them. The schema is
//! the same schema that the gearbox and rotor scenes use, so the viewer needs
//! no second path.

use serde::Serialize;

use nanocad_model::Part;
use nanocad_parts::{library, ParameterSet};

use crate::rotor_scene::zero_design;
use crate::scene::{
    AtomisticLayer, CoarseBody, CoarseCylinder, CoarseLayer, DeviceLayer, Scene, SceneAtom,
    SceneBody, SceneError, SCENE_SCHEMA, SCENE_VERSION,
};

/// The frame of a single-part scene.
pub const PART_FRAME: &str = "part";

/// The body role of the single part body.
pub const PART_ROLE: &str = "part";

/// The id of the machine scene of the sorting rotor.
pub const ROTOR_SCENE_ID: &str = "rotor";

/// The id of the machine scene of the planetary gearbox.
pub const GEARBOX_SCENE_ID: &str = "gearbox";

/// One entry of the scene catalog that the viewer shows.
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct SceneCatalogEntry {
    /// The stable scene id. The viewer sends it back to the server.
    pub id: String,
    /// The human-readable name.
    pub name: String,
    /// The scene kind: `machine` or `part`.
    pub kind: String,
    /// The library category, or `machine` for a whole mechanism.
    pub category: String,
}

/// Returns every scene that the viewer can show.
///
/// The two machines come first. The library parts follow, in registry order.
pub fn scene_catalog() -> Vec<SceneCatalogEntry> {
    let mut entries = vec![
        SceneCatalogEntry {
            id: ROTOR_SCENE_ID.to_string(),
            name: "Sorting rotor".to_string(),
            kind: "machine".to_string(),
            category: "machines".to_string(),
        },
        SceneCatalogEntry {
            id: GEARBOX_SCENE_ID.to_string(),
            name: "Planetary gearbox".to_string(),
            kind: "machine".to_string(),
            category: "machines".to_string(),
        },
    ];
    for entry in library() {
        entries.push(SceneCatalogEntry {
            id: entry.id.to_string(),
            name: entry.name.to_string(),
            kind: "part".to_string(),
            category: entry.category.to_string(),
        });
    }
    entries
}

/// Builds a one-body scene for the library part `id`.
///
/// An unknown id is [`SceneError::Part`].
pub fn build_part_scene(id: &str) -> Result<Scene, SceneError> {
    let part = nanocad_parts::generate(id, &ParameterSet::new())?;
    let name = nanocad_parts::generator(id)
        .map(|generator| generator.name().to_string())
        .unwrap_or_else(|| id.to_string());
    Ok(scene_of_part(id, &name, &part))
}

/// Wraps one part into a fixed one-body scene.
fn scene_of_part(id: &str, name: &str, part: &Part) -> Scene {
    let mut atoms = Vec::with_capacity(part.atom_count());
    for index in 0..part.atom_count() {
        let element = part.topology.element(index);
        let symbol = element
            .and_then(|element| element.symbol())
            .filter(|symbol| !symbol.is_empty())
            .map(str::to_string)
            .unwrap_or_else(|| {
                let number = element.map(|element| element.atomic_number()).unwrap_or(0);
                format!("Z{number}")
            });
        let atomic_number = element.map(|element| element.atomic_number()).unwrap_or(0);
        let position_m = part.topology.position_m(index).unwrap_or([0.0; 3]);
        atoms.push(SceneAtom {
            element: symbol,
            atomic_number,
            position_m,
            body: 0,
            name: name.to_string(),
            part_id: id.to_string(),
        });
    }
    let atom_count = atoms.len();
    let bounding_cylinder = wrap_along_z(&atoms);
    Scene {
        schema: SCENE_SCHEMA.to_string(),
        version: SCENE_VERSION,
        length_unit: "m".to_string(),
        angle_unit: "rad".to_string(),
        frame: PART_FRAME.to_string(),
        design: zero_design(),
        atomistic: AtomisticLayer { atom_count, atoms },
        device: DeviceLayer {
            body_count: 1,
            bodies: vec![SceneBody {
                index: 0,
                name: name.to_string(),
                role: PART_ROLE.to_string(),
                part_id: Some(id.to_string()),
                position_m: [0.0, 0.0, 0.0],
                orientation_wxyz: [1.0, 0.0, 0.0, 0.0],
                mass_kg: 0.0,
                inertia_diagonal_kg_m2: [0.0, 0.0, 0.0],
                fixed: true,
                linear_velocity_m_per_s: [0.0, 0.0, 0.0],
                angular_velocity_rad_per_s: [0.0, 0.0, 0.0],
            }],
            joints: Vec::new(),
            gear_couplings: Vec::new(),
            gear_constraints: Vec::new(),
        },
        coarse: CoarseLayer {
            body_count: 1,
            bodies: vec![CoarseBody {
                index: 0,
                name: name.to_string(),
                role: PART_ROLE.to_string(),
                pitch_circle: None,
                bounding_cylinder,
            }],
        },
    }
}

/// Returns a coarse cylinder that wraps the atoms along z.
fn wrap_along_z(atoms: &[SceneAtom]) -> Option<CoarseCylinder> {
    if atoms.is_empty() {
        return None;
    }
    let mut radius_m = 0.0_f64;
    let mut min_z_m = f64::INFINITY;
    let mut max_z_m = f64::NEG_INFINITY;
    for atom in atoms {
        let x = atom.position_m[0];
        let y = atom.position_m[1];
        let z = atom.position_m[2];
        radius_m = radius_m.max(x.hypot(y));
        min_z_m = min_z_m.min(z);
        max_z_m = max_z_m.max(z);
    }
    Some(CoarseCylinder {
        center_m: [0.0, 0.0, 0.5 * (min_z_m + max_z_m)],
        axis: [0.0, 0.0, 1.0],
        radius_m,
        half_length_m: 0.5 * (max_z_m - min_z_m),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_catalog_lists_the_machines_and_the_parts() {
        let catalog = scene_catalog();
        assert_eq!(catalog[0].id, ROTOR_SCENE_ID);
        assert_eq!(catalog[0].kind, "machine");
        assert_eq!(catalog[1].id, GEARBOX_SCENE_ID);
        assert_eq!(catalog.len(), library().len() + 2);
        assert!(catalog.iter().any(|entry| entry.id == "leaf_spring"));
    }

    #[test]
    fn a_part_scene_has_one_fixed_body() {
        let scene = build_part_scene("leaf_spring").expect("the leaf spring builds");
        assert_eq!(scene.frame, PART_FRAME);
        assert_eq!(scene.device.body_count, 1);
        assert!(scene.device.joints.is_empty());
        assert!(scene.device.bodies[0].fixed);
        assert_eq!(
            scene.device.bodies[0].part_id.as_deref(),
            Some("leaf_spring")
        );
        assert!(scene.atomistic.atom_count > 0);
        assert_eq!(scene.coarse.body_count, 1);
        assert!(scene.coarse.bodies[0].bounding_cylinder.is_some());
    }

    #[test]
    fn an_unknown_part_is_refused() {
        assert!(build_part_scene("no_such_part").is_err());
    }

    #[test]
    fn a_part_scene_round_trips_through_json() {
        let scene = build_part_scene("cam_hub").expect("the cam hub builds");
        let json = crate::scene::scene_to_json(&scene).expect("the scene writes");
        let restored = crate::scene::scene_from_json(&json).expect("the scene reads");
        assert_eq!(restored.atomistic.atom_count, scene.atomistic.atom_count);
        assert_eq!(restored.device.body_count, 1);
    }
}
