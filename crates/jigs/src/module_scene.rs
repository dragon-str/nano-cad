//! The respirocyte oxygen-loading module in a hull-wall section (ADR-0111,
//! ADR-0117).
//!
//! Body 0 is the hull wall (fixed): a diamond slab on the plane x = 0 with
//! blood on the +x side and the tank on the -x side, two end plates over the
//! rotor faces carrying the arm-cam grooves (ADR-0116), the two seal arcs
//! where the rim crosses the wall, and one hydroxyl-lined drain bore. Body 1
//! is the rotor, with a revolute joint about z (its axis lies in the wall
//! plane). Bodies 2..2+n are the heme-like sites; each has a fixed joint to
//! the rotor, and four two-carbon tethers covalently bond it to the slot wall
//! when the unit is relaxed (ADR-0113). Bodies 2+n..2+2n are the ligand arms,
//! each with a prismatic joint to the rotor along its pocket radius.
//!
//! The cam profile sets each pocket's state. The rotor turns +z: a pocket
//! loads in the blood with the arm in, crosses the seal at 90 degrees with
//! the arm in, releases in the tank after the arm comes out at 120 degrees,
//! and returns through the seal at 270 degrees, where the arm goes in again.
//! The repair donor and the discharge rotor are not built.

use nanocad_meter::{relax_valence, ValenceRelaxOptions, ValenceRelaxReport};
use nanocad_model::{Bond, BondType, Part};
use nanocad_parts::{
    hull_wall, water_rotor, HemeRotorGenerator, HemeRotorSet, HullWall, HullWallSpec, ParameterSet,
    WaterRotorSpec,
};

use crate::rotor_scene::zero_design;
use crate::scene::{
    AtomisticLayer, CoarseBody, CoarseCylinder, CoarseLayer, DeviceLayer, Scene, SceneAtom,
    SceneBody, SceneError, SceneJoint, SCENE_SCHEMA, SCENE_VERSION,
};

/// The scene id of the module.
pub const MODULE_SCENE_ID: &str = "respirocyte_module";

/// The scene id of the water-return rotor in its hull wall (ADR-0122).
pub const WATER_ROTOR_SCENE_ID: &str = "respirocyte_water_rotor";

const Z: [f64; 3] = [0.0, 0.0, 1.0];

/// The generated bodies of the module, in body order.
pub struct ModuleParts {
    /// The rotor set (rotor, sites, arms).
    pub set: HemeRotorSet,
    /// The hull wall that holds the rotor (body 0).
    pub housing: Part,
    /// The hull-wall geometry.
    pub wall: HullWall,
}

/// Generates the module parts, unrelaxed.
pub fn module_parts() -> Result<ModuleParts, SceneError> {
    module_parts_with(&ParameterSet::new())
}

/// Generates the module parts with rotor parameters, unrelaxed.
pub fn module_parts_with(parameters: &ParameterSet) -> Result<ModuleParts, SceneError> {
    let set = HemeRotorGenerator.build(parameters)?;
    let wall = hull_wall(&set, &HullWallSpec::default())?;
    let housing = wall.part.clone();
    Ok(ModuleParts { set, housing, wall })
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

/// Builds the module scene from its parts.
pub fn module_scene(parts: &ModuleParts) -> Scene {
    let set = &parts.set;
    let n = set.sites.len();
    let mut atoms = atoms_of(&parts.housing, 0, "hull_wall");
    atoms.extend(atoms_of(&set.rotor, 1, "heme_rotor"));
    let mut bodies = vec![
        body(0, "housing", "hull_wall", true),
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
            &format!("site_{k}_attachment_placeholder"),
            "fixed",
            1,
            2 + k,
            fe,
            Z,
        ));
    }
    for k in 0..set.arms.len() {
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
    let coarse = vec![CoarseBody {
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
    }];
    Scene {
        schema: SCENE_SCHEMA.to_string(),
        version: SCENE_VERSION,
        length_unit: "m".to_string(),
        angle_unit: "rad".to_string(),
        frame: MODULE_SCENE_ID.to_string(),
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

/// The relaxation settings of the module: run to convergence (about 1,900
/// iterations for the housing, 700-1,200 for the rotor unit).
pub fn module_relax_options() -> ValenceRelaxOptions {
    ValenceRelaxOptions {
        max_iterations: 6000,
        ..ValenceRelaxOptions::default()
    }
}

/// Generates the water-return rotor in its hull wall, unrelaxed (ADR-0122).
pub fn water_module_parts() -> Result<ModuleParts, SceneError> {
    let set = water_rotor(&WaterRotorSpec::default())?;
    let wall = hull_wall(&set, &HullWallSpec::default())?;
    let housing = wall.part.clone();
    Ok(ModuleParts { set, housing, wall })
}

/// Builds the relaxed water-return rotor scene.
pub fn build_relaxed_water_scene(
    options: &ValenceRelaxOptions,
) -> Result<(Scene, Vec<(String, ValenceRelaxReport)>), SceneError> {
    let mut parts = water_module_parts()?;
    let reports = relax_module_parts(&mut parts, options)?;
    let mut scene = module_scene(&parts);
    scene.frame = WATER_ROTOR_SCENE_ID.to_string();
    Ok((scene, reports))
}

/// Builds the unrelaxed module scene.
pub fn build_module_scene() -> Result<Scene, SceneError> {
    Ok(module_scene(&module_parts()?))
}

/// Builds the module scene with rotor parameters, relaxed or not.
pub fn build_module_scene_with(
    parameters: &ParameterSet,
    relax: Option<&ValenceRelaxOptions>,
) -> Result<(Scene, Vec<(String, ValenceRelaxReport)>), SceneError> {
    let mut parts = module_parts_with(parameters)?;
    let reports = match relax {
        Some(options) => relax_module_parts(&mut parts, options)?,
        None => Vec::new(),
    };
    Ok((module_scene(&parts), reports))
}

pub(crate) fn relax_part(
    part: &mut Part,
    options: &ValenceRelaxOptions,
) -> Result<ValenceRelaxReport, SceneError> {
    let report = relax_valence(&part.topology, 0..part.atom_count(), options)
        .map_err(|e| SceneError::Part(format!("{} did not relax: {e}", part.name)))?;
    for (i, p) in report.positions_m.iter().enumerate() {
        part.topology
            .set_position_m(i, *p)
            .map_err(|e| SceneError::Part(e.to_string()))?;
    }
    Ok(report)
}

/// Relaxes the module. The housing relaxes alone. The rotor, the sites and
/// the arms relax together as one system, with the Fe-N bond of each loaded
/// pocket included, so the pockets stay fitted to their sites. Returns the
/// scene and the reports (housing, then the rotor unit).
pub fn build_relaxed_module_scene(
    options: &ValenceRelaxOptions,
) -> Result<(Scene, Vec<(String, ValenceRelaxReport)>), SceneError> {
    let mut parts = module_parts()?;
    let reports = relax_module_parts(&mut parts, options)?;
    Ok((module_scene(&parts), reports))
}

/// The housing atoms within this distance of the rotor unit form the fixed
/// wall of the second relaxation stage, in metres.
pub const MODULE_WALL_REACH_M: f64 = 1.0e-9;

/// Relaxes the module parts in place, in two stages. The housing relaxes
/// alone. Then the rotor, the sites and the arms relax as one system, with
/// the Fe-N bond of each loaded pocket and every attachment tether, against
/// the relaxed housing atoms near them held fixed. The fixed wall matters
/// because the cam grooves hold each arm pin at contact (ADR-0116): relaxed
/// apart, the pins bent and the released arms drifted into the groove walls.
/// Returns the reports, housing first then the rotor unit.
pub fn relax_module_parts(
    parts: &mut ModuleParts,
    options: &ValenceRelaxOptions,
) -> Result<Vec<(String, ValenceRelaxReport)>, SceneError> {
    let mut reports = Vec::new();
    let name = parts.housing.name.clone();
    reports.push((name, relax_part(&mut parts.housing, options)?));

    // Merge rotor, sites and arms into one topology.
    let set = &mut parts.set;
    let mut unit = set.rotor.topology.clone();
    let mut offsets = Vec::new();
    for part in set.sites.iter().chain(set.arms.iter()) {
        let off = unit.atom_count() as u32;
        offsets.push(off);
        for i in 0..part.atom_count() {
            if let Some(atom) = part.topology.atom(i) {
                unit.add_atom(atom);
            }
        }
        for b in part.topology.bonds() {
            unit.add_bond(Bond::new(b.u + off, b.v + off, b.order, b.bond_type))
                .map_err(|e| SceneError::Part(e.to_string()))?;
        }
    }
    let n = set.sites.len();
    for k in 0..n {
        if set.loaded[k] {
            let fe = offsets[k] + set.site_iron[k] as u32;
            let tip = offsets[n + k] + set.arm_tip[k] as u32;
            unit.add_bond(Bond::new(fe, tip, 1, BondType::Single))
                .map_err(|e| SceneError::Part(e.to_string()))?;
        }
    }
    // The covalent attachment: bond each site tether end carbon to its rotor
    // wall carbon (ADR-0113).
    for (base, anchors) in offsets.iter().zip(set.site_anchor.iter()) {
        let base = *base as usize;
        for &(si, rj) in anchors {
            unit.add_bond(Bond::new(
                (base + si) as u32,
                rj as u32,
                1,
                BondType::Single,
            ))
            .map_err(|e| SceneError::Part(e.to_string()))?;
        }
    }

    // The fixed wall: the housing atoms near the unit, placed first.
    let wall = near_atoms(&parts.housing, &unit, MODULE_WALL_REACH_M);
    let mut system = nanocad_model::Topology::new();
    for &i in &wall {
        if let Some(atom) = parts.housing.topology.atom(i) {
            system.add_atom(atom);
        }
    }
    let off = system.atom_count() as u32;
    for i in 0..unit.atom_count() {
        if let Some(atom) = unit.atom(i) {
            system.add_atom(atom);
        }
    }
    for b in unit.bonds() {
        system
            .add_bond(Bond::new(b.u + off, b.v + off, b.order, b.bond_type))
            .map_err(|e| SceneError::Part(e.to_string()))?;
    }
    let staged = ValenceRelaxOptions {
        frozen_atoms: wall.len(),
        ..*options
    };
    let report = relax_valence(&system, 0..system.atom_count(), &staged)
        .map_err(|e| SceneError::Part(format!("the rotor unit did not relax: {e}")))?;
    let mut cursor = wall.len();
    let parts_mut: Vec<&mut Part> = std::iter::once(&mut set.rotor)
        .chain(set.sites.iter_mut())
        .chain(set.arms.iter_mut())
        .collect();
    for part in parts_mut {
        for i in 0..part.atom_count() {
            part.topology
                .set_position_m(i, report.positions_m[cursor])
                .map_err(|e| SceneError::Part(e.to_string()))?;
            cursor += 1;
        }
    }
    reports.push(("rotor_sites_arms".to_string(), report));
    Ok(reports)
}

/// The atoms of `part` within `reach_m` of any atom of `other`.
fn near_atoms(part: &Part, other: &nanocad_model::Topology, reach_m: f64) -> Vec<usize> {
    use std::collections::HashSet;
    let cell = |p: [f64; 3]| {
        (
            (p[0] / reach_m).floor() as i64,
            (p[1] / reach_m).floor() as i64,
            (p[2] / reach_m).floor() as i64,
        )
    };
    let mut grid: std::collections::HashMap<(i64, i64, i64), Vec<[f64; 3]>> =
        std::collections::HashMap::new();
    for i in 0..other.atom_count() {
        if let Some(p) = other.position_m(i) {
            grid.entry(cell(p)).or_default().push(p);
        }
    }
    let mut out = HashSet::new();
    for i in 0..part.atom_count() {
        let Some(p) = part.topology.position_m(i) else {
            continue;
        };
        let (x, y, z) = cell(p);
        'search: for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    for q in grid.get(&(x + dx, y + dy, z + dz)).into_iter().flatten() {
                        let d2 =
                            (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2);
                        if d2 <= reach_m * reach_m {
                            out.insert(i);
                            break 'search;
                        }
                    }
                }
            }
        }
    }
    let mut v: Vec<usize> = out.into_iter().collect();
    v.sort_unstable();
    v
}

/// The least distance between every pair of bodies, joints included, as
/// (body_a, body_b, distance_m) for pairs closer than `limit_m`.
pub fn body_contacts(scene: &Scene, limit_m: f64) -> Vec<(usize, usize, f64)> {
    use std::collections::HashMap;
    let atoms = &scene.atomistic.atoms;
    let cell = |p: [f64; 3]| {
        (
            (p[0] / limit_m).floor() as i64,
            (p[1] / limit_m).floor() as i64,
            (p[2] / limit_m).floor() as i64,
        )
    };
    let mut grid: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
    for (i, a) in atoms.iter().enumerate() {
        grid.entry(cell(a.position_m)).or_default().push(i);
    }
    let mut best: HashMap<(usize, usize), f64> = HashMap::new();
    for (i, a) in atoms.iter().enumerate() {
        let (x, y, z) = cell(a.position_m);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    for &j in grid.get(&(x + dx, y + dy, z + dz)).into_iter().flatten() {
                        let b = &atoms[j];
                        if j <= i || a.body == b.body {
                            continue;
                        }
                        let d = ((a.position_m[0] - b.position_m[0]).powi(2)
                            + (a.position_m[1] - b.position_m[1]).powi(2)
                            + (a.position_m[2] - b.position_m[2]).powi(2))
                        .sqrt();
                        if d < limit_m {
                            let key = (a.body.min(b.body), a.body.max(b.body));
                            let e = best.entry(key).or_insert(f64::INFINITY);
                            *e = e.min(d);
                        }
                    }
                }
            }
        }
    }
    let mut out: Vec<_> = best.into_iter().map(|((a, b), d)| (a, b, d)).collect();
    out.sort_by(|x, y| x.2.total_cmp(&y.2));
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sense::sense_report;

    #[test]
    fn the_module_has_its_bodies_and_joints() {
        let parts = module_parts().expect("the module builds");
        let scene = module_scene(&parts);
        assert_eq!(scene.device.body_count, 18);
        assert_eq!(scene.device.joints.len(), 17);
        let report = sense_report(&scene);
        assert!(report.floating_bodies.is_empty(), "{report:?}");
        // The loaded pockets (five of eight, set by the cam profile) carry the
        // Fe-NH2 ligand bond across bodies.
        let loaded = parts.set.loaded.iter().filter(|&&l| l).count();
        assert_eq!(loaded, 5);
        assert_eq!(report.cross_body_bonds, loaded);
    }

    #[test]
    fn every_site_has_four_attachment_tethers() {
        let set = HemeRotorGenerator
            .build(&ParameterSet::new())
            .expect("the rotor builds");
        assert_eq!(set.site_anchor.len(), set.sites.len());
        for (k, anchors) in set.site_anchor.iter().enumerate() {
            assert_eq!(anchors.len(), 4, "site {k} has tethers");
            for &(site_atom, rotor_atom) in anchors {
                assert!(site_atom < set.sites[k].atom_count());
                assert!(rotor_atom < set.rotor.atom_count());
            }
        }
    }

    /// The relaxed module is the delivered CAD model: no overlap, no loose
    /// piece, and no contact under 2.2 A between bodies except the tethers, the Fe-N
    /// ligand bonds and the seal rubs. A few minutes in a release build; run with
    /// `cargo test --release -p nanocad-jigs -- --ignored relaxed_module`.
    #[test]
    #[ignore]
    fn the_relaxed_module_passes_the_audit() {
        let mut parts = module_parts().expect("parts");
        let reports =
            relax_module_parts(&mut parts, &module_relax_options()).expect("the module relaxes");
        assert!(reports.iter().all(|(_, r)| r.converged));
        // The cam states survive the relaxation (ADR-0116): each loaded arm
        // stays on its iron, each released arm at 4.0 A or more.
        let set = &parts.set;
        for k in 0..set.sites.len() {
            let fe = set.sites[k]
                .topology
                .position_m(set.site_iron[k])
                .expect("fe");
            let n = set.arms[k].topology.position_m(set.arm_tip[k]).expect("n");
            let d =
                ((fe[0] - n[0]).powi(2) + (fe[1] - n[1]).powi(2) + (fe[2] - n[2]).powi(2)).sqrt();
            if set.loaded[k] {
                assert!((d - 2.1e-10).abs() < 0.05e-10, "pocket {k} Fe-N {d:e}");
            } else {
                assert!(d >= 4.0e-10, "pocket {k} Fe-N {d:e}");
            }
        }
        let scene = module_scene(&parts);
        let report = sense_report(&scene);
        assert_eq!(report.overlap_count, 0, "{report:?}");
        assert!(report.disjoint_bodies.is_empty(), "{report:?}");
        // One tether contact per site plus one Fe-N bond per loaded pocket.
        // The wall touches the rotor only at the seals (the two rim arcs and
        // the face seal rings), which are set at H-H contact (ADR-0117).
        let seal = parts.wall.seal_half_rad;
        let wall: Vec<[f64; 3]> = (0..parts.housing.atom_count())
            .filter_map(|i| parts.housing.topology.position_m(i))
            .collect();
        let mut rubs = 0usize;
        for i in 0..parts.set.rotor.atom_count() {
            let p = parts.set.rotor.topology.position_m(i).expect("p");
            for q in &wall {
                let d =
                    ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt();
                if d < 2.15e-10 {
                    let theta = p[1].atan2(p[0]);
                    let off = (theta.abs() - 0.5 * std::f64::consts::PI).abs();
                    eprintln!(
                        "wall-rotor {:.2} A at {:.1} deg, r {:.3} nm, z {:.3} nm",
                        d * 1e10,
                        theta.to_degrees(),
                        p[0].hypot(p[1]) * 1e9,
                        p[2] * 1e9
                    );
                    let at_rim_seal = off <= seal + 0.05;
                    let at_face_seal = p[0].hypot(p[1])
                        >= parts.wall.face_seal_inner_radius_m - 2.0e-10
                        && p[2].abs() >= parts.set.half_thickness_m - 1.0e-10;
                    assert!(at_rim_seal || at_face_seal, "a rub outside the seals");
                    assert!(d > 1.9e-10, "the seal is pressed: {d:e}");
                    rubs += 1;
                }
            }
        }
        let n = parts.set.sites.len();
        // An arm pin rides its groove at contact (ADR-0116): a wall-arm pair
        // at 2.0 A or more is the cam, not a clash.
        let cam = |a: usize, b: usize, d: f64| a == 0 && b >= 2 + n && d >= 2.0e-10;
        let close: Vec<_> = body_contacts(&scene, 2.15e-10)
            .into_iter()
            .filter(|&(a, b, _)| !(a == 0 && b == 1 && rubs > 0))
            .filter(|&(a, b, d)| !cam(a, b, d))
            .collect();
        let loaded = parts.set.loaded.iter().filter(|&&l| l).count();
        assert_eq!(close.len(), n + loaded, "{close:?}");
        for (a, b, d) in &close {
            if *a == 1 && (2..2 + n).contains(b) {
                assert!((1.35e-10..=1.65e-10).contains(d), "tether {a}-{b} {d}");
            }
        }
    }

    #[test]
    fn the_water_rotor_scene_has_two_bodies() {
        let parts = water_module_parts().expect("builds");
        let scene = module_scene(&parts);
        assert_eq!(scene.device.body_count, 8);
        assert_eq!(scene.device.joints.len(), 7);
    }

    /// The relaxed water rotor: no overlap, no loose body; the wall touches
    /// the rotor only at the seals.
    #[test]
    #[ignore]
    fn the_relaxed_water_rotor_passes_the_audit() {
        let (scene, reports) = build_relaxed_water_scene(&module_relax_options()).expect("relaxes");
        assert!(reports.iter().all(|(_, r)| r.converged), "{reports:?}");
        let report = sense_report(&scene);
        assert_eq!(report.overlap_count, 0, "{report:?}");
        assert!(report.disjoint_bodies.is_empty(), "{report:?}");
        assert!(report.floating_bodies.is_empty(), "{report:?}");
    }
}
