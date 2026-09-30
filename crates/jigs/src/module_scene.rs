//! The respirocyte oxygen-loading module (ADR-0111).
//!
//! Body 0 is the housing (fixed): a ring with end plates, an inlet opening at
//! 345 degrees and an outlet opening at 165 degrees. Body 1 is the rotor,
//! with a revolute joint about z. Bodies 2..2+n are the heme-like sites; each
//! has a fixed joint to the rotor, and four two-carbon tethers covalently
//! bond it to the slot wall when the unit is relaxed (ADR-0113). Bodies
//! 2+n..2+2n are the ligand arms,
//! each with a prismatic joint to the rotor along its pocket radius. Pockets
//! 0-5 show the loaded state (NH2 tip 2.1 A from the iron, O2 bound); pockets
//! 6-11 show the release state (tip withdrawn to 4.0 A, site empty). The
//! cam that moves the arms and the repair donor on the 180-330 degree arc
//! are not built.

use nanocad_meter::{relax_valence, ValenceRelaxOptions, ValenceRelaxReport};
use nanocad_model::{Bond, BondType, Part};
use nanocad_parts::{module_housing, HemeRotorGenerator, HemeRotorSet, ParameterSet};

use crate::rotor_scene::zero_design;
use crate::scene::{
    AtomisticLayer, CoarseBody, CoarseCylinder, CoarseLayer, DeviceLayer, Scene, SceneAtom,
    SceneBody, SceneError, SceneJoint, SCENE_SCHEMA, SCENE_VERSION,
};

/// The scene id of the module.
pub const MODULE_SCENE_ID: &str = "respirocyte_module";

/// The gap between the rotor rim carbons and the housing, in metres.
pub const MODULE_CLEARANCE_M: f64 = 4.6e-10;

const Z: [f64; 3] = [0.0, 0.0, 1.0];

/// The generated bodies of the module, in body order.
pub struct ModuleParts {
    /// The rotor set (rotor, sites, arms).
    pub set: HemeRotorSet,
    /// The housing.
    pub housing: Part,
}

/// Generates the module parts, unrelaxed.
pub fn module_parts() -> Result<ModuleParts, SceneError> {
    let set = HemeRotorGenerator.build(&ParameterSet::new())?;
    let half = set.half_thickness_m;
    let housing = module_housing(
        set.disc_radius_m + MODULE_CLEARANCE_M,
        0.8e-9,
        half + MODULE_CLEARANCE_M,
        0.6e-9,
        1.6e-9,
        -15.0_f64.to_radians(),
    )?;
    Ok(ModuleParts { set, housing })
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
    let mut atoms = atoms_of(&parts.housing, 0, "module_housing");
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
            &format!("site_{k}_attachment_placeholder"),
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
        max_iterations: 3000,
        ..ValenceRelaxOptions::default()
    }
}

/// Builds the unrelaxed module scene.
pub fn build_module_scene() -> Result<Scene, SceneError> {
    Ok(module_scene(&module_parts()?))
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

/// Relaxes the module parts in place. The housing relaxes alone. The rotor,
/// the sites and the arms relax together as one system, with the Fe-N bond of
/// each loaded pocket and every attachment tether included. Returns the
/// reports, housing first then the rotor unit.
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
    let report = relax_valence(&unit, 0..unit.atom_count(), options)
        .map_err(|e| SceneError::Part(format!("the rotor unit did not relax: {e}")))?;
    let mut cursor = 0usize;
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
        let scene = build_module_scene().expect("the module builds");
        assert_eq!(scene.device.body_count, 26);
        assert_eq!(scene.device.joints.len(), 25);
        let report = sense_report(&scene);
        assert!(report.floating_bodies.is_empty(), "{report:?}");
        // The six loaded pockets carry the Fe-NH2 ligand bond across bodies.
        assert_eq!(report.cross_body_bonds, 6);
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
    /// piece, and no contact under 2.2 A between bodies except the six Fe-N
    /// ligand bonds. About a minute in a release build; run with
    /// `cargo test --release -p nanocad-jigs -- --ignored relaxed_module`.
    #[test]
    #[ignore]
    fn the_relaxed_module_passes_the_audit() {
        let (scene, reports) =
            build_relaxed_module_scene(&module_relax_options()).expect("the module relaxes");
        assert!(reports.iter().all(|(_, r)| r.converged));
        let report = sense_report(&scene);
        assert_eq!(report.overlap_count, 0, "{report:?}");
        assert!(report.disjoint_bodies.is_empty(), "{report:?}");
        // Twelve tether bonds (one per site) plus the six Fe-N ligand bonds.
        let close: Vec<_> = body_contacts(&scene, 2.15e-10);
        assert_eq!(close.len(), 18, "{close:?}");
        for (a, b, d) in &close {
            if *a == 1 && (2..14).contains(b) {
                assert!((1.35e-10..=1.65e-10).contains(d), "tether {a}-{b} {d}");
            }
        }
    }
}
