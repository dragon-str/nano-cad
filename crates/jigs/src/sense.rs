//! The physical-sense audit of a scene.
//!
//! A scene must look like a real mechanism. This module checks the rules that
//! the viewer and the physics assume:
//!
//! 1. No two atoms overlap. Two atoms closer than `MIN_SEPARATION_M` are a
//!    defect. A non-bonded pair of two bodies closer than `OVERLAP_M` is also a
//!    defect.
//! 2. Every rigid body is one connected lattice. A body whose perceived bond
//!    graph has more than one component moves atoms that are not attached.
//! 3. Every non-fixed body is held by a joint. A body with no joint floats.
//! 4. Two different bodies do not share a perceived bond unless a joint joins
//!    them.
//!
//! Bonds are perceived from geometry: two atoms of the same body closer than
//! `BOND_CUTOFF_M` are bonded. The scene schema does not carry a bond list, so
//! the audit uses the same rule as the relaxation in `nanocad-meter`.

use std::collections::HashMap;

use serde::Serialize;

use crate::scene::Scene;

/// The distance below which two atoms overlap, in metres.
pub const OVERLAP_M: f64 = 1.6e-10;

/// The distance below which two atoms are too close even for a bond, in metres.
///
/// The shortest diamond bond is the C-H bond at 1.09e-10 m. A pair closer than
/// this cannot be a bond, so the pair is a defect.
pub const MIN_SEPARATION_M: f64 = 9.0e-11;

/// The distance below which two atoms of one body are bonded, in metres.
pub const BOND_CUTOFF_M: f64 = 1.62e-10;

/// The cell size of the neighbour grid, in metres.
const CELL_M: f64 = 1.0e-9;

/// One rigid body whose perceived bond graph is disconnected.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct DisjointBody {
    /// The body index.
    pub body: usize,
    /// The body name.
    pub name: String,
    /// The number of atoms of the body.
    pub atoms: usize,
    /// The number of connected components of the bond graph.
    pub components: usize,
}

/// The overlapping atom pairs between one pair of bodies.
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct OverlapBodies {
    /// The first body index.
    pub body_a: usize,
    /// The second body index.
    pub body_b: usize,
    /// The number of overlapping pairs between the two bodies.
    pub count: usize,
}

/// The result of the physical-sense audit.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct SenseReport {
    /// The number of atoms.
    pub atom_count: usize,
    /// The number of bodies.
    pub body_count: usize,
    /// The number of overlapping atom pairs.
    pub overlap_count: usize,
    /// The least distance of an overlapping pair, in metres. Zero when none.
    pub worst_overlap_m: f64,
    /// The overlapping pairs grouped by body pair, most first.
    pub overlap_bodies: Vec<OverlapBodies>,
    /// The non-fixed bodies with no joint.
    pub floating_bodies: Vec<usize>,
    /// The bodies whose bond graph is disconnected.
    pub disjoint_bodies: Vec<DisjointBody>,
    /// The number of perceived bonds between bodies that no joint joins.
    pub cross_body_bonds: usize,
    /// True when the scene obeys every rule.
    pub clean: bool,
}

impl SenseReport {
    /// Returns true when the scene obeys every rule.
    pub fn is_clean(&self) -> bool {
        self.clean
    }
}

fn cell_of(position_m: [f64; 3]) -> (i64, i64, i64) {
    (
        (position_m[0] / CELL_M).floor() as i64,
        (position_m[1] / CELL_M).floor() as i64,
        (position_m[2] / CELL_M).floor() as i64,
    )
}

fn distance_m(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

struct UnionFind {
    parent: Vec<usize>,
}

impl UnionFind {
    fn new(len: usize) -> Self {
        Self {
            parent: (0..len).collect(),
        }
    }

    fn find(&mut self, index: usize) -> usize {
        let mut root = index;
        while self.parent[root] != root {
            root = self.parent[root];
        }
        let mut walk = index;
        while self.parent[walk] != root {
            let next = self.parent[walk];
            self.parent[walk] = root;
            walk = next;
        }
        root
    }

    fn union(&mut self, a: usize, b: usize) {
        let root_a = self.find(a);
        let root_b = self.find(b);
        if root_a != root_b {
            self.parent[root_b] = root_a;
        }
    }
}

/// Audits a scene against the physical-sense rules.
pub fn sense_report(scene: &Scene) -> SenseReport {
    let atoms = &scene.atomistic.atoms;
    let body_count = scene.device.body_count;
    let mut grid: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
    for (index, atom) in atoms.iter().enumerate() {
        grid.entry(cell_of(atom.position_m))
            .or_default()
            .push(index);
    }

    let mut unions = UnionFind::new(atoms.len());
    let mut overlap_count = 0usize;
    let mut worst_overlap_m = f64::INFINITY;
    let mut cross_body_bonds = 0usize;
    let mut overlap_by_pair: HashMap<(usize, usize), usize> = HashMap::new();

    let mut joined: Vec<(usize, usize)> = Vec::new();
    for joint in &scene.device.joints {
        joined.push((
            joint.body_a.min(joint.body_b),
            joint.body_a.max(joint.body_b),
        ));
    }
    let is_joined = |a: usize, b: usize| {
        let pair = (a.min(b), a.max(b));
        joined.contains(&pair)
    };

    for (index, atom) in atoms.iter().enumerate() {
        let (cx, cy, cz) = cell_of(atom.position_m);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let Some(bucket) = grid.get(&(cx + dx, cy + dy, cz + dz)) else {
                        continue;
                    };
                    for &other in bucket {
                        if other <= index {
                            continue;
                        }
                        let distance = distance_m(atom.position_m, atoms[other].position_m);
                        if atom.body == atoms[other].body {
                            if distance < BOND_CUTOFF_M {
                                unions.union(index, other);
                            }
                        } else if distance < BOND_CUTOFF_M
                            && !is_joined(atom.body, atoms[other].body)
                        {
                            cross_body_bonds += 1;
                        }
                        let overlap = distance < MIN_SEPARATION_M
                            || (atom.body != atoms[other].body
                                && distance < OVERLAP_M
                                && !is_joined(atom.body, atoms[other].body));
                        if overlap {
                            overlap_count += 1;
                            if distance < worst_overlap_m {
                                worst_overlap_m = distance;
                            }
                            let pair = (
                                atom.body.min(atoms[other].body),
                                atom.body.max(atoms[other].body),
                            );
                            *overlap_by_pair.entry(pair).or_insert(0) += 1;
                        }
                    }
                }
            }
        }
    }

    let mut roots_by_body: Vec<HashMap<usize, usize>> = vec![HashMap::new(); body_count];
    for (index, atom) in atoms.iter().enumerate() {
        let body = atom.body;
        if body >= body_count {
            continue;
        }
        let root = unions.find(index);
        *roots_by_body[body].entry(root).or_insert(0) += 1;
    }

    let mut disjoint_bodies = Vec::new();
    for (body, roots) in roots_by_body.iter().enumerate() {
        let components = roots.len();
        if components > 1 {
            disjoint_bodies.push(DisjointBody {
                body,
                name: scene.device.bodies[body].name.clone(),
                atoms: atoms.iter().filter(|atom| atom.body == body).count(),
                components,
            });
        }
    }

    let mut floating_bodies = Vec::new();
    for body in 0..body_count {
        let fixed = scene.device.bodies[body].fixed;
        let held = scene
            .device
            .joints
            .iter()
            .any(|joint| joint.body_a == body || joint.body_b == body);
        if !fixed && !held {
            floating_bodies.push(body);
        }
    }

    let clean = overlap_count == 0
        && disjoint_bodies.is_empty()
        && floating_bodies.is_empty()
        && cross_body_bonds == 0;

    let mut overlap_bodies: Vec<OverlapBodies> = overlap_by_pair
        .into_iter()
        .map(|((body_a, body_b), count)| OverlapBodies {
            body_a,
            body_b,
            count,
        })
        .collect();
    overlap_bodies.sort_by_key(|a| std::cmp::Reverse(a.count));

    SenseReport {
        atom_count: atoms.len(),
        body_count,
        overlap_count,
        worst_overlap_m: if overlap_count == 0 {
            0.0
        } else {
            worst_overlap_m
        },
        overlap_bodies,
        floating_bodies,
        disjoint_bodies,
        cross_body_bonds,
        clean,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scene::{
        AtomisticLayer, CoarseLayer, DeviceLayer, SceneAtom, SceneBody, SceneDesign, SceneJoint,
        SCENE_SCHEMA, SCENE_VERSION,
    };

    const BOND_M: f64 = 1.5e-10;
    const CLOSE_M: f64 = 5.0e-11;

    fn atom(position_m: [f64; 3], body: usize) -> SceneAtom {
        SceneAtom {
            element: "C".to_owned(),
            atomic_number: 6,
            position_m,
            body,
            name: format!("body_{body}"),
            part_id: String::new(),
        }
    }

    fn body(index: usize, fixed: bool) -> SceneBody {
        SceneBody {
            index,
            name: format!("body_{index}"),
            role: "test".to_owned(),
            part_id: None,
            position_m: [0.0; 3],
            orientation_wxyz: [1.0, 0.0, 0.0, 0.0],
            mass_kg: 0.0,
            inertia_diagonal_kg_m2: [0.0; 3],
            fixed,
            linear_velocity_m_per_s: [0.0; 3],
            angular_velocity_rad_per_s: [0.0; 3],
        }
    }

    fn joint(body_a: usize, body_b: usize) -> SceneJoint {
        SceneJoint {
            name: format!("joint_{body_a}_{body_b}"),
            kind: "prismatic".to_owned(),
            body_a,
            body_b,
            anchor_world_m: None,
            axis_world: [0.0, 0.0, 1.0],
            axis_a_body: [0.0, 0.0, 1.0],
            axis_b_body: None,
        }
    }

    fn scene(atoms: Vec<SceneAtom>, bodies: Vec<SceneBody>, joints: Vec<SceneJoint>) -> Scene {
        let atom_count = atoms.len();
        let body_count = bodies.len();
        Scene {
            schema: SCENE_SCHEMA.to_owned(),
            version: SCENE_VERSION,
            length_unit: "m".to_owned(),
            angle_unit: "rad".to_owned(),
            frame: String::new(),
            design: SceneDesign {
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
            },
            atomistic: AtomisticLayer { atom_count, atoms },
            device: DeviceLayer {
                body_count,
                bodies,
                joints,
                gear_couplings: Vec::new(),
                gear_constraints: Vec::new(),
            },
            coarse: CoarseLayer {
                body_count,
                bodies: Vec::new(),
            },
        }
    }

    #[test]
    fn a_bonded_body_is_clean() {
        let report = sense_report(&scene(
            vec![atom([0.0; 3], 0), atom([BOND_M, 0.0, 0.0], 0)],
            vec![body(0, true)],
            Vec::new(),
        ));
        assert!(report.is_clean(), "a bonded body is a defect");
        assert_eq!(report.overlap_count, 0);
        assert_eq!(report.worst_overlap_m, 0.0);
    }

    #[test]
    fn two_atoms_too_close_overlap() {
        let report = sense_report(&scene(
            vec![atom([0.0; 3], 0), atom([CLOSE_M, 0.0, 0.0], 0)],
            vec![body(0, true)],
            Vec::new(),
        ));
        assert!(!report.is_clean(), "a close pair is clean");
        assert_eq!(report.overlap_count, 1);
        assert_eq!(report.worst_overlap_m, CLOSE_M);
    }

    #[test]
    fn a_split_body_is_disjoint() {
        let report = sense_report(&scene(
            vec![
                atom([0.0; 3], 0),
                atom([BOND_M, 0.0, 0.0], 0),
                atom([5.0e-9, 0.0, 0.0], 0),
                atom([5.0e-9 + BOND_M, 0.0, 0.0], 0),
            ],
            vec![body(0, true)],
            Vec::new(),
        ));
        assert_eq!(report.disjoint_bodies.len(), 1);
        assert_eq!(report.disjoint_bodies[0].components, 2);
        assert_eq!(report.disjoint_bodies[0].atoms, 4);
    }

    #[test]
    fn a_body_without_a_joint_floats() {
        let report = sense_report(&scene(
            vec![atom([0.0; 3], 0), atom([0.0; 3], 1)],
            vec![body(0, true), body(1, false)],
            Vec::new(),
        ));
        assert_eq!(report.floating_bodies, vec![1]);
    }

    #[test]
    fn a_cross_body_bond_is_a_defect() {
        let report = sense_report(&scene(
            vec![atom([0.0; 3], 0), atom([BOND_M, 0.0, 0.0], 1)],
            vec![body(0, true), body(1, true)],
            Vec::new(),
        ));
        assert_eq!(report.cross_body_bonds, 1);
        assert!(!report.is_clean(), "a cross-body bond is clean");
    }

    #[test]
    fn a_jointed_cross_body_pair_is_allowed() {
        let report = sense_report(&scene(
            vec![atom([0.0; 3], 0), atom([BOND_M, 0.0, 0.0], 1)],
            vec![body(0, true), body(1, false)],
            vec![joint(0, 1)],
        ));
        assert_eq!(report.cross_body_bonds, 0);
        assert_eq!(report.overlap_count, 0);
        assert!(report.floating_bodies.is_empty());
    }
}
