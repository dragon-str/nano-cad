//! One pocket of the site-blocking rotor (ADR-0090).
//!
//! The pocket is a short bore in a diamond block. The bore axis is `z` and
//! the mouth is the top face of the block. One iron atom replaces a bore-wall
//! carbon near the floor, so the iron keeps three lattice bonds and its open
//! side faces the bore. The iron is level with the carbon wall, one cap
//! length behind the hydrogen surface, so a passing rod does not strike it.
//!
//! A rod slides in the bore. Its end face is the pocket floor. Pulled back,
//! the rod leaves the iron open and the pocket loads. Pushed out, the rod
//! covers the iron: the rod-to-wall gap is too narrow for a bound guest, so
//! the rod drives the guest out of the mouth. The same narrow gap is the seal
//! against water along the rod.
//!
//! This part is a test cell, not the full rotor. It models one pocket, a
//! rigid lattice and no solvent.

use std::collections::HashMap;

use nanocad_model::{Atom, Bond, BondType, Element, Part, Topology};
use nanocad_units::Unit;

use crate::diamond_solid;
use crate::error::PartError;
use crate::gear_pump::keep_well_bonded;
use crate::generator::PartGenerator;
use crate::lattice_fill::fill_solid;
use crate::parameter::{ParameterSet, ParameterSpec};
use crate::shape::{Box3, Cylinder, Difference};

static ROD_POCKET_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new(
        "rod_radius_m",
        Some(Unit::Metre),
        2.6e-10,
        1.5e-10,
        2.0e-9,
        false,
        "radius of the rod lattice, carbon centres",
    ),
    ParameterSpec::new(
        "clearance_m",
        Some(Unit::Metre),
        3.0e-10,
        2.5e-10,
        1.0e-9,
        false,
        "carbon-to-carbon gap between the rod and the bore wall",
    ),
    ParameterSpec::new(
        "depth_m",
        Some(Unit::Metre),
        3.5e-10,
        2.0e-10,
        2.0e-9,
        false,
        "height of the bore above the iron, up to the mouth",
    ),
    ParameterSpec::new(
        "below_m",
        Some(Unit::Metre),
        1.0e-9,
        4.0e-10,
        3.0e-9,
        false,
        "height of the block below the iron, where the rod runs",
    ),
    ParameterSpec::new(
        "wall_m",
        Some(Unit::Metre),
        6.0e-10,
        3.0e-10,
        3.0e-9,
        false,
        "wall thickness of the block around the bore",
    ),
    ParameterSpec::new(
        "iron_rank",
        None,
        0.0,
        0.0,
        64.0,
        true,
        "which ranked bore-wall carbon becomes the iron; 0 points most into the bore",
    ),
    ParameterSpec::new(
        "rod_length_m",
        Some(Unit::Metre),
        1.4e-9,
        4.0e-10,
        5.0e-9,
        false,
        "length of the rod",
    ),
];

/// A generated pocket: the block with its iron, and the rod in its own frame.
#[derive(Clone, Debug, PartialEq)]
pub struct RodPocketSet {
    /// The block with the bore and the iron. The iron sits at `z` near zero.
    pub block: Part,
    /// The rod, with its end face at `z = 0` and its body below.
    pub rod: Part,
    /// The index of the iron atom in `block`.
    pub iron: usize,
    /// The bore radius, carbon centres, in metres.
    pub bore_radius_m: f64,
    /// The rod radius, carbon centres, in metres.
    pub rod_radius_m: f64,
    /// The `z` of the iron, in metres.
    pub iron_z_m: f64,
    /// The `z` of the pocket mouth (block top), in metres.
    pub mouth_z_m: f64,
    /// The `z` of the highest rod carbon above the rod frame origin. The
    /// face of the rod, in the rod frame, in metres.
    pub rod_face_offset_m: f64,
}

impl RodPocketSet {
    /// Returns the block and the rod in one part, with the rod face at
    /// `face_z_m`. The rod atoms follow the block atoms.
    pub fn assembled(&self, face_z_m: f64) -> Result<Part, PartError> {
        let mut topology = self.block.topology.clone();
        let offset = topology.atom_count() as u32;
        let shift_m = face_z_m - self.rod_face_offset_m;
        for index in 0..self.rod.atom_count() {
            let Some(mut atom) = self.rod.topology.atom(index) else {
                continue;
            };
            atom.position_m[2] += shift_m;
            topology.add_atom(atom);
        }
        for bond in self.rod.topology.bonds() {
            topology.add_bond(Bond::new(
                bond.u + offset,
                bond.v + offset,
                bond.order,
                bond.bond_type,
            ))?;
        }
        Ok(Part::new("rod_pocket_assembly", topology).with_material("diamond"))
    }

    /// The number of atoms in the block. Rod atoms start here in
    /// [`RodPocketSet::assembled`].
    pub fn block_atom_count(&self) -> usize {
        self.block.atom_count()
    }
}

/// Builds one pocket of the site-blocking rotor.
#[derive(Clone, Copy, Debug, Default)]
pub struct RodPocketGenerator;

impl RodPocketGenerator {
    /// Resolves the inputs and builds the block and the rod.
    pub fn build(&self, parameters: &ParameterSet) -> Result<RodPocketSet, PartError> {
        let resolved = self.resolve(parameters)?;
        let rod_radius_m = resolved.require("rod_radius_m")?;
        let clearance_m = resolved.require("clearance_m")?;
        let depth_m = resolved.require("depth_m")?;
        let below_m = resolved.require("below_m")?;
        let wall_m = resolved.require("wall_m")?;
        let rod_length_m = resolved.require("rod_length_m")?;
        let bore_radius_m = rod_radius_m + clearance_m;
        let half_width_m = bore_radius_m + wall_m;

        // The block, with the bore through it.
        let block_solid = Difference {
            outer: Box::new(Box3 {
                min_m: [-half_width_m, -half_width_m, -below_m],
                max_m: [half_width_m, half_width_m, depth_m],
            }),
            inner: Box::new(Cylinder {
                center_m: [0.0, 0.0, 0.5 * (depth_m - below_m)],
                radius_m: bore_radius_m,
                height_m: below_m + depth_m + 1.0e-9,
            }),
        };
        let mut block = Topology::new();
        let carbons: Vec<u32> = keep_well_bonded(fill_solid(&block_solid))
            .into_iter()
            .map(|site| block.add_atom(Atom::new(Element::CARBON, site, 0.0, "C")))
            .collect();
        if carbons.is_empty() {
            return Err(PartError::InvalidGeometry(
                "the pocket block has no lattice site".to_owned(),
            ));
        }
        diamond_solid::bond_and_cap(&mut block, &carbons)?;

        // The iron: a bore-wall carbon with exactly one cap.
        let mut caps: HashMap<usize, Vec<usize>> = HashMap::new();
        for bond in block.bonds() {
            let (u, v) = (bond.u as usize, bond.v as usize);
            if block.element(v) == Some(Element::HYDROGEN) {
                caps.entry(u).or_default().push(v);
            }
            if block.element(u) == Some(Element::HYDROGEN) {
                caps.entry(v).or_default().push(u);
            }
        }
        // The best wall carbon has its open side (the cap direction) pointing
        // most directly at the bore axis, so a bound guest points into the
        // bore and not along the wall. It must lie near z = 0.
        let score = |index: usize| -> Option<f64> {
            let p = block.position_m(index)?;
            let radius_m = p[0].hypot(p[1]);
            if radius_m > bore_radius_m + 1.0e-10 || p[2].abs() > 4.0e-10 || radius_m <= 0.0 {
                return None;
            }
            let h = block.position_m(*caps.get(&index)?.first()?)?;
            let d = [h[0] - p[0], h[1] - p[1], h[2] - p[2]];
            let length = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            let inward = -(d[0] * p[0] + d[1] * p[1]) / (length * radius_m);

            Some(-inward + 0.1 * p[2].abs() / 1.0e-10)
        };
        // The candidates are ranked by that score. A neighbouring cap can
        // still crowd the open side, so the caller measures the binding at a
        // few ranks and picks one with `iron_rank`.
        let mut ranked: Vec<(f64, usize)> = carbons
            .iter()
            .map(|&c| c as usize)
            .filter(|c| caps.get(c).is_some_and(|list| list.len() == 1))
            .filter_map(|c| score(c).map(|s| (s, c)))
            .collect();
        ranked.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        let rank = resolved.require("iron_rank")?.round() as usize;
        let (iron, cap) = ranked
            .get(rank)
            .map(|&(_, c)| (c, caps[&c][0]))
            .ok_or_else(|| {
                PartError::InvalidGeometry(format!(
                    "no bore-wall carbon at rank {rank} can take the iron ({} candidates)",
                    ranked.len()
                ))
            })?;
        let (block, iron) = without_atom(&block, cap, iron)?;
        let mut block = block;
        block.set_element(iron, Element::IRON)?;
        block.set_atom_type(iron, "Fe")?;
        let iron_z_m = block.position_m(iron).map_or(0.0, |p| p[2]);

        // The rod, face up, in its own frame.
        let rod_solid = Cylinder {
            center_m: [0.0, 0.0, -0.5 * rod_length_m],
            radius_m: rod_radius_m,
            height_m: rod_length_m,
        };
        let mut rod = Topology::new();
        let rod_carbons: Vec<u32> = keep_well_bonded(fill_solid(&rod_solid))
            .into_iter()
            .map(|site| rod.add_atom(Atom::new(Element::CARBON, site, 0.0, "C")))
            .collect();
        if rod_carbons.is_empty() {
            return Err(PartError::InvalidGeometry(
                "the rod has no lattice site".to_owned(),
            ));
        }
        diamond_solid::bond_and_cap(&mut rod, &rod_carbons)?;
        let rod_face_offset_m = rod_carbons
            .iter()
            .filter_map(|&c| rod.position_m(c as usize))
            .map(|p| p[2])
            .fold(f64::NEG_INFINITY, f64::max);

        Ok(RodPocketSet {
            block: Part::new("rod_pocket_block", block).with_material("diamond"),
            rod: Part::new("rod_pocket_rod", rod).with_material("diamond"),
            iron,
            bore_radius_m,
            rod_radius_m,
            iron_z_m,
            mouth_z_m: depth_m,
            rod_face_offset_m,
        })
    }
}

impl PartGenerator for RodPocketGenerator {
    fn id(&self) -> &'static str {
        "rod_pocket"
    }

    fn name(&self) -> &'static str {
        "Site-blocking pocket"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        ROD_POCKET_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        let set = self.build(parameters)?;
        set.assembled(set.iron_z_m - 2.5e-10)
    }
}

/// Copies a topology without one atom. Returns the copy and the new index of
/// `keep`.
fn without_atom(
    topology: &Topology,
    drop: usize,
    keep: usize,
) -> Result<(Topology, usize), PartError> {
    let mut out = Topology::new();
    let mut map = vec![u32::MAX; topology.atom_count()];
    for (index, slot) in map.iter_mut().enumerate() {
        if index == drop {
            continue;
        }
        if let Some(atom) = topology.atom(index) {
            *slot = out.add_atom(atom);
        }
    }
    for bond in topology.bonds() {
        let (u, v) = (map[bond.u as usize], map[bond.v as usize]);
        if u != u32::MAX && v != u32::MAX {
            out.add_bond(Bond::new(u, v, bond.order, BondType::Single))?;
        }
    }
    Ok((out, map[keep] as usize))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_iron_replaces_a_wall_carbon_and_keeps_three_bonds() {
        let set = RodPocketGenerator
            .build(&ParameterSet::new())
            .expect("the pocket builds");
        let topology = &set.block.topology;
        assert_eq!(topology.element(set.iron), Some(Element::IRON));
        let neighbours: Vec<Element> = topology
            .bonds()
            .filter_map(|bond| {
                if bond.u as usize == set.iron {
                    topology.element(bond.v as usize)
                } else if bond.v as usize == set.iron {
                    topology.element(bond.u as usize)
                } else {
                    None
                }
            })
            .collect();
        assert_eq!(neighbours, vec![Element::CARBON; 3]);
        assert!(set.iron_z_m.abs() < 4.0e-10, "{}", set.iron_z_m);
        assert!(set.mouth_z_m - set.iron_z_m > 3.0e-10);
    }

    #[test]
    fn the_rod_clears_the_bore_along_the_stroke() {
        let set = RodPocketGenerator
            .build(&ParameterSet::new())
            .expect("the pocket builds");
        let block = set.block.topology.positions_m().to_vec();
        let rod = set.rod.topology.positions_m().to_vec();
        let mut worst = f64::INFINITY;
        for step in 0..=10 {
            let face_m = set.iron_z_m - 4.0e-10 + 1.0e-10 * step as f64;
            let shift = face_m - set.rod_face_offset_m;
            for r in rod.chunks(3) {
                for b in block.chunks(3) {
                    let d = ((r[0] - b[0]).powi(2)
                        + (r[1] - b[1]).powi(2)
                        + (r[2] + shift - b[2]).powi(2))
                    .sqrt();
                    worst = worst.min(d);
                }
            }
        }
        // Hydrogen to hydrogen, the sense audit overlap limit is 1.6e-10 m.
        assert!(worst > 1.6e-10, "rod to block {worst:e} m");
    }
}
