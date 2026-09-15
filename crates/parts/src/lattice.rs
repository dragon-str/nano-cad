use nanocad_model::{Atom, Bond, BondType, Element, Part, Topology};

use crate::error::PartError;
use crate::generator::PartGenerator;
use crate::parameter::{ParameterSet, ParameterSpec};

/// The relative tolerance for a first-shell bond. Distances within this
/// fraction of the expected bond length become bonds.
const BOND_TOLERANCE_RELATIVE: f64 = 1.0e-3;

/// The diamond conventional cubic lattice constant at 300 K.
///
/// Source: CRC Handbook of Chemistry and Physics; N. W. Ashcroft and N. D.
/// Mermin, Solid State Physics (1976), chapter 4.
const DIAMOND_LATTICE_CONSTANT_M: f64 = 3.567e-10;

/// The graphite in-plane lattice constant `a` at 300 K.
///
/// Source: P. Trucano and R. Chen, "Structure of graphite by neutron
/// diffraction", Nature 258, 136 (1975).
const GRAPHITE_LATTICE_CONSTANT_M: f64 = 2.461e-10;

/// The graphite interlayer spacing `c/2` at 300 K.
///
/// Source: P. Trucano and R. Chen, Nature 258, 136 (1975).
const GRAPHITE_LAYER_SPACING_M: f64 = 3.354e-10;

/// The eight basis atoms of the diamond cubic cell, in fractional
/// coordinates. Two interpenetrating FCC lattices, offset by (1/4, 1/4, 1/4).
const DIAMOND_BASIS: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0],
    [0.0, 0.5, 0.5],
    [0.5, 0.0, 0.5],
    [0.5, 0.5, 0.0],
    [0.25, 0.25, 0.25],
    [0.25, 0.75, 0.75],
    [0.75, 0.25, 0.75],
    [0.75, 0.75, 0.25],
];

/// The two Cartesian basis offsets of the graphite honeycomb cell, in metres
/// along the in-plane `y` axis. Sublattice `A` sits at the lattice point.
fn graphite_basis_y_m() -> [f64; 2] {
    [0.0, GRAPHITE_LATTICE_CONSTANT_M / 3.0_f64.sqrt()]
}

static DIAMOND_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new("cells_x", None, 1.0, 1.0, 8.0, true, "unit cells along x"),
    ParameterSpec::new("cells_y", None, 1.0, 1.0, 8.0, true, "unit cells along y"),
    ParameterSpec::new("cells_z", None, 1.0, 1.0, 8.0, true, "unit cells along z"),
];

static GRAPHITE_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new("cells_a", None, 1.0, 1.0, 16.0, true, "unit cells along a1"),
    ParameterSpec::new("cells_b", None, 1.0, 1.0, 16.0, true, "unit cells along a2"),
    ParameterSpec::new("layers", None, 1.0, 1.0, 8.0, true, "graphene layers"),
];

/// Builds a diamond cubic lattice.
///
/// The requested size is a block of `cells_x` by `cells_y` by `cells_z`
/// conventional cubic cells. The block holds eight atoms per cell. Bonds join
/// carbon atoms at the first-shell distance `a * sqrt(3) / 4`.
#[derive(Clone, Copy, Debug, Default)]
pub struct DiamondGenerator;

impl PartGenerator for DiamondGenerator {
    fn id(&self) -> &'static str {
        "diamond"
    }

    fn name(&self) -> &'static str {
        "Diamond cubic lattice"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        DIAMOND_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        let resolved = self.resolve(parameters)?;
        let cells_x = resolved.require("cells_x")? as usize;
        let cells_y = resolved.require("cells_y")? as usize;
        let cells_z = resolved.require("cells_z")? as usize;

        let mut topology = Topology::new();
        for ix in 0..cells_x {
            for iy in 0..cells_y {
                for iz in 0..cells_z {
                    for basis in DIAMOND_BASIS {
                        let position_m = [
                            (ix as f64 + basis[0]) * DIAMOND_LATTICE_CONSTANT_M,
                            (iy as f64 + basis[1]) * DIAMOND_LATTICE_CONSTANT_M,
                            (iz as f64 + basis[2]) * DIAMOND_LATTICE_CONSTANT_M,
                        ];
                        topology.add_atom(Atom::new(Element::CARBON, position_m, 0.0, "C"));
                    }
                }
            }
        }

        let bond_length_m = DIAMOND_LATTICE_CONSTANT_M * 3.0_f64.sqrt() / 4.0;
        first_shell_bonds(&mut topology, bond_length_m)?;

        let name = format!("diamond-{cells_x}x{cells_y}x{cells_z}");
        Ok(Part::new(name, topology).with_material("diamond"))
    }
}

/// Builds an AB-stacked graphite (hexagonal) lattice.
///
/// The requested size is a block of `cells_a` by `cells_b` in-plane cells and
/// `layers` graphene layers. Each cell holds two atoms. Bonds join carbon atoms
/// at the in-plane first-shell distance `a / sqrt(3)`. The interlayer spacing
/// is larger, so the bonds stay in plane and each interior atom has three.
#[derive(Clone, Copy, Debug, Default)]
pub struct GraphiteGenerator;

impl PartGenerator for GraphiteGenerator {
    fn id(&self) -> &'static str {
        "graphite"
    }

    fn name(&self) -> &'static str {
        "Graphite hexagonal lattice"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        GRAPHITE_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        let resolved = self.resolve(parameters)?;
        let cells_a = resolved.require("cells_a")? as usize;
        let cells_b = resolved.require("cells_b")? as usize;
        let layers = resolved.require("layers")? as usize;

        let a_m = GRAPHITE_LATTICE_CONSTANT_M;
        let a1_m = [a_m, 0.0];
        let a2_m = [a_m * 0.5, a_m * 3.0_f64.sqrt() / 2.0];
        let basis_y_m = graphite_basis_y_m();

        let mut topology = Topology::new();
        for layer in 0..layers {
            let z_m = layer as f64 * GRAPHITE_LAYER_SPACING_M;
            let stacking_shift_m = if layer % 2 == 0 { 0.0 } else { basis_y_m[1] };
            for ia in 0..cells_a {
                for ib in 0..cells_b {
                    let origin_m = [
                        ia as f64 * a1_m[0] + ib as f64 * a2_m[0],
                        ia as f64 * a1_m[1] + ib as f64 * a2_m[1] + stacking_shift_m,
                    ];
                    for offset_y_m in graphite_basis_y_m() {
                        let position_m = [origin_m[0], origin_m[1] + offset_y_m, z_m];
                        topology.add_atom(Atom::new(Element::CARBON, position_m, 0.0, "C"));
                    }
                }
            }
        }

        let bond_length_m = a_m / 3.0_f64.sqrt();
        first_shell_bonds(&mut topology, bond_length_m)?;

        let name = format!("graphite-{cells_a}x{cells_b}x{layers}");
        Ok(Part::new(name, topology).with_material("graphite"))
    }
}

/// Adds one single bond for every atom pair at the first-shell distance.
///
/// The search visits each unordered pair once, so a bond can never duplicate.
fn first_shell_bonds(topology: &mut Topology, expected_m: f64) -> Result<(), PartError> {
    let low_m = expected_m * (1.0 - BOND_TOLERANCE_RELATIVE);
    let high_m = expected_m * (1.0 + BOND_TOLERANCE_RELATIVE);
    let atom_count = topology.atom_count();
    for u in 0..atom_count {
        for v in (u + 1)..atom_count {
            let position_u_m = match topology.position_m(u) {
                Some(position) => position,
                None => continue,
            };
            let position_v_m = match topology.position_m(v) {
                Some(position) => position,
                None => continue,
            };
            let distance_m = distance_m(position_u_m, position_v_m);
            if distance_m >= low_m && distance_m <= high_m {
                topology.add_bond(Bond::new(u as u32, v as u32, 1, BondType::Single))?;
            }
        }
    }
    Ok(())
}

fn distance_m(a_m: [f64; 3], b_m: [f64; 3]) -> f64 {
    let dx = a_m[0] - b_m[0];
    let dy = a_m[1] - b_m[1];
    let dz = a_m[2] - b_m[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use super::*;

    fn coordination(topology: &Topology) -> Vec<usize> {
        let mut counts = vec![0usize; topology.atom_count()];
        for bond in topology.bonds() {
            counts[bond.u as usize] += 1;
            counts[bond.v as usize] += 1;
        }
        counts
    }

    fn assert_bond_lengths(part: &Part, expected_m: f64) {
        assert!(part.bond_count() > 0, "expected at least one bond");
        for bond in part.topology.bonds() {
            let u_m = part.topology.position_m(bond.u as usize).expect("atom u");
            let v_m = part.topology.position_m(bond.v as usize).expect("atom v");
            let length_m = distance_m(u_m, v_m);
            let delta_m = (length_m - expected_m).abs();
            assert!(
                delta_m <= expected_m * BOND_TOLERANCE_RELATIVE,
                "bond length {length_m:e} differs from {expected_m:e}"
            );
        }
    }

    fn assert_no_duplicate_bonds(part: &Part) {
        let mut seen = BTreeSet::new();
        for bond in part.topology.bonds() {
            let key = if bond.u <= bond.v {
                (bond.u, bond.v)
            } else {
                (bond.v, bond.u)
            };
            assert!(seen.insert(key), "duplicate bond {key:?}");
        }
        assert_eq!(seen.len(), part.bond_count());
    }

    #[test]
    fn diamond_atom_count_is_eight_per_unit_cell() {
        let cells = ParameterSet::new()
            .with("cells_x", 2.0)
            .with("cells_y", 3.0)
            .with("cells_z", 1.0);
        let part = DiamondGenerator.generate(&cells).expect("generate");
        assert_eq!(part.atom_count(), 8 * 2 * 3);
        assert_eq!(part.material, "diamond");
        assert_eq!(part.name, "diamond-2x3x1");
    }

    #[test]
    fn diamond_bonds_match_the_first_shell_distance() {
        let cells = ParameterSet::new()
            .with("cells_x", 2.0)
            .with("cells_y", 2.0)
            .with("cells_z", 2.0);
        let part = DiamondGenerator.generate(&cells).expect("generate");
        let expected_m = DIAMOND_LATTICE_CONSTANT_M * 3.0_f64.sqrt() / 4.0;
        assert_bond_lengths(&part, expected_m);
    }

    #[test]
    fn diamond_has_no_duplicate_bonds() {
        let cells = ParameterSet::new()
            .with("cells_x", 2.0)
            .with("cells_y", 2.0)
            .with("cells_z", 2.0);
        let part = DiamondGenerator.generate(&cells).expect("generate");
        assert_no_duplicate_bonds(&part);
    }

    #[test]
    fn diamond_interior_atoms_have_four_neighbours() {
        let cells = ParameterSet::new()
            .with("cells_x", 2.0)
            .with("cells_y", 2.0)
            .with("cells_z", 2.0);
        let part = DiamondGenerator.generate(&cells).expect("generate");
        let counts = coordination(&part.topology);
        let interior = counts.iter().filter(|count| **count == 4).count();
        assert!(interior > 0, "expected an interior atom with four bonds");
        assert!(counts.iter().all(|count| *count <= 4));
    }

    #[test]
    fn invalid_diamond_size_is_an_error_not_a_panic() {
        let zero = ParameterSet::new().with("cells_x", 0.0);
        assert!(DiamondGenerator.generate(&zero).is_err());
        let fraction = ParameterSet::new().with("cells_x", 1.5);
        assert!(DiamondGenerator.generate(&fraction).is_err());
    }

    #[test]
    fn graphite_atom_count_is_two_per_cell_per_layer() {
        let cells = ParameterSet::new()
            .with("cells_a", 3.0)
            .with("cells_b", 2.0)
            .with("layers", 2.0);
        let part = GraphiteGenerator.generate(&cells).expect("generate");
        assert_eq!(part.atom_count(), 3 * 2 * 2 * 2);
        assert_eq!(part.material, "graphite");
        assert_eq!(part.name, "graphite-3x2x2");
    }

    #[test]
    fn graphite_bonds_match_the_first_shell_distance() {
        let cells = ParameterSet::new()
            .with("cells_a", 4.0)
            .with("cells_b", 4.0)
            .with("layers", 1.0);
        let part = GraphiteGenerator.generate(&cells).expect("generate");
        let expected_m = GRAPHITE_LATTICE_CONSTANT_M / 3.0_f64.sqrt();
        assert_bond_lengths(&part, expected_m);
    }

    #[test]
    fn graphite_has_no_duplicate_bonds() {
        let cells = ParameterSet::new()
            .with("cells_a", 4.0)
            .with("cells_b", 4.0)
            .with("layers", 1.0);
        let part = GraphiteGenerator.generate(&cells).expect("generate");
        assert_no_duplicate_bonds(&part);
    }

    #[test]
    fn graphite_in_plane_coordination_is_three() {
        let cells = ParameterSet::new()
            .with("cells_a", 4.0)
            .with("cells_b", 4.0)
            .with("layers", 1.0);
        let part = GraphiteGenerator.generate(&cells).expect("generate");
        let counts = coordination(&part.topology);
        assert!(counts.iter().all(|count| *count <= 3));
        let central_index = 2 * 4 * 2 + 2 * 2;
        assert_eq!(counts[central_index], 3);
    }

    #[test]
    fn invalid_graphite_size_is_an_error_not_a_panic() {
        let zero = ParameterSet::new().with("cells_a", 0.0);
        assert!(GraphiteGenerator.generate(&zero).is_err());
        let layers = ParameterSet::new().with("layers", 0.0);
        assert!(GraphiteGenerator.generate(&layers).is_err());
    }

    #[test]
    fn a_generator_reports_its_identity_and_defaults() {
        assert_eq!(DiamondGenerator.id(), "diamond");
        assert_eq!(DiamondGenerator.name(), "Diamond cubic lattice");
        assert_eq!(DiamondGenerator.parameters().len(), 3);
        assert_eq!(GraphiteGenerator.id(), "graphite");
        assert_eq!(GraphiteGenerator.parameters().len(), 3);
        assert!(GraphiteGenerator.generate_with_defaults().is_ok());
        assert_eq!(
            DIAMOND_PARAMETERS[0].unit, None,
            "cell counts are dimensionless"
        );
    }
}
