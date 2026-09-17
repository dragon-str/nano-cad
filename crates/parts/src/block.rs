//! Block part generators: plate, beam, and bracket.
//!
//! Each generator cuts a solid region from the diamond cubic lattice, bonds
//! every first-shell carbon pair, and caps each surface carbon with hydrogen.
//! The three generators share one fill-and-cap pass, [`part_from_solid`], so a
//! filled block is built the same way in each case. All lengths are SI metres.

use nanocad_model::{Atom, Element, Part, Topology};
use nanocad_units::Unit;

use crate::diamond_solid;
use crate::error::PartError;
use crate::generator::PartGenerator;
use crate::lattice_fill::fill_solid;
use crate::parameter::{ParameterSet, ParameterSpec};
use crate::shape::{Box3, Cylinder, Difference, Solid, Union};

/// Fills a solid with the capped diamond lattice and names the part.
///
/// The pass runs once for every block generator: enumerate the lattice sites
/// inside `solid`, add one carbon at each site, bond first-shell neighbours,
/// and cap each surface carbon with hydrogen. The result carries the
/// `"diamond"` material label.
fn part_from_solid(name: &str, solid: &dyn Solid) -> Result<Part, PartError> {
    let positions_m = fill_solid(solid);
    let mut topology = Topology::new();
    let mut indices = Vec::with_capacity(positions_m.len());
    for position_m in positions_m {
        indices.push(topology.add_atom(Atom::new(Element::CARBON, position_m, 0.0, "C")));
    }
    diamond_solid::bond_and_cap(&mut topology, &indices)?;
    Ok(Part::new(name, topology).with_material("diamond"))
}

static PLATE_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new(
        "size_x_m",
        Some(Unit::Metre),
        4.0e-9,
        1.5e-9,
        5.0e-8,
        false,
        "plate size along x in metres",
    ),
    ParameterSpec::new(
        "size_y_m",
        Some(Unit::Metre),
        4.0e-9,
        1.5e-9,
        5.0e-8,
        false,
        "plate size along y in metres",
    ),
    ParameterSpec::new(
        "thickness_m",
        Some(Unit::Metre),
        1.5e-9,
        7.0e-10,
        2.0e-8,
        false,
        "plate thickness along z in metres",
    ),
    ParameterSpec::new(
        "hole_radius_m",
        Some(Unit::Metre),
        0.0,
        0.0,
        5.0e-9,
        false,
        "centred through-hole radius in metres; zero means no hole",
    ),
];

static BEAM_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new(
        "length_m",
        Some(Unit::Metre),
        6.0e-9,
        1.5e-9,
        1.0e-7,
        false,
        "beam length along x in metres",
    ),
    ParameterSpec::new(
        "width_m",
        Some(Unit::Metre),
        1.5e-9,
        7.0e-10,
        2.0e-8,
        false,
        "beam width along y in metres",
    ),
    ParameterSpec::new(
        "height_m",
        Some(Unit::Metre),
        1.5e-9,
        7.0e-10,
        2.0e-8,
        false,
        "beam height along z in metres",
    ),
];

static BRACKET_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new(
        "leg_a_m",
        Some(Unit::Metre),
        3.0e-9,
        1.5e-9,
        5.0e-8,
        false,
        "length of the leg along +x in metres",
    ),
    ParameterSpec::new(
        "leg_b_m",
        Some(Unit::Metre),
        3.0e-9,
        1.5e-9,
        5.0e-8,
        false,
        "length of the leg along +y in metres",
    ),
    ParameterSpec::new(
        "thickness_m",
        Some(Unit::Metre),
        1.5e-9,
        7.0e-10,
        1.0e-8,
        false,
        "leg thickness in z and leg depth in the xy plane, in metres",
    ),
];

/// Builds a centred rectangular plate, with an optional centred through-hole.
///
/// The plate is a [`Box3`] of `size_x_m` by `size_y_m` by `thickness_m`. When
/// `hole_radius_m` is above zero, the generator removes a centred [`Cylinder`]
/// of that radius through the full thickness.
#[derive(Clone, Copy, Debug, Default)]
pub struct PlateGenerator;

impl PartGenerator for PlateGenerator {
    fn id(&self) -> &'static str {
        "plate"
    }

    fn name(&self) -> &'static str {
        "plate"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        PLATE_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        let resolved = self.resolve(parameters)?;
        let size_x_m = resolved.require("size_x_m")?;
        let size_y_m = resolved.require("size_y_m")?;
        let thickness_m = resolved.require("thickness_m")?;
        let hole_radius_m = resolved.require("hole_radius_m")?;

        let plate = Box3::from_center_m([0.0, 0.0, 0.0], [size_x_m, size_y_m, thickness_m]);
        if hole_radius_m > 0.0 {
            let hole = Cylinder {
                center_m: [0.0, 0.0, 0.0],
                radius_m: hole_radius_m,
                height_m: 4.0 * thickness_m,
            };
            let solid = Difference {
                outer: Box::new(plate),
                inner: Box::new(hole),
            };
            part_from_solid("plate", &solid)
        } else {
            part_from_solid("plate", &plate)
        }
    }
}

/// Builds a centred rectangular beam.
///
/// The beam is a [`Box3`] with `length_m` along x, `width_m` along y, and
/// `height_m` along z. Its centre is the origin.
#[derive(Clone, Copy, Debug, Default)]
pub struct BeamGenerator;

impl PartGenerator for BeamGenerator {
    fn id(&self) -> &'static str {
        "beam"
    }

    fn name(&self) -> &'static str {
        "beam"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        BEAM_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        let resolved = self.resolve(parameters)?;
        let length_m = resolved.require("length_m")?;
        let width_m = resolved.require("width_m")?;
        let height_m = resolved.require("height_m")?;

        let beam = Box3::from_center_m([0.0, 0.0, 0.0], [length_m, width_m, height_m]);
        part_from_solid("beam", &beam)
    }
}

/// Builds an L-shaped bracket from two legs.
///
/// The legs share the corner at the origin of the xy plane. Leg A runs along
/// +x and leg B runs along +y. Each leg is `thickness_m` deep in the plane and
/// `thickness_m` thick in z. The two legs form a [`Union`], so the L is one
/// solid.
#[derive(Clone, Copy, Debug, Default)]
pub struct BracketGenerator;

impl PartGenerator for BracketGenerator {
    fn id(&self) -> &'static str {
        "bracket"
    }

    fn name(&self) -> &'static str {
        "bracket"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        BRACKET_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        let resolved = self.resolve(parameters)?;
        let leg_a_m = resolved.require("leg_a_m")?;
        let leg_b_m = resolved.require("leg_b_m")?;
        let thickness_m = resolved.require("thickness_m")?;

        let half_m = 0.5 * thickness_m;
        let leg_a = Box3 {
            min_m: [0.0, 0.0, -half_m],
            max_m: [leg_a_m, thickness_m, half_m],
        };
        let leg_b = Box3 {
            min_m: [0.0, 0.0, -half_m],
            max_m: [thickness_m, leg_b_m, half_m],
        };
        let bracket = Union {
            a: Box::new(leg_a),
            b: Box::new(leg_b),
        };
        part_from_solid("bracket", &bracket)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn span_m(part: &Part, axis: usize) -> (f64, f64) {
        let mut low = f64::INFINITY;
        let mut high = f64::NEG_INFINITY;
        for atom in part.topology.atoms() {
            low = low.min(atom.position_m[axis]);
            high = high.max(atom.position_m[axis]);
        }
        (low, high)
    }

    #[test]
    fn the_plate_is_a_box() {
        let plate = PlateGenerator
            .generate(
                &ParameterSet::new()
                    .with("size_x_m", 1.5e-9)
                    .with("size_y_m", 1.5e-9)
                    .with("thickness_m", 7.0e-10),
            )
            .expect("valid plate");
        let (min_x_m, max_x_m) = span_m(&plate, 0);
        let (min_y_m, max_y_m) = span_m(&plate, 1);
        assert!(((max_x_m - min_x_m) - 1.5e-9).abs() <= 1.0e-9);
        assert!(((max_y_m - min_y_m) - 1.5e-9).abs() <= 1.0e-9);
    }

    #[test]
    fn a_plate_hole_removes_the_centre() {
        let hole_radius_m = 7.0e-10;
        let plate = PlateGenerator
            .generate(
                &ParameterSet::new()
                    .with("size_x_m", 1.5e-9)
                    .with("size_y_m", 1.5e-9)
                    .with("thickness_m", 7.0e-10)
                    .with("hole_radius_m", hole_radius_m),
            )
            .expect("valid plate");
        let inner_m = hole_radius_m * (1.0 - 1.0e-6);
        for atom in plate.topology.atoms() {
            if atom.element != Element::CARBON {
                continue;
            }
            let radius_m = (atom.position_m[0] * atom.position_m[0]
                + atom.position_m[1] * atom.position_m[1])
                .sqrt();
            assert!(
                radius_m >= inner_m,
                "a carbon sits inside the hole at radius {radius_m} m"
            );
        }
    }

    #[test]
    fn a_longer_beam_has_more_atoms() {
        let cross_m = 7.0e-10;
        let short = BeamGenerator
            .generate(
                &ParameterSet::new()
                    .with("length_m", 6.0e-9)
                    .with("width_m", cross_m)
                    .with("height_m", cross_m),
            )
            .expect("valid beam");
        let long = BeamGenerator
            .generate(
                &ParameterSet::new()
                    .with("length_m", 1.2e-8)
                    .with("width_m", cross_m)
                    .with("height_m", cross_m),
            )
            .expect("valid beam");
        assert!(short.atom_count() > 0);
        assert!(long.atom_count() > short.atom_count());
    }

    #[test]
    fn the_bracket_is_an_l() {
        let thickness_m = 7.0e-10;
        let limit_m = 2.0 * thickness_m;
        let bracket = BracketGenerator
            .generate(
                &ParameterSet::new()
                    .with("leg_a_m", 3.0e-9)
                    .with("leg_b_m", 3.0e-9)
                    .with("thickness_m", thickness_m),
            )
            .expect("valid bracket");
        let mut beyond_x = false;
        let mut beyond_y = false;
        for atom in bracket.topology.atoms() {
            let x_m = atom.position_m[0];
            let y_m = atom.position_m[1];
            if x_m > limit_m && y_m > limit_m {
                panic!("the corner region holds an atom at ({x_m}, {y_m}) m");
            }
            beyond_x |= x_m > limit_m;
            beyond_y |= y_m > limit_m;
        }
        assert!(beyond_x, "no atom lies beyond the leg A limit");
        assert!(beyond_y, "no atom lies beyond the leg B limit");
    }

    #[test]
    fn every_carbon_keeps_a_valid_valence() {
        let plate = PlateGenerator
            .generate_with_defaults()
            .expect("valid plate");
        let mut degree = vec![0_u32; plate.atom_count()];
        for bond in plate.topology.bonds() {
            degree[bond.u as usize] += 1;
            degree[bond.v as usize] += 1;
        }
        for (index, atom) in plate.topology.atoms().enumerate() {
            if atom.element == Element::CARBON {
                let valence = degree[index];
                assert!(
                    (1..=4).contains(&valence),
                    "a carbon at index {index} has valence {valence}"
                );
            }
        }
    }
}
