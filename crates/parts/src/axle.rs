//! Hexagonal axle and plain shaft generators.
//!
//! Both parts are solid diamond cubic carbon, filled by [`fill_solid`] and
//! capped with hydrogen by [`diamond_solid::bond_and_cap`]. The axis of each
//! part is z, and each part is centred on the origin. Every length is in SI
//! metres.
//!
//! The hex axle is a regular hexagonal prism. Its `radius_m` is the
//! circumradius, and the caller gives the across-flats distance. The
//! across-flats distance of a regular hexagon is twice the apothem, and the
//! apothem is `radius * cos(30 deg)`, so the circumradius is
//! `across_flats / (2 * cos(30 deg)) = across_flats / sqrt(3)`.
//!
//! The chamfer shortens the filled length. It is not a conical cut. A
//! `Cylinder` is not conical, and a `Placed` box is a blunt removal. The
//! honest, closable model is an axial shortening: `chamfer_m` removes
//! `2 * chamfer_m` from the filled length, in total. When the shortened height
//! is not positive, the generator uses the full length and treats the chamfer
//! as zero.

use nanocad_model::{Atom, Element, Part, Topology};
use nanocad_units::Unit;

use crate::diamond_solid;
use crate::error::PartError;
use crate::generator::PartGenerator;
use crate::lattice_fill::fill_solid;
use crate::parameter::{ParameterSet, ParameterSpec};
use crate::shape::{Cylinder, HexPrism};

static HEX_AXLE_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new(
        "across_flats_m",
        Some(Unit::Metre),
        3.0e-9,
        8.0e-10,
        2.0e-8,
        false,
        "across-flats distance of the hexagonal cross section",
    ),
    ParameterSpec::new(
        "length_m",
        Some(Unit::Metre),
        1.2e-8,
        1.5e-9,
        1.0e-7,
        false,
        "axial length of the axle",
    ),
    ParameterSpec::new(
        "chamfer_m",
        Some(Unit::Metre),
        2.5e-10,
        0.0,
        5.0e-10,
        false,
        "axial length removed at each end",
    ),
];

static PLAIN_SHAFT_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new(
        "radius_m",
        Some(Unit::Metre),
        1.0e-9,
        4.0e-10,
        1.0e-8,
        false,
        "radius of the circular cross section",
    ),
    ParameterSpec::new(
        "length_m",
        Some(Unit::Metre),
        1.0e-8,
        1.5e-9,
        1.0e-7,
        false,
        "axial length of the shaft",
    ),
];

/// Builds a hexagonal axle.
///
/// The hexagon is centred on the origin with its axis along z. The parameter
/// `across_flats_m` fixes the hexagon, and `length_m` fixes the axial length.
/// The parameter `chamfer_m` shortens the filled length. The result is a solid
/// diamond part with hydrogen-capped surface carbons.
#[derive(Clone, Copy, Debug, Default)]
pub struct HexAxleGenerator;

impl PartGenerator for HexAxleGenerator {
    fn id(&self) -> &'static str {
        "hex_axle"
    }

    fn name(&self) -> &'static str {
        "hex axle"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        HEX_AXLE_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        let resolved = self.resolve(parameters)?;
        let across_flats_m = resolved.require("across_flats_m")?;
        let length_m = resolved.require("length_m")?;
        let chamfer_m = resolved.require("chamfer_m")?;

        let radius_m = across_flats_m / 3.0_f64.sqrt();
        let height_m = chamfered_height_m(length_m, chamfer_m);
        let solid = HexPrism {
            center_m: [0.0; 3],
            radius_m,
            height_m,
            rotation_rad: 0.0,
        };

        let topology = filled_carbon_topology(&solid)?;
        let name = "hex-axle".to_owned();
        Ok(Part::new(name, topology).with_material("diamond"))
    }
}

/// Builds a plain cylindrical shaft.
///
/// The cylinder is centred on the origin with its axis along z. The parameter
/// `radius_m` fixes the circular cross section, and `length_m` fixes the axial
/// length. The result is a solid diamond part with hydrogen-capped surface
/// carbons.
#[derive(Clone, Copy, Debug, Default)]
pub struct PlainShaftGenerator;

impl PartGenerator for PlainShaftGenerator {
    fn id(&self) -> &'static str {
        "plain_shaft"
    }

    fn name(&self) -> &'static str {
        "plain shaft"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        PLAIN_SHAFT_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        let resolved = self.resolve(parameters)?;
        let radius_m = resolved.require("radius_m")?;
        let length_m = resolved.require("length_m")?;

        let solid = Cylinder {
            center_m: [0.0; 3],
            radius_m,
            height_m: length_m,
        };

        let topology = filled_carbon_topology(&solid)?;
        let name = "plain-shaft".to_owned();
        Ok(Part::new(name, topology).with_material("diamond"))
    }
}

/// Returns the filled height after the axial chamfer.
///
/// The chamfer removes `2 * chamfer_m` in total. When the result is not
/// positive, the function returns the full length, so the chamfer is zero.
fn chamfered_height_m(length_m: f64, chamfer_m: f64) -> f64 {
    let height_m = length_m - 2.0 * chamfer_m;
    if height_m > 0.0 {
        height_m
    } else {
        length_m
    }
}

/// Fills a solid with diamond carbon and caps the surface with hydrogen.
///
/// Every filled site becomes a carbon atom. The bond and cap pass gives each
/// carbon a valence of four, so no carbon is left with degree zero.
fn filled_carbon_topology(solid: &dyn crate::shape::Solid) -> Result<Topology, PartError> {
    let mut topology = Topology::new();
    let positions_m = fill_solid(solid);
    let mut indices = Vec::with_capacity(positions_m.len());
    for position_m in positions_m {
        indices.push(topology.add_atom(Atom::new(Element::CARBON, position_m, 0.0, "C")));
    }
    diamond_solid::bond_and_cap(&mut topology, &indices)?;
    Ok(topology)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Returns the radius in the xy plane for an atom.
    fn radius_m(position_m: [f64; 3]) -> f64 {
        (position_m[0] * position_m[0] + position_m[1] * position_m[1]).sqrt()
    }

    /// Returns the smallest radius over the filled carbon atoms.
    fn z_span_m(part: &Part) -> f64 {
        let mut min_m = f64::INFINITY;
        let mut max_m = f64::NEG_INFINITY;
        for atom in part.topology.atoms() {
            min_m = min_m.min(atom.position_m[2]);
            max_m = max_m.max(atom.position_m[2]);
        }
        max_m - min_m
    }

    /// Returns the bond degree of every carbon atom.
    fn carbon_degrees(part: &Part) -> Vec<u32> {
        let topology = &part.topology;
        let mut degrees = vec![0_u32; topology.atom_count()];
        for index in 0..topology.bond_count() {
            if let (Some(u), Some(v)) = (topology.bond_u(index), topology.bond_v(index)) {
                degrees[u as usize] += 1;
                degrees[v as usize] += 1;
            }
        }
        (0..topology.atom_count())
            .filter(|&index| topology.element(index) == Some(Element::CARBON))
            .map(|index| degrees[index])
            .collect()
    }

    #[test]
    fn the_hex_axle_has_six_flats() {
        let part = HexAxleGenerator.generate_with_defaults().expect("generate");
        assert!(!part.topology.is_empty(), "the axle filled with no atom");

        let circumradius_m = 3.0e-9 / 3.0_f64.sqrt();
        let max_radius_m = part
            .topology
            .atoms()
            .filter(|atom| atom.element == Element::CARBON)
            .map(|atom| radius_m(atom.position_m))
            .fold(0.0_f64, f64::max);
        assert!(max_radius_m > 0.0, "no carbon sits off the axis");
        assert!(
            max_radius_m <= circumradius_m * (1.0 + 1.0e-6),
            "a carbon sits at {max_radius_m} m, beyond the circumradius {circumradius_m} m"
        );
    }

    #[test]
    fn the_hex_axle_respects_its_length() {
        let part = HexAxleGenerator.generate_with_defaults().expect("generate");
        let length_m = 1.2e-8;
        let chamfer_m = 2.5e-10;
        let span_m = z_span_m(&part);
        assert!(span_m <= length_m, "the span {span_m} m exceeds the length");
        assert!(
            span_m >= length_m - 2.0 * chamfer_m - 1.0e-9,
            "the span {span_m} m is shorter than the chamfered length"
        );
    }

    #[test]
    fn a_bigger_across_flats_makes_more_atoms() {
        let small = HexAxleGenerator
            .generate(&ParameterSet::new().with("length_m", 1.5e-9))
            .expect("generate");
        let big = HexAxleGenerator
            .generate(
                &ParameterSet::new()
                    .with("length_m", 1.5e-9)
                    .with("across_flats_m", 6.0e-9),
            )
            .expect("generate");
        assert!(small.atom_count() > 0, "the small axle is empty");
        assert!(
            big.atom_count() > small.atom_count(),
            "the widened axle has {} atoms, not more than {}",
            big.atom_count(),
            small.atom_count()
        );
    }

    #[test]
    fn the_shaft_respects_its_radius() {
        let part = PlainShaftGenerator
            .generate_with_defaults()
            .expect("generate");
        let radius_limit_m = 1.0e-9 * (1.0 + 1.0e-6);
        for atom in part.topology.atoms() {
            if atom.element != Element::CARBON {
                continue;
            }
            let atom_radius_m = radius_m(atom.position_m);
            assert!(
                atom_radius_m <= radius_limit_m,
                "a carbon sits at {atom_radius_m} m, beyond {radius_limit_m} m"
            );
        }
    }

    #[test]
    fn both_parts_cap_their_surface_carbons() {
        let hex = HexAxleGenerator.generate_with_defaults().expect("generate");
        let shaft = PlainShaftGenerator
            .generate_with_defaults()
            .expect("generate");
        for part in [hex, shaft] {
            let degrees = carbon_degrees(&part);
            assert!(!degrees.is_empty(), "the part holds no carbon");
            for degree in degrees {
                assert!(
                    (1..=4).contains(&degree),
                    "a carbon has degree {degree}, not 1 to 4"
                );
            }
        }
    }
}
