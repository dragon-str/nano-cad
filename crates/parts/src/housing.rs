//! The sorting rotor housing.
//!
//! The housing holds the sorting rotor. It is an annular slab with a circular
//! chamber for the disk and two radial channels. The inlet channel lets the
//! outside solution reach the rim. The outlet channel leaves the chamber where
//! the pockets deliver their bound molecules.
//!
//! This generator builds the frame only. It does not build a binding site, and
//! it does not model a solvent or a solution.
//!
//! The housing is cut from the diamond cubic lattice, so every bond is the
//! diamond first-shell length and every carbon keeps a valence of four. All
//! lengths are SI metres and all angles are radians.

use nanocad_model::{Atom, Element, Part, Topology};
use nanocad_units::Unit;

use crate::diamond_solid;
use crate::error::PartError;
use crate::generator::PartGenerator;
use crate::lattice_fill::fill_solid;
use crate::parameter::{ParameterSet, ParameterSpec};
use crate::port::{Dof, Port};
use crate::shape::{Box3, Cylinder, Difference, Placed, Solid, Union};

/// The parameter specs of [`RotorHousingGenerator`], in a stable order.
static ROTOR_HOUSING_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new(
        "chamber_radius_m",
        Some(Unit::Metre),
        7.2e-9,
        2.0e-9,
        4.0e-8,
        false,
        "radius of the rotor chamber in metres",
    ),
    ParameterSpec::new(
        "wall_m",
        Some(Unit::Metre),
        1.5e-9,
        7.0e-10,
        1.0e-8,
        false,
        "radial wall thickness around the chamber in metres",
    ),
    ParameterSpec::new(
        "thickness_m",
        Some(Unit::Metre),
        2.4e-9,
        7.0e-10,
        2.0e-8,
        false,
        "axial thickness of the housing in metres",
    ),
    ParameterSpec::new(
        "channel_width_m",
        Some(Unit::Metre),
        2.0e-9,
        5.0e-10,
        1.0e-8,
        false,
        "width of each radial channel in metres",
    ),
    ParameterSpec::new(
        "inlet_angle_rad",
        None,
        0.0,
        -3.2,
        3.2,
        false,
        "angle of the inlet channel in radians",
    ),
    ParameterSpec::new(
        "outlet_angle_rad",
        None,
        std::f64::consts::PI,
        -3.2,
        3.2,
        false,
        "angle of the outlet channel in radians",
    ),
];

/// Generates a sorting rotor housing.
///
/// The part is the annulus between the chamber radius and the outside radius,
/// less one rectangular channel at each of the two stated angles. A channel runs
/// from inside the chamber wall to past the outside radius, so each one opens
/// the chamber to the outside.
///
/// The chamber is concentric about `z`, so the chamber port lies on that axis.
#[derive(Clone, Copy, Debug, Default)]
pub struct RotorHousingGenerator;

impl RotorHousingGenerator {
    /// Returns the chamber port.
    ///
    /// The port sits at the origin, its axis is `+z`, and it holds the housing
    /// rigidly. The chamber and the rotor share this axis.
    pub fn chamber_port(&self) -> Port {
        Port::named("chamber")
    }

    /// Returns the port of one channel mouth.
    ///
    /// The port sits on the chamber wall at `angle_rad`, and its axis points
    /// outward along the channel. It holds the housing rigidly.
    pub fn channel_port(&self, name: &str, angle_rad: f64, chamber_radius_m: f64) -> Port {
        let (sin_rad, cos_rad) = angle_rad.sin_cos();
        let mut port = Port::named(name);
        port.origin_m = [chamber_radius_m * cos_rad, chamber_radius_m * sin_rad, 0.0];
        port.axis_m = [cos_rad, sin_rad, 0.0];
        port.dof = Dof::Fixed;
        port
    }
}

impl PartGenerator for RotorHousingGenerator {
    fn id(&self) -> &'static str {
        "rotor_housing"
    }

    fn name(&self) -> &'static str {
        "rotor housing"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        ROTOR_HOUSING_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        let resolved = self.resolve(parameters)?;
        let chamber_radius_m = resolved.require("chamber_radius_m")?;
        let wall_m = resolved.require("wall_m")?;
        let thickness_m = resolved.require("thickness_m")?;
        let channel_width_m = resolved.require("channel_width_m")?;
        let inlet_angle_rad = resolved.require("inlet_angle_rad")?;
        let outlet_angle_rad = resolved.require("outlet_angle_rad")?;

        if wall_m <= 0.0 {
            return Err(PartError::InvalidGeometry(
                "the housing wall must be thicker than zero".to_string(),
            ));
        }
        if channel_width_m >= 2.0 * std::f64::consts::PI * chamber_radius_m {
            return Err(PartError::InvalidGeometry(format!(
                "the channel width {channel_width_m} m is wider than the chamber circumference, \
                 so the wall would vanish"
            )));
        }

        let outer_radius_m = chamber_radius_m + wall_m;
        let center_m = [0.0, 0.0, 0.0];
        let ring = Difference {
            outer: Box::new(Cylinder {
                center_m,
                radius_m: outer_radius_m,
                height_m: thickness_m,
            }),
            inner: Box::new(Cylinder {
                center_m,
                radius_m: chamber_radius_m,
                height_m: thickness_m,
            }),
        };

        // The channel spans the full wall and reaches past the outside radius,
        // so it cuts a clean slot instead of a trapped pocket.
        let reach_m = outer_radius_m + wall_m;
        let channel = |angle_rad: f64| -> Box<dyn Solid> {
            Box::new(Placed {
                inner: Box::new(Box3 {
                    min_m: [-reach_m, -0.5 * channel_width_m, -0.5 * thickness_m],
                    max_m: [reach_m, 0.5 * channel_width_m, 0.5 * thickness_m],
                }),
                rotation_rad: angle_rad,
                offset_m: center_m,
            })
        };
        let channels = Union {
            a: channel(inlet_angle_rad),
            b: channel(outlet_angle_rad),
        };

        let solid = Difference {
            outer: Box::new(ring),
            inner: Box::new(channels),
        };
        fill_part(&format!("rotor-housing-{chamber_radius_m:e}"), &solid)
    }
}

/// Fills a solid region with the diamond lattice and caps every surface atom.
///
/// The region must have a finite bounding box. The function returns an error
/// when the region cut no lattice site, so an empty body cannot reach a caller.
fn fill_part(name: &str, solid: &dyn Solid) -> Result<Part, PartError> {
    let sites_m = fill_solid(solid);
    let mut topology = Topology::new();
    let mut carbons = Vec::with_capacity(sites_m.len());
    for site_m in sites_m {
        carbons.push(topology.add_atom(Atom::new(Element::CARBON, site_m, 0.0, "C")));
    }
    if carbons.is_empty() {
        return Err(PartError::InvalidGeometry(format!(
            "{name} cut no atom from the diamond lattice"
        )));
    }
    diamond_solid::bond_and_cap(&mut topology, &carbons)?;
    Ok(Part::new(name, topology).with_material("diamond"))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn radial_xy_m(position_m: [f64; 3]) -> f64 {
        (position_m[0] * position_m[0] + position_m[1] * position_m[1]).sqrt()
    }

    #[test]
    fn the_housing_has_a_chamber() {
        let chamber_radius_m = 7.2e-9;
        let part = RotorHousingGenerator
            .generate_with_defaults()
            .expect("housing");
        let inside = part
            .topology
            .atoms()
            .filter(|atom| radial_xy_m(atom.position_m) < 0.5 * chamber_radius_m);
        assert_eq!(inside.count(), 0, "the chamber holds a carbon");
    }

    #[test]
    fn the_housing_has_two_channels() {
        let part = RotorHousingGenerator
            .generate_with_defaults()
            .expect("housing");
        let probe_m = 8.0e-9;
        for angle_rad in [0.0_f64, std::f64::consts::PI] {
            let (sin_rad, cos_rad) = angle_rad.sin_cos();
            let center_m = [probe_m * cos_rad, probe_m * sin_rad, 0.0];
            let half_m = 0.25e-9;
            let inside = part.topology.atoms().filter(|atom| {
                let dx = atom.position_m[0] - center_m[0];
                let dy = atom.position_m[1] - center_m[1];
                (dx * dx + dy * dy).sqrt() < half_m
            });
            assert_eq!(
                inside.count(),
                0,
                "the channel at {angle_rad} rad holds a carbon"
            );
        }
    }

    #[test]
    fn the_housing_wall_is_solid() {
        let part = RotorHousingGenerator
            .generate_with_defaults()
            .expect("housing");
        let count = part
            .topology
            .atoms()
            .filter(|atom| {
                let radial_m = radial_xy_m(atom.position_m);
                radial_m > 7.2e-9 && radial_m < 8.7e-9
            })
            .count();
        assert!(count > 0, "the housing wall holds no atom");
    }

    #[test]
    fn the_housing_outside_is_open() {
        let part = RotorHousingGenerator
            .generate_with_defaults()
            .expect("housing");
        let farthest_m = part
            .topology
            .atoms()
            .map(|atom| radial_xy_m(atom.position_m))
            .fold(0.0_f64, f64::max);
        assert!(
            farthest_m <= 8.7e-9 + 1.0e-10,
            "an atom sits at {farthest_m:e} m, past the outside radius"
        );
    }

    #[test]
    fn the_housing_ports_are_stated() {
        let chamber = RotorHousingGenerator.chamber_port();
        assert_eq!(chamber.name, "chamber");
        assert_eq!(chamber.dof, Dof::Fixed);
        assert_eq!(chamber.origin_m, [0.0, 0.0, 0.0]);
        assert_eq!(chamber.axis_m, [0.0, 0.0, 1.0]);

        let inlet = RotorHousingGenerator.channel_port("inlet", 0.0, 7.2e-9);
        assert_eq!(inlet.name, "inlet");
        assert_eq!(inlet.dof, Dof::Fixed);
        assert!((inlet.origin_m[0] - 7.2e-9).abs() < 1.0e-18);
        assert!((inlet.axis_m[0] - 1.0).abs() < 1.0e-15);

        let outlet = RotorHousingGenerator.channel_port("outlet", std::f64::consts::PI, 7.2e-9);
        assert!((outlet.origin_m[0] + 7.2e-9).abs() < 1.0e-18);
        assert!((outlet.axis_m[0] + 1.0).abs() < 1.0e-15);
    }

    #[test]
    fn a_zero_wall_is_refused() {
        let parameters = ParameterSet::new().with("wall_m", 0.0);
        assert!(RotorHousingGenerator.generate(&parameters).is_err());
    }

    #[test]
    fn every_carbon_keeps_a_valid_valence() {
        let part = RotorHousingGenerator
            .generate_with_defaults()
            .expect("housing");
        let mut degree = vec![0usize; part.atom_count()];
        for bond in part.topology.bonds() {
            degree[bond.u as usize] += 1;
            degree[bond.v as usize] += 1;
        }
        for (index, &count) in degree.iter().enumerate() {
            if part.topology.element(index) == Some(Element::CARBON) {
                assert!(
                    (1..=4).contains(&count),
                    "carbon {index} keeps a degree of {count}, outside 1 to 4"
                );
            }
        }
    }
}
