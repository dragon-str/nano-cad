//! Radial rolling bearing and plain bushing generators.
//!
//! A radial bearing is a concentric race pair with a stated count of cylindrical
//! rolling elements between the races. A bushing is a plain sleeve with a stated
//! running clearance to its shaft. Both parts are cut from the diamond cubic
//! lattice, so every bond is the diamond first-shell length and every carbon
//! keeps a valence of four.
//!
//! The generators follow [`PartGenerator`] and never panic. All lengths are SI
//! metres and all angles are radians.

use std::f64::consts::PI;

use nanocad_model::{Atom, Element, Part, Topology};
use nanocad_units::Unit;

use crate::diamond_solid;
use crate::error::PartError;
use crate::generator::PartGenerator;
use crate::lattice_fill::fill_solid;
use crate::parameter::{ParameterSet, ParameterSpec};
use crate::port::{Dof, Port};
use crate::shape::{Cylinder, Difference, Solid, Union};

/// The parameter specs of [`RadialBearingGenerator`], in a stable order.
static RADIAL_BEARING_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new(
        "inner_radius_m",
        Some(Unit::Metre),
        3.0e-9,
        1.0e-9,
        2.0e-8,
        false,
        "inner race bore radius in metres",
    ),
    ParameterSpec::new(
        "outer_radius_m",
        Some(Unit::Metre),
        5.0e-9,
        2.0e-9,
        4.0e-8,
        false,
        "outer race outside radius in metres",
    ),
    ParameterSpec::new(
        "height_m",
        Some(Unit::Metre),
        2.0e-9,
        7.0e-10,
        2.0e-8,
        false,
        "axial height of the bearing in metres",
    ),
    ParameterSpec::new(
        "roller_count",
        None,
        8.0,
        4.0,
        16.0,
        true,
        "number of cylindrical rolling elements",
    ),
    ParameterSpec::new(
        "radial_gap_m",
        Some(Unit::Metre),
        5.0e-10,
        1.0e-10,
        2.0e-9,
        false,
        "race-to-roller radial clearance in metres",
    ),
];

/// The parameter specs of [`BushingGenerator`], in a stable order.
static BUSHING_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new(
        "inner_radius_m",
        Some(Unit::Metre),
        3.0e-9,
        1.0e-9,
        2.0e-8,
        false,
        "sleeve bore radius in metres",
    ),
    ParameterSpec::new(
        "wall_m",
        Some(Unit::Metre),
        1.0e-9,
        5.0e-10,
        1.0e-8,
        false,
        "sleeve wall thickness in metres",
    ),
    ParameterSpec::new(
        "height_m",
        Some(Unit::Metre),
        2.0e-9,
        7.0e-10,
        2.0e-8,
        false,
        "sleeve axial length in metres",
    ),
    ParameterSpec::new(
        "clearance_m",
        Some(Unit::Metre),
        2.5e-10,
        0.0,
        1.0e-9,
        false,
        "running radial clearance to the shaft in metres",
    ),
];

/// Generates a radial rolling bearing.
///
/// The part is one diamondoid solid made from three regions. The inner race is
/// the annulus between the bore radius and the bore radius plus the radial gap.
/// The outer race is the annulus between the outside radius minus the radial gap
/// and the outside radius. The rolling elements are `roller_count` cylinders of
/// radius `radial_gap_m` on the pitch circle at the mean of the two race radii.
///
/// The two races and the rolling elements are concentric about `z`, so the inner
/// and outer ports share that axis.
#[derive(Clone, Copy, Debug, Default)]
pub struct RadialBearingGenerator;

impl RadialBearingGenerator {
    /// Returns the inner-ring port.
    ///
    /// The port sits at the origin, its axis is `+z`, and it allows one rotation
    /// about that axis. The gap is zero. The `inner_radius_m` argument is
    /// accepted for interface parity with the built geometry; the port position
    /// does not depend on it. The inner and outer ports share the `z` axis.
    pub fn inner_port(&self, inner_radius_m: f64) -> Port {
        let _ = inner_radius_m;
        let mut port = Port::named("inner");
        port.dof = Dof::Revolute;
        port
    }

    /// Returns the outer-ring port.
    ///
    /// The port sits at the origin, its axis is `+z`, and it holds both parts
    /// rigidly. The gap is zero. It shares the `z` axis with the inner port.
    pub fn outer_port(&self) -> Port {
        Port::named("outer")
    }
}

impl PartGenerator for RadialBearingGenerator {
    fn id(&self) -> &'static str {
        "radial_bearing"
    }

    fn name(&self) -> &'static str {
        "radial bearing"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        RADIAL_BEARING_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        let resolved = self.resolve(parameters)?;
        let inner_radius_m = resolved.require("inner_radius_m")?;
        let outer_radius_m = resolved.require("outer_radius_m")?;
        let height_m = resolved.require("height_m")?;
        let roller_count = resolved.require("roller_count")? as usize;
        let radial_gap_m = resolved.require("radial_gap_m")?;

        if outer_radius_m < inner_radius_m + 2.0 * radial_gap_m {
            return Err(PartError::InvalidGeometry(format!(
                "outer radius {outer_radius_m} m is below inner radius {inner_radius_m} m \
                 plus twice the radial gap {radial_gap_m} m, so the races would overlap"
            )));
        }

        let center_m = [0.0, 0.0, 0.0];
        let inner_race = Difference {
            outer: Box::new(Cylinder {
                center_m,
                radius_m: inner_radius_m + radial_gap_m,
                height_m,
            }),
            inner: Box::new(Cylinder {
                center_m,
                radius_m: inner_radius_m,
                height_m,
            }),
        };
        let outer_race = Difference {
            outer: Box::new(Cylinder {
                center_m,
                radius_m: outer_radius_m,
                height_m,
            }),
            inner: Box::new(Cylinder {
                center_m,
                radius_m: outer_radius_m - radial_gap_m,
                height_m,
            }),
        };

        let pitch_radius_m = 0.5 * (inner_radius_m + outer_radius_m);
        let mut solid: Box<dyn Solid> = Box::new(Union {
            a: Box::new(inner_race),
            b: Box::new(outer_race),
        });
        for index in 0..roller_count {
            let angle_rad = 2.0 * PI * index as f64 / roller_count as f64;
            let (sin_rad, cos_rad) = angle_rad.sin_cos();
            let roller = Cylinder {
                center_m: [pitch_radius_m * cos_rad, pitch_radius_m * sin_rad, 0.0],
                radius_m: radial_gap_m,
                height_m,
            };
            solid = Box::new(Union {
                a: solid,
                b: Box::new(roller),
            });
        }

        fill_part(&format!("radial-bearing-{roller_count}"), solid.as_ref())
    }
}

/// Generates a plain bushing.
///
/// The part is one diamondoid sleeve. It is the annulus between the bore radius
/// plus the wall and the bore radius plus the clearance. The clearance is the
/// stated running gap between the sleeve bore and its shaft.
#[derive(Clone, Copy, Debug, Default)]
pub struct BushingGenerator;

impl BushingGenerator {
    /// Returns the shaft port.
    ///
    /// The port sits on the axis at the origin, its axis is `+z`, and it allows
    /// one rotation about that axis. The gap equals the running clearance
    /// `clearance_m`, so a shaft mates at the same clearance the sleeve is cut
    /// with.
    pub fn shaft_port(&self, clearance_m: f64) -> Port {
        let mut port = Port::named("shaft");
        port.dof = Dof::Revolute;
        port.gap_m = clearance_m;
        port
    }
}

impl PartGenerator for BushingGenerator {
    fn id(&self) -> &'static str {
        "bushing"
    }

    fn name(&self) -> &'static str {
        "bushing"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        BUSHING_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        let resolved = self.resolve(parameters)?;
        let inner_radius_m = resolved.require("inner_radius_m")?;
        let wall_m = resolved.require("wall_m")?;
        let height_m = resolved.require("height_m")?;
        let clearance_m = resolved.require("clearance_m")?;

        let center_m = [0.0, 0.0, 0.0];
        let sleeve = Difference {
            outer: Box::new(Cylinder {
                center_m,
                radius_m: inner_radius_m + wall_m,
                height_m,
            }),
            inner: Box::new(Cylinder {
                center_m,
                radius_m: inner_radius_m + clearance_m,
                height_m,
            }),
        };

        fill_part(&format!("bushing-{height_m:e}"), &sleeve)
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
            "{name} cut no atoms from the diamond lattice"
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

    fn default_bushing_parameters(
        inner_radius_m: f64,
        wall_m: f64,
        height_m: f64,
        clearance_m: f64,
    ) -> ParameterSet {
        ParameterSet::new()
            .with("inner_radius_m", inner_radius_m)
            .with("wall_m", wall_m)
            .with("height_m", height_m)
            .with("clearance_m", clearance_m)
    }

    #[test]
    fn the_bearing_has_rollers() {
        let inner_radius_m = 1.0e-9;
        let outer_radius_m = 9.0e-9;
        let radial_gap_m = 2.0e-9;
        let parameters = ParameterSet::new()
            .with("inner_radius_m", inner_radius_m)
            .with("outer_radius_m", outer_radius_m)
            .with("height_m", 4.0e-9)
            .with("roller_count", 8.0)
            .with("radial_gap_m", radial_gap_m);
        let part = RadialBearingGenerator
            .generate(&parameters)
            .expect("bearing");
        assert!(part.atom_count() > 0, "the bearing is empty");
        let pitch_radius_m = 0.5 * (inner_radius_m + outer_radius_m);
        let nearest_error_m = part
            .topology
            .atoms()
            .map(|atom| (radial_xy_m(atom.position_m) - pitch_radius_m).abs())
            .fold(f64::INFINITY, f64::min);
        assert!(
            nearest_error_m <= radial_gap_m,
            "the nearest atom is {nearest_error_m:e} m from the pitch radius, \
             above the gap {radial_gap_m:e} m"
        );
    }

    #[test]
    fn a_bigger_gap_needs_a_bigger_outer_radius() {
        let inner_radius_m = 3.0e-9;
        let radial_gap_m = 5.0e-10;
        let parameters = ParameterSet::new()
            .with("inner_radius_m", inner_radius_m)
            .with("outer_radius_m", inner_radius_m + 1.5 * radial_gap_m)
            .with("radial_gap_m", radial_gap_m);
        assert!(RadialBearingGenerator.generate(&parameters).is_err());
    }

    #[test]
    fn the_bearing_inner_port_is_revolute() {
        let port = RadialBearingGenerator.inner_port(3.0e-9);
        assert_eq!(port.dof, Dof::Revolute);
        assert_eq!(port.name, "inner");
        assert_eq!(port.origin_m, [0.0, 0.0, 0.0]);
        assert_eq!(port.axis_m, [0.0, 0.0, 1.0]);
        assert_eq!(port.gap_m, 0.0);
        let outer = RadialBearingGenerator.outer_port();
        assert_eq!(outer.name, "outer");
        assert_eq!(outer.dof, Dof::Fixed);
        assert_eq!(outer.axis_m, [0.0, 0.0, 1.0]);
    }

    #[test]
    fn the_bushing_clearance_is_the_shaft_gap() {
        let clearance_m = 2.5e-10;
        let port = BushingGenerator.shaft_port(clearance_m);
        assert_eq!(port.gap_m, clearance_m);
        assert_eq!(port.dof, Dof::Revolute);
        assert_eq!(port.name, "shaft");
        assert_eq!(port.origin_m, [0.0, 0.0, 0.0]);
        assert_eq!(port.axis_m, [0.0, 0.0, 1.0]);
    }

    #[test]
    fn the_bushing_wall_is_solid() {
        let inner_radius_m = 3.0e-9;
        let wall_m = 1.0e-9;
        let parameters = default_bushing_parameters(inner_radius_m, wall_m, 4.0e-9, 0.0);
        let part = BushingGenerator.generate(&parameters).expect("bushing");
        let inner_band_m = inner_radius_m + 1.0e-10;
        let outer_band_m = inner_radius_m + wall_m + 1.0e-10;
        let count = part
            .topology
            .atoms()
            .filter(|atom| {
                let radial_m = radial_xy_m(atom.position_m);
                radial_m >= inner_band_m && radial_m <= outer_band_m
            })
            .count();
        assert!(count > 0, "the sleeve wall holds no atom");
    }

    #[test]
    fn every_carbon_keeps_a_valid_valence() {
        let part = BushingGenerator.generate_with_defaults().expect("bushing");
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
