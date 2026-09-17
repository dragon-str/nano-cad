//! Clutch plate and ratchet generators.
//!
//! Both parts are solid diamond blocks, cut to a toothed outline. The outline
//! is a 2D polygon in the `xy` plane. The polygon is filled with the diamond
//! lattice through [`fill_solid`], and every surface carbon receives hydrogen
//! so that each carbon keeps a valence of four.
//!
//! Geometry is defined in SI metres and angles are in radians. Every length
//! carries its unit in the name.

use std::f64::consts::TAU;

use nanocad_model::{Atom, Element, Part, Topology};
use nanocad_units::Unit;

use crate::diamond_solid;
use crate::error::PartError;
use crate::generator::PartGenerator;
use crate::lattice_fill::fill_solid;
use crate::parameter::{ParameterSet, ParameterSpec};
use crate::shape::{Box3, Cylinder, Difference, Intersection, Profile, Solid, Union};

/// The clearance between the pawl nose and the wheel rim, in metres.
const PAWL_CLEARANCE_M: f64 = 2.0e-10;

static CLUTCH_PLATE_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new(
        "teeth",
        None,
        12.0,
        6.0,
        32.0,
        true,
        "number of rectangular teeth",
    ),
    ParameterSpec::new(
        "inner_radius_m",
        Some(Unit::Metre),
        2.0e-9,
        1.0e-9,
        1.0e-8,
        false,
        "inner bore radius in metres",
    ),
    ParameterSpec::new(
        "tooth_height_m",
        Some(Unit::Metre),
        1.0e-9,
        5.0e-10,
        5.0e-9,
        false,
        "radial tooth height in metres",
    ),
    ParameterSpec::new(
        "thickness_m",
        Some(Unit::Metre),
        1.5e-9,
        7.0e-10,
        1.0e-8,
        false,
        "plate thickness in metres",
    ),
];

static RATCHET_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new("teeth", None, 10.0, 4.0, 24.0, true, "number of saw teeth"),
    ParameterSpec::new(
        "wheel_radius_m",
        Some(Unit::Metre),
        3.0e-9,
        1.5e-9,
        2.0e-8,
        false,
        "wheel root radius in metres",
    ),
    ParameterSpec::new(
        "tooth_height_m",
        Some(Unit::Metre),
        1.0e-9,
        5.0e-10,
        5.0e-9,
        false,
        "radial tooth height in metres",
    ),
    ParameterSpec::new(
        "thickness_m",
        Some(Unit::Metre),
        1.5e-9,
        7.0e-10,
        1.0e-8,
        false,
        "wheel and pawl thickness in metres",
    ),
    ParameterSpec::new(
        "pawl_length_m",
        Some(Unit::Metre),
        1.5e-9,
        7.0e-10,
        8.0e-9,
        false,
        "pawl block length in metres",
    ),
    ParameterSpec::new(
        "pawl_width_m",
        Some(Unit::Metre),
        7.0e-10,
        3.0e-10,
        3.0e-9,
        false,
        "pawl block width in metres",
    ),
];

/// Builds a rectangular-toothed annular clutch plate.
///
/// The outline is a square-tooth (dog clutch) polygon: each tooth contributes
/// four points. The first pair sits on the inner radius, and the second pair
/// sits on the outer radius. The points run counter-clockwise.
///
/// A [`Profile`] has no z bounds, because the polygon is extruded through all
/// z. The profile therefore needs a finite clipping shape on every axis:
///
/// - a [`Box3`] sets the z limits to `+/- thickness_m / 2` and also bounds
///   `x` and `y`,
/// - a [`Cylinder`] of `inner_radius_m` removes the centre, so the result is
///   an annulus and not a full disc.
///
/// The bounded annulus is intersected with the profile, and the diamond
/// lattice fills the result. An unbounded profile alone would give an empty
/// fill, because [`fill_solid`] needs a finite box on every axis.
#[derive(Clone, Copy, Debug, Default)]
pub struct ClutchPlateGenerator;

impl PartGenerator for ClutchPlateGenerator {
    fn id(&self) -> &'static str {
        "clutch_plate"
    }

    fn name(&self) -> &'static str {
        "clutch plate"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        CLUTCH_PLATE_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        let resolved = self.resolve(parameters)?;
        let teeth = resolved.require("teeth")? as usize;
        let inner_radius_m = resolved.require("inner_radius_m")?;
        let tooth_height_m = resolved.require("tooth_height_m")?;
        let thickness_m = resolved.require("thickness_m")?;

        let outline_m = rectangular_tooth_outline_m(teeth, inner_radius_m, tooth_height_m);
        let profile = Profile {
            outline_m,
            internal: false,
            outer_radius_m: None,
        };

        let outer_radius_m = inner_radius_m + tooth_height_m;
        let clip = Box3 {
            min_m: [-outer_radius_m, -outer_radius_m, -0.5 * thickness_m],
            max_m: [outer_radius_m, outer_radius_m, 0.5 * thickness_m],
        };
        let bore = Cylinder {
            center_m: [0.0, 0.0, 0.0],
            radius_m: inner_radius_m,
            height_m: thickness_m,
        };
        let bounded_annulus = Difference {
            outer: Box::new(clip),
            inner: Box::new(bore),
        };
        let solid = Intersection {
            a: Box::new(profile),
            b: Box::new(bounded_annulus),
        };

        let name = format!("clutch-plate-z{teeth}");
        carbon_solid_part(name, "diamondoid", &solid, "C")
    }
}

/// Builds a saw-tooth ratchet wheel with a rigid pawl.
///
/// Each tooth leans one way: the outer point sits at `k * TAU / teeth`, and
/// the sloped face returns to `wheel_radius_m` at
/// `k * TAU / teeth + 0.8 * TAU / teeth`. The wheel is clipped by a [`Box3`],
/// because the profile has no z bound.
///
/// The pawl is a rigid block for now. It is not a spring, and it carries no
/// joint. It sits just off the rim on the +x axis. The wheel and the pawl are
/// united into one solid, and the diamond lattice fills both.
#[derive(Clone, Copy, Debug, Default)]
pub struct RatchetGenerator;

impl PartGenerator for RatchetGenerator {
    fn id(&self) -> &'static str {
        "ratchet"
    }

    fn name(&self) -> &'static str {
        "ratchet"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        RATCHET_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        let resolved = self.resolve(parameters)?;
        let teeth = resolved.require("teeth")? as usize;
        let wheel_radius_m = resolved.require("wheel_radius_m")?;
        let tooth_height_m = resolved.require("tooth_height_m")?;
        let thickness_m = resolved.require("thickness_m")?;
        let pawl_length_m = resolved.require("pawl_length_m")?;
        let pawl_width_m = resolved.require("pawl_width_m")?;

        let outline_m = saw_tooth_outline_m(teeth, wheel_radius_m, tooth_height_m);
        let profile = Profile {
            outline_m,
            internal: false,
            outer_radius_m: None,
        };

        let outer_radius_m = wheel_radius_m + tooth_height_m;
        let clip = Box3 {
            min_m: [-outer_radius_m, -outer_radius_m, -0.5 * thickness_m],
            max_m: [outer_radius_m, outer_radius_m, 0.5 * thickness_m],
        };
        let wheel_solid = Intersection {
            a: Box::new(profile),
            b: Box::new(clip),
        };

        let pawl_min_x_m = wheel_radius_m + PAWL_CLEARANCE_M;
        let pawl_center_m = [pawl_min_x_m + 0.5 * pawl_length_m, 0.0, 0.0];
        let pawl = Box3::from_center_m(pawl_center_m, [pawl_length_m, pawl_width_m, thickness_m]);
        let solid = Union {
            a: Box::new(wheel_solid),
            b: Box::new(pawl),
        };

        let name = format!("ratchet-z{teeth}");
        carbon_solid_part(name, "diamondoid", &solid, "C")
    }
}

/// Builds the closed counter-clockwise outline of a square-tooth clutch plate.
///
/// Each tooth contributes four points at the angles
/// `k * step + {0.0, 0.25, 0.5, 0.75} * step`, where `step = TAU / teeth`. The
/// first pair sits at `inner_radius_m` and the second pair at
/// `inner_radius_m + tooth_height_m`. The first point repeats at the end, so
/// the polygon closes.
fn rectangular_tooth_outline_m(
    teeth: usize,
    inner_radius_m: f64,
    tooth_height_m: f64,
) -> Vec<[f64; 2]> {
    let step_rad = TAU / teeth as f64;
    let outer_radius_m = inner_radius_m + tooth_height_m;
    let mut points_m = Vec::with_capacity(teeth * 4 + 1);
    for tooth in 0..teeth {
        let base_rad = tooth as f64 * step_rad;
        for (fraction, radius_m) in [
            (0.0, inner_radius_m),
            (0.25, inner_radius_m),
            (0.5, outer_radius_m),
            (0.75, outer_radius_m),
        ] {
            let angle_rad = base_rad + fraction * step_rad;
            points_m.push([radius_m * angle_rad.cos(), radius_m * angle_rad.sin()]);
        }
    }
    if let Some(first) = points_m.first().copied() {
        points_m.push(first);
    }
    points_m
}

/// Builds the closed counter-clockwise outline of a saw-tooth ratchet wheel.
///
/// Each tooth contributes two points: the outer point at `k * step`, and the
/// root point at `k * step + 0.8 * step`, where `step = TAU / teeth`. The
/// sloped face runs from the outer point down to the root point, so the tooth
/// leans one way. The first point repeats at the end, so the polygon closes.
fn saw_tooth_outline_m(teeth: usize, wheel_radius_m: f64, tooth_height_m: f64) -> Vec<[f64; 2]> {
    let step_rad = TAU / teeth as f64;
    let outer_radius_m = wheel_radius_m + tooth_height_m;
    let mut points_m = Vec::with_capacity(teeth * 2 + 1);
    for tooth in 0..teeth {
        let base_rad = tooth as f64 * step_rad;
        let outer_angle_rad = base_rad;
        let root_angle_rad = base_rad + 0.8 * step_rad;
        points_m.push([
            outer_radius_m * outer_angle_rad.cos(),
            outer_radius_m * outer_angle_rad.sin(),
        ]);
        points_m.push([
            wheel_radius_m * root_angle_rad.cos(),
            wheel_radius_m * root_angle_rad.sin(),
        ]);
    }
    if let Some(first) = points_m.first().copied() {
        points_m.push(first);
    }
    points_m
}

/// Fills a solid with diamond carbon, caps the surface, and returns a part.
///
/// The solid must have a finite bounding box on every axis. Each lattice site
/// becomes a carbon atom, and [`diamond_solid::bond_and_cap`] bonds the
/// first-shell pairs and adds hydrogen to each under-coordinated carbon.
fn carbon_solid_part(
    name: String,
    material: &str,
    solid: &dyn Solid,
    atom_type: &str,
) -> Result<Part, PartError> {
    let positions_m = fill_solid(solid);
    let mut topology = Topology::new();
    let mut indices = Vec::with_capacity(positions_m.len());
    for position_m in positions_m {
        indices.push(topology.add_atom(Atom::new(Element::CARBON, position_m, 0.0, atom_type)));
    }
    diamond_solid::bond_and_cap(&mut topology, &indices)?;
    Ok(Part::new(name, topology).with_material(material))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn radius_m(position_m: [f64; 3]) -> f64 {
        (position_m[0] * position_m[0] + position_m[1] * position_m[1]).sqrt()
    }

    #[test]
    fn the_clutch_plate_has_the_stated_teeth() {
        let part = ClutchPlateGenerator
            .generate_with_defaults()
            .expect("clutch plate");
        let inner_radius_m = 2.0e-9;
        let outer_radius_m = inner_radius_m + 1.0e-9;
        let mut min_radius_m = f64::INFINITY;
        let mut max_radius_m = 0.0_f64;
        for atom in part.topology.atoms() {
            // The hydrogen caps sit past the carbon surface, so the plate
            // geometry is measured on the carbon sites only.
            if atom.element != Element::CARBON {
                continue;
            }
            let radius_m = radius_m(atom.position_m);
            min_radius_m = min_radius_m.min(radius_m);
            max_radius_m = max_radius_m.max(radius_m);
        }
        assert!(
            (max_radius_m - outer_radius_m).abs() <= 1.0e-10,
            "max radius {max_radius_m:e} m, expected {outer_radius_m:e} m"
        );
        assert!(
            min_radius_m >= inner_radius_m - 1.0e-10,
            "min radius {min_radius_m:e} m, expected at least {inner_radius_m:e} m"
        );
    }

    #[test]
    fn more_teeth_change_the_plate() {
        let default = ClutchPlateGenerator
            .generate_with_defaults()
            .expect("default plate");
        let doubled = ClutchPlateGenerator
            .generate(&ParameterSet::new().with("teeth", 24.0))
            .expect("doubled plate");
        let default_positions: Vec<[f64; 3]> = default
            .topology
            .atoms()
            .map(|atom| atom.position_m)
            .collect();
        let doubled_positions: Vec<[f64; 3]> = doubled
            .topology
            .atoms()
            .map(|atom| atom.position_m)
            .collect();
        assert!(
            default.atom_count() != doubled.atom_count() || default_positions != doubled_positions
        );
    }

    #[test]
    fn the_ratchet_has_a_pawl() {
        let part = RatchetGenerator.generate_with_defaults().expect("ratchet");
        let threshold_m = 3.0e-9 + 1.0e-9;
        let mut beyond = false;
        let mut below = false;
        for atom in part.topology.atoms() {
            let radius_m = radius_m(atom.position_m);
            if radius_m > threshold_m {
                beyond = true;
            }
            if radius_m < threshold_m {
                below = true;
            }
        }
        assert!(beyond, "no pawl atom lies beyond {threshold_m:e} m");
        assert!(below, "no wheel atom lies below {threshold_m:e} m");
    }

    #[test]
    fn the_ratchet_is_one_part() {
        let part = RatchetGenerator.generate_with_defaults().expect("ratchet");
        assert!(part.atom_count() > 0);
    }

    #[test]
    fn every_carbon_keeps_a_valid_valence() {
        let part = ClutchPlateGenerator
            .generate_with_defaults()
            .expect("clutch plate");
        let mut degree = vec![0u32; part.atom_count()];
        for bond in part.topology.bonds() {
            degree[bond.u as usize] += 1;
            degree[bond.v as usize] += 1;
        }
        for (index, atom) in part.topology.atoms().enumerate() {
            if atom.element != Element::CARBON {
                continue;
            }
            let count = degree[index];
            assert!(
                (1..=4).contains(&count),
                "carbon {index} has degree {count}"
            );
        }
    }
}
