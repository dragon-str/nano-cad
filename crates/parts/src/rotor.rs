//! The sorting rotor disk.
//!
//! The sorting rotor is the device in Freitas, *Nanomedicine* Volume I,
//! Section 3.4.2, after Drexler. A disk carries a row of binding pockets along
//! its rim. The disk turns about its axis, and each pocket carries its bound
//! molecule from the outside solution to an inner chamber.
//!
//! This generator builds the disk and its pockets. It does not build a binding
//! site, and it does not model a molecule. A pocket is a cylindrical void that
//! opens to the rim.
//!
//! The disk is cut from the diamond cubic lattice, so every bond is the diamond
//! first-shell length and every carbon keeps a valence of four. All lengths are
//! SI metres and all angles are radians.

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

/// The parameter specs of [`SortingRotorGenerator`], in a stable order.
static SORTING_ROTOR_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new(
        "disc_radius_m",
        Some(Unit::Metre),
        7.0e-9,
        2.0e-9,
        4.0e-8,
        false,
        "outside radius of the rotor disk in metres",
    ),
    ParameterSpec::new(
        "thickness_m",
        Some(Unit::Metre),
        2.0e-9,
        7.0e-10,
        2.0e-8,
        false,
        "axial thickness of the disk in metres",
    ),
    ParameterSpec::new(
        "pocket_count",
        None,
        12.0,
        3.0,
        24.0,
        true,
        "number of binding pockets on the rim",
    ),
    ParameterSpec::new(
        "pocket_radius_m",
        Some(Unit::Metre),
        1.0e-9,
        5.0e-10,
        5.0e-9,
        false,
        "radius of each pocket void in metres",
    ),
    ParameterSpec::new(
        "pocket_orbit_m",
        Some(Unit::Metre),
        6.5e-9,
        1.0e-9,
        4.0e-8,
        false,
        "radius of the pocket centre circle in metres",
    ),
    ParameterSpec::new(
        "bore_radius_m",
        Some(Unit::Metre),
        1.5e-9,
        0.0,
        2.0e-8,
        false,
        "radius of the central bore in metres",
    ),
];

/// Generates a sorting rotor disk.
///
/// The disk is the solid between the central bore and the outside radius, less
/// one cylindrical pocket for every slot. A pocket whose centre circle plus its
/// radius reaches past the outside radius opens to the rim, which is the
/// binding face the solution sees.
///
/// The disk turns about `z`, so the axis port lies on that axis.
#[derive(Clone, Copy, Debug, Default)]
pub struct SortingRotorGenerator;

impl SortingRotorGenerator {
    /// Returns the axis port.
    ///
    /// The port sits at the origin, its axis is `+z`, and it allows one
    /// rotation about that axis. The gap is zero.
    pub fn axis_port(&self) -> Port {
        let mut port = Port::named("axis");
        port.dof = Dof::Revolute;
        port
    }

    /// Returns the port of one pocket.
    ///
    /// The port sits at the pocket centre, its axis is `+z`, and it holds the
    /// pocket rigidly. The angle of pocket `index` is `2 pi index / count`.
    pub fn pocket_port(&self, index: usize, count: usize, pocket_orbit_m: f64) -> Port {
        let angle_rad = 2.0 * PI * index as f64 / count.max(1) as f64;
        let (sin_rad, cos_rad) = angle_rad.sin_cos();
        let mut port = Port::named(&format!("pocket_{index}"));
        port.origin_m = [pocket_orbit_m * cos_rad, pocket_orbit_m * sin_rad, 0.0];
        port
    }

    /// Returns the centre of one pocket, in metres.
    pub fn pocket_center_m(index: usize, count: usize, pocket_orbit_m: f64) -> [f64; 3] {
        let angle_rad = 2.0 * PI * index as f64 / count.max(1) as f64;
        let (sin_rad, cos_rad) = angle_rad.sin_cos();
        [pocket_orbit_m * cos_rad, pocket_orbit_m * sin_rad, 0.0]
    }
}

impl PartGenerator for SortingRotorGenerator {
    fn id(&self) -> &'static str {
        "sorting_rotor"
    }

    fn name(&self) -> &'static str {
        "sorting rotor"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        SORTING_ROTOR_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        let resolved = self.resolve(parameters)?;
        let disc_radius_m = resolved.require("disc_radius_m")?;
        let thickness_m = resolved.require("thickness_m")?;
        let pocket_count = resolved.require("pocket_count")? as usize;
        let pocket_radius_m = resolved.require("pocket_radius_m")?;
        let pocket_orbit_m = resolved.require("pocket_orbit_m")?;
        let bore_radius_m = resolved.require("bore_radius_m")?;

        if pocket_orbit_m < bore_radius_m + pocket_radius_m {
            return Err(PartError::InvalidGeometry(format!(
                "the pocket orbit {pocket_orbit_m} m is inside the bore {bore_radius_m} m \
                 plus the pocket radius {pocket_radius_m} m, so a pocket would breach the bore"
            )));
        }
        if pocket_orbit_m - pocket_radius_m >= disc_radius_m {
            return Err(PartError::InvalidGeometry(format!(
                "the pocket orbit {pocket_orbit_m} m minus the pocket radius {pocket_radius_m} m \
                 is at or outside the disc radius {disc_radius_m} m, so the pockets would not bite"
            )));
        }

        let center_m = [0.0, 0.0, 0.0];
        let disc = Cylinder {
            center_m,
            radius_m: disc_radius_m,
            height_m: thickness_m,
        };
        let mut void: Box<dyn Solid> = Box::new(Cylinder {
            center_m,
            radius_m: bore_radius_m,
            height_m: thickness_m,
        });
        for index in 0..pocket_count {
            let pocket = Cylinder {
                center_m: Self::pocket_center_m(index, pocket_count, pocket_orbit_m),
                radius_m: pocket_radius_m,
                height_m: thickness_m,
            };
            void = Box::new(Union {
                a: void,
                b: Box::new(pocket),
            });
        }

        let solid = Difference {
            outer: Box::new(disc),
            inner: void,
        };
        fill_part(
            &format!("sorting-rotor-{pocket_count}"),
            &solid,
            "the disc is thinner than one lattice layer",
        )
    }
}

/// Fills a solid region with the diamond lattice and caps every surface atom.
///
/// The region must have a finite bounding box. The function returns an error
/// when the region cut no lattice site, so an empty body cannot reach a caller.
fn fill_part(name: &str, solid: &dyn Solid, empty_note: &str) -> Result<Part, PartError> {
    let sites_m = fill_solid(solid);
    let mut topology = Topology::new();
    let mut carbons = Vec::with_capacity(sites_m.len());
    for site_m in sites_m {
        carbons.push(topology.add_atom(Atom::new(Element::CARBON, site_m, 0.0, "C")));
    }
    if carbons.is_empty() {
        return Err(PartError::InvalidGeometry(format!(
            "{name} cut no atom from the diamond lattice: {empty_note}"
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
    fn the_rotor_has_a_void_for_every_pocket() {
        let pocket_count = 12;
        let pocket_radius_m = 1.0e-9;
        let pocket_orbit_m = 6.5e-9;
        let part = SortingRotorGenerator
            .generate_with_defaults()
            .expect("rotor");
        for index in 0..pocket_count {
            let center_m =
                SortingRotorGenerator::pocket_center_m(index, pocket_count, pocket_orbit_m);
            let half_m = 0.5 * pocket_radius_m;
            let inside = part.topology.atoms().filter(|atom| {
                let dx = atom.position_m[0] - center_m[0];
                let dy = atom.position_m[1] - center_m[1];
                (dx * dx + dy * dy).sqrt() < half_m
            });
            assert_eq!(
                inside.count(),
                0,
                "pocket {index} holds a carbon at its centre"
            );
        }
    }

    #[test]
    fn the_rotor_respects_its_radius() {
        let disc_radius_m = 7.0e-9;
        let part = SortingRotorGenerator
            .generate_with_defaults()
            .expect("rotor");
        let farthest_m = part
            .topology
            .atoms()
            .map(|atom| radial_xy_m(atom.position_m))
            .fold(0.0_f64, f64::max);
        assert!(
            farthest_m <= disc_radius_m + 1.0e-10,
            "an atom sits at {farthest_m:e} m, past the disc radius {disc_radius_m:e} m"
        );
    }

    #[test]
    fn the_bore_is_open() {
        let bore_radius_m = 1.5e-9;
        let part = SortingRotorGenerator
            .generate_with_defaults()
            .expect("rotor");
        let inside = part
            .topology
            .atoms()
            .filter(|atom| radial_xy_m(atom.position_m) < 0.5 * bore_radius_m);
        assert_eq!(inside.count(), 0, "the bore holds a carbon");
    }

    #[test]
    fn more_pockets_change_the_rotor() {
        let parameters = ParameterSet::new()
            .with("pocket_count", 6.0)
            .with("pocket_radius_m", 1.3e-9)
            .with("pocket_orbit_m", 6.2e-9);
        let six = SortingRotorGenerator.generate(&parameters).expect("six");
        let parameters = ParameterSet::new()
            .with("pocket_count", 12.0)
            .with("pocket_radius_m", 1.3e-9)
            .with("pocket_orbit_m", 6.2e-9);
        let twelve = SortingRotorGenerator.generate(&parameters).expect("twelve");
        assert_ne!(six.atom_count(), twelve.atom_count());
    }

    #[test]
    fn a_pocket_that_breaches_the_bore_is_refused() {
        let parameters = ParameterSet::new()
            .with("pocket_orbit_m", 2.0e-9)
            .with("pocket_radius_m", 1.0e-9)
            .with("bore_radius_m", 1.5e-9);
        assert!(SortingRotorGenerator.generate(&parameters).is_err());
    }

    #[test]
    fn the_axis_port_is_revolute() {
        let port = SortingRotorGenerator.axis_port();
        assert_eq!(port.name, "axis");
        assert_eq!(port.dof, Dof::Revolute);
        assert_eq!(port.origin_m, [0.0, 0.0, 0.0]);
        assert_eq!(port.axis_m, [0.0, 0.0, 1.0]);
        assert_eq!(port.gap_m, 0.0);
    }

    #[test]
    fn a_pocket_port_sits_on_the_pocket_circle() {
        let port = SortingRotorGenerator.pocket_port(3, 12, 6.5e-9);
        assert_eq!(port.name, "pocket_3");
        assert_eq!(port.dof, Dof::Fixed);
        assert!((radial_xy_m(port.origin_m) - 6.5e-9).abs() < 1.0e-18);
        let quarter_turn = SortingRotorGenerator.pocket_port(3, 12, 6.5e-9);
        assert!(quarter_turn.origin_m[0].abs() < 1.0e-15);
        assert!((quarter_turn.origin_m[1] - 6.5e-9).abs() < 1.0e-18);
    }

    #[test]
    fn every_carbon_keeps_a_valid_valence() {
        let part = SortingRotorGenerator
            .generate_with_defaults()
            .expect("rotor");
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
