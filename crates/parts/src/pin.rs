//! The follower pin.
//!
//! The sorting rotor in Freitas, *Nanomedicine* Volume I, Section 3.4.2, drives
//! its ejection rods with a cam surface. A cam needs a follower: a small pin
//! that rides in a groove and takes the thrust of the track.
//!
//! This generator builds the follower pin. It is a stepped column along `z`: a
//! narrow neck that passes through the opening of a groove, and a wider head
//! that the walls of the groove push. The pin is cut from the diamond cubic
//! lattice, so every bond is the diamond first-shell length.
//!
//! The pin is a schematic. A diamond lattice cannot resolve a feature much
//! below 2e-10 m, so this part is a few tens of atoms and it is not a machined
//! bearing surface. This generator does not build the groove and it does not
//! model the sliding fit or the friction of the contact. All lengths are SI
//! metres and all angles are radians.

use nanocad_model::{Atom, Element, Part, Topology};
use nanocad_units::Unit;

use crate::diamond_solid;
use crate::error::PartError;
use crate::generator::PartGenerator;
use crate::lattice_fill::fill_solid;
use crate::parameter::{ParameterSet, ParameterSpec};
use crate::port::Port;
use crate::shape::{Cylinder, Union};

/// The parameter specs of [`FollowerPinGenerator`], in a stable order.
static FOLLOWER_PIN_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new(
        "neck_radius_m",
        Some(Unit::Metre),
        2.5e-10,
        1.0e-10,
        1.5e-9,
        false,
        "radius of the narrow neck in metres",
    ),
    ParameterSpec::new(
        "neck_length_m",
        Some(Unit::Metre),
        2.0e-10,
        1.0e-10,
        2.0e-9,
        false,
        "length of the narrow neck in metres",
    ),
    ParameterSpec::new(
        "head_radius_m",
        Some(Unit::Metre),
        3.5e-10,
        1.0e-10,
        2.0e-9,
        false,
        "radius of the wide head that the groove walls push, in metres",
    ),
    ParameterSpec::new(
        "head_length_m",
        Some(Unit::Metre),
        4.0e-10,
        1.0e-10,
        3.0e-9,
        false,
        "length of the wide head in metres",
    ),
];

/// The generator of the follower pin.
#[derive(Clone, Copy, Debug, Default)]
pub struct FollowerPinGenerator;

impl FollowerPinGenerator {
    /// Returns the port where the pin meets the rod.
    ///
    /// The port sits at the top of the pin, its axis points out along `+z`,
    /// and the joint is rigid. The top rests on the lower face of the rotor.
    pub fn rod_port(&self) -> Port {
        Port::named("rod")
    }

    /// Returns the port at the working surface of the head.
    ///
    /// The port sits at the bottom of the pin, its axis points out along `-z`,
    /// and the joint is rigid. The head rides in the cam groove here.
    pub fn groove_port(&self) -> Port {
        let mut port = Port::named("groove");
        port.origin_m = [0.0, 0.0, -self.length_m()];
        port.axis_m = [0.0, 0.0, -1.0];
        port
    }

    /// Returns the total length of the pin, in metres.
    pub fn length_m(&self) -> f64 {
        self.neck_length_m() + self.head_length_m()
    }

    /// Returns the radius of the head, in metres.
    pub fn head_radius_m(&self) -> f64 {
        self.default_m("head_radius_m")
    }

    /// Returns the radius of the neck, in metres.
    pub fn neck_radius_m(&self) -> f64 {
        self.default_m("neck_radius_m")
    }

    fn neck_length_m(&self) -> f64 {
        self.default_m("neck_length_m")
    }

    fn head_length_m(&self) -> f64 {
        self.default_m("head_length_m")
    }

    fn default_m(&self, name: &str) -> f64 {
        self.resolve(&ParameterSet::new())
            .ok()
            .and_then(|resolved| resolved.get(name))
            .unwrap_or(0.0)
    }
}

impl PartGenerator for FollowerPinGenerator {
    fn id(&self) -> &'static str {
        "follower_pin"
    }

    fn name(&self) -> &'static str {
        "follower pin"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        FOLLOWER_PIN_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        let resolved = self.resolve(parameters)?;
        let neck_radius_m = resolved.require("neck_radius_m")?;
        let neck_length_m = resolved.require("neck_length_m")?;
        let head_radius_m = resolved.require("head_radius_m")?;
        let head_length_m = resolved.require("head_length_m")?;

        if head_radius_m < neck_radius_m {
            return Err(PartError::InvalidGeometry(format!(
                "the head radius {head_radius_m} m is below the neck radius \
                 {neck_radius_m} m, so the pin has a waist and not a head"
            )));
        }

        let neck = Cylinder {
            center_m: [0.0, 0.0, -0.5 * neck_length_m],
            radius_m: neck_radius_m,
            height_m: neck_length_m,
        };
        let head = Cylinder {
            center_m: [0.0, 0.0, -neck_length_m - 0.5 * head_length_m],
            radius_m: head_radius_m,
            height_m: head_length_m,
        };
        let solid = Union {
            a: Box::new(neck),
            b: Box::new(head),
        };

        let sites_m = fill_solid(&solid);
        if sites_m.is_empty() {
            return Err(PartError::InvalidGeometry(
                "the pin parameters leave no lattice site inside the pin".to_string(),
            ));
        }

        let mut topology = Topology::new();
        let mut carbons = Vec::with_capacity(sites_m.len());
        for site_m in sites_m {
            carbons.push(topology.add_atom(Atom::new(Element::CARBON, site_m, 0.0, "C")));
        }
        diamond_solid::bond_and_cap(&mut topology, &carbons)?;

        Ok(Part::new("follower_pin", topology).with_material("diamond"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::port::Dof;

    fn pin() -> Part {
        FollowerPinGenerator
            .generate_with_defaults()
            .expect("the pin builds")
    }

    #[test]
    fn the_pin_builds_with_atoms() {
        let part = pin();
        assert!(part.atom_count() > 0);
        assert!(part.bond_count() > 0);
    }

    #[test]
    fn the_pin_is_stepped() {
        let part = pin();
        let head_radius_m = FollowerPinGenerator.head_radius_m();
        let neck_radius_m = FollowerPinGenerator.neck_radius_m();
        let mut largest_m: f64 = 0.0;
        for index in 0..part.atom_count() {
            let Some(position_m) = part.topology.position_m(index) else {
                continue;
            };
            let radius_m = position_m[0].hypot(position_m[1]);
            largest_m = largest_m.max(radius_m);
        }
        let lattice_step_m = crate::diamond_solid::diamond_bond_length_m();
        assert!(
            largest_m <= head_radius_m + lattice_step_m,
            "the widest atom sits at {largest_m} m, well outside the head"
        );
        assert!(largest_m > neck_radius_m);
    }

    #[test]
    fn the_pin_hangs_below_its_top_face() {
        let part = pin();
        let expected_m = FollowerPinGenerator.length_m();
        let mut low_m = f64::INFINITY;
        let mut high_m = f64::NEG_INFINITY;
        for index in 0..part.atom_count() {
            let Some(position_m) = part.topology.position_m(index) else {
                continue;
            };
            low_m = low_m.min(position_m[2]);
            high_m = high_m.max(position_m[2]);
        }
        assert!(high_m <= 1.0e-10, "the pin rises above its top face");
        assert!(low_m >= -expected_m - 1.0e-10);
        assert!(high_m - low_m > 0.8 * expected_m);
    }

    #[test]
    fn the_rod_port_sits_on_the_top_face() {
        let port = FollowerPinGenerator.rod_port();
        assert_eq!(port.origin_m, [0.0, 0.0, 0.0]);
        assert_eq!(port.axis_m, [0.0, 0.0, 1.0]);
        assert_eq!(port.dof, Dof::Fixed);
    }

    #[test]
    fn the_groove_port_points_along_negative_z() {
        let port = FollowerPinGenerator.groove_port();
        assert_eq!(port.axis_m, [0.0, 0.0, -1.0]);
        assert_eq!(port.dof, Dof::Fixed);
        assert!(port.origin_m[2] < 0.0);
    }

    #[test]
    fn a_head_narrower_than_the_neck_is_refused() {
        let mut parameters = ParameterSet::new();
        parameters.set("head_radius_m", 1.0e-10);
        parameters.set("neck_radius_m", 3.0e-10);
        assert!(FollowerPinGenerator.generate(&parameters).is_err());
    }

    #[test]
    fn every_carbon_keeps_a_valid_valence() {
        let part = pin();
        let mut degree = vec![0_usize; part.atom_count()];
        for bond in part.topology.bonds() {
            degree[bond.u as usize] += 1;
            degree[bond.v as usize] += 1;
        }
        for (index, count) in degree.iter().enumerate() {
            let element = part.topology.element(index).expect("an element");
            let expected = nanocad_model::valence(element).expect("a valence");
            assert!(
                *count <= expected as usize,
                "atom {index} has degree {count} above {expected}"
            );
        }
    }
}
