//! The drive shaft.
//!
//! The rotor must be turned by something. This generator builds the coaxial
//! drive shaft that carries the torque. The shaft lies on the rotor axis. It
//! passes through the central bore of the rotor and through the bore of each
//! housing end plate. A drive key stands on the shaft surface and enters the
//! keyway of the rotor, so the shaft drives the rotor without slip.
//!
//! The reference counts the drive among the parts that the rotor mass must
//! carry: "about 1e5 atoms including housing and a pro-rata share of the drive
//! system" (Freitas, *Nanomedicine* Volume I, Section 3.4.2). The external
//! motor or gear that turns the shaft is out of scope. This generator builds
//! the shaft and the key only. It does not model a bearing, a torque or a
//! sliding fit. The shaft is cut from the diamond cubic lattice. All lengths
//! are SI metres.

use nanocad_model::{Atom, Element, Part, Topology};
use nanocad_units::Unit;

use crate::diamond_solid;
use crate::error::PartError;
use crate::generator::PartGenerator;
use crate::lattice_fill::fill_solid;
use crate::parameter::{ParameterSet, ParameterSpec};
use crate::port::{Dof, Port};
use crate::shape::{Box3, Cylinder, Solid, Union};

/// The parameter specs of [`DriveShaftGenerator`], in a stable order.
static DRIVE_SHAFT_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new(
        "shaft_radius_m",
        Some(Unit::Metre),
        1.3e-9,
        2.0e-10,
        1.0e-8,
        false,
        "radius of the shaft in metres",
    ),
    ParameterSpec::new(
        "shaft_length_m",
        Some(Unit::Metre),
        4.0e-9,
        5.0e-10,
        4.0e-8,
        false,
        "length of the shaft along z in metres",
    ),
    ParameterSpec::new(
        "key_width_m",
        Some(Unit::Metre),
        2.4e-10,
        0.0,
        2.0e-9,
        false,
        "tangential width of the drive key in metres, or zero for no key",
    ),
    ParameterSpec::new(
        "key_height_m",
        Some(Unit::Metre),
        2.4e-10,
        1.0e-11,
        2.0e-9,
        false,
        "radial height of the drive key above the shaft surface in metres",
    ),
    ParameterSpec::new(
        "key_length_m",
        Some(Unit::Metre),
        6.0e-10,
        1.0e-11,
        2.0e-8,
        false,
        "length of the drive key along z in metres",
    ),
];

/// The generator of the drive shaft.
#[derive(Clone, Copy, Debug, Default)]
pub struct DriveShaftGenerator;

impl DriveShaftGenerator {
    /// Returns the port on the shaft axis.
    ///
    /// The shaft shares the rotor axis and turns with the rotor, so the joint
    /// that the scene builds here is revolute.
    pub fn axis_port(&self) -> Port {
        let mut port = Port::named("axis");
        port.dof = Dof::Revolute;
        port
    }

    /// Returns the port at the outer face of the drive key.
    ///
    /// The key seats in the keyway of the rotor, so the port points outward
    /// along `+x` at the shaft radius.
    pub fn key_port(&self) -> Port {
        let mut port = Port::named("key");
        port.origin_m = [self.shaft_radius_m(), 0.0, 0.0];
        port.axis_m = [1.0, 0.0, 0.0];
        port.dof = Dof::Fixed;
        port
    }

    /// Returns the radius of the shaft, in metres.
    pub fn shaft_radius_m(&self) -> f64 {
        self.default_m("shaft_radius_m")
    }

    /// Returns the length of the shaft along z, in metres.
    pub fn shaft_length_m(&self) -> f64 {
        self.default_m("shaft_length_m")
    }

    /// Returns the radial height of the key above the shaft, in metres.
    pub fn key_height_m(&self) -> f64 {
        self.default_m("key_height_m")
    }

    fn default_m(&self, name: &str) -> f64 {
        self.resolve(&ParameterSet::new())
            .ok()
            .and_then(|resolved| resolved.get(name))
            .unwrap_or(0.0)
    }
}

impl PartGenerator for DriveShaftGenerator {
    fn id(&self) -> &'static str {
        "drive_shaft"
    }

    fn name(&self) -> &'static str {
        "drive shaft"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        DRIVE_SHAFT_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        let resolved = self.resolve(parameters)?;
        let shaft_radius_m = resolved.require("shaft_radius_m")?;
        let shaft_length_m = resolved.require("shaft_length_m")?;
        let key_width_m = resolved.require("key_width_m")?;
        let key_height_m = resolved.require("key_height_m")?;
        let key_length_m = resolved.require("key_length_m")?;

        if shaft_radius_m <= 0.0 || shaft_length_m <= 0.0 {
            return Err(PartError::InvalidGeometry(
                "the shaft must have a radius and a length".to_string(),
            ));
        }
        if key_width_m > 0.0 && key_length_m > shaft_length_m {
            return Err(PartError::InvalidGeometry(format!(
                "the key length {key_length_m} m is longer than the shaft \
                 {shaft_length_m} m"
            )));
        }

        let shaft: Box<dyn Solid> = Box::new(Cylinder {
            center_m: [0.0, 0.0, 0.0],
            radius_m: shaft_radius_m,
            height_m: shaft_length_m,
        });
        let solid: Box<dyn Solid> = if key_width_m > 0.0 {
            let key = Box3 {
                min_m: [
                    shaft_radius_m - key_height_m,
                    -0.5 * key_width_m,
                    -0.5 * key_length_m,
                ],
                max_m: [
                    shaft_radius_m + key_height_m,
                    0.5 * key_width_m,
                    0.5 * key_length_m,
                ],
            };
            Box::new(Union {
                a: shaft,
                b: Box::new(key),
            })
        } else {
            shaft
        };

        let sites_m = fill_solid(&*solid);
        if sites_m.is_empty() {
            return Err(PartError::InvalidGeometry(
                "the shaft parameters leave no lattice site inside the shaft".to_string(),
            ));
        }

        let mut topology = Topology::new();
        let mut carbons = Vec::with_capacity(sites_m.len());
        for site_m in sites_m {
            carbons.push(topology.add_atom(Atom::new(Element::CARBON, site_m, 0.0, "C")));
        }
        diamond_solid::bond_and_cap(&mut topology, &carbons)?;

        Ok(Part::new("drive_shaft", topology).with_material("diamond"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::port::Dof;

    fn shaft() -> Part {
        DriveShaftGenerator
            .generate_with_defaults()
            .expect("the shaft builds")
    }

    fn carbon_count(part: &Part) -> usize {
        (0..part.atom_count())
            .filter(|index| part.topology.element(*index) == Some(Element::CARBON))
            .count()
    }

    #[test]
    fn the_shaft_builds_with_atoms() {
        let part = shaft();
        assert!(part.atom_count() > 0);
        assert!(part.bond_count() > 0);
    }

    #[test]
    fn the_shaft_stands_on_the_axis() {
        let part = shaft();
        let half_length_m = 0.5 * DriveShaftGenerator.shaft_length_m();
        for index in 0..part.atom_count() {
            if part.topology.element(index) != Some(Element::CARBON) {
                continue;
            }
            let Some(position_m) = part.topology.position_m(index) else {
                continue;
            };
            assert!(
                position_m[2].abs() <= half_length_m + 1.0e-10,
                "an atom sits at z {} m, past the shaft end {half_length_m} m",
                position_m[2]
            );
        }
    }

    #[test]
    fn the_key_stands_above_the_shaft_surface() {
        let part = shaft();
        let radius_m = DriveShaftGenerator.shaft_radius_m();
        let height_m = DriveShaftGenerator.key_height_m();
        let mut key_carbons = 0;
        for index in 0..part.atom_count() {
            if part.topology.element(index) != Some(Element::CARBON) {
                continue;
            }
            let Some(position_m) = part.topology.position_m(index) else {
                continue;
            };
            if position_m[0] > radius_m + 1.0e-11 {
                key_carbons += 1;
                assert!(
                    position_m[0] <= radius_m + height_m + 1.0e-10,
                    "a key atom reaches {} m, past the key height",
                    position_m[0]
                );
            }
        }
        assert!(key_carbons > 0, "the key adds no carbon above the shaft");
    }

    #[test]
    fn a_zero_key_removes_material() {
        let parameters = ParameterSet::new().with("key_width_m", 0.0);
        let plain = DriveShaftGenerator
            .generate(&parameters)
            .expect("the plain shaft builds");
        assert!(
            carbon_count(&shaft()) > carbon_count(&plain),
            "the key adds no carbon: {} plain, {} keyed",
            carbon_count(&plain),
            carbon_count(&shaft())
        );
    }

    #[test]
    fn the_axis_port_turns() {
        let port = DriveShaftGenerator.axis_port();
        assert_eq!(port.name, "axis");
        assert_eq!(port.axis_m, [0.0, 0.0, 1.0]);
        assert_eq!(port.dof, Dof::Revolute);
    }

    #[test]
    fn the_key_port_points_outward() {
        let port = DriveShaftGenerator.key_port();
        assert_eq!(port.name, "key");
        assert_eq!(port.axis_m, [1.0, 0.0, 0.0]);
        assert!((port.origin_m[0] - DriveShaftGenerator.shaft_radius_m()).abs() < 1.0e-18);
        assert_eq!(port.dof, Dof::Fixed);
    }

    #[test]
    fn a_key_longer_than_the_shaft_is_refused() {
        let parameters = ParameterSet::new().with("key_length_m", 5.0e-9);
        assert!(DriveShaftGenerator.generate(&parameters).is_err());
    }

    #[test]
    fn every_carbon_keeps_a_valid_valence() {
        let part = shaft();
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
