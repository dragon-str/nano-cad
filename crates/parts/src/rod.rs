//! The ejection rod.
//!
//! The sorting rotor in Freitas, *Nanomedicine* Volume I, Section 3.4.2, holds
//! its bound molecule in a pocket until the pocket reaches the inner chamber.
//! There the molecule is forced out: the reference says the bound molecules
//! "are forcibly ejected by rods thrust outward by the cam surface".
//!
//! This generator builds the rod. It is a stepped column: a wide shaft that
//! slides in a bore, and a narrow tip that reaches into the pocket and pushes
//! the guest. The rod is cut from the diamond cubic lattice, so every bond is
//! the diamond first-shell length and the stress path is a solid diamond body.
//!
//! This generator does not build the cam. It does not model the sliding fit,
//! the friction of the shaft, or the force on the rod. Those are measurements,
//! not geometry, and [`crate::ejection`] holds the one this milestone needs.
//! All lengths are SI metres and all angles are radians.

use nanocad_model::{Atom, Element, Part, Topology};
use nanocad_units::Unit;

use crate::diamond_solid;
use crate::error::PartError;
use crate::generator::PartGenerator;
use crate::lattice_fill::fill_solid;
use crate::parameter::{ParameterSet, ParameterSpec};
use crate::port::{Dof, Port};
use crate::shape::{Cylinder, Union};

/// The parameter specs of [`EjectionRodGenerator`], in a stable order.
static EJECTION_ROD_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new(
        "shaft_radius_m",
        Some(Unit::Metre),
        6.5e-10,
        1.5e-10,
        2.0e-9,
        false,
        "radius of the sliding shaft in metres",
    ),
    ParameterSpec::new(
        "shaft_length_m",
        Some(Unit::Metre),
        1.6e-9,
        3.0e-10,
        2.0e-8,
        false,
        "length of the sliding shaft in metres",
    ),
    ParameterSpec::new(
        "tip_radius_m",
        Some(Unit::Metre),
        4.0e-10,
        1.0e-10,
        1.0e-9,
        false,
        "radius of the ejection tip in metres",
    ),
    ParameterSpec::new(
        "tip_length_m",
        Some(Unit::Metre),
        4.0e-10,
        1.0e-10,
        3.0e-9,
        false,
        "length of the ejection tip in metres",
    ),
];

/// The generator of the ejection rod.
#[derive(Clone, Copy, Debug, Default)]
pub struct EjectionRodGenerator;

impl EjectionRodGenerator {
    /// Returns the port at the free end of the shaft.
    ///
    /// The port sits at the bottom of the shaft, its axis points out along
    /// `-z`, and it allows the slide along that axis. The cam pushes here.
    pub fn cam_port(&self) -> Port {
        let mut port = Port::named("cam");
        port.origin_m = [0.0, 0.0, -self.half_shaft_m()];
        port.axis_m = [0.0, 0.0, -1.0];
        port.dof = Dof::Prismatic;
        port
    }

    /// Returns the port at the working face of the tip.
    ///
    /// The port sits at the top of the tip, its axis points out along `+z`,
    /// and the joint is rigid.
    pub fn tip_port(&self) -> Port {
        let mut port = Port::named("tip");
        port.origin_m = [0.0, 0.0, self.half_shaft_m() + self.tip_length_m()];
        port
    }

    /// Returns the total length of the rod, in metres.
    pub fn length_m(&self) -> f64 {
        2.0 * self.half_shaft_m() + self.tip_length_m()
    }

    /// Returns the outer diameter of the shaft, in metres.
    pub fn shaft_diameter_m(&self) -> f64 {
        2.0 * self.shaft_radius_m()
    }

    fn half_shaft_m(&self) -> f64 {
        0.5 * self.shaft_length_m()
    }

    fn shaft_radius_m(&self) -> f64 {
        self.default_m("shaft_radius_m")
    }

    fn shaft_length_m(&self) -> f64 {
        self.default_m("shaft_length_m")
    }

    fn tip_length_m(&self) -> f64 {
        self.default_m("tip_length_m")
    }

    fn default_m(&self, name: &str) -> f64 {
        self.resolve(&ParameterSet::new())
            .ok()
            .and_then(|resolved| resolved.get(name))
            .unwrap_or(0.0)
    }
}

impl PartGenerator for EjectionRodGenerator {
    fn id(&self) -> &'static str {
        "ejection_rod"
    }

    fn name(&self) -> &'static str {
        "ejection rod"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        EJECTION_ROD_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        let resolved = self.resolve(parameters)?;
        let shaft_radius_m = resolved.require("shaft_radius_m")?;
        let shaft_length_m = resolved.require("shaft_length_m")?;
        let tip_radius_m = resolved.require("tip_radius_m")?;
        let tip_length_m = resolved.require("tip_length_m")?;

        if tip_radius_m >= shaft_radius_m {
            return Err(PartError::InvalidGeometry(format!(
                "the tip radius {tip_radius_m} m is at or above the shaft radius \
                 {shaft_radius_m} m, so the rod is not stepped and the tip cannot enter a bore"
            )));
        }

        let shaft = Cylinder {
            center_m: [0.0, 0.0, 0.0],
            radius_m: shaft_radius_m,
            height_m: shaft_length_m,
        };
        let tip = Cylinder {
            center_m: [0.0, 0.0, 0.5 * shaft_length_m + 0.5 * tip_length_m],
            radius_m: tip_radius_m,
            height_m: tip_length_m,
        };
        let solid = Union {
            a: Box::new(shaft),
            b: Box::new(tip),
        };

        let sites_m = fill_solid(&solid);
        if sites_m.is_empty() {
            return Err(PartError::InvalidGeometry(
                "the rod parameters leave no lattice site inside the rod".to_string(),
            ));
        }

        let mut topology = Topology::new();
        let mut carbons = Vec::with_capacity(sites_m.len());
        for site_m in sites_m {
            carbons.push(topology.add_atom(Atom::new(Element::CARBON, site_m, 0.0, "C")));
        }
        diamond_solid::bond_and_cap(&mut topology, &carbons)?;

        Ok(Part::new("ejection_rod", topology).with_material("diamond"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rod() -> Part {
        EjectionRodGenerator
            .generate_with_defaults()
            .expect("the rod builds")
    }

    #[test]
    fn the_rod_builds_with_atoms() {
        let part = rod();
        assert!(part.atom_count() > 0);
        assert!(part.bond_count() > 0);
    }

    #[test]
    fn the_rod_is_stepped() {
        let part = rod();
        let shaft_radius_m = EjectionRodGenerator.shaft_radius_m();
        let tip_radius_m = EjectionRodGenerator.default_m("tip_radius_m");
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
            largest_m <= shaft_radius_m + lattice_step_m,
            "the widest atom sits at {largest_m} m, well outside the shaft"
        );
        assert!(largest_m > tip_radius_m);
    }

    #[test]
    fn the_rod_spans_its_stated_length() {
        let part = rod();
        let expected_m = EjectionRodGenerator.length_m();
        let mut low_m = f64::INFINITY;
        let mut high_m = f64::NEG_INFINITY;
        for index in 0..part.atom_count() {
            let Some(position_m) = part.topology.position_m(index) else {
                continue;
            };
            low_m = low_m.min(position_m[2]);
            high_m = high_m.max(position_m[2]);
        }
        assert!(low_m >= -0.5 * EjectionRodGenerator.shaft_length_m() - 1.0e-10);
        assert!(high_m <= expected_m + 1.0e-10);
        assert!(high_m - low_m > 0.8 * expected_m);
    }

    #[test]
    fn the_tip_port_points_along_positive_z() {
        let port = EjectionRodGenerator.tip_port();
        assert_eq!(port.axis_m, [0.0, 0.0, 1.0]);
        assert_eq!(port.dof, Dof::Fixed);
        assert!(port.origin_m[2] > 0.0);
    }

    #[test]
    fn the_cam_port_points_along_negative_z() {
        let port = EjectionRodGenerator.cam_port();
        assert_eq!(port.axis_m, [0.0, 0.0, -1.0]);
        assert_eq!(port.dof, Dof::Prismatic);
        assert!(port.origin_m[2] < 0.0);
    }

    #[test]
    fn a_tip_as_wide_as_the_shaft_is_refused() {
        let mut parameters = ParameterSet::new();
        parameters.set("tip_radius_m", 6.5e-10);
        assert!(EjectionRodGenerator.generate(&parameters).is_err());
    }

    #[test]
    fn every_carbon_keeps_a_valid_valence() {
        let part = rod();
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
