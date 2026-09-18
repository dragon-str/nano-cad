//! The return leaf spring.
//!
//! A cam pushes the ejection rod out. Something must push the rod back in, or
//! the rod stays out and its pocket is not swept again. Freitas, *Nanomedicine*
//! Volume I, Section 3.4.2, returns the rods with a spring. This generator
//! builds that spring as a diamond leaf: a thin cantilever that bends when the
//! rod presses on its free tip.
//!
//! The leaf hangs along `-z` from its anchor at the origin. The width runs
//! along `x`, the thin direction runs along `y`, and the length runs along `z`.
//! The load is along `y`, the thin direction, so the leaf is compliant in the
//! radial direction and stiff in the other two. The rotor carries the anchor
//! and the rod forces the free tip, so the leaf turns with the rotor.
//!
//! The stiffness is an estimate from the Euler-Bernoulli cantilever formula,
//! `k = 3*E*I/L^3` with `I = w*t^3/12`. It is an order-of-magnitude number, not
//! a measured value. The leaf is a schematic: the diamond lattice cannot
//! resolve a feature much below 2e-10 m, so this part is a few tens of atoms and
//! it is not a machined spring. This generator does not build the anchor boss,
//! the rod or the cam, and it does not model the contact at the tip or the
//! bond strain under load. All lengths are SI metres and all angles are radians.

use nanocad_model::{Atom, Element, Part, Topology};
use nanocad_units::Unit;

use crate::diamond_solid;
use crate::error::PartError;
use crate::generator::PartGenerator;
use crate::lattice_fill::fill_solid;
use crate::parameter::{ParameterSet, ParameterSpec};
use crate::port::Port;
use crate::shape::Box3;

/// The Young's modulus of diamond, in pascals.
///
/// The value is a typical polycrystalline figure. It sets the stiffness
/// estimate only.
pub const DIAMOND_YOUNG_MODULUS_PA: f64 = 1.05e12;

/// The parameter specs of [`LeafSpringGenerator`], in a stable order.
static LEAF_SPRING_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new(
        "leaf_length_m",
        Some(Unit::Metre),
        1.6e-9,
        2.0e-10,
        1.0e-8,
        false,
        "length of the leaf from the anchor to the free tip in metres",
    ),
    ParameterSpec::new(
        "leaf_width_m",
        Some(Unit::Metre),
        6.0e-10,
        1.0e-10,
        5.0e-9,
        false,
        "width of the leaf across the load direction in metres",
    ),
    ParameterSpec::new(
        "leaf_thickness_m",
        Some(Unit::Metre),
        3.6e-10,
        1.0e-10,
        2.0e-9,
        false,
        "thickness of the leaf along the load direction in metres",
    ),
];

/// The generator of the return leaf spring.
#[derive(Clone, Copy, Debug, Default)]
pub struct LeafSpringGenerator;

impl LeafSpringGenerator {
    /// Returns the port where the leaf meets the rotor.
    ///
    /// The port sits at the top of the leaf, its axis points out along `+z`,
    /// and the joint is rigid.
    pub fn anchor_port(&self) -> Port {
        Port::named("anchor")
    }

    /// Returns the port at the free tip of the leaf.
    ///
    /// The port sits at the bottom of the leaf, its axis points out along `+y`
    /// toward the rod, and the joint is a contact, not a bond.
    pub fn load_port(&self) -> Port {
        let mut port = Port::named("load");
        port.origin_m = [0.0, 0.0, -self.length_m()];
        port.axis_m = [0.0, 1.0, 0.0];
        port
    }

    /// Returns the length of the leaf, in metres.
    pub fn length_m(&self) -> f64 {
        self.default_m("leaf_length_m")
    }

    /// Returns the width of the leaf, in metres.
    pub fn width_m(&self) -> f64 {
        self.default_m("leaf_width_m")
    }

    /// Returns the thickness of the leaf, in metres.
    pub fn thickness_m(&self) -> f64 {
        self.default_m("leaf_thickness_m")
    }

    /// Returns the second moment of area about the bending axis, in metres^4.
    ///
    /// The load is along the thin direction, so the width is the other
    /// cross-section side.
    pub fn second_moment_m4(&self) -> f64 {
        self.width_m() * self.thickness_m().powi(3) / 12.0
    }

    /// Returns the estimated tip stiffness, in newtons per metre.
    ///
    /// This is the Euler-Bernoulli cantilever value, an estimate.
    pub fn tip_stiffness_n_per_m(&self) -> f64 {
        3.0 * DIAMOND_YOUNG_MODULUS_PA * self.second_moment_m4() / self.length_m().powi(3)
    }

    /// Returns the estimated tip load for a deflection, in newtons.
    pub fn tip_load_n(&self, deflection_m: f64) -> f64 {
        self.tip_stiffness_n_per_m() * deflection_m
    }

    fn default_m(&self, name: &str) -> f64 {
        self.resolve(&ParameterSet::new())
            .ok()
            .and_then(|resolved| resolved.get(name))
            .unwrap_or(0.0)
    }
}

impl PartGenerator for LeafSpringGenerator {
    fn id(&self) -> &'static str {
        "leaf_spring"
    }

    fn name(&self) -> &'static str {
        "leaf spring"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        LEAF_SPRING_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        let resolved = self.resolve(parameters)?;
        let leaf_length_m = resolved.require("leaf_length_m")?;
        let leaf_width_m = resolved.require("leaf_width_m")?;
        let leaf_thickness_m = resolved.require("leaf_thickness_m")?;

        if leaf_thickness_m >= leaf_width_m {
            return Err(PartError::InvalidGeometry(format!(
                "the leaf thickness {leaf_thickness_m} m is at or above the width \
                 {leaf_width_m} m, so the leaf is a block and not a beam"
            )));
        }

        let solid = Box3::from_center_m(
            [0.0, 0.0, -0.5 * leaf_length_m],
            [leaf_width_m, leaf_thickness_m, leaf_length_m],
        );

        let sites_m = fill_solid(&solid);
        if sites_m.is_empty() {
            return Err(PartError::InvalidGeometry(
                "the leaf parameters leave no lattice site inside the leaf".to_string(),
            ));
        }

        let mut topology = Topology::new();
        let mut carbons = Vec::with_capacity(sites_m.len());
        for site_m in sites_m {
            carbons.push(topology.add_atom(Atom::new(Element::CARBON, site_m, 0.0, "C")));
        }
        diamond_solid::bond_and_cap(&mut topology, &carbons)?;

        Ok(Part::new("leaf_spring", topology).with_material("diamond"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::port::Dof;

    fn spring() -> Part {
        LeafSpringGenerator
            .generate_with_defaults()
            .expect("the spring builds")
    }

    #[test]
    fn the_spring_builds_with_atoms() {
        let part = spring();
        assert!(part.atom_count() > 0);
        assert!(part.bond_count() > 0);
    }

    #[test]
    fn the_spring_hangs_below_its_anchor() {
        let part = spring();
        let expected_m = LeafSpringGenerator.length_m();
        let mut low_m = f64::INFINITY;
        let mut high_m = f64::NEG_INFINITY;
        for index in 0..part.atom_count() {
            let Some(position_m) = part.topology.position_m(index) else {
                continue;
            };
            low_m = low_m.min(position_m[2]);
            high_m = high_m.max(position_m[2]);
        }
        let lattice_step_m = crate::diamond_solid::diamond_bond_length_m();
        assert!(high_m <= lattice_step_m, "the leaf rises above its anchor");
        assert!(low_m >= -expected_m - lattice_step_m);
        assert!(high_m - low_m > 0.8 * expected_m);
    }

    #[test]
    fn the_leaf_is_thin_along_the_load_axis() {
        let part = spring();
        let mut widest_y_m: f64 = 0.0;
        let mut widest_x_m: f64 = 0.0;
        for index in 0..part.atom_count() {
            let Some(position_m) = part.topology.position_m(index) else {
                continue;
            };
            widest_y_m = widest_y_m.max(position_m[1].abs());
            widest_x_m = widest_x_m.max(position_m[0].abs());
        }
        let lattice_step_m = crate::diamond_solid::diamond_bond_length_m();
        assert!(widest_y_m <= 0.5 * LeafSpringGenerator.thickness_m() + lattice_step_m);
        assert!(widest_x_m <= 0.5 * LeafSpringGenerator.width_m() + lattice_step_m);
    }

    #[test]
    fn the_anchor_port_sits_on_the_top_face() {
        let port = LeafSpringGenerator.anchor_port();
        assert_eq!(port.origin_m, [0.0, 0.0, 0.0]);
        assert_eq!(port.axis_m, [0.0, 0.0, 1.0]);
        assert_eq!(port.dof, Dof::Fixed);
    }

    #[test]
    fn the_load_port_points_along_positive_y() {
        let port = LeafSpringGenerator.load_port();
        assert_eq!(port.axis_m, [0.0, 1.0, 0.0]);
        assert_eq!(port.dof, Dof::Fixed);
        assert!((port.origin_m[2] + LeafSpringGenerator.length_m()).abs() < 1.0e-18);
    }

    #[test]
    fn the_stiffness_estimate_is_positive_and_thickness_driven() {
        let stiffness = LeafSpringGenerator.tip_stiffness_n_per_m();
        assert!(stiffness.is_finite() && stiffness > 0.0);
        let mut parameters = ParameterSet::new();
        parameters.set("leaf_thickness_m", 2.0 * LeafSpringGenerator.thickness_m());
        let thicker = LeafSpringGenerator.resolve(&parameters).expect("resolve");
        let thick_stiffness = 3.0
            * DIAMOND_YOUNG_MODULUS_PA
            * (thicker.require("leaf_width_m").expect("width")
                * thicker
                    .require("leaf_thickness_m")
                    .expect("thickness")
                    .powi(3)
                / 12.0)
            / thicker.require("leaf_length_m").expect("length").powi(3);
        assert!(thick_stiffness > 7.0 * stiffness);
        assert!((LeafSpringGenerator.tip_load_n(1.5e-9) - 1.5e-9 * stiffness).abs() < 1.0e-30);
    }

    #[test]
    fn a_block_leaf_is_refused() {
        let mut parameters = ParameterSet::new();
        parameters.set("leaf_thickness_m", 6.0e-10);
        parameters.set("leaf_width_m", 6.0e-10);
        assert!(LeafSpringGenerator.generate(&parameters).is_err());
    }

    #[test]
    fn every_carbon_keeps_a_valid_valence() {
        let part = spring();
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
