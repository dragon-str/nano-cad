//! The cam ring.
//!
//! The sorting rotor in Freitas, *Nanomedicine* Volume I, Section 3.4.2, holds
//! its bound molecule in a pocket until the pocket reaches the inner chamber.
//! There the molecule is forced out: the reference says the bound molecules
//! "are forcibly ejected by rods thrust outward by the cam surface".
//!
//! This generator builds the cam ring. It is a flat annulus that shares the
//! rotor axis and lies in the rotor plane. Its outer edge is the cam surface:
//! the base radius carries each rod retracted, and one lobe stands above the
//! base so that the rod which reaches the lobe is thrust outward. The rotor
//! turns, so a fixed lobe thrusts each rod once for each turn. The ring is cut
//! from the diamond cubic lattice.
//!
//! This generator does not build the housing, the rods, the sliding fit, or a
//! return spring. The lobe is a radial key, not a tuned motion law, so the
//! acceleration of the rod is not designed. The profile is stated here; the
//! force on the rod is a measurement, not geometry. All lengths are SI metres
//! and all angles are radians.

use std::f64::consts::PI;

use nanocad_model::{Atom, Element, Part, Topology};
use nanocad_units::Unit;

use crate::diamond_solid;
use crate::error::PartError;
use crate::generator::PartGenerator;
use crate::lattice_fill::fill_solid;
use crate::parameter::{ParameterSet, ParameterSpec};
use crate::port::{Dof, Port};
use crate::shape::{Box3, Cylinder, Difference, Placed, Solid, Union};

/// The parameter specs of [`CamRingGenerator`], in a stable order.
static CAM_RING_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new(
        "inner_radius_m",
        Some(Unit::Metre),
        5.0e-10,
        0.0,
        5.0e-9,
        false,
        "radius of the central hole in the cam ring in metres",
    ),
    ParameterSpec::new(
        "base_radius_m",
        Some(Unit::Metre),
        1.5e-9,
        1.0e-10,
        1.0e-8,
        false,
        "base radius of the cam ring, where a retracted rod rests, in metres",
    ),
    ParameterSpec::new(
        "lobe_height_m",
        Some(Unit::Metre),
        2.0e-9,
        0.0,
        1.0e-8,
        false,
        "radial height of the lobe above the base radius in metres, or zero for no lobe",
    ),
    ParameterSpec::new(
        "lobe_width_m",
        Some(Unit::Metre),
        3.0e-10,
        1.0e-11,
        5.0e-9,
        false,
        "full tangential width of the lobe in metres",
    ),
    ParameterSpec::new(
        "lobe_angle_rad",
        None,
        PI,
        -PI,
        PI,
        false,
        "azimuth of the lobe centre in radians",
    ),
    ParameterSpec::new(
        "thickness_m",
        Some(Unit::Metre),
        2.0e-9,
        1.0e-10,
        1.0e-8,
        false,
        "thickness of the cam ring along z in metres",
    ),
];

/// The generator of the cam ring.
#[derive(Clone, Copy, Debug, Default)]
pub struct CamRingGenerator;

impl CamRingGenerator {
    /// Returns the port at the rotor axis.
    ///
    /// The ring shares the rotor axis and does not turn, so the joint is rigid.
    pub fn axis_port(&self) -> Port {
        Port::named("axis")
    }

    /// Returns the port at the peak of the cam lobe.
    ///
    /// The port sits on the cam surface above the base and points outward along
    /// the lobe azimuth. The rod rests here while the lobe carries it.
    pub fn cam_port(&self) -> Port {
        let angle_rad = self.lobe_angle_rad();
        let (sin_rad, cos_rad) = angle_rad.sin_cos();
        let peak_m = self.base_radius_m() + self.lobe_height_m();
        let mut port = Port::named("cam");
        port.origin_m = [peak_m * cos_rad, peak_m * sin_rad, 0.0];
        port.axis_m = [cos_rad, sin_rad, 0.0];
        port.dof = Dof::Fixed;
        port
    }

    /// Returns the base radius of the ring, in metres.
    pub fn base_radius_m(&self) -> f64 {
        self.default_m("base_radius_m")
    }

    /// Returns the height of the lobe above the base, in metres.
    pub fn lobe_height_m(&self) -> f64 {
        self.default_m("lobe_height_m")
    }

    /// Returns the azimuth of the lobe centre, in radians.
    pub fn lobe_angle_rad(&self) -> f64 {
        self.default_m("lobe_angle_rad")
    }

    /// Returns the thickness of the ring along z, in metres.
    pub fn thickness_m(&self) -> f64 {
        self.default_m("thickness_m")
    }

    /// Returns the full tangential width of the lobe, in metres.
    pub fn lobe_width_m(&self) -> f64 {
        self.default_m("lobe_width_m")
    }

    /// Returns the radius of the central hole in the ring, in metres.
    pub fn inner_radius_m(&self) -> f64 {
        self.default_m("inner_radius_m")
    }

    #[allow(clippy::too_many_arguments)]
    fn solid_m(
        inner_radius_m: f64,
        base_radius_m: f64,
        lobe_height_m: f64,
        lobe_width_m: f64,
        thickness_m: f64,
        lobe_angle_rad: f64,
    ) -> Box<dyn Solid> {
        let inner = Cylinder {
            center_m: [0.0, 0.0, 0.0],
            radius_m: inner_radius_m,
            height_m: thickness_m,
        };
        let outer = Cylinder {
            center_m: [0.0, 0.0, 0.0],
            radius_m: base_radius_m,
            height_m: thickness_m,
        };
        let ring = Difference {
            outer: Box::new(outer),
            inner: Box::new(inner),
        };
        let mut solid: Box<dyn Solid> = Box::new(ring);

        if lobe_height_m > 0.0 {
            let center_radius_m = base_radius_m + 0.5 * (lobe_height_m - lobe_width_m);
            let lobe = Placed {
                inner: Box::new(Box3::from_center_m(
                    [center_radius_m, 0.0, 0.0],
                    [lobe_height_m + lobe_width_m, lobe_width_m, thickness_m],
                )),
                rotation_rad: lobe_angle_rad,
                offset_m: [0.0, 0.0, 0.0],
            };
            solid = Box::new(Union {
                a: solid,
                b: Box::new(lobe),
            });
        }

        solid
    }

    fn default_m(&self, name: &str) -> f64 {
        self.resolve(&ParameterSet::new())
            .ok()
            .and_then(|resolved| resolved.get(name))
            .unwrap_or(0.0)
    }
}

impl PartGenerator for CamRingGenerator {
    fn id(&self) -> &'static str {
        "cam_ring"
    }

    fn name(&self) -> &'static str {
        "cam ring"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        CAM_RING_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        let resolved = self.resolve(parameters)?;
        let inner_radius_m = resolved.require("inner_radius_m")?;
        let base_radius_m = resolved.require("base_radius_m")?;
        let lobe_height_m = resolved.require("lobe_height_m")?;
        let lobe_width_m = resolved.require("lobe_width_m")?;
        let lobe_angle_rad = resolved.require("lobe_angle_rad")?;
        let thickness_m = resolved.require("thickness_m")?;

        if base_radius_m <= inner_radius_m {
            return Err(PartError::InvalidGeometry(format!(
                "the base radius {base_radius_m} m is at or below the hole radius \
                 {inner_radius_m} m, so the ring has no material"
            )));
        }
        if lobe_height_m > 0.0 && lobe_width_m <= 0.0 {
            return Err(PartError::InvalidGeometry(
                "the lobe has a height but no width, so it cannot carry a rod".to_string(),
            ));
        }

        let solid = Self::solid_m(
            inner_radius_m,
            base_radius_m,
            lobe_height_m,
            lobe_width_m,
            thickness_m,
            lobe_angle_rad,
        );
        let sites_m = fill_solid(&*solid);
        if sites_m.is_empty() {
            return Err(PartError::InvalidGeometry(
                "the cam ring parameters leave no lattice site inside the ring".to_string(),
            ));
        }

        let mut topology = Topology::new();
        let mut carbons = Vec::with_capacity(sites_m.len());
        for site_m in sites_m {
            carbons.push(topology.add_atom(Atom::new(Element::CARBON, site_m, 0.0, "C")));
        }
        diamond_solid::bond_and_cap(&mut topology, &carbons)?;

        Ok(Part::new("cam_ring", topology).with_material("diamond"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cam() -> Part {
        CamRingGenerator
            .generate_with_defaults()
            .expect("the cam ring builds")
    }

    fn carbon_count(part: &Part) -> usize {
        (0..part.atom_count())
            .filter(|index| part.topology.element(*index) == Some(Element::CARBON))
            .count()
    }

    fn carbon_radii_m(part: &Part) -> Vec<f64> {
        (0..part.atom_count())
            .filter(|index| part.topology.element(*index) == Some(Element::CARBON))
            .filter_map(|index| part.topology.position_m(index))
            .map(|position_m| position_m[0].hypot(position_m[1]))
            .collect()
    }

    #[test]
    fn the_cam_ring_builds_with_atoms() {
        let part = cam();
        assert!(part.atom_count() > 0);
        assert!(part.bond_count() > 0);
    }

    #[test]
    fn the_cam_ring_has_a_hole() {
        let part = cam();
        let hole_m = CamRingGenerator.inner_radius_m();
        for radius_m in carbon_radii_m(&part) {
            assert!(
                radius_m >= hole_m - 1.0e-10,
                "an atom sits at {radius_m} m inside the hole radius {hole_m} m"
            );
        }
    }

    #[test]
    fn the_lobe_stands_above_the_base() {
        let part = cam();
        let base_m = CamRingGenerator.base_radius_m();
        let height_m = CamRingGenerator.lobe_height_m();
        let peak_m = base_m + height_m;
        let highest_m = carbon_radii_m(&part).into_iter().fold(0.0_f64, f64::max);
        assert!(
            highest_m > base_m + 0.5 * height_m,
            "the highest atom sits at {highest_m} m, not well above the base {base_m} m"
        );
        assert!(highest_m <= peak_m + 1.0e-10);
    }

    #[test]
    fn every_atom_above_the_base_lies_in_the_lobe() {
        let part = cam();
        let base_m = CamRingGenerator.base_radius_m();
        let angle_rad = CamRingGenerator.lobe_angle_rad();
        let (sin_rad, cos_rad) = angle_rad.sin_cos();
        let half_width_m = 0.5 * CamRingGenerator.lobe_width_m();
        let tolerance_m = diamond_solid::diamond_bond_length_m();
        for index in 0..part.atom_count() {
            if part.topology.element(index) != Some(Element::CARBON) {
                continue;
            }
            let Some(position_m) = part.topology.position_m(index) else {
                continue;
            };
            let radius_m = position_m[0].hypot(position_m[1]);
            if radius_m <= base_m {
                continue;
            }
            let across_m = -position_m[0] * sin_rad + position_m[1] * cos_rad;
            assert!(
                across_m.abs() <= half_width_m + tolerance_m,
                "an atom at radius {radius_m} m stands {across_m} m across the lobe axis"
            );
        }
    }

    #[test]
    fn a_lobe_of_zero_height_removes_nothing() {
        let with_lobe = cam();
        let mut parameters = ParameterSet::new();
        parameters.set("lobe_height_m", 0.0);
        let without_lobe = CamRingGenerator
            .generate(&parameters)
            .expect("the ring builds");
        assert!(without_lobe.atom_count() > 0);
        assert!(
            carbon_count(&without_lobe) < carbon_count(&with_lobe),
            "the lobe adds no carbon: {} without, {} with",
            carbon_count(&without_lobe),
            carbon_count(&with_lobe)
        );
    }

    #[test]
    fn the_cam_port_points_outward_along_the_lobe() {
        let port = CamRingGenerator.cam_port();
        let angle_rad = CamRingGenerator.lobe_angle_rad();
        let (sin_rad, cos_rad) = angle_rad.sin_cos();
        assert!((port.axis_m[0] - cos_rad).abs() < 1.0e-12);
        assert!((port.axis_m[1] - sin_rad).abs() < 1.0e-12);
        assert_eq!(port.axis_m[2], 0.0);
        assert_eq!(port.dof, Dof::Fixed);
    }

    #[test]
    fn the_axis_port_points_along_positive_z() {
        let port = CamRingGenerator.axis_port();
        assert_eq!(port.name, "axis");
        assert_eq!(port.axis_m, [0.0, 0.0, 1.0]);
        assert_eq!(port.dof, Dof::Fixed);
    }

    #[test]
    fn a_base_radius_at_or_below_the_hole_is_refused() {
        let mut parameters = ParameterSet::new();
        parameters.set("inner_radius_m", 1.5e-9);
        parameters.set("base_radius_m", 1.5e-9);
        assert!(CamRingGenerator.generate(&parameters).is_err());
    }

    #[test]
    fn every_carbon_keeps_a_valid_valence() {
        let part = cam();
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
