//! The one-sided cam hub.
//!
//! The cam plate of the earlier design carried a closed groove, so it pushed a
//! rod out and pulled it back. A spring now gives the return. The cam is
//! therefore one-sided: a stationary hub below the rotor whose outer surface
//! rises from `base_radius_m` to `base_radius_m + rise_m` at the azimuth
//! `angle_rad`. The rise is a raised cosine over the half-angle
//! `ramp_half_angle_rad`, so the hub dwells over most of the turn and lifts a
//! rod only as its pocket meets the housing outlet. Each follower pin rests on
//! the hub surface, and a leaf spring on the rotor holds the pin against it.
//!
//! The hub does not turn, so the surface meets each pin once for each turn of
//! the rotor. The surface is a radial key, not a tuned motion law. An optional
//! base flange below the hub lets the housing hold the hub, so the hub does not
//! float. This generator builds no rotor, rod, pin, spring or drive. All lengths
//! are SI metres and all angles are radians.

use std::f64::consts::PI;

use nanocad_model::{Atom, Element, Part, Topology};
use nanocad_units::Unit;

use crate::diamond_solid;
use crate::error::PartError;
use crate::generator::PartGenerator;
use crate::lattice_fill::fill_solid;
use crate::parameter::{ParameterSet, ParameterSpec};
use crate::port::{Dof, Port};
use crate::shape::{Bounds, Cylinder, Difference, Solid, Union};

/// The parameter specs of [`CamHubGenerator`], in a stable order.
static CAM_HUB_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new(
        "bore_radius_m",
        Some(Unit::Metre),
        1.5e-9,
        0.0,
        5.0e-8,
        false,
        "radius of the central bore that clears the drive shaft in metres",
    ),
    ParameterSpec::new(
        "thickness_m",
        Some(Unit::Metre),
        9.0e-10,
        1.0e-10,
        1.0e-8,
        false,
        "thickness of the hub along z in metres",
    ),
    ParameterSpec::new(
        "base_radius_m",
        Some(Unit::Metre),
        2.0e-9,
        1.0e-10,
        5.0e-8,
        false,
        "outer radius of the hub over the dwell in metres",
    ),
    ParameterSpec::new(
        "rise_m",
        Some(Unit::Metre),
        1.5e-9,
        0.0,
        5.0e-8,
        false,
        "rise of the hub surface above the dwell in metres; the stroke is this",
    ),
    ParameterSpec::new(
        "ramp_half_angle_rad",
        None,
        0.21,
        1.0e-3,
        PI,
        false,
        "half-angle of the lobe ramp in radians; the hub dwells outside it",
    ),
    ParameterSpec::new(
        "angle_rad",
        None,
        PI,
        -PI,
        PI,
        false,
        "azimuth of the lobe in radians; a rod is fully out here",
    ),
    ParameterSpec::new(
        "base_outer_radius_m",
        Some(Unit::Metre),
        0.0,
        0.0,
        5.0e-8,
        false,
        "outer radius of the base flange in metres; zero means no flange",
    ),
    ParameterSpec::new(
        "base_thickness_m",
        Some(Unit::Metre),
        0.0,
        0.0,
        1.0e-8,
        false,
        "thickness of the base flange below the hub in metres",
    ),
];

/// The generator of the one-sided cam hub.
#[derive(Clone, Copy, Debug, Default)]
pub struct CamHubGenerator;

impl CamHubGenerator {
    /// Returns the port at the rotor axis.
    ///
    /// The hub shares the rotor axis and does not turn, so the joint is rigid.
    pub fn axis_port(&self) -> Port {
        let mut port = Port::named("axis");
        port.origin_m = [0.0, 0.0, 0.0];
        port.axis_m = [0.0, 0.0, 1.0];
        port.dof = Dof::Fixed;
        port
    }

    /// Returns the port at the top of the lobe.
    ///
    /// The point lies on the hub surface at the lobe azimuth, which is the
    /// farthest point of the hub from the axis.
    pub fn lobe_port(&self) -> Port {
        let angle_rad = self.angle_rad();
        let (sin_rad, cos_rad) = angle_rad.sin_cos();
        let reach_m = self.base_radius_m() + self.rise_m();
        let mut port = Port::named("lobe");
        port.origin_m = [reach_m * cos_rad, reach_m * sin_rad, 0.0];
        port.axis_m = [cos_rad, sin_rad, 0.0];
        port.dof = Dof::Fixed;
        port
    }

    /// Returns the bore radius, in metres.
    pub fn bore_radius_m(&self) -> f64 {
        self.default_m("bore_radius_m")
    }

    /// Returns the thickness along z, in metres.
    pub fn thickness_m(&self) -> f64 {
        self.default_m("thickness_m")
    }

    /// Returns the dwell radius of the hub surface, in metres.
    pub fn base_radius_m(&self) -> f64 {
        self.default_m("base_radius_m")
    }

    /// Returns the rise of the hub surface at the lobe, in metres.
    pub fn rise_m(&self) -> f64 {
        self.default_m("rise_m")
    }

    /// Returns the half-angle of the lobe ramp, in radians.
    pub fn ramp_half_angle_rad(&self) -> f64 {
        self.default_m("ramp_half_angle_rad")
    }

    /// Returns the azimuth of the lobe, in radians.
    pub fn angle_rad(&self) -> f64 {
        self.default_m("angle_rad")
    }

    /// Returns the outer radius of the base flange, in metres.
    ///
    /// A value of zero means the hub has no flange.
    pub fn base_outer_radius_m(&self) -> f64 {
        self.default_m("base_outer_radius_m")
    }

    /// Returns the thickness of the base flange, in metres.
    pub fn base_thickness_m(&self) -> f64 {
        self.default_m("base_thickness_m")
    }

    /// Returns the hub surface radius at `theta_rad`, in metres.
    ///
    /// `theta_rad` is an absolute azimuth. The surface dwells at the base radius
    /// outside the ramp. Inside the ramp it rises by a raised cosine, so the pin
    /// reaches full lift at the lobe azimuth with no step.
    pub fn surface_radius_m(&self, theta_rad: f64) -> f64 {
        Self::surface_radius_with(
            self.base_radius_m(),
            self.rise_m(),
            self.ramp_half_angle_rad(),
            theta_rad - self.angle_rad(),
        )
    }

    /// Returns the centre radius of a follower pin of radius `pin_radius_m`
    /// that rests on the hub surface at `theta_rad`, in metres.
    pub fn pin_center_radius_m(&self, theta_rad: f64, pin_radius_m: f64) -> f64 {
        self.surface_radius_m(theta_rad) + pin_radius_m
    }

    /// Returns the dwell-and-ramp surface radius for the given parameters.
    ///
    /// `delta_rad` is the angle from the lobe azimuth. It is wrapped to
    /// `(-PI, PI]`, so a caller may pass an unwrapped angle.
    fn surface_radius_with(base_radius_m: f64, rise_m: f64, ramp_rad: f64, delta_rad: f64) -> f64 {
        let delta_rad = (delta_rad + PI).rem_euclid(2.0 * PI) - PI;
        if delta_rad.abs() >= ramp_rad {
            return base_radius_m;
        }
        base_radius_m + rise_m * 0.5 * (1.0 + (PI * delta_rad / ramp_rad).cos())
    }

    fn default_m(&self, name: &str) -> f64 {
        self.resolve(&ParameterSet::new())
            .ok()
            .and_then(|resolved| resolved.get(name))
            .unwrap_or(0.0)
    }
}

/// A disc whose outer radius follows the cam profile.
///
/// The surface test uses the exact profile, not a set of samples, so the hub
/// and the follower share one radius law.
struct PolarDisc {
    bore_radius_m: f64,
    thickness_m: f64,
    base_radius_m: f64,
    rise_m: f64,
    ramp_rad: f64,
    angle_rad: f64,
}

impl Solid for PolarDisc {
    fn contains_m(&self, point_m: [f64; 3]) -> bool {
        if point_m[2].abs() > 0.5 * self.thickness_m {
            return false;
        }
        let radius_m = point_m[0].hypot(point_m[1]);
        if radius_m < self.bore_radius_m {
            return false;
        }
        let theta_rad = point_m[1].atan2(point_m[0]);
        let surface_m = CamHubGenerator::surface_radius_with(
            self.base_radius_m,
            self.rise_m,
            self.ramp_rad,
            theta_rad - self.angle_rad,
        );
        radius_m <= surface_m
    }

    fn bounds_m(&self) -> Bounds {
        let reach_m = self.base_radius_m + self.rise_m;
        Bounds {
            min_m: [-reach_m, -reach_m, -0.5 * self.thickness_m],
            max_m: [reach_m, reach_m, 0.5 * self.thickness_m],
        }
    }

    fn clone_box(&self) -> Box<dyn Solid> {
        Box::new(PolarDisc {
            bore_radius_m: self.bore_radius_m,
            thickness_m: self.thickness_m,
            base_radius_m: self.base_radius_m,
            rise_m: self.rise_m,
            ramp_rad: self.ramp_rad,
            angle_rad: self.angle_rad,
        })
    }
}

impl PartGenerator for CamHubGenerator {
    fn id(&self) -> &'static str {
        "cam_hub"
    }

    fn name(&self) -> &'static str {
        "cam hub"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        CAM_HUB_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        let resolved = self.resolve(parameters)?;
        let bore_radius_m = resolved.require("bore_radius_m")?;
        let thickness_m = resolved.require("thickness_m")?;
        let base_radius_m = resolved.require("base_radius_m")?;
        let rise_m = resolved.require("rise_m")?;
        let ramp_half_angle_rad = resolved.require("ramp_half_angle_rad")?;
        let angle_rad = resolved.require("angle_rad")?;

        let base_outer_radius_m = resolved.require("base_outer_radius_m")?;
        let base_thickness_m = resolved.require("base_thickness_m")?;

        if base_radius_m <= bore_radius_m {
            return Err(PartError::InvalidGeometry(format!(
                "the base radius {base_radius_m} m is at or below the bore radius \
                 {bore_radius_m} m, so the hub has no material"
            )));
        }
        let flanged = base_outer_radius_m > 0.0 && base_thickness_m > 0.0;
        if flanged && base_outer_radius_m <= bore_radius_m {
            return Err(PartError::InvalidGeometry(format!(
                "the base outer radius {base_outer_radius_m} m is at or below the bore radius \
                 {bore_radius_m} m, so the flange has no material"
            )));
        }
        if flanged && base_outer_radius_m < base_radius_m {
            return Err(PartError::InvalidGeometry(format!(
                "the base outer radius {base_outer_radius_m} m is below the hub base radius \
                 {base_radius_m} m, so the flange does not reach the hub"
            )));
        }

        let hub = PolarDisc {
            bore_radius_m,
            thickness_m,
            base_radius_m,
            rise_m,
            ramp_rad: ramp_half_angle_rad,
            angle_rad,
        };
        let solid: Box<dyn Solid> = if flanged {
            let center_z_m = -0.5 * thickness_m - 0.5 * base_thickness_m;
            Box::new(Union {
                a: Box::new(hub),
                b: Box::new(Difference {
                    outer: Box::new(Cylinder {
                        center_m: [0.0, 0.0, center_z_m],
                        radius_m: base_outer_radius_m,
                        height_m: base_thickness_m,
                    }),
                    inner: Box::new(Cylinder {
                        center_m: [0.0, 0.0, center_z_m],
                        radius_m: bore_radius_m,
                        height_m: base_thickness_m,
                    }),
                }),
            })
        } else {
            Box::new(hub)
        };

        let sites_m = fill_solid(solid.as_ref());
        if sites_m.is_empty() {
            return Err(PartError::InvalidGeometry(
                "the hub has no lattice sites".to_string(),
            ));
        }

        let mut topology = Topology::new();
        let mut carbons = Vec::with_capacity(sites_m.len());
        for site_m in sites_m {
            let index = topology.add_atom(Atom::new(Element::CARBON, site_m, 0.0, "C"));
            carbons.push(index);
        }
        diamond_solid::bond_and_cap(&mut topology, &carbons)?;

        Ok(Part::new(self.id(), topology).with_material("diamond"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hub() -> Part {
        CamHubGenerator
            .generate_with_defaults()
            .expect("the cam hub builds")
    }

    fn carbon_radius_m(part: &Part, index: usize) -> Option<f64> {
        if part.topology.element(index)? != Element::CARBON {
            return None;
        }
        let position_m = part.topology.position_m(index)?;
        Some(position_m[0].hypot(position_m[1]))
    }

    #[test]
    fn the_hub_builds_with_atoms() {
        let part = hub();
        assert!(part.atom_count() > 0);
        assert!(part.bond_count() > 0);
    }

    #[test]
    fn the_hub_has_a_central_bore() {
        let part = hub();
        let bore_m = CamHubGenerator.bore_radius_m();
        for index in 0..part.atom_count() {
            if let Some(radius_m) = carbon_radius_m(&part, index) {
                assert!(
                    radius_m >= bore_m - 1.0e-12,
                    "a carbon at {radius_m} m sits inside the bore"
                );
            }
        }
    }

    #[test]
    fn the_surface_follows_the_profile() {
        let base_m = CamHubGenerator.base_radius_m();
        let rise_m = CamHubGenerator.rise_m();
        let ramp_rad = CamHubGenerator.ramp_half_angle_rad();
        let lobe_rad = CamHubGenerator.angle_rad();
        assert!((CamHubGenerator.surface_radius_m(lobe_rad) - (base_m + rise_m)).abs() < 1.0e-18);
        assert!((CamHubGenerator.surface_radius_m(lobe_rad + ramp_rad) - base_m).abs() < 1.0e-18);
        assert!((CamHubGenerator.surface_radius_m(lobe_rad + PI) - base_m).abs() < 1.0e-18);
        let mid_m = CamHubGenerator.surface_radius_m(lobe_rad + 0.5 * ramp_rad);
        assert!(mid_m > base_m && mid_m < base_m + rise_m);
    }

    #[test]
    fn the_pin_center_radius_adds_the_pin_radius() {
        let pin_radius_m = 2.5e-10;
        let lobe_rad = CamHubGenerator.angle_rad();
        let center_m = CamHubGenerator.pin_center_radius_m(lobe_rad, pin_radius_m);
        let surface_m = CamHubGenerator.surface_radius_m(lobe_rad);
        assert!((center_m - surface_m - pin_radius_m).abs() < 1.0e-18);
    }

    #[test]
    fn no_carbon_sits_past_the_lobe() {
        let part = hub();
        let reach_m = CamHubGenerator.base_radius_m() + CamHubGenerator.rise_m();
        for index in 0..part.atom_count() {
            if let Some(radius_m) = carbon_radius_m(&part, index) {
                assert!(
                    radius_m <= reach_m + 1.0e-12,
                    "a carbon at {radius_m} m sits past the lobe"
                );
            }
        }
    }

    #[test]
    fn the_lobe_port_sits_at_the_farthest_point() {
        let port = CamHubGenerator.lobe_port();
        let angle_rad = CamHubGenerator.angle_rad();
        let (sin_rad, cos_rad) = angle_rad.sin_cos();
        let reach_m = CamHubGenerator.base_radius_m() + CamHubGenerator.rise_m();
        assert!((port.origin_m[0] - reach_m * cos_rad).abs() < 1.0e-18);
        assert!((port.origin_m[1] - reach_m * sin_rad).abs() < 1.0e-18);
        assert_eq!(port.dof, Dof::Fixed);
    }

    #[test]
    fn the_axis_port_points_along_positive_z() {
        let port = CamHubGenerator.axis_port();
        assert_eq!(port.axis_m, [0.0, 0.0, 1.0]);
        assert_eq!(port.dof, Dof::Fixed);
    }

    #[test]
    fn a_bore_wider_than_the_base_is_refused() {
        let mut parameters = ParameterSet::new();
        parameters.set("bore_radius_m", 3.0e-9);
        assert!(CamHubGenerator.generate(&parameters).is_err());
    }

    #[test]
    fn the_flange_reaches_below_the_hub() {
        let mut parameters = ParameterSet::new();
        parameters.set("base_outer_radius_m", 8.7e-9);
        parameters.set("base_thickness_m", 1.0e-9);
        let part = CamHubGenerator
            .generate(&parameters)
            .expect("the flanged hub builds");
        let thickness_m = CamHubGenerator.thickness_m();
        let mut deepest_m: f64 = 0.0;
        for index in 0..part.atom_count() {
            if let Some(position_m) = part.topology.position_m(index) {
                deepest_m = deepest_m.min(position_m[2]);
            }
        }
        assert!(
            deepest_m < -0.5 * thickness_m,
            "the flange did not reach below the hub"
        );
    }

    #[test]
    fn a_flange_that_does_not_reach_the_hub_is_refused() {
        let mut parameters = ParameterSet::new();
        parameters.set("base_outer_radius_m", 1.8e-9);
        parameters.set("base_thickness_m", 1.0e-9);
        assert!(CamHubGenerator.generate(&parameters).is_err());
    }

    #[test]
    fn every_carbon_keeps_a_valid_valence() {
        let part = hub();
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
