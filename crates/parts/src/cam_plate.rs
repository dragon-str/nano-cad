//! The cam plate.
//!
//! The cam ring of the earlier design lay in the rotor plane, so a fixed lobe
//! could not stay clear of the turning rods. This generator fixes that. It
//! builds a stationary plate that lies *below* the rotor, out of its plane. The
//! plate carries one closed groove. The centreline of the groove is a circle of
//! radius `groove_radius_m`, centred at a point offset `groove_eccentricity_m`
//! from the rotor axis, at the azimuth `groove_angle_rad`.
//!
//! Each rod carries a follower pin that reaches down into the groove. The pin
//! must lie on the groove and on its own radial guide, so the groove fixes the
//! rod radius. The groove is profiled, not a circle. Its centreline dwells at
//! `groove_base_radius_m` for most of the turn, and rises to
//! `groove_base_radius_m + groove_rise_m` at the azimuth `groove_angle_rad`.
//! The rise is a raised cosine over the half-angle `groove_ramp_half_angle_rad`,
//! so the rod stays retracted until the ramp and returns on the far side. The
//! two walls of the groove drive the rod in both directions, so no return
//! spring is needed and the motion is reversible. This is the off-centre cam of
//! the reference: the rods "are forcibly ejected by rods thrust outward by the
//! cam surface" (Freitas, *Nanomedicine* Volume I, Section 3.4.2).
//!
//! The ramps hold the stroke to a short arc, so the rod pushes the guest out
//! only as the pocket meets the housing outlet. The rise and the rod length are
//! chosen together so a fully extended rod stops at the rotor rim and does not
//! touch the housing. The groove is a radial key, not a tuned motion law. This
//! generator builds no housing, rod, pin, sliding fit or drive. All lengths are
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

/// The number of follower positions that shape the groove centreline. Each
/// position cuts a disk of the groove width, and the disks overlap into one
/// slot. More samples make the ramp smoother.
pub const GROOVE_SAMPLES: usize = 240;

/// The parameter specs of [`CamPlateGenerator`], in a stable order.
static CAM_PLATE_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new(
        "outer_radius_m",
        Some(Unit::Metre),
        8.7e-9,
        1.0e-9,
        5.0e-8,
        false,
        "outer radius of the cam plate in metres",
    ),
    ParameterSpec::new(
        "inner_radius_m",
        Some(Unit::Metre),
        1.5e-9,
        0.0,
        5.0e-8,
        false,
        "radius of the central hole in the cam plate in metres",
    ),
    ParameterSpec::new(
        "thickness_m",
        Some(Unit::Metre),
        9.0e-10,
        1.0e-10,
        1.0e-8,
        false,
        "thickness of the cam plate along z in metres",
    ),
    ParameterSpec::new(
        "groove_base_radius_m",
        Some(Unit::Metre),
        2.0e-9,
        1.0e-10,
        5.0e-8,
        false,
        "radius of the groove centreline over the dwell in metres",
    ),
    ParameterSpec::new(
        "groove_rise_m",
        Some(Unit::Metre),
        1.5e-9,
        0.0,
        5.0e-8,
        false,
        "rise of the groove centreline above the dwell in metres; the stroke is this",
    ),
    ParameterSpec::new(
        "groove_ramp_half_angle_rad",
        None,
        0.21,
        1.0e-3,
        PI,
        false,
        "half-angle of the lobe ramp in radians; the groove dwells outside it",
    ),
    ParameterSpec::new(
        "groove_width_m",
        Some(Unit::Metre),
        8.0e-10,
        1.0e-11,
        1.0e-8,
        false,
        "full width of the groove across its centreline in metres",
    ),
    ParameterSpec::new(
        "groove_depth_m",
        Some(Unit::Metre),
        6.0e-10,
        1.0e-11,
        1.0e-8,
        false,
        "depth of the groove below the plate face that meets the rotor in metres",
    ),
    ParameterSpec::new(
        "groove_angle_rad",
        None,
        PI,
        -PI,
        PI,
        false,
        "azimuth of the groove centre in radians; the rod is fully out here",
    ),
];

/// The generator of the cam plate.
#[derive(Clone, Copy, Debug, Default)]
pub struct CamPlateGenerator;

impl CamPlateGenerator {
    /// Returns the port at the rotor axis.
    ///
    /// The plate shares the rotor axis and does not turn, so the joint is rigid.
    pub fn axis_port(&self) -> Port {
        Port::named("axis")
    }

    /// Returns the port at the point where a rod is fully extended.
    ///
    /// The point lies on the groove centreline at the groove azimuth, which is
    /// the farthest point of the groove from the rotor axis.
    pub fn groove_port(&self) -> Port {
        let angle_rad = self.groove_angle_rad();
        let (sin_rad, cos_rad) = angle_rad.sin_cos();
        let reach_m = self.groove_base_radius_m() + self.groove_rise_m();
        let mut port = Port::named("groove");
        port.origin_m = [
            reach_m * cos_rad,
            reach_m * sin_rad,
            0.5 * self.thickness_m(),
        ];
        port.axis_m = [cos_rad, sin_rad, 0.0];
        port.dof = Dof::Fixed;
        port
    }

    /// Returns the outer radius of the plate, in metres.
    pub fn outer_radius_m(&self) -> f64 {
        self.default_m("outer_radius_m")
    }

    /// Returns the radius of the central hole, in metres.
    pub fn inner_radius_m(&self) -> f64 {
        self.default_m("inner_radius_m")
    }

    /// Returns the thickness of the plate along z, in metres.
    pub fn thickness_m(&self) -> f64 {
        self.default_m("thickness_m")
    }

    /// Returns the groove centreline radius over the dwell, in metres.
    pub fn groove_base_radius_m(&self) -> f64 {
        self.default_m("groove_base_radius_m")
    }

    /// Returns the groove centreline rise at the lobe, in metres.
    pub fn groove_rise_m(&self) -> f64 {
        self.default_m("groove_rise_m")
    }

    /// Returns the half-angle of the lobe ramp, in radians.
    pub fn groove_ramp_half_angle_rad(&self) -> f64 {
        self.default_m("groove_ramp_half_angle_rad")
    }

    /// Returns the full width of the groove, in metres.
    pub fn groove_width_m(&self) -> f64 {
        self.default_m("groove_width_m")
    }

    /// Returns the depth of the groove, in metres.
    pub fn groove_depth_m(&self) -> f64 {
        self.default_m("groove_depth_m")
    }

    /// Returns the azimuth of the groove centre, in radians.
    pub fn groove_angle_rad(&self) -> f64 {
        self.default_m("groove_angle_rad")
    }

    /// Returns the groove centreline radius at `theta_rad`, in metres.
    ///
    /// `theta_rad` is measured from the groove azimuth. The centreline dwells at
    /// the base radius outside the ramp. Inside the ramp it rises by a raised
    /// cosine, so the rod reaches full extension at the groove azimuth and
    /// returns with no step.
    pub fn rod_radius_m(&self, theta_rad: f64) -> f64 {
        Self::rod_radius_with(
            self.groove_base_radius_m(),
            self.groove_rise_m(),
            self.groove_ramp_half_angle_rad(),
            theta_rad,
        )
    }

    /// Returns the swept slot that the groove removes from the plate face.
    ///
    /// One disk of the groove width is placed at each of [`GROOVE_SAMPLES`]
    /// follower positions around the centreline, and the disks union into one
    /// slot. The disks overlap, so the slot keeps the groove width through the
    /// ramp.
    fn groove_slot_m(
        base_radius_m: f64,
        rise_m: f64,
        ramp_rad: f64,
        groove_width_m: f64,
        groove_depth_m: f64,
        groove_angle_rad: f64,
        height_m: f64,
    ) -> Box<dyn Solid> {
        let half_width_m = 0.5 * groove_width_m;
        let center_z_m = 0.5 * height_m - 0.5 * groove_depth_m;
        let disk = |index: usize| {
            let delta_rad = 2.0 * PI * index as f64 / GROOVE_SAMPLES as f64;
            let profile = Self::rod_radius_with(base_radius_m, rise_m, ramp_rad, delta_rad);
            let angle_rad = groove_angle_rad + delta_rad;
            let (sin_rad, cos_rad) = angle_rad.sin_cos();
            Cylinder {
                center_m: [profile * cos_rad, profile * sin_rad, center_z_m],
                radius_m: half_width_m,
                height_m: groove_depth_m,
            }
        };
        let mut slot: Box<dyn Solid> = Box::new(disk(0));
        for index in 1..GROOVE_SAMPLES {
            slot = Box::new(Union {
                a: slot,
                b: Box::new(disk(index)),
            });
        }
        slot
    }

    /// Returns the dwell-and-ramp centreline radius for the given parameters.
    ///
    /// `delta_rad` is the angle from the groove azimuth. It is wrapped to
    /// `(-PI, PI]`, so a caller may pass an unwrapped sample angle.
    fn rod_radius_with(base_radius_m: f64, rise_m: f64, ramp_rad: f64, delta_rad: f64) -> f64 {
        let delta_rad = (delta_rad + PI).rem_euclid(2.0 * PI) - PI;
        if delta_rad.abs() >= ramp_rad {
            return base_radius_m;
        }
        base_radius_m + rise_m * 0.5 * (1.0 + (PI * delta_rad / ramp_rad).cos())
    }

    #[allow(clippy::too_many_arguments)]
    fn solid_m(
        outer_radius_m: f64,
        inner_radius_m: f64,
        thickness_m: f64,
        base_radius_m: f64,
        rise_m: f64,
        ramp_rad: f64,
        groove_width_m: f64,
        groove_depth_m: f64,
        groove_angle_rad: f64,
    ) -> Box<dyn Solid> {
        let plate = Difference {
            outer: Box::new(Cylinder {
                center_m: [0.0, 0.0, 0.0],
                radius_m: outer_radius_m,
                height_m: thickness_m,
            }),
            inner: Box::new(Cylinder {
                center_m: [0.0, 0.0, 0.0],
                radius_m: inner_radius_m,
                height_m: thickness_m,
            }),
        };
        let groove = Self::groove_slot_m(
            base_radius_m,
            rise_m,
            ramp_rad,
            groove_width_m,
            groove_depth_m,
            groove_angle_rad,
            thickness_m,
        );
        Box::new(Difference {
            outer: Box::new(plate),
            inner: groove,
        })
    }

    fn default_m(&self, name: &str) -> f64 {
        self.resolve(&ParameterSet::new())
            .ok()
            .and_then(|resolved| resolved.get(name))
            .unwrap_or(0.0)
    }
}

impl PartGenerator for CamPlateGenerator {
    fn id(&self) -> &'static str {
        "cam_plate"
    }

    fn name(&self) -> &'static str {
        "cam plate"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        CAM_PLATE_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        let resolved = self.resolve(parameters)?;
        let outer_radius_m = resolved.require("outer_radius_m")?;
        let inner_radius_m = resolved.require("inner_radius_m")?;
        let thickness_m = resolved.require("thickness_m")?;
        let groove_base_radius_m = resolved.require("groove_base_radius_m")?;
        let groove_rise_m = resolved.require("groove_rise_m")?;
        let groove_ramp_half_angle_rad = resolved.require("groove_ramp_half_angle_rad")?;
        let groove_width_m = resolved.require("groove_width_m")?;
        let groove_depth_m = resolved.require("groove_depth_m")?;
        let groove_angle_rad = resolved.require("groove_angle_rad")?;

        if outer_radius_m <= inner_radius_m {
            return Err(PartError::InvalidGeometry(format!(
                "the outer radius {outer_radius_m} m is at or below the hole radius \
                 {inner_radius_m} m, so the plate has no material"
            )));
        }
        if groove_width_m >= 2.0 * groove_base_radius_m {
            return Err(PartError::InvalidGeometry(format!(
                "the groove width {groove_width_m} m is at or above twice the base radius \
                 {groove_base_radius_m} m, so the groove has no centreline"
            )));
        }
        if groove_depth_m > thickness_m {
            return Err(PartError::InvalidGeometry(format!(
                "the groove depth {groove_depth_m} m is deeper than the plate \
                 {thickness_m} m"
            )));
        }
        let reach_m = groove_base_radius_m + groove_rise_m + 0.5 * groove_width_m;
        if reach_m > outer_radius_m {
            return Err(PartError::InvalidGeometry(format!(
                "the groove reaches {reach_m} m and breaks the outer radius \
                 {outer_radius_m} m"
            )));
        }
        let near_m = groove_base_radius_m - 0.5 * groove_width_m;
        if near_m <= inner_radius_m {
            return Err(PartError::InvalidGeometry(format!(
                "the groove comes within {near_m} m of the axis and opens into the \
                 hole radius {inner_radius_m} m"
            )));
        }

        let solid = Self::solid_m(
            outer_radius_m,
            inner_radius_m,
            thickness_m,
            groove_base_radius_m,
            groove_rise_m,
            groove_ramp_half_angle_rad,
            groove_width_m,
            groove_depth_m,
            groove_angle_rad,
        );
        let sites_m = fill_solid(&*solid);
        if sites_m.is_empty() {
            return Err(PartError::InvalidGeometry(
                "the cam plate parameters leave no lattice site inside the plate".to_string(),
            ));
        }

        let mut topology = Topology::new();
        let mut carbons = Vec::with_capacity(sites_m.len());
        for site_m in sites_m {
            carbons.push(topology.add_atom(Atom::new(Element::CARBON, site_m, 0.0, "C")));
        }
        diamond_solid::bond_and_cap(&mut topology, &carbons)?;

        Ok(Part::new("cam_plate", topology).with_material("diamond"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::port::Dof;

    fn plate() -> Part {
        CamPlateGenerator
            .generate_with_defaults()
            .expect("the cam plate builds")
    }

    fn carbon_count(part: &Part) -> usize {
        (0..part.atom_count())
            .filter(|index| part.topology.element(*index) == Some(Element::CARBON))
            .count()
    }

    fn groove_centers_m() -> Vec<[f64; 2]> {
        let angle_rad = CamPlateGenerator.groove_angle_rad();
        (0..GROOVE_SAMPLES)
            .map(|index| {
                let delta_rad = 2.0 * PI * index as f64 / GROOVE_SAMPLES as f64;
                let profile_m = CamPlateGenerator.rod_radius_m(delta_rad);
                let theta_rad = angle_rad + delta_rad;
                let (sin_rad, cos_rad) = theta_rad.sin_cos();
                [profile_m * cos_rad, profile_m * sin_rad]
            })
            .collect()
    }

    fn in_groove_m(position_m: [f64; 3], centers_m: &[[f64; 2]]) -> bool {
        let half_width_m = 0.5 * CamPlateGenerator.groove_width_m();
        let top_m = 0.5 * CamPlateGenerator.thickness_m();
        let bottom_m = top_m - CamPlateGenerator.groove_depth_m();
        if position_m[2] >= top_m - 1.0e-12 || position_m[2] <= bottom_m + 1.0e-12 {
            return false;
        }
        centers_m.iter().any(|center_m| {
            (position_m[0] - center_m[0]).hypot(position_m[1] - center_m[1])
                < half_width_m - 1.0e-12
        })
    }

    #[test]
    fn the_plate_builds_with_atoms() {
        let part = plate();
        assert!(part.atom_count() > 0);
        assert!(part.bond_count() > 0);
    }

    #[test]
    fn the_plate_has_a_central_hole() {
        let part = plate();
        let hole_m = CamPlateGenerator.inner_radius_m();
        for index in 0..part.atom_count() {
            if part.topology.element(index) != Some(Element::CARBON) {
                continue;
            }
            let Some(position_m) = part.topology.position_m(index) else {
                continue;
            };
            let radius_m = position_m[0].hypot(position_m[1]);
            assert!(
                radius_m >= hole_m - 1.0e-10,
                "an atom sits at {radius_m} m inside the hole radius {hole_m} m"
            );
        }
    }

    #[test]
    fn the_groove_is_void() {
        let part = plate();
        let centers_m = groove_centers_m();
        for index in 0..part.atom_count() {
            if part.topology.element(index) != Some(Element::CARBON) {
                continue;
            }
            let Some(position_m) = part.topology.position_m(index) else {
                continue;
            };
            assert!(
                !in_groove_m(position_m, &centers_m),
                "an atom at {position_m:?} m sits inside the groove"
            );
        }
    }

    #[test]
    fn the_groove_removes_material() {
        let mut parameters = ParameterSet::new();
        parameters.set("groove_depth_m", 1.0e-11);
        let shallow = CamPlateGenerator
            .generate(&parameters)
            .expect("the shallow plate builds");
        assert!(
            carbon_count(&shallow) > carbon_count(&plate()),
            "the groove removes no carbon: {} shallow, {} default",
            carbon_count(&shallow),
            carbon_count(&plate())
        );
    }

    #[test]
    fn the_groove_dwells_and_rises_at_the_lobe() {
        let base_m = CamPlateGenerator.groove_base_radius_m();
        let rise_m = CamPlateGenerator.groove_rise_m();
        let ramp_rad = CamPlateGenerator.groove_ramp_half_angle_rad();
        assert!((CamPlateGenerator.rod_radius_m(0.0) - (base_m + rise_m)).abs() < 1.0e-18);
        assert!((CamPlateGenerator.rod_radius_m(ramp_rad) - base_m).abs() < 1.0e-18);
        assert!((CamPlateGenerator.rod_radius_m(PI) - base_m).abs() < 1.0e-18);
        let mid_m = CamPlateGenerator.rod_radius_m(0.5 * ramp_rad);
        assert!(mid_m > base_m && mid_m < base_m + rise_m);
    }

    #[test]
    fn the_groove_port_sits_at_the_farthest_point() {
        let port = CamPlateGenerator.groove_port();
        let angle_rad = CamPlateGenerator.groove_angle_rad();
        let (sin_rad, cos_rad) = angle_rad.sin_cos();
        let reach_m = CamPlateGenerator.groove_base_radius_m() + CamPlateGenerator.groove_rise_m();
        assert!((port.origin_m[0] - reach_m * cos_rad).abs() < 1.0e-18);
        assert!((port.origin_m[1] - reach_m * sin_rad).abs() < 1.0e-18);
        assert_eq!(port.dof, Dof::Fixed);
    }

    #[test]
    fn the_axis_port_points_along_positive_z() {
        let port = CamPlateGenerator.axis_port();
        assert_eq!(port.name, "axis");
        assert_eq!(port.axis_m, [0.0, 0.0, 1.0]);
        assert_eq!(port.dof, Dof::Fixed);
    }

    #[test]
    fn a_groove_that_breaks_the_rim_is_refused() {
        let mut parameters = ParameterSet::new();
        parameters.set("groove_rise_m", 10.0e-9);
        assert!(CamPlateGenerator.generate(&parameters).is_err());
    }

    #[test]
    fn a_groove_that_opens_into_the_hole_is_refused() {
        let mut parameters = ParameterSet::new();
        parameters.set("groove_base_radius_m", 1.0e-9);
        assert!(CamPlateGenerator.generate(&parameters).is_err());
    }

    #[test]
    fn a_groove_wider_than_its_radius_is_refused() {
        let mut parameters = ParameterSet::new();
        parameters.set("groove_width_m", 4.0e-9);
        assert!(CamPlateGenerator.generate(&parameters).is_err());
    }

    #[test]
    fn every_carbon_keeps_a_valid_valence() {
        let part = plate();
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
