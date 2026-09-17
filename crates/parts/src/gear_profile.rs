use std::f64::consts::{FRAC_PI_2, PI};

use nanocad_model::{Atom, Element, Part, Topology};
use nanocad_units::Unit;

use crate::error::PartError;
use crate::generator::PartGenerator;
use crate::parameter::{ParameterSet, ParameterSpec};

/// The full-depth addendum coefficient `h_a*`. The addendum is `h_a* m`.
///
/// Source: ISO 21771 / AGMA 908-B89 standard full-depth tooth form.
const ADDENDUM_COEFFICIENT: f64 = 0.8;

/// The full-depth dedendum coefficient `h_f*`. The dedendum is `h_f* m`.
///
/// Source: ISO 21771 standard full-depth tooth form.
const DEDENDUM_COEFFICIENT: f64 = 1.25;

/// The undercut factor. The minimum tooth count without undercut is
/// `2 h_a* / sin^2(alpha)`.
///
/// Source: J. E. Shigley, Mechanical Engineering Design, involute undercut.
const UNDERCUT_FACTOR: f64 = 2.0;

/// The default pressure angle, 20 degrees, in radians.
///
/// Source: ISO 53 standard basic rack tooth profile.
const DEFAULT_PRESSURE_ANGLE_RAD: f64 = 20.0 * PI / 180.0;

/// The involute function `inv(alpha) = tan(alpha) - alpha`.
fn involute_function(alpha_rad: f64) -> f64 {
    alpha_rad.tan() - alpha_rad
}

/// The 2D geometry of a standard full-depth involute spur gear.
///
/// The profile is driven by the module `m`, the tooth count `z`, and the
/// pressure angle `alpha`. All lengths are SI metres. The involute flank is
/// placed so that the tooth is `pi m / 2` thick on the pitch circle.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GearProfile {
    module_m: f64,
    teeth: usize,
    pressure_angle_rad: f64,
    backlash_m: f64,
}

impl GearProfile {
    /// Builds a gear profile and rejects non-physical inputs.
    pub fn new(module_m: f64, teeth: usize, pressure_angle_rad: f64) -> Result<Self, PartError> {
        if !module_m.is_finite() || module_m <= 0.0 {
            return Err(PartError::InvalidGeometry(format!(
                "module {module_m} m must be finite and positive"
            )));
        }
        if teeth < 3 {
            return Err(PartError::InvalidGeometry(format!(
                "tooth count {teeth} must be at least 3"
            )));
        }
        if !pressure_angle_rad.is_finite()
            || pressure_angle_rad <= 0.0
            || pressure_angle_rad >= FRAC_PI_2
        {
            return Err(PartError::InvalidGeometry(format!(
                "pressure angle {pressure_angle_rad} rad must be in (0, pi/2)"
            )));
        }
        Ok(Self {
            module_m,
            teeth,
            pressure_angle_rad,
            backlash_m: 0.0,
        })
    }

    /// Returns the profile with a circumferential backlash, in metres.
    ///
    /// The tooth is thinner by `backlash_m`. Two mating gears that each use
    /// this profile then clear each other by twice `backlash_m` at the flanks.
    pub fn with_backlash(mut self, backlash_m: f64) -> Self {
        self.backlash_m = backlash_m.max(0.0);
        self
    }

    /// The circumferential backlash, in metres.
    pub fn backlash_m(&self) -> f64 {
        self.backlash_m
    }

    /// The angular half-thickness removed from the tooth by the backlash.
    fn backlash_half_angle_rad(&self) -> f64 {
        self.backlash_m / (4.0 * self.pitch_radius_m())
    }

    /// The module, in metres.
    pub fn module_m(&self) -> f64 {
        self.module_m
    }

    /// The tooth count.
    pub fn teeth(&self) -> usize {
        self.teeth
    }

    /// The pressure angle, in radians.
    pub fn pressure_angle_rad(&self) -> f64 {
        self.pressure_angle_rad
    }

    /// The minimum tooth count without undercut, for a pressure angle.
    pub fn minimum_tooth_count(pressure_angle_rad: f64) -> usize {
        let sine = pressure_angle_rad.sin();
        let raw = UNDERCUT_FACTOR / (sine * sine);
        if !raw.is_finite() || raw <= 0.0 {
            return usize::MAX;
        }
        raw.ceil() as usize
    }

    /// Rejects a tooth count below the undercut limit.
    pub fn check_undercut(&self) -> Result<(), PartError> {
        let min_teeth = Self::minimum_tooth_count(self.pressure_angle_rad);
        if self.teeth < min_teeth {
            return Err(PartError::ToothCountBelowUndercut {
                teeth: self.teeth,
                min_teeth,
            });
        }
        Ok(())
    }

    /// The pitch radius `d / 2 = m z / 2`, in metres.
    pub fn pitch_radius_m(&self) -> f64 {
        self.module_m * self.teeth as f64 / 2.0
    }

    /// The pitch diameter `d = m z`, in metres.
    pub fn pitch_diameter_m(&self) -> f64 {
        self.module_m * self.teeth as f64
    }

    /// The base radius `d_b / 2 = (d / 2) cos(alpha)`, in metres.
    pub fn base_radius_m(&self) -> f64 {
        self.pitch_radius_m() * self.pressure_angle_rad.cos()
    }

    /// The addendum `a = h_a* m`, in metres.
    pub fn addendum_m(&self) -> f64 {
        ADDENDUM_COEFFICIENT * self.module_m
    }

    /// The dedendum `b = h_f* m`, in metres.
    pub fn dedendum_m(&self) -> f64 {
        DEDENDUM_COEFFICIENT * self.module_m
    }

    /// The outside (addendum) radius `r_a = d / 2 + a`, in metres.
    pub fn outer_radius_m(&self) -> f64 {
        self.pitch_radius_m() + self.addendum_m()
    }

    /// The root (dedendum) radius `r_f = d / 2 - b`, in metres.
    pub fn root_radius_m(&self) -> f64 {
        self.pitch_radius_m() - self.dedendum_m()
    }

    /// The polar angle of one involute flank at a radius.
    ///
    /// `sign` is `-1.0` for the right flank and `+1.0` for the left flank. The
    /// point at the pitch radius sits at `sign pi / (2 z)`.
    pub fn flank_angle_rad(&self, radius_m: f64, sign: f64) -> Result<f64, PartError> {
        let base_m = self.base_radius_m();
        if radius_m < base_m {
            return Err(PartError::InvalidGeometry(format!(
                "radius {radius_m} m is inside the base circle"
            )));
        }
        let cosine = (base_m / radius_m).clamp(-1.0, 1.0);
        let alpha_r = cosine.acos();
        let half_tooth_rad = PI / (2.0 * self.teeth as f64) - self.backlash_half_angle_rad();
        Ok(sign
            * (half_tooth_rad + involute_function(self.pressure_angle_rad)
                - involute_function(alpha_r)))
    }

    /// Builds the closed 2D outline of the gear, in metres.
    ///
    /// The outline runs counter-clockwise through every tooth. The first and
    /// last points are equal, so the curve closes.
    pub fn outline_points(
        &self,
        samples_per_flank: usize,
        samples_per_arc: usize,
    ) -> Result<Vec<[f64; 2]>, PartError> {
        let flank_samples = samples_per_flank.max(1);
        let arc_samples = samples_per_arc.max(1);
        let teeth_f = self.teeth as f64;
        let tooth_pitch_rad = 2.0 * PI / teeth_f;
        let half_pitch_rad = PI / teeth_f;
        let outer_m = self.outer_radius_m();
        let root_m = self.root_radius_m();
        let base_m = self.base_radius_m();
        if root_m <= 0.0 {
            return Err(PartError::InvalidGeometry(format!(
                "root radius {root_m} m must be positive"
            )));
        }

        // The tip arc must join the two flanks. The half-angle at the outer
        // radius is the flank angle there, not the half pitch: the half pitch
        // is twice the tooth half-thickness and would overhang the flanks.
        let tip_half_angle_rad = self.flank_angle_rad(outer_m, 1.0)?;
        if tip_half_angle_rad <= 0.0 {
            return Err(PartError::InvalidGeometry(
                "tooth is pointed at the tip circle".to_owned(),
            ));
        }

        let flank_start_m = root_m.max(base_m);
        let flank_start_angle_rad = self.flank_angle_rad(flank_start_m, -1.0)?;

        let mut local: Vec<[f64; 2]> = Vec::new();
        if root_m < base_m {
            for i in 0..=arc_samples {
                let t = i as f64 / arc_samples as f64;
                let angle_rad = -half_pitch_rad + (half_pitch_rad + flank_start_angle_rad) * t;
                local.push([root_m * angle_rad.cos(), root_m * angle_rad.sin()]);
            }
        }
        local.push([
            flank_start_m * flank_start_angle_rad.cos(),
            flank_start_m * flank_start_angle_rad.sin(),
        ]);
        for i in 1..=flank_samples {
            let t = i as f64 / flank_samples as f64;
            let radius_m = flank_start_m + (outer_m - flank_start_m) * t;
            let angle_rad = self.flank_angle_rad(radius_m, -1.0)?;
            local.push([radius_m * angle_rad.cos(), radius_m * angle_rad.sin()]);
        }
        for i in 0..=arc_samples {
            let t = i as f64 / arc_samples as f64;
            let angle_rad = -tip_half_angle_rad + 2.0 * tip_half_angle_rad * t;
            local.push([outer_m * angle_rad.cos(), outer_m * angle_rad.sin()]);
        }
        for i in (0..flank_samples).rev() {
            let t = i as f64 / flank_samples as f64;
            let radius_m = flank_start_m + (outer_m - flank_start_m) * t;
            let angle_rad = self.flank_angle_rad(radius_m, 1.0)?;
            local.push([radius_m * angle_rad.cos(), radius_m * angle_rad.sin()]);
        }
        if root_m < base_m {
            for i in 0..=arc_samples {
                let t = i as f64 / arc_samples as f64;
                let angle_rad =
                    -flank_start_angle_rad + (half_pitch_rad + flank_start_angle_rad) * t;
                local.push([root_m * angle_rad.cos(), root_m * angle_rad.sin()]);
            }
        }

        let mut points: Vec<[f64; 2]> = Vec::with_capacity(local.len() * self.teeth);
        for tooth in 0..self.teeth {
            let rotation_rad = tooth as f64 * tooth_pitch_rad;
            let (sin_rot, cos_rot) = rotation_rad.sin_cos();
            for point in &local {
                points.push([
                    point[0] * cos_rot - point[1] * sin_rot,
                    point[0] * sin_rot + point[1] * cos_rot,
                ]);
            }
        }

        points.dedup_by(|a, b| (a[0] - b[0]).abs() < 1.0e-15 && (a[1] - b[1]).abs() < 1.0e-15);
        let Some(first) = points.first().copied() else {
            return Err(PartError::InvalidGeometry(
                "profile produced no points".to_owned(),
            ));
        };
        if let Some(last) = points.last().copied() {
            if (last[0] - first[0]).abs() < 1.0e-12 && (last[1] - first[1]).abs() < 1.0e-12 {
                points.pop();
            }
        }
        points.push(first);
        Ok(points)
    }
}

static GEAR_PROFILE_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new(
        "module_m",
        Some(Unit::Metre),
        0.5e-9,
        0.05e-9,
        5.0e-9,
        false,
        "gear module in metres",
    ),
    ParameterSpec::new(
        "tooth_count",
        None,
        20.0,
        3.0,
        200.0,
        true,
        "number of teeth",
    ),
    ParameterSpec::new(
        "pressure_angle_rad",
        None,
        DEFAULT_PRESSURE_ANGLE_RAD,
        10.0 * PI / 180.0,
        25.0 * PI / 180.0,
        false,
        "pressure angle in radians",
    ),
    ParameterSpec::new(
        "samples_per_flank",
        None,
        8.0,
        1.0,
        64.0,
        true,
        "samples along each involute flank",
    ),
    ParameterSpec::new(
        "samples_per_arc",
        None,
        4.0,
        1.0,
        64.0,
        true,
        "samples along each tip and root arc",
    ),
];

/// Builds the 2D involute tooth profile as a point set.
///
/// The result is a [`Part`] whose atoms are the outline points in the `z = 0`
/// plane. The part carries no bonds; it is a geometry carrier for the profile,
/// not a molecule.
#[derive(Clone, Copy, Debug, Default)]
pub struct GearProfileGenerator;

impl GearProfileGenerator {
    /// Resolves the inputs and returns the profile geometry.
    pub fn profile(&self, parameters: &ParameterSet) -> Result<GearProfile, PartError> {
        let resolved = self.resolve(parameters)?;
        let module_m = resolved.require("module_m")?;
        let teeth = resolved.require("tooth_count")? as usize;
        let pressure_angle_rad = resolved.require("pressure_angle_rad")?;
        let profile = GearProfile::new(module_m, teeth, pressure_angle_rad)?;
        profile.check_undercut()?;
        Ok(profile)
    }
}

impl PartGenerator for GearProfileGenerator {
    fn id(&self) -> &'static str {
        "gear_profile"
    }

    fn name(&self) -> &'static str {
        "Involute gear tooth profile"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        GEAR_PROFILE_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        let resolved = self.resolve(parameters)?;
        let profile = self.profile(parameters)?;
        let samples_per_flank = resolved.require("samples_per_flank")? as usize;
        let samples_per_arc = resolved.require("samples_per_arc")? as usize;
        let points = profile.outline_points(samples_per_flank, samples_per_arc)?;

        let mut topology = Topology::new();
        for point in &points {
            topology.add_atom(Atom::new(
                Element::CARBON,
                [point[0], point[1], 0.0],
                0.0,
                "gear_profile",
            ));
        }

        let name = format!("gear-profile-z{}", profile.teeth());
        Ok(Part::new(name, topology).with_material("geometry"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_profile() -> GearProfile {
        GearProfile::new(0.5e-9, 20, DEFAULT_PRESSURE_ANGLE_RAD).expect("valid profile")
    }

    fn count_teeth(points: &[[f64; 2]], pitch_radius_m: f64) -> usize {
        let mut teeth = 0;
        let mut in_tooth = false;
        for point in points {
            let radius_m = (point[0] * point[0] + point[1] * point[1]).sqrt();
            let inside = radius_m > pitch_radius_m;
            if inside && !in_tooth {
                teeth += 1;
            }
            in_tooth = inside;
        }
        teeth
    }

    #[test]
    fn pitch_diameter_is_module_times_tooth_count() {
        let profile = default_profile();
        assert_eq!(profile.pitch_diameter_m(), 0.5e-9 * 20.0);
        assert_eq!(profile.pitch_radius_m(), 0.5e-9 * 20.0 / 2.0);
    }

    #[test]
    fn the_profile_has_the_designed_tooth_count() {
        let profile = default_profile();
        let points = profile.outline_points(8, 4).expect("outline");
        assert_eq!(count_teeth(&points, profile.pitch_radius_m()), 20);
    }

    #[test]
    fn the_profile_closes() {
        let profile = default_profile();
        let points = profile.outline_points(8, 4).expect("outline");
        assert!(points.len() > 3);
        let first = points.first().expect("first");
        let last = points.last().expect("last");
        assert_eq!(first, last);
    }

    #[test]
    fn the_outline_stays_between_the_root_and_tip_circles() {
        let profile = default_profile();
        let points = profile.outline_points(8, 4).expect("outline");
        let root_m = profile.root_radius_m();
        let outer_m = profile.outer_radius_m();
        for point in &points {
            let radius_m = (point[0] * point[0] + point[1] * point[1]).sqrt();
            assert!(radius_m >= root_m * (1.0 - 1.0e-12));
            assert!(radius_m <= outer_m * (1.0 + 1.0e-12));
        }
        let max_radius_m = points
            .iter()
            .map(|p| (p[0] * p[0] + p[1] * p[1]).sqrt())
            .fold(0.0_f64, f64::max);
        assert!((max_radius_m - outer_m).abs() <= outer_m * 1.0e-12);
    }
    #[test]
    fn a_tooth_count_below_undercut_is_rejected() {
        let too_few = GearProfile::new(0.5e-9, 5, DEFAULT_PRESSURE_ANGLE_RAD).expect("shape");
        assert!(matches!(
            too_few.check_undercut(),
            Err(PartError::ToothCountBelowUndercut { .. })
        ));
        assert_eq!(
            GearProfile::minimum_tooth_count(DEFAULT_PRESSURE_ANGLE_RAD),
            18
        );
        assert!(GearProfileGenerator
            .generate(&ParameterSet::new().with("tooth_count", 5.0))
            .is_err());
    }

    #[test]
    fn invalid_inputs_are_errors_not_panics() {
        assert!(GearProfile::new(0.0, 20, DEFAULT_PRESSURE_ANGLE_RAD).is_err());
        assert!(GearProfile::new(0.5e-9, 2, DEFAULT_PRESSURE_ANGLE_RAD).is_err());
        assert!(GearProfile::new(0.5e-9, 20, 0.0).is_err());
        assert!(GearProfile::new(0.5e-9, 20, FRAC_PI_2).is_err());
    }

    #[test]
    fn the_generator_emits_one_atom_per_outline_point() {
        let part = GearProfileGenerator
            .generate_with_defaults()
            .expect("generate");
        let profile = default_profile();
        let points = profile.outline_points(8, 4).expect("outline");
        assert_eq!(part.atom_count(), points.len());
        assert_eq!(part.bond_count(), 0);
        assert_eq!(part.name, "gear-profile-z20");
    }

    #[test]
    fn a_generator_reports_its_identity() {
        assert_eq!(GearProfileGenerator.id(), "gear_profile");
        assert_eq!(GearProfileGenerator.parameters().len(), 5);
    }
}
