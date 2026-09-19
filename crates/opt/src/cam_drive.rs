//! The cam-drive search: a CMA-ES over the radial cam profile.

use std::f64::consts::PI;

use nanocad_meter::{cam_drive_report, CamDriveTarget};

use crate::cmaes::{Objective, ParameterBounds};

/// The stroke that the rotor needs, in metres. The default rise gives it.
pub const STROKE_TARGET_M: f64 = 1.5e-9;

/// The largest ramp half-angle that keeps one pin on the ramp, in radians.
///
/// The rotor has 12 pockets, so two pins are `2*PI/12 = 30` degrees apart. A
/// ramp half-angle above half that spacing puts two pins on the ramp at once.
pub const POCKET_RAMP_LIMIT_RAD: f64 = PI / 12.0;

/// The largest base radius that still fits under the rotor, in metres.
///
/// The lobe outer radius is `base + rise`. It must clear the pocket inner
/// radius, and the pin rides on the lobe, so a larger base needs a larger
/// rotor. This value keeps the hub inside the present rotor.
pub const HUB_BASE_LIMIT_M: f64 = 3.0e-9;

/// The jerk that costs one cost unit, in metres per cubic radian.
const JERK_SCALE_M_PER_RAD3: f64 = 1.0e-7;

/// The weight of the jerk term, relative to the pressure angle in radians.
const JERK_WEIGHT: f64 = 0.05;

/// The result of one cam-drive design.
#[derive(Clone, Debug, PartialEq)]
pub struct CamDriveScore {
    /// The follower stroke, in metres.
    pub stroke_m: f64,
    /// The largest pressure angle over one turn, in radians.
    pub max_pressure_angle_rad: f64,
    /// The largest jerk over one turn, in metres per cubic radian.
    pub max_jerk_m_per_rad3: f64,
    /// The parameter vector that produced the score.
    pub parameters: Vec<f64>,
    /// True when the stroke meets or exceeds its target.
    pub feasible: bool,
    /// The distance of the stroke from the target, in metres. It is zero when
    /// the stroke equals the target.
    pub stroke_error_m: f64,
}

/// A CMA-ES search over the base radius, the ramp and the pin radius.
///
/// The rise holds the stroke target, so every design is feasible. The cost is
/// the largest pressure angle over one turn, plus a small jerk term. The ramp
/// stays inside one pocket spacing.
#[derive(Clone, Debug, PartialEq)]
pub struct CamDriveSearch {
    /// The searched parameter names, in vector order.
    pub names: Vec<&'static str>,
    /// The range of each searched parameter.
    pub bounds: Vec<ParameterBounds>,
    /// The stroke target, in metres.
    pub stroke_target_m: f64,
    /// The ramp half-angle at which the cost is measured, in radians.
    pub ramp_limit_rad: f64,
}

impl Default for CamDriveSearch {
    fn default() -> Self {
        let ramp_limit_rad = POCKET_RAMP_LIMIT_RAD;
        Self {
            names: vec!["base_radius_m", "ramp_half_angle_rad", "pin_radius_m"],
            bounds: vec![
                ParameterBounds::new(1.0e-9, HUB_BASE_LIMIT_M, false),
                ParameterBounds::new(0.05, ramp_limit_rad, false),
                ParameterBounds::new(0.15e-9, 0.80e-9, false),
            ],
            stroke_target_m: STROKE_TARGET_M,
            ramp_limit_rad,
        }
    }
}

impl CamDriveSearch {
    /// Builds a search with the default ranges.
    pub fn new() -> Self {
        Self::default()
    }

    /// Returns the starting vector, at the cam defaults, clamped to the ranges.
    pub fn initial(&self) -> Vec<f64> {
        let target = CamDriveTarget::default();
        let defaults = [
            target.base_radius_m,
            target.ramp_half_angle_rad,
            target.pin_radius_m,
        ];
        defaults
            .iter()
            .zip(self.bounds.iter())
            .map(|(value, bounds)| bounds.clamp(*value))
            .collect()
    }

    /// Builds a cam target from a parameter vector.
    pub fn target_of(&self, vector: &[f64]) -> CamDriveTarget {
        let mut target = CamDriveTarget {
            rise_m: self.stroke_target_m,
            ..CamDriveTarget::default()
        };
        for (name, value) in self.names.iter().zip(vector.iter()) {
            match *name {
                "base_radius_m" => target.base_radius_m = *value,
                "rise_m" => target.rise_m = *value,
                "ramp_half_angle_rad" => target.ramp_half_angle_rad = *value,
                "pin_radius_m" => target.pin_radius_m = *value,
                _ => {}
            }
        }
        target
    }

    /// Scores one design. It is always valid, because the metric is analytic.
    pub fn score(&self, vector: &[f64]) -> CamDriveScore {
        let target = self.target_of(vector);
        let report = cam_drive_report(&target);
        let stroke_m = report.stroke_m;
        let feasible = stroke_m >= self.stroke_target_m;
        let stroke_error_m = (stroke_m - self.stroke_target_m).abs();
        CamDriveScore {
            stroke_m,
            max_pressure_angle_rad: report.max_pressure_angle_rad,
            max_jerk_m_per_rad3: report.max_jerk_m_per_rad3,
            parameters: vector.to_vec(),
            feasible,
            stroke_error_m,
        }
    }
}

impl Objective for CamDriveSearch {
    fn evaluate(&self, parameters: &[f64]) -> Option<f64> {
        let score = self.score(parameters);
        let jerk_term = score.max_jerk_m_per_rad3 / JERK_SCALE_M_PER_RAD3 * JERK_WEIGHT;
        Some(score.max_pressure_angle_rad + jerk_term)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cmaes::{minimize, CmaEsOptions};

    #[test]
    fn the_starting_vector_is_inside_every_range() {
        let search = CamDriveSearch::default();
        let initial = search.initial();
        assert_eq!(initial.len(), search.bounds.len());
        for (value, bounds) in initial.iter().zip(search.bounds.iter()) {
            assert!(*value >= bounds.min && *value <= bounds.max);
        }
    }

    #[test]
    fn the_default_cam_is_not_smooth() {
        let search = CamDriveSearch::default();
        let score = search.score(&search.initial());
        assert!(score.feasible);
        assert!(score.max_pressure_angle_rad > PI / 6.0);
    }

    #[test]
    fn the_search_lowers_the_pressure_angle() {
        let search = CamDriveSearch::default();
        let initial = search.initial();
        let start = search.score(&initial);
        let start_cost = search.evaluate(&initial).unwrap();
        let report = minimize(
            &initial,
            &search.bounds,
            &CmaEsOptions {
                generations: 60,
                ..CmaEsOptions::default()
            },
            &search,
        );
        assert!(report.best_cost <= start_cost);
        let best = search.score(&report.best_parameters);
        assert!(best.feasible);
        assert!(best.max_pressure_angle_rad < start.max_pressure_angle_rad);
    }

    #[test]
    fn a_smooth_angle_needs_a_larger_hub() {
        let search = CamDriveSearch::default();
        let baseline = cam_drive_report(&CamDriveTarget {
            rise_m: STROKE_TARGET_M,
            ramp_half_angle_rad: POCKET_RAMP_LIMIT_RAD,
            ..CamDriveTarget::default()
        });
        assert!(!baseline.is_smooth());
        let large = cam_drive_report(&CamDriveTarget {
            base_radius_m: 2.0e-8,
            rise_m: STROKE_TARGET_M,
            ramp_half_angle_rad: POCKET_RAMP_LIMIT_RAD,
            ..CamDriveTarget::default()
        });
        assert!(large.is_smooth());
        let _ = search;
    }
}
