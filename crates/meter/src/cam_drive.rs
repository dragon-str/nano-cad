//! The cam drive smoothness.
//!
//! The sorting rotor turns a radial cam. Each rod carries a follower pin. The
//! pin rides the cam outer surface, so the rod slides along its own axis. The
//! push is radial, and it matches the rod axis.
//!
//! The **pressure angle** is the angle between the contact normal and the rod
//! axis. It is the one number that states whether the cam drives the rod
//! smoothly. A small angle pushes the rod along its axis. A large angle loads
//! the prismatic joint sideways and wastes force on friction.
//!
//! For a translating radial follower the angle obeys
//! `tan(angle) = |dr/dtheta| / r`, where `r` is the follower-centre radius,
//! `surface radius + pin radius`. The profile is a raised cosine over the ramp,
//! so the derivatives have a closed form.
//!
//! This metric samples one turn. It needs no forces and no time, so it carries
//! [`Fidelity::Geometric`]. Every value is a model value. All lengths are SI
//! metres. See ADR-0078.

use std::f64::consts::PI;

use nanocad_parts::{CamHubGenerator, FollowerPinGenerator};

use crate::{Fidelity, MetricValue};

/// The pressure angle at or below which the drive is judged smooth, in radians.
///
/// Thirty degrees is the usual upper limit for a cam with a translating
/// follower. Above it the side load and the friction rise quickly.
pub const SMOOTH_PRESSURE_ANGLE_RAD: f64 = PI / 6.0;

/// The cam profile under study.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CamDriveTarget {
    /// The radius of the dwell surface, in metres.
    pub base_radius_m: f64,
    /// The full lift of the ramp, in metres.
    pub rise_m: f64,
    /// The half-angle of the ramp, in radians.
    pub ramp_half_angle_rad: f64,
    /// The radius of the follower pin, in metres.
    pub pin_radius_m: f64,
    /// The number of samples over one turn.
    pub samples: usize,
}

impl Default for CamDriveTarget {
    fn default() -> Self {
        Self {
            base_radius_m: CamHubGenerator.base_radius_m(),
            rise_m: CamHubGenerator.rise_m(),
            ramp_half_angle_rad: CamHubGenerator.ramp_half_angle_rad(),
            pin_radius_m: FollowerPinGenerator.head_radius_m(),
            samples: 4096,
        }
    }
}

impl CamDriveTarget {
    /// Returns the target for the current cam and follower defaults.
    pub fn from_defaults() -> Self {
        Self::default()
    }
}

/// The measured motion of the cam drive over one turn.
#[derive(Clone, Debug, PartialEq)]
pub struct CamDriveReport {
    /// The follower lift, in metres.
    pub stroke_m: f64,
    /// The least follower-centre radius, in metres.
    pub follower_radius_min_m: f64,
    /// The greatest follower-centre radius, in metres.
    pub follower_radius_max_m: f64,
    /// The greatest pressure angle over the turn, in radians.
    pub max_pressure_angle_rad: f64,
    /// The azimuth where the greatest pressure angle occurs, in radians.
    pub max_pressure_angle_at_rad: f64,
    /// The least mechanical advantage, `cos(pressure angle)`.
    pub min_mechanical_advantage: f64,
    /// The greatest side-load ratio, `tan(pressure angle)`.
    pub max_side_load_ratio: f64,
    /// The greatest follower speed `|ds/dtheta|`, in metres per radian.
    pub max_velocity_m_per_rad: f64,
    /// The greatest follower acceleration `|d2s/dtheta2|`, in metres per
    /// radian squared.
    pub max_acceleration_m_per_rad2: f64,
    /// The greatest follower jerk `|d3s/dtheta3|`, in metres per radian cubed.
    pub max_jerk_m_per_rad3: f64,
    /// The number of samples.
    pub samples: usize,
    /// The fidelity of the measurement.
    pub fidelity: Fidelity,
}

impl CamDriveReport {
    /// Returns true when the pressure angle stays within the smooth limit.
    pub fn is_smooth(&self) -> bool {
        self.max_pressure_angle_rad <= SMOOTH_PRESSURE_ANGLE_RAD
    }

    /// Returns the report as named metric values.
    pub fn metric_values(&self) -> Vec<MetricValue> {
        vec![
            MetricValue {
                name: "cam_stroke".to_string(),
                value: self.stroke_m,
                unit: "m".to_string(),
                fidelity: self.fidelity,
                note: format!("{:.3} nm rise", self.stroke_m * 1e9),
            },
            MetricValue {
                name: "cam_max_pressure_angle".to_string(),
                value: self.max_pressure_angle_rad,
                unit: "rad".to_string(),
                fidelity: self.fidelity,
                note: format!(
                    "{:.1} deg at {:.1} deg; {:.2} mechanical advantage",
                    self.max_pressure_angle_rad.to_degrees(),
                    self.max_pressure_angle_at_rad.to_degrees(),
                    self.min_mechanical_advantage
                ),
            },
            MetricValue {
                name: "cam_max_jerk".to_string(),
                value: self.max_jerk_m_per_rad3,
                unit: "m/rad^3".to_string(),
                fidelity: self.fidelity,
                note: "greatest follower jerk over the turn".to_string(),
            },
        ]
    }
}

/// Samples the cam profile over one turn and measures the drive.
pub fn cam_drive_report(target: &CamDriveTarget) -> CamDriveReport {
    let samples = target.samples.max(8);
    let step_rad = 2.0 * PI / samples as f64;
    let ramp_rad = target.ramp_half_angle_rad;
    let amplitude_m = target.rise_m;

    let mut max_pressure_angle_rad = 0.0_f64;
    let mut max_pressure_angle_at_rad = 0.0_f64;
    let mut max_velocity_m_per_rad = 0.0_f64;
    let mut max_acceleration_m_per_rad2 = 0.0_f64;
    let mut max_jerk_m_per_rad3 = 0.0_f64;

    for index in 0..samples {
        let delta_rad = -PI + (index as f64 + 0.5) * step_rad;
        let surface_m = CamHubGenerator::surface_radius_with(
            target.base_radius_m,
            amplitude_m,
            ramp_rad,
            delta_rad,
        );
        let follower_m = surface_m + target.pin_radius_m;

        let (velocity, acceleration, jerk) = profile_derivatives(amplitude_m, ramp_rad, delta_rad);
        let pressure_angle_rad = (velocity.abs() / follower_m).atan();
        if pressure_angle_rad > max_pressure_angle_rad {
            max_pressure_angle_rad = pressure_angle_rad;
            max_pressure_angle_at_rad = delta_rad;
        }
        max_velocity_m_per_rad = max_velocity_m_per_rad.max(velocity.abs());
        max_acceleration_m_per_rad2 = max_acceleration_m_per_rad2.max(acceleration.abs());
        max_jerk_m_per_rad3 = max_jerk_m_per_rad3.max(jerk.abs());
    }

    let follower_radius_min_m = target.base_radius_m + target.pin_radius_m;
    let follower_radius_max_m = target.base_radius_m + target.rise_m + target.pin_radius_m;

    let min_mechanical_advantage = max_pressure_angle_rad.cos();
    let max_side_load_ratio = max_pressure_angle_rad.tan();

    CamDriveReport {
        stroke_m: amplitude_m,
        follower_radius_min_m,
        follower_radius_max_m,
        max_pressure_angle_rad,
        max_pressure_angle_at_rad,
        min_mechanical_advantage,
        max_side_load_ratio,
        max_velocity_m_per_rad,
        max_acceleration_m_per_rad2,
        max_jerk_m_per_rad3,
        samples,
        fidelity: Fidelity::Geometric,
    }
}

/// Returns `(ds/dtheta, d2s/dtheta2, d3s/dtheta3)` for the raised-cosine ramp.
///
/// The derivatives are zero on the dwell. Inside the ramp they follow the
/// closed form of `rise * 0.5 * (1 + cos(PI * delta / ramp))`.
fn profile_derivatives(rise_m: f64, ramp_rad: f64, delta_rad: f64) -> (f64, f64, f64) {
    if ramp_rad <= 0.0 || delta_rad.abs() >= ramp_rad {
        return (0.0, 0.0, 0.0);
    }
    let k = PI / ramp_rad;
    let phase = k * delta_rad;
    let velocity = -0.5 * rise_m * k * phase.sin();
    let acceleration = -0.5 * rise_m * k * k * phase.cos();
    let jerk = 0.5 * rise_m * k * k * k * phase.sin();
    (velocity, acceleration, jerk)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_default_cam_has_a_steep_pressure_angle() {
        let report = cam_drive_report(&CamDriveTarget::default());
        assert!(
            (report.max_pressure_angle_rad - 1.30).abs() < 0.02,
            "default pressure angle {:.4} rad",
            report.max_pressure_angle_rad
        );
        assert!(!report.is_smooth());
        assert!(report.min_mechanical_advantage < 0.3);
    }

    #[test]
    fn the_stroke_equals_the_rise() {
        let target = CamDriveTarget::default();
        let report = cam_drive_report(&target);
        assert_eq!(report.stroke_m, target.rise_m);
        assert!(
            (report.follower_radius_min_m - (target.base_radius_m + target.pin_radius_m)).abs()
                < 1e-15
        );
        assert!(
            (report.follower_radius_max_m
                - (target.base_radius_m + target.rise_m + target.pin_radius_m))
                .abs()
                < 1e-15
        );
    }

    #[test]
    fn the_velocity_matches_the_closed_form() {
        let target = CamDriveTarget::default();
        let report = cam_drive_report(&target);
        let expected = target.rise_m * PI / (2.0 * target.ramp_half_angle_rad);
        assert!(
            (report.max_velocity_m_per_rad - expected).abs() / expected < 1e-3,
            "velocity {:.6e} vs {:.6e}",
            report.max_velocity_m_per_rad,
            expected
        );
    }

    #[test]
    fn the_velocity_matches_a_finite_difference() {
        let target = CamDriveTarget::default();
        let step = 1e-6;
        let delta = 0.5 * target.ramp_half_angle_rad;
        let surface_plus = CamHubGenerator::surface_radius_with(
            target.base_radius_m,
            target.rise_m,
            target.ramp_half_angle_rad,
            delta + step,
        );
        let surface_minus = CamHubGenerator::surface_radius_with(
            target.base_radius_m,
            target.rise_m,
            target.ramp_half_angle_rad,
            delta - step,
        );
        let finite_difference = (surface_plus - surface_minus) / (2.0 * step);
        let (velocity, _, _) =
            profile_derivatives(target.rise_m, target.ramp_half_angle_rad, delta);
        assert!(
            (velocity - finite_difference).abs() / velocity.abs() < 1e-6,
            "velocity {:.6e} vs finite difference {:.6e}",
            velocity,
            finite_difference
        );
    }

    #[test]
    fn a_gentle_ramp_is_smooth() {
        let target = CamDriveTarget {
            ramp_half_angle_rad: 1.5,
            ..CamDriveTarget::default()
        };
        let report = cam_drive_report(&target);
        assert!(
            report.is_smooth(),
            "angle {:.4} rad",
            report.max_pressure_angle_rad
        );
    }

    #[test]
    fn a_larger_base_lowers_the_pressure_angle() {
        let base = cam_drive_report(&CamDriveTarget::default());
        let wide = cam_drive_report(&CamDriveTarget {
            base_radius_m: 4.0e-9,
            ..CamDriveTarget::default()
        });
        assert!(wide.max_pressure_angle_rad < base.max_pressure_angle_rad);
    }
}
