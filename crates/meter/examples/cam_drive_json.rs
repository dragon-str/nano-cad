//! Measure the pressure angle and the motion of the radial cam drive, and print
//! the result as one line of JSON.
//!
//! The pressure angle is the angle between the contact normal and the direction
//! that the rod moves. A small angle transmits the cam force into the rod well.
//! A large angle sends most of the force into the rod guide as a side load.
//!
//! The measurement is a model measurement. It reads the cam profile and the pin
//! radius only, and it models no friction, no spring preload and no contact
//! force. See ADR-0078.

use nanocad_meter::{cam_drive_report, CamDriveTarget};

/// Format a number for JSON. A non-finite value becomes null.
fn number(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.6e}")
    } else {
        "null".to_string()
    }
}

/// Format one cam target and its report as a JSON object.
fn report_json(target: &CamDriveTarget) -> String {
    let report = cam_drive_report(target);
    format!(
        "{{\"base_radius_m\":{},\"rise_m\":{},\"ramp_half_angle_rad\":{},\
         \"pin_radius_m\":{},\"samples\":{},\
         \"stroke_m\":{},\"follower_radius_min_m\":{},\"follower_radius_max_m\":{},\
         \"max_pressure_angle_rad\":{},\"max_pressure_angle_deg\":{:.3},\
         \"max_pressure_angle_at_rad\":{},\"min_mechanical_advantage\":{},\
         \"max_side_load_ratio\":{},\"max_velocity_m_per_rad\":{},\
         \"max_acceleration_m_per_rad2\":{},\"max_jerk_m_per_rad3\":{},\
         \"smooth\":{}}}",
        number(target.base_radius_m),
        number(target.rise_m),
        number(target.ramp_half_angle_rad),
        number(target.pin_radius_m),
        target.samples,
        number(report.stroke_m),
        number(report.follower_radius_min_m),
        number(report.follower_radius_max_m),
        number(report.max_pressure_angle_rad),
        report.max_pressure_angle_rad.to_degrees(),
        number(report.max_pressure_angle_at_rad),
        number(report.min_mechanical_advantage),
        number(report.max_side_load_ratio),
        number(report.max_velocity_m_per_rad),
        number(report.max_acceleration_m_per_rad2),
        number(report.max_jerk_m_per_rad3),
        report.is_smooth()
    )
}

fn main() {
    let defaults = CamDriveTarget::default();
    let gentle = CamDriveTarget {
        ramp_half_angle_rad: 1.5,
        ..CamDriveTarget::default()
    };
    println!(
        "{{\"default\":{},\"gentle\":{}}}",
        report_json(&defaults),
        report_json(&gentle)
    );
}
