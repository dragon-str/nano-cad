//! Estimate the loading of one pocket from the oxygen and carbon dioxide
//! partial pressures, and print the result as one line of JSON.
//!
//! The model is a single site in a gas. The arrival rate is the gas kinetic
//! flux onto the pocket mouth. The site reaches the target occupancy when that
//! many molecules have arrived. The model states no diffusion in a liquid, no
//! crowding and a rigid pocket. Treat the result as an estimate. See ADR-0081.

use nanocad_meter::{loading_report, GasLoadingReport, LoadingTarget};

/// Format a number for JSON. A non-finite value becomes null.
fn number(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.6e}")
    } else {
        "null".to_string()
    }
}

/// Format one gas estimate as a JSON object.
fn gas_json(report: &GasLoadingReport) -> String {
    format!(
        "{{\"gas\":\"{}\",\"number_density_per_m3\":{},\
         \"mean_speed_m_per_s\":{},\"capture_rate_m3_per_s\":{},\
         \"saturation_time_s\":{},\"loading_angle_rad\":{},\
         \"loading_angle_deg\":{:.3}}}",
        report.gas,
        number(report.number_density_per_m3),
        number(report.mean_speed_m_per_s),
        number(report.capture_rate_m3_per_s),
        number(report.saturation_time_s),
        number(report.loading_angle_rad),
        report.loading_angle_rad.to_degrees()
    )
}

/// Format one loading target and its report as a JSON object.
fn report_json(target: &LoadingTarget) -> String {
    let report = loading_report(target);
    format!(
        "{{\"o2_pressure_pa\":{},\"co2_pressure_pa\":{},\"temperature_k\":{},\
         \"opening_radius_m\":{},\"sticking_probability\":{},\
         \"angular_rate_rad_per_s\":{},\"target_occupancy\":{},\
         \"opening_area_m2\":{},\"o2\":{},\"co2\":{},\
         \"loading_angle_rad\":{},\"loading_angle_deg\":{:.3}}}",
        number(target.o2_pressure_pa),
        number(target.co2_pressure_pa),
        number(target.temperature_k),
        number(target.opening_radius_m),
        number(target.sticking_probability),
        number(target.angular_rate_rad_per_s),
        number(target.target_occupancy),
        number(report.opening_area_m2),
        gas_json(&report.o2),
        gas_json(&report.co2),
        number(report.loading_angle_rad),
        report.loading_angle_rad.to_degrees()
    )
}

fn main() {
    println!("{}", report_json(&LoadingTarget::default()));
}
