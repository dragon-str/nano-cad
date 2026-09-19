//! Search the radial cam profile for the least pressure angle, and print the
//! result as one line of JSON.
//!
//! The search keeps the ramp inside one pocket spacing, so only one pin rides
//! the ramp at a time. The stroke stays at the rotor target.
//!
//! The search is a model search. It reads the cam profile and the pin radius
//! only. It models no friction, no spring preload and no contact force. See
//! ADR-0078.

use nanocad_opt::{minimize, CamDriveSearch, CmaEsOptions, Objective};

/// Format a number for JSON. A non-finite value becomes null.
fn number(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.6e}")
    } else {
        "null".to_string()
    }
}

fn main() {
    let search = CamDriveSearch::default();
    let initial = search.initial();
    let options = CmaEsOptions {
        generations: 200,
        population: 12,
        seed: 1,
        ..CmaEsOptions::default()
    };
    let report = minimize(&initial, &search.bounds, &options, &search);
    let start = search.score(&initial);
    let best = search.score(&report.best_parameters);
    let start_cost = search.evaluate(&initial).unwrap();
    let best_cost = search.evaluate(&report.best_parameters).unwrap();

    let mut names = Vec::new();
    for (name, value) in search.names.iter().zip(best.parameters.iter()) {
        names.push(format!("\"{name}\":{}", number(*value)));
    }

    println!(
        "{{\"start_pressure_angle_deg\":{:.3},\"best_pressure_angle_deg\":{:.3},\
         \"best_stroke_m\":{},\"best_jerk_m_per_rad3\":{},\"feasible\":{},\
         \"smooth\":{},\"start_cost\":{},\"best_cost\":{},\"generations\":{},\
         \"evaluations\":{},\"converged\":{},\"parameters\":{{{}}}}}",
        start.max_pressure_angle_rad.to_degrees(),
        best.max_pressure_angle_rad.to_degrees(),
        number(best.stroke_m),
        number(best.max_jerk_m_per_rad3),
        best.feasible,
        best.max_pressure_angle_rad <= std::f64::consts::PI / 6.0 && best.feasible,
        number(start_cost),
        number(best_cost),
        report.generations,
        report.evaluations,
        report.converged,
        names.join(",")
    );
}
