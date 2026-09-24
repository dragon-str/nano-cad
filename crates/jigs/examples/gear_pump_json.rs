//! Builds the relaxed gear-pump scene, audits it, and prints the facts.
//!
//! With one argument, the example also writes the scene JSON to that path.
//! The last stdout line is one JSON object: the relaxation reports of gear A
//! and the housing, and the sense audit of the relaxed scene.

use nanocad_jigs::{build_relaxed_gear_pump_scene, sense_report, write_scene_json};

fn main() {
    let (scene, reports) = match build_relaxed_gear_pump_scene() {
        Ok(result) => result,
        Err(error) => {
            eprintln!("the pump scene failed: {error}");
            std::process::exit(1);
        }
    };
    if let Some(path) = std::env::args().nth(1) {
        if let Err(error) = write_scene_json(&scene, &path, false) {
            eprintln!("the scene was not written: {error}");
            std::process::exit(1);
        }
    }
    let relax: Vec<String> = reports
        .iter()
        .zip(["gear", "housing"])
        .map(|(report, name)| {
            format!(
                "{{\"body\":\"{name}\",\"atoms\":{},\"initial_energy_j\":{:e},\
                 \"final_energy_j\":{:e},\"iterations\":{},\"converged\":{},\
                 \"gradient_n\":{:e},\"initial_max_c_c_m\":{:e},\"max_c_c_m\":{:e},\
                 \"max_displacement_m\":{:e}}}",
                report.atom_count,
                report.initial_energy_j,
                report.final_energy_j,
                report.iterations,
                report.converged,
                report.gradient_norm_n,
                report.initial_max_c_c_m,
                report.max_c_c_m,
                report.max_displacement_m,
            )
        })
        .collect();
    let audit = match serde_json::to_string(&sense_report(&scene)) {
        Ok(json) => json,
        Err(error) => {
            eprintln!("the audit did not write: {error}");
            std::process::exit(1);
        }
    };
    println!("{{\"relax\":[{}],\"sense\":{audit}}}", relax.join(","));
}
