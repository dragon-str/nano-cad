//! Audits the rotor scene against the physical-sense rules and prints the report.
//!
//! The rules are: no overlapping atoms, every rigid body is one connected
//! lattice, every non-fixed body is held by a joint, and no bond crosses between
//! bodies that no joint joins.

use nanocad_jigs::{build_part_scene, build_rotor_scene, sense_report};

fn main() {
    let target = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "rotor".to_string());
    let scene = if target == "rotor" {
        build_rotor_scene()
    } else {
        build_part_scene(&target)
    };
    let scene = match scene {
        Ok(scene) => scene,
        Err(error) => {
            eprintln!("the scene failed: {error}");
            std::process::exit(1);
        }
    };
    let report = sense_report(&scene);
    match serde_json::to_string(&report) {
        Ok(json) => println!("{json}"),
        Err(error) => {
            eprintln!("the report did not write: {error}");
            std::process::exit(1);
        }
    }
}
