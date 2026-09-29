//! Builds the respirocyte oxygen module, relaxed or not, audits it, and
//! optionally writes the scene JSON.
//!
//! Usage: `module_json [relaxed|raw] [path]`. The last stdout line is a JSON
//! object with the per-body relaxation summary and the sense audit.

use nanocad_jigs::{
    body_contacts, build_module_scene, build_relaxed_module_scene, module_relax_options,
    sense_report, write_scene_json,
};
use nanocad_meter::ValenceRelaxOptions;

fn main() {
    let mut args = std::env::args().skip(1);
    let mode = args.next().unwrap_or_else(|| "relaxed".to_string());
    let path = args.next();
    let (scene, reports) = if mode == "raw" {
        match build_module_scene() {
            Ok(s) => (s, Vec::new()),
            Err(e) => {
                eprintln!("{e}");
                std::process::exit(1);
            }
        }
    } else {
        match build_relaxed_module_scene(&ValenceRelaxOptions {
            max_iterations: std::env::var("RELAX_ITERATIONS")
                .ok()
                .and_then(|v| v.parse().ok())
                .unwrap_or(module_relax_options().max_iterations),
            ..module_relax_options()
        }) {
            Ok(r) => r,
            Err(e) => {
                eprintln!("{e}");
                std::process::exit(1);
            }
        }
    };
    if let Some(path) = path {
        if let Err(e) = write_scene_json(&scene, &path, false) {
            eprintln!("{e}");
            std::process::exit(1);
        }
    }
    let relax: Vec<String> = reports
        .iter()
        .filter(|(name, _)| !name.starts_with("heme_site_") && !name.starts_with("ligand_arm_") || name.ends_with("_0") || name.ends_with("_6"))
        .map(|(name, r)| {
            format!(
                "{{\"body\":\"{name}\",\"atoms\":{},\"converged\":{},\"iterations\":{},\"max_c_c_m\":{:.3e},\"rms_displacement_m\":{:.3e},\"max_displacement_m\":{:.3e},\"min_nonbonded_m\":{:.3e}}}",
                r.atom_count, r.converged, r.iterations, r.max_c_c_m, r.rms_displacement_m, r.max_displacement_m, r.min_nonbonded_m
            )
        })
        .collect();
    let audit = serde_json::to_string(&sense_report(&scene)).unwrap_or_default();
    // Contacts between every pair of bodies, joints included. Fe-N at 2.1 A
    // in a loaded pocket is the intended ligand bond.
    let contacts: Vec<String> = body_contacts(&scene, 2.2e-10)
        .iter()
        .map(|(a, b, d)| format!("[{a},{b},{d:.3e}]"))
        .collect();
    println!(
        "{{\"contacts_below_2p2A\":[{}],\"atoms\":{},\"bodies\":{},\"relax\":[{}],\"sense\":{audit}}}",
        contacts.join(","),
        scene.atomistic.atom_count,
        scene.device.body_count,
        relax.join(",")
    );
}
