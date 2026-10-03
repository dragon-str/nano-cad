//! Builds the respirocyte oxygen module, relaxed or not, audits it, and
//! optionally writes the scene JSON.
//!
//! Usage: `module_json [relaxed|raw] [path] [key=value ...]`; the pairs set
//! rotor parameters. The last stdout line is a JSON
//! object with the per-body relaxation summary and the sense audit.

use nanocad_jigs::{
    body_contacts, build_module_scene_with, module_relax_options, sense_report, write_scene_json,
};
use nanocad_meter::ValenceRelaxOptions;
use nanocad_parts::ParameterSet;

fn main() {
    let mut args = std::env::args().skip(1);
    let mode = args.next().unwrap_or_else(|| "relaxed".to_string());
    let path = args.next();
    // Optional rotor parameters as key=value after the path.
    let mut parameters = ParameterSet::new();
    for pair in args {
        if let Some((key, value)) = pair.split_once('=') {
            if let Ok(v) = value.parse::<f64>() {
                parameters = parameters.with(key, v);
            }
        }
    }
    let options = ValenceRelaxOptions {
        max_iterations: std::env::var("RELAX_ITERATIONS")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(module_relax_options().max_iterations),
        ..module_relax_options()
    };
    let relax = if mode == "raw" { None } else { Some(&options) };
    let (scene, reports) = match build_module_scene_with(&parameters, relax) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(1);
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
