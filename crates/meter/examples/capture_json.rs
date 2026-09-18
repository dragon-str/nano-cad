//! Measures the sidewall capture of the sorting rotor.
//!
//! The example reads the default ejection rod and the default sorting rotor.
//! It then reports, for each guest, the fraction of the pocket that the rod
//! can reach. It compares the narrow tip face with the wider shaft face, and
//! it states the face radius that closes the dead zone.
//!
//! Run it with `cargo run -p nanocad-meter --example capture_json`.

use nanocad_parts::{
    guest, EjectionRodGenerator, ParameterSet, PartGenerator, SortingRotorGenerator,
};

use nanocad_meter::{capture_report, CaptureTarget};

fn resolved(generator: &impl PartGenerator, name: &str) -> f64 {
    generator
        .resolve(&ParameterSet::new())
        .ok()
        .and_then(|values| values.get(name))
        .unwrap_or(0.0)
}

fn report_json(label: &str, target: &CaptureTarget) -> String {
    let report = capture_report(target);
    format!(
        "{{\"face\":\"{label}\",\"pocket_radius_m\":{:e},\"guest_radius_m\":{:e},\
         \"piston_radius_m\":{:e},\"allowed_offset_m\":{:e},\"reach_m\":{:e},\
         \"dead_zone_m\":{:e},\"coverage\":{:.6},\"complete\":{},\
         \"required_piston_radius_m\":{:e}}}",
        report.pocket_radius_m,
        report.guest_radius_m,
        report.piston_radius_m,
        report.allowed_offset_m,
        report.reach_m,
        report.dead_zone_m,
        report.coverage,
        report.complete,
        report.required_piston_radius_m,
    )
}

fn main() {
    let pocket_radius_m = resolved(&SortingRotorGenerator, "pocket_radius_m");
    let shaft_radius_m = resolved(&EjectionRodGenerator, "shaft_radius_m");
    let tip_radius_m = resolved(&EjectionRodGenerator, "tip_radius_m");

    for id in ["ethanol", "dimethyl_ether"] {
        let guest = match guest(id) {
            Some(guest) => guest,
            None => {
                eprintln!("the guest {id} is not in the library");
                std::process::exit(1);
            }
        };
        let guest_radius_m = guest.extent_m();
        let tip = CaptureTarget {
            pocket_radius_m,
            guest_radius_m,
            piston_radius_m: tip_radius_m,
        };
        let shaft = CaptureTarget {
            pocket_radius_m,
            guest_radius_m,
            piston_radius_m: shaft_radius_m,
        };
        println!(
            "{{\"guest\":\"{id}\",\"extent_m\":{guest_radius_m:e},\"faces\":[{},{}]}}",
            report_json("tip", &tip),
            report_json("shaft", &shaft),
        );
    }
}
