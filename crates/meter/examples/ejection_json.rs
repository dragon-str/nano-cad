//! Measure the work that an ejection rod must supply to free each guest from
//! one pocket, for each wall group, and print the result as one line of JSON.
//!
//! The work is the height of the barrier above the bound state. The reference
//! figure of 10 to 40 zJ per molecule at 310 K is the reversible work at a
//! stated concentration ratio, so it is a lower bound and not this work.
//!
//! The measurement is a model measurement. The guest charges, the wall charges
//! and the nonbonded parameters come from stated sources, and no value is
//! validated against experiment. See ADR-0063, ADR-0064, ADR-0066 and
//! ADR-0067.

use nanocad_meter::{binding_energy_j, ejection_work_j, BindingTarget, EjectionTarget};
use nanocad_parts::{
    guests, BindingPocketGenerator, ParameterSet, PartGenerator, FUNCTIONAL_GROUPS,
};

/// The reference work per molecule for one sorting stage, in joules.
const REFERENCE_WORK_LOW_J: f64 = 1.0e-20;

/// Format a number for JSON. A non-finite value becomes null.
fn number(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.6e}")
    } else {
        "null".to_string()
    }
}

fn main() {
    let radius_m = std::env::args()
        .nth(1)
        .and_then(|value| value.parse::<f64>().ok())
        .unwrap_or(7.0e-10);
    let mut parameters = ParameterSet::new();
    parameters.set("pocket_radius_m", radius_m);
    parameters.set("pocket_depth_m", 1.2e-9);
    parameters.set("wall_m", 0.8e-9);
    parameters.set("thickness_m", 2.0e-9);

    let generator = BindingPocketGenerator;
    let centre_m = generator.cavity_center_m();
    let binding = BindingTarget::default();
    let target = EjectionTarget::default();

    let mut groups = Vec::new();
    for group in FUNCTIONAL_GROUPS {
        parameters.set("wall_group", group.index() as f64);
        let pocket = match generator.generate(&parameters) {
            Ok(pocket) => pocket,
            Err(error) => {
                println!("{{\"error\":\"{}\"}}", error);
                return;
            }
        };
        let mut rows = Vec::new();
        for guest in guests() {
            let report = ejection_work_j(&pocket, centre_m, guest, &binding, &target);
            let over_reference = report.work_j / REFERENCE_WORK_LOW_J;
            rows.push(format!(
                "{{\"id\":\"{}\",\"work_j\":{},\"over_kt\":{},\"over_reference\":{},\
                 \"bound_energy_j\":{},\"peak_energy_j\":{},\"released_energy_j\":{},\
                 \"peak_travel_m\":{:.6e},\"blocked_positions\":{}}}",
                report.guest,
                number(report.work_j),
                number(report.over_kt()),
                number(over_reference),
                number(report.bound_energy_j),
                number(report.peak_energy_j),
                number(report.released_energy_j),
                report.peak_travel_m,
                report.blocked_positions
            ));
        }
        let bound = binding_energy_j(&pocket, centre_m, &guests()[1], &binding);
        groups.push(format!(
            "{{\"wall_group\":\"{}\",\"pocket_atoms\":{},\
             \"ethanol_bound_j\":{},\"guests\":[{}]}}",
            group.label(),
            pocket.atom_count(),
            number(bound.minimum_energy_j),
            rows.join(",")
        ));
    }

    println!(
        "{{\"pocket_radius_m\":{:.6e},\"cavity_centre_m\":[{:.6e},{:.6e},{:.6e}],\
         \"travel_m\":{:.6e},\"steps\":{},\"reference_work_low_j\":{:.6e},\
         \"groups\":[{}]}}",
        radius_m,
        centre_m[0],
        centre_m[1],
        centre_m[2],
        target.travel_m,
        target.steps,
        REFERENCE_WORK_LOW_J,
        groups.join(",")
    );
}
