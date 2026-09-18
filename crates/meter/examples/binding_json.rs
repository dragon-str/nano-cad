//! Measure the binding energy of every guest in one pocket, for each wall
//! group, and print the result as one line of JSON.
//!
//! The measurement is a model measurement. The guest charges, the wall charges
//! and the nonbonded parameters come from stated sources, and no value is
//! validated against experiment. See ADR-0063, ADR-0064 and ADR-0066.

use nanocad_meter::{binding_energy_j, selectivity_ratio, BindingTarget};

/// Format a number for JSON. A non-finite value becomes null.
fn number(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.6e}")
    } else {
        "null".to_string()
    }
}
use nanocad_parts::{
    guests, BindingPocketGenerator, ParameterSet, PartGenerator, FUNCTIONAL_GROUPS,
};

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
    let target = BindingTarget::default();

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
            let report = binding_energy_j(&pocket, centre_m, guest, &target);
            rows.push(format!(
                "{{\"id\":\"{}\",\"minimum_energy_j\":{},\"over_kt\":{},\
                 \"displacement_m\":{:.6e},\"azimuth_rad\":{:.6},\
                 \"centre_energy_j\":{},\"blocked_positions\":{}}}",
                report.guest,
                number(report.minimum_energy_j),
                number(report.over_kt()),
                report.displacement_m,
                report.azimuth_rad,
                number(report.centre_energy_j),
                report.blocked_positions
            ));
        }
        let ethanol_j = binding_energy_j(&pocket, centre_m, &guests()[1], &target).minimum_energy_j;
        let ether_j = binding_energy_j(&pocket, centre_m, &guests()[2], &target).minimum_energy_j;
        groups.push(format!(
            "{{\"wall_group\":\"{}\",\"pocket_atoms\":{},\"isomer_ratio\":{:.4},\
             \"guests\":[{}]}}",
            group.label(),
            pocket.atom_count(),
            selectivity_ratio(ethanol_j, ether_j),
            rows.join(",")
        ));
    }

    println!(
        "{{\"pocket_radius_m\":{:.6e},\"cavity_centre_m\":[{:.6e},{:.6e},{:.6e}],\
         \"displacements\":{},\"azimuths\":{},\"cutoff_m\":{:.6e},\"groups\":[{}]}}",
        radius_m,
        centre_m[0],
        centre_m[1],
        centre_m[2],
        target.displacements,
        target.azimuths,
        target.cutoff_m,
        groups.join(",")
    );
}
