//! Relaxes each body of the gear pump with the valence force field and
//! prints one JSON line per body: the energy drop, the longest C-C bond
//! before and after, the displacements, the closest nonbonded pair and the
//! run time. Arguments: `gears` or `all`, then an optional iteration budget.

use std::time::Instant;

use nanocad_meter::{relax_valence, ValenceRelaxOptions};
use nanocad_parts::{GearPumpGenerator, ParameterSet};

fn main() {
    let set = match GearPumpGenerator.build(&ParameterSet::new()) {
        Ok(set) => set,
        Err(error) => {
            eprintln!("the pump failed: {error}");
            std::process::exit(1);
        }
    };
    let which = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "gears".to_string());
    let mut bodies = vec![("gear_a", set.gear_a_atoms.clone())];
    if which == "all" {
        bodies.push(("housing", set.housing_atoms.clone()));
    }
    for (name, range) in bodies {
        let start = Instant::now();
        let iterations = std::env::args()
            .nth(2)
            .and_then(|text| text.parse().ok())
            .unwrap_or(ValenceRelaxOptions::default().max_iterations);
        let options = ValenceRelaxOptions {
            max_iterations: iterations,
            ..ValenceRelaxOptions::default()
        };
        match relax_valence(&set.part.topology, range, &options) {
            Ok(report) => println!(
                "{{\"body\":\"{name}\",\"atoms\":{},\"bonds\":{},\"angles\":{},\
                 \"initial_energy_j\":{:e},\"final_energy_j\":{:e},\"iterations\":{},\
                 \"converged\":{},\"gradient_n\":{:e},\"initial_max_c_c_m\":{:e},\
                 \"max_c_c_m\":{:e},\"max_displacement_m\":{:e},\"max_carbon_displacement_m\":{:e},\
                 \"rms_displacement_m\":{:e},\"min_nonbonded_m\":{:e},\"seconds\":{:.1}}}",
                report.atom_count,
                report.bond_count,
                report.angle_count,
                report.initial_energy_j,
                report.final_energy_j,
                report.iterations,
                report.converged,
                report.gradient_norm_n,
                report.initial_max_c_c_m,
                report.max_c_c_m,
                report.max_displacement_m,
                report.max_carbon_displacement_m,
                report.rms_displacement_m,
                report.min_nonbonded_m,
                start.elapsed().as_secs_f64(),
            ),
            Err(error) => {
                eprintln!("{name} did not relax: {error}");
                std::process::exit(1);
            }
        }
    }
}
