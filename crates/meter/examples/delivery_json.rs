//! Compare the delivered fraction of several passive sites over the oxygen and
//! carbon dioxide swings of the body, and print the result as one line of JSON.
//!
//! A site follows the Hill curve. The delivered fraction of one cycle is the
//! occupancy at the loading pressure minus the occupancy at the unloading
//! pressure. The oxygen swing is large and the carbon dioxide swing is small.
//!
//! The model states no kinetics, no hysteresis and no competing ligand. Treat
//! every result as an estimate. See ADR-0082.

use nanocad_meter::{delivery_report, optimal_p50_pa, DeliveryReport, DeliveryTarget};

/// Format a number for JSON. A non-finite value becomes null.
fn number(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.6e}")
    } else {
        "null".to_string()
    }
}

/// Format one delivery report as a JSON object.
fn report_json(name: &str, report: &DeliveryReport, hill_coefficient: f64) -> String {
    format!(
        "{{\"design\":\"{name}\",\"load_pressure_pa\":{},\"unload_pressure_pa\":{},\
         \"p50_pa\":{},\"binding_energy_j\":{},\"hill_coefficient\":{},\
         \"occupancy_loaded\":{:.6},\"occupancy_unloaded\":{:.6},\
         \"delivered_fraction\":{:.6},\"pressure_ratio\":{:.3}}}",
        number(report.load_pressure_pa),
        number(report.unload_pressure_pa),
        number(report.p50_pa),
        number(report.binding_energy_j),
        hill_coefficient,
        report.occupancy_loaded,
        report.occupancy_unloaded,
        report.delivered_fraction,
        report.pressure_ratio()
    )
}

/// Build a site that is tuned to one swing, with the stated cooperativity.
fn tuned(load_pa: f64, unload_pa: f64, hill_coefficient: f64) -> DeliveryTarget {
    DeliveryTarget {
        load_pressure_pa: load_pa,
        unload_pressure_pa: unload_pa,
        p50_pa: optimal_p50_pa(load_pa, unload_pa),
        hill_coefficient,
        ..DeliveryTarget::default()
    }
}

fn main() {
    let designs = [
        ("o2_flat_blood", tuned(13.3e3, 5.3e3, 1.0)),
        ("o2_flat_tissue", tuned(13.3e3, 3.0e3, 1.0)),
        ("o2_hill_tissue", tuned(13.3e3, 3.0e3, 2.8)),
        ("co2_blood", tuned(6.1e3, 5.3e3, 1.0)),
        ("co2_tissue", tuned(6.7e3, 5.3e3, 1.0)),
    ];
    let rows: Vec<String> = designs
        .iter()
        .map(|(name, target)| report_json(name, &delivery_report(target), target.hill_coefficient))
        .collect();
    println!("{{\"designs\":[{}]}}", rows.join(","));
}
