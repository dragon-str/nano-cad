//! Print the water kill test for the gear-pump site as JSON (ADR-0086).
//!
//! Two sites compete oxygen, water and carbon dioxide at alveolar partial
//! pressures and body temperature. Water enters at its vapour pressure, which
//! is exact for a site in the liquid. The enthalpies are estimates: the metal
//! site uses the Fe(NH3)4 quantum-chemistry proxy, and the bare pocket uses the
//! classical methyl wall at a 4.5e-10 m radius.

use nanocad_meter::{competition_report, CompetitionReport, Competitor, AVOGADRO_PER_MOL};

const BODY_TEMPERATURE_K: f64 = 310.15;
const OXYGEN_PA: f64 = 13_300.0;
const WATER_VAPOUR_PA: f64 = 6_280.0;
const CARBON_DIOXIDE_PA: f64 = 5_330.0;
// Standard gas entropies at 298 K, one bar, J/(mol K).
const OXYGEN_S: f64 = 205.2;
const WATER_S: f64 = 188.8;
const CARBON_DIOXIDE_S: f64 = 213.8;
// kJ/mol for one kT at body temperature.
const KT_KJ: f64 = 1.380_649e-23 * BODY_TEMPERATURE_K * AVOGADRO_PER_MOL / 1.0e3;

fn report_json(name: &str, report: &CompetitionReport) -> String {
    let rows: Vec<String> = report
        .guests
        .iter()
        .map(|g| {
            format!(
                "{{\"id\":\"{}\",\"free_energy_kj_per_mol\":{:.2},\"occupancy\":{:.6e}}}",
                g.id,
                g.binding_free_energy_j * AVOGADRO_PER_MOL / 1.0e3,
                g.occupancy
            )
        })
        .collect();
    let ratio = report
        .ratio("o2", "water")
        .map_or("null".to_string(), |r| format!("{r:.3e}"));
    format!(
        "{{\"site\":\"{name}\",\"empty\":{:.6e},\"o2_over_water\":{ratio},\"guests\":[{}]}}",
        report.empty,
        rows.join(",")
    )
}

fn main() {
    let metal = [
        Competitor::from_molar("o2", OXYGEN_PA, -95.6, OXYGEN_S),
        Competitor::from_molar("water", WATER_VAPOUR_PA, -10.3, WATER_S),
        Competitor::from_molar("co2", CARBON_DIOXIDE_PA, 0.6, CARBON_DIOXIDE_S),
    ];
    // Classical methyl pocket: O2 6.01 kT, water 7.78 kT; CO2 does not fit.
    let bare = [
        Competitor::from_molar("o2", OXYGEN_PA, -6.01 * KT_KJ, OXYGEN_S),
        Competitor::from_molar("water", WATER_VAPOUR_PA, -7.78 * KT_KJ, WATER_S),
    ];
    let sites = [
        report_json(
            "metal_fe_nh3_4_proxy",
            &competition_report(&metal, BODY_TEMPERATURE_K),
        ),
        report_json(
            "bare_methyl_pocket_4p5a",
            &competition_report(&bare, BODY_TEMPERATURE_K),
        ),
    ];
    println!(
        "{{\"temperature_k\":{BODY_TEMPERATURE_K},\"note\":\"estimates; gas entropy fully lost on binding\",\"sites\":[{}]}}",
        sites.join(",")
    );
}
