//! Prints the steered-drive report of a planetary gear set as JSON.
//!
//! Each argument is `key=value`, for example `sun_teeth=12`. Parameters that
//! are not listed take their defaults. The keys `steps`, `dt_s`, `cluster_m`,
//! `max_atoms`, `cutoff_m`, `contact_m` and `temperature_k` set the drive
//! target. Unknown keys are refused.
//!
//! ```sh
//! cargo run -q -p nanocad-meter --example drive_json -- steps=600
//! ```

use std::error::Error;

use nanocad_meter::{DriveTarget, Fidelity, SteeredDrive};
use nanocad_parts::{ParameterSet, PlanetaryGenerator};

/// The sun rate that the viewer and the scene use, in radians per second.
const SUN_RAD_PER_S: f64 = -0.34;

fn main() -> Result<(), Box<dyn Error>> {
    let mut parameters = ParameterSet::new();
    let mut target = DriveTarget::default();
    for argument in std::env::args().skip(1) {
        let (key, value) = argument
            .split_once('=')
            .ok_or_else(|| format!("expected key=value, got {argument}"))?;
        match key {
            "steps" => target.steps = value.parse()?,
            "dt_s" => target.dt_s = value.parse()?,
            "cluster_m" => target.cluster_m = value.parse()?,
            "max_atoms" => target.max_atoms = value.parse()?,
            "cutoff_m" => target.cutoff_m = value.parse()?,
            "contact_m" => target.contact_m = value.parse()?,
            "temperature_k" => target.temperature_k = value.parse()?,
            _ => parameters = parameters.with(key, value.parse()?),
        }
    }
    let set = PlanetaryGenerator.build(&parameters)?;
    let report = SteeredDrive::with_target(target).measure(&set, SUN_RAD_PER_S);
    let metric = report.to_metric_value();
    println!(
        "{{\"name\":\"{}\",\"value\":{},\"unit\":\"{}\",\"fidelity\":\"{}\",\"note\":\"{}\"}}",
        metric.name,
        metric.value,
        metric.unit,
        fidelity_name(metric.fidelity),
        metric.note.replace('\\', "\\\\").replace('"', "\\\"")
    );
    println!(
        "{{\"steps\":{},\"atom_count\":{},\"free_atom_count\":{},\"torque_n_m\":{},\"peak_torque_n_m\":{},\"work_j\":{},\"energy_loss_j\":{},\"temperature_rise_k\":{},\"final_temperature_k\":{}}}",
        report.steps,
        report.atom_count,
        report.free_atom_count,
        report.torque_n_m,
        report.peak_torque_n_m,
        report.work_j,
        report.energy_loss_j,
        report.temperature_rise_k,
        report.final_temperature_k,
    );
    Ok(())
}

fn fidelity_name(fidelity: Fidelity) -> &'static str {
    match fidelity {
        Fidelity::Geometric => "geometric",
        Fidelity::QuasiStatic => "quasi_static",
        Fidelity::Harmonic => "harmonic",
        Fidelity::Dynamics => "dynamics",
    }
}
