//! Prints the metric score of a planetary gear set as JSON.
//!
//! Each argument is `key=value`, for example `module_m=5e-10`. Parameters that
//! are not listed take their defaults. Unknown keys are refused.
//!
//! ```sh
//! cargo run -q -p nanocad-meter --example score_json -- sun_teeth=12
//! ```

use std::error::Error;

use nanocad_meter::{score_planetary, Fidelity, Score};
use nanocad_parts::{ParameterSet, PlanetaryGenerator};

/// The sun rate that the viewer and the scene use, in radians per second.
const SUN_RAD_PER_S: f64 = -0.34;

fn main() -> Result<(), Box<dyn Error>> {
    let mut parameters = ParameterSet::new();
    for argument in std::env::args().skip(1) {
        let (key, value) = argument
            .split_once('=')
            .ok_or_else(|| format!("expected key=value, got {argument}"))?;
        parameters = parameters.with(key, value.parse()?);
    }
    let set = PlanetaryGenerator.build(&parameters)?;
    let score = score_planetary(&set, SUN_RAD_PER_S)?;
    println!("{}", to_json(&score));
    Ok(())
}

fn to_json(score: &Score) -> String {
    let items: Vec<String> = score
        .values
        .iter()
        .map(|value| {
            format!(
                "{{\"name\":\"{}\",\"value\":{},\"unit\":\"{}\",\"fidelity\":\"{}\",\"note\":\"{}\"}}",
                escape(&value.name),
                value.value,
                escape(&value.unit),
                fidelity_name(value.fidelity),
                escape(&value.note),
            )
        })
        .collect();
    format!("{{\"metrics\":[{}]}}", items.join(","))
}

fn fidelity_name(fidelity: Fidelity) -> &'static str {
    match fidelity {
        Fidelity::Geometric => "geometric",
        Fidelity::QuasiStatic => "quasi_static",
        Fidelity::Harmonic => "harmonic",
        Fidelity::Dynamics => "dynamics",
    }
}

fn escape(text: &str) -> String {
    text.replace('\\', "\\\\").replace('"', "\\\"")
}
