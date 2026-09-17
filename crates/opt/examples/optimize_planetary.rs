//! Searches the planetary gear geometry and prints the result as JSON.
//!
//! Each argument is `key=value`. The keys are `population`, `generations`,
//! `steps` and `seed`. The search runs the real generator and the real
//! metrics for every candidate, so it is slow. Run it in release.
//!
//! ```sh
//! cargo run --release -p nanocad-opt --example optimize_planetary -- generations=6
//! ```
//!
//! The objective is the slip barrier in units of the thermal energy at 300 K.
//! The module and the layer count stay at their defaults, so the candidates
//! have a comparable size. A candidate whose clearance fails is rejected.

use std::error::Error;

use nanocad_opt::{minimize, CmaEsOptions, PlanetarySearch};

fn main() -> Result<(), Box<dyn Error>> {
    let mut options = CmaEsOptions::default();
    let mut search = PlanetarySearch::default();
    for argument in std::env::args().skip(1) {
        let (key, value) = argument
            .split_once('=')
            .ok_or_else(|| format!("expected key=value, got {argument}"))?;
        match key {
            "population" => options.population = value.parse()?,
            "generations" => options.generations = value.parse()?,
            "seed" => options.seed = value.parse()?,
            "steps" => search.slip_steps = value.parse()?,
            other => return Err(format!("unknown key {other}").into()),
        }
    }

    let initial = search.initial();
    let start = search.score(&initial);
    let report = minimize(&initial, &search.bounds, &options, &search);
    let best = search.score(&report.best_parameters);

    println!("{}", to_json(&search, &initial, &start, &report, &best));
    Ok(())
}

fn to_json(
    search: &PlanetarySearch,
    initial: &[f64],
    start: &Option<nanocad_opt::PlanetaryScore>,
    report: &nanocad_opt::CmaEsReport,
    best: &Option<nanocad_opt::PlanetaryScore>,
) -> String {
    let names: Vec<String> = search.names.iter().map(|n| format!("\"{n}\"")).collect();
    let initial_values: Vec<String> = initial.iter().map(|v| v.to_string()).collect();
    let best_values: Vec<String> = report
        .best_parameters
        .iter()
        .map(|v| v.to_string())
        .collect();
    let history: Vec<String> = report.history.iter().map(|v| v.to_string()).collect();
    let start_cost = start
        .as_ref()
        .map(|s| s.barrier_over_kt)
        .unwrap_or(f64::INFINITY);
    let (barrier_j, barrier_over_kt, clearance_m, atom_count) = match best {
        Some(score) => (
            score.barrier_j,
            score.barrier_over_kt,
            score.clearance_m,
            score.atom_count,
        ),
        None => (f64::INFINITY, f64::INFINITY, 0.0, 0),
    };
    format!(
        "{{\"names\":[{}],\"initial\":[{}],\"best_parameters\":[{}],\"cost\":{},\"start_cost\":{},\
\"barrier_j\":{},\"barrier_over_kt\":{},\"clearance_m\":{},\"atom_count\":{},\
\"evaluations\":{},\"generations\":{},\"rejected\":{},\"converged\":{},\"history\":[{}]}}",
        names.join(","),
        initial_values.join(","),
        best_values.join(","),
        report.best_cost,
        start_cost,
        barrier_j,
        barrier_over_kt,
        clearance_m,
        atom_count,
        report.evaluations,
        report.generations,
        report.rejected,
        report.converged,
        history.join(","),
    )
}
