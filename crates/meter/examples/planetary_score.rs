//! Prints the clearance score of the default planetary gear set.
//!
//! Run it from the repository root:
//!
//! ```sh
//! cargo run -p nanocad-meter --example planetary_score
//! ```

use std::error::Error;

use nanocad_meter::{Clearance, MovingAtoms};
use nanocad_parts::{ParameterSet, PlanetaryGenerator};

fn main() -> Result<(), Box<dyn Error>> {
    let parameters = ParameterSet::new().with("planet_count", 3.0);
    let set = PlanetaryGenerator.build(&parameters)?;
    let moving = MovingAtoms::planetary(&set, &set.design, -0.34);
    let clearance = Clearance::default();
    let report = clearance.measure(&moving);
    let value = report.to_metric_value(&clearance.target);

    println!("atoms {}", moving.positions_m.len());
    println!("period_s {:.6}", moving.period_s);
    println!(
        "least distance {:.4e} m between {} and {} at {:.4} s over {} samples",
        report.minimum_m, report.body_a, report.body_b, report.time_s, report.samples
    );
    println!("{} [{:.2e} m]: {}", value.name, value.value, value.note);
    Ok(())
}
