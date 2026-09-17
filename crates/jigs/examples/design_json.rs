//! Prints a design document with a two-step edit, an undo and a redo.
//!
//! The example is small. It shows the history rules: a commit appends, an
//! identical commit is refused, undo and redo move the cursor, and a new commit
//! after an undo drops the redo tail.

use nanocad_jigs::{DesignDocument, DesignSnapshot, Measurement, PartRecord};

fn snapshot(sun_teeth: f64, barrier_kt: f64) -> DesignSnapshot {
    DesignSnapshot::new()
        .with_parameter("sun_teeth", sun_teeth)
        .with_parameter("planet_teeth", 9.0)
        .with_part(PartRecord {
            name: "planetary".to_string(),
            generator_id: "planetary".to_string(),
            atom_count: 142091,
            bond_count: 173341,
            material: "diamondoid".to_string(),
        })
        .with_measurement(Measurement::new(
            "slip_barrier",
            barrier_kt * 4.142e-21,
            "J",
            "quasi_static",
            format!("{barrier_kt:.1} kT at 300 K"),
        ))
}

fn main() {
    let mut document = DesignDocument::new(snapshot(12.0, 13.8));
    let appended = document.commit(snapshot(13.0, 12.1));
    let refused = document.commit(snapshot(13.0, 12.1));
    let undone = document.undo();
    let redone = document.redo();
    println!(
        "appended={appended} refused={refused} undone={undone} redone={redone} depth={} cursor={}",
        document.depth(),
        document.cursor()
    );
    match document.to_json_pretty() {
        Ok(text) => println!("{text}"),
        Err(error) => eprintln!("write failed: {error}"),
    }
}
