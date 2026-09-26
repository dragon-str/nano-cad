//! Writes two parallel hydrogen-capped diamond plates as JSON, for the
//! gas-space test of option A (M28-05). The plates are `PlateGenerator`
//! blocks with reconstructed (001) faces. The gap is the distance between the
//! two facing carbon planes.
//!
//! Usage: `slot_json <gap_m> <path>`.

use nanocad_model::Element;
use nanocad_parts::{ParameterSet, PartGenerator, PlateGenerator};

fn main() {
    let mut args = std::env::args().skip(1);
    let gap_m: f64 = args.next().and_then(|v| v.parse().ok()).unwrap_or(1.5e-9);
    let path = args.next().unwrap_or_else(|| "slot.json".to_string());
    let plate = match PlateGenerator.generate(
        &ParameterSet::new()
            .with("size_x_m", 3.2e-9)
            .with("size_y_m", 3.2e-9)
            .with("thickness_m", 0.7e-9)
            .with("hole_radius_m", 0.0),
    ) {
        Ok(plate) => plate,
        Err(error) => {
            eprintln!("the plate failed: {error}");
            std::process::exit(1);
        }
    };
    let topology = &plate.topology;
    let carbons_z: Vec<f64> = (0..plate.atom_count())
        .filter(|&i| topology.element(i) == Some(Element::CARBON))
        .filter_map(|i| topology.position_m(i).map(|p| p[2]))
        .collect();
    let top_m = carbons_z.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let bottom_m = carbons_z.iter().copied().fold(f64::INFINITY, f64::min);
    let mut atoms = Vec::new();
    // Lower plate: its top carbon plane at z = -gap/2. Upper plate: its bottom
    // carbon plane at z = +gap/2.
    for (shift_m, name) in [(-0.5 * gap_m - top_m, 0), (0.5 * gap_m - bottom_m, 1)] {
        for i in 0..plate.atom_count() {
            let (Some(p), Some(e)) = (topology.position_m(i), topology.element(i)) else {
                continue;
            };
            atoms.push(format!(
                "[\"{}\",{:.5},{:.5},{:.5},{name}]",
                e.symbol().unwrap_or("X"),
                p[0] * 1e10,
                p[1] * 1e10,
                (p[2] + shift_m) * 1e10
            ));
        }
    }
    let json = format!(
        "{{\"units\":\"angstrom\",\"gap\":{:.4},\"atoms\":[{}]}}",
        gap_m * 1e10,
        atoms.join(",")
    );
    if let Err(error) = std::fs::write(&path, json) {
        eprintln!("the slot was not written: {error}");
        std::process::exit(1);
    }
    println!("{} atoms, gap {:.2} A", atoms.len(), gap_m * 1e10);
}
