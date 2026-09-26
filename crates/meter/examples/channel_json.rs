//! Writes a lined through-channel in a diamond slab as JSON, with the atom
//! charges, for the drain-channel MD test (M28-09).
//!
//! Usage: `channel_json <wall_group> <radius_m> <thickness_m> <path>`.

use nanocad_parts::{BindingPocketGenerator, ParameterSet, PartGenerator};

const ELEMENTARY_CHARGE_C: f64 = 1.602_176_634e-19;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let group: f64 = args.first().and_then(|v| v.parse().ok()).unwrap_or(0.0);
    let radius_m: f64 = args.get(1).and_then(|v| v.parse().ok()).unwrap_or(3.8e-10);
    let thickness_m: f64 = args.get(2).and_then(|v| v.parse().ok()).unwrap_or(1.0e-9);
    let path = args.get(3).cloned().unwrap_or_else(|| "channel.json".to_string());
    let part = match BindingPocketGenerator.generate(
        &ParameterSet::new()
            .with("pocket_radius_m", radius_m)
            .with("pocket_depth_m", thickness_m)
            .with("thickness_m", thickness_m)
            .with("wall_m", 1.0e-9)
            .with("wall_group", group),
    ) {
        Ok(part) => part,
        Err(error) => {
            eprintln!("the channel failed: {error}");
            std::process::exit(1);
        }
    };
    let atoms: Vec<String> = (0..part.atom_count())
        .filter_map(|i| {
            let p = part.topology.position_m(i)?;
            let e = part.topology.element(i)?;
            let q = part.topology.charge_c(i).unwrap_or(0.0) / ELEMENTARY_CHARGE_C;
            Some(format!(
                "[\"{}\",{:.5},{:.5},{:.5},{:.5}]",
                e.symbol().unwrap_or("X"),
                p[0] * 1e10,
                p[1] * 1e10,
                p[2] * 1e10,
                q
            ))
        })
        .collect();
    let json = format!(
        "{{\"units\":\"angstrom\",\"radius\":{:.4},\"thickness\":{:.4},\"atoms\":[{}]}}",
        radius_m * 1e10,
        thickness_m * 1e10,
        atoms.join(",")
    );
    if let Err(error) = std::fs::write(&path, json) {
        eprintln!("the channel was not written: {error}");
        std::process::exit(1);
    }
    println!("{} atoms", atoms.len());
}
