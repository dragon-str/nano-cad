//! Writes the as-built drain bore of the hull wall (M28-35) for the drain
//! transit MD, `scripts/drain_transit_md.py`.
//!
//! The wall is relaxed as in the module scene, or left as built with `raw`
//! (the frozen as-built channel of ADR-0098; a vacuum relaxation folds the
//! hydroxyl lining into the empty bore). A square patch of the wall
//! around the drain axis is cut out, with each atom's charge, and turned so
//! the bore runs along z: (x, y, z) -> (y - y_drain, z, x).
//!
//! Usage: `drain_patch_json [half_width_m] [path] [raw|relaxed]`.

use nanocad_jigs::{module_parts, module_relax_options};
use nanocad_meter::relax_valence;

const ELEMENTARY_CHARGE_C: f64 = 1.602_176_634e-19;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let half_m: f64 = args.first().and_then(|v| v.parse().ok()).unwrap_or(0.6e-9);
    let path = args
        .get(1)
        .cloned()
        .unwrap_or_else(|| "drain.json".to_string());
    let raw = args.get(2).is_some_and(|mode| mode == "raw");
    let mut parts = match module_parts() {
        Ok(parts) => parts,
        Err(error) => {
            eprintln!("the module failed: {error}");
            std::process::exit(1);
        }
    };
    let wall = &mut parts.housing;
    let (positions, converged) = if raw {
        let positions: Vec<[f64; 3]> = (0..wall.atom_count())
            .map(|i| wall.topology.position_m(i).unwrap_or([0.0; 3]))
            .collect();
        (positions, true)
    } else {
        match relax_valence(
            &wall.topology,
            0..wall.atom_count(),
            &module_relax_options(),
        ) {
            Ok(report) => (report.positions_m, report.converged),
            Err(error) => {
                eprintln!("the wall did not relax: {error}");
                std::process::exit(1);
            }
        }
    };
    let drain_y = parts.wall.drain_y_m;
    let mut atoms = Vec::new();
    let mut thickness_m: f64 = 0.0;
    for (i, p) in positions.iter().enumerate() {
        if (p[1] - drain_y).abs() >= half_m || p[2].abs() >= half_m {
            continue;
        }
        let Some(element) = wall.topology.element(i) else {
            continue;
        };
        thickness_m = thickness_m.max(2.0 * p[0].abs());
        let q = wall.topology.charge_c(i).unwrap_or(0.0) / ELEMENTARY_CHARGE_C;
        atoms.push(format!(
            "[\"{}\",{:.5},{:.5},{:.5},{:.5}]",
            element.symbol().unwrap_or("X"),
            (p[1] - drain_y) * 1e10,
            p[2] * 1e10,
            p[0] * 1e10,
            q
        ));
    }
    let json = format!(
        "{{\"units\":\"angstrom\",\"radius\":{:.4},\"thickness\":{:.4},\"converged\":{},\"atoms\":[{}]}}",
        nanocad_parts::HullWallSpec::default().drain_radius_m * 1e10,
        thickness_m * 1e10,
        converged,
        atoms.join(",")
    );
    if let Err(error) = std::fs::write(&path, json) {
        eprintln!("the patch was not written: {error}");
        std::process::exit(1);
    }
    println!("{} atoms, relaxed {}", atoms.len(), converged);
}
