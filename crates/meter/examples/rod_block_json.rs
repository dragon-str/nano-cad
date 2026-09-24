//! The site-block test of the rod pocket (ADR-0090, test 1).
//!
//! The rod face moves from below the iron to past the mouth. At each face
//! position the example searches oxygen and water at the iron, with the
//! coordination term, the wall and the rod. It prints one JSON object with
//! the energy profile, the blocked fraction of poses, the least rod-to-block
//! distance and the steepest rise of the oxygen energy, which is the peak
//! push force of a quasi-static stroke.
//!
//! With one argument it also writes the loaded-pocket geometry (rod pulled
//! back) as JSON for the explicit-water test.

use nanocad_meter::{iron_coordination, metal_site_binding, MetalSiteTarget, AVOGADRO_PER_MOL};
use nanocad_parts::{guest, ParameterSet, RodPocketGenerator};

/// The retracted rod face, below the iron, in metres.
const RETRACT_M: f64 = 1.5e-10;

fn main() {
    // POCKET_CLEARANCE_M and POCKET_ROD_RADIUS_M override the defaults, for
    // a scan.
    let mut parameters = ParameterSet::new();
    for (variable, name) in [
        ("POCKET_CLEARANCE_M", "clearance_m"),
        ("POCKET_ROD_RADIUS_M", "rod_radius_m"),
        ("POCKET_IRON_RANK", "iron_rank"),
        ("POCKET_DEPTH_M", "depth_m"),
    ] {
        if let Some(value) = std::env::var(variable).ok().and_then(|v| v.parse().ok()) {
            parameters = parameters.with(name, value);
        }
    }
    let set = match RodPocketGenerator.build(&parameters) {
        Ok(set) => set,
        Err(error) => {
            eprintln!("the pocket failed: {error}");
            std::process::exit(1);
        }
    };
    let target = MetalSiteTarget {
        directions: 32,
        spins: 6,
        donor_tilts: 2,
        ..MetalSiteTarget::default()
    };
    let per_mol = |j: f64| j * AVOGADRO_PER_MOL / 1.0e3;
    let start_m = set.iron_z_m - 4.0e-10;
    let end_m = set.mouth_z_m + 1.0e-10;
    let step_m = 0.5e-10;
    let steps = ((end_m - start_m) / step_m).round() as usize;
    let mut rows = Vec::new();
    let mut profile = Vec::new();
    for step in 0..=steps {
        let face_m = start_m + step_m * step as f64;
        let Ok(part) = set.assembled(face_m) else {
            eprintln!("the assembly failed");
            std::process::exit(1);
        };
        let mut cells = Vec::new();
        let mut oxygen_j = 0.0;
        for id in ["o2", "water"] {
            let coordination = iron_coordination(id);
            let Some(report) = metal_site_binding(
                &part,
                set.iron,
                guest(id).expect("guest"),
                coordination.as_ref(),
                &target,
            ) else {
                eprintln!("the search failed");
                std::process::exit(1);
            };
            // A guest with no allowed pose, or no bound pose, is free: zero.
            let energy_j = report.minimum_energy_j.min(0.0);
            if id == "o2" {
                oxygen_j = energy_j;
            }
            cells.push(format!(
                "\"{id}_kj_per_mol\":{:.2},\"{id}_blocked\":{:.3},\"{id}_donor_m\":{:.3e}",
                per_mol(energy_j),
                report.blocked as f64 / report.poses.max(1) as f64,
                report.donor_distance_m
            ));
        }
        // Least distance from a rod atom to a block atom.
        let block = set.block_atom_count();
        let positions = part.topology.positions_m();
        let mut nearest = f64::INFINITY;
        let mut nearest_iron = f64::INFINITY;
        for r in positions[3 * block..].chunks(3) {
            for (index, b) in positions[..3 * block].chunks(3).enumerate() {
                let d =
                    ((r[0] - b[0]).powi(2) + (r[1] - b[1]).powi(2) + (r[2] - b[2]).powi(2)).sqrt();
                nearest = nearest.min(d);
                if index == set.iron {
                    nearest_iron = nearest_iron.min(d);
                }
            }
        }
        profile.push((face_m, oxygen_j));
        rows.push(format!(
            "{{\"face_above_iron_m\":{:.2e},{},\"rod_to_block_m\":{:.3e},\"rod_to_iron_m\":{:.3e}}}",
            face_m - set.iron_z_m,
            cells.join(","),
            nearest,
            nearest_iron
        ));
    }
    let peak_force_n = profile
        .windows(2)
        .map(|w| (w[1].1 - w[0].1) / (w[1].0 - w[0].0))
        .fold(0.0_f64, f64::max);

    if let Some(path) = std::env::args().nth(1) {
        let retract_m = std::env::var("POCKET_RETRACT_M")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(RETRACT_M);
        let Ok(part) = set.assembled(set.iron_z_m - retract_m) else {
            std::process::exit(1);
        };
        let atoms: Vec<String> = (0..part.atom_count())
            .filter_map(|i| {
                let p = part.topology.position_m(i)?;
                let e = part.topology.element(i)?;
                Some(format!(
                    "[\"{}\",{:.5},{:.5},{:.5}]",
                    e.symbol().unwrap_or("X"),
                    p[0] * 1e10,
                    p[1] * 1e10,
                    p[2] * 1e10
                ))
            })
            .collect();
        let json = format!(
            "{{\"units\":\"angstrom\",\"iron\":{},\"bore_radius\":{:.4},\"rod_radius\":{:.4},\
             \"iron_z\":{:.4},\"mouth_z\":{:.4},\"face_z\":{:.4},\"block_atoms\":{},\"atoms\":[{}]}}",
            set.iron,
            set.bore_radius_m * 1e10,
            set.rod_radius_m * 1e10,
            set.iron_z_m * 1e10,
            set.mouth_z_m * 1e10,
            (set.iron_z_m - retract_m) * 1e10,
            set.block_atom_count(),
            atoms.join(",")
        );
        if let Err(error) = std::fs::write(&path, json) {
            eprintln!("the geometry was not written: {error}");
            std::process::exit(1);
        }
    }

    println!(
        "{{\"bore_radius_m\":{:.3e},\"rod_radius_m\":{:.3e},\"iron_z_m\":{:.3e},\"mouth_z_m\":{:.3e},\
         \"peak_push_force_n\":{:.3e},\"profile\":[{}]}}",
        set.bore_radius_m,
        set.rod_radius_m,
        set.iron_z_m,
        set.mouth_z_m,
        peak_force_n,
        rows.join(",")
    );
}
