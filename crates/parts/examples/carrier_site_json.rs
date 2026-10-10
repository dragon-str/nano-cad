//! Writes the carrier-site geometry as JSON for the explicit-water test.
//!
//! The pocket axis is +z: the porphine plane holds the iron, and the O2 pocket
//! and the mouth sit above it. The scaffold atoms come first, then the site
//! atoms, so `block_atoms` separates them. `axis` is written because the iron
//! has four nitrogen neighbours, not three carbons, so the water test cannot
//! infer the open side from the bonds.
//!
//! Usage:
//!
//!     cargo run --release -p nanocad-parts --example carrier_site_json -- site.json

use nanocad_parts::{CarrierSiteGenerator, CarrierSiteSpec};

fn main() {
    // Optional second argument: the ring slot radius in angstrom (ADR-0137).
    let mut spec = CarrierSiteSpec::default();
    if let Some(r) = std::env::args().nth(2).and_then(|a| a.parse::<f64>().ok()) {
        spec.slot_radius_m = r * 1.0e-10;
    }
    // Optional third argument: the slot half-thickness in angstrom.
    if let Some(h) = std::env::args().nth(3).and_then(|a| a.parse::<f64>().ok()) {
        spec.slot_half_m = h * 1.0e-10;
    }
    let site = match CarrierSiteGenerator.build_spec(&spec) {
        Ok(site) => site,
        Err(error) => {
            eprintln!("the carrier site failed: {error}");
            std::process::exit(1);
        }
    };
    let part = match site.merged_part() {
        Ok(part) => part,
        Err(error) => {
            eprintln!("the merge failed: {error}");
            std::process::exit(1);
        }
    };
    let block_atoms = site.scaffold.atom_count();
    let tethers: Vec<String> = site
        .tethers
        .iter()
        .map(|t| {
            format!(
                "{{\"bridge_atoms\":{},\"strain_kj_mol\":{:.2}}}",
                t.bridge_atoms, t.strain_kj_mol
            )
        })
        .collect();
    let iron = site.merged_iron();
    let face_z = spec.slot_half_m * 1.0e10;
    let mouth_z = (spec.slot_half_m + spec.pocket_height_m + spec.mouth_length_m) * 1.0e10;

    let atoms: Vec<String> = (0..part.atom_count())
        .filter_map(|i| {
            let p = part.topology.position_m(i)?;
            let e = part.topology.element(i)?;
            Some(format!(
                "[\"{}\",{:.5},{:.5},{:.5}]",
                e.symbol().unwrap_or("X"),
                p[0] * 1.0e10,
                p[1] * 1.0e10,
                p[2] * 1.0e10
            ))
        })
        .collect();

    // The topology bonds, so a consumer need not guess them from distances.
    let bonds: Vec<String> = (0..part.topology.bond_count())
        .filter_map(|b| {
            Some(format!(
                "[{},{},{}]",
                part.topology.bond_u(b)?,
                part.topology.bond_v(b)?,
                part.topology.bond_order(b)?
            ))
        })
        .collect();

    let donor = match &site.donor {
        Some(d) => format!(
            "{{\"oxygen\":{},\"hydrogen\":{},\"fe_distance\":{:.3}}}",
            d.oxygen,
            d.hydrogen,
            d.fe_distance_m * 1.0e10
        ),
        None => "null".to_owned(),
    };
    let catalyst = match &site.catalyst {
        Some(c) => format!(
            "{{\"manganese\":{},\"fe_distance\":{:.3}}}",
            c.manganese,
            c.fe_distance_m * 1.0e10
        ),
        None => "null".to_owned(),
    };

    let json = format!(
        "{{\"units\":\"angstrom\",\"iron\":{iron},\"bore_radius\":{:.4},\
         \"face_z\":{face_z:.4},\"mouth_z\":{mouth_z:.4},\"block_atoms\":{block_atoms},\
         \"axis\":[0.0,0.0,1.0],\"slot_radius\":{:.4},\"tethers\":[{}],\
         \"donor\":{donor},\"catalyst\":{catalyst},\"atoms\":[{}],\"bonds\":[{}]}}",
        spec.pocket_radius_m * 1.0e10,
        spec.slot_radius_m * 1.0e10,
        tethers.join(","),
        atoms.join(","),
        bonds.join(",")
    );
    if let Some(path) = std::env::args().nth(1) {
        if let Err(error) = std::fs::write(&path, &json) {
            eprintln!("cannot write {path}: {error}");
            std::process::exit(1);
        }
        let strain: f64 = site.tethers.iter().map(|t| t.strain_kj_mol).sum();
        let steric: f64 = site.tethers.iter().map(|t| t.steric_kj_mol).sum();
        println!(
            "wrote {path}: {block_atoms} scaffold atoms, iron {iron}, slot {:.2} A half {:.2} A, tethers {:?}, valence strain {strain:.1}, steric {steric:.1} kJ/mol, donor {donor}, catalyst {catalyst}",
            spec.slot_radius_m * 1.0e10,
            spec.slot_half_m * 1.0e10,
            site.tethers.iter().map(|t| t.bridge_atoms).collect::<Vec<_>>()
        );
    } else {
        println!("{json}");
    }
}
