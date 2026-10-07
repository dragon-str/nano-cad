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
    let spec = CarrierSiteSpec::default();
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

    let json = format!(
        "{{\"units\":\"angstrom\",\"iron\":{iron},\"bore_radius\":{:.4},\
         \"face_z\":{face_z:.4},\"mouth_z\":{mouth_z:.4},\"block_atoms\":{block_atoms},\
         \"axis\":[0.0,0.0,1.0],\"atoms\":[{}]}}",
        spec.pocket_radius_m * 1.0e10,
        atoms.join(",")
    );
    if let Some(path) = std::env::args().nth(1) {
        if let Err(error) = std::fs::write(&path, &json) {
            eprintln!("cannot write {path}: {error}");
            std::process::exit(1);
        }
        println!("wrote {path}: {block_atoms} scaffold atoms, iron {iron}");
    } else {
        println!("{json}");
    }
}
