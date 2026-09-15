//! Golden round-trip tests for the format crate.
//!
//! These tests are wired as a `[[test]]` target of `nanocad-format`, so
//! `cargo test -p nanocad-format` runs them. They compare exact `f64` bit
//! patterns, not approximate distances.
//!
//! The NanoEngineer sample lives outside this repository. It is GPL code, so we
//! never copy it in. The test reads it in place and skips when it is absent.

use std::path::{Path, PathBuf};

use nanocad_format::{export_mmp, from_ncz_bytes, import_mmp, to_ncz_bytes};
use nanocad_model::Document;

const REAL_MMP: &str = "/tmp/nanoengineer/cad/partlib/couplings/Universal Joint.mmp";

fn synthetic_mmp_path() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/synthetic_chain.mmp")
}

fn assert_same_coordinates_bits(left: &Document, right: &Document, label: &str) {
    assert_eq!(left.part_count(), right.part_count(), "{label}: part count");
    for (part_index, (left_part, right_part)) in left.parts.iter().zip(&right.parts).enumerate() {
        assert_eq!(
            left_part.atom_count(),
            right_part.atom_count(),
            "{label}: part {part_index} atom count"
        );
        assert_eq!(
            left_part.bond_count(),
            right_part.bond_count(),
            "{label}: part {part_index} bond count"
        );
        assert_eq!(
            left_part.topology.bonds().collect::<Vec<_>>(),
            right_part.topology.bonds().collect::<Vec<_>>(),
            "{label}: part {part_index} topology"
        );
        for atom_index in 0..left_part.atom_count() {
            let left_position = left_part
                .topology
                .position_m(atom_index)
                .expect("left position");
            let right_position = right_part
                .topology
                .position_m(atom_index)
                .expect("right position");
            for axis in 0..3 {
                assert_eq!(
                    left_position[axis].to_bits(),
                    right_position[axis].to_bits(),
                    "{label}: part {part_index} atom {atom_index} axis {axis}"
                );
            }
        }
    }
}

#[test]
fn synthetic_mmp_round_trips_through_mmp() {
    let path = synthetic_mmp_path();
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|error| panic!("cannot read {}: {error}", path.display()));
    let document = import_mmp(&text).expect("import synthetic mmp");
    assert!(document.part_count() >= 1, "the fixture has a part");
    let part = document.part(0).expect("part");
    assert_eq!(part.atom_count(), 4, "the fixture has four atoms");
    assert_eq!(part.bond_count(), 3, "the fixture has three bonds");

    let exported = export_mmp(&document).expect("export mmp");
    let reimported = import_mmp(&exported).expect("reimport mmp");
    assert_same_coordinates_bits(&document, &reimported, "mmp export");
    assert_eq!(document, reimported, "mmp export keeps every field");
}

#[test]
fn synthetic_mmp_round_trips_through_ncz() {
    let path = synthetic_mmp_path();
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|error| panic!("cannot read {}: {error}", path.display()));
    let document = import_mmp(&text).expect("import synthetic mmp");

    let bytes = to_ncz_bytes(&document).expect("write ncz");
    let reloaded = from_ncz_bytes(&bytes).expect("read ncz");
    assert_same_coordinates_bits(&document, &reloaded, "ncz reload");
    assert_eq!(document, reloaded, "ncz reload keeps every field");
}

#[test]
fn real_mmp_round_trips_through_ncz_bit_for_bit() {
    let path = Path::new(REAL_MMP);
    if !path.exists() {
        eprintln!(
            "SKIP: the real MMP file is absent at {REAL_MMP}. \
             Put the NanoEngineer tree at /tmp/nanoengineer to run this test."
        );
        return;
    }

    let text = std::fs::read_to_string(path).expect("read real mmp");
    let document = import_mmp(&text).expect("import real mmp");
    assert!(document.part_count() >= 1, "the real file has a part");
    let part = document.part(0).expect("part");
    assert!(part.atom_count() > 1000, "the real file has many atoms");
    assert!(part.bond_count() > 1000, "the real file has many bonds");

    let bytes = to_ncz_bytes(&document).expect("write ncz");
    let reloaded = from_ncz_bytes(&bytes).expect("read ncz");
    assert_same_coordinates_bits(&document, &reloaded, "real mmp ncz reload");
    assert_eq!(document, reloaded, "the real file reloads exactly");
    eprintln!(
        "RAN: {} atoms and {} bonds reproduced bit for bit via NCZ.",
        part.atom_count(),
        part.bond_count()
    );
}
