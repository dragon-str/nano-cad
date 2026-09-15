mod common;

use common::good_record;
use nanocad_params::{
    ParamError, ParameterLibrary, PartRecord, Validation, LIBRARY_SCHEMA, LIBRARY_VERSION,
};

/// A good record with a chosen id and geometry hash.
fn record_with(part_id: &str, geometry_ref: &str) -> PartRecord {
    let mut record = good_record();
    record.part_id = part_id.to_string();
    record.geometry_ref = geometry_ref.to_string();
    record
}

fn library_with_two_parts() -> ParameterLibrary {
    let mut library = ParameterLibrary::new();
    library
        .insert(record_with("gear.sun", "sha256:sun1"), 1)
        .unwrap();
    library
        .insert(record_with("gear.sun", "sha256:sun2"), 2)
        .unwrap();
    library
        .insert(record_with("gear.planet", "sha256:planet1"), 1)
        .unwrap();
    library
}

#[test]
fn insert_and_retrieve() {
    let library = library_with_two_parts();
    assert_eq!(library.len(), 3);
    assert_eq!(library.revision(), 3);
    assert!(!library.is_empty());

    let latest = library.get("gear.sun").unwrap();
    assert_eq!(latest.geometry_ref, "sha256:sun2");
    assert_eq!(latest.part_id, "gear.sun");

    assert_eq!(
        library.get("gear.planet").unwrap().geometry_ref,
        "sha256:planet1"
    );
    assert!(library.get("missing").is_none());

    let ids: Vec<&str> = library.part_ids().collect();
    assert_eq!(ids, vec!["gear.planet", "gear.sun"]);
    assert_eq!(library.list().len(), 3);
}

#[test]
fn lookup_by_version() {
    let library = library_with_two_parts();
    assert_eq!(
        library.get_version("gear.sun", 1).unwrap().geometry_ref,
        "sha256:sun1"
    );
    assert_eq!(
        library.get_version("gear.sun", 2).unwrap().geometry_ref,
        "sha256:sun2"
    );
    assert!(library.get_version("gear.sun", 3).is_none());
    assert_eq!(library.latest_version("gear.sun"), Some(2));
    assert_eq!(library.latest_version("gear.planet"), Some(1));
    assert_eq!(library.latest_version("missing"), None);
}

#[test]
fn duplicate_record_is_rejected() {
    let mut library = ParameterLibrary::new();
    library
        .insert(record_with("gear.sun", "sha256:sun1"), 1)
        .unwrap();

    let error = library
        .insert(record_with("gear.sun", "sha256:sun2"), 1)
        .unwrap_err();
    assert!(matches!(
        error,
        ParamError::DuplicateRecord {
            part_id,
            version: 1
        } if part_id == "gear.sun"
    ));
    // The failed insert changed nothing.
    assert_eq!(library.len(), 1);
    assert_eq!(library.revision(), 1);
}

#[test]
fn geometry_mismatch_is_rejected_on_insert() {
    let mut library = ParameterLibrary::new();
    let error = library
        .insert_with_geometry(record_with("gear.sun", "sha256:sun1"), 1, "sha256:other")
        .unwrap_err();
    assert!(matches!(
        error,
        ParamError::GeometryMismatch { version: 1, .. }
    ));
    assert!(library.is_empty());
    assert_eq!(library.revision(), 0);

    library
        .insert_with_geometry(record_with("gear.sun", "sha256:sun1"), 1, "sha256:sun1")
        .unwrap();
    assert_eq!(library.revision(), 1);
}

#[test]
fn geometry_mismatch_is_flagged_by_lookup() {
    let library = library_with_two_parts();
    assert!(library
        .verify_geometry("gear.sun", 1, "sha256:sun1")
        .is_ok());
    assert!(matches!(
        library.verify_geometry("gear.sun", 1, "sha256:wrong"),
        Err(ParamError::GeometryMismatch { .. })
    ));
    assert!(matches!(
        library.verify_geometry("missing", 1, "sha256:sun1"),
        Err(ParamError::UnknownPart { .. })
    ));
}

#[test]
fn persist_to_bytes_and_reload() {
    let library = library_with_two_parts();
    let bytes = library.to_bytes().unwrap();
    let reloaded = ParameterLibrary::from_bytes(&bytes).unwrap();
    assert_eq!(reloaded, library);
    assert_eq!(reloaded.revision(), 3);
    assert_eq!(reloaded.len(), 3);

    // Provenance and validation status survive the reload.
    let record = reloaded.get_version("gear.sun", 2).unwrap();
    assert_eq!(record.provenance.source, "nanocad-engine");
    assert_eq!(record.provenance.code_version, "0.1.0");
    assert_eq!(record.validation, Validation::Unverified);
}

#[test]
fn persist_to_file_and_reload() {
    let library = library_with_two_parts();
    let path = std::env::temp_dir().join(format!(
        "nanocad-params-library-{}-{}.json",
        std::process::id(),
        library.revision()
    ));
    library.save_to_file(&path).unwrap();
    let reloaded = ParameterLibrary::load_from_file(&path).unwrap();
    assert_eq!(reloaded, library);
    let _ = std::fs::remove_file(&path);
}

#[test]
fn reload_rejects_a_foreign_library_schema() {
    let library = library_with_two_parts();
    let json = library
        .to_json()
        .unwrap()
        .replace(LIBRARY_SCHEMA, "nanocad.param.other");
    assert!(matches!(
        ParameterLibrary::from_json(&json),
        Err(ParamError::UnsupportedLibrarySchema { .. })
    ));
}

#[test]
fn reload_rejects_a_foreign_library_version() {
    let library = library_with_two_parts();
    let mut value: serde_json::Value = serde_json::from_str(&library.to_json().unwrap()).unwrap();
    value["version"] = serde_json::Value::from(LIBRARY_VERSION + 1);
    let json = serde_json::to_string(&value).unwrap();
    assert!(matches!(
        ParameterLibrary::from_json(&json),
        Err(ParamError::UnsupportedLibraryVersion { .. })
    ));
}

#[test]
fn reload_rejects_a_bad_record_header() {
    let library = library_with_two_parts();
    let mut value: serde_json::Value = serde_json::from_str(&library.to_json().unwrap()).unwrap();
    value["records"]["gear.planet"][0]["record"]["version"] = serde_json::Value::from(99);
    let json = serde_json::to_string(&value).unwrap();
    assert!(matches!(
        ParameterLibrary::from_json(&json),
        Err(ParamError::UnsupportedVersion { .. })
    ));
}

#[test]
fn remove_part_and_version() {
    let mut library = library_with_two_parts();
    assert!(library.remove("gear.planet"));
    assert!(!library.remove("gear.planet"));
    assert_eq!(library.revision(), 4);
    assert_eq!(library.len(), 2);

    let removed = library.remove_version("gear.sun", 1).unwrap();
    assert_eq!(removed.geometry_ref, "sha256:sun1");
    assert!(library.get_version("gear.sun", 1).is_none());
    assert_eq!(library.get("gear.sun").unwrap().geometry_ref, "sha256:sun2");

    assert!(library.remove_version("gear.sun", 2).is_some());
    assert!(library.is_empty());
    assert_eq!(library.part_ids().count(), 0);
}

#[test]
fn empty_library_json_round_trip() {
    let library = ParameterLibrary::new();
    let json = library.to_json().unwrap();
    assert!(json.contains(LIBRARY_SCHEMA));
    let back = ParameterLibrary::from_json(&json).unwrap();
    assert_eq!(back, library);
    assert_eq!(back.revision(), 0);
}
