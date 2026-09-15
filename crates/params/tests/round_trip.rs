mod common;

use common::{good_record, quantity};
use nanocad_params::{Method, PartRecord, Quantity, Validation, SCHEMA, VERSION};

#[test]
fn json_round_trip_asserts_every_field() {
    let record = good_record();
    let json = record.to_json().unwrap();
    let back = PartRecord::from_json(&json).unwrap();

    assert_eq!(back, record);
    assert_eq!(back.schema, SCHEMA);
    assert_eq!(back.version, VERSION);
    assert_eq!(back.part_id, "gear.sun.diamond.2nm.12t");
    assert_eq!(back.geometry_ref, "sha256:deadbeef");
    assert_eq!(back.material, "diamond");
    assert_eq!(back.atoms, 1234);

    assert_eq!(back.mass_kg.value_si(), 1.5e-21);
    assert_eq!(back.mass_kg.unit(), "kg");
    assert_eq!(back.mass_kg.uncertainty_si(), Some(0.05e-21));

    assert_eq!(back.inertia_kg_m2[0].value_si(), 1.0e-40);
    assert_eq!(back.inertia_kg_m2[1].value_si(), 1.1e-40);
    assert_eq!(back.inertia_kg_m2[2].value_si(), 1.2e-40);
    assert_eq!(back.inertia_kg_m2[2].uncertainty_si(), None);

    assert_eq!(back.elastic_modulus_pa.value_si(), 1.05e12);
    assert_eq!(back.elastic_modulus_pa.uncertainty_si(), Some(0.08e12));
    assert_eq!(back.shear_modulus_pa.value_si(), 4.2e11);
    assert_eq!(back.poisson_ratio.value_si(), 0.2);
    assert_eq!(back.failure_stress_pa.value_si(), 2.0e10);
    assert_eq!(back.friction_coefficient.value_si(), 0.05);
    assert_eq!(back.thermal_conductivity_w_m_k.value_si(), 1000.0);
    assert_eq!(back.specific_heat_j_kg_k.value_si(), 500.0);

    assert_eq!(back.method, Method::Md);
    assert_eq!(back.validation, Validation::Unverified);

    assert_eq!(back.provenance.source, "nanocad-engine");
    assert_eq!(back.provenance.method, Method::Md);
    assert_eq!(back.provenance.code_version, "0.1.0");
    assert_eq!(back.provenance.force_field.as_deref(), Some("MM4"));
    assert_eq!(back.provenance.timestamp, "2026-09-13T00:00:00Z");
    assert_eq!(back.provenance.uncertainty, None);
    assert_eq!(back.provenance.validation, Validation::Unverified);
    assert_eq!(
        back.provenance.notes,
        "synthetic test record; not a measurement"
    );
}

#[test]
fn byte_round_trip_asserts_every_field() {
    let record = good_record();
    let bytes = record.to_bytes().unwrap();
    let back = PartRecord::from_bytes(&bytes).unwrap();
    assert_eq!(back, record);
}

#[test]
fn json_uses_the_schema_field_names() {
    let json = good_record().to_json().unwrap();
    assert!(json.contains("\"elastic_modulus_Pa\""));
    assert!(json.contains("\"shear_modulus_Pa\""));
    assert!(json.contains("\"failure_stress_Pa\""));
    assert!(json.contains("\"thermal_conductivity_W_m_K\""));
    assert!(json.contains("\"specific_heat_J_kg_K\""));
    assert!(json.contains("\"unit\":\"Pa\""));
    assert!(json.contains("\"uncertainty\":null"));
    assert!(json.contains("\"validation\":\"unverified\""));
}

#[test]
fn unknown_uncertainty_stays_null() {
    let json = serde_json::to_string(&quantity(1.2e-40, "kg*m^2", None)).unwrap();
    assert_eq!(
        json,
        r#"{"value":1.2e-40,"unit":"kg*m^2","uncertainty":null}"#
    );
}

#[test]
fn deserialization_rejects_a_foreign_version() {
    let mut record = good_record();
    record.version = 2;
    let json = record.to_json().unwrap();
    assert!(PartRecord::from_json(&json).is_err());
}

#[test]
fn deserialization_rejects_a_foreign_schema() {
    let mut record = good_record();
    record.schema = "nanocad.param.body".to_string();
    let json = record.to_json().unwrap();
    assert!(PartRecord::from_json(&json).is_err());
}

#[test]
fn deserialization_rejects_a_missing_required_field() {
    let json = r#"{"schema":"nanocad.param.part","version":1}"#;
    assert!(PartRecord::from_json(json).is_err());
}

#[test]
fn quantity_round_trips_a_non_si_unit() {
    let q: Quantity = quantity(2.5e-9, "nm", Some(1.0e-10));
    let json = serde_json::to_string(&q).unwrap();
    let back: Quantity = serde_json::from_str(&json).unwrap();
    assert_eq!(back.unit(), "nm");
    // A non-SI unit converts on the boundary, so allow one rounding step.
    assert!((back.value_si() - q.value_si()).abs() <= 1.0e-23);
    assert!((back.uncertainty_si().unwrap() - q.uncertainty_si().unwrap()).abs() <= 1.0e-25);
}
