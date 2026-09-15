mod common;

use common::{good_record, quantity};
use nanocad_params::{check, is_consistent, Method, Problem, Validation};

#[test]
fn good_record_has_no_problems() {
    let problems = check(&good_record());
    assert!(problems.is_empty(), "unexpected problems: {problems:?}");
    assert!(is_consistent(&good_record()));
}

#[test]
fn bad_schema_is_reported() {
    let mut record = good_record();
    record.schema = "nanocad.param.body".to_string();
    assert!(check(&record)
        .iter()
        .any(|p| matches!(p, Problem::SchemaMismatch { .. })));
}

#[test]
fn bad_version_is_reported() {
    let mut record = good_record();
    record.version = 7;
    assert!(check(&record).iter().any(|p| matches!(
        p,
        Problem::UnsupportedVersion {
            found: 7,
            supported: 1
        }
    )));
}

#[test]
fn missing_required_fields_are_reported() {
    let mut record = good_record();
    record.part_id = String::new();
    record.geometry_ref = "   ".to_string();
    record.material = String::new();
    record.atoms = 0;
    let problems = check(&record);
    assert!(problems
        .iter()
        .any(|p| matches!(p, Problem::MissingField { field: "part_id" })));
    assert!(problems.iter().any(|p| matches!(
        p,
        Problem::MissingField {
            field: "geometry_ref"
        }
    )));
    assert!(problems
        .iter()
        .any(|p| matches!(p, Problem::MissingField { field: "material" })));
    assert!(problems
        .iter()
        .any(|p| matches!(p, Problem::Implausible { field: "atoms", .. })));
}

#[test]
fn wrong_dimension_is_reported() {
    let mut record = good_record();
    record.elastic_modulus_pa = quantity(1.05e12, "s", Some(0.08e12));
    assert!(check(&record).iter().any(|p| matches!(
        p,
        Problem::DimensionMismatch {
            field: "elastic_modulus_Pa",
            ..
        }
    )));
}

#[test]
fn negative_uncertainty_is_reported() {
    let mut record = good_record();
    record.mass_kg = quantity(1.5e-21, "kg", Some(-0.05e-21));
    assert!(check(&record).iter().any(|p| matches!(
        p,
        Problem::NegativeUncertainty {
            field: "mass_kg",
            ..
        }
    )));
}

#[test]
fn null_uncertainty_is_allowed() {
    let mut record = good_record();
    record.mass_kg = quantity(1.5e-21, "kg", None);
    assert!(is_consistent(&record));
}

#[test]
fn poisson_ratio_out_of_range_is_reported() {
    for bad in [0.5, 0.9, -1.0, -1.5] {
        let mut record = good_record();
        record.poisson_ratio = quantity(bad, "1", None);
        assert!(
            check(&record).iter().any(|p| matches!(
                p,
                Problem::Implausible {
                    field: "poisson_ratio",
                    ..
                }
            )),
            "poisson_ratio {bad} must be reported"
        );
    }
}

#[test]
fn nonpositive_modulus_is_reported() {
    for bad in [0.0, -1.0] {
        let mut record = good_record();
        record.elastic_modulus_pa = quantity(bad, "Pa", None);
        assert!(
            check(&record).iter().any(|p| matches!(
                p,
                Problem::Implausible {
                    field: "elastic_modulus_Pa",
                    ..
                }
            )),
            "elastic modulus {bad} must be reported"
        );
    }
}

#[test]
fn nonpositive_mass_and_inertia_are_reported() {
    let mut record = good_record();
    record.mass_kg = quantity(-1.0, "kg", None);
    record.inertia_kg_m2[0] = quantity(0.0, "kg*m^2", None);
    let problems = check(&record);
    assert!(problems.iter().any(|p| matches!(
        p,
        Problem::Implausible {
            field: "mass_kg",
            ..
        }
    )));
    assert!(problems.iter().any(|p| matches!(
        p,
        Problem::Implausible {
            field: "inertia_kg_m2[0]",
            ..
        }
    )));
}

#[test]
fn empty_provenance_is_reported() {
    let mut record = good_record();
    record.provenance.source = String::new();
    record.provenance.code_version = String::new();
    record.provenance.timestamp = "yesterday".to_string();
    let problems = check(&record);
    assert!(problems.iter().any(|p| matches!(
        p,
        Problem::MissingField {
            field: "provenance.source"
        }
    )));
    assert!(problems.iter().any(|p| matches!(
        p,
        Problem::MissingField {
            field: "provenance.code_version"
        }
    )));
    assert!(problems.iter().any(|p| matches!(
        p,
        Problem::Implausible {
            field: "provenance.timestamp",
            ..
        }
    )));
}

#[test]
fn provenance_method_and_validation_must_match_the_record() {
    let mut record = good_record();
    record.method = Method::Dft;
    record.validation = Validation::Certified;
    let problems = check(&record);
    assert!(problems.iter().any(|p| matches!(
        p,
        Problem::Implausible {
            field: "provenance.method",
            ..
        }
    )));
    assert!(problems.iter().any(|p| matches!(
        p,
        Problem::Implausible {
            field: "provenance.validation",
            ..
        }
    )));
}

#[test]
fn check_returns_every_problem_at_once() {
    let mut record = good_record();
    record.part_id = String::new();
    record.poisson_ratio = quantity(0.9, "1", None);
    record.shear_modulus_pa = quantity(-1.0, "Pa", None);
    assert!(check(&record).len() >= 3);
}
