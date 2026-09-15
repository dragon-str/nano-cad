use core::fmt;

use nanocad_units::Dimension;

use crate::quantity::Quantity;
use crate::record::{PartRecord, SCHEMA, VERSION};

const DIMENSIONLESS: Dimension = Dimension::new(0, 0, 0, 0, 0, 0, 0);
const MASS: Dimension = Dimension::new(0, 1, 0, 0, 0, 0, 0);
const MOMENT_OF_INERTIA: Dimension = Dimension::new(2, 1, 0, 0, 0, 0, 0);
const PRESSURE: Dimension = Dimension::new(-1, 1, -2, 0, 0, 0, 0);
const THERMAL_CONDUCTIVITY: Dimension = Dimension::new(1, 1, -3, 0, -1, 0, 0);
const SPECIFIC_HEAT: Dimension = Dimension::new(2, 0, -2, 0, -1, 0, 0);

const INERTIA_FIELDS: [&str; 3] = ["inertia_kg_m2[0]", "inertia_kg_m2[1]", "inertia_kg_m2[2]"];

/// One problem that a consistency check found.
///
/// [`check`] returns a list of these. It never returns a single boolean, so
/// the caller sees every problem in one pass.
#[derive(Clone, Debug, PartialEq)]
pub enum Problem {
    /// The schema name is not [`SCHEMA`].
    SchemaMismatch { found: String },
    /// The schema version is not [`VERSION`].
    UnsupportedVersion { found: u32, supported: u32 },
    /// A required string field is empty.
    MissingField { field: &'static str },
    /// A unit symbol is not in the registry.
    UnknownUnit { field: &'static str, unit: String },
    /// The unit has the wrong dimension for the field.
    DimensionMismatch {
        field: &'static str,
        unit: String,
        expected: Dimension,
    },
    /// An uncertainty is negative or not a number. Rule 3 allows a null, not a
    /// negative value.
    NegativeUncertainty { field: &'static str, value_si: f64 },
    /// A value is outside the physically plausible range.
    Implausible { field: &'static str, reason: String },
}

impl fmt::Display for Problem {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Problem::SchemaMismatch { found } => {
                write!(f, "schema {found:?} does not match {SCHEMA:?}")
            }
            Problem::UnsupportedVersion { found, supported } => {
                write!(f, "version {found} is not supported; expected {supported}")
            }
            Problem::MissingField { field } => write!(f, "required field {field} is empty"),
            Problem::UnknownUnit { field, unit } => {
                write!(f, "field {field} has unknown unit {unit:?}")
            }
            Problem::DimensionMismatch {
                field,
                unit,
                expected,
            } => write!(
                f,
                "field {field} has unit {unit:?} with dimension {expected:?}, which is wrong"
            ),
            Problem::NegativeUncertainty { field, value_si } => {
                write!(
                    f,
                    "field {field} has uncertainty {value_si}, which is not nonnegative"
                )
            }
            Problem::Implausible { field, reason } => {
                write!(f, "field {field} is not plausible: {reason}")
            }
        }
    }
}

/// Check one parameter record. Return every problem that was found.
///
/// The checks are:
///
/// - the schema name and the version match this crate,
/// - the required fields are present,
/// - each unit has the correct dimension for its field,
/// - each uncertainty is nonnegative or null,
/// - each value is physically plausible.
pub fn check(record: &PartRecord) -> Vec<Problem> {
    let mut problems = Vec::new();

    if record.schema != SCHEMA {
        problems.push(Problem::SchemaMismatch {
            found: record.schema.clone(),
        });
    }
    if record.version != VERSION {
        problems.push(Problem::UnsupportedVersion {
            found: record.version,
            supported: VERSION,
        });
    }

    require_text("part_id", &record.part_id, &mut problems);
    require_text("geometry_ref", &record.geometry_ref, &mut problems);
    require_text("material", &record.material, &mut problems);
    if record.atoms == 0 {
        problems.push(Problem::Implausible {
            field: "atoms",
            reason: "the atom count must be positive".to_string(),
        });
    }

    check_dimension("mass_kg", &record.mass_kg, MASS, &mut problems);
    for (index, quantity) in record.inertia_kg_m2.iter().enumerate() {
        let field = INERTIA_FIELDS
            .get(index)
            .copied()
            .unwrap_or("inertia_kg_m2");
        check_dimension(field, quantity, MOMENT_OF_INERTIA, &mut problems);
    }
    check_dimension(
        "elastic_modulus_Pa",
        &record.elastic_modulus_pa,
        PRESSURE,
        &mut problems,
    );
    check_dimension(
        "shear_modulus_Pa",
        &record.shear_modulus_pa,
        PRESSURE,
        &mut problems,
    );
    check_dimension(
        "poisson_ratio",
        &record.poisson_ratio,
        DIMENSIONLESS,
        &mut problems,
    );
    check_dimension(
        "failure_stress_Pa",
        &record.failure_stress_pa,
        PRESSURE,
        &mut problems,
    );
    check_dimension(
        "friction_coefficient",
        &record.friction_coefficient,
        DIMENSIONLESS,
        &mut problems,
    );
    check_dimension(
        "thermal_conductivity_W_m_K",
        &record.thermal_conductivity_w_m_k,
        THERMAL_CONDUCTIVITY,
        &mut problems,
    );
    check_dimension(
        "specific_heat_J_kg_K",
        &record.specific_heat_j_kg_k,
        SPECIFIC_HEAT,
        &mut problems,
    );

    let quantities = [
        ("mass_kg", &record.mass_kg),
        ("inertia_kg_m2[0]", &record.inertia_kg_m2[0]),
        ("inertia_kg_m2[1]", &record.inertia_kg_m2[1]),
        ("inertia_kg_m2[2]", &record.inertia_kg_m2[2]),
        ("elastic_modulus_Pa", &record.elastic_modulus_pa),
        ("shear_modulus_Pa", &record.shear_modulus_pa),
        ("poisson_ratio", &record.poisson_ratio),
        ("failure_stress_Pa", &record.failure_stress_pa),
        ("friction_coefficient", &record.friction_coefficient),
        (
            "thermal_conductivity_W_m_K",
            &record.thermal_conductivity_w_m_k,
        ),
        ("specific_heat_J_kg_K", &record.specific_heat_j_kg_k),
    ];
    for (field, quantity) in quantities {
        check_uncertainty(field, quantity, &mut problems);
    }
    if let Some(uncertainty) = &record.provenance.uncertainty {
        check_uncertainty("provenance.uncertainty", uncertainty, &mut problems);
    }

    check_positive("mass_kg", &record.mass_kg, &mut problems);
    for (index, quantity) in record.inertia_kg_m2.iter().enumerate() {
        let field = INERTIA_FIELDS
            .get(index)
            .copied()
            .unwrap_or("inertia_kg_m2");
        check_positive(field, quantity, &mut problems);
    }
    check_positive(
        "elastic_modulus_Pa",
        &record.elastic_modulus_pa,
        &mut problems,
    );
    check_positive("shear_modulus_Pa", &record.shear_modulus_pa, &mut problems);
    check_positive(
        "failure_stress_Pa",
        &record.failure_stress_pa,
        &mut problems,
    );
    check_positive(
        "thermal_conductivity_W_m_K",
        &record.thermal_conductivity_w_m_k,
        &mut problems,
    );
    check_positive(
        "specific_heat_J_kg_K",
        &record.specific_heat_j_kg_k,
        &mut problems,
    );
    check_nonnegative(
        "friction_coefficient",
        &record.friction_coefficient,
        &mut problems,
    );

    let poisson = record.poisson_ratio.value_si();
    if !poisson.is_finite() || poisson <= -1.0 || poisson >= 0.5 {
        problems.push(Problem::Implausible {
            field: "poisson_ratio",
            reason: format!("the Poisson ratio must be in (-1, 0.5), got {poisson}"),
        });
    }

    require_text(
        "provenance.source",
        &record.provenance.source,
        &mut problems,
    );
    require_text(
        "provenance.code_version",
        &record.provenance.code_version,
        &mut problems,
    );
    if !looks_like_iso8601(&record.provenance.timestamp) {
        problems.push(Problem::Implausible {
            field: "provenance.timestamp",
            reason: format!(
                "the timestamp {:?} is not an ISO-8601 instant",
                record.provenance.timestamp
            ),
        });
    }
    if record.provenance.method != record.method {
        problems.push(Problem::Implausible {
            field: "provenance.method",
            reason: "the provenance method must match the record method".to_string(),
        });
    }
    if record.provenance.validation != record.validation {
        problems.push(Problem::Implausible {
            field: "provenance.validation",
            reason: "the provenance validation must match the record validation".to_string(),
        });
    }

    problems
}

/// True when [`check`] finds no problem.
pub fn is_consistent(record: &PartRecord) -> bool {
    check(record).is_empty()
}

fn require_text(field: &'static str, value: &str, problems: &mut Vec<Problem>) {
    if value.trim().is_empty() {
        problems.push(Problem::MissingField { field });
    }
}

fn check_dimension(
    field: &'static str,
    quantity: &Quantity,
    expected: Dimension,
    problems: &mut Vec<Problem>,
) {
    match quantity.dimension() {
        None => problems.push(Problem::UnknownUnit {
            field,
            unit: quantity.unit().to_string(),
        }),
        Some(found) if found != expected => problems.push(Problem::DimensionMismatch {
            field,
            unit: quantity.unit().to_string(),
            expected,
        }),
        Some(_) => {}
    }
}

fn check_uncertainty(field: &'static str, quantity: &Quantity, problems: &mut Vec<Problem>) {
    if let Some(value_si) = quantity.uncertainty_si() {
        if !value_si.is_finite() || value_si < 0.0 {
            problems.push(Problem::NegativeUncertainty { field, value_si });
        }
    }
}

fn check_positive(field: &'static str, quantity: &Quantity, problems: &mut Vec<Problem>) {
    let value_si = quantity.value_si();
    if !value_si.is_finite() || value_si <= 0.0 {
        problems.push(Problem::Implausible {
            field,
            reason: format!("the value must be finite and positive, got {value_si}"),
        });
    }
}

fn check_nonnegative(field: &'static str, quantity: &Quantity, problems: &mut Vec<Problem>) {
    let value_si = quantity.value_si();
    if !value_si.is_finite() || value_si < 0.0 {
        problems.push(Problem::Implausible {
            field,
            reason: format!("the value must be finite and nonnegative, got {value_si}"),
        });
    }
}

/// A light ISO-8601 check. It pins the `YYYY-MM-DDThh:mm:ss` prefix.
fn looks_like_iso8601(value: &str) -> bool {
    let bytes = value.as_bytes();
    if bytes.len() < 19 {
        return false;
    }
    let digit = |index: usize| bytes.get(index).is_some_and(u8::is_ascii_digit);
    digit(0)
        && digit(1)
        && digit(2)
        && digit(3)
        && bytes.get(4) == Some(&b'-')
        && digit(5)
        && digit(6)
        && bytes.get(7) == Some(&b'-')
        && digit(8)
        && digit(9)
        && bytes.get(10) == Some(&b'T')
        && digit(11)
        && digit(12)
        && bytes.get(13) == Some(&b':')
        && digit(14)
        && digit(15)
        && bytes.get(16) == Some(&b':')
        && digit(17)
        && digit(18)
}
