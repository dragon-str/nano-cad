//! Conversion helpers shared by the binding modules.
//!
//! Every Rust error becomes a Python `ValueError`. Library code never panics,
//! so a caller sees a typed Python exception instead of an abort.

use nanocad_model::BondType;
use nanocad_params::{Method, Validation};
use nanocad_units::Unit;
use pyo3::exceptions::PyValueError;
use pyo3::PyErr;

/// Converts any displayable error into a Python `ValueError`.
pub fn err<E: std::fmt::Display>(error: E) -> PyErr {
    PyValueError::new_err(error.to_string())
}

/// Converts a 3-vector argument, which is length checked.
pub fn vec3(values: Vec<f64>) -> Result<[f64; 3], PyErr> {
    if values.len() != 3 {
        return Err(PyValueError::new_err(format!(
            "expected 3 values, got {}",
            values.len()
        )));
    }
    Ok([values[0], values[1], values[2]])
}

/// Converts a 4-vector argument, which is length checked.
pub fn vec4(values: Vec<f64>) -> Result<[f64; 4], PyErr> {
    if values.len() != 4 {
        return Err(PyValueError::new_err(format!(
            "expected 4 values, got {}",
            values.len()
        )));
    }
    Ok([values[0], values[1], values[2], values[3]])
}

/// Converts a 3x3 matrix argument, which is shape checked.
pub fn mat3(rows: Vec<Vec<f64>>) -> Result<[[f64; 3]; 3], PyErr> {
    if rows.len() != 3 || rows.iter().any(|row| row.len() != 3) {
        return Err(PyValueError::new_err(
            "expected 3 rows of 3 values for the inertia tensor",
        ));
    }
    Ok([
        [rows[0][0], rows[0][1], rows[0][2]],
        [rows[1][0], rows[1][1], rows[1][2]],
        [rows[2][0], rows[2][1], rows[2][2]],
    ])
}

/// Every unit symbol this crate accepts.
pub const UNIT_SYMBOLS: [&str; 14] = [
    "m", "kg", "s", "A", "K", "mol", "cd", "Å", "nm", "kcal/mol", "fs", "ps", "u", "e",
];

/// Parses a unit symbol.
pub fn unit_from_symbol(symbol: &str) -> Result<Unit, PyErr> {
    Ok(match symbol {
        "m" | "metre" | "meter" => Unit::Metre,
        "kg" | "kilogram" => Unit::Kilogram,
        "s" | "second" => Unit::Second,
        "A" | "ampere" => Unit::Ampere,
        "K" | "kelvin" => Unit::Kelvin,
        "mol" | "mole" => Unit::Mole,
        "cd" | "candela" => Unit::Candela,
        "Å" | "Aangstrom" | "angstrom" => Unit::Angstrom,
        "nm" | "nanometre" | "nanometer" => Unit::Nanometre,
        "kcal/mol" | "kcal_per_mol" => Unit::KilocaloriePerMole,
        "fs" | "femtosecond" => Unit::Femtosecond,
        "ps" | "picosecond" => Unit::Picosecond,
        "u" | "amu" => Unit::AtomicMassUnit,
        "e" | "electron_charge" => Unit::ElectronCharge,
        other => {
            return Err(PyValueError::new_err(format!(
                "unknown unit symbol {other:?}"
            )))
        }
    })
}

/// Parses a bond-type name.
pub fn bond_type_from_str(name: &str) -> Result<BondType, PyErr> {
    Ok(match name.to_ascii_lowercase().as_str() {
        "single" => BondType::Single,
        "double" => BondType::Double,
        "triple" => BondType::Triple,
        "aromatic" => BondType::Aromatic,
        "unknown" => BondType::Unknown,
        other => {
            return Err(PyValueError::new_err(format!(
                "unknown bond type {other:?}"
            )))
        }
    })
}

/// Returns the bond-type name used on the Python side.
pub fn bond_type_name(bond_type: BondType) -> &'static str {
    match bond_type {
        BondType::Single => "single",
        BondType::Double => "double",
        BondType::Triple => "triple",
        BondType::Aromatic => "aromatic",
        BondType::Unknown => "unknown",
    }
}

/// Parses a provenance method name.
pub fn method_from_str(name: &str) -> Result<Method, PyErr> {
    Ok(match name.to_ascii_lowercase().as_str() {
        "md" => Method::Md,
        "dft" => Method::Dft,
        "experiment" => Method::Experiment,
        "fit" => Method::Fit,
        "estimate" => Method::Estimate,
        other => {
            return Err(PyValueError::new_err(format!(
                "unknown method {other:?}; expected md, dft, experiment, fit, or estimate"
            )))
        }
    })
}

/// Returns the method name used on the Python side.
pub fn method_name(method: Method) -> &'static str {
    match method {
        Method::Md => "md",
        Method::Dft => "dft",
        Method::Experiment => "experiment",
        Method::Fit => "fit",
        Method::Estimate => "estimate",
    }
}

/// Parses a validation-status name.
pub fn validation_from_str(name: &str) -> Result<Validation, PyErr> {
    Ok(match name.to_ascii_lowercase().as_str() {
        "unverified" => Validation::Unverified,
        "cross-checked" | "cross_checked" => Validation::CrossChecked,
        "certified" => Validation::Certified,
        other => {
            return Err(PyValueError::new_err(format!(
            "unknown validation status {other:?}; expected unverified, cross-checked, or certified"
        )))
        }
    })
}

/// Returns the validation-status name used on the Python side.
pub fn validation_name(validation: Validation) -> &'static str {
    match validation {
        Validation::Unverified => "unverified",
        Validation::CrossChecked => "cross-checked",
        Validation::Certified => "certified",
    }
}
