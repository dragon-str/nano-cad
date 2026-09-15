use serde::{Deserialize, Serialize};

use crate::quantity::Quantity;

/// How a value or a record was produced.
///
/// This follows Rule 5 of `PARAMETERS.md`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Method {
    /// Molecular dynamics.
    Md,
    /// Density functional theory.
    Dft,
    /// A physical experiment.
    Experiment,
    /// A fit to other data.
    Fit,
    /// An estimate.
    Estimate,
}

/// The validation status of a record.
///
/// This follows Rule 6 of `PARAMETERS.md`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Validation {
    /// Not yet checked.
    #[serde(rename = "unverified")]
    Unverified,
    /// Checked against an independent source.
    #[serde(rename = "cross-checked")]
    CrossChecked,
    /// Formally certified.
    #[serde(rename = "certified")]
    Certified,
}

/// Who or what produced a value, how, and when.
///
/// This follows the `Provenance` block in `PARAMETERS.md`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Provenance {
    /// A tool, a person, or a dataset.
    pub source: String,
    /// The method that produced the value.
    pub method: Method,
    /// The version of the code that produced the value.
    pub code_version: String,
    /// The force field, when the method used one.
    #[serde(default)]
    pub force_field: Option<String>,
    /// The ISO-8601 time of production.
    pub timestamp: String,
    /// The uncertainty of the record, or [`None`] when it is unknown.
    #[serde(default)]
    pub uncertainty: Option<Quantity>,
    /// The validation status.
    pub validation: Validation,
    /// Free-form notes.
    #[serde(default)]
    pub notes: String,
}

/// Builds the provenance block for one extraction run.
///
/// The source is the engine, the code version is this crate, and the validation
/// status starts unverified. An extraction never certifies its own output.
pub(crate) fn extraction_provenance(method: Method, notes: String) -> Provenance {
    Provenance {
        source: "nanocad-engine".to_string(),
        method,
        code_version: env!("CARGO_PKG_VERSION").to_string(),
        force_field: None,
        timestamp: utc_timestamp_now(),
        uncertainty: None,
        validation: Validation::Unverified,
        notes,
    }
}

/// The current UTC time as an ISO-8601 instant with second resolution.
///
/// This avoids a date-library dependency. The civil-date arithmetic is the
/// public-domain algorithm of Howard Hinnant, transcribed here.
pub(crate) fn utc_timestamp_now() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let seconds_since_epoch = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|duration| duration.as_secs() as i64)
        .unwrap_or(0);
    let days = seconds_since_epoch.div_euclid(86_400);
    let second_of_day = seconds_since_epoch.rem_euclid(86_400);
    let (year, month, day) = civil_from_days(days);
    format!(
        "{year:04}-{month:02}-{day:02}T{:02}:{:02}:{:02}Z",
        second_of_day / 3_600,
        (second_of_day % 3_600) / 60,
        second_of_day % 60
    )
}

fn civil_from_days(days_since_epoch: i64) -> (i64, u32, u32) {
    let shifted = days_since_epoch + 719_468;
    let era = shifted.div_euclid(146_097);
    let day_of_era = (shifted - era * 146_097) as u64;
    let year_of_era =
        (day_of_era - day_of_era / 1_460 + day_of_era / 36_524 - day_of_era / 146_096) / 365;
    let year = year_of_era as i64 + era * 400;
    let day_of_year = day_of_era - (365 * year_of_era + year_of_era / 4 - year_of_era / 100);
    let month_prime = (5 * day_of_year + 2) / 153;
    let day = (day_of_year - (153 * month_prime + 2) / 5 + 1) as u32;
    let month = if month_prime < 10 {
        month_prime + 3
    } else {
        month_prime - 9
    } as u32;
    (year + i64::from(month <= 2), month, day)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn method_serializes_lowercase() {
        assert_eq!(serde_json::to_string(&Method::Dft).unwrap(), "\"dft\"");
        assert_eq!(
            serde_json::to_string(&Method::Estimate).unwrap(),
            "\"estimate\""
        );
    }

    #[test]
    fn validation_serializes_with_hyphen() {
        assert_eq!(
            serde_json::to_string(&Validation::CrossChecked).unwrap(),
            "\"cross-checked\""
        );
        assert_eq!(
            serde_json::to_string(&Validation::Unverified).unwrap(),
            "\"unverified\""
        );
        assert_eq!(
            serde_json::to_string(&Validation::Certified).unwrap(),
            "\"certified\""
        );
    }

    #[test]
    fn validation_round_trips() {
        for status in [
            Validation::Unverified,
            Validation::CrossChecked,
            Validation::Certified,
        ] {
            let json = serde_json::to_string(&status).unwrap();
            let back: Validation = serde_json::from_str(&json).unwrap();
            assert_eq!(back, status);
        }
    }

    #[test]
    fn the_timestamp_has_the_iso8601_shape() {
        let timestamp = utc_timestamp_now();
        assert_eq!(timestamp.len(), 20);
        assert!(timestamp.ends_with('Z'));
        assert_eq!(timestamp.as_bytes()[10], b'T');
        assert!(timestamp[..4].chars().all(|c| c.is_ascii_digit()));
    }
}
