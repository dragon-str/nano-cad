use serde::{Deserialize, Serialize};

use crate::error::ParamError;
use crate::provenance::{Method, Provenance, Validation};
use crate::quantity::Quantity;

/// The schema name of a [`PartRecord`]. This follows `PARAMETERS.md`.
pub const SCHEMA: &str = "nanocad.param.part";

/// The current `PartRecord` schema version.
pub const VERSION: u32 = 1;

/// The measured or computed properties of one part. It is the L1 output and the
/// L2 input.
///
/// Every quantity carries a unit and an uncertainty. The uncertainty is `null`
/// when it is unknown. The record carries one method, one validation status,
/// and one [`Provenance`] block. This follows the schema in `PARAMETERS.md`.
///
/// Field names in JSON follow the schema, for example `elastic_modulus_Pa`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PartRecord {
    /// The schema name. It must equal [`SCHEMA`].
    pub schema: String,
    /// The schema version. It must equal [`VERSION`].
    pub version: u32,
    /// A stable identifier for the part.
    pub part_id: String,
    /// A hash of the NCZ part geometry.
    pub geometry_ref: String,
    /// The material name.
    pub material: String,
    /// The number of atoms.
    pub atoms: u64,
    /// The total mass.
    pub mass_kg: Quantity,
    /// The three principal moments of inertia.
    pub inertia_kg_m2: [Quantity; 3],
    /// The elastic modulus.
    #[serde(rename = "elastic_modulus_Pa")]
    pub elastic_modulus_pa: Quantity,
    /// The shear modulus.
    #[serde(rename = "shear_modulus_Pa")]
    pub shear_modulus_pa: Quantity,
    /// The Poisson ratio. It is dimensionless.
    pub poisson_ratio: Quantity,
    /// The failure stress.
    #[serde(rename = "failure_stress_Pa")]
    pub failure_stress_pa: Quantity,
    /// The friction coefficient. It is dimensionless.
    pub friction_coefficient: Quantity,
    /// The thermal conductivity.
    #[serde(rename = "thermal_conductivity_W_m_K")]
    pub thermal_conductivity_w_m_k: Quantity,
    /// The specific heat.
    #[serde(rename = "specific_heat_J_kg_K")]
    pub specific_heat_j_kg_k: Quantity,
    /// The method that produced the record.
    pub method: Method,
    /// The validation status of the record.
    pub validation: Validation,
    /// The provenance of the record.
    pub provenance: Provenance,
}

impl PartRecord {
    /// Serialize the record to compact JSON.
    pub fn to_json(&self) -> Result<String, ParamError> {
        Ok(serde_json::to_string(self)?)
    }

    /// Serialize the record to indented JSON.
    pub fn to_json_pretty(&self) -> Result<String, ParamError> {
        Ok(serde_json::to_string_pretty(self)?)
    }

    /// Serialize the record to a JSON byte buffer.
    pub fn to_bytes(&self) -> Result<Vec<u8>, ParamError> {
        Ok(serde_json::to_vec(self)?)
    }

    /// Deserialize a record from JSON and check its versioned header.
    pub fn from_json(json: &str) -> Result<Self, ParamError> {
        let record: Self = serde_json::from_str(json)?;
        record.check_header()?;
        Ok(record)
    }

    /// Deserialize a record from a byte buffer and check its versioned header.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, ParamError> {
        let record: Self = serde_json::from_slice(bytes)?;
        record.check_header()?;
        Ok(record)
    }

    /// Reject a foreign schema name or an unsupported version.
    pub fn check_header(&self) -> Result<(), ParamError> {
        if self.schema != SCHEMA {
            return Err(ParamError::UnsupportedSchema {
                found: self.schema.clone(),
                expected: SCHEMA.to_string(),
            });
        }
        if self.version != VERSION {
            return Err(ParamError::UnsupportedVersion {
                found: self.version,
                expected: VERSION,
            });
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A minimal record for header tests. It is not physically plausible.
    fn header_only_record() -> PartRecord {
        let q = Quantity::derived(1.0, "1").unwrap();
        PartRecord {
            schema: SCHEMA.to_string(),
            version: VERSION,
            part_id: "p".to_string(),
            geometry_ref: "g".to_string(),
            material: "m".to_string(),
            atoms: 1,
            mass_kg: q.clone(),
            inertia_kg_m2: [q.clone(), q.clone(), q.clone()],
            elastic_modulus_pa: q.clone(),
            shear_modulus_pa: q.clone(),
            poisson_ratio: q.clone(),
            failure_stress_pa: q.clone(),
            friction_coefficient: q.clone(),
            thermal_conductivity_w_m_k: q.clone(),
            specific_heat_j_kg_k: q.clone(),
            method: Method::Estimate,
            validation: Validation::Unverified,
            provenance: Provenance {
                source: "test".to_string(),
                method: Method::Estimate,
                code_version: "0.1.0".to_string(),
                force_field: None,
                timestamp: "2026-09-13T00:00:00Z".to_string(),
                uncertainty: None,
                validation: Validation::Unverified,
                notes: String::new(),
            },
        }
    }

    #[test]
    fn good_header_passes() {
        assert!(header_only_record().check_header().is_ok());
    }

    #[test]
    fn foreign_schema_is_rejected() {
        let mut record = header_only_record();
        record.schema = "nanocad.param.joint".to_string();
        assert!(matches!(
            record.check_header(),
            Err(ParamError::UnsupportedSchema { .. })
        ));
    }

    #[test]
    fn foreign_version_is_rejected() {
        let mut record = header_only_record();
        record.version = 99;
        assert!(matches!(
            record.check_header(),
            Err(ParamError::UnsupportedVersion {
                found: 99,
                expected: 1
            })
        ));
    }

    #[test]
    fn json_round_trip_preserves_the_record() {
        let record = header_only_record();
        let json = record.to_json().unwrap();
        let back = PartRecord::from_json(&json).unwrap();
        assert_eq!(back, record);
    }

    #[test]
    fn byte_round_trip_preserves_the_record() {
        let record = header_only_record();
        let bytes = record.to_bytes().unwrap();
        let back = PartRecord::from_bytes(&bytes).unwrap();
        assert_eq!(back, record);
    }
}
