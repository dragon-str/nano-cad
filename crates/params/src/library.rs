use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::ParamError;
use crate::record::PartRecord;

/// The schema name of a [`ParameterLibrary`]. This is separate from the
/// [`PartRecord`](crate::PartRecord) schema, because the library header is a
/// different contract.
pub const LIBRARY_SCHEMA: &str = "nanocad.param.library";

/// The current [`ParameterLibrary`] schema version.
pub const LIBRARY_VERSION: u32 = 1;

/// One stored record together with its part revision.
///
/// A part id can hold several revisions. Each revision has its own version
/// number. This is how the library keeps history for one part.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct LibraryEntry {
    /// The part revision number. It is unique within one part id.
    pub version: u32,
    /// The stored record. It keeps its own provenance and validation status.
    pub record: PartRecord,
}

/// A versioned store of many [`PartRecord`] values.
///
/// The library is keyed by `part_id`. Each part id can hold several revisions,
/// keyed by a `u32` version. The library itself has a schema version and a
/// revision counter. The counter increases on every change.
///
/// The library does not change a record. Every stored record keeps its
/// provenance and its validation status.
///
/// # Errors
///
/// Library methods return [`ParamError`] instead of panicking. An insert
/// rejects a duplicate `(part_id, version)` pair. An insert with a geometry
/// check rejects a `geometry_ref` that does not match the supplied hash.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ParameterLibrary {
    schema: String,
    version: u32,
    revision: u64,
    records: BTreeMap<String, Vec<LibraryEntry>>,
}

impl Default for ParameterLibrary {
    fn default() -> Self {
        Self::new()
    }
}

impl ParameterLibrary {
    /// Create an empty library at revision zero.
    pub fn new() -> Self {
        Self {
            schema: LIBRARY_SCHEMA.to_string(),
            version: LIBRARY_VERSION,
            revision: 0,
            records: BTreeMap::new(),
        }
    }

    /// The library schema name. It equals [`LIBRARY_SCHEMA`].
    pub fn schema(&self) -> &str {
        &self.schema
    }

    /// The library schema version. It equals [`LIBRARY_VERSION`].
    pub fn version(&self) -> u32 {
        self.version
    }

    /// The library revision counter. It increases on every change.
    pub fn revision(&self) -> u64 {
        self.revision
    }

    /// The number of stored records, across every part id and revision.
    pub fn len(&self) -> usize {
        self.records.values().map(Vec::len).sum()
    }

    /// True when the library holds no record.
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// The distinct part ids, in sorted order.
    pub fn part_ids(&self) -> impl Iterator<Item = &str> + '_ {
        self.records.keys().map(String::as_str)
    }

    /// Store a record at a part revision.
    ///
    /// The key is the `part_id` of the record. The key pair
    /// `(part_id, version)` must be new. A repeat is
    /// [`ParamError::DuplicateRecord`].
    pub fn insert(&mut self, record: PartRecord, version: u32) -> Result<(), ParamError> {
        let part_id = record.part_id.clone();
        let entries = self.records.entry(part_id.clone()).or_default();
        if entries.iter().any(|entry| entry.version == version) {
            return Err(ParamError::DuplicateRecord { part_id, version });
        }
        entries.push(LibraryEntry { version, record });
        entries.sort_by_key(|entry| entry.version);
        self.revision += 1;
        Ok(())
    }

    /// Store a record only when its `geometry_ref` equals the supplied hash.
    ///
    /// A mismatch is [`ParamError::GeometryMismatch`]. A mismatch changes
    /// nothing, so the revision does not increase.
    pub fn insert_with_geometry(
        &mut self,
        record: PartRecord,
        version: u32,
        expected_geometry_ref: &str,
    ) -> Result<(), ParamError> {
        if record.geometry_ref != expected_geometry_ref {
            return Err(ParamError::GeometryMismatch {
                part_id: record.part_id,
                version,
                expected: expected_geometry_ref.to_string(),
                found: record.geometry_ref,
            });
        }
        self.insert(record, version)
    }

    /// Get the newest revision of one part.
    pub fn get(&self, part_id: &str) -> Option<&PartRecord> {
        let version = self.latest_version(part_id)?;
        self.get_version(part_id, version)
    }

    /// Get one exact part revision.
    pub fn get_version(&self, part_id: &str, version: u32) -> Option<&PartRecord> {
        self.records
            .get(part_id)?
            .iter()
            .find(|entry| entry.version == version)
            .map(|entry| &entry.record)
    }

    /// The newest revision number of one part, or [`None`] when absent.
    pub fn latest_version(&self, part_id: &str) -> Option<u32> {
        self.records
            .get(part_id)?
            .iter()
            .map(|entry| entry.version)
            .max()
    }

    /// Every stored record, ordered by part id and then by version.
    pub fn list(&self) -> Vec<&PartRecord> {
        self.records
            .values()
            .flat_map(|entries| entries.iter())
            .map(|entry| &entry.record)
            .collect()
    }

    /// Remove every revision of one part. Return true when something was
    /// removed.
    pub fn remove(&mut self, part_id: &str) -> bool {
        let removed = self.records.remove(part_id).is_some();
        if removed {
            self.revision += 1;
        }
        removed
    }

    /// Remove one exact revision of one part.
    pub fn remove_version(&mut self, part_id: &str, version: u32) -> Option<PartRecord> {
        let entries = self.records.get_mut(part_id)?;
        let index = entries.iter().position(|entry| entry.version == version)?;
        let removed = entries.remove(index);
        if entries.is_empty() {
            self.records.remove(part_id);
        }
        self.revision += 1;
        Some(removed.record)
    }

    /// Check a stored record against an expected geometry hash.
    ///
    /// A mismatch is [`ParamError::GeometryMismatch`]. A missing part is
    /// [`ParamError::UnknownPart`]. This is the lookup-side geometry check. It
    /// flags a bad record without a write.
    pub fn verify_geometry(
        &self,
        part_id: &str,
        version: u32,
        expected_geometry_ref: &str,
    ) -> Result<(), ParamError> {
        let record = self
            .get_version(part_id, version)
            .ok_or_else(|| ParamError::UnknownPart {
                part_id: part_id.to_string(),
                version,
            })?;
        if record.geometry_ref != expected_geometry_ref {
            return Err(ParamError::GeometryMismatch {
                part_id: part_id.to_string(),
                version,
                expected: expected_geometry_ref.to_string(),
                found: record.geometry_ref.clone(),
            });
        }
        Ok(())
    }

    /// Reject a foreign library schema or an unsupported library version.
    pub fn check_header(&self) -> Result<(), ParamError> {
        if self.schema != LIBRARY_SCHEMA {
            return Err(ParamError::UnsupportedLibrarySchema {
                found: self.schema.clone(),
                expected: LIBRARY_SCHEMA.to_string(),
            });
        }
        if self.version != LIBRARY_VERSION {
            return Err(ParamError::UnsupportedLibraryVersion {
                found: self.version,
                expected: LIBRARY_VERSION,
            });
        }
        Ok(())
    }

    /// Serialize the library to compact JSON.
    pub fn to_json(&self) -> Result<String, ParamError> {
        Ok(serde_json::to_string(self)?)
    }

    /// Serialize the library to indented JSON.
    pub fn to_json_pretty(&self) -> Result<String, ParamError> {
        Ok(serde_json::to_string_pretty(self)?)
    }

    /// Serialize the library to a JSON byte buffer.
    pub fn to_bytes(&self) -> Result<Vec<u8>, ParamError> {
        Ok(serde_json::to_vec(self)?)
    }

    /// Deserialize a library from JSON. Check the header and the contents.
    pub fn from_json(json: &str) -> Result<Self, ParamError> {
        let library: Self = serde_json::from_str(json)?;
        library.validate()?;
        Ok(library)
    }

    /// Deserialize a library from a byte buffer. Check the header and the
    /// contents.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, ParamError> {
        let library: Self = serde_json::from_slice(bytes)?;
        library.validate()?;
        Ok(library)
    }

    /// Write the library as JSON to a file.
    pub fn save_to_file(&self, path: impl AsRef<Path>) -> Result<(), ParamError> {
        std::fs::write(path, self.to_bytes()?)?;
        Ok(())
    }

    /// Read a library as JSON from a file.
    pub fn load_from_file(path: impl AsRef<Path>) -> Result<Self, ParamError> {
        let bytes = std::fs::read(path)?;
        Self::from_bytes(&bytes)
    }

    /// Check the library header, every record header, and every key.
    fn validate(&self) -> Result<(), ParamError> {
        self.check_header()?;
        for (part_id, entries) in &self.records {
            let mut seen = BTreeSet::new();
            for entry in entries {
                if entry.record.part_id != *part_id {
                    return Err(ParamError::LibraryIntegrity {
                        reason: format!(
                            "key {part_id:?} does not match record part id {:?}",
                            entry.record.part_id
                        ),
                    });
                }
                if !seen.insert(entry.version) {
                    return Err(ParamError::LibraryIntegrity {
                        reason: format!("part {part_id:?} holds version {} twice", entry.version),
                    });
                }
                entry.record.check_header()?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::provenance::{Method, Provenance, Validation};
    use crate::quantity::Quantity;

    fn record(part_id: &str, geometry_ref: &str) -> PartRecord {
        let q = Quantity::derived(1.0, "1").unwrap();
        PartRecord {
            schema: crate::record::SCHEMA.to_string(),
            version: crate::record::VERSION,
            part_id: part_id.to_string(),
            geometry_ref: geometry_ref.to_string(),
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
    fn new_library_has_the_header_and_no_records() {
        let library = ParameterLibrary::new();
        assert_eq!(library.schema(), LIBRARY_SCHEMA);
        assert_eq!(library.version(), LIBRARY_VERSION);
        assert_eq!(library.revision(), 0);
        assert!(library.is_empty());
        assert!(library.check_header().is_ok());
    }

    #[test]
    fn duplicate_version_is_rejected_without_a_revision_bump() {
        let mut library = ParameterLibrary::new();
        library.insert(record("a", "g1"), 1).unwrap();
        let err = library.insert(record("a", "g2"), 1).unwrap_err();
        assert!(matches!(
            err,
            ParamError::DuplicateRecord {
                part_id,
                version: 1
            } if part_id == "a"
        ));
        assert_eq!(library.revision(), 1);
    }

    #[test]
    fn geometry_checked_insert_rejects_a_mismatch() {
        let mut library = ParameterLibrary::new();
        let err = library
            .insert_with_geometry(record("a", "sha256:aa"), 1, "sha256:bb")
            .unwrap_err();
        assert!(matches!(err, ParamError::GeometryMismatch { .. }));
        assert!(library.is_empty());
        assert_eq!(library.revision(), 0);

        library
            .insert_with_geometry(record("a", "sha256:aa"), 1, "sha256:aa")
            .unwrap();
        assert_eq!(library.revision(), 1);
    }

    #[test]
    fn lookup_geometry_flags_a_mismatch() {
        let mut library = ParameterLibrary::new();
        library.insert(record("a", "sha256:aa"), 1).unwrap();
        assert!(library.verify_geometry("a", 1, "sha256:aa").is_ok());
        assert!(matches!(
            library.verify_geometry("a", 1, "sha256:zz"),
            Err(ParamError::GeometryMismatch { .. })
        ));
        assert!(matches!(
            library.verify_geometry("missing", 1, "sha256:aa"),
            Err(ParamError::UnknownPart { .. })
        ));
    }

    #[test]
    fn json_round_trip_preserves_the_library() {
        let mut library = ParameterLibrary::new();
        library.insert(record("a", "g1"), 1).unwrap();
        library.insert(record("a", "g2"), 2).unwrap();
        library.insert(record("b", "g3"), 1).unwrap();
        let json = library.to_json().unwrap();
        let back = ParameterLibrary::from_json(&json).unwrap();
        assert_eq!(back, library);
    }

    #[test]
    fn foreign_library_schema_is_rejected() {
        let mut library = ParameterLibrary::new();
        library.schema = "nanocad.param.other".to_string();
        assert!(matches!(
            library.check_header(),
            Err(ParamError::UnsupportedLibrarySchema { .. })
        ));
    }

    #[test]
    fn foreign_library_version_is_rejected() {
        let mut library = ParameterLibrary::new();
        library.version = 99;
        assert!(matches!(
            library.check_header(),
            Err(ParamError::UnsupportedLibraryVersion {
                found: 99,
                expected: 1
            })
        ));
    }

    #[test]
    fn key_mismatch_is_rejected_on_load() {
        let mut library = ParameterLibrary::new();
        library.insert(record("a", "g1"), 1).unwrap();
        let mut value: serde_json::Value =
            serde_json::from_str(&library.to_json().unwrap()).unwrap();
        let mut entry = value["records"]["a"][0].clone();
        entry["record"]["part_id"] = serde_json::Value::String("b".to_string());
        value["records"]["a"] = serde_json::Value::Array(vec![entry]);
        let json = serde_json::to_string(&value).unwrap();
        assert!(matches!(
            ParameterLibrary::from_json(&json),
            Err(ParamError::LibraryIntegrity { .. })
        ));
    }
}
