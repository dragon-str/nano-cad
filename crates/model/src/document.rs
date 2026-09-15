use std::collections::BTreeMap;

use crate::encoding::{Reader, Writer};
use crate::error::ModelError;
use crate::part::Part;

/// The document schema version written into every buffer header.
pub const SCHEMA_VERSION: u32 = 1;

const MAGIC: [u8; 4] = *b"NCAD";

/// A versioned, serializable design graph.
///
/// The document holds parts and metadata only. It is separate from simulation
/// state, so it can be saved, loaded, and diffed on its own.
#[derive(Clone, Debug, PartialEq)]
pub struct Document {
    pub schema_version: u32,
    pub name: String,
    pub parts: Vec<Part>,
    pub metadata: BTreeMap<String, String>,
}

impl Document {
    /// Creates an empty document at the current schema version.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            schema_version: SCHEMA_VERSION,
            name: name.into(),
            parts: Vec::new(),
            metadata: BTreeMap::new(),
        }
    }

    /// Appends a part and returns its index.
    pub fn add_part(&mut self, part: Part) -> usize {
        let index = self.parts.len();
        self.parts.push(part);
        index
    }

    /// Returns the number of parts.
    pub fn part_count(&self) -> usize {
        self.parts.len()
    }

    /// Returns a part by index.
    pub fn part(&self, index: usize) -> Option<&Part> {
        self.parts.get(index)
    }

    /// Returns a mutable part by index.
    pub fn part_mut(&mut self, index: usize) -> Option<&mut Part> {
        self.parts.get_mut(index)
    }

    /// Returns the first part with a given name.
    pub fn part_by_name(&self, name: &str) -> Option<&Part> {
        self.parts.iter().find(|part| part.name == name)
    }

    /// Encodes the document. The header carries the magic bytes and the schema
    /// version. All integers and floats are little-endian.
    pub fn to_bytes(&self) -> Result<Vec<u8>, ModelError> {
        let mut writer = Writer::new();
        writer.write_bytes(&MAGIC);
        writer.write_u32(self.schema_version);
        writer.write_str(&self.name)?;
        writer.write_map(&self.metadata)?;
        writer.write_count(self.parts.len())?;
        for part in &self.parts {
            part.encode(&mut writer)?;
        }
        Ok(writer.into_bytes())
    }

    /// Decodes a document. It rejects a bad header, a foreign schema version,
    /// a truncated buffer, and trailing bytes.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, ModelError> {
        let mut reader = Reader::new(bytes);
        let magic = reader.read_array::<4>()?;
        if magic != MAGIC {
            return Err(ModelError::BadMagic(magic));
        }
        let schema_version = reader.read_u32()?;
        if schema_version != SCHEMA_VERSION {
            return Err(ModelError::UnsupportedSchemaVersion {
                found: schema_version,
                expected: SCHEMA_VERSION,
            });
        }
        let name = reader.read_str()?;
        let metadata = reader.read_map()?;
        let part_count = reader.read_count()?;
        let mut parts = Vec::new();
        for _ in 0..part_count {
            parts.push(Part::decode(&mut reader)?);
        }
        if reader.remaining() != 0 {
            return Err(ModelError::TrailingBytes {
                remaining: reader.remaining(),
            });
        }
        Ok(Self {
            schema_version,
            name,
            parts,
            metadata,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::bond::{Bond, BondType};
    use crate::element::Element;
    use crate::topology::Topology;

    fn sample_part(name: &str, atom_type: &str) -> Part {
        let mut topology = Topology::new();
        topology.add_atom(Atom::new(Element::SILICON, [0.0, 0.0, 0.0], 0.5, atom_type));
        topology.add_atom(Atom::new(
            Element::SILICON,
            [2.35e-10, 0.0, 0.0],
            -0.5,
            atom_type,
        ));
        topology
            .add_bond(Bond::new(0, 1, 1, BondType::Single))
            .expect("valid bond");
        Part::new(name, topology).with_material("diamond")
    }

    fn sample_document() -> Document {
        let mut document = Document::new("planetary-gearbox");
        document
            .metadata
            .insert("author".to_owned(), "nanocad".to_owned());
        document
            .metadata
            .insert("revision".to_owned(), "1".to_owned());
        document.add_part(sample_part("sun", "Si"));
        document.add_part(sample_part("ring", "Si"));
        document
    }

    #[test]
    fn the_schema_version_constant_is_one() {
        assert_eq!(SCHEMA_VERSION, 1);
    }

    #[test]
    fn a_new_document_carries_the_current_schema_version() {
        let document = Document::new("empty");
        assert_eq!(document.schema_version, SCHEMA_VERSION);
        assert_eq!(document.part_count(), 0);
    }

    #[test]
    fn a_document_round_trips_through_a_byte_buffer() {
        let document = sample_document();
        let bytes = document.to_bytes().expect("encode");
        let restored = Document::from_bytes(&bytes).expect("decode");
        assert_eq!(restored, document);
        assert_eq!(restored.part_count(), 2);
        assert_eq!(
            restored.part_by_name("sun").map(|part| part.bond_count()),
            Some(1)
        );
    }

    #[test]
    fn every_field_survives_the_round_trip() {
        let document = sample_document();
        let restored = Document::from_bytes(&document.to_bytes().expect("encode")).expect("decode");
        assert_eq!(restored.schema_version, document.schema_version);
        assert_eq!(restored.name, document.name);
        assert_eq!(restored.metadata, document.metadata);
        assert_eq!(restored.parts, document.parts);
        let part = restored.part(0).expect("part");
        assert_eq!(part.name, "sun");
        assert_eq!(part.material, "diamond");
        assert_eq!(part.topology.atom_count(), 2);
        assert_eq!(
            part.topology.atom(1).map(|atom| atom.position_m),
            Some([2.35e-10, 0.0, 0.0])
        );
    }

    #[test]
    fn a_document_with_no_parts_round_trips() {
        let document = Document::new("empty");
        let restored = Document::from_bytes(&document.to_bytes().expect("encode")).expect("decode");
        assert_eq!(restored, document);
    }

    #[test]
    fn a_bad_magic_is_rejected() {
        let mut bytes = sample_document().to_bytes().expect("encode");
        bytes[0] = b'X';
        assert_eq!(
            Document::from_bytes(&bytes),
            Err(ModelError::BadMagic(*b"XCAD"))
        );
    }

    #[test]
    fn a_foreign_schema_version_is_rejected() {
        let mut bytes = sample_document().to_bytes().expect("encode");
        bytes[4..8].copy_from_slice(&999u32.to_le_bytes());
        assert_eq!(
            Document::from_bytes(&bytes),
            Err(ModelError::UnsupportedSchemaVersion {
                found: 999,
                expected: SCHEMA_VERSION,
            })
        );
    }

    #[test]
    fn trailing_bytes_are_rejected() {
        let mut bytes = sample_document().to_bytes().expect("encode");
        bytes.extend_from_slice(&[0, 0]);
        assert_eq!(
            Document::from_bytes(&bytes),
            Err(ModelError::TrailingBytes { remaining: 2 })
        );
    }

    #[test]
    fn every_truncated_prefix_is_an_error_not_a_panic() {
        let bytes = sample_document().to_bytes().expect("encode");
        for len in 0..bytes.len() {
            assert!(
                Document::from_bytes(&bytes[..len]).is_err(),
                "prefix of length {len} decoded"
            );
        }
    }

    #[test]
    fn parts_can_be_reached_mutably_and_by_name() {
        let mut document = sample_document();
        if let Some(part) = document.part_mut(0) {
            part.material = "silicon-carbide".to_owned();
        }
        assert_eq!(
            document
                .part_by_name("sun")
                .map(|part| part.material.as_str()),
            Some("silicon-carbide")
        );
        assert_eq!(document.part(9), None);
        assert_eq!(document.part_by_name("missing"), None);
    }
}
