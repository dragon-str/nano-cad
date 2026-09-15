//! NCZ: the native, versioned document format.
//!
//! An NCZ file is a ZIP archive. A header entry carries the format tag and the
//! schema version. A document entry carries the JSON metadata. Coordinates,
//! charges, elements, and bonds live in per-part binary entries as little-endian
//! arrays. Coordinates are `f64` so a round trip is bit for bit.
//!
//! Archive layout:
//!
//! ```text
//! header.json            { format, schema_version, generator }
//! document.json          { name, metadata, parts: [ { name, material,
//!                          metadata, atom_count, bond_count, *_file } ] }
//! part_<i>_types.json    JSON array of atom type strings
//! part_<i>_elements.bin  u8[atom_count]
//! part_<i>_positions.bin f64 little-endian[3 * atom_count]
//! part_<i>_charges.bin   f64 little-endian[atom_count]
//! part_<i>_bonds.bin     per bond: u32 u, u32 v, u8 order, u8 bond_type
//! ```

use std::collections::BTreeMap;
use std::io::{Cursor, Read, Seek, Write};

use nanocad_model::{Atom, Bond, BondType, Document, Element, Part, Topology};
use serde::{Deserialize, Serialize};
use zip::write::SimpleFileOptions;
use zip::{CompressionMethod, ZipArchive, ZipWriter};

use crate::error::FormatError;

/// The format tag in every NCZ header.
pub const NCZ_FORMAT: &str = "nanocad.ncz";
/// The NCZ schema version written by this crate.
pub const NCZ_SCHEMA_VERSION: u32 = 1;

const HEADER_ENTRY: &str = "header.json";
const DOCUMENT_ENTRY: &str = "document.json";

#[derive(Serialize, Deserialize)]
struct NczHeader {
    format: String,
    schema_version: u32,
    generator: String,
}

#[derive(Serialize, Deserialize)]
struct NczDocument {
    name: String,
    metadata: BTreeMap<String, String>,
    parts: Vec<NczPart>,
}

#[derive(Serialize, Deserialize)]
struct NczPart {
    name: String,
    material: String,
    metadata: BTreeMap<String, String>,
    atom_count: u64,
    bond_count: u64,
    types_file: String,
    elements_file: String,
    positions_file: String,
    charges_file: String,
    bonds_file: String,
}

/// Writes a document to an NCZ archive and returns the inner writer.
pub fn write_ncz<W: Write + Seek>(writer: W, document: &Document) -> Result<W, FormatError> {
    let mut zip = ZipWriter::new(writer);
    let options = SimpleFileOptions::default().compression_method(CompressionMethod::Stored);

    let header = NczHeader {
        format: NCZ_FORMAT.to_owned(),
        schema_version: NCZ_SCHEMA_VERSION,
        generator: format!("nanocad-format {}", env!("CARGO_PKG_VERSION")),
    };
    write_json_entry(&mut zip, HEADER_ENTRY, &header, options)?;

    let mut parts = Vec::with_capacity(document.parts.len());
    for (index, part) in document.parts.iter().enumerate() {
        let types_file = format!("part_{index}_types.json");
        let elements_file = format!("part_{index}_elements.bin");
        let positions_file = format!("part_{index}_positions.bin");
        let charges_file = format!("part_{index}_charges.bin");
        let bonds_file = format!("part_{index}_bonds.bin");

        let atom_types: Vec<&str> = part
            .topology
            .atom_types()
            .iter()
            .map(String::as_str)
            .collect();
        write_json_entry(&mut zip, &types_file, &atom_types, options)?;

        let elements: Vec<u8> = part
            .topology
            .atoms()
            .map(|atom| atom.element.atomic_number())
            .collect();
        write_bytes_entry(&mut zip, &elements_file, &elements, options)?;

        let mut positions = Vec::with_capacity(part.topology.positions_m().len() * 8);
        for value in part.topology.positions_m() {
            positions.extend_from_slice(&value.to_le_bytes());
        }
        write_bytes_entry(&mut zip, &positions_file, &positions, options)?;

        let mut charges = Vec::with_capacity(part.topology.charges_c().len() * 8);
        for value in part.topology.charges_c() {
            charges.extend_from_slice(&value.to_le_bytes());
        }
        write_bytes_entry(&mut zip, &charges_file, &charges, options)?;

        let mut bonds = Vec::with_capacity(part.topology.bond_count() * 10);
        for bond in part.topology.bonds() {
            bonds.extend_from_slice(&bond.u.to_le_bytes());
            bonds.extend_from_slice(&bond.v.to_le_bytes());
            bonds.push(bond.order);
            bonds.push(bond.bond_type.to_u8());
        }
        write_bytes_entry(&mut zip, &bonds_file, &bonds, options)?;

        parts.push(NczPart {
            name: part.name.clone(),
            material: part.material.clone(),
            metadata: part.metadata.clone(),
            atom_count: part.atom_count() as u64,
            bond_count: part.bond_count() as u64,
            types_file,
            elements_file,
            positions_file,
            charges_file,
            bonds_file,
        });
    }

    let document_entry = NczDocument {
        name: document.name.clone(),
        metadata: document.metadata.clone(),
        parts,
    };
    write_json_entry(&mut zip, DOCUMENT_ENTRY, &document_entry, options)?;

    Ok(zip.finish()?)
}

/// Reads a document from an NCZ archive.
pub fn read_ncz<R: Read + Seek>(reader: R) -> Result<Document, FormatError> {
    let mut archive = ZipArchive::new(reader)?;
    read_archive(&mut archive)
}

/// Writes a document to an in-memory NCZ archive.
pub fn to_ncz_bytes(document: &Document) -> Result<Vec<u8>, FormatError> {
    let cursor = write_ncz(Cursor::new(Vec::new()), document)?;
    Ok(cursor.into_inner())
}

/// Reads a document from an in-memory NCZ archive.
pub fn from_ncz_bytes(bytes: &[u8]) -> Result<Document, FormatError> {
    read_ncz(Cursor::new(bytes))
}

fn read_archive<R: Read + Seek>(archive: &mut ZipArchive<R>) -> Result<Document, FormatError> {
    let header: NczHeader = read_json_entry(archive, HEADER_ENTRY)?;
    if header.format != NCZ_FORMAT {
        return Err(FormatError::UnsupportedFormat {
            found: header.format,
            expected: NCZ_FORMAT.to_owned(),
        });
    }
    if header.schema_version != NCZ_SCHEMA_VERSION {
        return Err(FormatError::UnsupportedSchemaVersion {
            found: header.schema_version,
            expected: NCZ_SCHEMA_VERSION,
        });
    }

    let document: NczDocument = read_json_entry(archive, DOCUMENT_ENTRY)?;
    let mut parts = Vec::with_capacity(document.parts.len());
    for entry in &document.parts {
        parts.push(read_part(archive, entry)?);
    }

    let mut restored = Document::new(document.name);
    restored.metadata = document.metadata;
    for part in parts {
        restored.add_part(part);
    }
    Ok(restored)
}

fn read_part<R: Read + Seek>(
    archive: &mut ZipArchive<R>,
    entry: &NczPart,
) -> Result<Part, FormatError> {
    let atom_count = usize::try_from(entry.atom_count).map_err(|_| FormatError::IntegerOverflow)?;
    let bond_count = usize::try_from(entry.bond_count).map_err(|_| FormatError::IntegerOverflow)?;

    let atom_types: Vec<String> = read_json_entry(archive, &entry.types_file)?;
    let elements = read_bytes_entry(archive, &entry.elements_file)?;
    let positions = read_f64_entry(archive, &entry.positions_file, atom_count * 3)?;
    let charges = read_f64_entry(archive, &entry.charges_file, atom_count)?;
    let bonds = read_bytes_entry(archive, &entry.bonds_file)?;

    check_len(&entry.types_file, atom_types.len(), atom_count)?;
    check_len(&entry.elements_file, elements.len(), atom_count)?;
    check_len(&entry.bonds_file, bonds.len(), bond_count * 10)?;

    let mut topology = Topology::new();
    for index in 0..atom_count {
        let element = Element::from_atomic_number(elements[index])
            .ok_or(FormatError::UnknownElement(elements[index]))?;
        let position_m = [
            positions[index * 3],
            positions[index * 3 + 1],
            positions[index * 3 + 2],
        ];
        topology.add_atom(Atom::new(
            element,
            position_m,
            charges[index],
            atom_types[index].clone(),
        ));
    }

    for index in 0..bond_count {
        let base = index * 10;
        let u = u32::from_le_bytes([
            bonds[base],
            bonds[base + 1],
            bonds[base + 2],
            bonds[base + 3],
        ]);
        let v = u32::from_le_bytes([
            bonds[base + 4],
            bonds[base + 5],
            bonds[base + 6],
            bonds[base + 7],
        ]);
        let order = bonds[base + 8];
        let tag = bonds[base + 9];
        let bond_type =
            BondType::from_u8(tag).ok_or(nanocad_model::ModelError::InvalidBondType(tag))?;
        topology.add_bond(Bond::new(u, v, order, bond_type))?;
    }

    let mut part = Part::new(entry.name.clone(), topology).with_material(entry.material.clone());
    part.metadata = entry.metadata.clone();
    Ok(part)
}

fn check_len(name: &str, found: usize, expected: usize) -> Result<(), FormatError> {
    if found == expected {
        Ok(())
    } else {
        Err(FormatError::BadArrayLength {
            name: name.to_owned(),
            found,
            expected,
        })
    }
}

fn write_json_entry<W: Write + Seek, T: Serialize>(
    zip: &mut ZipWriter<W>,
    name: &str,
    value: &T,
    options: SimpleFileOptions,
) -> Result<(), FormatError> {
    let bytes = serde_json::to_vec_pretty(value)?;
    zip.start_file(name, options)?;
    zip.write_all(&bytes)?;
    Ok(())
}

fn write_bytes_entry<W: Write + Seek>(
    zip: &mut ZipWriter<W>,
    name: &str,
    bytes: &[u8],
    options: SimpleFileOptions,
) -> Result<(), FormatError> {
    zip.start_file(name, options)?;
    zip.write_all(bytes)?;
    Ok(())
}

fn read_entry<R: Read + Seek>(
    archive: &mut ZipArchive<R>,
    name: &str,
) -> Result<Vec<u8>, FormatError> {
    let mut file = match archive.by_name(name) {
        Ok(file) => file,
        Err(zip::result::ZipError::FileNotFound) => {
            return Err(FormatError::MissingEntry(name.to_owned()));
        }
        Err(error) => return Err(error.into()),
    };
    let mut bytes = Vec::new();
    file.read_to_end(&mut bytes)?;
    Ok(bytes)
}

fn read_json_entry<R: Read + Seek, T: for<'de> Deserialize<'de>>(
    archive: &mut ZipArchive<R>,
    name: &str,
) -> Result<T, FormatError> {
    let bytes = read_entry(archive, name)?;
    let text = String::from_utf8(bytes).map_err(|_| FormatError::InvalidUtf8)?;
    Ok(serde_json::from_str(&text)?)
}

fn read_bytes_entry<R: Read + Seek>(
    archive: &mut ZipArchive<R>,
    name: &str,
) -> Result<Vec<u8>, FormatError> {
    read_entry(archive, name)
}

fn read_f64_entry<R: Read + Seek>(
    archive: &mut ZipArchive<R>,
    name: &str,
    expected: usize,
) -> Result<Vec<f64>, FormatError> {
    let bytes = read_entry(archive, name)?;
    if bytes.len() % 8 != 0 {
        return Err(FormatError::BadArrayLength {
            name: name.to_owned(),
            found: bytes.len(),
            expected: expected * 8,
        });
    }
    let values: Vec<f64> = bytes
        .chunks(8)
        .map(|chunk| {
            f64::from_le_bytes([
                chunk[0], chunk[1], chunk[2], chunk[3], chunk[4], chunk[5], chunk[6], chunk[7],
            ])
        })
        .collect();
    check_len(name, values.len(), expected)?;
    Ok(values)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nanocad_model::{BondType, Element};

    fn sample_document() -> Document {
        let mut topology = Topology::new();
        topology.add_atom(Atom::new(
            Element::CARBON,
            [1.0e-10, 2.0e-10, 3.0e-10],
            -0.2,
            "C3",
        ));
        topology.add_atom(Atom::new(
            Element::OXYGEN,
            [-1.0e-10, 0.5e-10, 0.0],
            -0.4,
            "O2",
        ));
        topology
            .add_bond(Bond::new(0, 1, 2, BondType::Double))
            .expect("valid bond");
        let mut part = Part::new("methanol", topology).with_material("diamond");
        part.metadata
            .insert("generator".to_owned(), "manual".to_owned());
        let mut document = Document::new("water");
        document
            .metadata
            .insert("author".to_owned(), "nanocad".to_owned());
        document.add_part(part);
        document
    }

    #[test]
    fn a_document_round_trips_through_ncz_bytes() {
        let document = sample_document();
        let bytes = to_ncz_bytes(&document).expect("write");
        let restored = from_ncz_bytes(&bytes).expect("read");
        assert_eq!(restored, document);
    }

    #[test]
    fn coordinates_survive_bit_for_bit() {
        let document = sample_document();
        let bytes = to_ncz_bytes(&document).expect("write");
        let restored = from_ncz_bytes(&bytes).expect("read");
        let original = document.part(0).expect("part");
        let result = restored.part(0).expect("part");
        for index in 0..original.atom_count() {
            assert_eq!(
                result.topology.position_m(index),
                original.topology.position_m(index)
            );
        }
    }

    #[test]
    fn a_foreign_schema_version_is_rejected() {
        let bytes = archive_with_header(
            b"{\"format\":\"nanocad.ncz\",\"schema_version\":999,\"generator\":\"x\"}",
        );
        assert!(matches!(
            from_ncz_bytes(&bytes),
            Err(FormatError::UnsupportedSchemaVersion { found: 999, .. })
        ));
    }

    #[test]
    fn a_missing_entry_is_rejected() {
        let bytes = archive_with_header(
            b"{\"format\":\"nanocad.ncz\",\"schema_version\":1,\"generator\":\"x\"}",
        );
        assert!(matches!(
            from_ncz_bytes(&bytes),
            Err(FormatError::MissingEntry(_))
        ));
    }

    fn archive_with_header(header: &[u8]) -> Vec<u8> {
        let mut cursor = Cursor::new(Vec::new());
        {
            let mut writer = ZipWriter::new(&mut cursor);
            let options =
                SimpleFileOptions::default().compression_method(CompressionMethod::Stored);
            writer.start_file(HEADER_ENTRY, options).expect("start");
            writer.write_all(header).expect("write");
            writer.finish().expect("finish");
        }
        cursor.into_inner()
    }
}
