//! PDB import and export.
//!
//! PDB is a fixed-column text format. This clean-room reader handles the
//! records a molecular model needs:
//!
//! - `ATOM` and `HETATM` give the serial number, the element, and the x, y, z
//!   coordinates in Angstrom. The model stores SI metres, so this module
//!   converts at the boundary.
//! - `CONECT` gives the bonds. A repeated index raises the bond order, so a
//!   partner that appears twice is a double bond.
//! - `CRYST1` gives the unit cell. The reader stores the cell in the document
//!   metadata under `pdb_cryst1_*` keys in SI units.
//! - `TITLE` gives the document name.
//! - `MODEL` and `ENDMDL` start a new part per model. A file without `MODEL`
//!   becomes one part.
//!
//! This is a subset, not full PDB. The reader ignores ANISOU, TER, SEQRES,
//! HELIX, and every display record. The reader keeps no partial charge,
//! because PDB has no charge field that maps to the model. The writer emits
//! element symbols and three decimal places of Angstrom.

use std::collections::{BTreeMap, HashMap};

use nanocad_model::{Atom, Bond, BondType, Document, Element, Part, Topology};

use crate::error::FormatError;
use crate::periodic;

const METRES_PER_ANGSTROM: f64 = 1.0e-10;
const ANGSTROMS_PER_METRE: f64 = 1.0e10;

const CRYST1_LENGTH_KEYS: [&str; 3] = ["pdb_cryst1_a_m", "pdb_cryst1_b_m", "pdb_cryst1_c_m"];
const CRYST1_ANGLE_KEYS: [&str; 3] = [
    "pdb_cryst1_alpha_deg",
    "pdb_cryst1_beta_deg",
    "pdb_cryst1_gamma_deg",
];
const CRYST1_SPACE_GROUP_KEY: &str = "pdb_cryst1_space_group";

/// The parsed `CRYST1` record: three lengths in metres, three angles in
/// degrees, and the space group when it is present.
type Cryst1Record = ([f64; 3], [f64; 3], Option<String>);

/// Imports PDB text into a document.
pub fn import_pdb(text: &str) -> Result<Document, FormatError> {
    let mut parts: Vec<Part> = Vec::new();
    let mut part_serials: Vec<HashMap<u32, u32>> = Vec::new();
    let mut directed_counts: Vec<HashMap<(u32, u32), u8>> = Vec::new();
    let mut current_part: Option<usize> = None;
    let mut metadata: BTreeMap<String, String> = BTreeMap::new();
    let mut document_name: Option<String> = None;
    let mut title_seen = false;

    for (line_index, raw_line) in text.lines().enumerate() {
        let line_number = line_index + 1;
        let record = field(raw_line, 0, 6).trim();

        match record {
            "ATOM" | "HETATM" => {
                let part_index = ensure_part(
                    &mut parts,
                    &mut part_serials,
                    &mut directed_counts,
                    &mut current_part,
                    document_name.as_deref(),
                );
                let (serial, atom) = parse_atom_line(raw_line, line_number)?;
                let atom_index = parts[part_index].topology.add_atom(atom);
                if part_serials[part_index]
                    .insert(serial, atom_index)
                    .is_some()
                {
                    return Err(FormatError::PdbParse {
                        line: line_number,
                        message: format!("atom serial {serial} appears twice"),
                    });
                }
            }
            "CONECT" => {
                let part_index = ensure_part(
                    &mut parts,
                    &mut part_serials,
                    &mut directed_counts,
                    &mut current_part,
                    document_name.as_deref(),
                );
                let (serial, partners) = parse_conect_line(raw_line, line_number)?;
                for partner in partners {
                    if partner == serial {
                        continue;
                    }
                    let entry = directed_counts[part_index]
                        .entry((serial, partner))
                        .or_insert(0);
                    *entry = entry.saturating_add(1).min(Bond::MAX_ORDER);
                }
            }
            "CRYST1" => {
                let (lengths_m, angles_deg, space_group) = parse_cryst1_line(raw_line)?;
                for (index, key) in CRYST1_LENGTH_KEYS.iter().enumerate() {
                    metadata.insert((*key).to_owned(), lengths_m[index].to_string());
                }
                for (index, key) in CRYST1_ANGLE_KEYS.iter().enumerate() {
                    metadata.insert((*key).to_owned(), angles_deg[index].to_string());
                }
                if let Some(space_group) = space_group {
                    metadata.insert(CRYST1_SPACE_GROUP_KEY.to_owned(), space_group);
                }
            }
            "TITLE" => {
                if !title_seen {
                    if let Some(title) = optional_text(raw_line, 10) {
                        document_name = Some(title);
                    }
                    title_seen = true;
                }
            }
            "MODEL" if current_part.is_some_and(|index| parts[index].atom_count() > 0) => {
                current_part = None;
            }
            _ => {}
        }
    }

    if parts.iter().all(|part| part.atom_count() == 0) {
        return Err(FormatError::PdbParse {
            line: 1,
            message: "the file has no ATOM or HETATM records".to_owned(),
        });
    }

    for (part_index, part) in parts.iter_mut().enumerate() {
        let mut orders: BTreeMap<(u32, u32), u8> = BTreeMap::new();
        for (&(start, end), &count) in &directed_counts[part_index] {
            if start == end {
                continue;
            }
            let key = if start < end {
                (start, end)
            } else {
                (end, start)
            };
            let entry = orders.entry(key).or_insert(0);
            *entry = (*entry).max(count);
        }
        for (&(u, v), &order) in &orders {
            let serials = &part_serials[part_index];
            let (u_index, v_index) = match (serials.get(&u), serials.get(&v)) {
                (Some(&u_index), Some(&v_index)) => (u_index, v_index),
                _ => {
                    return Err(FormatError::PdbParse {
                        line: 1,
                        message: format!("CONECT names unknown atom serial {u} or {v}"),
                    })
                }
            };
            let order = order.clamp(1, Bond::MAX_ORDER);
            let bond_type = match order {
                2 => BondType::Double,
                3 => BondType::Triple,
                _ => BondType::Single,
            };
            part.topology
                .add_bond(Bond::new(u_index, v_index, order, bond_type))?;
        }
    }

    let name = document_name.unwrap_or_else(|| "pdb".to_owned());
    let mut document = Document::new(name);
    document.metadata = metadata;
    for part in parts {
        if part.atom_count() > 0 {
            document.add_part(part);
        }
    }
    Ok(document)
}

/// Exports a document as PDB text.
///
/// Every atom of every part is written in order with a fresh serial number.
/// Bonds are written once from each end, as the format expects. Partial
/// charges and atom types beyond four characters do not survive.
pub fn export_pdb(document: &Document) -> Result<String, FormatError> {
    let mut out = String::new();

    if !document.name.is_empty() {
        out.push_str(&format!("TITLE     {}\n", document.name));
    }

    if let Some(cryst1) = format_cryst1(&document.metadata) {
        out.push_str(&cryst1);
    }

    let mut serials_per_part: Vec<Vec<u32>> = Vec::with_capacity(document.parts.len());
    let mut next_serial: u32 = 1;
    for part in &document.parts {
        let mut serials = Vec::with_capacity(part.atom_count());
        for atom in part.topology.atoms() {
            let serial = next_serial;
            next_serial = if next_serial >= 99_999 {
                1
            } else {
                next_serial + 1
            };
            serials.push(serial);
            out.push_str(&format_atom(serial, &atom));
        }
        serials_per_part.push(serials);
    }

    for (part_index, part) in document.parts.iter().enumerate() {
        let serials = &serials_per_part[part_index];
        let mut partners: Vec<BTreeMap<u32, u8>> = vec![BTreeMap::new(); part.atom_count()];
        for bond in part.topology.bonds() {
            let order = bond.order.clamp(1, Bond::MAX_ORDER);
            if let Some(slot) = partners.get_mut(bond.u as usize) {
                let entry = slot.entry(bond.v).or_insert(0);
                *entry = (*entry).max(order);
            }
            if let Some(slot) = partners.get_mut(bond.v as usize) {
                let entry = slot.entry(bond.u).or_insert(0);
                *entry = (*entry).max(order);
            }
        }
        for (atom_index, slots) in partners.iter().enumerate() {
            if slots.is_empty() {
                continue;
            }
            out.push_str(&format!("CONECT{:>5}", serials[atom_index]));
            for (&partner, &order) in slots {
                for _ in 0..order {
                    out.push_str(&format!("{:>5}", serials[partner as usize]));
                }
            }
            out.push('\n');
        }
    }

    out.push_str("END\n");
    Ok(out)
}

fn ensure_part(
    parts: &mut Vec<Part>,
    part_serials: &mut Vec<HashMap<u32, u32>>,
    directed_counts: &mut Vec<HashMap<(u32, u32), u8>>,
    current_part: &mut Option<usize>,
    document_name: Option<&str>,
) -> usize {
    if let Some(index) = *current_part {
        return index;
    }
    let name = document_name.unwrap_or("pdb").to_owned();
    parts.push(Part::new(name, Topology::new()));
    part_serials.push(HashMap::new());
    directed_counts.push(HashMap::new());
    let index = parts.len() - 1;
    *current_part = Some(index);
    index
}

fn parse_atom_line(line: &str, line_number: usize) -> Result<(u32, Atom), FormatError> {
    let serial = parse_u32(field(line, 6, 11), line_number, "atom serial")?;
    let x = parse_f64(field(line, 30, 38), line_number, "x coordinate")?;
    let y = parse_f64(field(line, 38, 46), line_number, "y coordinate")?;
    let z = parse_f64(field(line, 46, 54), line_number, "z coordinate")?;
    let element_field = field(line, 76, 78).trim();
    let atom_name = field(line, 12, 16).trim();
    let symbol = if element_field.is_empty() {
        element_from_name(atom_name)
    } else {
        element_field.to_owned()
    };
    let atomic_number =
        periodic::atomic_number_for(&symbol).ok_or_else(|| FormatError::PdbParse {
            line: line_number,
            message: format!("unknown element {symbol:?}"),
        })?;
    let element = Element::from_atomic_number(atomic_number)
        .ok_or(FormatError::UnknownElement(atomic_number))?;
    let atom_type = if atom_name.is_empty() {
        periodic::symbol_for(atomic_number)
            .unwrap_or("X")
            .to_owned()
    } else {
        atom_name.to_owned()
    };
    Ok((
        serial,
        Atom::new(
            element,
            [
                x * METRES_PER_ANGSTROM,
                y * METRES_PER_ANGSTROM,
                z * METRES_PER_ANGSTROM,
            ],
            0.0,
            atom_type,
        ),
    ))
}

fn parse_conect_line(line: &str, line_number: usize) -> Result<(u32, Vec<u32>), FormatError> {
    let serial = parse_u32(field(line, 6, 11), line_number, "CONECT serial")?;
    let mut partners = Vec::new();
    let mut start = 11;
    while start < line.len().min(31) {
        let text = field(line, start, start + 5).trim();
        if !text.is_empty() {
            partners.push(text.parse().map_err(|_| FormatError::PdbParse {
                line: line_number,
                message: format!("CONECT partner {text:?} is not a number"),
            })?);
        }
        start += 5;
    }
    Ok((serial, partners))
}

fn parse_cryst1_line(line: &str) -> Result<Cryst1Record, FormatError> {
    let mut lengths_m = [0.0_f64; 3];
    for (index, (start, end)) in [(6, 15), (15, 24), (24, 33)].iter().enumerate() {
        lengths_m[index] =
            parse_f64(field(line, *start, *end), 1, "CRYST1 length")? * METRES_PER_ANGSTROM;
    }
    let mut angles_deg = [0.0_f64; 3];
    for (index, (start, end)) in [(33, 40), (40, 47), (47, 54)].iter().enumerate() {
        angles_deg[index] = parse_f64(field(line, *start, *end), 1, "CRYST1 angle")?;
    }
    let space_group = {
        let text = field(line, 55, 66).trim();
        if text.is_empty() {
            None
        } else {
            Some(text.to_owned())
        }
    };
    Ok((lengths_m, angles_deg, space_group))
}

fn format_atom(serial: u32, atom: &Atom) -> String {
    let symbol = periodic::symbol_for(atom.element.atomic_number())
        .map(str::to_owned)
        .unwrap_or_else(|| format!("Z{}", atom.element.atomic_number()));
    let name: String = if atom.atom_type.is_empty() {
        symbol.chars().take(4).collect()
    } else {
        atom.atom_type.chars().take(4).collect()
    };
    let x = atom.position_m[0] * ANGSTROMS_PER_METRE;
    let y = atom.position_m[1] * ANGSTROMS_PER_METRE;
    let z = atom.position_m[2] * ANGSTROMS_PER_METRE;
    format!(
        "ATOM  {serial:>5} {name:<4} MOL A   1    {x:>8.3}{y:>8.3}{z:>8.3}  1.00  0.00          {element:>2}\n",
        element = symbol
    )
}

fn format_cryst1(metadata: &BTreeMap<String, String>) -> Option<String> {
    let mut lengths = [0.0_f64; 3];
    for (index, key) in CRYST1_LENGTH_KEYS.iter().enumerate() {
        let text = metadata.get(*key)?;
        lengths[index] = text.parse::<f64>().ok()? * ANGSTROMS_PER_METRE;
    }
    let mut angles = [90.0_f64; 3];
    for (index, key) in CRYST1_ANGLE_KEYS.iter().enumerate() {
        if let Some(text) = metadata.get(*key) {
            angles[index] = text.parse().ok()?;
        }
    }
    let space_group = metadata
        .get(CRYST1_SPACE_GROUP_KEY)
        .map(String::as_str)
        .unwrap_or("P 1");
    Some(format!(
        "CRYST1{:>9.3}{:>9.3}{:>9.3}{:>7.2}{:>7.2}{:>7.2} {:<11}\n",
        lengths[0], lengths[1], lengths[2], angles[0], angles[1], angles[2], space_group
    ))
}

fn field(line: &str, start: usize, end: usize) -> &str {
    let length = line.len();
    if start >= end || start >= length {
        return "";
    }
    line.get(start..end.min(length)).unwrap_or("")
}

fn element_from_name(atom_name: &str) -> String {
    let letters: String = atom_name
        .chars()
        .take_while(|character| character.is_ascii_alphabetic())
        .collect();
    if letters.is_empty() {
        atom_name.to_owned()
    } else {
        letters
    }
}

fn optional_text(line: &str, start: usize) -> Option<String> {
    let text = line.get(start..).unwrap_or("").trim();
    if text.is_empty() {
        None
    } else {
        Some(text.to_owned())
    }
}

fn parse_u32(text: &str, line: usize, what: &str) -> Result<u32, FormatError> {
    text.trim().parse().map_err(|_| FormatError::PdbParse {
        line,
        message: format!("{what} {:?} is not a number", text.trim()),
    })
}

fn parse_f64(text: &str, line: usize, what: &str) -> Result<f64, FormatError> {
    text.trim().parse().map_err(|_| FormatError::PdbParse {
        line,
        message: format!("{what} {:?} is not a number", text.trim()),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ncz::{from_ncz_bytes, to_ncz_bytes};

    const SAMPLE: &str = "\
TITLE     WATER LIKE
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
ATOM      1  O   MOL A   1       0.000   0.000   0.000  1.00  0.00           O
ATOM      2  H   MOL A   1       1.000   0.000   0.000  1.00  0.00           H
ATOM      3  H   MOL A   1      -0.500   1.000   0.000  1.00  0.00           H
CONECT    1    2    2    3
CONECT    2    1
CONECT    3    1
END
";

    #[test]
    fn a_pdb_text_imports_the_expected_topology() {
        let document = import_pdb(SAMPLE).expect("import");
        assert_eq!(document.name, "WATER LIKE");
        let part = document.part(0).expect("part");
        assert_eq!(part.atom_count(), 3);
        assert_eq!(part.bond_count(), 2);
        assert_eq!(part.topology.element(0), Some(Element::OXYGEN));
        assert_eq!(part.topology.element(1), Some(Element::HYDROGEN));
        assert_eq!(part.topology.bond_order(0), Some(2));
        assert_eq!(part.topology.bond_type(0), Some(BondType::Double));
    }

    #[test]
    fn coordinates_convert_from_angstrom_to_metres() {
        let document = import_pdb(SAMPLE).expect("import");
        let position = document
            .part(0)
            .and_then(|part| part.topology.position_m(1))
            .expect("position");
        assert!((position[0] - 1.0e-10).abs() < 1.0e-20);
    }

    #[test]
    fn the_cryst1_box_is_kept_in_si_metres() {
        let document = import_pdb(SAMPLE).expect("import");
        let length = document
            .metadata
            .get("pdb_cryst1_a_m")
            .expect("length")
            .parse::<f64>()
            .expect("number");
        assert!((length - 20.0e-10).abs() < 1.0e-20);
        assert_eq!(
            document.metadata.get("pdb_cryst1_space_group"),
            Some(&"P 1".to_owned())
        );
    }

    #[test]
    fn a_pdb_text_round_trips_through_export() {
        let document = import_pdb(SAMPLE).expect("import");
        let text = export_pdb(&document).expect("export");
        let again = import_pdb(&text).expect("reimport");
        assert_eq!(again.name, document.name);
        assert_eq!(again.part_count(), document.part_count());
        assert_eq!(again.metadata, document.metadata);
        let first = document.part(0).expect("part");
        let second = again.part(0).expect("part");
        assert_eq!(first.bond_count(), second.bond_count());
        assert_eq!(first.bond_count(), 2);
        for index in 0..first.atom_count() {
            for axis in 0..3 {
                let left = first.topology.position_m(index).expect("left")[axis];
                let right = second.topology.position_m(index).expect("right")[axis];
                assert!((left - right).abs() < 1.0e-20);
            }
        }
    }

    #[test]
    fn a_pdb_text_round_trips_through_ncz() {
        let document = import_pdb(SAMPLE).expect("import");
        let bytes = to_ncz_bytes(&document).expect("write ncz");
        let reloaded = from_ncz_bytes(&bytes).expect("read ncz");
        assert_eq!(reloaded, document);
        assert_eq!(
            export_pdb(&reloaded).expect("export"),
            export_pdb(&document).expect("export")
        );
    }

    #[test]
    fn a_file_with_no_atoms_is_rejected() {
        assert!(matches!(import_pdb(""), Err(FormatError::PdbParse { .. })));
    }

    #[test]
    fn a_malformed_coordinate_is_rejected_not_panicking() {
        let text =
            "ATOM      1  O   MOL A   1       oops   0.000   0.000  1.00  0.00           O\n";
        assert!(matches!(
            import_pdb(text),
            Err(FormatError::PdbParse { .. })
        ));
    }

    #[test]
    fn a_single_part_is_written_without_an_empty_second_part() {
        let document = import_pdb(SAMPLE).expect("import");
        let text = export_pdb(&document).expect("export");
        assert_eq!(text.matches("ATOM  ").count(), 3);
        assert_eq!(text.matches("CONECT").count(), 3);
    }
}
