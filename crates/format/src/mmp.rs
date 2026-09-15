//! MMP import and export for NanoEngineer interoperability.
//!
//! The format is line based. An `atom` line gives a serial number, an atomic
//! number in parentheses, three integer coordinates, and a type token. A
//! `bond<N>` line that follows lists the serials bonded to the current atom.
//! The integer coordinate unit is one thousandth of an Angstrom, so one unit is
//! 1e-13 metres. An `info atom atomtype = <name>` line overrides the type.
//!
//! This is a clean-room reader. It handles the records that a molecular part
//! uses: `atom`, `bond1`, `bond2`, `bond3`, `bonda`, `bondg`, `bondc`,
//! `info atom`, `mol`, and `group`. It ignores display and view records.

use std::collections::{BTreeMap, HashMap};

use nanocad_model::{Atom, Bond, BondType, Document, Element, Part, Topology};

use crate::error::FormatError;

const METRES_PER_MMP_UNIT: f64 = 1.0e-13;

/// Imports MMP text into a document with one part for each `mol` record.
pub fn import_mmp(text: &str) -> Result<Document, FormatError> {
    let mut parts: Vec<Part> = Vec::new();
    let mut current_part: Option<usize> = None;
    let mut current_atom: Option<(usize, u32)> = None;
    let mut serial_to_atom: HashMap<u32, (usize, u32)> = HashMap::new();
    let mut document_name: Option<String> = None;

    for (line_index, raw_line) in text.lines().enumerate() {
        let line_number = line_index + 1;
        let line = raw_line.trim();
        if line.is_empty() {
            continue;
        }
        let keyword = line.split_whitespace().next().unwrap_or("");

        if keyword == "mol" {
            let name = parenthesized_name(line).unwrap_or_else(|| "mmp".to_owned());
            parts.push(Part::new(name.clone(), Topology::new()));
            current_part = Some(parts.len() - 1);
            current_atom = None;
            if document_name.is_none() {
                document_name = Some(name);
            }
        } else if keyword == "atom" {
            let part_index = match current_part {
                Some(index) => index,
                None => {
                    let name = document_name.clone().unwrap_or_else(|| "mmp".to_owned());
                    parts.push(Part::new(name, Topology::new()));
                    current_part = Some(parts.len() - 1);
                    parts.len() - 1
                }
            };
            let (serial, atom) = parse_atom_line(line, line_number)?;
            let atom_index = parts[part_index].topology.add_atom(atom);
            serial_to_atom.insert(serial, (part_index, atom_index));
            current_atom = Some((part_index, atom_index));
        } else if let Some(tag) = keyword.strip_prefix("bond") {
            let (part_index, atom_index) = current_atom.ok_or_else(|| FormatError::MmpParse {
                line: line_number,
                message: "a bond line appears before any atom".to_owned(),
            })?;
            let (order, bond_type) = bond_kind(tag).ok_or_else(|| FormatError::MmpParse {
                line: line_number,
                message: format!("unknown bond record {keyword:?}"),
            })?;
            for partner in line.split_whitespace().skip(1) {
                let partner_serial: u32 = partner.parse().map_err(|_| FormatError::MmpParse {
                    line: line_number,
                    message: format!("bond partner {partner:?} is not a number"),
                })?;
                let (partner_part, partner_atom) = serial_to_atom
                    .get(&partner_serial)
                    .copied()
                    .ok_or_else(|| FormatError::MmpParse {
                        line: line_number,
                        message: format!("bond names unknown atom {partner_serial}"),
                    })?;
                if partner_part != part_index {
                    return Err(FormatError::MmpParse {
                        line: line_number,
                        message: "a bond crosses two mol records".to_owned(),
                    });
                }
                parts[part_index].topology.add_bond(Bond::new(
                    atom_index,
                    partner_atom,
                    order,
                    bond_type,
                ))?;
            }
        } else if keyword == "info" && line.contains("atomtype") {
            if let Some(value) = line.split('=').nth(1) {
                let value = value.trim();
                if !value.is_empty() {
                    if let Some((part_index, atom_index)) = current_atom {
                        parts[part_index]
                            .topology
                            .set_atom_type(atom_index as usize, value)?;
                    }
                }
            }
        }
    }

    if parts.iter().all(|part| part.atom_count() == 0) {
        return Err(FormatError::MmpEmpty);
    }

    let name = document_name.unwrap_or_else(|| "mmp".to_owned());
    let mut document = Document::new(name);
    for part in parts {
        document.add_part(part);
    }
    Ok(document)
}

/// Exports a document as MMP text. Records use the classic integer coordinates.
pub fn export_mmp(document: &Document) -> Result<String, FormatError> {
    let mut out = String::new();
    out.push_str("mmpformat 050920 required; 051103 preferred\n");
    let mut serial: u32 = 0;

    for part in &document.parts {
        out.push_str(&format!("mol ({}) def\n", part.name));
        let atom_count = part.atom_count();
        let mut local_to_serial: Vec<u32> = Vec::with_capacity(atom_count);
        for _ in 0..atom_count {
            serial = serial.checked_add(1).ok_or(FormatError::IntegerOverflow)?;
            local_to_serial.push(serial);
        }

        let mut write_after: Vec<Vec<(u32, u8, BondType)>> = vec![Vec::new(); atom_count];
        for bond in part.topology.bonds() {
            let after = bond.u.max(bond.v) as usize;
            let other = bond.u.min(bond.v);
            if after < atom_count {
                write_after[after].push((other, bond.order, bond.bond_type));
            }
        }

        for index in 0..atom_count {
            let atom = part
                .topology
                .atom(index)
                .ok_or(FormatError::IntegerOverflow)?;
            let x = round_coordinate(atom.position_m[0]);
            let y = round_coordinate(atom.position_m[1]);
            let z = round_coordinate(atom.position_m[2]);
            out.push_str(&format!(
                "atom {} ({}) ({}, {}, {}) def\n",
                local_to_serial[index],
                atom.element.atomic_number(),
                x,
                y,
                z
            ));
            let atom_type = atom.atom_type.as_str();
            if !atom_type.is_empty() && atom_type != "def" {
                out.push_str(&format!("info atom atomtype = {atom_type}\n"));
            }

            let mut groups: BTreeMap<&'static str, Vec<u32>> = BTreeMap::new();
            for &(other, order, bond_type) in &write_after[index] {
                let tag = bond_tag(bond_type, order);
                groups
                    .entry(tag)
                    .or_default()
                    .push(local_to_serial[other as usize]);
            }
            for (tag, mut partners) in groups {
                partners.sort_unstable();
                let joined = partners
                    .iter()
                    .map(u32::to_string)
                    .collect::<Vec<_>>()
                    .join(" ");
                out.push_str(&format!("bond{tag} {joined}\n"));
            }
        }
        out.push_str(&format!("egroup ({})\n", part.name));
    }

    out.push_str(&format!("end molecular machine part {}\n", document.name));
    Ok(out)
}

fn round_coordinate(position_m: f64) -> i64 {
    (position_m / METRES_PER_MMP_UNIT).round() as i64
}

fn parenthesized_name(line: &str) -> Option<String> {
    let start = line.find('(')?;
    let end = line[start + 1..].find(')')? + start + 1;
    Some(line[start + 1..end].trim().to_owned())
}

fn trim_wrappers(token: &str) -> &str {
    token.trim_matches(|character: char| {
        character == '(' || character == ')' || character == ',' || character.is_whitespace()
    })
}

fn parse_atom_line(line: &str, line_number: usize) -> Result<(u32, Atom), FormatError> {
    let mut tokens = line.split_whitespace();
    tokens.next();
    let serial: u32 = tokens
        .next()
        .and_then(|token| token.parse().ok())
        .ok_or_else(|| FormatError::MmpParse {
            line: line_number,
            message: "atom serial is not a number".to_owned(),
        })?;
    let atomic_number: u8 = tokens
        .next()
        .and_then(|token| trim_wrappers(token).parse().ok())
        .ok_or_else(|| FormatError::MmpParse {
            line: line_number,
            message: "atom element is not a number".to_owned(),
        })?;
    let element = Element::from_atomic_number(atomic_number)
        .ok_or(FormatError::UnknownElement(atomic_number))?;
    let mut position_m = [0.0_f64; 3];
    for (axis, slot) in position_m.iter_mut().enumerate() {
        let token = tokens.next().ok_or_else(|| FormatError::MmpParse {
            line: line_number,
            message: format!("atom is missing coordinate {axis}"),
        })?;
        let value: i64 = trim_wrappers(token)
            .parse()
            .map_err(|_| FormatError::MmpParse {
                line: line_number,
                message: format!("atom coordinate {axis} is not a number"),
            })?;
        *slot = value as f64 * METRES_PER_MMP_UNIT;
    }
    let type_token = tokens.next().unwrap_or("def");
    let atom_type = if type_token.is_empty() {
        "def".to_owned()
    } else {
        type_token.to_owned()
    };
    Ok((serial, Atom::new(element, position_m, 0.0, atom_type)))
}

fn bond_kind(tag: &str) -> Option<(u8, BondType)> {
    match tag {
        "1" => Some((1, BondType::Single)),
        "2" => Some((2, BondType::Double)),
        "3" => Some((3, BondType::Triple)),
        "a" => Some((1, BondType::Aromatic)),
        "g" => Some((0, BondType::Unknown)),
        "c" => Some((0, BondType::Unknown)),
        _ => None,
    }
}

fn bond_tag(bond_type: BondType, order: u8) -> &'static str {
    match bond_type {
        BondType::Single => "1",
        BondType::Double => "2",
        BondType::Triple => "3",
        BondType::Aromatic => "a",
        BondType::Unknown => match order {
            2 => "2",
            3 => "3",
            _ => "1",
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "\
mmpformat 050920 required; 051103 preferred
kelvin 300
group (Part)
info opengroup open = True
mol (Synthetic) def
atom 1 (6) (0, 0, 0) def
atom 2 (8) (1200, 0, 0) def
info atom atomtype = sp2
bond2 1
atom 3 (6) (2400, 1000, 0) def
bond1 2
atom 4 (1) (3600, 1000, 0) def
bond1 3
egroup (Part)
end molecular machine part Synthetic
";

    #[test]
    fn an_mmp_text_imports_the_expected_topology() {
        let document = import_mmp(SAMPLE).expect("import");
        assert_eq!(document.name, "Synthetic");
        let part = document.part_by_name("Synthetic").expect("part");
        assert_eq!(part.atom_count(), 4);
        assert_eq!(part.bond_count(), 3);
        assert_eq!(part.topology.element(0), Some(Element::CARBON));
        assert_eq!(part.topology.element(1), Some(Element::OXYGEN));
        assert_eq!(part.topology.atom_type(1), Some("sp2"));
        assert_eq!(part.topology.atom_type(0), Some("def"));
        assert_eq!(part.topology.bond_order(0), Some(2));
        assert_eq!(part.topology.bond_type(0), Some(BondType::Double));
    }

    #[test]
    fn the_coordinate_unit_is_one_hundredth_of_a_picometre() {
        let document = import_mmp(SAMPLE).expect("import");
        let position = document
            .part(0)
            .and_then(|part| part.topology.position_m(1))
            .expect("position");
        assert_eq!(position[0], 1200.0 * 1.0e-13);
    }

    #[test]
    fn an_mmp_text_round_trips_through_export() {
        let document = import_mmp(SAMPLE).expect("import");
        let text = export_mmp(&document).expect("export");
        let again = import_mmp(&text).expect("reimport");
        assert_eq!(again, document);
    }

    #[test]
    fn an_empty_text_is_rejected() {
        assert!(matches!(import_mmp(""), Err(FormatError::MmpEmpty)));
    }

    #[test]
    fn a_bond_to_an_unknown_atom_is_rejected() {
        let text = "mol (bad) def\natom 1 (6) (0, 0, 0) def\nbond1 99\n";
        assert!(matches!(
            import_mmp(text),
            Err(FormatError::MmpParse { .. })
        ));
    }

    #[test]
    fn a_malformed_atom_is_rejected_not_panicking() {
        let text = "mol (bad) def\natom 1 (6) (oops, 0, 0) def\n";
        assert!(matches!(
            import_mmp(text),
            Err(FormatError::MmpParse { .. })
        ));
    }
}
