//! XYZ import and export.
//!
//! XYZ is a plain-text format. Line one gives the atom count, line two is a
//! comment, and each later line gives an element symbol and three coordinates
//! in Angstroms. The model stores SI metres, so this module converts at the
//! boundary.

use nanocad_model::{Atom, Document, Element, Part, Topology};

use crate::error::FormatError;
use crate::periodic;

const METRES_PER_ANGSTROM: f64 = 1.0e-10;
const ANGSTROMS_PER_METRE: f64 = 1.0e10;

/// Imports an XYZ text into a one-part document.
pub fn import_xyz(text: &str) -> Result<Document, FormatError> {
    let mut lines = text.lines().enumerate();
    let (_, count_line) = lines.next().ok_or_else(|| FormatError::XyzParse {
        line: 1,
        message: "missing count line".to_owned(),
    })?;
    let atom_count: usize = count_line
        .trim()
        .parse()
        .map_err(|_| FormatError::XyzParse {
            line: 1,
            message: "count line is not a number".to_owned(),
        })?;

    let (_, comment) = lines.next().unwrap_or((1, ""));
    let part_name = if comment.trim().is_empty() {
        "xyz".to_owned()
    } else {
        comment.trim().to_owned()
    };

    let mut topology = Topology::new();
    for (index, raw_line) in lines {
        if topology.atom_count() == atom_count {
            break;
        }
        let line_number = index + 1;
        let line = raw_line.trim();
        if line.is_empty() {
            continue;
        }
        let mut tokens = line.split_whitespace();
        let symbol = tokens.next().ok_or_else(|| FormatError::XyzParse {
            line: line_number,
            message: "missing element symbol".to_owned(),
        })?;
        let atomic_number = periodic::atomic_number_for(symbol)
            .ok_or_else(|| FormatError::UnknownElementSymbol(symbol.to_owned()))?;
        let element = Element::from_atomic_number(atomic_number)
            .ok_or(FormatError::UnknownElement(atomic_number))?;
        let mut coordinate = [0.0_f64; 3];
        for (axis, slot) in coordinate.iter_mut().enumerate() {
            let token = tokens.next().ok_or_else(|| FormatError::XyzParse {
                line: line_number,
                message: format!("missing coordinate {axis}"),
            })?;
            *slot = token.parse::<f64>().map_err(|_| FormatError::XyzParse {
                line: line_number,
                message: format!("coordinate {axis} is not a number"),
            })? * METRES_PER_ANGSTROM;
        }
        topology.add_atom(Atom::new(element, coordinate, 0.0, symbol.to_owned()));
    }

    if topology.atom_count() != atom_count {
        return Err(FormatError::XyzCountMismatch);
    }

    let part = Part::new(part_name.clone(), topology);
    let mut document = Document::new(part_name);
    document.add_part(part);
    Ok(document)
}

/// Exports every atom of a document as one XYZ block.
pub fn export_xyz(document: &Document) -> Result<String, FormatError> {
    let total = document.parts.iter().map(Part::atom_count).sum::<usize>();
    let comment = if document.name.is_empty() {
        "nanocad"
    } else {
        document.name.as_str()
    };
    let mut out = String::new();
    out.push_str(&format!("{total}\n{comment}\n"));
    for part in &document.parts {
        for atom in part.topology.atoms() {
            let symbol = periodic::symbol_for(atom.element.atomic_number())
                .map(str::to_owned)
                .unwrap_or_else(|| format!("Z{}", atom.element.atomic_number()));
            let x = atom.position_m[0] * ANGSTROMS_PER_METRE;
            let y = atom.position_m[1] * ANGSTROMS_PER_METRE;
            let z = atom.position_m[2] * ANGSTROMS_PER_METRE;
            out.push_str(&format!("{symbol} {x} {y} {z}\n"));
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nanocad_model::Element;

    #[test]
    fn an_xyz_text_round_trips_through_a_document() {
        let text = "2\nwater\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\n";
        let document = import_xyz(text).expect("import");
        assert_eq!(document.part_count(), 1);
        let part = document.part(0).expect("part");
        assert_eq!(part.atom_count(), 2);
        assert_eq!(part.topology.element(0), Some(Element::OXYGEN));
        assert_eq!(part.topology.element(1), Some(Element::HYDROGEN));
        let restored = export_xyz(&document).expect("export");
        let again = import_xyz(&restored).expect("reimport");
        let original_position = document
            .part(0)
            .and_then(|part| part.topology.position_m(1))
            .expect("original");
        let restored_position = again
            .part(0)
            .and_then(|part| part.topology.position_m(1))
            .expect("restored");
        for axis in 0..3 {
            assert!((original_position[axis] - restored_position[axis]).abs() < 1.0e-15);
        }
    }

    #[test]
    fn a_count_that_does_not_match_is_rejected() {
        let text = "3\nbad\nH 0.0 0.0 0.0\n";
        assert!(matches!(
            import_xyz(text),
            Err(FormatError::XyzCountMismatch)
        ));
    }

    #[test]
    fn an_unknown_symbol_is_rejected() {
        let text = "1\nbad\nXx 0.0 0.0 0.0\n";
        assert!(matches!(
            import_xyz(text),
            Err(FormatError::UnknownElementSymbol(_))
        ));
    }

    #[test]
    fn a_broken_coordinate_is_rejected() {
        let text = "1\nbad\nH nope 0.0 0.0\n";
        assert!(matches!(
            import_xyz(text),
            Err(FormatError::XyzParse { .. })
        ));
    }

    #[test]
    fn exported_coordinates_use_angstroms() {
        let mut topology = Topology::new();
        topology.add_atom(Atom::new(Element::HYDROGEN, [1.0e-10, 0.0, 0.0], 0.0, "H"));
        let mut document = Document::new("one");
        document.add_part(Part::new("one", topology));
        let text = export_xyz(&document).expect("export");
        assert!(text.contains("H 1 0 0"), "text was {text:?}");
    }
}
