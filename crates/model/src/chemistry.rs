//! Element chemistry: valence, atomic mass and covalent radius.
//!
//! A generator that builds a part from more than one element needs three
//! facts about each element: the number of bonds it forms, the length of
//! those bonds, and the mass it carries. This module holds those facts in
//! one table, so that the lattice code and the mass code read the same
//! numbers.
//!
//! The valence is the usual valence of the neutral element in an organic
//! molecule. The mass is the standard atomic weight of IUPAC. The covalent
//! radius is the single-bond radius of Cordero and co-workers (2008). All
//! values are data, not fitted constants.

use crate::element::Element;

/// The chemistry of one element.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElementChemistry {
    /// The number of single bonds this element usually forms.
    pub valence: u8,
    /// The standard atomic weight in kilograms.
    pub atomic_mass_kg: f64,
    /// The single-bond covalent radius in metres.
    pub covalent_radius_m: f64,
}

/// The table entry for one atomic number.
struct Entry {
    atomic_number: u8,
    valence: u8,
    atomic_mass_u: f64,
    covalent_radius_m: f64,
}

/// The unified atomic mass unit in kilograms.
pub const ATOMIC_MASS_UNIT_KG: f64 = 1.660_539_066_60e-27;

const ENTRIES: &[Entry] = &[
    Entry {
        atomic_number: 1,
        valence: 1,
        atomic_mass_u: 1.008,
        covalent_radius_m: 0.31e-10,
    },
    Entry {
        atomic_number: 6,
        valence: 4,
        atomic_mass_u: 12.011,
        covalent_radius_m: 0.76e-10,
    },
    Entry {
        atomic_number: 7,
        valence: 3,
        atomic_mass_u: 14.007,
        covalent_radius_m: 0.71e-10,
    },
    Entry {
        atomic_number: 8,
        valence: 2,
        atomic_mass_u: 15.999,
        covalent_radius_m: 0.66e-10,
    },
    Entry {
        atomic_number: 9,
        valence: 1,
        atomic_mass_u: 18.998,
        covalent_radius_m: 0.57e-10,
    },
    Entry {
        atomic_number: 15,
        valence: 3,
        atomic_mass_u: 30.974,
        covalent_radius_m: 1.07e-10,
    },
    Entry {
        atomic_number: 16,
        valence: 2,
        atomic_mass_u: 32.06,
        covalent_radius_m: 1.05e-10,
    },
    Entry {
        atomic_number: 17,
        valence: 1,
        atomic_mass_u: 35.45,
        covalent_radius_m: 1.02e-10,
    },
    Entry {
        atomic_number: 35,
        valence: 1,
        atomic_mass_u: 79.904,
        covalent_radius_m: 1.20e-10,
    },
];

/// Returns the chemistry of an element, when the table holds it.
///
/// The table holds the elements that an organic molecule needs: hydrogen,
/// carbon, nitrogen, oxygen, fluorine, phosphorus, sulfur, chlorine and
/// bromine. Any other element returns `None`.
pub fn chemistry(element: Element) -> Option<ElementChemistry> {
    let number = element.atomic_number();
    let entry = ENTRIES.iter().find(|entry| entry.atomic_number == number)?;
    Some(ElementChemistry {
        valence: entry.valence,
        atomic_mass_kg: entry.atomic_mass_u * ATOMIC_MASS_UNIT_KG,
        covalent_radius_m: entry.covalent_radius_m,
    })
}

/// Returns the usual valence of an element, when the table holds it.
pub fn valence(element: Element) -> Option<u8> {
    chemistry(element).map(|chemistry| chemistry.valence)
}

/// Returns the standard atomic weight in kilograms, when the table holds it.
pub fn atomic_mass_kg(element: Element) -> Option<f64> {
    chemistry(element).map(|chemistry| chemistry.atomic_mass_kg)
}

/// Returns the single-bond covalent radius in metres, when the table holds it.
pub fn covalent_radius_m(element: Element) -> Option<f64> {
    chemistry(element).map(|chemistry| chemistry.covalent_radius_m)
}

/// Returns the length of a single bond between two elements in metres.
///
/// The length is the sum of the two covalent radii. Two atoms of the same
/// element therefore give twice the radius of that element.
pub fn bond_length_m(left: Element, right: Element) -> Option<f64> {
    Some(covalent_radius_m(left)? + covalent_radius_m(right)?)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_table_holds_the_organic_elements() {
        for element in [
            Element::HYDROGEN,
            Element::CARBON,
            Element::NITROGEN,
            Element::OXYGEN,
            Element::FLUORINE,
            Element::PHOSPHORUS,
            Element::SULFUR,
            Element::CHLORINE,
        ] {
            assert!(chemistry(element).is_some(), "{element} is missing");
        }
    }

    #[test]
    fn the_valences_match_the_usual_chemistry() {
        assert_eq!(valence(Element::HYDROGEN), Some(1));
        assert_eq!(valence(Element::CARBON), Some(4));
        assert_eq!(valence(Element::NITROGEN), Some(3));
        assert_eq!(valence(Element::OXYGEN), Some(2));
        assert_eq!(valence(Element::SULFUR), Some(2));
        assert_eq!(valence(Element::PHOSPHORUS), Some(3));
    }

    #[test]
    fn an_unlisted_element_has_no_chemistry() {
        assert_eq!(valence(Element::HELIUM), None);
        assert_eq!(atomic_mass_kg(Element::GOLD), None);
        assert_eq!(covalent_radius_m(Element::GOLD), None);
    }

    #[test]
    fn the_carbon_mass_is_twelve_units() {
        let mass_kg = atomic_mass_kg(Element::CARBON).expect("carbon is in the table");
        let ratio = mass_kg / ATOMIC_MASS_UNIT_KG;
        assert!((ratio - 12.011).abs() < 1.0e-9, "ratio {ratio}");
    }

    #[test]
    fn a_bond_between_two_carbons_is_twice_the_radius() {
        let radius_m = covalent_radius_m(Element::CARBON).expect("carbon is in the table");
        let length_m = bond_length_m(Element::CARBON, Element::CARBON).expect("both are known");
        assert!((length_m - 2.0 * radius_m).abs() < 1.0e-18);
    }

    #[test]
    fn a_hydrogen_carbon_bond_is_near_the_known_length() {
        let length_m = bond_length_m(Element::CARBON, Element::HYDROGEN).expect("both are known");
        assert!(
            (length_m - 1.09e-10).abs() < 0.1e-10,
            "length {length_m} m is far from 1.09e-10 m"
        );
    }

    #[test]
    fn a_bond_with_an_unknown_element_is_refused() {
        assert_eq!(bond_length_m(Element::CARBON, Element::GOLD), None);
    }

    #[test]
    fn every_entry_is_valid() {
        for entry in ENTRIES {
            assert!(
                entry.valence >= 1,
                "element {} has no valence",
                entry.atomic_number
            );
            assert!(entry.atomic_mass_u > 0.0);
            assert!(entry.covalent_radius_m > 0.0);
        }
    }

    #[test]
    fn no_element_appears_twice() {
        for (index, entry) in ENTRIES.iter().enumerate() {
            for other in &ENTRIES[index + 1..] {
                assert_ne!(
                    entry.atomic_number, other.atomic_number,
                    "element {} appears twice",
                    entry.atomic_number
                );
            }
        }
    }
}
