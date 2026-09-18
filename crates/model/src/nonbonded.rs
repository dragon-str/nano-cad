//! Element nonbonded parameters: well depth and contact distance.
//!
//! A pocket that binds one molecule and rejects another needs a repulsion
//! term with a size for every element. This module holds that size in one
//! table, so that the pocket geometry and the binding metric read the same
//! numbers.
//!
//! The values are the nonbonded parameters of the Universal Force Field
//! (UFF), from Rappe, Casewit, Colwell, Goddard and Skiff, *Journal of the
//! American Chemical Society* 1992, 114, 10024. They were transcribed
//! through the RDKit implementation, which is BSD-3 licensed. The two
//! numbers per element are `x1`, the distance at the energy minimum, and
//! `D1`, the depth of the well. The values are data, not fitted constants.
//!
//! UFF gives every atom type of an element the same `x1` and `D1`, so the
//! table reduces to one row for each element.

use crate::element::Element;

/// The nonbonded parameters of one element.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NonbondedParams {
    /// The well depth in joules.
    pub well_depth_j: f64,
    /// The distance at the minimum of the well in metres.
    pub vdw_distance_m: f64,
}

/// The table entry for one atomic number.
struct Entry {
    atomic_number: u8,
    /// The UFF `x1` value in Angstrom.
    vdw_distance_a: f64,
    /// The UFF `D1` value in kilocalorie per mole.
    well_depth_kcal_per_mol: f64,
}

/// One kilocalorie per mole in joules.
pub const KILOCALORIE_PER_MOL_J: f64 = 4184.0 / 6.022_140_76e23;

/// The Lennard-Jones length conversion, `2^(1/6)`.
///
/// UFF states its distance at the minimum. The Lennard-Jones `sigma` is the
/// distance at which the potential crosses zero, which is smaller by this
/// factor.
pub const LENNARD_JONES_LENGTH_FACTOR: f64 = 1.122_462_048_309_373;

const ENTRIES: &[Entry] = &[
    Entry {
        atomic_number: 1,
        vdw_distance_a: 2.886,
        well_depth_kcal_per_mol: 0.044,
    },
    Entry {
        atomic_number: 6,
        vdw_distance_a: 3.851,
        well_depth_kcal_per_mol: 0.105,
    },
    Entry {
        atomic_number: 7,
        vdw_distance_a: 3.660,
        well_depth_kcal_per_mol: 0.069,
    },
    Entry {
        atomic_number: 8,
        vdw_distance_a: 3.500,
        well_depth_kcal_per_mol: 0.060,
    },
    Entry {
        atomic_number: 9,
        vdw_distance_a: 3.364,
        well_depth_kcal_per_mol: 0.050,
    },
    Entry {
        atomic_number: 15,
        vdw_distance_a: 4.147,
        well_depth_kcal_per_mol: 0.305,
    },
    Entry {
        atomic_number: 16,
        vdw_distance_a: 4.035,
        well_depth_kcal_per_mol: 0.274,
    },
    Entry {
        atomic_number: 17,
        vdw_distance_a: 3.947,
        well_depth_kcal_per_mol: 0.227,
    },
];

/// Returns the nonbonded parameters of an element, when the table holds it.
///
/// The table holds the elements that an organic molecule needs: hydrogen,
/// carbon, nitrogen, oxygen, fluorine, phosphorus, sulfur and chlorine. Any
/// other element returns `None`, so a caller cannot use a silent default.
pub fn nonbonded(element: Element) -> Option<NonbondedParams> {
    let number = element.atomic_number();
    let entry = ENTRIES.iter().find(|entry| entry.atomic_number == number)?;
    Some(NonbondedParams {
        well_depth_j: entry.well_depth_kcal_per_mol * KILOCALORIE_PER_MOL_J,
        vdw_distance_m: entry.vdw_distance_a * 1.0e-10,
    })
}

/// Returns the well depth of an element in joules, when the table holds it.
pub fn well_depth_j(element: Element) -> Option<f64> {
    nonbonded(element).map(|params| params.well_depth_j)
}

/// Returns the distance at the well minimum in metres, when the table holds it.
pub fn vdw_distance_m(element: Element) -> Option<f64> {
    nonbonded(element).map(|params| params.vdw_distance_m)
}

/// Returns the Lennard-Jones pair parameters for two elements.
///
/// The result is `(epsilon_j, sigma_m)`. UFF mixes the two elements with
/// `x_ij = sqrt(x_i x_j)` and `D_ij = sqrt(D_i D_j)`. The Lennard-Jones form
/// `U(r) = 4 epsilon ((sigma/r)^12 - (sigma/r)^6)` then needs
/// `epsilon = D_ij / 4` and `sigma = x_ij / 2^(1/6)`, because UFF states the
/// distance at the minimum of the well.
pub fn pair_params(left: Element, right: Element) -> Option<(f64, f64)> {
    let left = nonbonded(left)?;
    let right = nonbonded(right)?;
    let well_depth_j = (left.well_depth_j * right.well_depth_j).sqrt();
    let vdw_distance_m = (left.vdw_distance_m * right.vdw_distance_m).sqrt();
    Some((
        well_depth_j / 4.0,
        vdw_distance_m / LENNARD_JONES_LENGTH_FACTOR,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(left: f64, right: f64, tolerance: f64) -> bool {
        (left - right).abs() <= tolerance * right.abs().max(1.0)
    }

    #[test]
    fn the_table_holds_the_group_elements() {
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
            assert!(nonbonded(element).is_some(), "{element} is missing");
        }
    }

    #[test]
    fn an_unlisted_element_has_no_parameters() {
        assert!(nonbonded(Element::GOLD).is_none());
        assert!(pair_params(Element::CARBON, Element::GOLD).is_none());
        assert!(well_depth_j(Element::GOLD).is_none());
        assert!(vdw_distance_m(Element::GOLD).is_none());
    }

    #[test]
    fn the_carbon_values_match_the_source() {
        let carbon = nonbonded(Element::CARBON).unwrap();
        assert!(close(carbon.vdw_distance_m, 3.851e-10, 1.0e-12));
        assert!(close(
            carbon.well_depth_j,
            0.105 * KILOCALORIE_PER_MOL_J,
            1.0e-12
        ));
    }

    #[test]
    fn the_carbon_well_is_about_a_tenth_of_a_kilocalorie() {
        let carbon = well_depth_j(Element::CARBON).unwrap();
        assert!(close(carbon, 7.3e-22, 5.0e-2), "{carbon}");
    }

    #[test]
    fn the_table_holds_positive_values() {
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
            let params = nonbonded(element).unwrap();
            assert!(params.well_depth_j > 0.0);
            assert!(params.vdw_distance_m > 0.0);
        }
    }

    #[test]
    fn a_pair_of_one_element_recovers_its_values() {
        let (epsilon_j, sigma_m) = pair_params(Element::CARBON, Element::CARBON).unwrap();
        let carbon = nonbonded(Element::CARBON).unwrap();
        assert!(close(epsilon_j, carbon.well_depth_j / 4.0, 1.0e-12));
        assert!(close(
            sigma_m * LENNARD_JONES_LENGTH_FACTOR,
            carbon.vdw_distance_m,
            1.0e-12
        ));
    }

    #[test]
    fn the_mixing_rule_uses_the_geometric_mean() {
        let (epsilon_j, sigma_m) = pair_params(Element::HYDROGEN, Element::OXYGEN).unwrap();
        let hydrogen = nonbonded(Element::HYDROGEN).unwrap();
        let oxygen = nonbonded(Element::OXYGEN).unwrap();
        let expected_well = (hydrogen.well_depth_j * oxygen.well_depth_j).sqrt() / 4.0;
        let expected_distance =
            (hydrogen.vdw_distance_m * oxygen.vdw_distance_m).sqrt() / LENNARD_JONES_LENGTH_FACTOR;
        assert!(close(epsilon_j, expected_well, 1.0e-12));
        assert!(close(sigma_m, expected_distance, 1.0e-12));
    }

    #[test]
    fn the_mixing_rule_does_not_depend_on_the_order() {
        let forward = pair_params(Element::CARBON, Element::SULFUR).unwrap();
        let backward = pair_params(Element::SULFUR, Element::CARBON).unwrap();
        assert_eq!(forward, backward);
    }

    #[test]
    fn the_sigma_is_below_the_uff_distance() {
        let (_, sigma_m) = pair_params(Element::CARBON, Element::CARBON).unwrap();
        let carbon = nonbonded(Element::CARBON).unwrap();
        assert!(sigma_m < carbon.vdw_distance_m);
        assert!(sigma_m > 0.8 * carbon.vdw_distance_m);
    }

    #[test]
    fn a_heavy_pair_has_a_deeper_well_than_a_light_pair() {
        let (light, _) = pair_params(Element::HYDROGEN, Element::HYDROGEN).unwrap();
        let (heavy, _) = pair_params(Element::CHLORINE, Element::CHLORINE).unwrap();
        assert!(heavy > light);
    }
}
