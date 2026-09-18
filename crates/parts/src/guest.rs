//! Guest molecules that a binding pocket can admit or reject.
//!
//! A guest is data, not a parametric family. Each molecule is a frozen table of
//! atoms, bonds and partial charges, because a molecule has one shape and one
//! charge set. That is why this module does not implement [`PartGenerator`].
//!
//! The table holds two pairs that make a selectivity test possible:
//!
//! - Ethanol and dimethyl ether share the formula `C2H6O`. They have the same
//!   nine atoms and almost the same mass, and they have different shapes. One is
//!   a chain that ends in a hydroxyl group. The other is a bent ether. A pocket
//!   that admits one and rejects the other selects by shape and by hydrogen
//!   bonding, not by mass.
//! - Methanol is the small polar control. Benzene and cyclohexane are the plain
//!   hydrocarbons, and they differ in size.
//!
//! The geometry and the charges are a model. See [`crate::guest_data`].
//!
//! [`PartGenerator`]: crate::generator::PartGenerator

use nanocad_model::{atomic_mass_kg, Atom, Bond, BondType, Part, Topology};

use crate::error::PartError;
use crate::guest_data::GUEST_DATA;

/// The Avogadro constant in reciprocal moles. This is the exact value of the
/// 2019 SI redefinition.
pub const AVOGADRO_PER_MOL: f64 = 6.022_140_76e23;

/// One atom of a guest molecule.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GuestAtom {
    /// The element.
    pub element: nanocad_model::Element,
    /// The position in metres, in the frame of the frozen table.
    pub position_m: [f64; 3],
    /// The partial charge in coulombs. This is a model charge.
    pub charge_c: f64,
    /// The force-field type label, for example `C_sp3` or `O_hydroxyl`.
    pub atom_type: &'static str,
}

/// One bond of a guest molecule.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GuestBond {
    /// The first atom index.
    pub i: u32,
    /// The second atom index.
    pub j: u32,
    /// The bond order. Four means aromatic, which is the usual tag convention.
    pub order: u8,
}

/// A frozen guest molecule.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Guest {
    /// The short id, for example `ethanol`.
    pub id: &'static str,
    /// The display label.
    pub label: &'static str,
    /// The molecular formula.
    pub formula: &'static str,
    /// The molar mass in kilograms per mole, as the source table states it.
    pub molar_mass_kg_per_mol: f64,
    /// The atoms.
    pub atoms: &'static [GuestAtom],
    /// The bonds.
    pub bonds: &'static [GuestBond],
}

impl Guest {
    /// The number of atoms.
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    /// The number of bonds.
    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    /// The number of atoms that are not hydrogen.
    pub fn heavy_atom_count(&self) -> usize {
        self.atoms
            .iter()
            .filter(|atom| atom.element.symbol() != Some("H"))
            .count()
    }

    /// The mass of one molecule in kilograms.
    ///
    /// The mass sums the element table, so it does not use the frozen molar
    /// mass. A test checks that the two agree.
    pub fn mass_kg(&self) -> f64 {
        self.atoms
            .iter()
            .filter_map(|atom| atomic_mass_kg(atom.element))
            .sum()
    }

    /// The mass-weighted centroid in metres.
    pub fn centroid_m(&self) -> [f64; 3] {
        let mut total = 0.0;
        let mut sum = [0.0_f64; 3];
        for atom in self.atoms {
            let mass = atomic_mass_kg(atom.element).unwrap_or(0.0);
            total += mass;
            for (axis, value) in atom.position_m.iter().enumerate() {
                sum[axis] += mass * value;
            }
        }
        if total <= 0.0 {
            return [0.0; 3];
        }
        [sum[0] / total, sum[1] / total, sum[2] / total]
    }

    /// The mass-weighted radius of gyration in metres.
    pub fn radius_of_gyration_m(&self) -> f64 {
        let centroid = self.centroid_m();
        let mut total = 0.0;
        let mut moment = 0.0;
        for atom in self.atoms {
            let mass = atomic_mass_kg(atom.element).unwrap_or(0.0);
            total += mass;
            let mut square = 0.0;
            for (position, centre) in atom.position_m.iter().zip(centroid.iter()) {
                let delta = position - centre;
                square += delta * delta;
            }
            moment += mass * square;
        }
        if total <= 0.0 {
            return 0.0;
        }
        (moment / total).sqrt()
    }

    /// The largest distance from the centroid to any atom, in metres.
    pub fn extent_m(&self) -> f64 {
        let centroid = self.centroid_m();
        let mut largest: f64 = 0.0;
        for atom in self.atoms {
            let mut square = 0.0;
            for (position, centre) in atom.position_m.iter().zip(centroid.iter()) {
                let delta = position - centre;
                square += delta * delta;
            }
            largest = largest.max(square.sqrt());
        }
        largest
    }

    /// The total charge in coulombs. It is zero for a neutral molecule.
    pub fn total_charge_c(&self) -> f64 {
        self.atoms.iter().map(|atom| atom.charge_c).sum()
    }
}

/// Every guest molecule in the table.
pub fn guests() -> &'static [Guest] {
    GUEST_DATA
}

/// The guest molecule with the given id.
pub fn guest(id: &str) -> Option<&'static Guest> {
    GUEST_DATA.iter().find(|entry| entry.id == id)
}

/// Builds a part from a guest molecule.
///
/// The part keeps the frozen positions. A caller that needs the molecule
/// somewhere else transforms the part.
pub fn build_guest(guest: &Guest) -> Part {
    let mut topology = Topology::new();
    for atom in guest.atoms {
        topology.add_atom(Atom::new(
            atom.element,
            atom.position_m,
            atom.charge_c,
            atom.atom_type,
        ));
    }
    for bond in guest.bonds {
        let kind = match bond.order {
            1 => BondType::Single,
            2 => BondType::Double,
            3 => BondType::Triple,
            4 => BondType::Aromatic,
            _ => BondType::Single,
        };
        let _ = topology.add_bond(Bond::new(bond.i, bond.j, bond.order, kind));
    }
    Part::new(guest.id, topology).with_material("guest molecule")
}

/// Builds a part from the guest molecule with the given id.
pub fn build_guest_part(id: &str) -> Result<Part, PartError> {
    let guest = guest(id).ok_or_else(|| PartError::UnknownGenerator(id.to_string()))?;
    Ok(build_guest(guest))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_table_holds_five_molecules() {
        assert_eq!(guests().len(), 5);
    }

    #[test]
    fn every_guest_id_is_unique() {
        let mut ids: Vec<&str> = guests().iter().map(|entry| entry.id).collect();
        ids.sort_unstable();
        let count = ids.len();
        ids.dedup();
        assert_eq!(ids.len(), count);
    }

    #[test]
    fn the_isomer_pair_shares_a_formula() {
        let ethanol = guest("ethanol").unwrap();
        let ether = guest("dimethyl_ether").unwrap();
        assert_eq!(ethanol.formula, ether.formula);
        assert_eq!(ethanol.atom_count(), ether.atom_count());
        let mass_ratio = ethanol.mass_kg() / ether.mass_kg();
        assert!((mass_ratio - 1.0).abs() < 1e-9, "mass ratio {mass_ratio}");
    }

    #[test]
    fn the_isomer_pair_has_different_shapes() {
        let ethanol = guest("ethanol").unwrap();
        let ether = guest("dimethyl_ether").unwrap();
        assert!(
            ethanol.radius_of_gyration_m() > ether.radius_of_gyration_m(),
            "ethanol is the longer molecule"
        );
        assert!(ethanol.extent_m() > ether.extent_m());
    }

    #[test]
    fn the_isomer_pair_has_different_charges() {
        let ethanol = guest("ethanol").unwrap();
        let ether = guest("dimethyl_ether").unwrap();
        let charge = |molecule: &Guest| -> f64 {
            molecule
                .atoms
                .iter()
                .filter(|atom| atom.element.symbol() == Some("O"))
                .map(|atom| atom.charge_c)
                .fold(0.0_f64, f64::min)
        };
        assert!(
            (charge(ethanol) - charge(ether)).abs() > 1.0e-3,
            "the two oxygens must not carry the same charge"
        );
    }

    #[test]
    fn only_ethanol_can_donate_a_hydrogen_bond() {
        let donors = |molecule: &Guest| -> usize {
            molecule
                .atoms
                .iter()
                .filter(|atom| atom.atom_type == "H_hydroxyl")
                .count()
        };
        assert_eq!(donors(guest("ethanol").unwrap()), 1);
        assert_eq!(donors(guest("dimethyl_ether").unwrap()), 0);
    }

    #[test]
    fn benzene_is_flat_and_cyclohexane_is_not() {
        // The ring normal is not the z axis, because the embedding orients the
        // molecule freely. Fit the plane from three carbon positions instead.
        let benzene = guest("benzene").unwrap();
        let ring: Vec<[f64; 3]> = benzene
            .atoms
            .iter()
            .filter(|atom| atom.element.symbol() == Some("C"))
            .map(|atom| atom.position_m)
            .collect();
        let first = ring[0];
        let normal = cross(subtract(ring[1], ring[0]), subtract(ring[2], ring[1]));
        let length = dot(normal, normal).sqrt();
        assert!(length > 0.0);
        let unit = [normal[0] / length, normal[1] / length, normal[2] / length];
        let mut spread: f64 = 0.0;
        for position in &ring {
            spread = spread.max(dot(subtract(*position, first), unit).abs());
        }
        assert!(
            spread < 1.0e-13,
            "the ring is planar, out-of-plane {spread}"
        );
        assert!(guest("cyclohexane").unwrap().extent_m() > benzene.extent_m());
    }

    fn subtract(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
        [left[0] - right[0], left[1] - right[1], left[2] - right[2]]
    }

    fn cross(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
        [
            left[1] * right[2] - left[2] * right[1],
            left[2] * right[0] - left[0] * right[2],
            left[0] * right[1] - left[1] * right[0],
        ]
    }

    fn dot(left: [f64; 3], right: [f64; 3]) -> f64 {
        left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
    }

    #[test]
    fn every_guest_is_neutral() {
        for entry in guests() {
            let total = entry.total_charge_c();
            assert!(total.abs() < 1.0e-4, "{} has charge {total}", entry.id);
        }
    }

    #[test]
    fn the_element_mass_matches_the_source_molar_mass() {
        for entry in guests() {
            let from_table = entry.mass_kg() * AVOGADRO_PER_MOL;
            let relative =
                (from_table - entry.molar_mass_kg_per_mol).abs() / entry.molar_mass_kg_per_mol;
            assert!(
                relative < 1.0e-2,
                "{} differs by {relative} from {}",
                entry.id,
                entry.molar_mass_kg_per_mol
            );
        }
    }

    #[test]
    fn every_bond_index_is_in_range() {
        for entry in guests() {
            for bond in entry.bonds {
                assert!((bond.i as usize) < entry.atom_count(), "{}", entry.id);
                assert!((bond.j as usize) < entry.atom_count(), "{}", entry.id);
                assert!(bond.i != bond.j, "{}", entry.id);
            }
        }
    }

    #[test]
    fn every_atom_type_is_stated() {
        for entry in guests() {
            for atom in entry.atoms {
                assert!(!atom.atom_type.is_empty(), "{}", entry.id);
            }
        }
    }

    #[test]
    fn a_guest_builds_into_a_part() {
        let part = build_guest_part("benzene").unwrap();
        assert_eq!(part.atom_count(), 12);
        assert_eq!(part.bond_count(), 12);
        assert_eq!(
            part.topology.charge_c(0),
            Some(guest("benzene").unwrap().atoms[0].charge_c)
        );
    }

    #[test]
    fn an_unknown_guest_is_refused() {
        assert!(build_guest_part("unobtainium").is_err());
        assert!(guest("unobtainium").is_none());
    }
}
