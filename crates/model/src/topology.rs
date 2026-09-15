use crate::atom::Atom;
use crate::bond::{Bond, BondType};
use crate::element::Element;
use crate::error::ModelError;

const COORDS_PER_ATOM: usize = 3;

/// Atoms and bonds as struct-of-arrays, addressed by integer index.
///
/// Positions live in one flat `f64` buffer with three values per atom. Bonds
/// hold atom indices, never pointers. Every accessor and mutator checks its
/// index and returns `Option` or `Result`; none of them panic.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Topology {
    elements: Vec<Element>,
    positions_m: Vec<f64>,
    charges_c: Vec<f64>,
    atom_types: Vec<String>,
    bond_u: Vec<u32>,
    bond_v: Vec<u32>,
    bond_orders: Vec<u8>,
    bond_types: Vec<BondType>,
}

impl Topology {
    /// Creates an empty topology.
    pub fn new() -> Self {
        Self::default()
    }

    /// Returns the number of atoms.
    pub fn atom_count(&self) -> usize {
        self.elements.len()
    }

    /// Returns the number of bonds.
    pub fn bond_count(&self) -> usize {
        self.bond_u.len()
    }

    /// Returns the number of atoms. An alias of [`Topology::atom_count`].
    pub fn len(&self) -> usize {
        self.atom_count()
    }

    /// Reports whether the topology has no atoms.
    pub fn is_empty(&self) -> bool {
        self.atom_count() == 0
    }

    /// Appends an atom and returns its index.
    pub fn add_atom(&mut self, atom: Atom) -> u32 {
        let index = self.elements.len() as u32;
        self.elements.push(atom.element);
        self.positions_m.extend_from_slice(&atom.position_m);
        self.charges_c.push(atom.charge_c);
        self.atom_types.push(atom.atom_type);
        index
    }

    /// Returns the element at an index.
    pub fn element(&self, index: usize) -> Option<Element> {
        self.elements.get(index).copied()
    }

    /// Returns the position at an index, in SI metres.
    pub fn position_m(&self, index: usize) -> Option<[f64; 3]> {
        let slice = self.position_slice(index)?;
        Some([slice[0], slice[1], slice[2]])
    }

    /// Returns the charge at an index, in SI coulombs.
    pub fn charge_c(&self, index: usize) -> Option<f64> {
        self.charges_c.get(index).copied()
    }

    /// Returns the atom type name at an index.
    pub fn atom_type(&self, index: usize) -> Option<&str> {
        self.atom_types.get(index).map(String::as_str)
    }

    /// Returns a full atom value at an index.
    pub fn atom(&self, index: usize) -> Option<Atom> {
        let element = self.element(index)?;
        let position_m = self.position_m(index)?;
        let charge_c = self.charge_c(index)?;
        let atom_type = self.atom_types.get(index)?.clone();
        Some(Atom::new(element, position_m, charge_c, atom_type))
    }

    /// Returns the flat element buffer.
    pub fn elements(&self) -> &[Element] {
        &self.elements
    }

    /// Returns the flat position buffer. Its length is three times the atom
    /// count.
    pub fn positions_m(&self) -> &[f64] {
        &self.positions_m
    }

    /// Returns the flat position buffer for in-place mutation.
    pub fn positions_m_mut(&mut self) -> &mut [f64] {
        &mut self.positions_m
    }

    /// Returns the flat charge buffer.
    pub fn charges_c(&self) -> &[f64] {
        &self.charges_c
    }

    /// Returns the flat atom-type-name buffer.
    pub fn atom_types(&self) -> &[String] {
        &self.atom_types
    }

    /// Replaces the element at an index.
    pub fn set_element(&mut self, index: usize, element: Element) -> Result<(), ModelError> {
        self.check_atom(index)?;
        if let Some(slot) = self.elements.get_mut(index) {
            *slot = element;
        }
        Ok(())
    }

    /// Replaces the position at an index, in SI metres.
    pub fn set_position_m(&mut self, index: usize, position_m: [f64; 3]) -> Result<(), ModelError> {
        self.check_atom(index)?;
        let start = index * COORDS_PER_ATOM;
        if let Some(slot) = self.positions_m.get_mut(start..start + COORDS_PER_ATOM) {
            slot.copy_from_slice(&position_m);
        }
        Ok(())
    }

    /// Replaces the charge at an index, in SI coulombs.
    pub fn set_charge_c(&mut self, index: usize, charge_c: f64) -> Result<(), ModelError> {
        self.check_atom(index)?;
        if let Some(slot) = self.charges_c.get_mut(index) {
            *slot = charge_c;
        }
        Ok(())
    }

    /// Replaces the atom type name at an index.
    pub fn set_atom_type(
        &mut self,
        index: usize,
        atom_type: impl Into<String>,
    ) -> Result<(), ModelError> {
        self.check_atom(index)?;
        if let Some(slot) = self.atom_types.get_mut(index) {
            *slot = atom_type.into();
        }
        Ok(())
    }

    /// Appends a bond. Both endpoints must name existing atoms. Returns the
    /// bond index.
    pub fn add_bond(&mut self, bond: Bond) -> Result<u32, ModelError> {
        let atom_count = self.atom_count();
        if bond.u as usize >= atom_count {
            return Err(ModelError::BondEndpointOutOfBounds {
                endpoint: bond.u,
                atom_count,
            });
        }
        if bond.v as usize >= atom_count {
            return Err(ModelError::BondEndpointOutOfBounds {
                endpoint: bond.v,
                atom_count,
            });
        }
        if !Bond::is_order_valid(bond.order) {
            return Err(ModelError::InvalidBondOrder(bond.order));
        }
        let index = self.bond_u.len() as u32;
        self.bond_u.push(bond.u);
        self.bond_v.push(bond.v);
        self.bond_orders.push(bond.order);
        self.bond_types.push(bond.bond_type);
        Ok(index)
    }

    /// Returns the bond at an index.
    pub fn bond(&self, index: usize) -> Option<Bond> {
        Some(Bond::new(
            *self.bond_u.get(index)?,
            *self.bond_v.get(index)?,
            *self.bond_orders.get(index)?,
            *self.bond_types.get(index)?,
        ))
    }

    /// Returns the first endpoint of a bond.
    pub fn bond_u(&self, index: usize) -> Option<u32> {
        self.bond_u.get(index).copied()
    }

    /// Returns the second endpoint of a bond.
    pub fn bond_v(&self, index: usize) -> Option<u32> {
        self.bond_v.get(index).copied()
    }

    /// Returns the order of a bond.
    pub fn bond_order(&self, index: usize) -> Option<u8> {
        self.bond_orders.get(index).copied()
    }

    /// Returns the kind of a bond.
    pub fn bond_type(&self, index: usize) -> Option<BondType> {
        self.bond_types.get(index).copied()
    }

    /// Replaces the order of a bond.
    pub fn set_bond_order(&mut self, index: usize, order: u8) -> Result<(), ModelError> {
        self.check_bond(index)?;
        if !Bond::is_order_valid(order) {
            return Err(ModelError::InvalidBondOrder(order));
        }
        if let Some(slot) = self.bond_orders.get_mut(index) {
            *slot = order;
        }
        Ok(())
    }

    /// Replaces the kind of a bond.
    pub fn set_bond_type(&mut self, index: usize, bond_type: BondType) -> Result<(), ModelError> {
        self.check_bond(index)?;
        if let Some(slot) = self.bond_types.get_mut(index) {
            *slot = bond_type;
        }
        Ok(())
    }

    /// Iterates over atom values in index order.
    pub fn atoms(&self) -> impl Iterator<Item = Atom> + '_ {
        (0..self.atom_count()).filter_map(|index| self.atom(index))
    }

    /// Iterates over bond values in index order.
    pub fn bonds(&self) -> impl Iterator<Item = Bond> + '_ {
        (0..self.bond_count()).filter_map(|index| self.bond(index))
    }

    fn position_slice(&self, index: usize) -> Option<&[f64]> {
        if index >= self.atom_count() {
            return None;
        }
        let start = index * COORDS_PER_ATOM;
        self.positions_m.get(start..start + COORDS_PER_ATOM)
    }

    fn check_atom(&self, index: usize) -> Result<(), ModelError> {
        if index < self.atom_count() {
            Ok(())
        } else {
            Err(ModelError::AtomIndexOutOfBounds {
                index,
                len: self.atom_count(),
            })
        }
    }

    fn check_bond(&self, index: usize) -> Result<(), ModelError> {
        if index < self.bond_count() {
            Ok(())
        } else {
            Err(ModelError::BondIndexOutOfBounds {
                index,
                len: self.bond_count(),
            })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn atom(element: Element, x: f64, atom_type: &str) -> Atom {
        Atom::new(element, [x, 0.0, 0.0], 0.0, atom_type)
    }

    #[test]
    fn a_new_topology_is_empty() {
        let topology = Topology::new();
        assert!(topology.is_empty());
        assert_eq!(topology.len(), 0);
        assert_eq!(topology.atom_count(), 0);
        assert_eq!(topology.bond_count(), 0);
    }

    #[test]
    fn add_atom_returns_sequential_indices_and_stores_every_field() {
        let mut topology = Topology::new();
        assert_eq!(topology.add_atom(atom(Element::CARBON, 1.0e-10, "C3")), 0);
        assert_eq!(topology.add_atom(atom(Element::OXYGEN, 2.0e-10, "O2")), 1);
        assert_eq!(topology.atom_count(), 2);
        assert_eq!(topology.element(0), Some(Element::CARBON));
        assert_eq!(topology.element(1), Some(Element::OXYGEN));
        assert_eq!(topology.position_m(1), Some([2.0e-10, 0.0, 0.0]));
        assert_eq!(topology.atom_type(0), Some("C3"));
    }

    #[test]
    fn positions_are_stored_flat_with_three_values_per_atom() {
        let mut topology = Topology::new();
        topology.add_atom(Atom::new(Element::HYDROGEN, [1.0, 2.0, 3.0], 0.0, "H"));
        topology.add_atom(Atom::new(Element::HYDROGEN, [4.0, 5.0, 6.0], 0.0, "H"));
        assert_eq!(topology.positions_m(), &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        assert_eq!(topology.positions_m().len(), topology.atom_count() * 3);
    }

    #[test]
    fn atom_access_out_of_bounds_is_none_not_a_panic() {
        let topology = Topology::new();
        assert_eq!(topology.atom(0), None);
        assert_eq!(topology.element(7), None);
        assert_eq!(topology.position_m(7), None);
        assert_eq!(topology.charge_c(7), None);
        assert_eq!(topology.atom_type(7), None);
        assert_eq!(topology.bond(7), None);
    }

    #[test]
    fn mutators_reject_an_out_of_bounds_index() {
        let mut topology = Topology::new();
        assert!(topology.set_position_m(0, [0.0, 0.0, 0.0]).is_err());
        assert!(topology.set_charge_c(0, 1.0).is_err());
        assert!(topology.set_element(0, Element::HYDROGEN).is_err());
        assert!(topology.set_atom_type(0, "H").is_err());
        assert!(topology.set_bond_order(0, 1).is_err());
        assert!(topology.set_bond_type(0, BondType::Single).is_err());
    }

    #[test]
    fn a_in_bounds_mutator_changes_only_its_target() {
        let mut topology = Topology::new();
        topology.add_atom(atom(Element::CARBON, 1.0, "C3"));
        topology.add_atom(atom(Element::CARBON, 2.0, "C3"));
        topology
            .set_position_m(1, [9.0, 9.0, 9.0])
            .expect("in bounds");
        topology.set_charge_c(1, -0.5).expect("in bounds");
        assert_eq!(topology.position_m(0), Some([1.0, 0.0, 0.0]));
        assert_eq!(topology.position_m(1), Some([9.0, 9.0, 9.0]));
        assert_eq!(topology.charge_c(1), Some(-0.5));
    }

    #[test]
    fn atom_values_round_trip_through_the_topology() {
        let mut topology = Topology::new();
        let carbon = Atom::new(Element::CARBON, [1.0e-10, 2.0e-10, 3.0e-10], -0.2, "C3");
        topology.add_atom(carbon.clone());
        assert_eq!(topology.atom(0), Some(carbon));
    }

    #[test]
    fn add_bond_returns_an_index_and_stores_every_field() {
        let mut topology = Topology::new();
        topology.add_atom(atom(Element::CARBON, 1.0, "C3"));
        topology.add_atom(atom(Element::CARBON, 2.0, "C3"));
        let index = topology
            .add_bond(Bond::new(0, 1, 1, BondType::Single))
            .expect("valid bond");
        assert_eq!(index, 0);
        assert_eq!(topology.bond_count(), 1);
        assert_eq!(topology.bond(0), Some(Bond::new(0, 1, 1, BondType::Single)));
        assert_eq!(topology.bond_u(0), Some(0));
        assert_eq!(topology.bond_v(0), Some(1));
        assert_eq!(topology.bond_order(0), Some(1));
        assert_eq!(topology.bond_type(0), Some(BondType::Single));
    }

    #[test]
    fn add_bond_rejects_an_endpoint_that_names_no_atom() {
        let mut topology = Topology::new();
        topology.add_atom(atom(Element::CARBON, 1.0, "C3"));
        assert_eq!(
            topology.add_bond(Bond::new(0, 1, 1, BondType::Single)),
            Err(ModelError::BondEndpointOutOfBounds {
                endpoint: 1,
                atom_count: 1,
            })
        );
        assert_eq!(
            topology.add_bond(Bond::new(5, 0, 1, BondType::Single)),
            Err(ModelError::BondEndpointOutOfBounds {
                endpoint: 5,
                atom_count: 1,
            })
        );
        assert_eq!(topology.bond_count(), 0);
    }

    #[test]
    fn add_bond_rejects_an_invalid_order() {
        let mut topology = Topology::new();
        topology.add_atom(atom(Element::CARBON, 1.0, "C3"));
        topology.add_atom(atom(Element::CARBON, 2.0, "C3"));
        assert_eq!(
            topology.add_bond(Bond::new(0, 1, 9, BondType::Single)),
            Err(ModelError::InvalidBondOrder(9))
        );
    }

    #[test]
    fn set_bond_order_and_type_validate_their_input() {
        let mut topology = Topology::new();
        topology.add_atom(atom(Element::CARBON, 1.0, "C3"));
        topology.add_atom(atom(Element::CARBON, 2.0, "C3"));
        topology
            .add_bond(Bond::new(0, 1, 1, BondType::Single))
            .expect("valid bond");
        topology.set_bond_order(0, 2).expect("valid order");
        topology.set_bond_type(0, BondType::Double).expect("valid");
        assert_eq!(topology.bond_order(0), Some(2));
        assert_eq!(topology.bond_type(0), Some(BondType::Double));
        assert!(topology.set_bond_order(0, 200).is_err());
    }

    #[test]
    fn iterators_visit_every_atom_and_bond_in_order() {
        let mut topology = Topology::new();
        topology.add_atom(atom(Element::CARBON, 1.0, "C3"));
        topology.add_atom(atom(Element::OXYGEN, 2.0, "O2"));
        topology
            .add_bond(Bond::new(0, 1, 1, BondType::Single))
            .expect("valid bond");
        let elements: Vec<Element> = topology.atoms().map(|a| a.element).collect();
        assert_eq!(elements, vec![Element::CARBON, Element::OXYGEN]);
        assert_eq!(topology.bonds().count(), 1);
    }

    #[test]
    fn the_mutable_position_buffer_has_the_expected_layout() {
        let mut topology = Topology::new();
        topology.add_atom(atom(Element::HYDROGEN, 0.0, "H"));
        topology.positions_m_mut().copy_from_slice(&[7.0, 8.0, 9.0]);
        assert_eq!(topology.position_m(0), Some([7.0, 8.0, 9.0]));
    }
}
