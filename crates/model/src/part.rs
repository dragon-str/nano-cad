use std::collections::BTreeMap;

use crate::atom::Atom;
use crate::bond::{Bond, BondType};
use crate::element::Element;
use crate::encoding::{Reader, Writer};
use crate::error::ModelError;
use crate::topology::Topology;

#[derive(Clone, Debug, PartialEq)]
pub struct Part {
    pub name: String,
    pub topology: Topology,
    pub material: String,
    pub metadata: BTreeMap<String, String>,
}

impl Part {
    pub fn new(name: impl Into<String>, topology: Topology) -> Self {
        Self {
            name: name.into(),
            topology,
            material: String::new(),
            metadata: BTreeMap::new(),
        }
    }

    pub fn with_material(mut self, material: impl Into<String>) -> Self {
        self.material = material.into();
        self
    }

    pub fn atom_count(&self) -> usize {
        self.topology.atom_count()
    }

    pub fn bond_count(&self) -> usize {
        self.topology.bond_count()
    }

    pub fn to_bytes(&self) -> Result<Vec<u8>, ModelError> {
        let mut writer = Writer::new();
        self.encode(&mut writer)?;
        Ok(writer.into_bytes())
    }

    pub fn from_bytes(bytes: &[u8]) -> Result<Self, ModelError> {
        let mut reader = Reader::new(bytes);
        let part = Self::decode(&mut reader)?;
        if reader.remaining() != 0 {
            return Err(ModelError::TrailingBytes {
                remaining: reader.remaining(),
            });
        }
        Ok(part)
    }

    pub(crate) fn encode(&self, writer: &mut Writer) -> Result<(), ModelError> {
        writer.write_str(&self.name)?;
        writer.write_str(&self.material)?;
        writer.write_map(&self.metadata)?;

        writer.write_count(self.topology.atom_count())?;
        for atom in self.topology.atoms() {
            writer.write_u8(atom.element.atomic_number());
            writer.write_f64(atom.position_m[0]);
            writer.write_f64(atom.position_m[1]);
            writer.write_f64(atom.position_m[2]);
            writer.write_f64(atom.charge_c);
            writer.write_str(&atom.atom_type)?;
        }

        writer.write_count(self.topology.bond_count())?;
        for bond in self.topology.bonds() {
            writer.write_u32(bond.u);
            writer.write_u32(bond.v);
            writer.write_u8(bond.order);
            writer.write_u8(bond.bond_type.to_u8());
        }
        Ok(())
    }

    pub(crate) fn decode(reader: &mut Reader<'_>) -> Result<Self, ModelError> {
        let name = reader.read_str()?;
        let material = reader.read_str()?;
        let metadata = reader.read_map()?;

        let atom_count = reader.read_count()?;
        let mut topology = Topology::new();
        for _ in 0..atom_count {
            let atomic_number = reader.read_u8()?;
            let element = Element::from_atomic_number(atomic_number)
                .ok_or(ModelError::UnknownElement(atomic_number))?;
            let position_m = [reader.read_f64()?, reader.read_f64()?, reader.read_f64()?];
            let charge_c = reader.read_f64()?;
            let atom_type = reader.read_str()?;
            topology.add_atom(Atom::new(element, position_m, charge_c, atom_type));
        }

        let bond_count = reader.read_count()?;
        for _ in 0..bond_count {
            let u = reader.read_u32()?;
            let v = reader.read_u32()?;
            let order = reader.read_u8()?;
            let tag = reader.read_u8()?;
            let bond_type = BondType::from_u8(tag).ok_or(ModelError::InvalidBondType(tag))?;
            topology.add_bond(Bond::new(u, v, order, bond_type))?;
        }

        Ok(Self {
            name,
            topology,
            material,
            metadata,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_part() -> Part {
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
        topology.add_atom(Atom::new(Element::HYDROGEN, [0.0, 0.0, 1.0e-10], 0.2, "H1"));
        topology
            .add_bond(Bond::new(0, 1, 2, BondType::Double))
            .expect("valid bond");
        topology
            .add_bond(Bond::new(0, 2, 1, BondType::Single))
            .expect("valid bond");

        let mut part = Part::new("methanol-fragment", topology).with_material("diamond");
        part.metadata
            .insert("generator".to_owned(), "manual".to_owned());
        part.metadata.insert("units".to_owned(), "si".to_owned());
        part
    }

    #[test]
    fn a_part_round_trips_through_a_byte_buffer() {
        let part = sample_part();
        let bytes = part.to_bytes().expect("encode");
        let restored = Part::from_bytes(&bytes).expect("decode");
        assert_eq!(restored, part);
        assert_eq!(restored.atom_count(), 3);
        assert_eq!(restored.bond_count(), 2);
    }

    #[test]
    fn an_empty_part_round_trips() {
        let part = Part::new("empty", Topology::new());
        let bytes = part.to_bytes().expect("encode");
        assert_eq!(Part::from_bytes(&bytes).expect("decode"), part);
    }

    #[test]
    fn trailing_bytes_are_rejected() {
        let part = sample_part();
        let mut bytes = part.to_bytes().expect("encode");
        bytes.push(0);
        assert_eq!(
            Part::from_bytes(&bytes),
            Err(ModelError::TrailingBytes { remaining: 1 })
        );
    }

    #[test]
    fn every_truncated_prefix_is_an_error_not_a_panic() {
        let part = sample_part();
        let bytes = part.to_bytes().expect("encode");
        for len in 0..bytes.len() {
            assert!(
                Part::from_bytes(&bytes[..len]).is_err(),
                "prefix of length {len} decoded"
            );
        }
    }

    #[test]
    fn an_unknown_element_is_rejected() {
        let mut writer = Writer::new();
        writer.write_str("bad").expect("fits");
        writer.write_str("").expect("fits");
        writer.write_count(0).expect("fits");
        writer.write_count(1).expect("fits");
        writer.write_u8(200);
        let bytes = writer.into_bytes();
        assert_eq!(
            Part::from_bytes(&bytes),
            Err(ModelError::UnknownElement(200))
        );
    }

    #[test]
    fn an_unknown_bond_type_is_rejected() {
        let mut writer = Writer::new();
        writer.write_str("bad").expect("fits");
        writer.write_str("").expect("fits");
        writer.write_count(0).expect("fits");
        writer.write_count(1).expect("fits");
        writer.write_u8(6);
        writer.write_f64(0.0);
        writer.write_f64(0.0);
        writer.write_f64(0.0);
        writer.write_f64(0.0);
        writer.write_str("C3").expect("fits");
        writer.write_count(1).expect("fits");
        writer.write_u32(0);
        writer.write_u32(0);
        writer.write_u8(1);
        writer.write_u8(99);
        let bytes = writer.into_bytes();
        assert_eq!(
            Part::from_bytes(&bytes),
            Err(ModelError::InvalidBondType(99))
        );
    }
}
