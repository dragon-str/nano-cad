/// The kind of a chemical bond.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum BondType {
    /// A single bond.
    Single,
    /// A double bond.
    Double,
    /// A triple bond.
    Triple,
    /// An aromatic bond.
    Aromatic,
    /// A bond of unknown kind.
    #[default]
    Unknown,
}

impl BondType {
    /// Returns the stable tag for the on-disk encoding.
    pub const fn to_u8(self) -> u8 {
        match self {
            Self::Single => 0,
            Self::Double => 1,
            Self::Triple => 2,
            Self::Aromatic => 3,
            Self::Unknown => 4,
        }
    }

    /// Reads a bond type from its on-disk tag.
    pub const fn from_u8(tag: u8) -> Option<Self> {
        match tag {
            0 => Some(Self::Single),
            1 => Some(Self::Double),
            2 => Some(Self::Triple),
            3 => Some(Self::Aromatic),
            4 => Some(Self::Unknown),
            _ => None,
        }
    }
}

/// A bond between two atoms, held as integer indices.
///
/// `u` and `v` index the atom arrays of the same topology. They are never
/// pointers.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Bond {
    /// Index of the first atom.
    pub u: u32,
    /// Index of the second atom.
    pub v: u32,
    /// Bond order. Zero means unspecified.
    pub order: u8,
    /// The kind of the bond.
    pub bond_type: BondType,
}

impl Bond {
    /// The largest supported bond order.
    pub const MAX_ORDER: u8 = 4;

    /// Builds a bond.
    pub const fn new(u: u32, v: u32, order: u8, bond_type: BondType) -> Self {
        Self {
            u,
            v,
            order,
            bond_type,
        }
    }

    /// Reports whether a bond order is supported.
    pub const fn is_order_valid(order: u8) -> bool {
        order <= Self::MAX_ORDER
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bond_type_tags_round_trip() {
        for tag in 0..=4u8 {
            let bond_type = BondType::from_u8(tag).expect("known tag");
            assert_eq!(bond_type.to_u8(), tag);
        }
    }

    #[test]
    fn an_unknown_bond_type_tag_is_rejected() {
        assert_eq!(BondType::from_u8(5), None);
    }

    #[test]
    fn bond_orders_above_four_are_invalid() {
        assert!(Bond::is_order_valid(0));
        assert!(Bond::is_order_valid(4));
        assert!(!Bond::is_order_valid(5));
    }
}
