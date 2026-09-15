use core::fmt;

/// A chemical element, identified by atomic number.
///
/// The internal representation is one byte, so an element array stays dense.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Element(u8);

impl Element {
    /// Hydrogen.
    pub const HYDROGEN: Self = Self(1);
    /// Helium.
    pub const HELIUM: Self = Self(2);
    /// Lithium.
    pub const LITHIUM: Self = Self(3);
    /// Beryllium.
    pub const BERYLLIUM: Self = Self(4);
    /// Boron.
    pub const BORON: Self = Self(5);
    /// Carbon.
    pub const CARBON: Self = Self(6);
    /// Nitrogen.
    pub const NITROGEN: Self = Self(7);
    /// Oxygen.
    pub const OXYGEN: Self = Self(8);
    /// Fluorine.
    pub const FLUORINE: Self = Self(9);
    /// Neon.
    pub const NEON: Self = Self(10);
    /// Sodium.
    pub const SODIUM: Self = Self(11);
    /// Magnesium.
    pub const MAGNESIUM: Self = Self(12);
    /// Aluminium.
    pub const ALUMINIUM: Self = Self(13);
    /// Silicon.
    pub const SILICON: Self = Self(14);
    /// Phosphorus.
    pub const PHOSPHORUS: Self = Self(15);
    /// Sulfur.
    pub const SULFUR: Self = Self(16);
    /// Chlorine.
    pub const CHLORINE: Self = Self(17);
    /// Argon.
    pub const ARGON: Self = Self(18);
    /// Potassium.
    pub const POTASSIUM: Self = Self(19);
    /// Calcium.
    pub const CALCIUM: Self = Self(20);
    /// Iron.
    pub const IRON: Self = Self(26);
    /// Copper.
    pub const COPPER: Self = Self(29);
    /// Zinc.
    pub const ZINC: Self = Self(30);
    /// Silver.
    pub const SILVER: Self = Self(47);
    /// Gold.
    pub const GOLD: Self = Self(79);
    /// Platinum.
    pub const PLATINUM: Self = Self(78);

    /// The largest supported atomic number, oganesson.
    pub const MAX_ATOMIC_NUMBER: u8 = 118;

    /// Builds an element from an atomic number. Zero is allowed as a dummy.
    pub const fn from_atomic_number(atomic_number: u8) -> Option<Self> {
        if atomic_number <= Self::MAX_ATOMIC_NUMBER {
            Some(Self(atomic_number))
        } else {
            None
        }
    }

    /// Returns the atomic number.
    pub const fn atomic_number(self) -> u8 {
        self.0
    }

    /// Returns the chemical symbol for the elements this build supports.
    pub const fn symbol(self) -> Option<&'static str> {
        match self.0 {
            1 => Some("H"),
            2 => Some("He"),
            3 => Some("Li"),
            4 => Some("Be"),
            5 => Some("B"),
            6 => Some("C"),
            7 => Some("N"),
            8 => Some("O"),
            9 => Some("F"),
            10 => Some("Ne"),
            11 => Some("Na"),
            12 => Some("Mg"),
            13 => Some("Al"),
            14 => Some("Si"),
            15 => Some("P"),
            16 => Some("S"),
            17 => Some("Cl"),
            18 => Some("Ar"),
            19 => Some("K"),
            20 => Some("Ca"),
            21 => Some("Sc"),
            22 => Some("Ti"),
            23 => Some("V"),
            24 => Some("Cr"),
            25 => Some("Mn"),
            26 => Some("Fe"),
            27 => Some("Co"),
            28 => Some("Ni"),
            29 => Some("Cu"),
            30 => Some("Zn"),
            31 => Some("Ga"),
            32 => Some("Ge"),
            33 => Some("As"),
            34 => Some("Se"),
            35 => Some("Br"),
            36 => Some("Kr"),
            47 => Some("Ag"),
            78 => Some("Pt"),
            79 => Some("Au"),
            _ => None,
        }
    }
}

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.symbol() {
            Some(symbol) => f.write_str(symbol),
            None => write!(f, "Z{}", self.0),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_known_element_round_trips_through_its_atomic_number() {
        let carbon = Element::from_atomic_number(6);
        assert_eq!(carbon, Some(Element::CARBON));
        assert_eq!(carbon.map(Element::atomic_number), Some(6));
        assert_eq!(carbon.and_then(Element::symbol), Some("C"));
    }

    #[test]
    fn an_out_of_range_atomic_number_is_rejected() {
        assert_eq!(Element::from_atomic_number(119), None);
        assert_eq!(Element::from_atomic_number(u8::MAX), None);
    }

    #[test]
    fn zero_is_a_supported_dummy_element() {
        let dummy = Element::from_atomic_number(0);
        assert!(dummy.is_some());
        assert_eq!(dummy.and_then(Element::symbol), None);
    }

    #[test]
    fn display_prefers_the_symbol_and_falls_back_to_the_number() {
        assert_eq!(Element::OXYGEN.to_string(), "O");
        let unknown = Element::from_atomic_number(100).expect("known");
        assert_eq!(unknown.to_string(), "Z100");
    }
}
