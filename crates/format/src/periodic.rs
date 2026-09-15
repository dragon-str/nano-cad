//! Element symbol lookup.
//!
//! The model crate can build an element from an atomic number but exposes
//! symbols for only a common subset. Text formats such as XYZ name elements by
//! symbol, so this crate carries the full table.

const SYMBOLS: [&str; 119] = [
    "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
    "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
];

/// Returns the symbol for an atomic number, or `None` when out of range.
pub(crate) fn symbol_for(atomic_number: u8) -> Option<&'static str> {
    SYMBOLS
        .get(atomic_number as usize)
        .copied()
        .filter(|symbol| !symbol.is_empty())
}

/// Returns the atomic number for a symbol. The match ignores case.
pub(crate) fn atomic_number_for(symbol: &str) -> Option<u8> {
    let symbol = symbol.trim();
    if symbol.is_empty() {
        return None;
    }
    let first = symbol[..1].to_ascii_uppercase();
    let rest = symbol[1..].to_ascii_lowercase();
    let normalized = format!("{first}{rest}");
    SYMBOLS
        .iter()
        .position(|known| *known == normalized)
        .map(|position| position as u8)
        .filter(|number| *number != 0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn common_symbols_round_trip() {
        for (number, symbol) in [(1u8, "H"), (6, "C"), (14, "Si"), (79, "Au")] {
            assert_eq!(symbol_for(number), Some(symbol));
            assert_eq!(atomic_number_for(symbol), Some(number));
        }
    }

    #[test]
    fn a_lowercase_symbol_is_accepted() {
        assert_eq!(atomic_number_for("si"), Some(14));
        assert_eq!(atomic_number_for("HE"), Some(2));
    }

    #[test]
    fn an_unknown_symbol_or_number_is_rejected() {
        assert_eq!(atomic_number_for("Xx"), None);
        assert_eq!(atomic_number_for(""), None);
        assert_eq!(symbol_for(0), None);
        assert_eq!(symbol_for(119), None);
    }
}
