use crate::element::Element;

/// One atom as a value: element, position, charge, and type.
///
/// `Topology` does not store `Atom` values. It stores the same fields in
/// struct-of-arrays form. This type carries one atom across an accessor.
#[derive(Clone, Debug, PartialEq)]
pub struct Atom {
    /// The chemical element.
    pub element: Element,
    /// Position in SI metres.
    pub position_m: [f64; 3],
    /// Partial charge in SI coulombs.
    pub charge_c: f64,
    /// The force-field atom type name.
    pub atom_type: String,
}

impl Atom {
    /// Builds an atom.
    pub fn new(
        element: Element,
        position_m: [f64; 3],
        charge_c: f64,
        atom_type: impl Into<String>,
    ) -> Self {
        Self {
            element,
            position_m,
            charge_c,
            atom_type: atom_type.into(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_atom_carries_every_field() {
        let atom = Atom::new(Element::CARBON, [1.0e-10, 2.0e-10, 3.0e-10], -0.2, "C3");
        assert_eq!(atom.element, Element::CARBON);
        assert_eq!(atom.position_m, [1.0e-10, 2.0e-10, 3.0e-10]);
        assert_eq!(atom.charge_c, -0.2);
        assert_eq!(atom.atom_type, "C3");
    }
}
