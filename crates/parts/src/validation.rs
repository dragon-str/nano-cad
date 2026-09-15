use nanocad_model::{Element, Part};

/// The default bond-strain tolerance. A bond within this relative deviation
/// from the reference length passes.
pub const DEFAULT_STRAIN_TOLERANCE_RELATIVE: f64 = 0.05;

/// One validation failure. The report collects these; it does not stop at the
/// first one.
#[derive(Clone, Debug, PartialEq)]
pub enum Violation {
    /// An atom has more bonds than its element valence allows.
    OverCoordinated {
        /// The atom index in the part topology.
        atom_index: usize,
        /// The element of the atom.
        element: Element,
        /// The element valence.
        valence: u32,
        /// The measured bond count.
        bond_count: u32,
    },
    /// A bond length is outside the strain tolerance of the reference length.
    BondStrain {
        /// The bond index in the part topology.
        bond_index: usize,
        /// The measured bond length, in metres.
        length_m: f64,
        /// The reference bond length, in metres.
        reference_m: f64,
        /// The relative deviation `|length - reference| / reference`.
        deviation_relative: f64,
    },
    /// The element has no tabulated valence, so over-coordination is not
    /// checked for this atom.
    UnknownValence {
        /// The atom index in the part topology.
        atom_index: usize,
        /// The element of the atom.
        element: Element,
    },
}

/// The result of a validation pass. An empty list means the part passed.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ValidationReport {
    /// Every violation found, in atom order then bond order.
    pub violations: Vec<Violation>,
}

impl ValidationReport {
    /// Creates an empty report.
    pub fn new() -> Self {
        Self::default()
    }

    /// Reports whether the part passed with no violations.
    pub fn is_ok(&self) -> bool {
        self.violations.is_empty()
    }

    /// The number of violations.
    pub fn len(&self) -> usize {
        self.violations.len()
    }

    /// Reports whether the report has no violations.
    pub fn is_empty(&self) -> bool {
        self.violations.is_empty()
    }

    /// Returns the violations as a slice.
    pub fn violations(&self) -> &[Violation] {
        &self.violations
    }

    /// Appends one violation.
    pub fn push(&mut self, violation: Violation) {
        self.violations.push(violation);
    }
}

/// Returns the valence of an element, or `None` when it is not tabulated.
///
/// The value is the maximum number of single bonds the element forms. Carbon
/// is 4, oxygen is 2, and so on. A caller with an unusual element can extend
/// this table.
pub fn element_valence(element: Element) -> Option<u32> {
    match element.atomic_number() {
        1 => Some(1),
        3 => Some(1),
        4 => Some(2),
        5 => Some(3),
        6 => Some(4),
        7 => Some(3),
        8 => Some(2),
        9 => Some(1),
        11 => Some(1),
        12 => Some(2),
        13 => Some(3),
        14 => Some(4),
        15 => Some(5),
        16 => Some(6),
        17 => Some(1),
        19 => Some(1),
        20 => Some(2),
        26 => Some(3),
        29 => Some(2),
        30 => Some(2),
        47 => Some(1),
        78 => Some(4),
        79 => Some(3),
        _ => None,
    }
}

/// Validates a part against its element valences and a reference bond length.
///
/// The check returns every violation in one pass:
/// - an atom with more bonds than its valence,
/// - a bond outside `strain_tolerance_relative` of `reference_m`,
/// - an atom with an untabulated element.
pub fn validate_part(
    part: &Part,
    reference_m: f64,
    strain_tolerance_relative: f64,
) -> ValidationReport {
    let mut report = ValidationReport::new();

    let mut bond_counts = vec![0u32; part.atom_count()];
    for bond in part.topology.bonds() {
        bond_counts[bond.u as usize] += 1;
        bond_counts[bond.v as usize] += 1;
    }
    for (index, count) in bond_counts.iter().enumerate() {
        let Some(element) = part.topology.element(index) else {
            continue;
        };
        match element_valence(element) {
            Some(valence) if *count > valence => report.push(Violation::OverCoordinated {
                atom_index: index,
                element,
                valence,
                bond_count: *count,
            }),
            None => report.push(Violation::UnknownValence {
                atom_index: index,
                element,
            }),
            _ => {}
        }
    }

    for (index, bond) in part.topology.bonds().enumerate() {
        let Some(u_m) = part.topology.position_m(bond.u as usize) else {
            continue;
        };
        let Some(v_m) = part.topology.position_m(bond.v as usize) else {
            continue;
        };
        let length_m = crate::geometry::distance_m(u_m, v_m);
        if reference_m <= 0.0 {
            continue;
        }
        let deviation_relative = (length_m - reference_m).abs() / reference_m;
        if deviation_relative > strain_tolerance_relative {
            report.push(Violation::BondStrain {
                bond_index: index,
                length_m,
                reference_m,
                deviation_relative,
            });
        }
    }

    report
}

#[cfg(test)]
mod tests {
    use super::*;
    use nanocad_model::{Atom, Bond, BondType, Topology};

    const CARBON_BOND_M: f64 = 1.544e-10;

    fn carbon_atom(x_m: f64) -> Atom {
        Atom::new(Element::CARBON, [x_m, 0.0, 0.0], 0.0, "C")
    }

    fn krypton_element() -> Element {
        Element::from_atomic_number(36).expect("krypton is supported")
    }

    fn star_atom(index: u32, bond_m: f64) -> Atom {
        let angle_rad = std::f64::consts::TAU * index as f64 / 5.0;
        Atom::new(
            Element::CARBON,
            [bond_m * angle_rad.cos(), bond_m * angle_rad.sin(), 0.0],
            0.0,
            "C",
        )
    }

    fn good_part() -> Part {
        let mut topology = Topology::new();
        topology.add_atom(carbon_atom(0.0));
        topology.add_atom(carbon_atom(CARBON_BOND_M));
        topology
            .add_bond(Bond::new(0, 1, 1, BondType::Single))
            .expect("valid bond");
        Part::new("good", topology).with_material("diamondoid")
    }

    #[test]
    fn a_good_part_has_an_empty_report() {
        let report = validate_part(&good_part(), CARBON_BOND_M, 1.0e-3);
        assert!(report.is_ok());
        assert_eq!(report.len(), 0);
        assert!(report.violations().is_empty());
    }

    #[test]
    fn an_over_coordinated_atom_is_reported() {
        let mut topology = Topology::new();
        topology.add_atom(carbon_atom(0.0));
        for index in 1..=5u32 {
            topology.add_atom(star_atom(index, CARBON_BOND_M));
        }
        for index in 1..=5u32 {
            topology
                .add_bond(Bond::new(0, index, 1, BondType::Single))
                .expect("valid bond");
        }
        let part = Part::new("over", topology);
        let report = validate_part(&part, CARBON_BOND_M, 10.0);
        assert!(!report.is_ok());
        let over: Vec<&Violation> = report
            .violations()
            .iter()
            .filter(|violation| matches!(violation, Violation::OverCoordinated { .. }))
            .collect();
        assert_eq!(over.len(), 1);
        assert_eq!(
            over[0],
            &Violation::OverCoordinated {
                atom_index: 0,
                element: Element::CARBON,
                valence: 4,
                bond_count: 5,
            }
        );
    }

    #[test]
    fn a_stretched_bond_is_reported() {
        let mut topology = Topology::new();
        topology.add_atom(carbon_atom(0.0));
        topology.add_atom(carbon_atom(CARBON_BOND_M * 1.2));
        topology
            .add_bond(Bond::new(0, 1, 1, BondType::Single))
            .expect("valid bond");
        let part = Part::new("stretched", topology);
        let report = validate_part(&part, CARBON_BOND_M, 0.01);
        assert!(!report.is_ok());
        let strain: Vec<&Violation> = report
            .violations()
            .iter()
            .filter(|violation| matches!(violation, Violation::BondStrain { .. }))
            .collect();
        assert_eq!(strain.len(), 1);
        match strain[0] {
            Violation::BondStrain {
                bond_index,
                deviation_relative,
                ..
            } => {
                assert_eq!(*bond_index, 0);
                assert!((deviation_relative - 0.2).abs() < 1.0e-12);
            }
            other => panic!("unexpected violation {other:?}"),
        }
    }

    #[test]
    fn a_part_with_two_faults_reports_both() {
        let mut topology = Topology::new();
        topology.add_atom(carbon_atom(0.0));
        for index in 1..=4u32 {
            topology.add_atom(star_atom(index, CARBON_BOND_M));
        }
        topology.add_atom(carbon_atom(CARBON_BOND_M * 2.0));
        for index in 1..=5u32 {
            topology
                .add_bond(Bond::new(0, index, 1, BondType::Single))
                .expect("valid bond");
        }
        let part = Part::new("two-faults", topology);
        let report = validate_part(&part, CARBON_BOND_M, 0.01);
        assert!(!report.is_ok());
        assert!(report
            .violations()
            .iter()
            .any(|violation| matches!(violation, Violation::OverCoordinated { .. })));
        assert!(report
            .violations()
            .iter()
            .any(|violation| matches!(violation, Violation::BondStrain { .. })));
    }

    #[test]
    fn an_untabulated_element_is_reported() {
        let krypton = krypton_element();
        let mut topology = Topology::new();
        topology.add_atom(Atom::new(krypton, [0.0, 0.0, 0.0], 0.0, "Kr"));
        let part = Part::new("krypton", topology);
        let report = validate_part(&part, CARBON_BOND_M, 0.01);
        assert_eq!(
            report.violations(),
            &[Violation::UnknownValence {
                atom_index: 0,
                element: krypton,
            }]
        );
    }

    #[test]
    fn the_valence_table_covers_the_common_elements() {
        assert_eq!(element_valence(Element::HYDROGEN), Some(1));
        assert_eq!(element_valence(Element::CARBON), Some(4));
        assert_eq!(element_valence(Element::NITROGEN), Some(3));
        assert_eq!(element_valence(Element::OXYGEN), Some(2));
        assert_eq!(element_valence(Element::SILICON), Some(4));
        assert_eq!(element_valence(krypton_element()), None);
    }

    #[test]
    fn a_generated_diamond_part_passes() {
        use crate::generator::PartGenerator;
        use crate::lattice::DiamondGenerator;
        use crate::parameter::ParameterSet;
        let part = DiamondGenerator
            .generate(
                &ParameterSet::new()
                    .with("cells_x", 2.0)
                    .with("cells_y", 2.0)
                    .with("cells_z", 2.0),
            )
            .expect("generate");
        let reference_m = 3.567e-10 * 3.0_f64.sqrt() / 4.0;
        let report = validate_part(&part, reference_m, 1.0e-3);
        assert!(report.is_ok(), "violations: {:?}", report.violations());
    }
}
