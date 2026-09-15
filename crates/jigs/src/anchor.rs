use nanocad_model::Topology;
use nanocad_units::{Quantity, Unit};

use crate::error::JigError;
use crate::jig::{atom_position_m, validate_positions, Jig, JigKind};

/// An anchor jig: it holds one or more atoms at fixed positions.
///
/// The anchor is a stiff spring between each held atom and its hold point. The
/// potential is
///
/// `U = 0.5 * k_n_per_m * |r_i - r_hold_i|^2`,
///
/// summed over the held atoms. The gradient with respect to atom `i` is
/// `k_n_per_m * (r_i - r_hold_i)`. A coincident atom and hold point give zero
/// energy and zero force.
#[derive(Clone, Debug, PartialEq)]
pub struct AnchorJig {
    atom_count: usize,
    k_n_per_m: f64,
    atom_indices: Vec<u32>,
    hold_positions_m: Vec<[f64; 3]>,
}

impl AnchorJig {
    /// Creates an anchor with no held atoms over a topology.
    ///
    /// The stiffness must be finite and non-negative.
    pub fn new(topology: &Topology, k_n_per_m: f64) -> Result<Self, JigError> {
        if !k_n_per_m.is_finite() || k_n_per_m < 0.0 {
            return Err(JigError::InvalidStiffness { k_n_per_m });
        }
        Ok(Self {
            atom_count: topology.atom_count(),
            k_n_per_m,
            atom_indices: Vec::new(),
            hold_positions_m: Vec::new(),
        })
    }

    /// Holds one atom at a position in SI metres.
    ///
    /// The atom index must name an atom of the topology. The hold position must
    /// be finite.
    pub fn hold_atom(
        &mut self,
        topology: &Topology,
        atom: u32,
        hold_position_m: [f64; 3],
    ) -> Result<(), JigError> {
        let atom_count = topology.atom_count();
        if atom as usize >= atom_count {
            return Err(JigError::AtomIndexOutOfBounds {
                atom: atom as usize,
                atom_count,
            });
        }
        for value in &hold_position_m {
            if !value.is_finite() {
                return Err(JigError::NonFiniteHoldPosition {
                    atom: atom as usize,
                });
            }
        }
        self.atom_count = atom_count;
        self.atom_indices.push(atom);
        self.hold_positions_m.push(hold_position_m);
        Ok(())
    }

    /// Holds one atom at its current position. This reads the topology.
    pub fn hold_atom_current(&mut self, topology: &Topology, atom: u32) -> Result<(), JigError> {
        let position_m =
            topology
                .position_m(atom as usize)
                .ok_or(JigError::AtomIndexOutOfBounds {
                    atom: atom as usize,
                    atom_count: topology.atom_count(),
                })?;
        self.hold_atom(topology, atom, position_m)
    }

    /// Returns the spring stiffness in newtons per metre.
    pub fn k_n_per_m(&self) -> f64 {
        self.k_n_per_m
    }

    /// Returns the number of held atoms.
    pub fn held_count(&self) -> usize {
        self.atom_indices.len()
    }

    /// Returns the index of a held atom.
    pub fn held_atom(&self, index: usize) -> Option<u32> {
        self.atom_indices.get(index).copied()
    }

    /// Returns the hold position of a held atom in SI metres.
    pub fn hold_position_m(&self, index: usize) -> Option<[f64; 3]> {
        self.hold_positions_m.get(index).copied()
    }

    /// Returns the hold position of a held atom as a quantity per axis.
    pub fn hold_position_quantity_m(&self, index: usize) -> Option<[Quantity; 3]> {
        let position_m = self.hold_position_m(index)?;
        Some([
            Quantity::from_si(position_m[0], Unit::Metre),
            Quantity::from_si(position_m[1], Unit::Metre),
            Quantity::from_si(position_m[2], Unit::Metre),
        ])
    }
}

impl Jig for AnchorJig {
    fn kind(&self) -> JigKind {
        JigKind::Anchor
    }

    fn atom_count(&self) -> usize {
        self.atom_count
    }

    fn energy_j(&self, positions_m: &[f64]) -> Result<f64, JigError> {
        validate_positions(positions_m, self.atom_count)?;
        let mut energy_j = 0.0;
        for (slot, atom) in self.atom_indices.iter().enumerate() {
            let position_m = atom_position_m(positions_m, *atom);
            let hold_m = self.hold_positions_m[slot];
            let dx = position_m[0] - hold_m[0];
            let dy = position_m[1] - hold_m[1];
            let dz = position_m[2] - hold_m[2];
            energy_j += 0.5 * self.k_n_per_m * (dx * dx + dy * dy + dz * dz);
        }
        Ok(energy_j)
    }

    fn gradient_j_per_m(&self, positions_m: &[f64]) -> Result<Vec<f64>, JigError> {
        validate_positions(positions_m, self.atom_count)?;
        let mut gradient = vec![0.0; 3 * self.atom_count];
        for (slot, atom) in self.atom_indices.iter().enumerate() {
            let base = *atom as usize * 3;
            let hold_m = self.hold_positions_m[slot];
            gradient[base] = self.k_n_per_m * (positions_m[base] - hold_m[0]);
            gradient[base + 1] = self.k_n_per_m * (positions_m[base + 1] - hold_m[1]);
            gradient[base + 2] = self.k_n_per_m * (positions_m[base + 2] - hold_m[2]);
        }
        Ok(gradient)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nanocad_model::{Atom, Element};

    struct Rng(u64);

    impl Rng {
        fn next_f64(&mut self) -> f64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            (x >> 11) as f64 / (1u64 << 53) as f64
        }

        fn symmetric(&mut self, scale_m: f64) -> f64 {
            (self.next_f64() - 0.5) * 2.0 * scale_m
        }
    }

    fn random_topology(atom_count: usize, seed: u64) -> Topology {
        let mut rng = Rng(seed);
        let mut topology = Topology::new();
        for _ in 0..atom_count {
            let position_m = [
                rng.symmetric(1.5e-10),
                rng.symmetric(1.5e-10),
                rng.symmetric(1.5e-10),
            ];
            topology.add_atom(Atom::new(Element::CARBON, position_m, 0.0, "C3"));
        }
        topology
    }

    #[test]
    fn an_empty_anchor_has_zero_energy_and_an_empty_gradient() {
        let topology = random_topology(3, 0x01);
        let anchor = AnchorJig::new(&topology, 300.0).expect("valid stiffness");
        assert_eq!(anchor.kind(), JigKind::Anchor);
        assert_eq!(anchor.held_count(), 0);
        let positions_m = topology.positions_m().to_vec();
        assert_eq!(anchor.energy_j(&positions_m).expect("valid"), 0.0);
        assert!(anchor
            .gradient_j_per_m(&positions_m)
            .expect("valid")
            .iter()
            .all(|value| *value == 0.0));
    }

    #[test]
    fn an_atom_at_its_hold_point_has_zero_energy_and_zero_force() {
        let topology = random_topology(2, 0x02);
        let mut anchor = AnchorJig::new(&topology, 300.0).expect("valid stiffness");
        let hold_m = topology.position_m(1).expect("atom exists");
        anchor.hold_atom(&topology, 1, hold_m).expect("valid hold");
        let positions_m = topology.positions_m().to_vec();
        assert!(anchor.energy_j(&positions_m).expect("valid").abs() < 1.0e-30);
        for force in anchor.forces_n(&positions_m).expect("valid") {
            assert!(force.abs() < 1.0e-20);
        }
    }

    #[test]
    fn a_displaced_atom_has_the_expected_energy() {
        let topology = random_topology(1, 0x03);
        let mut anchor = AnchorJig::new(&topology, 300.0).expect("valid stiffness");
        let hold_m = [0.0, 0.0, 0.0];
        anchor.hold_atom(&topology, 0, hold_m).expect("valid hold");
        let delta_m = 1.0e-11;
        let positions_m = vec![delta_m, 0.0, 0.0];
        let expected_j = 0.5 * 300.0 * delta_m * delta_m;
        let energy_j = anchor.energy_j(&positions_m).expect("valid");
        assert!((energy_j - expected_j).abs() / expected_j < 1.0e-12);
    }

    #[test]
    fn hold_atom_current_captures_the_topology_position() {
        let topology = random_topology(1, 0x04);
        let mut anchor = AnchorJig::new(&topology, 100.0).expect("valid stiffness");
        anchor.hold_atom_current(&topology, 0).expect("valid hold");
        assert_eq!(anchor.hold_position_m(0), topology.position_m(0));
    }

    #[test]
    fn the_hold_position_is_available_as_a_quantity() {
        let topology = random_topology(1, 0x05);
        let mut anchor = AnchorJig::new(&topology, 100.0).expect("valid stiffness");
        anchor
            .hold_atom(&topology, 0, [1.0e-10, 2.0e-10, 3.0e-10])
            .expect("valid hold");
        let quantity = anchor.hold_position_quantity_m(0).expect("held atom");
        assert_eq!(quantity[0].value_in(Unit::Metre).expect("length"), 1.0e-10);
        assert_eq!(quantity[1].value_in(Unit::Angstrom).expect("length"), 2.0);
    }

    #[test]
    fn invalid_anchor_input_is_rejected() {
        let topology = random_topology(2, 0x06);
        assert!(matches!(
            AnchorJig::new(&topology, -1.0),
            Err(JigError::InvalidStiffness { .. })
        ));
        let mut anchor = AnchorJig::new(&topology, 300.0).expect("valid stiffness");
        assert!(matches!(
            anchor.hold_atom(&topology, 7, [0.0, 0.0, 0.0]),
            Err(JigError::AtomIndexOutOfBounds { .. })
        ));
        assert!(matches!(
            anchor.hold_atom(&topology, 0, [f64::NAN, 0.0, 0.0]),
            Err(JigError::NonFiniteHoldPosition { .. })
        ));
    }

    #[test]
    fn a_position_buffer_of_the_wrong_size_is_rejected() {
        let topology = random_topology(3, 0x07);
        let anchor = AnchorJig::new(&topology, 300.0).expect("valid stiffness");
        assert!(matches!(
            anchor.energy_j(&[0.0, 0.0, 0.0]),
            Err(JigError::PositionBufferSizeMismatch { .. })
        ));
    }

    #[test]
    fn the_analytic_gradient_matches_a_central_finite_difference() {
        let topology = random_topology(5, 0x08);
        let mut anchor = AnchorJig::new(&topology, 300.0).expect("valid stiffness");
        let mut rng = Rng(0x09);
        for atom in 0..topology.atom_count() as u32 {
            let position_m = topology.position_m(atom as usize).expect("atom exists");
            let hold_m = [
                position_m[0] + rng.symmetric(0.02e-10),
                position_m[1] + rng.symmetric(0.02e-10),
                position_m[2] + rng.symmetric(0.02e-10),
            ];
            anchor
                .hold_atom(&topology, atom, hold_m)
                .expect("valid hold");
        }
        let positions_m = topology.positions_m().to_vec();
        let analytic = anchor.gradient_j_per_m(&positions_m).expect("valid");

        let step_m = 1.0e-14;
        let atol_n = 1.0e-15;
        let rtol = 1.0e-7;
        let mut max_error_n = 0.0_f64;
        for index in 0..positions_m.len() {
            let mut forward_m = positions_m.clone();
            forward_m[index] += step_m;
            let mut backward_m = positions_m.clone();
            backward_m[index] -= step_m;
            let forward_j = anchor.energy_j(&forward_m).expect("valid");
            let backward_j = anchor.energy_j(&backward_m).expect("valid");
            let finite_difference_n = (forward_j - backward_j) / (2.0 * step_m);
            let error_n = (analytic[index] - finite_difference_n).abs();
            max_error_n = max_error_n.max(error_n);
            assert!(
                error_n <= atol_n + rtol * finite_difference_n.abs(),
                "coordinate {index}: analytic {} finite difference {} error {}",
                analytic[index],
                finite_difference_n,
                error_n
            );
        }
        assert!(
            max_error_n < 1.0e-15,
            "maximum gradient error {max_error_n} exceeds the absolute tolerance"
        );
    }

    #[test]
    fn the_force_is_the_negative_gradient() {
        let topology = random_topology(4, 0x0a);
        let mut anchor = AnchorJig::new(&topology, 250.0).expect("valid stiffness");
        anchor
            .hold_atom(&topology, 0, [0.0, 0.0, 0.0])
            .expect("valid hold");
        anchor
            .hold_atom(&topology, 3, [1.0e-10, 1.0e-10, 1.0e-10])
            .expect("valid hold");
        let positions_m = topology.positions_m();
        let gradient = anchor.gradient_j_per_m(positions_m).expect("valid");
        let forces = anchor.forces_n(positions_m).expect("valid");
        for (force, grad) in forces.iter().zip(gradient.iter()) {
            assert!((force + grad).abs() < 1.0e-30);
        }
    }
}
