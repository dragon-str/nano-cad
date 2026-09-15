use nanocad_model::Topology;

use crate::error::EngineError;

/// Harmonic bond-stretch parameters for one bond.
///
/// The potential is `U(r) = 0.5 * k_n_per_m * (r - r0_m)^2`, with `r` the
/// current distance in metres. The force constant is in newtons per metre.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BondStretchParams {
    /// Force constant in newtons per metre.
    pub k_n_per_m: f64,
    /// Equilibrium bond length in metres.
    pub r0_m: f64,
}

impl BondStretchParams {
    /// Builds parameters for one bond.
    pub const fn new(k_n_per_m: f64, r0_m: f64) -> Self {
        Self { k_n_per_m, r0_m }
    }
}

/// A harmonic bond-stretch term over a set of bonds.
///
/// The term reads its endpoints from the bond list of a [`Topology`] and holds
/// one parameter set per bond. Energy is in joules. The gradient is the
/// derivative of the energy with respect to each Cartesian coordinate, in
/// newtons.
#[derive(Clone, Debug, PartialEq)]
pub struct BondStretchTerm {
    atom_count: usize,
    endpoints: Vec<[u32; 2]>,
    k_n_per_m: Vec<f64>,
    r0_m: Vec<f64>,
}

impl Default for BondStretchTerm {
    fn default() -> Self {
        Self::new()
    }
}

impl BondStretchTerm {
    /// Creates an empty term with no bonds and no atoms.
    pub const fn new() -> Self {
        Self {
            atom_count: 0,
            endpoints: Vec::new(),
            k_n_per_m: Vec::new(),
            r0_m: Vec::new(),
        }
    }

    /// Builds a term over the bonds of a topology.
    ///
    /// `params` must hold one entry per bond, in bond index order. A mismatch
    /// returns an error.
    pub fn from_topology(
        topology: &Topology,
        params: &[BondStretchParams],
    ) -> Result<Self, EngineError> {
        if params.len() != topology.bond_count() {
            return Err(EngineError::BondParameterCountMismatch {
                provided: params.len(),
                expected: topology.bond_count(),
            });
        }
        let mut term = Self {
            atom_count: topology.atom_count(),
            endpoints: Vec::with_capacity(params.len()),
            k_n_per_m: Vec::with_capacity(params.len()),
            r0_m: Vec::with_capacity(params.len()),
        };
        for (bond, param) in params.iter().enumerate() {
            let u = topology
                .bond_u(bond)
                .ok_or(EngineError::BondEndpointOutOfBounds {
                    bond,
                    atom_count: topology.atom_count(),
                })?;
            let v = topology
                .bond_v(bond)
                .ok_or(EngineError::BondEndpointOutOfBounds {
                    bond,
                    atom_count: topology.atom_count(),
                })?;
            term.add_bond(u, v, param.k_n_per_m, param.r0_m)?;
        }
        Ok(term)
    }

    /// Appends one bond with its parameters.
    ///
    /// The term grows its atom count to cover the larger endpoint. The force
    /// constant must be finite and non-negative. The equilibrium length must
    /// be finite and positive. The endpoints must differ.
    pub fn add_bond(
        &mut self,
        u: u32,
        v: u32,
        k_n_per_m: f64,
        r0_m: f64,
    ) -> Result<(), EngineError> {
        let bond = self.endpoints.len();
        if u == v {
            return Err(EngineError::InvalidBondParameters { bond });
        }
        if !k_n_per_m.is_finite() || k_n_per_m < 0.0 {
            return Err(EngineError::InvalidBondParameters { bond });
        }
        if !r0_m.is_finite() || r0_m <= 0.0 {
            return Err(EngineError::InvalidBondParameters { bond });
        }
        self.endpoints.push([u, v]);
        self.k_n_per_m.push(k_n_per_m);
        self.r0_m.push(r0_m);
        let span = (u.max(v) as usize).saturating_add(1);
        if span > self.atom_count {
            self.atom_count = span;
        }
        Ok(())
    }

    /// Returns the number of bonds.
    pub fn bond_count(&self) -> usize {
        self.endpoints.len()
    }

    /// Returns the number of atoms the term covers.
    pub fn atom_count(&self) -> usize {
        self.atom_count
    }

    /// Returns the total bond energy in joules for a flat position buffer.
    pub fn energy_j(&self, positions_m: &[f64]) -> Result<f64, EngineError> {
        self.validate_positions(positions_m)?;
        let mut energy_j = 0.0;
        for bond in 0..self.bond_count() {
            let [u, v] = self.endpoints[bond];
            let distance_m = self.distance_m(positions_m, u, v);
            let delta_m = distance_m - self.r0_m[bond];
            energy_j += 0.5 * self.k_n_per_m[bond] * delta_m * delta_m;
        }
        Ok(energy_j)
    }

    /// Returns the energy gradient in newtons for a flat position buffer.
    ///
    /// Entry `3 * atom + axis` is `dU/dx`. A bond with coincident atoms has no
    /// defined direction and returns an error.
    pub fn gradient_j_per_m(&self, positions_m: &[f64]) -> Result<Vec<f64>, EngineError> {
        self.validate_positions(positions_m)?;
        let mut gradient = vec![0.0; 3 * self.atom_count];
        for bond in 0..self.bond_count() {
            let [u, v] = self.endpoints[bond];
            let base_u = u as usize * 3;
            let base_v = v as usize * 3;
            let dx = positions_m[base_u] - positions_m[base_v];
            let dy = positions_m[base_u + 1] - positions_m[base_v + 1];
            let dz = positions_m[base_u + 2] - positions_m[base_v + 2];
            let distance_m = (dx * dx + dy * dy + dz * dz).sqrt();
            if distance_m <= 0.0 {
                return Err(EngineError::CoincidentBondAtoms { bond });
            }
            let scale = self.k_n_per_m[bond] * (distance_m - self.r0_m[bond]) / distance_m;
            gradient[base_u] += scale * dx;
            gradient[base_u + 1] += scale * dy;
            gradient[base_u + 2] += scale * dz;
            gradient[base_v] -= scale * dx;
            gradient[base_v + 1] -= scale * dy;
            gradient[base_v + 2] -= scale * dz;
        }
        Ok(gradient)
    }

    /// Returns the force in newtons for a flat position buffer. The force is
    /// the negative of the energy gradient.
    pub fn forces_n(&self, positions_m: &[f64]) -> Result<Vec<f64>, EngineError> {
        let mut forces = self.gradient_j_per_m(positions_m)?;
        for force in &mut forces {
            *force = -*force;
        }
        Ok(forces)
    }

    /// Returns the total energy and the gradient together.
    pub fn energy_and_gradient_j(
        &self,
        positions_m: &[f64],
    ) -> Result<(f64, Vec<f64>), EngineError> {
        let energy_j = self.energy_j(positions_m)?;
        let gradient = self.gradient_j_per_m(positions_m)?;
        Ok((energy_j, gradient))
    }

    fn validate_positions(&self, positions_m: &[f64]) -> Result<(), EngineError> {
        if positions_m.len() != 3 * self.atom_count {
            return Err(EngineError::PositionBufferSizeMismatch {
                len: positions_m.len(),
                atom_count: self.atom_count,
            });
        }
        for (index, value) in positions_m.iter().enumerate() {
            if !value.is_finite() {
                return Err(EngineError::NonFinitePosition { index });
            }
        }
        Ok(())
    }

    fn distance_m(&self, positions_m: &[f64], u: u32, v: u32) -> f64 {
        let base_u = u as usize * 3;
        let base_v = v as usize * 3;
        let dx = positions_m[base_u] - positions_m[base_v];
        let dy = positions_m[base_u + 1] - positions_m[base_v + 1];
        let dz = positions_m[base_u + 2] - positions_m[base_v + 2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nanocad_model::{Atom, Bond, BondType, Element};

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

    fn chain_topology(atom_count: usize, r0_m: f64, seed: u64) -> Topology {
        let mut rng = Rng(seed);
        let mut topology = Topology::new();
        for atom in 0..atom_count {
            let x_m = atom as f64 * r0_m + rng.symmetric(0.02 * r0_m);
            let y_m = rng.symmetric(0.05 * r0_m);
            let z_m = rng.symmetric(0.05 * r0_m);
            topology.add_atom(Atom::new(Element::CARBON, [x_m, y_m, z_m], 0.0, "C3"));
        }
        for atom in 0..atom_count.saturating_sub(1) {
            topology
                .add_bond(Bond::new(atom as u32, atom as u32 + 1, 1, BondType::Single))
                .expect("valid bond");
        }
        topology
    }

    fn chain_term(topology: &Topology, k_n_per_m: f64, r0_m: f64) -> BondStretchTerm {
        let params: Vec<BondStretchParams> = (0..topology.bond_count())
            .map(|_| BondStretchParams::new(k_n_per_m, r0_m))
            .collect();
        BondStretchTerm::from_topology(topology, &params).expect("matching parameters")
    }

    #[test]
    fn an_empty_term_has_zero_energy_and_an_empty_gradient() {
        let term = BondStretchTerm::new();
        assert_eq!(term.bond_count(), 0);
        assert_eq!(term.atom_count(), 0);
        assert_eq!(term.energy_j(&[]).expect("valid"), 0.0);
        assert!(term.gradient_j_per_m(&[]).expect("valid").is_empty());
    }

    #[test]
    fn a_diatomic_at_its_equilibrium_length_has_zero_energy_and_zero_force() {
        let r0_m = 1.5e-10;
        let topology = chain_topology(2, r0_m, 0x1);
        let term = chain_term(&topology, 300.0, r0_m);
        let positions_m = vec![-0.5 * r0_m, 0.0, 0.0, 0.5 * r0_m, 0.0, 0.0];
        assert!(term.energy_j(&positions_m).expect("valid").abs() < 1.0e-30);
        for force in term.forces_n(&positions_m).expect("valid") {
            assert!(force.abs() < 1.0e-20);
        }
    }

    #[test]
    fn a_stretched_bond_has_positive_energy() {
        let r0_m = 1.5e-10;
        let topology = chain_topology(2, r0_m, 0x2);
        let term = chain_term(&topology, 300.0, r0_m);
        let positions_m = vec![0.0, 0.0, 0.0, 2.0 * r0_m, 0.0, 0.0];
        let expected_j = 0.5 * 300.0 * r0_m * r0_m;
        let energy_j = term.energy_j(&positions_m).expect("valid");
        assert!((energy_j - expected_j).abs() / expected_j < 1.0e-12);
    }

    #[test]
    fn a_parameter_count_mismatch_is_rejected() {
        let topology = chain_topology(3, 1.5e-10, 0x3);
        assert!(matches!(
            BondStretchTerm::from_topology(&topology, &[]),
            Err(EngineError::BondParameterCountMismatch { .. })
        ));
    }

    #[test]
    fn invalid_bond_parameters_are_rejected() {
        let mut term = BondStretchTerm::new();
        assert!(matches!(
            term.add_bond(0, 1, -1.0, 1.5e-10),
            Err(EngineError::InvalidBondParameters { .. })
        ));
        assert!(matches!(
            term.add_bond(0, 1, 300.0, 0.0),
            Err(EngineError::InvalidBondParameters { .. })
        ));
        assert!(matches!(
            term.add_bond(0, 0, 300.0, 1.5e-10),
            Err(EngineError::InvalidBondParameters { .. })
        ));
    }

    #[test]
    fn a_position_buffer_of_the_wrong_size_is_rejected() {
        let topology = chain_topology(3, 1.5e-10, 0x4);
        let term = chain_term(&topology, 300.0, 1.5e-10);
        assert!(matches!(
            term.energy_j(&[0.0, 0.0, 0.0]),
            Err(EngineError::PositionBufferSizeMismatch { .. })
        ));
    }

    #[test]
    fn the_analytic_gradient_matches_a_central_finite_difference() {
        let r0_m = 1.5e-10;
        let k_n_per_m = 300.0;
        let topology = chain_topology(6, r0_m, 0x5);
        let term = chain_term(&topology, k_n_per_m, r0_m);
        let positions_m = topology.positions_m().to_vec();
        let analytic = term.gradient_j_per_m(&positions_m).expect("valid");

        let step_m = 1.0e-14;
        let atol_n = 1.0e-15;
        let rtol = 1.0e-7;
        let mut max_error_n = 0.0_f64;
        for index in 0..positions_m.len() {
            let mut forward_m = positions_m.clone();
            forward_m[index] += step_m;
            let mut backward_m = positions_m.clone();
            backward_m[index] -= step_m;
            let forward_j = term.energy_j(&forward_m).expect("valid");
            let backward_j = term.energy_j(&backward_m).expect("valid");
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
        let topology = chain_topology(4, 1.5e-10, 0x6);
        let term = chain_term(&topology, 250.0, 1.4e-10);
        let positions_m = topology.positions_m();
        let gradient = term.gradient_j_per_m(positions_m).expect("valid");
        let forces = term.forces_n(positions_m).expect("valid");
        for (force, grad) in forces.iter().zip(gradient.iter()) {
            assert!((force + grad).abs() < 1.0e-30);
        }
    }
}
