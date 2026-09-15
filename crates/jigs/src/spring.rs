use nanocad_model::Topology;
use nanocad_units::{Quantity, Unit};

use crate::error::JigError;
use crate::jig::{atom_position_m, distance_m, validate_positions, Jig, JigKind};

/// One endpoint of a spring: a topology atom or a fixed point in space.
#[derive(Clone, Copy, Debug, PartialEq)]
enum SpringEndpoint {
    Atom(u32),
    Fixed([f64; 3]),
}

/// One linear spring with a stiffness and a rest length.
#[derive(Clone, Copy, Debug, PartialEq)]
struct Spring {
    u: SpringEndpoint,
    v: SpringEndpoint,
    k_n_per_m: f64,
    rest_length_m: f64,
}

/// A spring jig: it connects two atoms, or an atom and a fixed point, with a
/// linear spring.
///
/// The potential of one spring is
///
/// `U = 0.5 * k_n_per_m * (d - rest_length_m)^2`,
///
/// with `d` the current endpoint distance in metres. A fixed endpoint carries
/// no gradient. A coincident pair has no defined direction and returns an
/// error.
#[derive(Clone, Debug, PartialEq)]
pub struct SpringJig {
    atom_count: usize,
    springs: Vec<Spring>,
}

impl SpringJig {
    /// Creates a spring jig with no springs over a topology.
    pub fn new(topology: &Topology) -> Self {
        Self {
            atom_count: topology.atom_count(),
            springs: Vec::new(),
        }
    }

    /// Adds a spring between two atoms.
    ///
    /// The endpoints must differ. The stiffness must be finite and
    /// non-negative. The rest length must be finite and positive.
    pub fn add_spring(
        &mut self,
        topology: &Topology,
        u: u32,
        v: u32,
        k_n_per_m: f64,
        rest_length_m: f64,
    ) -> Result<(), JigError> {
        self.check_atom(topology, u)?;
        self.check_atom(topology, v)?;
        let spring = self.springs.len();
        if u == v {
            return Err(JigError::SelfConnectedSpring {
                spring,
                atom: u as usize,
            });
        }
        self.push_spring(
            SpringEndpoint::Atom(u),
            SpringEndpoint::Atom(v),
            k_n_per_m,
            rest_length_m,
        )
    }

    /// Adds a spring between one atom and a fixed point in SI metres.
    pub fn add_spring_to_point(
        &mut self,
        topology: &Topology,
        u: u32,
        point_m: [f64; 3],
        k_n_per_m: f64,
        rest_length_m: f64,
    ) -> Result<(), JigError> {
        self.check_atom(topology, u)?;
        let spring = self.springs.len();
        if point_m.iter().any(|value| !value.is_finite()) {
            return Err(JigError::NonFiniteFixedPoint { spring });
        }
        self.push_spring(
            SpringEndpoint::Atom(u),
            SpringEndpoint::Fixed(point_m),
            k_n_per_m,
            rest_length_m,
        )
    }

    /// Adds a spring between two atoms, with the rest length as a quantity.
    ///
    /// The quantity must measure a length, or the call returns
    /// [`JigError::UnitMismatch`].
    pub fn add_spring_with_rest_length(
        &mut self,
        topology: &Topology,
        u: u32,
        v: u32,
        k_n_per_m: f64,
        rest_length: Quantity,
    ) -> Result<(), JigError> {
        let rest_length_m =
            rest_length
                .value_in(Unit::Metre)
                .map_err(|_| JigError::UnitMismatch {
                    from: rest_length.unit(),
                    to: Unit::Metre,
                })?;
        self.add_spring(topology, u, v, k_n_per_m, rest_length_m)
    }

    fn push_spring(
        &mut self,
        u: SpringEndpoint,
        v: SpringEndpoint,
        k_n_per_m: f64,
        rest_length_m: f64,
    ) -> Result<(), JigError> {
        if !k_n_per_m.is_finite() || k_n_per_m < 0.0 {
            return Err(JigError::InvalidStiffness { k_n_per_m });
        }
        if !rest_length_m.is_finite() || rest_length_m <= 0.0 {
            return Err(JigError::InvalidRestLength { rest_length_m });
        }
        self.springs.push(Spring {
            u,
            v,
            k_n_per_m,
            rest_length_m,
        });
        Ok(())
    }

    fn check_atom(&self, topology: &Topology, atom: u32) -> Result<(), JigError> {
        let atom_count = topology.atom_count();
        if atom as usize >= atom_count {
            return Err(JigError::AtomIndexOutOfBounds {
                atom: atom as usize,
                atom_count,
            });
        }
        Ok(())
    }

    /// Returns the number of springs.
    pub fn spring_count(&self) -> usize {
        self.springs.len()
    }

    /// Returns the stiffness of a spring in newtons per metre.
    pub fn k_n_per_m(&self, spring: usize) -> Option<f64> {
        self.springs.get(spring).map(|spring| spring.k_n_per_m)
    }

    /// Returns the rest length of a spring in SI metres.
    pub fn rest_length_m(&self, spring: usize) -> Option<f64> {
        self.springs.get(spring).map(|spring| spring.rest_length_m)
    }

    /// Returns the rest length of a spring as a quantity.
    pub fn rest_length_quantity_m(&self, spring: usize) -> Option<Quantity> {
        self.rest_length_m(spring)
            .map(|rest_length_m| Quantity::from_si(rest_length_m, Unit::Metre))
    }

    /// Reads one endpoint position from a validated flat buffer.
    fn endpoint_position_m(&self, endpoint: SpringEndpoint, positions_m: &[f64]) -> [f64; 3] {
        match endpoint {
            SpringEndpoint::Atom(atom) => atom_position_m(positions_m, atom),
            SpringEndpoint::Fixed(point_m) => point_m,
        }
    }
}

impl Jig for SpringJig {
    fn kind(&self) -> JigKind {
        JigKind::Spring
    }

    fn atom_count(&self) -> usize {
        self.atom_count
    }

    fn energy_j(&self, positions_m: &[f64]) -> Result<f64, JigError> {
        validate_positions(positions_m, self.atom_count)?;
        let mut energy_j = 0.0;
        for spring in &self.springs {
            let u_m = self.endpoint_position_m(spring.u, positions_m);
            let v_m = self.endpoint_position_m(spring.v, positions_m);
            let delta_m = distance_m(u_m, v_m) - spring.rest_length_m;
            energy_j += 0.5 * spring.k_n_per_m * delta_m * delta_m;
        }
        Ok(energy_j)
    }

    fn gradient_j_per_m(&self, positions_m: &[f64]) -> Result<Vec<f64>, JigError> {
        validate_positions(positions_m, self.atom_count)?;
        let mut gradient = vec![0.0; 3 * self.atom_count];
        for (index, spring) in self.springs.iter().enumerate() {
            let u_m = self.endpoint_position_m(spring.u, positions_m);
            let v_m = self.endpoint_position_m(spring.v, positions_m);
            let dx = u_m[0] - v_m[0];
            let dy = u_m[1] - v_m[1];
            let dz = u_m[2] - v_m[2];
            let distance = (dx * dx + dy * dy + dz * dz).sqrt();
            if distance <= 0.0 {
                return Err(JigError::CoincidentSpringAtoms { spring: index });
            }
            let scale = spring.k_n_per_m * (distance - spring.rest_length_m) / distance;
            if let SpringEndpoint::Atom(atom) = spring.u {
                let base = atom as usize * 3;
                gradient[base] += scale * dx;
                gradient[base + 1] += scale * dy;
                gradient[base + 2] += scale * dz;
            }
            if let SpringEndpoint::Atom(atom) = spring.v {
                let base = atom as usize * 3;
                gradient[base] -= scale * dx;
                gradient[base + 1] -= scale * dy;
                gradient[base + 2] -= scale * dz;
            }
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
    fn an_empty_spring_jig_has_zero_energy_and_an_empty_gradient() {
        let topology = random_topology(3, 0x11);
        let jig = SpringJig::new(&topology);
        assert_eq!(jig.kind(), JigKind::Spring);
        assert_eq!(jig.spring_count(), 0);
        let positions_m = topology.positions_m().to_vec();
        assert_eq!(jig.energy_j(&positions_m).expect("valid"), 0.0);
        assert!(jig
            .gradient_j_per_m(&positions_m)
            .expect("valid")
            .iter()
            .all(|value| *value == 0.0));
    }

    #[test]
    fn a_spring_at_its_rest_length_has_zero_energy() {
        let mut topology = Topology::new();
        topology.add_atom(Atom::new(Element::CARBON, [0.0, 0.0, 0.0], 0.0, "C3"));
        topology.add_atom(Atom::new(Element::CARBON, [1.5e-10, 0.0, 0.0], 0.0, "C3"));
        let mut jig = SpringJig::new(&topology);
        jig.add_spring(&topology, 0, 1, 300.0, 1.5e-10)
            .expect("valid spring");
        let positions_m = topology.positions_m().to_vec();
        assert!(jig.energy_j(&positions_m).expect("valid").abs() < 1.0e-30);
    }

    #[test]
    fn a_stretched_spring_has_the_expected_energy() {
        let mut topology = Topology::new();
        topology.add_atom(Atom::new(Element::CARBON, [0.0, 0.0, 0.0], 0.0, "C3"));
        topology.add_atom(Atom::new(Element::CARBON, [2.0e-10, 0.0, 0.0], 0.0, "C3"));
        let mut jig = SpringJig::new(&topology);
        let rest_length_m = 1.5e-10;
        jig.add_spring(&topology, 0, 1, 300.0, rest_length_m)
            .expect("valid spring");
        let positions_m = topology.positions_m().to_vec();
        let expected_j = 0.5 * 300.0 * (0.5e-10_f64).powi(2);
        let energy_j = jig.energy_j(&positions_m).expect("valid");
        assert!((energy_j - expected_j).abs() / expected_j < 1.0e-12);
    }

    #[test]
    fn an_atom_to_point_spring_has_the_expected_energy() {
        let mut topology = Topology::new();
        topology.add_atom(Atom::new(Element::CARBON, [0.0, 0.0, 0.0], 0.0, "C3"));
        let mut jig = SpringJig::new(&topology);
        let point_m = [1.0e-10, 0.0, 0.0];
        let rest_length_m = 0.5e-10;
        jig.add_spring_to_point(&topology, 0, point_m, 200.0, rest_length_m)
            .expect("valid spring");
        let positions_m = vec![0.0, 0.0, 0.0];
        let expected_j = 0.5 * 200.0 * (0.5e-10_f64).powi(2);
        let energy_j = jig.energy_j(&positions_m).expect("valid");
        assert!((energy_j - expected_j).abs() / expected_j < 1.0e-12);
    }

    #[test]
    fn a_rest_length_quantity_is_accepted() {
        let mut topology = Topology::new();
        topology.add_atom(Atom::new(Element::CARBON, [0.0, 0.0, 0.0], 0.0, "C3"));
        topology.add_atom(Atom::new(Element::CARBON, [1.5e-10, 0.0, 0.0], 0.0, "C3"));
        let mut jig = SpringJig::new(&topology);
        jig.add_spring_with_rest_length(&topology, 0, 1, 300.0, Quantity::new(1.5, Unit::Angstrom))
            .expect("length quantity");
        assert_eq!(jig.rest_length_m(0), Some(1.5e-10));
        assert_eq!(
            jig.rest_length_quantity_m(0)
                .expect("spring exists")
                .value_in(Unit::Angstrom)
                .expect("length"),
            1.5
        );
    }

    #[test]
    fn a_rest_length_quantity_of_the_wrong_dimension_is_rejected() {
        let mut topology = Topology::new();
        topology.add_atom(Atom::new(Element::CARBON, [0.0, 0.0, 0.0], 0.0, "C3"));
        topology.add_atom(Atom::new(Element::CARBON, [1.5e-10, 0.0, 0.0], 0.0, "C3"));
        let mut jig = SpringJig::new(&topology);
        assert!(matches!(
            jig.add_spring_with_rest_length(
                &topology,
                0,
                1,
                300.0,
                Quantity::new(1.0, Unit::Second)
            ),
            Err(JigError::UnitMismatch { .. })
        ));
    }

    #[test]
    fn invalid_spring_input_is_rejected() {
        let mut topology = Topology::new();
        topology.add_atom(Atom::new(Element::CARBON, [0.0, 0.0, 0.0], 0.0, "C3"));
        topology.add_atom(Atom::new(Element::CARBON, [1.5e-10, 0.0, 0.0], 0.0, "C3"));
        let mut jig = SpringJig::new(&topology);
        assert!(matches!(
            jig.add_spring(&topology, 0, 0, 300.0, 1.5e-10),
            Err(JigError::SelfConnectedSpring { .. })
        ));
        assert!(matches!(
            jig.add_spring(&topology, 0, 1, -1.0, 1.5e-10),
            Err(JigError::InvalidStiffness { .. })
        ));
        assert!(matches!(
            jig.add_spring(&topology, 0, 1, 300.0, 0.0),
            Err(JigError::InvalidRestLength { .. })
        ));
        assert!(matches!(
            jig.add_spring(&topology, 0, 7, 300.0, 1.5e-10),
            Err(JigError::AtomIndexOutOfBounds { .. })
        ));
        assert!(matches!(
            jig.add_spring_to_point(&topology, 0, [f64::NAN, 0.0, 0.0], 300.0, 1.5e-10),
            Err(JigError::NonFiniteFixedPoint { .. })
        ));
    }

    #[test]
    fn coincident_spring_atoms_are_rejected_during_gradient() {
        let mut topology = Topology::new();
        topology.add_atom(Atom::new(Element::CARBON, [0.0, 0.0, 0.0], 0.0, "C3"));
        topology.add_atom(Atom::new(Element::CARBON, [0.0, 0.0, 0.0], 0.0, "C3"));
        let mut jig = SpringJig::new(&topology);
        jig.add_spring(&topology, 0, 1, 300.0, 1.5e-10)
            .expect("valid spring");
        let positions_m = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        assert!(matches!(
            jig.gradient_j_per_m(&positions_m),
            Err(JigError::CoincidentSpringAtoms { .. })
        ));
    }

    #[test]
    fn the_analytic_gradient_matches_a_central_finite_difference() {
        let topology = random_topology(6, 0x12);
        let mut jig = SpringJig::new(&topology);
        jig.add_spring(&topology, 0, 1, 300.0, 1.5e-10)
            .expect("valid spring");
        jig.add_spring(&topology, 1, 2, 250.0, 1.4e-10)
            .expect("valid spring");
        jig.add_spring(&topology, 3, 4, 200.0, 1.6e-10)
            .expect("valid spring");
        jig.add_spring_to_point(&topology, 5, [0.0, 0.0, 0.0], 180.0, 1.3e-10)
            .expect("valid spring");
        let positions_m = topology.positions_m().to_vec();
        let analytic = jig.gradient_j_per_m(&positions_m).expect("valid");

        let step_m = 1.0e-14;
        let atol_n = 1.0e-15;
        let rtol = 1.0e-7;
        let mut max_error_n = 0.0_f64;
        for index in 0..positions_m.len() {
            let mut forward_m = positions_m.clone();
            forward_m[index] += step_m;
            let mut backward_m = positions_m.clone();
            backward_m[index] -= step_m;
            let forward_j = jig.energy_j(&forward_m).expect("valid");
            let backward_j = jig.energy_j(&backward_m).expect("valid");
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
        let topology = random_topology(3, 0x13);
        let mut jig = SpringJig::new(&topology);
        jig.add_spring(&topology, 0, 1, 250.0, 1.4e-10)
            .expect("valid spring");
        jig.add_spring_to_point(&topology, 2, [0.0, 0.0, 0.0], 190.0, 1.2e-10)
            .expect("valid spring");
        let positions_m = topology.positions_m();
        let gradient = jig.gradient_j_per_m(positions_m).expect("valid");
        let forces = jig.forces_n(positions_m).expect("valid");
        for (force, grad) in forces.iter().zip(gradient.iter()) {
            assert!((force + grad).abs() < 1.0e-30);
        }
    }
}
