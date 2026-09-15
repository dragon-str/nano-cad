use nanocad_model::Topology;

use crate::error::EngineError;
use crate::geometry::{atom_position_m, dot, sub, validate_positions};
use crate::nonbonded::{all_pairs, Cutoff, PeriodicBox};

/// The Coulomb constant in newton square metres per coulomb squared.
///
/// This is `1 / (4 pi epsilon_0)` from the 2018 CODATA recommended value.
pub const COULOMB_CONSTANT_N_M2_PER_C2: f64 = 8.987_551_792_3e9;

/// An electrostatic term over all atom pairs, with a cutoff.
///
/// The potential is `U(r) = k_e q_i q_j / r`. It is multiplied by the switching
/// function of [`Cutoff`] and reduced by the minimum-image convention of
/// [`PeriodicBox`]. Energy is in joules. The gradient is in newtons.
#[derive(Clone, Debug, PartialEq)]
pub struct ElectrostaticTerm {
    atom_count: usize,
    charges_c: Vec<f64>,
    cutoff: Cutoff,
    periodic_box: PeriodicBox,
}

impl ElectrostaticTerm {
    /// Creates an empty term with no atoms.
    pub fn new(cutoff: Cutoff, periodic_box: PeriodicBox) -> Result<Self, EngineError> {
        cutoff.validate_box(&periodic_box)?;
        Ok(Self {
            atom_count: 0,
            charges_c: Vec::new(),
            cutoff,
            periodic_box,
        })
    }

    /// Builds a term from one charge per atom, in SI coulombs.
    ///
    /// Every charge must be finite.
    pub fn from_charges(
        charges_c: &[f64],
        cutoff: Cutoff,
        periodic_box: PeriodicBox,
    ) -> Result<Self, EngineError> {
        let mut term = Self::new(cutoff, periodic_box)?;
        for (atom, charge_c) in charges_c.iter().enumerate() {
            if !charge_c.is_finite() {
                return Err(EngineError::InvalidCharge { atom });
            }
            term.charges_c.push(*charge_c);
            term.atom_count += 1;
        }
        Ok(term)
    }

    /// Builds a term from the charges of a topology.
    pub fn from_topology(
        topology: &Topology,
        cutoff: Cutoff,
        periodic_box: PeriodicBox,
    ) -> Result<Self, EngineError> {
        Self::from_charges(topology.charges_c(), cutoff, periodic_box)
    }

    /// Returns the number of atoms.
    pub fn atom_count(&self) -> usize {
        self.atom_count
    }

    /// Returns the charge at an atom index, in SI coulombs.
    pub fn charge_c(&self, atom: usize) -> Option<f64> {
        self.charges_c.get(atom).copied()
    }

    /// Returns the cutoff.
    pub fn cutoff(&self) -> &Cutoff {
        &self.cutoff
    }

    /// Returns the per-atom charges, in SI coulombs.
    pub(crate) fn charges_c(&self) -> &[f64] {
        &self.charges_c
    }

    /// Returns the periodic box.
    pub fn periodic_box(&self) -> &PeriodicBox {
        &self.periodic_box
    }

    /// Replaces the periodic box. The cutoff must not exceed half of any
    /// periodic length.
    pub fn set_periodic_box(&mut self, periodic_box: PeriodicBox) -> Result<(), EngineError> {
        self.cutoff.validate_box(&periodic_box)?;
        self.periodic_box = periodic_box;
        Ok(())
    }

    /// Returns the total electrostatic energy in joules.
    pub fn energy_j(&self, positions_m: &[f64]) -> Result<f64, EngineError> {
        Ok(self.energy_and_gradient_j(positions_m)?.0)
    }

    /// Returns the energy gradient in newtons.
    pub fn gradient_j_per_m(&self, positions_m: &[f64]) -> Result<Vec<f64>, EngineError> {
        Ok(self.energy_and_gradient_j(positions_m)?.1)
    }

    /// Returns the force in newtons.
    pub fn forces_n(&self, positions_m: &[f64]) -> Result<Vec<f64>, EngineError> {
        let mut forces = self.gradient_j_per_m(positions_m)?;
        for force in &mut forces {
            *force = -*force;
        }
        Ok(forces)
    }

    /// Returns the total energy and the gradient together over every atom pair.
    ///
    /// This is the reference path. It scans all unordered pairs. The
    /// pair-list path of
    /// [`ElectrostaticTerm::energy_and_gradient_from_pairs_j`] gives the same
    /// result for a pair list that holds every pair within the cutoff.
    pub fn energy_and_gradient_j(
        &self,
        positions_m: &[f64],
    ) -> Result<(f64, Vec<f64>), EngineError> {
        let pairs = all_pairs(self.atom_count);
        self.energy_and_gradient_from_pairs_j(positions_m, &pairs)
    }

    /// Returns the energy and the gradient over an explicit pair list.
    ///
    /// `pairs` is flat. Each consecutive pair of entries, at offsets `2k` and
    /// `2k + 1`, is one pair. A pair list from a [`VerletList`] is accepted
    /// directly. Pairs beyond the cutoff contribute nothing. A pair list of odd
    /// length, or an index outside the atom count, returns an error.
    pub fn energy_and_gradient_from_pairs_j(
        &self,
        positions_m: &[f64],
        pairs: &[u32],
    ) -> Result<(f64, Vec<f64>), EngineError> {
        validate_positions(positions_m, self.atom_count)?;
        if pairs.len() & 1 == 1 {
            return Err(EngineError::PairListSizeMismatch { len: pairs.len() });
        }
        let mut energy_j = 0.0;
        let mut gradient = vec![0.0; 3 * self.atom_count];
        let mut entry = 0;
        while entry + 1 < pairs.len() {
            let i = pairs[entry];
            let j = pairs[entry + 1];
            self.check_pair(entry, i, j)?;
            self.accumulate_pair(positions_m, i, j, &mut energy_j, &mut gradient)?;
            entry += 2;
        }
        Ok((energy_j, gradient))
    }

    fn check_pair(&self, entry: usize, i: u32, j: u32) -> Result<(), EngineError> {
        let atom_count = self.atom_count;
        if i as usize >= atom_count {
            return Err(EngineError::PairIndexOutOfBounds {
                entry,
                index: i,
                atom_count,
            });
        }
        if j as usize >= atom_count {
            return Err(EngineError::PairIndexOutOfBounds {
                entry,
                index: j,
                atom_count,
            });
        }
        Ok(())
    }

    fn accumulate_pair(
        &self,
        positions_m: &[f64],
        i: u32,
        j: u32,
        energy_j: &mut f64,
        gradient: &mut [f64],
    ) -> Result<(), EngineError> {
        let Some((value_j, pair_gradient)) = self.pair_energy_and_gradient_j(positions_m, i, j)?
        else {
            return Ok(());
        };
        *energy_j += value_j;
        let base_i = i as usize * 3;
        let base_j = j as usize * 3;
        for axis in 0..3 {
            gradient[base_i + axis] += pair_gradient[axis];
            gradient[base_j + axis] += pair_gradient[3 + axis];
        }
        Ok(())
    }

    /// Returns the energy and the gradient of one unordered pair.
    ///
    /// The six gradient entries are `dU/dx_i` then `dU/dx_j`, in newtons. A
    /// pair at or beyond the cutoff returns `None`.
    pub fn pair_energy_and_gradient_j(
        &self,
        positions_m: &[f64],
        i: u32,
        j: u32,
    ) -> Result<Option<(f64, [f64; 6])>, EngineError> {
        let delta_m = self.periodic_box.minimum_image(sub(
            atom_position_m(positions_m, i),
            atom_position_m(positions_m, j),
        ));
        let r_sq_m2 = dot(delta_m, delta_m);
        let cutoff_sq_m2 = self.cutoff.cutoff_m() * self.cutoff.cutoff_m();
        if r_sq_m2 >= cutoff_sq_m2 {
            return Ok(None);
        }
        let r_m = r_sq_m2.sqrt();
        if r_m <= 0.0 {
            return Err(EngineError::CoincidentNonbondedAtoms { i, j });
        }
        let coulomb_j =
            COULOMB_CONSTANT_N_M2_PER_C2 * self.charges_c[i as usize] * self.charges_c[j as usize];
        let value_j = coulomb_j / r_m;
        let value_prime_n = -coulomb_j / (r_m * r_m);
        let (switch_value, switch_derivative) = self.cutoff.switch_value(r_m);
        let du_dr_n = switch_value * value_prime_n + switch_derivative * value_j;
        let factor_n_per_m = du_dr_n / r_m;
        let mut pair_gradient = [0.0; 6];
        for axis in 0..3 {
            let contribution_n = factor_n_per_m * delta_m[axis];
            pair_gradient[axis] = contribution_n;
            pair_gradient[3 + axis] = -contribution_n;
        }
        Ok(Some((switch_value * value_j, pair_gradient)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::central_difference;
    use nanocad_model::{Atom, Element};

    fn cutoff() -> Cutoff {
        Cutoff::new(1.0e-9, 0.8e-9).expect("valid cutoff")
    }

    fn charges_c() -> Vec<f64> {
        vec![1.0e-19, -1.2e-19, 0.8e-19, -0.6e-19]
    }

    #[test]
    fn an_empty_term_has_zero_energy_and_an_empty_gradient() {
        let term = ElectrostaticTerm::new(cutoff(), PeriodicBox::non_periodic()).expect("valid");
        assert_eq!(term.atom_count(), 0);
        assert_eq!(term.energy_j(&[]).expect("valid"), 0.0);
        assert!(term.gradient_j_per_m(&[]).expect("valid").is_empty());
    }

    #[test]
    fn a_non_finite_charge_is_rejected() {
        assert!(matches!(
            ElectrostaticTerm::from_charges(&[f64::NAN], cutoff(), PeriodicBox::non_periodic()),
            Err(EngineError::InvalidCharge { .. })
        ));
    }

    #[test]
    fn a_topology_supplies_its_charges() {
        let mut topology = Topology::new();
        topology.add_atom(Atom::new(Element::HYDROGEN, [0.0, 0.0, 0.0], 0.5e-19, "H"));
        topology.add_atom(Atom::new(
            Element::HYDROGEN,
            [1.0e-10, 0.0, 0.0],
            -0.5e-19,
            "H",
        ));
        let term =
            ElectrostaticTerm::from_topology(&topology, cutoff(), PeriodicBox::non_periodic())
                .expect("valid term");
        assert_eq!(term.atom_count(), 2);
        assert_eq!(term.charge_c(0), Some(0.5e-19));
    }

    #[test]
    fn a_pair_beyond_the_cutoff_contributes_nothing() {
        let term = ElectrostaticTerm::from_charges(
            &[1.0e-19, -1.0e-19],
            cutoff(),
            PeriodicBox::non_periodic(),
        )
        .expect("valid term");
        let far_m = vec![0.0, 0.0, 0.0, 5.0e-9, 0.0, 0.0];
        assert_eq!(term.energy_j(&far_m).expect("valid"), 0.0);
    }

    #[test]
    fn a_pair_at_the_cutoff_is_switched_off() {
        let term = ElectrostaticTerm::from_charges(
            &[1.0e-19, -1.0e-19],
            cutoff(),
            PeriodicBox::non_periodic(),
        )
        .expect("valid term");
        let positions_m = vec![0.0, 0.0, 0.0, 1.0e-9, 0.0, 0.0];
        assert_eq!(term.energy_j(&positions_m).expect("valid"), 0.0);
        assert!(
            term.energy_j(&[0.0, 0.0, 0.0, 3.0e-10, 0.0, 0.0])
                .expect("valid")
                < 0.0
        );
    }

    #[test]
    fn the_pair_list_path_rejects_a_bad_list() {
        let term = ElectrostaticTerm::from_charges(
            &[1.0e-19, -1.0e-19],
            cutoff(),
            PeriodicBox::non_periodic(),
        )
        .expect("valid term");
        let positions_m = vec![0.0, 0.0, 0.0, 3.0e-10, 0.0, 0.0];
        assert!(matches!(
            term.energy_and_gradient_from_pairs_j(&positions_m, &[0, 1, 0]),
            Err(EngineError::PairListSizeMismatch { .. })
        ));
        assert!(matches!(
            term.energy_and_gradient_from_pairs_j(&positions_m, &[0, 2]),
            Err(EngineError::PairIndexOutOfBounds { .. })
        ));
    }

    #[test]
    fn the_analytic_gradient_matches_a_central_finite_difference() {
        let term =
            ElectrostaticTerm::from_charges(&charges_c(), cutoff(), PeriodicBox::non_periodic())
                .expect("valid term");
        let positions_m = vec![
            0.0, 0.0, 0.0, 3.4e-10, 0.0, 0.0, 0.0, 3.6e-10, 0.0, 0.0, 0.0, 9.0e-10,
        ];
        let analytic = term.gradient_j_per_m(&positions_m).expect("valid");
        let step_m = 1.0e-14;
        let finite_difference =
            central_difference(|p| term.energy_j(p).expect("valid"), &positions_m, step_m);
        let mut max_error_n = 0.0_f64;
        let rtol = 1.0e-6;
        let atol_n = 1.0e-18;
        for index in 0..positions_m.len() {
            let error_n = (analytic[index] - finite_difference[index]).abs();
            max_error_n = max_error_n.max(error_n);
            assert!(
                error_n <= atol_n + rtol * finite_difference[index].abs(),
                "coordinate {index}: analytic {} finite difference {} error {}",
                analytic[index],
                finite_difference[index],
                error_n
            );
        }
        println!("electrostatic max gradient error: {max_error_n:e} N");
    }
}
