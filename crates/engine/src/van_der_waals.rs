use crate::error::EngineError;
use crate::geometry::{atom_position_m, dot, sub, validate_positions};
use crate::nonbonded::{all_pairs, Cutoff, PeriodicBox};

/// Buckingham parameters for one atom.
///
/// An unlike pair uses the combining rules
/// `A_ij = sqrt(A_i A_j)`, `B_ij = 0.5 (B_i + B_j)`, and
/// `C_ij = sqrt(C_i C_j)`. The pair potential is
/// `U(r) = A_ij exp(-B_ij r) - C_ij / r^6`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VdwParams {
    /// Repulsive prefactor in joules.
    pub a_j: f64,
    /// Repulsive decay constant in reciprocal metres.
    pub b_per_m: f64,
    /// Dispersion coefficient in joule cubic metres to the sixth power.
    pub c_j_m6: f64,
}

impl VdwParams {
    /// Builds Buckingham parameters for one atom.
    pub const fn new(a_j: f64, b_per_m: f64, c_j_m6: f64) -> Self {
        Self {
            a_j,
            b_per_m,
            c_j_m6,
        }
    }
}

/// A van der Waals term over all atom pairs, with a cutoff.
///
/// The potential is multiplied by the switching function of [`Cutoff`] and
/// reduced by the minimum-image convention of [`PeriodicBox`]. Energy is in
/// joules. The gradient is in newtons.
#[derive(Clone, Debug, PartialEq)]
pub struct VanDerWaalsTerm {
    atom_count: usize,
    a_j: Vec<f64>,
    b_per_m: Vec<f64>,
    c_j_m6: Vec<f64>,
    cutoff: Cutoff,
    periodic_box: PeriodicBox,
}

impl VanDerWaalsTerm {
    /// Creates an empty term with no atoms.
    pub fn new(cutoff: Cutoff, periodic_box: PeriodicBox) -> Result<Self, EngineError> {
        cutoff.validate_box(&periodic_box)?;
        Ok(Self {
            atom_count: 0,
            a_j: Vec::new(),
            b_per_m: Vec::new(),
            c_j_m6: Vec::new(),
            cutoff,
            periodic_box,
        })
    }

    /// Builds a term from one parameter set per atom.
    pub fn from_params(
        params: &[VdwParams],
        cutoff: Cutoff,
        periodic_box: PeriodicBox,
    ) -> Result<Self, EngineError> {
        let mut term = Self::new(cutoff, periodic_box)?;
        for param in params {
            term.add_atom(param.a_j, param.b_per_m, param.c_j_m6)?;
        }
        Ok(term)
    }

    /// Appends one atom with its parameters.
    ///
    /// Every value must be finite. The repulsive prefactor and the dispersion
    /// coefficient must be non-negative. The decay constant must be
    /// non-negative, and positive when the repulsive prefactor is positive.
    pub fn add_atom(&mut self, a_j: f64, b_per_m: f64, c_j_m6: f64) -> Result<(), EngineError> {
        let atom = self.atom_count;
        if !a_j.is_finite() || a_j < 0.0 {
            return Err(EngineError::InvalidVdwParameters { atom });
        }
        if !b_per_m.is_finite() || b_per_m < 0.0 {
            return Err(EngineError::InvalidVdwParameters { atom });
        }
        if a_j > 0.0 && b_per_m <= 0.0 {
            return Err(EngineError::InvalidVdwParameters { atom });
        }
        if !c_j_m6.is_finite() || c_j_m6 < 0.0 {
            return Err(EngineError::InvalidVdwParameters { atom });
        }
        self.a_j.push(a_j);
        self.b_per_m.push(b_per_m);
        self.c_j_m6.push(c_j_m6);
        self.atom_count += 1;
        Ok(())
    }

    /// Returns the number of atoms.
    pub fn atom_count(&self) -> usize {
        self.atom_count
    }

    /// Returns the cutoff.
    pub fn cutoff(&self) -> &Cutoff {
        &self.cutoff
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

    /// Returns the total van der Waals energy in joules.
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
    /// pair-list path of [`VanDerWaalsTerm::energy_and_gradient_from_pairs_j`]
    /// gives the same result for a pair list that holds every pair within the
    /// cutoff.
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
        let a_ij = (self.a_j[i as usize] * self.a_j[j as usize]).sqrt();
        let b_ij = 0.5 * (self.b_per_m[i as usize] + self.b_per_m[j as usize]);
        let c_ij = (self.c_j_m6[i as usize] * self.c_j_m6[j as usize]).sqrt();
        let exp_term = (-b_ij * r_m).exp();
        let value_j = a_ij * exp_term - c_ij / (r_m * r_m * r_m * r_m * r_m * r_m);
        let value_prime_n =
            -a_ij * b_ij * exp_term + 6.0 * c_ij / (r_m * r_m * r_m * r_m * r_m * r_m * r_m);
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
    use crate::test_support::{central_difference, Rng};

    fn vdw_params() -> Vec<VdwParams> {
        vec![
            VdwParams::new(2.5e-19, 4.0e10, 1.0e-78),
            VdwParams::new(1.8e-19, 4.5e10, 6.0e-79),
            VdwParams::new(2.2e-19, 4.2e10, 8.0e-79),
            VdwParams::new(2.0e-19, 4.1e10, 9.0e-79),
        ]
    }

    fn cutoff() -> Cutoff {
        Cutoff::new(1.0e-9, 0.8e-9).expect("valid cutoff")
    }

    fn two_atom_params() -> Vec<VdwParams> {
        vec![
            VdwParams::new(2.5e-19, 4.0e10, 1.0e-78),
            VdwParams::new(1.8e-19, 4.5e10, 6.0e-79),
        ]
    }

    fn small_molecule(seed: u64) -> Vec<f64> {
        let mut rng = Rng::new(seed);
        let mut positions_m = Vec::new();
        let base_m = [
            [0.0, 0.0, 0.0],
            [3.4e-10, 0.0, 0.0],
            [0.0, 3.6e-10, 0.0],
            [0.0, 0.0, 9.0e-10],
        ];
        for atom in base_m {
            positions_m.push(atom[0] + rng.symmetric(1.0e-12));
            positions_m.push(atom[1] + rng.symmetric(1.0e-12));
            positions_m.push(atom[2] + rng.symmetric(1.0e-12));
        }
        positions_m
    }

    #[test]
    fn an_empty_term_has_zero_energy_and_an_empty_gradient() {
        let term = VanDerWaalsTerm::new(cutoff(), PeriodicBox::non_periodic()).expect("valid");
        assert_eq!(term.atom_count(), 0);
        assert_eq!(term.energy_j(&[]).expect("valid"), 0.0);
        assert!(term.gradient_j_per_m(&[]).expect("valid").is_empty());
    }

    #[test]
    fn invalid_parameters_are_rejected() {
        let mut term = VanDerWaalsTerm::new(cutoff(), PeriodicBox::non_periodic()).expect("valid");
        assert!(matches!(
            term.add_atom(-1.0, 4.0e10, 1.0e-78),
            Err(EngineError::InvalidVdwParameters { .. })
        ));
        assert!(matches!(
            term.add_atom(1.0e-19, 0.0, 1.0e-78),
            Err(EngineError::InvalidVdwParameters { .. })
        ));
        assert!(matches!(
            term.add_atom(1.0e-19, 4.0e10, -1.0e-78),
            Err(EngineError::InvalidVdwParameters { .. })
        ));
    }

    #[test]
    fn a_pair_beyond_the_cutoff_contributes_nothing() {
        let term =
            VanDerWaalsTerm::from_params(&two_atom_params(), cutoff(), PeriodicBox::non_periodic())
                .expect("valid term");
        let far_m = vec![0.0, 0.0, 0.0, 5.0e-9, 0.0, 0.0];
        assert_eq!(term.energy_j(&far_m).expect("valid"), 0.0);
    }

    #[test]
    fn a_pair_at_the_cutoff_is_switched_off() {
        let term =
            VanDerWaalsTerm::from_params(&two_atom_params(), cutoff(), PeriodicBox::non_periodic())
                .expect("valid term");
        let positions_m = vec![0.0, 0.0, 0.0, 1.0e-9, 0.0, 0.0];
        assert_eq!(term.energy_j(&positions_m).expect("valid"), 0.0);
        assert!(
            term.energy_j(&[0.0, 0.0, 0.0, 3.4e-10, 0.0, 0.0])
                .expect("valid")
                < 0.0
        );
    }

    #[test]
    fn the_minimum_image_gives_the_near_pair_distance() {
        let local_cutoff = Cutoff::new(0.15e-9, 0.1e-9).expect("valid cutoff");
        let periodic_box = PeriodicBox::new([3.4e-10, 0.0, 0.0]).expect("valid box");
        let term = VanDerWaalsTerm::from_params(&two_atom_params(), local_cutoff, periodic_box)
            .expect("valid term");
        let non_periodic = VanDerWaalsTerm::from_params(
            &two_atom_params(),
            local_cutoff,
            PeriodicBox::non_periodic(),
        )
        .expect("valid term");
        let positions_m = vec![0.05e-10, 0.0, 0.0, 3.35e-10, 0.0, 0.0];
        let image_energy_j = term.energy_j(&positions_m).expect("valid");
        let direct_energy_j = non_periodic.energy_j(&positions_m).expect("valid");
        assert_eq!(direct_energy_j, 0.0);
        assert!(image_energy_j.abs() > 1.0e-30);
    }

    #[test]
    fn the_pair_list_path_rejects_a_bad_list() {
        let term =
            VanDerWaalsTerm::from_params(&two_atom_params(), cutoff(), PeriodicBox::non_periodic())
                .expect("valid term");
        let positions_m = vec![0.0, 0.0, 0.0, 3.4e-10, 0.0, 0.0];
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
    fn the_pair_list_path_matches_the_all_pairs_path() {
        let term =
            VanDerWaalsTerm::from_params(&vdw_params(), cutoff(), PeriodicBox::non_periodic())
                .expect("valid term");
        let positions_m = small_molecule(0x2);
        let energy_and_gradient_j = term.energy_and_gradient_j(&positions_m).expect("valid");
        let from_pairs = term
            .energy_and_gradient_from_pairs_j(&positions_m, &[0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3])
            .expect("valid");
        assert!((energy_and_gradient_j.0 - from_pairs.0).abs() < 1.0e-30);
    }

    #[test]
    fn the_analytic_gradient_matches_a_central_finite_difference() {
        let term =
            VanDerWaalsTerm::from_params(&vdw_params(), cutoff(), PeriodicBox::non_periodic())
                .expect("valid term");
        let positions_m = small_molecule(0x2);
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
        println!("van der Waals max gradient error: {max_error_n:e} N");
    }
}
