//! A Lennard-Jones 12-6 term over all atom pairs, with a cutoff.
//!
//! The pair potential is
//! `U(r) = 4 epsilon_ij [ (sigma_ij / r)^12 - (sigma_ij / r)^6 ]`.
//!
//! A pocket that must bind one molecule and reject another needs a size and
//! a well depth for every element. The table in
//! [`nanocad_model::nonbonded`] supplies them. The Buckingham term of
//! [`crate::van_der_waals`] cannot hold those values, because a Lennard-Jones
//! curve has no exact Buckingham form, so this term exists beside it.
//!
//! Energy is in joules. The gradient is in newtons. The potential is
//! multiplied by the switching function of [`Cutoff`] and reduced by the
//! minimum-image convention of [`PeriodicBox`].

use crate::error::EngineError;
use crate::geometry::{atom_position_m, dot, sub, validate_positions};
use crate::nonbonded::{all_pairs, Cutoff, PeriodicBox};

/// Lennard-Jones parameters for one atom.
///
/// An unlike pair uses the mixing rules of the Universal Force Field,
/// `sigma_ij = sqrt(sigma_i sigma_j)` and `epsilon_ij = sqrt(epsilon_i
/// epsilon_j)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LjParams {
    /// The depth of the well in joules.
    pub epsilon_j: f64,
    /// The distance at which the potential crosses zero, in metres.
    pub sigma_m: f64,
}

impl LjParams {
    /// Builds Lennard-Jones parameters for one atom.
    pub const fn new(epsilon_j: f64, sigma_m: f64) -> Self {
        Self { epsilon_j, sigma_m }
    }
}

/// A Lennard-Jones term over all atom pairs, with a cutoff.
#[derive(Clone, Debug, PartialEq)]
pub struct LennardJonesTerm {
    atom_count: usize,
    epsilon_j: Vec<f64>,
    sigma_m: Vec<f64>,
    cutoff: Cutoff,
    periodic_box: PeriodicBox,
}

impl LennardJonesTerm {
    /// Creates an empty term with no atoms.
    pub fn new(cutoff: Cutoff, periodic_box: PeriodicBox) -> Result<Self, EngineError> {
        cutoff.validate_box(&periodic_box)?;
        Ok(Self {
            atom_count: 0,
            epsilon_j: Vec::new(),
            sigma_m: Vec::new(),
            cutoff,
            periodic_box,
        })
    }

    /// Builds a term from one parameter set per atom.
    pub fn from_params(
        params: &[LjParams],
        cutoff: Cutoff,
        periodic_box: PeriodicBox,
    ) -> Result<Self, EngineError> {
        let mut term = Self::new(cutoff, periodic_box)?;
        for param in params {
            term.add_atom(param.epsilon_j, param.sigma_m)?;
        }
        Ok(term)
    }

    /// Appends one atom with its parameters.
    ///
    /// The well depth must be finite and non-negative. The distance must be
    /// finite and positive.
    pub fn add_atom(&mut self, epsilon_j: f64, sigma_m: f64) -> Result<(), EngineError> {
        let atom = self.atom_count;
        if !epsilon_j.is_finite() || epsilon_j < 0.0 {
            return Err(EngineError::InvalidVdwParameters { atom });
        }
        if !sigma_m.is_finite() || sigma_m <= 0.0 {
            return Err(EngineError::InvalidVdwParameters { atom });
        }
        self.epsilon_j.push(epsilon_j);
        self.sigma_m.push(sigma_m);
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

    /// Returns the per-atom well depths in joules.
    pub(crate) fn epsilon_j(&self) -> &[f64] {
        &self.epsilon_j
    }

    /// Returns the per-atom distances in metres.
    pub(crate) fn sigma_m(&self) -> &[f64] {
        &self.sigma_m
    }

    /// Replaces the periodic box. The cutoff must not exceed half of any
    /// periodic length.
    pub fn set_periodic_box(&mut self, periodic_box: PeriodicBox) -> Result<(), EngineError> {
        self.cutoff.validate_box(&periodic_box)?;
        self.periodic_box = periodic_box;
        Ok(())
    }

    /// Returns the total Lennard-Jones energy in joules.
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
    /// This is the reference path. It scans all unordered pairs.
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
    /// `2k + 1`, is one pair. A pair list from a [`crate::VerletList`] is
    /// accepted directly. Pairs beyond the cutoff contribute nothing.
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
        let epsilon_ij = (self.epsilon_j[i as usize] * self.epsilon_j[j as usize]).sqrt();
        let sigma_ij = (self.sigma_m[i as usize] * self.sigma_m[j as usize]).sqrt();
        let (value_j, value_prime_n) = lj_energy_and_slope_j(epsilon_ij, sigma_ij, r_sq_m2, r_m);
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

/// Returns the Lennard-Jones energy in joules and its slope `dU/dr` in
/// newtons for one pair.
///
/// The ratio `sigma / r` is formed from the squared quantities, because
/// `(sigma/r)^12` overflows far sooner than `(sigma^2 / r^2)^6` and the
/// squared form is exact for the same values.
pub(crate) fn lj_energy_and_slope_j(
    epsilon_ij: f64,
    sigma_ij: f64,
    r_sq_m2: f64,
    r_m: f64,
) -> (f64, f64) {
    let ratio_sq = (sigma_ij * sigma_ij) / r_sq_m2;
    let inverse_6 = ratio_sq * ratio_sq * ratio_sq;
    let inverse_12 = inverse_6 * inverse_6;
    let four_epsilon = 4.0 * epsilon_ij;
    let value_j = four_epsilon * (inverse_12 - inverse_6);
    let value_prime_n = four_epsilon * (6.0 * inverse_6 - 12.0 * inverse_12) / r_m;
    (value_j, value_prime_n)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::{central_difference, Rng};

    fn lj_params() -> Vec<LjParams> {
        vec![
            LjParams::new(1.5e-22, 3.4e-10),
            LjParams::new(2.5e-22, 3.0e-10),
            LjParams::new(1.0e-22, 3.8e-10),
        ]
    }

    fn cutoff() -> Cutoff {
        Cutoff::new(2.0e-9, 1.8e-9).expect("valid cutoff")
    }

    fn positions() -> Vec<f64> {
        vec![0.0, 0.0, 0.0, 4.0e-10, 0.0, 0.0, 0.0, 5.0e-10, 0.0]
    }

    #[test]
    fn a_pair_of_one_element_matches_the_analytic_curve() {
        let term = LennardJonesTerm::from_params(
            &[
                LjParams::new(1.0e-21, 3.0e-10),
                LjParams::new(1.0e-21, 3.0e-10),
            ],
            cutoff(),
            PeriodicBox::non_periodic(),
        )
        .expect("valid term");
        let positions_m = vec![0.0, 0.0, 0.0, 3.0e-10, 0.0, 0.0];
        let energy_j = term.energy_j(&positions_m).expect("valid positions");
        assert!(close(energy_j, -1.0e-21, 1.0e-9), "{energy_j}");
    }

    #[test]
    fn the_well_minimum_sits_at_the_sigma_scaled_distance() {
        let sigma_m = 3.0e-10;
        let epsilon_j = 1.0e-21;
        let term = LennardJonesTerm::from_params(
            &[
                LjParams::new(epsilon_j, sigma_m),
                LjParams::new(epsilon_j, sigma_m),
            ],
            Cutoff::new(2.0e-9, 1.8e-9).expect("valid cutoff"),
            PeriodicBox::non_periodic(),
        )
        .expect("valid term");
        let minimum_m = sigma_m * 2.0_f64.powf(1.0 / 6.0);
        let mut best_energy_j = f64::INFINITY;
        let mut best_r_m = 0.0;
        for step in 1..400 {
            let r_m = 2.0e-10 + 1.0e-12 * step as f64;
            let positions_m = vec![0.0, 0.0, 0.0, r_m, 0.0, 0.0];
            let energy_j = term.energy_j(&positions_m).expect("valid positions");
            if energy_j < best_energy_j {
                best_energy_j = energy_j;
                best_r_m = r_m;
            }
        }
        assert!(close(best_energy_j, -epsilon_j, 1.0e-2), "{best_energy_j}");
        assert!(close(best_r_m, minimum_m, 5.0e-3), "{best_r_m}");
    }

    #[test]
    fn the_gradient_matches_a_central_difference() {
        let term =
            LennardJonesTerm::from_params(&lj_params(), cutoff(), PeriodicBox::non_periodic())
                .expect("valid term");
        let positions_m = positions();
        let step_m = 1.0e-14;
        let numeric =
            central_difference(|p| term.energy_j(p).expect("valid"), &positions_m, step_m);
        let analytic = term.gradient_j_per_m(&positions_m).expect("valid");
        let error = relative_error(&numeric, &analytic);
        assert!(error < 1.0e-4, "{error}");
    }

    #[test]
    fn the_pair_list_path_matches_the_all_pairs_path() {
        let term =
            LennardJonesTerm::from_params(&lj_params(), cutoff(), PeriodicBox::non_periodic())
                .expect("valid term");
        let positions_m = positions();
        let pairs = all_pairs(3);
        let (list_energy_j, list_gradient) = term
            .energy_and_gradient_from_pairs_j(&positions_m, &pairs)
            .expect("valid");
        let (all_energy_j, all_gradient) = term.energy_and_gradient_j(&positions_m).expect("valid");
        assert_eq!(list_energy_j, all_energy_j);
        assert_eq!(list_gradient, all_gradient);
    }

    #[test]
    fn a_pair_beyond_the_cutoff_contributes_nothing() {
        let term = LennardJonesTerm::from_params(
            &[
                LjParams::new(1.0e-21, 3.0e-10),
                LjParams::new(1.0e-21, 3.0e-10),
            ],
            cutoff(),
            PeriodicBox::non_periodic(),
        )
        .expect("valid term");
        let positions_m = vec![0.0, 0.0, 0.0, 5.0e-9, 0.0, 0.0];
        let pair = term
            .pair_energy_and_gradient_j(&positions_m, 0, 1)
            .expect("valid pair");
        assert!(pair.is_none());
        assert_eq!(term.energy_j(&positions_m).expect("valid"), 0.0);
    }

    #[test]
    fn a_pair_of_one_element_recovers_its_parameters() {
        let term = LennardJonesTerm::from_params(
            &[
                LjParams::new(1.0e-21, 3.0e-10),
                LjParams::new(1.0e-21, 3.0e-10),
            ],
            cutoff(),
            PeriodicBox::non_periodic(),
        )
        .expect("valid term");
        let (value_j, _) = term
            .pair_energy_and_gradient_j(&positions(), 0, 1)
            .expect("valid pair")
            .expect("inside cutoff");
        assert!(value_j.is_finite());
    }

    #[test]
    fn a_negative_well_depth_is_refused() {
        let mut term =
            LennardJonesTerm::new(cutoff(), PeriodicBox::non_periodic()).expect("valid term");
        assert!(term.add_atom(-1.0, 3.0e-10).is_err());
    }

    #[test]
    fn a_zero_distance_is_refused() {
        let mut term =
            LennardJonesTerm::new(cutoff(), PeriodicBox::non_periodic()).expect("valid term");
        assert!(term.add_atom(1.0e-21, 0.0).is_err());
    }

    #[test]
    fn coincident_atoms_are_refused() {
        let term = LennardJonesTerm::from_params(
            &[
                LjParams::new(1.0e-21, 3.0e-10),
                LjParams::new(1.0e-21, 3.0e-10),
            ],
            cutoff(),
            PeriodicBox::non_periodic(),
        )
        .expect("valid term");
        let positions_m = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let pair = term.pair_energy_and_gradient_j(&positions_m, 0, 1);
        assert!(matches!(
            pair,
            Err(EngineError::CoincidentNonbondedAtoms { .. })
        ));
    }

    #[test]
    fn the_repulsion_grows_as_the_atoms_approach() {
        let term = LennardJonesTerm::from_params(
            &[
                LjParams::new(1.0e-21, 3.0e-10),
                LjParams::new(1.0e-21, 3.0e-10),
            ],
            Cutoff::new(2.0e-9, 1.8e-9).expect("valid cutoff"),
            PeriodicBox::non_periodic(),
        )
        .expect("valid term");
        let far_j = term
            .energy_j(&[0.0, 0.0, 0.0, 3.5e-10, 0.0, 0.0])
            .expect("valid");
        let near_j = term
            .energy_j(&[0.0, 0.0, 0.0, 2.5e-10, 0.0, 0.0])
            .expect("valid");
        assert!(near_j > far_j);
    }

    #[test]
    fn an_odd_pair_list_is_refused() {
        let term =
            LennardJonesTerm::from_params(&lj_params(), cutoff(), PeriodicBox::non_periodic())
                .expect("valid term");
        let result = term.energy_and_gradient_from_pairs_j(&positions(), &[0, 1, 2]);
        assert!(matches!(
            result,
            Err(EngineError::PairListSizeMismatch { .. })
        ));
    }

    #[test]
    fn a_pair_index_outside_the_atom_count_is_refused() {
        let term =
            LennardJonesTerm::from_params(&lj_params(), cutoff(), PeriodicBox::non_periodic())
                .expect("valid term");
        let result = term.energy_and_gradient_from_pairs_j(&positions(), &[0, 7]);
        assert!(matches!(
            result,
            Err(EngineError::PairIndexOutOfBounds { .. })
        ));
    }

    #[test]
    fn a_random_configuration_matches_a_central_difference() {
        let term = LennardJonesTerm::from_params(
            &lj_params(),
            Cutoff::new(2.0e-9, 1.8e-9).expect("valid cutoff"),
            PeriodicBox::non_periodic(),
        )
        .expect("valid term");
        let mut rng = Rng::new(7);
        let mut positions_m = Vec::with_capacity(9);
        for _ in 0..3 {
            for _ in 0..3 {
                positions_m.push(3.0e-10 + 4.0e-10 * rng.next_f64());
            }
        }
        let numeric =
            central_difference(|p| term.energy_j(p).expect("valid"), &positions_m, 1.0e-14);
        let analytic = term.gradient_j_per_m(&positions_m).expect("valid");
        let error = relative_error(&numeric, &analytic);
        assert!(error < 1.0e-4, "{error}");
    }

    fn relative_error(numeric: &[f64], analytic: &[f64]) -> f64 {
        let mut worst: f64 = 0.0;
        for (left, right) in numeric.iter().zip(analytic.iter()) {
            let scale = left.abs().max(right.abs()).max(1.0e-30);
            worst = worst.max((left - right).abs() / scale);
        }
        worst
    }

    fn close(left: f64, right: f64, tolerance: f64) -> bool {
        (left - right).abs() <= tolerance * right.abs().max(1.0)
    }
}
