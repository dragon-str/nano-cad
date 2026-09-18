//! A four-lane non-bonded pair kernel for Buckingham plus electrostatic terms.
//!
//! The kernel evaluates one flat pair list with two paths: a scalar fallback
//! and a manually unrolled four-lane path. The four-lane path loads four pairs
//! into fixed-size `[f64; 4]` arrays and computes the separation vectors and
//! squared distances with constant-trip-count loops, so the optimizer can
//! auto-vectorize the arithmetic. The per-pair finishing step, the switching
//! function and the scatter stay scalar, because `f64::exp` and `round` have
//! no portable vector form in this workspace.
//!
//! # Why `std::simd` is not used yet
//!
//! Portable SIMD (`std::simd`) is nightly-only. This workspace is stable:
//! rustc 1.98 with a declared MSRV of 1.85. `std::simd` is not stable on
//! 1.85, so the crate cannot depend on it. The task also forbids a new
//! dependency, so a third-party SIMD wrapper is out of scope. When `std::simd`
//! stabilizes, the lane arrays below are the natural place to swap in
//! `f64x4` with a masked gather, an explicit `simd_sqrt`, and a vector
//! exponential approximation.
//!
//! The two paths add each pair's contribution in the same pair order, so they
//! agree bit for bit. The test below still compares them with a tolerance, as
//! the task requires.

use crate::electrostatic::COULOMB_CONSTANT_N_M2_PER_C2;
use crate::error::EngineError;
use crate::geometry::{atom_position_m, dot, sub};
use crate::lennard_jones::lj_energy_and_slope_j;
use crate::nonbonded::{Cutoff, PeriodicBox};

/// The number of lanes in the unrolled kernel.
pub(crate) const LANES: usize = 4;

/// The van der Waals view of a kernel.
pub(crate) struct KernelVanDerWaals<'a> {
    pub(crate) a_j: &'a [f64],
    pub(crate) b_per_m: &'a [f64],
    pub(crate) c_j_m6: &'a [f64],
    pub(crate) cutoff: Cutoff,
}

/// The Lennard-Jones view of a kernel.
pub(crate) struct KernelLennardJones<'a> {
    pub(crate) epsilon_j: &'a [f64],
    pub(crate) sigma_m: &'a [f64],
    pub(crate) cutoff: Cutoff,
}

/// The electrostatic view of a kernel.
pub(crate) struct KernelElectrostatic<'a> {
    pub(crate) charges_c: &'a [f64],
    pub(crate) cutoff: Cutoff,
}

/// A kernel over one position buffer and the active non-bonded terms.
pub(crate) struct NonbondedKernel<'a> {
    pub(crate) positions_m: &'a [f64],
    pub(crate) periodic_box: PeriodicBox,
    pub(crate) van_der_waals: Option<KernelVanDerWaals<'a>>,
    pub(crate) lennard_jones: Option<KernelLennardJones<'a>>,
    pub(crate) electrostatic: Option<KernelElectrostatic<'a>>,
}

impl NonbondedKernel<'_> {
    /// Accumulates a flat pair list through the scalar fallback.
    #[cfg(test)]
    pub(crate) fn accumulate_scalar(
        &self,
        pairs: &[u32],
        energy_j: &mut f64,
        gradient: &mut [f64],
    ) -> Result<(), EngineError> {
        let pair_count = pairs.len() / 2;
        for pair in 0..pair_count {
            self.accumulate_pair(pairs[2 * pair], pairs[2 * pair + 1], energy_j, gradient)?;
        }
        Ok(())
    }

    /// Accumulates a flat pair list through the four-lane path.
    ///
    /// Groups of four pairs run through the lane arrays. A remainder of one to
    /// three pairs runs through the scalar fallback.
    pub(crate) fn accumulate_lanes4(
        &self,
        pairs: &[u32],
        energy_j: &mut f64,
        gradient: &mut [f64],
    ) -> Result<(), EngineError> {
        let pair_count = pairs.len() / 2;
        let mut pair = 0;
        while pair + LANES <= pair_count {
            self.accumulate_group4(&pairs[2 * pair..2 * (pair + LANES)], energy_j, gradient)?;
            pair += LANES;
        }
        while pair < pair_count {
            self.accumulate_pair(pairs[2 * pair], pairs[2 * pair + 1], energy_j, gradient)?;
            pair += 1;
        }
        Ok(())
    }

    fn accumulate_group4(
        &self,
        pairs: &[u32],
        energy_j: &mut f64,
        gradient: &mut [f64],
    ) -> Result<(), EngineError> {
        let mut index_i = [0_u32; LANES];
        let mut index_j = [0_u32; LANES];
        let mut delta_x = [0.0_f64; LANES];
        let mut delta_y = [0.0_f64; LANES];
        let mut delta_z = [0.0_f64; LANES];
        for lane in 0..LANES {
            let i = pairs[2 * lane];
            let j = pairs[2 * lane + 1];
            index_i[lane] = i;
            index_j[lane] = j;
            let delta_m = self.periodic_box.minimum_image(sub(
                atom_position_m(self.positions_m, i),
                atom_position_m(self.positions_m, j),
            ));
            delta_x[lane] = delta_m[0];
            delta_y[lane] = delta_m[1];
            delta_z[lane] = delta_m[2];
        }

        let mut r_sq_m2 = [0.0_f64; LANES];
        for lane in 0..LANES {
            r_sq_m2[lane] = delta_x[lane] * delta_x[lane]
                + delta_y[lane] * delta_y[lane]
                + delta_z[lane] * delta_z[lane];
        }

        for lane in 0..LANES {
            let delta_m = [delta_x[lane], delta_y[lane], delta_z[lane]];
            self.finish_pair(
                r_sq_m2[lane],
                delta_m,
                index_i[lane],
                index_j[lane],
                energy_j,
                gradient,
            )?;
        }
        Ok(())
    }

    fn accumulate_pair(
        &self,
        i: u32,
        j: u32,
        energy_j: &mut f64,
        gradient: &mut [f64],
    ) -> Result<(), EngineError> {
        let delta_m = self.periodic_box.minimum_image(sub(
            atom_position_m(self.positions_m, i),
            atom_position_m(self.positions_m, j),
        ));
        let r_sq_m2 = dot(delta_m, delta_m);
        self.finish_pair(r_sq_m2, delta_m, i, j, energy_j, gradient)
    }

    fn finish_pair(
        &self,
        r_sq_m2: f64,
        delta_m: [f64; 3],
        i: u32,
        j: u32,
        energy_j: &mut f64,
        gradient: &mut [f64],
    ) -> Result<(), EngineError> {
        let van_der_waals = self
            .van_der_waals
            .as_ref()
            .filter(|term| (i as usize) < term.a_j.len() && (j as usize) < term.a_j.len());
        let lennard_jones = self.lennard_jones.as_ref().filter(|term| {
            (i as usize) < term.epsilon_j.len() && (j as usize) < term.epsilon_j.len()
        });
        let electrostatic = self.electrostatic.as_ref().filter(|term| {
            (i as usize) < term.charges_c.len() && (j as usize) < term.charges_c.len()
        });
        if van_der_waals.is_none() && lennard_jones.is_none() && electrostatic.is_none() {
            return Ok(());
        }
        if r_sq_m2 <= 0.0 {
            return Err(EngineError::CoincidentNonbondedAtoms { i, j });
        }
        let r_m = r_sq_m2.sqrt();

        let mut switched_energy_j = 0.0;
        let mut du_dr_n = 0.0;
        if let Some(term) = van_der_waals {
            let cutoff_sq_m2 = term.cutoff.cutoff_m() * term.cutoff.cutoff_m();
            if r_sq_m2 < cutoff_sq_m2 {
                let a_ij = (term.a_j[i as usize] * term.a_j[j as usize]).sqrt();
                let b_ij = 0.5 * (term.b_per_m[i as usize] + term.b_per_m[j as usize]);
                let c_ij = (term.c_j_m6[i as usize] * term.c_j_m6[j as usize]).sqrt();
                let r2 = r_m * r_m;
                let r6 = r2 * r2 * r2;
                let r7 = r6 * r_m;
                let exp_term = (-b_ij * r_m).exp();
                let value_j = a_ij * exp_term - c_ij / r6;
                let value_prime_n = -a_ij * b_ij * exp_term + 6.0 * c_ij / r7;
                let (switch_value, switch_derivative) = term.cutoff.switch_value(r_m);
                switched_energy_j += switch_value * value_j;
                du_dr_n += switch_value * value_prime_n + switch_derivative * value_j;
            }
        }
        if let Some(term) = lennard_jones {
            let cutoff_sq_m2 = term.cutoff.cutoff_m() * term.cutoff.cutoff_m();
            if r_sq_m2 < cutoff_sq_m2 {
                let epsilon_ij = (term.epsilon_j[i as usize] * term.epsilon_j[j as usize]).sqrt();
                let sigma_ij = (term.sigma_m[i as usize] * term.sigma_m[j as usize]).sqrt();
                let (value_j, value_prime_n) =
                    lj_energy_and_slope_j(epsilon_ij, sigma_ij, r_sq_m2, r_m);
                let (switch_value, switch_derivative) = term.cutoff.switch_value(r_m);
                switched_energy_j += switch_value * value_j;
                du_dr_n += switch_value * value_prime_n + switch_derivative * value_j;
            }
        }
        if let Some(term) = electrostatic {
            let cutoff_sq_m2 = term.cutoff.cutoff_m() * term.cutoff.cutoff_m();
            if r_sq_m2 < cutoff_sq_m2 {
                let coulomb_j = COULOMB_CONSTANT_N_M2_PER_C2
                    * term.charges_c[i as usize]
                    * term.charges_c[j as usize];
                let value_j = coulomb_j / r_m;
                let value_prime_n = -coulomb_j / (r_m * r_m);
                let (switch_value, switch_derivative) = term.cutoff.switch_value(r_m);
                switched_energy_j += switch_value * value_j;
                du_dr_n += switch_value * value_prime_n + switch_derivative * value_j;
            }
        }
        if switched_energy_j == 0.0 && du_dr_n == 0.0 {
            return Ok(());
        }

        let factor_n_per_m = du_dr_n / r_m;
        let base_i = i as usize * 3;
        let base_j = j as usize * 3;
        for axis in 0..3 {
            let contribution_n = factor_n_per_m * delta_m[axis];
            gradient[base_i + axis] += contribution_n;
            gradient[base_j + axis] -= contribution_n;
        }
        *energy_j += switched_energy_j;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::Rng;

    fn kernel_terms(
        atom_count: usize,
        seed: u64,
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Cutoff, PeriodicBox) {
        let mut rng = Rng::new(seed);
        let box_m = 3.0e-9;
        let mut positions_m = Vec::with_capacity(3 * atom_count);
        for _ in 0..atom_count {
            positions_m.push(rng.next_f64() * box_m);
            positions_m.push(rng.next_f64() * box_m);
            positions_m.push(rng.next_f64() * box_m);
        }
        let a_j: Vec<f64> = (0..atom_count)
            .map(|atom| 2.2e-19 + 1.0e-21 * atom as f64)
            .collect();
        let b_per_m: Vec<f64> = (0..atom_count).map(|_| 4.2e10).collect();
        let c_j_m6: Vec<f64> = (0..atom_count)
            .map(|atom| 8.0e-79 + 1.0e-80 * atom as f64)
            .collect();
        let cutoff = Cutoff::new(1.2e-9, 0.9e-9).expect("valid cutoff");
        let periodic_box = PeriodicBox::new([box_m, box_m, box_m]).expect("valid box");
        (positions_m, a_j, b_per_m, c_j_m6, cutoff, periodic_box)
    }

    #[test]
    fn the_four_lane_path_matches_the_scalar_path_on_random_pairs() {
        let atom_count = 40;
        let (positions_m, a_j, b_per_m, c_j_m6, cutoff, periodic_box) =
            kernel_terms(atom_count, 0x51D);
        let charges_c: Vec<f64> = (0..atom_count)
            .map(|atom| {
                let magnitude = 0.5e-19 + 0.1e-19 * atom as f64;
                if atom % 2 == 0 {
                    magnitude
                } else {
                    -magnitude
                }
            })
            .collect();
        let kernel = NonbondedKernel {
            positions_m: &positions_m,
            periodic_box,
            van_der_waals: Some(KernelVanDerWaals {
                a_j: &a_j,
                b_per_m: &b_per_m,
                c_j_m6: &c_j_m6,
                cutoff,
            }),
            lennard_jones: None,
            electrostatic: Some(KernelElectrostatic {
                charges_c: &charges_c,
                cutoff,
            }),
        };

        let mut rng = Rng::new(0xBEEF);
        let mut pairs = Vec::new();
        for _ in 0..143 {
            let i = (rng.next_f64() * atom_count as f64) as u32;
            let mut j = (rng.next_f64() * atom_count as f64) as u32;
            if j == i {
                j = (j + 1) % atom_count as u32;
            }
            pairs.push(i);
            pairs.push(j);
        }

        let mut scalar_energy_j = 0.0;
        let mut scalar_gradient = vec![0.0; 3 * atom_count];
        kernel
            .accumulate_scalar(&pairs, &mut scalar_energy_j, &mut scalar_gradient)
            .expect("valid scalar");

        let mut lane_energy_j = 0.0;
        let mut lane_gradient = vec![0.0; 3 * atom_count];
        kernel
            .accumulate_lanes4(&pairs, &mut lane_energy_j, &mut lane_gradient)
            .expect("valid lanes");

        let energy_error_j = (scalar_energy_j - lane_energy_j).abs();
        assert!(energy_error_j <= 1.0e-30 + 1.0e-13 * scalar_energy_j.abs());
        let mut max_error_n = 0.0_f64;
        for coordinate in 0..scalar_gradient.len() {
            let error_n = (scalar_gradient[coordinate] - lane_gradient[coordinate]).abs();
            max_error_n = max_error_n.max(error_n);
            assert!(
                error_n <= 1.0e-24 + 1.0e-13 * scalar_gradient[coordinate].abs(),
                "coordinate {coordinate}: scalar {} lane {} error {}",
                scalar_gradient[coordinate],
                lane_gradient[coordinate],
                error_n
            );
        }
        println!(
            "four-lane kernel: max gradient error {max_error_n:e} N, energy error {energy_error_j:e} J over {} pairs",
            pairs.len() / 2
        );
    }
}
