//! The Hessian of the potential energy and its normal modes.
//!
//! The Hessian is the matrix of second derivatives of the energy. Its
//! eigenvalues, divided by the atomic masses, give the squares of the normal
//! mode frequencies. A negative eigenvalue is an unstable mode: the
//! configuration is not a minimum.
//!
//! The Hessian here is a finite difference of the analytic gradient. It needs
//! two gradient evaluations for each coordinate, so it suits a small cluster,
//! not a whole device.

use std::f64::consts::TAU;

use crate::error::EngineError;
use crate::system::System;

/// A normal-mode spectrum.
///
/// The eigenvalues are mass-weighted, in reciprocal seconds squared. The
/// frequencies are in hertz.
#[derive(Clone, Debug, PartialEq)]
pub struct Spectrum {
    /// The mass-weighted eigenvalues, sorted from the smallest.
    pub eigenvalues: Vec<f64>,
    /// The frequency of each eigenvalue in hertz. It is zero for an unstable
    /// mode.
    pub frequencies_hz: Vec<f64>,
    /// The number of eigenvalues below the negative tolerance.
    pub unstable_count: usize,
}

/// The eigenvalue below which a mode counts as unstable, in reciprocal
/// seconds squared.
pub const UNSTABLE_TOLERANCE_PER_S2: f64 = -1.0e20;

/// Returns the eigenvalues of a symmetric matrix, sorted from the smallest.
///
/// The matrix is flat and row-major with `n` rows. The method is the cyclic
/// Jacobi rotation. It converges for a symmetric matrix. The input is not
/// modified.
pub fn symmetric_eigenvalues(matrix: &[f64], n: usize) -> Result<Vec<f64>, EngineError> {
    if n == 0 {
        return Err(EngineError::EmptySpectrum);
    }
    if matrix.len() != n * n {
        return Err(EngineError::BufferSizeMismatch {
            len: matrix.len(),
            expected: n * n,
        });
    }
    let mut a = matrix.to_vec();
    let max_sweeps = 100;
    for _ in 0..max_sweeps {
        let mut off_diagonal = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                off_diagonal += a[i * n + j] * a[i * n + j];
            }
        }
        if off_diagonal <= 1.0e-300 {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = a[p * n + q];
                if apq == 0.0 {
                    continue;
                }
                let app = a[p * n + p];
                let aqq = a[q * n + q];
                let theta = (aqq - app) / (2.0 * apq);
                let sign = if theta >= 0.0 { 1.0 } else { -1.0 };
                let t = sign / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let akp = a[k * n + p];
                    let akq = a[k * n + q];
                    a[k * n + p] = c * akp - s * akq;
                    a[k * n + q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let apk = a[p * n + k];
                    let aqk = a[q * n + k];
                    a[p * n + k] = c * apk - s * aqk;
                    a[q * n + k] = s * apk + c * aqk;
                }
            }
        }
    }
    let mut eigenvalues: Vec<f64> = (0..n).map(|i| a[i * n + i]).collect();
    eigenvalues.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    Ok(eigenvalues)
}

/// Returns the Hessian of the energy, as a flat row-major matrix.
///
/// Entry `(i, j)` is `d2U/dx_i dx_j`. The method is a central finite
/// difference of the gradient with the given step in metres. The step must be
/// much smaller than a bond, or the perturbation crosses a bond minimum. The
/// matrix is symmetrised, because the finite difference is not exactly
/// symmetric.
pub fn hessian_finite_difference(
    system: &mut System,
    positions_m: &[f64],
    step_m: f64,
) -> Result<Vec<f64>, EngineError> {
    if !step_m.is_finite() || step_m <= 0.0 {
        return Err(EngineError::NonPositiveDifferenceStep { step_m });
    }
    let coordinates = 3 * system.atom_count();
    let mut hessian = vec![0.0; coordinates * coordinates];
    let mut plus = positions_m.to_vec();
    let mut minus = positions_m.to_vec();
    system.rebuild_neighbors(positions_m)?;
    for i in 0..coordinates {
        plus[i] = positions_m[i] + step_m;
        minus[i] = positions_m[i] - step_m;
        let gradient_plus = system.gradient_j_per_m(&plus)?;
        let gradient_minus = system.gradient_j_per_m(&minus)?;
        plus[i] = positions_m[i];
        minus[i] = positions_m[i];
        let scale = 1.0 / (2.0 * step_m);
        for j in 0..coordinates {
            hessian[i * coordinates + j] = (gradient_plus[j] - gradient_minus[j]) * scale;
        }
    }
    for i in 0..coordinates {
        for j in (i + 1)..coordinates {
            let average = 0.5 * (hessian[i * coordinates + j] + hessian[j * coordinates + i]);
            hessian[i * coordinates + j] = average;
            hessian[j * coordinates + i] = average;
        }
    }
    Ok(hessian)
}

/// Returns the mass-weighted Hessian.
///
/// Entry `(i, j)` is divided by the square root of the product of the two
/// atom masses. The eigenvalues are then in reciprocal seconds squared.
pub fn mass_weighted(hessian: &[f64], masses_kg: &[f64]) -> Result<Vec<f64>, EngineError> {
    let coordinates = 3 * masses_kg.len();
    if hessian.len() != coordinates * coordinates {
        return Err(EngineError::BufferSizeMismatch {
            len: hessian.len(),
            expected: coordinates * coordinates,
        });
    }
    let mut weighted = vec![0.0; hessian.len()];
    for i in 0..coordinates {
        for j in 0..coordinates {
            let mi = masses_kg[i / 3];
            let mj = masses_kg[j / 3];
            weighted[i * coordinates + j] = hessian[i * coordinates + j] / (mi * mj).sqrt();
        }
    }
    Ok(weighted)
}

/// Builds a spectrum from mass-weighted eigenvalues.
///
/// The frequency of a positive eigenvalue is `sqrt(eigenvalue) / 2 pi`. A
/// negative eigenvalue is an unstable mode, and its frequency is zero.
pub fn spectrum(eigenvalues: Vec<f64>) -> Spectrum {
    let mut unstable_count = 0;
    let frequencies_hz = eigenvalues
        .iter()
        .map(|&eigenvalue| {
            if eigenvalue < UNSTABLE_TOLERANCE_PER_S2 {
                unstable_count += 1;
                0.0
            } else if eigenvalue > 0.0 {
                eigenvalue.sqrt() / TAU
            } else {
                0.0
            }
        })
        .collect();
    Spectrum {
        eigenvalues,
        frequencies_hz,
        unstable_count,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bond_stretch::BondStretchTerm;

    const CARBON_MASS_KG: f64 = 1.992_646_879_92e-26;
    const BOND_K_N_PER_M: f64 = 300.0;
    const BOND_LENGTH_M: f64 = 1.5e-10;

    /// Jacobi returns the known eigenvalues of a symmetric two by two matrix.
    #[test]
    fn jacobi_solves_a_two_by_two_matrix() {
        let matrix = vec![2.0, 1.0, 1.0, 2.0];
        let eigenvalues = symmetric_eigenvalues(&matrix, 2).unwrap();
        assert!((eigenvalues[0] - 1.0).abs() < 1.0e-12);
        assert!((eigenvalues[1] - 3.0).abs() < 1.0e-12);
    }

    /// A diatomic bond has one nonzero mode at sqrt(2 k / m) / 2 pi.
    #[test]
    fn a_diatomic_bond_has_one_mode() {
        let mut system = System::with_masses_kg(2, &[CARBON_MASS_KG; 2]).unwrap();
        let mut bonds = BondStretchTerm::new();
        bonds.add_bond(0, 1, BOND_K_N_PER_M, BOND_LENGTH_M).unwrap();
        system.set_bond_stretch(bonds).unwrap();

        let positions_m = vec![0.0, 0.0, 0.0, BOND_LENGTH_M, 0.0, 0.0];
        let step_m = 1.0e-12;
        let hessian = hessian_finite_difference(&mut system, &positions_m, step_m).unwrap();
        let weighted = mass_weighted(&hessian, &[CARBON_MASS_KG; 2]).unwrap();
        let eigenvalues = symmetric_eigenvalues(&weighted, 6).unwrap();
        let found = spectrum(eigenvalues);
        assert_eq!(found.unstable_count, 0);

        let expected = (2.0 * BOND_K_N_PER_M / CARBON_MASS_KG).sqrt() / TAU;
        let largest = found.frequencies_hz[found.frequencies_hz.len() - 1];
        assert!(
            (largest - expected).abs() / expected < 1.0e-6,
            "the mode is {largest} Hz, not {expected} Hz"
        );
        for frequency in &found.frequencies_hz[..5] {
            assert!(
                *frequency < 0.01 * expected,
                "a free mode is not small: {frequency} Hz against {expected} Hz"
            );
        }
    }

    /// A zero step is an error.
    #[test]
    fn a_zero_step_is_rejected() {
        let mut system = System::with_masses_kg(2, &[CARBON_MASS_KG; 2]).unwrap();
        let positions_m = vec![0.0; 6];
        let error = hessian_finite_difference(&mut system, &positions_m, 0.0).unwrap_err();
        assert_eq!(
            error,
            EngineError::NonPositiveDifferenceStep { step_m: 0.0 }
        );
    }
}
