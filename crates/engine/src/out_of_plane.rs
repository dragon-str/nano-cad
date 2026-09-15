use crate::error::EngineError;
use crate::geometry::{
    add_assign, atom_position_m, cross, dot, norm, scale, sub, validate_positions,
};

/// Harmonic out-of-plane parameters for one improper.
///
/// The potential is `U(chi) = 0.5 * k_j_per_rad2 * (chi - chi0_rad)^2`. The
/// signed angle `chi` is between the central-to-out-of-plane bond and the plane
/// of the other three atoms. The force constant is in joules per radian squared.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OutOfPlaneParams {
    /// Force constant in joules per radian squared.
    pub k_j_per_rad2: f64,
    /// Equilibrium out-of-plane angle in radians.
    pub chi0_rad: f64,
}

impl OutOfPlaneParams {
    /// Builds parameters for one improper.
    pub const fn new(k_j_per_rad2: f64, chi0_rad: f64) -> Self {
        Self {
            k_j_per_rad2,
            chi0_rad,
        }
    }
}

/// A harmonic out-of-plane term over a set of atom quads.
///
/// Each quad is `[i, j, k, l]`, where `j` is the central atom and `k` is the
/// out-of-plane atom. Atoms `i`, `j`, and `l` define the reference plane.
/// Energy is in joules. The gradient is the derivative of the energy with
/// respect to each Cartesian coordinate, in newtons.
#[derive(Clone, Debug, PartialEq)]
pub struct OutOfPlaneTerm {
    atom_count: usize,
    quads: Vec<[u32; 4]>,
    k_j_per_rad2: Vec<f64>,
    chi0_rad: Vec<f64>,
}

impl Default for OutOfPlaneTerm {
    fn default() -> Self {
        Self::new()
    }
}

impl OutOfPlaneTerm {
    /// Creates an empty term with no impropers and no atoms.
    pub const fn new() -> Self {
        Self {
            atom_count: 0,
            quads: Vec::new(),
            k_j_per_rad2: Vec::new(),
            chi0_rad: Vec::new(),
        }
    }

    /// Builds a term from a list of quads and one parameter set per quad.
    pub fn from_quads(
        quads: &[[u32; 4]],
        params: &[OutOfPlaneParams],
    ) -> Result<Self, EngineError> {
        if params.len() != quads.len() {
            return Err(EngineError::OutOfPlaneParameterCountMismatch {
                provided: params.len(),
                expected: quads.len(),
            });
        }
        let mut term = Self::new();
        for (quad, param) in quads.iter().zip(params.iter()) {
            term.add_improper(
                quad[0],
                quad[1],
                quad[2],
                quad[3],
                param.k_j_per_rad2,
                param.chi0_rad,
            )?;
        }
        Ok(term)
    }

    /// Appends one improper with its parameters.
    ///
    /// The four atoms must be distinct. The force constant must be finite and
    /// non-negative. The equilibrium angle must be finite and within
    /// `[-pi/2, pi/2]`.
    #[allow(clippy::too_many_arguments)]
    pub fn add_improper(
        &mut self,
        i: u32,
        j: u32,
        k: u32,
        l: u32,
        k_j_per_rad2: f64,
        chi0_rad: f64,
    ) -> Result<(), EngineError> {
        let improper = self.quads.len();
        if i == j || j == k || k == l || i == k || i == l || j == l {
            return Err(EngineError::InvalidOutOfPlaneParameters { improper });
        }
        if !k_j_per_rad2.is_finite() || k_j_per_rad2 < 0.0 {
            return Err(EngineError::InvalidOutOfPlaneParameters { improper });
        }
        if !chi0_rad.is_finite()
            || chi0_rad < -core::f64::consts::FRAC_PI_2
            || chi0_rad > core::f64::consts::FRAC_PI_2
        {
            return Err(EngineError::InvalidOutOfPlaneParameters { improper });
        }
        self.quads.push([i, j, k, l]);
        self.k_j_per_rad2.push(k_j_per_rad2);
        self.chi0_rad.push(chi0_rad);
        let span = (i.max(j).max(k).max(l) as usize).saturating_add(1);
        if span > self.atom_count {
            self.atom_count = span;
        }
        Ok(())
    }

    /// Returns the number of impropers.
    pub fn improper_count(&self) -> usize {
        self.quads.len()
    }

    /// Returns the number of atoms the term covers.
    pub fn atom_count(&self) -> usize {
        self.atom_count
    }

    /// Returns the total out-of-plane energy in joules for a flat buffer.
    pub fn energy_j(&self, positions_m: &[f64]) -> Result<f64, EngineError> {
        validate_positions(positions_m, self.atom_count)?;
        let mut energy_j = 0.0;
        for (improper, quad) in self.quads.iter().enumerate() {
            let chi_rad = self.chi_rad(positions_m, *quad, improper)?;
            let delta_rad = chi_rad - self.chi0_rad[improper];
            energy_j += 0.5 * self.k_j_per_rad2[improper] * delta_rad * delta_rad;
        }
        Ok(energy_j)
    }

    /// Returns the energy gradient in newtons for a flat position buffer.
    pub fn gradient_j_per_m(&self, positions_m: &[f64]) -> Result<Vec<f64>, EngineError> {
        validate_positions(positions_m, self.atom_count)?;
        let mut gradient = vec![0.0; 3 * self.atom_count];
        for (improper, quad) in self.quads.iter().enumerate() {
            let [i, j, k, l] = *quad;
            let a = sub(
                atom_position_m(positions_m, i),
                atom_position_m(positions_m, j),
            );
            let b = sub(
                atom_position_m(positions_m, l),
                atom_position_m(positions_m, j),
            );
            let c = sub(
                atom_position_m(positions_m, k),
                atom_position_m(positions_m, j),
            );
            let n = cross(a, b);
            let n_len = norm(n);
            let c_len = norm(c);
            if n_len <= 0.0 || c_len <= 0.0 {
                return Err(EngineError::OutOfPlaneGeometryUndefined { improper });
            }
            let g = dot(c, n);
            let sin_chi = (g / (c_len * n_len)).clamp(-1.0, 1.0);
            let cos_chi = (1.0 - sin_chi * sin_chi).sqrt();
            if cos_chi < 1.0e-12 {
                return Err(EngineError::OutOfPlaneGeometryUndefined { improper });
            }
            let chi_rad = sin_chi.asin();
            let du_dchi = self.k_j_per_rad2[improper] * (chi_rad - self.chi0_rad[improper]);
            let inv_cn = 1.0 / (c_len * n_len);
            let inv_cn3 = inv_cn / (n_len * n_len);
            let inv_c3n = 1.0 / (c_len * c_len * c_len * n_len);
            let a_a = scale(
                sub(scale(cross(b, c), inv_cn), scale(cross(b, n), g * inv_cn3)),
                du_dchi / cos_chi,
            );
            let a_b = scale(
                sub(scale(cross(c, a), inv_cn), scale(cross(n, a), g * inv_cn3)),
                du_dchi / cos_chi,
            );
            let a_c = scale(
                sub(scale(n, inv_cn), scale(c, g * inv_c3n)),
                du_dchi / cos_chi,
            );
            let g_j = scale(
                [
                    a_a[0] + a_b[0] + a_c[0],
                    a_a[1] + a_b[1] + a_c[1],
                    a_a[2] + a_b[2] + a_c[2],
                ],
                -1.0,
            );
            add_assign(&mut gradient, i, a_a);
            add_assign(&mut gradient, j, g_j);
            add_assign(&mut gradient, k, a_c);
            add_assign(&mut gradient, l, a_b);
        }
        Ok(gradient)
    }

    /// Returns the force in newtons for a flat position buffer.
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

    fn chi_rad(
        &self,
        positions_m: &[f64],
        quad: [u32; 4],
        improper: usize,
    ) -> Result<f64, EngineError> {
        let [i, j, k, l] = quad;
        let a = sub(
            atom_position_m(positions_m, i),
            atom_position_m(positions_m, j),
        );
        let b = sub(
            atom_position_m(positions_m, l),
            atom_position_m(positions_m, j),
        );
        let c = sub(
            atom_position_m(positions_m, k),
            atom_position_m(positions_m, j),
        );
        let n = cross(a, b);
        let n_len = norm(n);
        let c_len = norm(c);
        if n_len <= 0.0 || c_len <= 0.0 {
            return Err(EngineError::OutOfPlaneGeometryUndefined { improper });
        }
        Ok((dot(c, n) / (c_len * n_len)).clamp(-1.0, 1.0).asin())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::{central_difference, Rng};

    fn pyramidal_quad(seed: u64) -> ([u32; 4], Vec<f64>) {
        let mut rng = Rng::new(seed);
        let r_m = 1.4e-10;
        let base_m = [
            [r_m, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.3 * r_m, 0.2 * r_m, 0.15 * r_m],
            [0.0, r_m, 0.0],
        ];
        let mut positions_m = Vec::with_capacity(12);
        for atom in base_m {
            positions_m.push(atom[0] + rng.symmetric(0.02 * r_m));
            positions_m.push(atom[1] + rng.symmetric(0.02 * r_m));
            positions_m.push(atom[2] + rng.symmetric(0.02 * r_m));
        }
        ([0, 1, 2, 3], positions_m)
    }

    #[test]
    fn an_empty_term_has_zero_energy_and_an_empty_gradient() {
        let term = OutOfPlaneTerm::new();
        assert_eq!(term.improper_count(), 0);
        assert_eq!(term.atom_count(), 0);
        assert_eq!(term.energy_j(&[]).expect("valid"), 0.0);
        assert!(term.gradient_j_per_m(&[]).expect("valid").is_empty());
    }

    #[test]
    fn a_planar_improper_at_equilibrium_has_zero_energy() {
        let term =
            OutOfPlaneTerm::from_quads(&[[0, 1, 2, 3]], &[OutOfPlaneParams::new(1.0e-18, 0.0)])
                .expect("valid term");
        let r_m = 1.4e-10;
        let positions_m = vec![
            r_m, 0.0, 0.0, // i
            0.0, 0.0, 0.0, // j
            r_m, r_m, 0.0, // k, in the plane
            0.0, r_m, 0.0, // l
        ];
        assert!(term.energy_j(&positions_m).expect("valid").abs() < 1.0e-30);
    }

    #[test]
    fn invalid_out_of_plane_parameters_are_rejected() {
        let mut term = OutOfPlaneTerm::new();
        assert!(matches!(
            term.add_improper(0, 1, 2, 3, -1.0, 0.0),
            Err(EngineError::InvalidOutOfPlaneParameters { .. })
        ));
        assert!(matches!(
            term.add_improper(0, 1, 2, 3, 1.0, 2.0),
            Err(EngineError::InvalidOutOfPlaneParameters { .. })
        ));
        assert!(matches!(
            term.add_improper(0, 1, 2, 0, 1.0, 0.0),
            Err(EngineError::InvalidOutOfPlaneParameters { .. })
        ));
    }

    #[test]
    fn a_parameter_count_mismatch_is_rejected() {
        assert!(matches!(
            OutOfPlaneTerm::from_quads(&[[0, 1, 2, 3]], &[]),
            Err(EngineError::OutOfPlaneParameterCountMismatch { .. })
        ));
    }

    #[test]
    fn the_analytic_gradient_matches_a_central_finite_difference() {
        let (quad, positions_m) = pyramidal_quad(0x0F);
        let term = OutOfPlaneTerm::from_quads(&[quad], &[OutOfPlaneParams::new(4.0e-19, 0.0)])
            .expect("valid term");
        let analytic = term.gradient_j_per_m(&positions_m).expect("valid");
        let step_m = 1.0e-14;
        let finite_difference =
            central_difference(|p| term.energy_j(p).expect("valid"), &positions_m, step_m);
        let mut max_error_n = 0.0_f64;
        let rtol = 1.0e-6;
        let atol_n = 1.0e-20;
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
        println!("out-of-plane max gradient error: {max_error_n:e} N");
    }
}
