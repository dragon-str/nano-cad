use crate::error::EngineError;
use crate::geometry::{
    add, add_assign, atom_position_m, dot, norm, scale, sub, validate_positions,
};

/// Harmonic angle-bend parameters for one angle.
///
/// The potential is `U(theta) = 0.5 * k_j_per_rad2 * (theta - theta0_rad)^2`.
/// The force constant is in joules per radian squared. The equilibrium angle is
/// in radians.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngleBendParams {
    /// Force constant in joules per radian squared.
    pub k_j_per_rad2: f64,
    /// Equilibrium angle in radians.
    pub theta0_rad: f64,
}

impl AngleBendParams {
    /// Builds parameters for one angle.
    pub const fn new(k_j_per_rad2: f64, theta0_rad: f64) -> Self {
        Self {
            k_j_per_rad2,
            theta0_rad,
        }
    }
}

/// A harmonic angle-bend term over a set of atom triples.
///
/// Each triple is `[i, j, k]`, where `j` is the central atom. Energy is in
/// joules. The gradient is the derivative of the energy with respect to each
/// Cartesian coordinate, in newtons.
#[derive(Clone, Debug, PartialEq)]
pub struct AngleBendTerm {
    atom_count: usize,
    triples: Vec<[u32; 3]>,
    k_j_per_rad2: Vec<f64>,
    theta0_rad: Vec<f64>,
}

impl Default for AngleBendTerm {
    fn default() -> Self {
        Self::new()
    }
}

impl AngleBendTerm {
    /// Creates an empty term with no angles and no atoms.
    pub const fn new() -> Self {
        Self {
            atom_count: 0,
            triples: Vec::new(),
            k_j_per_rad2: Vec::new(),
            theta0_rad: Vec::new(),
        }
    }

    /// Builds a term from a list of triples and one parameter set per triple.
    pub fn from_triples(
        triples: &[[u32; 3]],
        params: &[AngleBendParams],
    ) -> Result<Self, EngineError> {
        if params.len() != triples.len() {
            return Err(EngineError::AngleParameterCountMismatch {
                provided: params.len(),
                expected: triples.len(),
            });
        }
        let mut term = Self::new();
        for (triple, param) in triples.iter().zip(params.iter()) {
            term.add_angle(
                triple[0],
                triple[1],
                triple[2],
                param.k_j_per_rad2,
                param.theta0_rad,
            )?;
        }
        Ok(term)
    }

    /// Appends one angle with its parameters.
    ///
    /// The three atoms must be distinct. The force constant must be finite and
    /// non-negative. The equilibrium angle must be finite and within
    /// `[0, pi]`.
    pub fn add_angle(
        &mut self,
        i: u32,
        j: u32,
        k: u32,
        k_j_per_rad2: f64,
        theta0_rad: f64,
    ) -> Result<(), EngineError> {
        let angle = self.triples.len();
        if i == j || j == k || i == k {
            return Err(EngineError::InvalidAngleParameters { angle });
        }
        if !k_j_per_rad2.is_finite() || k_j_per_rad2 < 0.0 {
            return Err(EngineError::InvalidAngleParameters { angle });
        }
        if !theta0_rad.is_finite() || !(0.0..=core::f64::consts::PI).contains(&theta0_rad) {
            return Err(EngineError::InvalidAngleParameters { angle });
        }
        self.triples.push([i, j, k]);
        self.k_j_per_rad2.push(k_j_per_rad2);
        self.theta0_rad.push(theta0_rad);
        let span = (i.max(j).max(k) as usize).saturating_add(1);
        if span > self.atom_count {
            self.atom_count = span;
        }
        Ok(())
    }

    /// Returns the number of angles.
    pub fn angle_count(&self) -> usize {
        self.triples.len()
    }

    /// Returns the number of atoms the term covers.
    pub fn atom_count(&self) -> usize {
        self.atom_count
    }

    /// Returns the total angle energy in joules for a flat position buffer.
    pub fn energy_j(&self, positions_m: &[f64]) -> Result<f64, EngineError> {
        validate_positions(positions_m, self.atom_count)?;
        let mut energy_j = 0.0;
        for (angle, triple) in self.triples.iter().enumerate() {
            let theta_rad = self.theta_rad(positions_m, *triple, angle)?;
            let delta_rad = theta_rad - self.theta0_rad[angle];
            energy_j += 0.5 * self.k_j_per_rad2[angle] * delta_rad * delta_rad;
        }
        Ok(energy_j)
    }

    /// Returns the energy gradient in newtons for a flat position buffer.
    pub fn gradient_j_per_m(&self, positions_m: &[f64]) -> Result<Vec<f64>, EngineError> {
        validate_positions(positions_m, self.atom_count)?;
        let mut gradient = vec![0.0; 3 * self.atom_count];
        for (angle, triple) in self.triples.iter().enumerate() {
            let [i, j, k] = *triple;
            let u = sub(
                atom_position_m(positions_m, i),
                atom_position_m(positions_m, j),
            );
            let w = sub(
                atom_position_m(positions_m, k),
                atom_position_m(positions_m, j),
            );
            let a = norm(u);
            let b = norm(w);
            if a <= 0.0 || b <= 0.0 {
                return Err(EngineError::AngleGeometryUndefined { angle });
            }
            let cos_theta = (dot(u, w) / (a * b)).clamp(-1.0, 1.0);
            let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
            if sin_theta < 1.0e-12 {
                return Err(EngineError::AngleGeometryUndefined { angle });
            }
            let theta_rad = cos_theta.acos();
            let du_dtheta = self.k_j_per_rad2[angle] * (theta_rad - self.theta0_rad[angle]);
            let e_u = scale(u, 1.0 / a);
            let e_w = scale(w, 1.0 / b);
            let g_i = scale(
                sub(e_w, scale(e_u, cos_theta)),
                -du_dtheta / (a * sin_theta),
            );
            let g_k = scale(
                sub(e_u, scale(e_w, cos_theta)),
                -du_dtheta / (b * sin_theta),
            );
            let g_j = scale(add(g_i, g_k), -1.0);
            add_assign(&mut gradient, i, g_i);
            add_assign(&mut gradient, j, g_j);
            add_assign(&mut gradient, k, g_k);
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

    fn theta_rad(
        &self,
        positions_m: &[f64],
        triple: [u32; 3],
        angle: usize,
    ) -> Result<f64, EngineError> {
        let [i, j, k] = triple;
        let u = sub(
            atom_position_m(positions_m, i),
            atom_position_m(positions_m, j),
        );
        let w = sub(
            atom_position_m(positions_m, k),
            atom_position_m(positions_m, j),
        );
        let a = norm(u);
        let b = norm(w);
        if a <= 0.0 || b <= 0.0 {
            return Err(EngineError::AngleGeometryUndefined { angle });
        }
        Ok((dot(u, w) / (a * b)).clamp(-1.0, 1.0).acos())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::{central_difference, Rng};

    fn bent_triple(seed: u64) -> ([u32; 3], Vec<f64>) {
        let mut rng = Rng::new(seed);
        let r0_m = 1.5e-10;
        let positions_m = vec![0.0, r0_m * 0.5, 0.0, 0.0, 0.0, 0.0, r0_m, -r0_m * 0.3, 0.0];
        let mut jittered = positions_m;
        for value in &mut jittered {
            *value += rng.symmetric(0.02 * r0_m);
        }
        ([0, 1, 2], jittered)
    }

    #[test]
    fn an_empty_term_has_zero_energy_and_an_empty_gradient() {
        let term = AngleBendTerm::new();
        assert_eq!(term.angle_count(), 0);
        assert_eq!(term.atom_count(), 0);
        assert_eq!(term.energy_j(&[]).expect("valid"), 0.0);
        assert!(term.gradient_j_per_m(&[]).expect("valid").is_empty());
    }

    #[test]
    fn a_straight_angle_at_equilibrium_has_zero_energy() {
        let term = AngleBendTerm::from_triples(
            &[[0, 1, 2]],
            &[AngleBendParams::new(1.0e-18, core::f64::consts::PI)],
        )
        .expect("valid term");
        let positions_m = vec![0.0, 0.0, 0.0, 1.0e-10, 0.0, 0.0, 2.0e-10, 0.0, 0.0];
        assert!(term.energy_j(&positions_m).expect("valid").abs() < 1.0e-30);
    }

    #[test]
    fn invalid_angle_parameters_are_rejected() {
        let mut term = AngleBendTerm::new();
        assert!(matches!(
            term.add_angle(0, 1, 2, -1.0, 1.0),
            Err(EngineError::InvalidAngleParameters { .. })
        ));
        assert!(matches!(
            term.add_angle(0, 1, 2, 1.0, 4.0),
            Err(EngineError::InvalidAngleParameters { .. })
        ));
        assert!(matches!(
            term.add_angle(0, 1, 1, 1.0, 1.0),
            Err(EngineError::InvalidAngleParameters { .. })
        ));
    }

    #[test]
    fn a_parameter_count_mismatch_is_rejected() {
        assert!(matches!(
            AngleBendTerm::from_triples(&[[0, 1, 2]], &[]),
            Err(EngineError::AngleParameterCountMismatch { .. })
        ));
    }

    #[test]
    fn the_analytic_gradient_matches_a_central_finite_difference() {
        let (triple, positions_m) = bent_triple(0xA1);
        let term = AngleBendTerm::from_triples(&[triple], &[AngleBendParams::new(5.0e-19, 1.9)])
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
        println!("angle-bend max gradient error: {max_error_n:e} N");
    }
}
