use crate::error::EngineError;
use crate::geometry::{
    add, add_assign, atom_position_m, cross, dot, norm, scale, sub, validate_positions,
};

/// Periodic torsion parameters for one dihedral.
///
/// The potential is the three-term cosine series
/// `U(phi) = 0.5 V1 (1 + cos(phi)) + 0.5 V2 (1 - cos(2 phi)) + 0.5 V3 (1 + cos(3 phi))`.
/// The amplitudes are in joules. A zero amplitude disables its term.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TorsionParams {
    /// Onefold barrier amplitude in joules.
    pub v1_j: f64,
    /// Twofold barrier amplitude in joules.
    pub v2_j: f64,
    /// Threefold barrier amplitude in joules.
    pub v3_j: f64,
}

impl TorsionParams {
    /// Builds parameters for one dihedral.
    pub const fn new(v1_j: f64, v2_j: f64, v3_j: f64) -> Self {
        Self { v1_j, v2_j, v3_j }
    }
}

/// A periodic torsion term over a set of atom quads.
///
/// Each quad is `[i, j, k, l]`, where the dihedral turns about the `j-k` bond.
/// Energy is in joules. The gradient is the derivative of the energy with
/// respect to each Cartesian coordinate, in newtons.
#[derive(Clone, Debug, PartialEq)]
pub struct TorsionTerm {
    atom_count: usize,
    quads: Vec<[u32; 4]>,
    v1_j: Vec<f64>,
    v2_j: Vec<f64>,
    v3_j: Vec<f64>,
}

impl Default for TorsionTerm {
    fn default() -> Self {
        Self::new()
    }
}

impl TorsionTerm {
    /// Creates an empty term with no torsions and no atoms.
    pub const fn new() -> Self {
        Self {
            atom_count: 0,
            quads: Vec::new(),
            v1_j: Vec::new(),
            v2_j: Vec::new(),
            v3_j: Vec::new(),
        }
    }

    /// Builds a term from a list of quads and one parameter set per quad.
    pub fn from_quads(quads: &[[u32; 4]], params: &[TorsionParams]) -> Result<Self, EngineError> {
        if params.len() != quads.len() {
            return Err(EngineError::TorsionParameterCountMismatch {
                provided: params.len(),
                expected: quads.len(),
            });
        }
        let mut term = Self::new();
        for (quad, param) in quads.iter().zip(params.iter()) {
            term.add_torsion(
                quad[0], quad[1], quad[2], quad[3], param.v1_j, param.v2_j, param.v3_j,
            )?;
        }
        Ok(term)
    }

    /// Appends one torsion with its parameters.
    ///
    /// The four atoms must be distinct. Every amplitude must be finite.
    #[allow(clippy::too_many_arguments)]
    pub fn add_torsion(
        &mut self,
        i: u32,
        j: u32,
        k: u32,
        l: u32,
        v1_j: f64,
        v2_j: f64,
        v3_j: f64,
    ) -> Result<(), EngineError> {
        let torsion = self.quads.len();
        if i == j || j == k || k == l || i == k || i == l || j == l {
            return Err(EngineError::InvalidTorsionParameters { torsion });
        }
        if !v1_j.is_finite() || !v2_j.is_finite() || !v3_j.is_finite() {
            return Err(EngineError::InvalidTorsionParameters { torsion });
        }
        self.quads.push([i, j, k, l]);
        self.v1_j.push(v1_j);
        self.v2_j.push(v2_j);
        self.v3_j.push(v3_j);
        let span = (i.max(j).max(k).max(l) as usize).saturating_add(1);
        if span > self.atom_count {
            self.atom_count = span;
        }
        Ok(())
    }

    /// Returns the number of torsions.
    pub fn torsion_count(&self) -> usize {
        self.quads.len()
    }

    /// Returns the number of atoms the term covers.
    pub fn atom_count(&self) -> usize {
        self.atom_count
    }

    /// Returns the total torsion energy in joules for a flat position buffer.
    pub fn energy_j(&self, positions_m: &[f64]) -> Result<f64, EngineError> {
        validate_positions(positions_m, self.atom_count)?;
        let mut energy_j = 0.0;
        for (torsion, quad) in self.quads.iter().enumerate() {
            let phi_rad = self.phi_rad(positions_m, *quad, torsion)?;
            energy_j += 0.5 * self.v1_j[torsion] * (1.0 + phi_rad.cos());
            energy_j += 0.5 * self.v2_j[torsion] * (1.0 - (2.0 * phi_rad).cos());
            energy_j += 0.5 * self.v3_j[torsion] * (1.0 + (3.0 * phi_rad).cos());
        }
        Ok(energy_j)
    }

    /// Returns the energy gradient in newtons for a flat position buffer.
    pub fn gradient_j_per_m(&self, positions_m: &[f64]) -> Result<Vec<f64>, EngineError> {
        validate_positions(positions_m, self.atom_count)?;
        let mut gradient = vec![0.0; 3 * self.atom_count];
        for (torsion, quad) in self.quads.iter().enumerate() {
            let [i, j, k, l] = *quad;
            let b1 = sub(
                atom_position_m(positions_m, j),
                atom_position_m(positions_m, i),
            );
            let b2 = sub(
                atom_position_m(positions_m, k),
                atom_position_m(positions_m, j),
            );
            let b3 = sub(
                atom_position_m(positions_m, l),
                atom_position_m(positions_m, k),
            );
            let n1 = cross(b1, b2);
            let n2 = cross(b2, b3);
            let b2_sq = dot(b2, b2);
            let n1_sq = dot(n1, n1);
            let n2_sq = dot(n2, n2);
            if b2_sq <= 0.0 || n1_sq <= 0.0 || n2_sq <= 0.0 {
                return Err(EngineError::TorsionGeometryUndefined { torsion });
            }
            let b2_len = b2_sq.sqrt();
            let d1 = b2_len / n1_sq;
            let d2 = b2_len / n2_sq;
            let c1 = dot(b1, b2) / b2_sq;
            let c2 = dot(b3, b2) / b2_sq;
            let x = dot(n1, n2);
            let y = dot(cross(n1, n2), b2) / b2_len;
            let phi_rad = y.atan2(x);
            let v1 = self.v1_j[torsion];
            let v2 = self.v2_j[torsion];
            let v3 = self.v3_j[torsion];
            let du_dphi = -0.5 * v1 * phi_rad.sin() + v2 * (2.0 * phi_rad).sin()
                - 1.5 * v3 * (3.0 * phi_rad).sin();
            let g_i = scale(n1, -d1 * du_dphi);
            let g_j = add(
                scale(n1, (1.0 + c1) * d1 * du_dphi),
                scale(n2, c2 * d2 * du_dphi),
            );
            let g_k = add(
                scale(n1, -c1 * d1 * du_dphi),
                scale(n2, -(1.0 + c2) * d2 * du_dphi),
            );
            let g_l = scale(n2, d2 * du_dphi);
            add_assign(&mut gradient, i, g_i);
            add_assign(&mut gradient, j, g_j);
            add_assign(&mut gradient, k, g_k);
            add_assign(&mut gradient, l, g_l);
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

    fn phi_rad(
        &self,
        positions_m: &[f64],
        quad: [u32; 4],
        torsion: usize,
    ) -> Result<f64, EngineError> {
        let [i, j, k, l] = quad;
        let b1 = sub(
            atom_position_m(positions_m, j),
            atom_position_m(positions_m, i),
        );
        let b2 = sub(
            atom_position_m(positions_m, k),
            atom_position_m(positions_m, j),
        );
        let b3 = sub(
            atom_position_m(positions_m, l),
            atom_position_m(positions_m, k),
        );
        let n1 = cross(b1, b2);
        let n2 = cross(b2, b3);
        let b2_len = norm(b2);
        if b2_len <= 0.0 || norm(n1) <= 0.0 || norm(n2) <= 0.0 {
            return Err(EngineError::TorsionGeometryUndefined { torsion });
        }
        let x = dot(n1, n2);
        let y = dot(cross(n1, n2), b2) / b2_len;
        Ok(y.atan2(x))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::{central_difference, Rng};

    fn butane_quad(seed: u64) -> ([u32; 4], Vec<f64>) {
        let mut rng = Rng::new(seed);
        let r_m = 1.5e-10;
        let base = [
            [-1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.5, 1.0, 0.5],
        ];
        let mut positions_m = Vec::with_capacity(12);
        for atom in base {
            positions_m.push(atom[0] * r_m + rng.symmetric(0.02 * r_m));
            positions_m.push(atom[1] * r_m + rng.symmetric(0.02 * r_m));
            positions_m.push(atom[2] * r_m + rng.symmetric(0.02 * r_m));
        }
        ([0, 1, 2, 3], positions_m)
    }

    #[test]
    fn an_empty_term_has_zero_energy_and_an_empty_gradient() {
        let term = TorsionTerm::new();
        assert_eq!(term.torsion_count(), 0);
        assert_eq!(term.atom_count(), 0);
        assert_eq!(term.energy_j(&[]).expect("valid"), 0.0);
        assert!(term.gradient_j_per_m(&[]).expect("valid").is_empty());
    }

    #[test]
    fn invalid_torsion_parameters_are_rejected() {
        let mut term = TorsionTerm::new();
        assert!(matches!(
            term.add_torsion(0, 1, 2, 3, f64::NAN, 0.0, 0.0),
            Err(EngineError::InvalidTorsionParameters { .. })
        ));
        assert!(matches!(
            term.add_torsion(0, 1, 2, 0, 1.0, 0.0, 0.0),
            Err(EngineError::InvalidTorsionParameters { .. })
        ));
    }

    #[test]
    fn a_parameter_count_mismatch_is_rejected() {
        assert!(matches!(
            TorsionTerm::from_quads(&[[0, 1, 2, 3]], &[]),
            Err(EngineError::TorsionParameterCountMismatch { .. })
        ));
    }

    #[test]
    fn the_analytic_gradient_matches_a_central_finite_difference() {
        let (quad, positions_m) = butane_quad(0x70);
        let term =
            TorsionTerm::from_quads(&[quad], &[TorsionParams::new(2.0e-20, -3.0e-20, 5.0e-20)])
                .expect("valid term");
        let analytic = term.gradient_j_per_m(&positions_m).expect("valid");
        let step_m = 1.0e-14;
        let finite_difference =
            central_difference(|p| term.energy_j(p).expect("valid"), &positions_m, step_m);
        let mut max_error_n = 0.0_f64;
        let rtol = 1.0e-6;
        let atol_n = 1.0e-21;
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
        println!("torsion max gradient error: {max_error_n:e} N");
    }
}
