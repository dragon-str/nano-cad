//! Metal-ligand coordination, modelled with a Morse bond.
//!
//! A transition metal binds a gas such as dioxygen by an electronic
//! interaction: sigma donation from the ligand and pi back-donation from the
//! metal. That interaction has no representation in a Lennard-Jones plus
//! point-charge force field. This term supplies it.
//!
//! Each pair holds a Morse potential
//! `U(r) = D * (1 - exp(-a * (r - r0)))^2`, with `D` the well depth in joules,
//! `r0` the equilibrium length in metres, and `a` the steepness in reciprocal
//! metres. The energy at `r0` is zero, and the energy at infinite separation is
//! `D`. The well depth is the binding energy of the contact.
//!
//! The term states no angle. Real coordination is directional. An angular term
//! is a later task. See ADR-0083.

use std::f64::consts::LN_2;

use crate::error::EngineError;

/// Morse parameters for one metal-ligand contact.
///
/// The potential is `U(r) = D * (1 - exp(-a * (r - r0)))^2`. The binding energy
/// is `D` in joules. The equilibrium length is `r0` in metres. The steepness is
/// `a` in reciprocal metres.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MetalLigandParams {
    /// The well depth in joules. The binding energy of the contact.
    pub well_depth_j: f64,
    /// The equilibrium metal-ligand length in metres.
    pub equilibrium_m: f64,
    /// The Morse steepness in reciprocal metres.
    pub steepness_per_m: f64,
}

impl MetalLigandParams {
    /// Builds parameters from a well depth, an equilibrium length and a
    /// steepness.
    pub const fn new(well_depth_j: f64, equilibrium_m: f64, steepness_per_m: f64) -> Self {
        Self {
            well_depth_j,
            equilibrium_m,
            steepness_per_m,
        }
    }

    /// Builds parameters from a well depth, an equilibrium length and a
    /// harmonic force constant.
    ///
    /// Near the minimum the Morse potential equals a harmonic spring of force
    /// constant `k = 2 * D * a^2`, so `a = sqrt(k / (2 * D))`. This lets a
    /// measured stretching force constant set the steepness. A zero or negative
    /// well depth gives a straight line (`a` is zero), which the term rejects.
    pub fn from_force_constant(well_depth_j: f64, equilibrium_m: f64, k_n_per_m: f64) -> Self {
        let steepness_per_m = (k_n_per_m / (2.0 * well_depth_j)).sqrt();
        Self::new(well_depth_j, equilibrium_m, steepness_per_m)
    }

    /// Builds parameters whose well depth is the magnitude of a binding free
    /// energy.
    ///
    /// The Morse well depth is a zero-kelvin potential energy. A measured
    /// binding free energy includes entropy and the choice of standard state.
    /// This constructor sets the well depth to `-delta_g_j` as a first
    /// approximation, and takes the steepness from a force constant. State the
    /// limitation where the result is used. A non-negative free energy gives a
    /// zero depth, which the term rejects.
    pub fn from_binding_energy_j(delta_g_j: f64, equilibrium_m: f64, k_n_per_m: f64) -> Self {
        Self::from_force_constant(-delta_g_j, equilibrium_m, k_n_per_m)
    }

    /// Returns the harmonic force constant at the minimum, in newtons per
    /// metre. It is `2 * D * a^2`.
    pub fn force_constant_n_per_m(&self) -> f64 {
        2.0 * self.well_depth_j * self.steepness_per_m * self.steepness_per_m
    }

    /// Returns the Morse energy in joules at a metal-ligand distance.
    ///
    /// The value at the equilibrium length is zero. The value at infinite
    /// separation is the well depth.
    pub fn energy_j(&self, distance_m: f64) -> f64 {
        let decay = (-self.steepness_per_m * (distance_m - self.equilibrium_m)).exp();
        let shifted = 1.0 - decay;
        self.well_depth_j * shifted * shifted
    }

    /// Returns the distance at which the Morse energy is zero on the repulsive
    /// side, in metres. It is `r0 - ln(2) / a`.
    pub fn crossing_m(&self) -> f64 {
        self.equilibrium_m - LN_2 / self.steepness_per_m
    }
}

/// A read-only view of one pair in a [`MetalLigandTerm`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MetalLigandInfo {
    /// Index of the metal atom.
    pub metal: u32,
    /// Index of the ligand donor atom.
    pub ligand: u32,
    /// The well depth in joules.
    pub well_depth_j: f64,
    /// The equilibrium length in metres.
    pub equilibrium_m: f64,
    /// The Morse steepness in reciprocal metres.
    pub steepness_per_m: f64,
}

/// A Morse coordination term over a set of metal-ligand pairs.
///
/// Energy is in joules. The gradient is the derivative of the energy with
/// respect to each Cartesian coordinate, in newtons.
#[derive(Clone, Debug, PartialEq)]
pub struct MetalLigandTerm {
    atom_count: usize,
    metals: Vec<u32>,
    ligands: Vec<u32>,
    params: Vec<MetalLigandParams>,
}

impl Default for MetalLigandTerm {
    fn default() -> Self {
        Self::new()
    }
}

impl MetalLigandTerm {
    /// Creates an empty term with no pairs and no atoms.
    pub const fn new() -> Self {
        Self {
            atom_count: 0,
            metals: Vec::new(),
            ligands: Vec::new(),
            params: Vec::new(),
        }
    }

    /// Appends one metal-ligand pair.
    ///
    /// The metal and the ligand must differ. The well depth must be finite and
    /// positive. The equilibrium length must be finite and positive. The
    /// steepness must be finite and positive. The term grows its atom count to
    /// cover the larger index.
    pub fn add_pair(
        &mut self,
        metal: u32,
        ligand: u32,
        params: MetalLigandParams,
    ) -> Result<(), EngineError> {
        let pair = self.metals.len();
        if metal == ligand {
            return Err(EngineError::InvalidMetalLigandParameters { pair });
        }
        if !params.well_depth_j.is_finite() || params.well_depth_j <= 0.0 {
            return Err(EngineError::InvalidMetalLigandParameters { pair });
        }
        if !params.equilibrium_m.is_finite() || params.equilibrium_m <= 0.0 {
            return Err(EngineError::InvalidMetalLigandParameters { pair });
        }
        if !params.steepness_per_m.is_finite() || params.steepness_per_m <= 0.0 {
            return Err(EngineError::InvalidMetalLigandParameters { pair });
        }
        self.metals.push(metal);
        self.ligands.push(ligand);
        self.params.push(params);
        let span = (metal.max(ligand) as usize).saturating_add(1);
        if span > self.atom_count {
            self.atom_count = span;
        }
        Ok(())
    }

    /// Returns the number of metal-ligand pairs.
    pub fn pair_count(&self) -> usize {
        self.metals.len()
    }

    /// Returns the number of atoms the term covers.
    pub fn atom_count(&self) -> usize {
        self.atom_count
    }

    /// Returns a read-only view of one pair, or `None` when the index is out
    /// of range.
    pub fn pair_info(&self, pair: usize) -> Option<MetalLigandInfo> {
        Some(MetalLigandInfo {
            metal: *self.metals.get(pair)?,
            ligand: *self.ligands.get(pair)?,
            well_depth_j: self.params.get(pair)?.well_depth_j,
            equilibrium_m: self.params.get(pair)?.equilibrium_m,
            steepness_per_m: self.params.get(pair)?.steepness_per_m,
        })
    }

    /// Returns the coordination energy in joules for a flat position buffer.
    ///
    /// A pair at its equilibrium length contributes zero. A pair at infinite
    /// separation contributes the well depth. A separated pair cannot be
    /// evaluated from positions, so the value is always below the well depth.
    pub fn energy_j(&self, positions_m: &[f64]) -> Result<f64, EngineError> {
        self.validate_positions(positions_m)?;
        let mut energy_j = 0.0;
        for pair in 0..self.pair_count() {
            let distance_m = self.distance_m(positions_m, pair);
            energy_j += self.params[pair].energy_j(distance_m);
        }
        Ok(energy_j)
    }

    /// Returns the energy gradient in newtons for a flat position buffer.
    ///
    /// Entry `3 * atom + axis` is `dU/dx`. A pair with coincident atoms has no
    /// defined direction and returns an error.
    pub fn gradient_j_per_m(&self, positions_m: &[f64]) -> Result<Vec<f64>, EngineError> {
        self.validate_positions(positions_m)?;
        let mut gradient = vec![0.0; 3 * self.atom_count];
        for pair in 0..self.pair_count() {
            let metal = self.metals[pair];
            let ligand = self.ligands[pair];
            let base_metal = metal as usize * 3;
            let base_ligand = ligand as usize * 3;
            let dx = positions_m[base_metal] - positions_m[base_ligand];
            let dy = positions_m[base_metal + 1] - positions_m[base_ligand + 1];
            let dz = positions_m[base_metal + 2] - positions_m[base_ligand + 2];
            let distance_m = (dx * dx + dy * dy + dz * dz).sqrt();
            if distance_m <= 0.0 {
                return Err(EngineError::CoincidentMetalLigandAtoms { pair });
            }
            let scale = self.radial_gradient_n(pair, distance_m) / distance_m;
            gradient[base_metal] += scale * dx;
            gradient[base_metal + 1] += scale * dy;
            gradient[base_metal + 2] += scale * dz;
            gradient[base_ligand] -= scale * dx;
            gradient[base_ligand + 1] -= scale * dy;
            gradient[base_ligand + 2] -= scale * dz;
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

    /// Returns the radial derivative `dU/dr` in newtons for one pair.
    ///
    /// The value is zero at the equilibrium length, negative inside it (the
    /// energy falls as the atoms separate) and positive outside it.
    pub fn radial_gradient_n(&self, pair: usize, distance_m: f64) -> f64 {
        let params = self.params[pair];
        let exponent = -params.steepness_per_m * (distance_m - params.equilibrium_m);
        let decay = exponent.exp();
        2.0 * params.well_depth_j * params.steepness_per_m * decay * (1.0 - decay)
    }

    /// Returns the radial force in newtons for one pair at a distance.
    ///
    /// The force is `-dU/dr`. A positive value pushes the atoms apart. The
    /// force is zero at the equilibrium length, positive inside it and negative
    /// outside it.
    pub fn radial_force_n(&self, pair: usize, distance_m: f64) -> f64 {
        -self.radial_gradient_n(pair, distance_m)
    }

    /// Returns the coordination distance of each pair in metres.
    pub fn distances_m(&self, positions_m: &[f64]) -> Result<Vec<f64>, EngineError> {
        self.validate_positions(positions_m)?;
        Ok((0..self.pair_count())
            .map(|pair| self.distance_m(positions_m, pair))
            .collect())
    }

    /// Returns the sum of the well depths in joules. It is the binding energy
    /// of a fully separated contact.
    pub fn dissociation_energy_j(&self) -> f64 {
        self.params.iter().map(|params| params.well_depth_j).sum()
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

    fn distance_m(&self, positions_m: &[f64], pair: usize) -> f64 {
        let base_metal = self.metals[pair] as usize * 3;
        let base_ligand = self.ligands[pair] as usize * 3;
        let dx = positions_m[base_metal] - positions_m[base_ligand];
        let dy = positions_m[base_metal + 1] - positions_m[base_ligand + 1];
        let dz = positions_m[base_metal + 2] - positions_m[base_ligand + 2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A plausible iron-oxygen contact. Depth about 20 kT at body temperature.
    fn iron_oxygen() -> MetalLigandParams {
        MetalLigandParams::new(1.5e-20, 2.0e-10, 2.0e10)
    }

    fn diatomic(separation_m: f64) -> Vec<f64> {
        vec![0.0, 0.0, 0.0, separation_m, 0.0, 0.0]
    }

    #[test]
    fn an_empty_term_has_zero_energy_and_an_empty_gradient() {
        let term = MetalLigandTerm::new();
        assert_eq!(term.pair_count(), 0);
        assert_eq!(term.atom_count(), 0);
        assert_eq!(term.energy_j(&[]).expect("valid"), 0.0);
        assert!(term.gradient_j_per_m(&[]).expect("valid").is_empty());
        assert_eq!(term.dissociation_energy_j(), 0.0);
    }

    #[test]
    fn a_pair_at_its_equilibrium_length_has_zero_energy_and_zero_force() {
        let params = iron_oxygen();
        let mut term = MetalLigandTerm::new();
        term.add_pair(0, 1, params).expect("valid pair");
        let positions_m = diatomic(params.equilibrium_m);
        assert!(term.energy_j(&positions_m).expect("valid").abs() < 1.0e-30);
        for force in term.forces_n(&positions_m).expect("valid") {
            assert!(force.abs() < 1.0e-24);
        }
        assert!(term.radial_force_n(0, params.equilibrium_m).abs() < 1.0e-24);
    }

    #[test]
    fn the_energy_approaches_the_well_depth_at_long_range() {
        let params = iron_oxygen();
        let mut term = MetalLigandTerm::new();
        term.add_pair(0, 1, params).expect("valid pair");
        let far_m = params.equilibrium_m + 20.0 / params.steepness_per_m;
        let energy_j = term.energy_j(&diatomic(far_m)).expect("valid");
        assert!(
            (energy_j - params.well_depth_j).abs() / params.well_depth_j < 1.0e-8,
            "the long-range energy {energy_j} must approach the well depth {}",
            params.well_depth_j
        );
        assert!((term.dissociation_energy_j() - params.well_depth_j).abs() < 1.0e-30);
    }

    #[test]
    fn the_energy_is_repulsive_inside_the_equilibrium_length() {
        let params = iron_oxygen();
        let mut term = MetalLigandTerm::new();
        term.add_pair(0, 1, params).expect("valid pair");
        let inside_m = 0.7 * params.equilibrium_m;
        let energy_j = term.energy_j(&diatomic(inside_m)).expect("valid");
        assert!(
            energy_j > params.well_depth_j,
            "the energy {energy_j} inside the well must exceed the depth"
        );
        assert!(
            term.radial_force_n(0, inside_m) > 0.0,
            "the force must repel"
        );
        assert!(term.radial_force_n(0, 1.3 * params.equilibrium_m) < 0.0);
    }

    #[test]
    fn a_force_constant_sets_the_steepness() {
        let k_n_per_m = 5.0e2;
        let params = MetalLigandParams::from_force_constant(1.2e-20, 2.0e-10, k_n_per_m);
        let restored = params.force_constant_n_per_m();
        assert!(
            (restored - k_n_per_m).abs() / k_n_per_m < 1.0e-12,
            "the force constant must round trip: {restored}"
        );
        assert!(params.crossing_m() < params.equilibrium_m);
    }

    #[test]
    fn invalid_parameters_are_rejected() {
        let mut term = MetalLigandTerm::new();
        assert!(matches!(
            term.add_pair(0, 0, iron_oxygen()),
            Err(EngineError::InvalidMetalLigandParameters { .. })
        ));
        assert!(matches!(
            term.add_pair(0, 1, MetalLigandParams::new(0.0, 2.0e-10, 2.0e10)),
            Err(EngineError::InvalidMetalLigandParameters { .. })
        ));
        assert!(matches!(
            term.add_pair(0, 1, MetalLigandParams::new(1.5e-20, 0.0, 2.0e10)),
            Err(EngineError::InvalidMetalLigandParameters { .. })
        ));
        assert!(matches!(
            term.add_pair(0, 1, MetalLigandParams::new(1.5e-20, 2.0e-10, 0.0)),
            Err(EngineError::InvalidMetalLigandParameters { .. })
        ));
    }

    #[test]
    fn a_position_buffer_of_the_wrong_size_is_rejected() {
        let mut term = MetalLigandTerm::new();
        term.add_pair(0, 1, iron_oxygen()).expect("valid pair");
        assert!(matches!(
            term.energy_j(&[0.0, 0.0, 0.0]),
            Err(EngineError::PositionBufferSizeMismatch { .. })
        ));
    }

    #[test]
    fn a_coincident_pair_is_rejected() {
        let mut term = MetalLigandTerm::new();
        term.add_pair(0, 1, iron_oxygen()).expect("valid pair");
        assert!(matches!(
            term.gradient_j_per_m(&diatomic(0.0)),
            Err(EngineError::CoincidentMetalLigandAtoms { .. })
        ));
    }

    #[test]
    fn pair_info_and_distances_match_the_geometry() {
        let params = iron_oxygen();
        let mut term = MetalLigandTerm::new();
        term.add_pair(1, 0, params).expect("valid pair");
        let info = term.pair_info(0).expect("pair 0");
        assert_eq!(info.metal, 1);
        assert_eq!(info.ligand, 0);
        assert_eq!(info.well_depth_j, params.well_depth_j);
        assert_eq!(term.pair_info(9), None);
        assert_eq!(term.atom_count(), 2);

        let separation_m = 1.2 * params.equilibrium_m;
        let distances = term.distances_m(&diatomic(separation_m)).expect("valid");
        assert!((distances[0] - separation_m).abs() < 1.0e-24);
    }

    #[test]
    fn the_analytic_gradient_matches_a_central_finite_difference() {
        let params = iron_oxygen();
        let mut term = MetalLigandTerm::new();
        term.add_pair(0, 1, params).expect("first pair");
        term.add_pair(1, 2, iron_oxygen()).expect("second pair");
        let positions_m = vec![
            0.0, 0.0, 0.0, 1.31e-10, 0.4e-10, -0.2e-10, 2.62e-10, 0.0, 0.3e-10,
        ];
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
            "maximum gradient error {max_error_n}"
        );
    }

    #[test]
    fn the_force_is_the_negative_gradient() {
        let mut term = MetalLigandTerm::new();
        term.add_pair(0, 1, iron_oxygen()).expect("valid pair");
        let positions_m = diatomic(1.4e-10);
        let gradient = term.gradient_j_per_m(&positions_m).expect("valid");
        let forces = term.forces_n(&positions_m).expect("valid");
        for (force, grad) in forces.iter().zip(gradient.iter()) {
            assert!((force + grad).abs() < 1.0e-30);
        }
    }
}
