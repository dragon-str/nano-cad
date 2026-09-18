use std::collections::HashSet;

use crate::angle_bend::AngleBendTerm;
use crate::bond_stretch::{BondInfo, BondStretchTerm};
use crate::electrostatic::ElectrostaticTerm;
use crate::error::EngineError;
use crate::geometry::validate_positions;
use crate::lennard_jones::LennardJonesTerm;
use crate::neighbor::VerletList;
use crate::nonbonded::PeriodicBox;
use crate::out_of_plane::OutOfPlaneTerm;
use crate::pair_kernel::{
    KernelElectrostatic, KernelLennardJones, KernelVanDerWaals, NonbondedKernel,
};
use crate::torsion::TorsionTerm;
use crate::van_der_waals::VanDerWaalsTerm;

const DEFAULT_SKIN_M: f64 = 2.0e-10;

/// The pair count in one reduction block.
///
/// The block partition is fixed by the pair count, not by the thread count, so
/// the reduction order stays bit-identical for any number of threads. The
/// serial path uses the same blocks, so the parallel result is also
/// bit-identical to the serial result.
const PARALLEL_BLOCK_PAIRS: usize = 2048;

/// The Boltzmann constant in joules per kelvin.
pub const BOLTZMANN_J_PER_K: f64 = 1.380_649e-23;

/// The total potential of a molecular system.
///
/// A `System` owns the bonded terms and the non-bonded terms. It evaluates the
/// total energy and the total gradient. The non-bonded terms read a
/// [`VerletList`] pair list, not an all-pairs scan. The list rebuilds when any
/// atom moves more than half the skin from the last build.
///
/// Units are SI. Positions are in metres, energy is in joules, and the gradient
/// is in newtons.
///
/// The neighbor list uses the same periodic box as the non-bonded terms, so
/// the list and the energy are consistent. The box is non-periodic by default.
#[derive(Clone, Debug, PartialEq)]
pub struct System {
    atom_count: usize,
    masses_kg: Vec<f64>,
    bond_stretch: Option<BondStretchTerm>,
    angle_bend: Option<AngleBendTerm>,
    torsion: Option<TorsionTerm>,
    out_of_plane: Option<OutOfPlaneTerm>,
    van_der_waals: Option<VanDerWaalsTerm>,
    lennard_jones: Option<LennardJonesTerm>,
    electrostatic: Option<ElectrostaticTerm>,
    exclusions: HashSet<u64>,
    skin_m: f64,
    periodic_box: PeriodicBox,
    verlet: Option<VerletList>,
    nonbonded_pairs: Vec<u32>,
}

impl System {
    /// Creates an empty system for an atom count.
    ///
    /// No term is set. No mass is set. The default skin is 2e-10 m.
    pub fn new(atom_count: usize) -> Self {
        Self {
            atom_count,
            masses_kg: Vec::new(),
            bond_stretch: None,
            angle_bend: None,
            torsion: None,
            out_of_plane: None,
            van_der_waals: None,
            lennard_jones: None,
            electrostatic: None,
            exclusions: HashSet::new(),
            skin_m: DEFAULT_SKIN_M,
            periodic_box: PeriodicBox::non_periodic(),
            verlet: None,
            nonbonded_pairs: Vec::new(),
        }
    }

    /// Creates a system with one mass per atom, in kilograms.
    pub fn with_masses_kg(atom_count: usize, masses_kg: &[f64]) -> Result<Self, EngineError> {
        let mut system = Self::new(atom_count);
        system.set_masses_kg(masses_kg)?;
        Ok(system)
    }

    /// Returns the number of atoms.
    pub fn atom_count(&self) -> usize {
        self.atom_count
    }

    /// Replaces the mass list, in kilograms. The length must match the atom
    /// count. Every mass must be finite and positive.
    pub fn set_masses_kg(&mut self, masses_kg: &[f64]) -> Result<(), EngineError> {
        if masses_kg.len() != self.atom_count {
            return Err(EngineError::MassCountMismatch {
                provided: masses_kg.len(),
                expected: self.atom_count,
            });
        }
        for (atom, mass_kg) in masses_kg.iter().enumerate() {
            if !mass_kg.is_finite() || *mass_kg <= 0.0 {
                return Err(EngineError::InvalidMass { atom });
            }
        }
        self.masses_kg = masses_kg.to_vec();
        Ok(())
    }

    /// Returns the mass list, in kilograms. It is empty when no mass is set.
    pub fn masses_kg(&self) -> &[f64] {
        &self.masses_kg
    }

    /// Reports whether a mass is set for every atom.
    pub fn has_masses(&self) -> bool {
        self.masses_kg.len() == self.atom_count
    }

    /// Returns the neighbor skin, in metres.
    pub fn skin_m(&self) -> f64 {
        self.skin_m
    }

    /// Replaces the neighbor skin, in metres. The skin must be finite and
    /// non-negative. This invalidates the neighbor list.
    pub fn set_skin_m(&mut self, skin_m: f64) -> Result<(), EngineError> {
        if !skin_m.is_finite() || skin_m < 0.0 {
            return Err(EngineError::InvalidSkin { skin_m });
        }
        self.skin_m = skin_m;
        self.invalidate_neighbors();
        Ok(())
    }

    /// Returns the periodic box. It is non-periodic by default.
    pub fn periodic_box(&self) -> &PeriodicBox {
        &self.periodic_box
    }

    /// Replaces the periodic box. The neighbor list and the non-bonded terms
    /// all apply the same minimum image. The cutoff plus the skin must not
    /// exceed half of any periodic length. This invalidates the neighbor list.
    pub fn set_periodic_box(&mut self, periodic_box: PeriodicBox) -> Result<(), EngineError> {
        let reach_m = self.nonbonded_cutoff_m().unwrap_or(0.0) + self.skin_m;
        for length_m in periodic_box.lengths_m() {
            if length_m > 0.0 && reach_m > 0.5 * length_m {
                return Err(EngineError::CutoffExceedsHalfBox {
                    cutoff_m: reach_m,
                    length_m,
                });
            }
        }
        if let Some(term) = &mut self.van_der_waals {
            term.set_periodic_box(periodic_box)?;
        }
        if let Some(term) = &mut self.lennard_jones {
            term.set_periodic_box(periodic_box)?;
        }
        if let Some(term) = &mut self.electrostatic {
            term.set_periodic_box(periodic_box)?;
        }
        self.periodic_box = periodic_box;
        self.invalidate_neighbors();
        Ok(())
    }

    /// Sets the bond-stretch term. The term must cover at most the system atom
    /// count.
    pub fn set_bond_stretch(&mut self, term: BondStretchTerm) -> Result<(), EngineError> {
        self.check_term_atom_count("bond-stretch", term.atom_count())?;
        self.bond_stretch = Some(term);
        Ok(())
    }

    /// Sets the angle-bend term.
    pub fn set_angle_bend(&mut self, term: AngleBendTerm) -> Result<(), EngineError> {
        self.check_term_atom_count("angle-bend", term.atom_count())?;
        self.angle_bend = Some(term);
        Ok(())
    }

    /// Sets the torsion term.
    pub fn set_torsion(&mut self, term: TorsionTerm) -> Result<(), EngineError> {
        self.check_term_atom_count("torsion", term.atom_count())?;
        self.torsion = Some(term);
        Ok(())
    }

    /// Sets the out-of-plane term.
    pub fn set_out_of_plane(&mut self, term: OutOfPlaneTerm) -> Result<(), EngineError> {
        self.check_term_atom_count("out-of-plane", term.atom_count())?;
        self.out_of_plane = Some(term);
        Ok(())
    }

    /// Sets the van der Waals term. The term must cover at most the system atom
    /// count. This invalidates the neighbor list.
    pub fn set_van_der_waals(&mut self, term: VanDerWaalsTerm) -> Result<(), EngineError> {
        self.check_term_atom_count("van der Waals", term.atom_count())?;
        self.van_der_waals = Some(term);
        self.invalidate_neighbors();
        Ok(())
    }

    /// Sets the Lennard-Jones term. The term must cover at most the system
    /// atom count. This invalidates the neighbor list.
    pub fn set_lennard_jones(&mut self, term: LennardJonesTerm) -> Result<(), EngineError> {
        self.check_term_atom_count("Lennard-Jones", term.atom_count())?;
        self.lennard_jones = Some(term);
        self.invalidate_neighbors();
        Ok(())
    }

    /// Sets the electrostatic term. The term must cover at most the system atom
    /// count. This invalidates the neighbor list.
    pub fn set_electrostatic(&mut self, term: ElectrostaticTerm) -> Result<(), EngineError> {
        self.check_term_atom_count("electrostatic", term.atom_count())?;
        self.electrostatic = Some(term);
        self.invalidate_neighbors();
        Ok(())
    }

    /// Excludes an unordered atom pair from the non-bonded terms.
    ///
    /// Use this for bonded pairs and their near neighbors. The pair order does
    /// not matter. A self-pair is rejected. This invalidates the neighbor list.
    pub fn add_exclusion(&mut self, i: u32, j: u32) -> Result<(), EngineError> {
        if i == j {
            return Err(EngineError::PairIndexOutOfBounds {
                entry: 0,
                index: i,
                atom_count: self.atom_count,
            });
        }
        self.exclusions.insert(exclusion_key(i, j));
        self.invalidate_neighbors();
        Ok(())
    }

    /// Returns the number of excluded pairs.
    pub fn exclusion_count(&self) -> usize {
        self.exclusions.len()
    }

    /// Returns the number of pairs in the current neighbor list.
    pub fn neighbor_pair_count(&self) -> usize {
        self.verlet.as_ref().map_or(0, VerletList::pair_count)
    }

    /// Returns the number of bonds in the bond-stretch term.
    ///
    /// It is zero when no bond-stretch term is set.
    pub fn bond_count(&self) -> usize {
        self.bond_stretch
            .as_ref()
            .map_or(0, BondStretchTerm::bond_count)
    }

    /// Returns a read-only view of one bond, or `None` when the index is out
    /// of range or no bond-stretch term is set.
    pub fn bond(&self, bond: usize) -> Option<BondInfo> {
        self.bond_stretch.as_ref()?.bond_info(bond)
    }

    /// Returns the magnitude of the harmonic force on each bond, in newtons.
    ///
    /// The buffer holds one value per bond, in bond index order. The magnitude
    /// is `|k_n_per_m * (r - r0_m)|`. It is never negative. See
    /// [`System::bond_strains`] for the local strain.
    pub fn bond_forces_n(&self, positions_m: &[f64]) -> Result<Vec<f64>, EngineError> {
        validate_positions(positions_m, self.atom_count)?;
        match &self.bond_stretch {
            Some(term) => {
                let span = term.atom_count();
                term.bond_force_magnitudes_n(&positions_m[..3 * span])
            }
            None => Ok(Vec::new()),
        }
    }

    /// Returns the local strain of each bond, dimensionless.
    ///
    /// The buffer holds one value per bond, in bond index order. The strain is
    /// `(r - r0_m) / r0_m`.
    pub fn bond_strains(&self, positions_m: &[f64]) -> Result<Vec<f64>, EngineError> {
        validate_positions(positions_m, self.atom_count)?;
        match &self.bond_stretch {
            Some(term) => {
                let span = term.atom_count();
                term.bond_strains(&positions_m[..3 * span])
            }
            None => Ok(Vec::new()),
        }
    }

    /// Rebuilds the neighbor list if a non-bonded term is set and the atoms
    /// moved more than half the skin. Does nothing when no non-bonded term is
    /// set.
    pub fn rebuild_neighbors(&mut self, positions_m: &[f64]) -> Result<(), EngineError> {
        validate_positions(positions_m, self.atom_count)?;
        self.ensure_neighbors(positions_m)
    }

    /// Returns the total potential energy in joules.
    pub fn energy_j(&mut self, positions_m: &[f64]) -> Result<f64, EngineError> {
        Ok(self.energy_and_gradient_j(positions_m)?.0)
    }

    /// Returns the total energy gradient in newtons.
    ///
    /// Entry `3 * atom + axis` is `dU/dx`.
    pub fn gradient_j_per_m(&mut self, positions_m: &[f64]) -> Result<Vec<f64>, EngineError> {
        Ok(self.energy_and_gradient_j(positions_m)?.1)
    }

    /// Returns the total force in newtons, the negative of the gradient.
    pub fn forces_n(&mut self, positions_m: &[f64]) -> Result<Vec<f64>, EngineError> {
        let mut forces = self.gradient_j_per_m(positions_m)?;
        for force in &mut forces {
            *force = -*force;
        }
        Ok(forces)
    }

    /// Returns the total potential energy and the total gradient together.
    pub fn energy_and_gradient_j(
        &mut self,
        positions_m: &[f64],
    ) -> Result<(f64, Vec<f64>), EngineError> {
        validate_positions(positions_m, self.atom_count)?;
        self.ensure_neighbors(positions_m)?;

        let mut energy_j = 0.0;
        let mut gradient = vec![0.0; 3 * self.atom_count];

        self.accumulate_bonded(positions_m, &mut energy_j, &mut gradient)?;
        self.accumulate_nonbonded_blocked(positions_m, &mut energy_j, &mut gradient)?;
        Ok((energy_j, gradient))
    }

    /// Builds the four-lane non-bonded kernel over the active terms.
    fn nonbonded_kernel<'a>(&'a self, positions_m: &'a [f64]) -> NonbondedKernel<'a> {
        NonbondedKernel {
            positions_m,
            periodic_box: self.periodic_box,
            van_der_waals: self.van_der_waals.as_ref().map(|term| KernelVanDerWaals {
                a_j: term.a_j(),
                b_per_m: term.b_per_m(),
                c_j_m6: term.c_j_m6(),
                cutoff: *term.cutoff(),
            }),
            lennard_jones: self.lennard_jones.as_ref().map(|term| KernelLennardJones {
                epsilon_j: term.epsilon_j(),
                sigma_m: term.sigma_m(),
                cutoff: *term.cutoff(),
            }),
            electrostatic: self.electrostatic.as_ref().map(|term| KernelElectrostatic {
                charges_c: term.charges_c(),
                cutoff: *term.cutoff(),
            }),
        }
    }

    /// Accumulates the non-bonded terms over the canonical pair blocks.
    ///
    /// The canonical block partition depends only on the pair count. The
    /// per-block partial gradients are summed in block order. The parallel
    /// path uses the same partition and the same order, so the two paths agree
    /// bit for bit.
    fn accumulate_nonbonded_blocked(
        &self,
        positions_m: &[f64],
        energy_j: &mut f64,
        gradient: &mut [f64],
    ) -> Result<(), EngineError> {
        let pair_count = self.nonbonded_pairs.len() / 2;
        if pair_count == 0
            || (self.van_der_waals.is_none()
                && self.lennard_jones.is_none()
                && self.electrostatic.is_none())
        {
            return Ok(());
        }
        let kernel = self.nonbonded_kernel(positions_m);
        let block_count = pair_count.div_ceil(PARALLEL_BLOCK_PAIRS);
        let mut block_gradient = vec![0.0; gradient.len()];
        for block in 0..block_count {
            let start = block * PARALLEL_BLOCK_PAIRS;
            let end = ((block + 1) * PARALLEL_BLOCK_PAIRS).min(pair_count);
            let pairs = &self.nonbonded_pairs[2 * start..2 * end];
            block_gradient.fill(0.0);
            let mut block_energy_j = 0.0;
            kernel.accumulate_lanes4(pairs, &mut block_energy_j, &mut block_gradient)?;
            *energy_j += block_energy_j;
            for (slot, value) in gradient.iter_mut().zip(block_gradient.iter()) {
                *slot += *value;
            }
        }
        Ok(())
    }

    fn accumulate_bonded(
        &self,
        positions_m: &[f64],
        energy_j: &mut f64,
        gradient: &mut [f64],
    ) -> Result<(), EngineError> {
        if let Some(term) = &self.bond_stretch {
            let span = term.atom_count();
            if span > 0 {
                let (e, g) = term.energy_and_gradient_j(&positions_m[..3 * span])?;
                accumulate(energy_j, gradient, e, &g);
            }
        }
        if let Some(term) = &self.angle_bend {
            let span = term.atom_count();
            if span > 0 {
                let (e, g) = term.energy_and_gradient_j(&positions_m[..3 * span])?;
                accumulate(energy_j, gradient, e, &g);
            }
        }
        if let Some(term) = &self.torsion {
            let span = term.atom_count();
            if span > 0 {
                let (e, g) = term.energy_and_gradient_j(&positions_m[..3 * span])?;
                accumulate(energy_j, gradient, e, &g);
            }
        }
        if let Some(term) = &self.out_of_plane {
            let span = term.atom_count();
            if span > 0 {
                let (e, g) = term.energy_and_gradient_j(&positions_m[..3 * span])?;
                accumulate(energy_j, gradient, e, &g);
            }
        }
        Ok(())
    }

    /// Returns the total force in newtons with a deterministic parallel
    /// reduction.
    ///
    /// The result is bit-identical run to run and independent of
    /// `thread_count`. The reduction partition does not depend on the thread
    /// count, so the floating-point sum order is fixed. The default serial
    /// path is unchanged. `thread_count` must be at least one.
    pub fn forces_n_parallel(
        &mut self,
        positions_m: &[f64],
        thread_count: usize,
    ) -> Result<Vec<f64>, EngineError> {
        let mut forces = self.gradient_j_per_m_parallel(positions_m, thread_count)?;
        for force in &mut forces {
            *force = -*force;
        }
        Ok(forces)
    }

    /// Returns the total energy gradient in newtons with a deterministic
    /// parallel reduction over the neighbor pair list.
    ///
    /// See [`System::forces_n_parallel`] for the determinism contract.
    pub fn gradient_j_per_m_parallel(
        &mut self,
        positions_m: &[f64],
        thread_count: usize,
    ) -> Result<Vec<f64>, EngineError> {
        if thread_count == 0 {
            return Err(EngineError::InvalidThreadCount { thread_count });
        }
        validate_positions(positions_m, self.atom_count)?;
        self.ensure_neighbors(positions_m)?;

        let mut energy_j = 0.0;
        let mut gradient = vec![0.0; 3 * self.atom_count];
        self.accumulate_bonded(positions_m, &mut energy_j, &mut gradient)?;

        let atom_count = self.atom_count;
        let pair_count = self.nonbonded_pairs.len() / 2;
        let has_nonbonded = self.van_der_waals.is_some()
            || self.lennard_jones.is_some()
            || self.electrostatic.is_some();
        if pair_count == 0 || !has_nonbonded {
            return Ok(gradient);
        }

        let block_count = pair_count.div_ceil(PARALLEL_BLOCK_PAIRS);
        let worker_count = thread_count.min(block_count);
        let block_len = 3 * atom_count;
        let kernel = self.nonbonded_kernel(positions_m);
        let kernel = &kernel;
        let pairs: &[u32] = &self.nonbonded_pairs;

        std::thread::scope(|scope| -> Result<Vec<Vec<f64>>, EngineError> {
            let mut handles = Vec::with_capacity(worker_count);
            for worker in 0..worker_count {
                let start_block = worker * block_count / worker_count;
                let end_block = (worker + 1) * block_count / worker_count;
                handles.push(scope.spawn(move || -> Result<Vec<f64>, EngineError> {
                    let mut local = vec![0.0; (end_block - start_block) * block_len];
                    for block in start_block..end_block {
                        let start = block * PARALLEL_BLOCK_PAIRS;
                        let end = ((block + 1) * PARALLEL_BLOCK_PAIRS).min(pair_count);
                        let offset = (block - start_block) * block_len;
                        let block_gradient = &mut local[offset..offset + block_len];
                        let mut block_energy_j = 0.0;
                        kernel.accumulate_lanes4(
                            &pairs[2 * start..2 * end],
                            &mut block_energy_j,
                            block_gradient,
                        )?;
                    }
                    Ok(local)
                }));
            }
            let mut partials = Vec::with_capacity(worker_count);
            for handle in handles {
                let local = handle.join().map_err(|_| EngineError::WorkerPanicked)?;
                partials.push(local?);
            }
            Ok(partials)
        })
        .map(|partials| {
            for (worker, local) in partials.into_iter().enumerate() {
                let start_block = worker * block_count / worker_count;
                for block in start_block..(start_block + local.len() / block_len) {
                    let offset = (block - start_block) * block_len;
                    for (slot, value) in gradient
                        .iter_mut()
                        .zip(local[offset..offset + block_len].iter())
                    {
                        *slot += *value;
                    }
                }
            }
            gradient
        })
    }

    /// Returns the kinetic energy in joules for a velocity buffer.
    ///
    /// Velocities are in metres per second. A mass must be set.
    pub fn kinetic_energy_j(&self, velocities_m_per_s: &[f64]) -> Result<f64, EngineError> {
        if velocities_m_per_s.len() != 3 * self.atom_count {
            return Err(EngineError::BufferSizeMismatch {
                len: velocities_m_per_s.len(),
                expected: 3 * self.atom_count,
            });
        }
        if !self.has_masses() {
            return Err(EngineError::MassesNotSet);
        }
        let mut kinetic_j = 0.0;
        for atom in 0..self.atom_count {
            let base = atom * 3;
            let vx = velocities_m_per_s[base];
            let vy = velocities_m_per_s[base + 1];
            let vz = velocities_m_per_s[base + 2];
            kinetic_j += 0.5 * self.masses_kg[atom] * (vx * vx + vy * vy + vz * vz);
        }
        Ok(kinetic_j)
    }

    /// Returns the total energy in joules, potential plus kinetic.
    pub fn total_energy_j(
        &mut self,
        positions_m: &[f64],
        velocities_m_per_s: &[f64],
    ) -> Result<f64, EngineError> {
        let potential_j = self.energy_j(positions_m)?;
        let kinetic_j = self.kinetic_energy_j(velocities_m_per_s)?;
        Ok(potential_j + kinetic_j)
    }

    /// Returns the instantaneous kinetic temperature in kelvin.
    ///
    /// The value is `2 K / (3 N k_B)` for `N` atoms. An empty system returns
    /// zero. A mass must be set.
    pub fn temperature_k(&self, velocities_m_per_s: &[f64]) -> Result<f64, EngineError> {
        if self.atom_count == 0 {
            return Ok(0.0);
        }
        let kinetic_j = self.kinetic_energy_j(velocities_m_per_s)?;
        let degrees_of_freedom = 3 * self.atom_count;
        Ok(2.0 * kinetic_j / (degrees_of_freedom as f64 * BOLTZMANN_J_PER_K))
    }

    fn check_term_atom_count(
        &self,
        term: &'static str,
        provided: usize,
    ) -> Result<(), EngineError> {
        if provided > self.atom_count {
            return Err(EngineError::TermAtomCountMismatch {
                term,
                provided,
                expected: self.atom_count,
            });
        }
        Ok(())
    }

    fn invalidate_neighbors(&mut self) {
        self.verlet = None;
        self.nonbonded_pairs.clear();
    }

    fn nonbonded_cutoff_m(&self) -> Option<f64> {
        let mut cutoff_m: Option<f64> = None;
        if let Some(term) = &self.van_der_waals {
            if term.atom_count() > 1 {
                cutoff_m = Some(max_option(cutoff_m, term.cutoff().cutoff_m()));
            }
        }
        if let Some(term) = &self.lennard_jones {
            if term.atom_count() > 1 {
                cutoff_m = Some(max_option(cutoff_m, term.cutoff().cutoff_m()));
            }
        }
        if let Some(term) = &self.electrostatic {
            if term.atom_count() > 1 {
                cutoff_m = Some(max_option(cutoff_m, term.cutoff().cutoff_m()));
            }
        }
        cutoff_m
    }

    fn ensure_neighbors(&mut self, positions_m: &[f64]) -> Result<(), EngineError> {
        let Some(cutoff_m) = self.nonbonded_cutoff_m() else {
            return Ok(());
        };
        let rebuild = match &self.verlet {
            None => true,
            Some(list) => list.needs_rebuild(positions_m),
        };
        if rebuild {
            let mut list = VerletList::with_periodic_box(cutoff_m, self.skin_m, self.periodic_box)?;
            list.build(positions_m)?;
            self.verlet = Some(list);
            self.nonbonded_pairs = self.filtered_pairs();
        }
        Ok(())
    }

    fn filtered_pairs(&self) -> Vec<u32> {
        let Some(list) = &self.verlet else {
            return Vec::new();
        };
        if self.exclusions.is_empty() {
            return list.pairs().to_vec();
        }
        let mut kept = Vec::with_capacity(list.pairs().len());
        for chunk in list.pairs().as_chunks::<2>().0 {
            if !self.exclusions.contains(&exclusion_key(chunk[0], chunk[1])) {
                kept.push(chunk[0]);
                kept.push(chunk[1]);
            }
        }
        kept
    }
}

fn accumulate(energy_j: &mut f64, gradient: &mut [f64], term_energy_j: f64, term_gradient: &[f64]) {
    *energy_j += term_energy_j;
    for (slot, value) in gradient.iter_mut().zip(term_gradient.iter()) {
        *slot += value;
    }
}

fn exclusion_key(i: u32, j: u32) -> u64 {
    let low = i.min(j) as u64;
    let high = i.max(j) as u64;
    (low << 32) | high
}

fn max_option(current: Option<f64>, value: f64) -> f64 {
    match current {
        Some(existing) => existing.max(value),
        None => value,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nonbonded::{all_pairs, Cutoff, PeriodicBox};
    use crate::test_support::{
        central_difference, full_system, small_molecule_positions_m, Rng, CARBON_MASS_KG,
    };
    use crate::van_der_waals::VdwParams;

    fn carbon_masses(atom_count: usize) -> Vec<f64> {
        vec![CARBON_MASS_KG; atom_count]
    }

    #[test]
    fn an_empty_system_has_zero_energy_and_an_empty_gradient() {
        let mut system = System::new(0);
        let (energy_j, gradient) = system.energy_and_gradient_j(&[]).expect("valid");
        assert_eq!(energy_j, 0.0);
        assert!(gradient.is_empty());
    }

    #[test]
    fn a_term_that_covers_too_many_atoms_is_rejected() {
        let mut system = System::new(2);
        let mut term = BondStretchTerm::new();
        term.add_bond(0, 1, 300.0, 1.5e-10).expect("valid");
        term.add_bond(1, 2, 300.0, 1.5e-10).expect("valid");
        assert!(matches!(
            system.set_bond_stretch(term),
            Err(EngineError::TermAtomCountMismatch { .. })
        ));
    }

    #[test]
    fn a_mass_count_mismatch_and_an_invalid_mass_are_rejected() {
        let mut system = System::new(2);
        assert!(matches!(
            system.set_masses_kg(&[1.0e-26]),
            Err(EngineError::MassCountMismatch { .. })
        ));
        assert!(matches!(
            system.set_masses_kg(&[1.0e-26, 0.0]),
            Err(EngineError::InvalidMass { .. })
        ));
    }

    #[test]
    fn the_analytic_gradient_matches_a_central_finite_difference() {
        let positions_m = small_molecule_positions_m();
        let mut system = full_system(&positions_m, &carbon_masses(6));
        let analytic = system.gradient_j_per_m(&positions_m).expect("valid");

        let step_m = 1.0e-14;
        let finite_difference =
            central_difference(|p| system.energy_j(p).expect("valid"), &positions_m, step_m);
        let mut max_error_n = 0.0_f64;
        let rtol = 1.0e-5;
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
        println!("system max gradient error: {max_error_n:e} N");
    }

    #[test]
    fn every_term_contributes_to_the_total_energy() {
        let positions_m = small_molecule_positions_m();
        let masses = carbon_masses(6);
        let mut bonded = crate::test_support::bonded_system(&positions_m, &masses);
        let mut full = full_system(&positions_m, &masses);
        let bonded_energy_j = bonded.energy_j(&positions_m).expect("valid");
        let full_energy_j = full.energy_j(&positions_m).expect("valid");
        assert!(bonded_energy_j > 0.0);
        assert!(full_energy_j != bonded_energy_j);
    }

    #[test]
    fn the_nonbonded_pair_list_path_matches_the_all_pairs_reference() {
        let mut rng = Rng::new(0x0FF1CE);
        let mut positions_m = Vec::with_capacity(3 * 60);
        for _ in 0..60 {
            positions_m.push(rng.next_f64() * 3.0e-9);
            positions_m.push(rng.next_f64() * 3.0e-9);
            positions_m.push(rng.next_f64() * 3.0e-9);
        }
        let cutoff = Cutoff::new(1.0e-9, 0.8e-9).expect("valid cutoff");
        let params: Vec<VdwParams> = (0..60)
            .map(|atom| VdwParams::new(2.2e-19 + 1.0e-21 * atom as f64, 4.2e10, 8.0e-79))
            .collect();
        let vdw = VanDerWaalsTerm::from_params(&params, cutoff, PeriodicBox::non_periodic())
            .expect("valid term");
        let mut system = System::with_masses_kg(60, &carbon_masses(60)).expect("valid");
        system.set_van_der_waals(vdw.clone()).expect("fits");

        let system_energy_j = system.energy_j(&positions_m).expect("valid");
        let reference_energy_j = vdw.energy_j(&positions_m).expect("valid");
        assert!(system.neighbor_pair_count() < all_pairs(60).len() / 2);
        assert!(
            (system_energy_j - reference_energy_j).abs()
                <= 1.0e-24 + 1.0e-9 * reference_energy_j.abs()
        );
        assert!(system_energy_j.abs() > 1.0e-26);
    }

    #[test]
    fn an_excluded_pair_is_removed_from_the_nonbonded_energy() {
        let positions_m = vec![0.0, 0.0, 0.0, 3.0e-10, 0.0, 0.0];
        let cutoff = Cutoff::new(1.0e-9, 0.8e-9).expect("valid cutoff");
        let params = vec![
            VdwParams::new(2.2e-19, 4.2e10, 8.0e-79),
            VdwParams::new(2.2e-19, 4.2e10, 8.0e-79),
        ];
        let vdw = VanDerWaalsTerm::from_params(&params, cutoff, PeriodicBox::non_periodic())
            .expect("valid term");
        let mut with_pair = System::with_masses_kg(2, &carbon_masses(2)).expect("valid");
        with_pair.set_van_der_waals(vdw).expect("fits");
        let with_energy_j = with_pair.energy_j(&positions_m).expect("valid");

        let vdw = VanDerWaalsTerm::from_params(&params, cutoff, PeriodicBox::non_periodic())
            .expect("valid term");
        let mut excluded = System::with_masses_kg(2, &carbon_masses(2)).expect("valid");
        excluded.set_van_der_waals(vdw).expect("fits");
        excluded.add_exclusion(0, 1).expect("valid exclusion");
        let excluded_energy_j = excluded.energy_j(&positions_m).expect("valid");

        assert!(with_energy_j.abs() > 1.0e-26);
        assert_eq!(excluded_energy_j, 0.0);
    }

    #[test]
    fn the_kinetic_energy_matches_the_analytic_value() {
        let system = System::with_masses_kg(1, &[2.0]).expect("valid mass");
        let kinetic_j = system
            .kinetic_energy_j(&[3.0, 4.0, 0.0])
            .expect("valid velocity");
        assert!((kinetic_j - 25.0).abs() < 1.0e-12);
        assert!(matches!(
            system.kinetic_energy_j(&[1.0]),
            Err(EngineError::BufferSizeMismatch { .. })
        ));
    }

    fn periodic_positions_m() -> Vec<f64> {
        let mut rng = Rng::new(0x9E3779B9);
        let box_m = 2.0e-9;
        let mut positions_m = Vec::with_capacity(24);
        positions_m.extend_from_slice(&[0.05e-9, 0.05e-9, 0.05e-9]);
        positions_m.extend_from_slice(&[1.95e-9, 0.05e-9, 0.05e-9]);
        for _ in 0..6 {
            positions_m.push(rng.next_f64() * box_m);
            positions_m.push(rng.next_f64() * box_m);
            positions_m.push(rng.next_f64() * box_m);
        }
        positions_m
    }

    fn periodic_cutoff() -> Cutoff {
        Cutoff::new(0.5e-9, 0.4e-9).expect("valid cutoff")
    }

    fn periodic_vdw_term(atom_count: usize, periodic_box: PeriodicBox) -> VanDerWaalsTerm {
        let params: Vec<VdwParams> = (0..atom_count)
            .map(|atom| VdwParams::new(2.2e-19 + 1.0e-21 * atom as f64, 4.2e10, 8.0e-79))
            .collect();
        VanDerWaalsTerm::from_params(&params, periodic_cutoff(), periodic_box).expect("valid term")
    }

    fn periodic_electrostatic_term(
        atom_count: usize,
        periodic_box: PeriodicBox,
    ) -> ElectrostaticTerm {
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
        ElectrostaticTerm::from_charges(&charges_c, periodic_cutoff(), periodic_box)
            .expect("valid term")
    }

    #[test]
    fn the_periodic_pair_list_energy_matches_the_all_pairs_reference() {
        let positions_m = periodic_positions_m();
        let atom_count = positions_m.len() / 3;
        let periodic_box = PeriodicBox::new([2.0e-9, 2.0e-9, 2.0e-9]).expect("valid box");
        let vdw = periodic_vdw_term(atom_count, periodic_box);
        let electrostatic = periodic_electrostatic_term(atom_count, periodic_box);

        let mut system =
            System::with_masses_kg(atom_count, &carbon_masses(atom_count)).expect("valid masses");
        system.set_periodic_box(periodic_box).expect("valid box");
        system.set_van_der_waals(vdw.clone()).expect("fits");
        system
            .set_electrostatic(electrostatic.clone())
            .expect("fits");

        let system_energy_j = system.energy_j(&positions_m).expect("valid");
        let reference_energy_j = vdw.energy_j(&positions_m).expect("valid")
            + electrostatic.energy_j(&positions_m).expect("valid");
        assert!(system_energy_j.abs() > 1.0e-26);
        assert!(
            (system_energy_j - reference_energy_j).abs()
                <= 1.0e-24 + 1.0e-9 * reference_energy_j.abs(),
            "system {system_energy_j:e} reference {reference_energy_j:e}"
        );
    }

    #[test]
    fn the_periodic_analytic_gradient_matches_a_central_finite_difference() {
        let positions_m = periodic_positions_m();
        let atom_count = positions_m.len() / 3;
        let periodic_box = PeriodicBox::new([2.0e-9, 2.0e-9, 2.0e-9]).expect("valid box");
        let vdw = periodic_vdw_term(atom_count, periodic_box);
        let electrostatic = periodic_electrostatic_term(atom_count, periodic_box);
        let mut system =
            System::with_masses_kg(atom_count, &carbon_masses(atom_count)).expect("valid masses");
        system.set_periodic_box(periodic_box).expect("valid box");
        system.set_van_der_waals(vdw).expect("fits");
        system.set_electrostatic(electrostatic).expect("fits");

        let analytic = system.gradient_j_per_m(&positions_m).expect("valid");
        let step_m = 1.0e-15;
        let finite_difference =
            central_difference(|p| system.energy_j(p).expect("valid"), &positions_m, step_m);
        let mut max_error_n = 0.0_f64;
        let rtol = 1.0e-4;
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
        println!("periodic system max gradient error: {max_error_n:e} N");
    }

    #[test]
    fn the_deterministic_parallel_forces_match_for_one_and_many_threads() {
        let positions_m = periodic_positions_m();
        let atom_count = positions_m.len() / 3;
        let periodic_box = PeriodicBox::new([2.0e-9, 2.0e-9, 2.0e-9]).expect("valid box");
        let mut system =
            System::with_masses_kg(atom_count, &carbon_masses(atom_count)).expect("valid masses");
        system.set_periodic_box(periodic_box).expect("valid box");
        system
            .set_van_der_waals(periodic_vdw_term(atom_count, periodic_box))
            .expect("fits");
        system
            .set_electrostatic(periodic_electrostatic_term(atom_count, periodic_box))
            .expect("fits");

        let serial = system.forces_n(&positions_m).expect("valid");
        let one_thread = system.forces_n_parallel(&positions_m, 1).expect("valid");
        let four_threads = system.forces_n_parallel(&positions_m, 4).expect("valid");
        let many_threads = system.forces_n_parallel(&positions_m, 17).expect("valid");

        for index in 0..serial.len() {
            assert_eq!(one_thread[index].to_bits(), four_threads[index].to_bits());
            assert_eq!(one_thread[index].to_bits(), many_threads[index].to_bits());
        }
        assert!(matches!(
            system.forces_n_parallel(&positions_m, 0),
            Err(EngineError::InvalidThreadCount { .. })
        ));
    }

    #[test]
    fn parallel_forces_are_bit_identical_to_serial_for_one_two_four_and_eight_workers() {
        // A dense periodic system gives many canonical blocks, so the worker
        // chunks are real. A sparse system would fit in one block and would
        // not exercise the partition.
        let atom_count = 400;
        let box_m = 3.0e-9;
        let mut rng = Rng::new(0x5EED);
        let mut positions_m = Vec::with_capacity(3 * atom_count);
        for _ in 0..atom_count {
            positions_m.push(rng.next_f64() * box_m);
            positions_m.push(rng.next_f64() * box_m);
            positions_m.push(rng.next_f64() * box_m);
        }
        let periodic_box = PeriodicBox::new([box_m, box_m, box_m]).expect("valid box");
        let cutoff = Cutoff::new(1.0e-9, 0.8e-9).expect("valid cutoff");
        let params: Vec<VdwParams> = (0..atom_count)
            .map(|atom| VdwParams::new(2.2e-19 + 1.0e-21 * atom as f64, 4.2e10, 8.0e-79))
            .collect();
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
        let mut system =
            System::with_masses_kg(atom_count, &carbon_masses(atom_count)).expect("valid masses");
        system.set_periodic_box(periodic_box).expect("valid box");
        system
            .set_van_der_waals(
                VanDerWaalsTerm::from_params(&params, cutoff, periodic_box).expect("valid term"),
            )
            .expect("fits");
        system
            .set_electrostatic(
                ElectrostaticTerm::from_charges(&charges_c, cutoff, periodic_box)
                    .expect("valid term"),
            )
            .expect("fits");

        let serial = system.forces_n(&positions_m).expect("valid");
        let block_count = system.neighbor_pair_count().div_ceil(PARALLEL_BLOCK_PAIRS);
        assert!(
            block_count > 1,
            "expected more than one block, got {block_count}"
        );
        println!(
            "parallel bit-identity: {} atoms, {} pairs, {} canonical blocks",
            atom_count,
            system.neighbor_pair_count(),
            block_count
        );
        for thread_count in [1_usize, 2, 4, 8] {
            let parallel = system
                .forces_n_parallel(&positions_m, thread_count)
                .expect("valid");
            for index in 0..serial.len() {
                assert_eq!(
                    serial[index].to_bits(),
                    parallel[index].to_bits(),
                    "thread count {thread_count}, coordinate {index}"
                );
            }
        }
    }

    #[test]
    fn setting_the_system_box_moves_the_nonbonded_terms_to_it() {
        let positions_m = periodic_positions_m();
        let atom_count = positions_m.len() / 3;
        let periodic_box = PeriodicBox::new([2.0e-9, 2.0e-9, 2.0e-9]).expect("valid box");

        let mut system =
            System::with_masses_kg(atom_count, &carbon_masses(atom_count)).expect("valid masses");
        system
            .set_van_der_waals(periodic_vdw_term(atom_count, PeriodicBox::non_periodic()))
            .expect("fits");
        system
            .set_electrostatic(periodic_electrostatic_term(
                atom_count,
                PeriodicBox::non_periodic(),
            ))
            .expect("fits");
        system.set_periodic_box(periodic_box).expect("valid box");

        assert_eq!(
            system.van_der_waals.as_ref().unwrap().periodic_box(),
            &periodic_box
        );
        assert_eq!(
            system.electrostatic.as_ref().unwrap().periodic_box(),
            &periodic_box
        );
        let system_energy_j = system.energy_j(&positions_m).expect("valid");
        let reference_energy_j = periodic_vdw_term(atom_count, periodic_box)
            .energy_j(&positions_m)
            .expect("valid")
            + periodic_electrostatic_term(atom_count, periodic_box)
                .energy_j(&positions_m)
                .expect("valid");
        assert!(
            (system_energy_j - reference_energy_j).abs()
                <= 1.0e-24 + 1.0e-9 * reference_energy_j.abs()
        );
    }

    #[test]
    fn the_parallel_forces_match_the_serial_forces_within_tolerance() {
        let positions_m = periodic_positions_m();
        let atom_count = positions_m.len() / 3;
        let periodic_box = PeriodicBox::new([2.0e-9, 2.0e-9, 2.0e-9]).expect("valid box");
        let mut system =
            System::with_masses_kg(atom_count, &carbon_masses(atom_count)).expect("valid masses");
        system.set_periodic_box(periodic_box).expect("valid box");
        system
            .set_van_der_waals(periodic_vdw_term(atom_count, periodic_box))
            .expect("fits");
        system
            .set_electrostatic(periodic_electrostatic_term(atom_count, periodic_box))
            .expect("fits");

        let serial = system.forces_n(&positions_m).expect("valid");
        let parallel = system.forces_n_parallel(&positions_m, 8).expect("valid");
        for index in 0..serial.len() {
            let error_n = (serial[index] - parallel[index]).abs();
            assert!(
                error_n <= 1.0e-18 + 1.0e-10 * serial[index].abs(),
                "coordinate {index}: serial {} parallel {} error {}",
                serial[index],
                parallel[index],
                error_n
            );
        }
    }

    #[test]
    fn the_bond_accessors_read_the_bond_stretch_term() {
        let positions_m = small_molecule_positions_m();
        let system = crate::test_support::bonded_system(&positions_m, &carbon_masses(6));
        assert_eq!(system.bond_count(), 5);
        let bond = system.bond(0).expect("bond 0");
        assert_eq!((bond.u, bond.v), (0, 1));
        assert_eq!(bond.order, 1);
        assert_eq!(bond.r0_m, 1.5e-10);
        assert!(system.bond(99).is_none());

        let magnitudes_n = system.bond_forces_n(&positions_m).expect("valid");
        let strains = system.bond_strains(&positions_m).expect("valid");
        assert_eq!(magnitudes_n.len(), 5);
        assert_eq!(strains.len(), 5);
        for (magnitude_n, strain) in magnitudes_n.iter().zip(strains.iter()) {
            assert!((magnitude_n - 300.0 * strain.abs() * 1.5e-10).abs() < 1.0e-20);
        }
    }

    #[test]
    fn an_empty_bond_term_reports_no_bonds() {
        let system = System::new(1);
        assert_eq!(system.bond_count(), 0);
        assert!(system.bond(0).is_none());
        assert!(system
            .bond_forces_n(&[0.0, 0.0, 0.0])
            .expect("valid")
            .is_empty());
        assert!(system
            .bond_strains(&[0.0, 0.0, 0.0])
            .expect("valid")
            .is_empty());
        assert!(matches!(
            system.bond_forces_n(&[0.0]),
            Err(EngineError::PositionBufferSizeMismatch { .. })
        ));
    }

    #[test]
    fn the_bond_force_magnitude_matches_a_finite_difference_of_the_bond_energy() {
        let positions_m = small_molecule_positions_m();
        let atom_count = positions_m.len() / 3;
        let mut system =
            System::with_masses_kg(atom_count, &carbon_masses(atom_count)).expect("valid");
        let mut term = BondStretchTerm::new();
        for [u, v] in [[0u32, 1u32], [1, 2], [2, 3], [3, 4], [4, 5]] {
            term.add_bond(u, v, 300.0, 1.5e-10).expect("valid bond");
        }
        system.set_bond_stretch(term).expect("bond term fits");

        let magnitudes_n = system.bond_forces_n(&positions_m).expect("valid");
        let step_m = 1.0e-14;
        for (bond, magnitude_n) in magnitudes_n.iter().enumerate() {
            let info = system.bond(bond).expect("bond");
            let mut bond_energy = BondStretchTerm::new();
            bond_energy
                .add_bond(info.u, info.v, info.k_n_per_m, info.r0_m)
                .expect("valid bond");
            let u = info.u as usize * 3;
            let v = info.v as usize * 3;
            let mut dx = positions_m[u] - positions_m[v];
            let mut dy = positions_m[u + 1] - positions_m[v + 1];
            let mut dz = positions_m[u + 2] - positions_m[v + 2];
            let distance_m = (dx * dx + dy * dy + dz * dz).sqrt();
            dx /= distance_m;
            dy /= distance_m;
            dz /= distance_m;

            let mut forward_m = positions_m.clone();
            forward_m[u] += step_m * dx;
            forward_m[u + 1] += step_m * dy;
            forward_m[u + 2] += step_m * dz;
            let mut backward_m = positions_m.clone();
            backward_m[u] -= step_m * dx;
            backward_m[u + 1] -= step_m * dy;
            backward_m[u + 2] -= step_m * dz;
            let span = bond_energy.atom_count();
            let forward_j = bond_energy.energy_j(&forward_m[..3 * span]).expect("valid");
            let backward_j = bond_energy
                .energy_j(&backward_m[..3 * span])
                .expect("valid");
            let finite_difference_n = ((forward_j - backward_j) / (2.0 * step_m)).abs();
            let error_n = (magnitude_n - finite_difference_n).abs();
            assert!(
                error_n <= 1.0e-20 + 1.0e-6 * finite_difference_n,
                "bond {bond}: magnitude {} finite difference {} error {}",
                magnitude_n,
                finite_difference_n,
                error_n
            );
        }
    }
}
