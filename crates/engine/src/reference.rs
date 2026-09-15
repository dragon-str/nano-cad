//! A flexible SPC-like water box for periodic reference runs.
//!
//! The model is close to SPC: a three-site water with point charges, a
//! Buckingham oxygen-oxygen repulsion-dispersion term, harmonic bonds, and a
//! harmonic angle. It is not SPC exactly. The bonded terms are flexible, and
//! the dispersion term is Buckingham rather than Lennard-Jones.
//!
//! The box is periodic. The reference case is self-consistent: the neighbor
//! list and the non-bonded terms use the same minimum image.
//!
//! A stable NVE run for 200 molecules needs a time step of 1e-17 s. A step of
//! 1e-16 s is too large for the flexible O-H mode of this model.

use crate::angle_bend::{AngleBendParams, AngleBendTerm};
use crate::bond_stretch::BondStretchTerm;
use crate::electrostatic::ElectrostaticTerm;
use crate::error::EngineError;
use crate::nonbonded::{Cutoff, PeriodicBox};
use crate::rng::Rng;
use crate::system::{System, BOLTZMANN_J_PER_K};
use crate::van_der_waals::{VanDerWaalsTerm, VdwParams};

/// The mass of an oxygen-16 atom in kilograms.
pub const OXYGEN_MASS_KG: f64 = 2.656_696_2e-26;
/// The mass of a hydrogen-1 atom in kilograms.
pub const HYDROGEN_MASS_KG: f64 = 1.673_557_5e-27;
/// The elementary charge in coulombs.
pub const ELEMENTARY_CHARGE_C: f64 = 1.602_176_634e-19;
/// The SPC oxygen partial charge in coulombs.
pub const SPC_OXYGEN_CHARGE_C: f64 = -0.82 * ELEMENTARY_CHARGE_C;
/// The SPC hydrogen partial charge in coulombs.
pub const SPC_HYDROGEN_CHARGE_C: f64 = 0.41 * ELEMENTARY_CHARGE_C;
/// The SPC oxygen-hydrogen bond length in metres.
pub const SPC_BOND_LENGTH_M: f64 = 1.0e-10;
/// The SPC hydrogen-oxygen-hydrogen angle in radians.
pub const SPC_ANGLE_RAD: f64 = 109.47 * std::f64::consts::PI / 180.0;
/// The harmonic oxygen-hydrogen force constant in newtons per metre.
pub const SPC_BOND_K_N_PER_M: f64 = 5.0e2;
/// The harmonic angle force constant in joules per radian squared.
pub const SPC_ANGLE_K_J_PER_RAD2: f64 = 6.0e-19;
/// The Buckingham oxygen-oxygen repulsive prefactor in joules.
pub const SPC_VDW_A_J: f64 = 1.151e-15;
/// The Buckingham oxygen-oxygen repulsive decay constant in reciprocal metres.
pub const SPC_VDW_B_PER_M: f64 = 4.2e10;
/// The Buckingham oxygen-oxygen dispersion coefficient in joule metre to the sixth power.
pub const SPC_VDW_C_J_M6: f64 = 2.6645e-78;

/// A periodic SPC-like water box with positions and velocities.
#[derive(Clone, Debug, PartialEq)]
pub struct WaterBox {
    system: System,
    positions_m: Vec<f64>,
    velocities_m_per_s: Vec<f64>,
    periodic_box: PeriodicBox,
    molecule_count: usize,
}

impl WaterBox {
    /// Returns the number of water molecules.
    pub fn molecule_count(&self) -> usize {
        self.molecule_count
    }

    /// Returns the number of atoms.
    pub fn atom_count(&self) -> usize {
        self.system.atom_count()
    }

    /// Returns the periodic box.
    pub fn periodic_box(&self) -> &PeriodicBox {
        &self.periodic_box
    }

    /// Returns the flat positions in metres.
    pub fn positions_m(&self) -> &[f64] {
        &self.positions_m
    }

    /// Returns the flat velocities in metres per second.
    pub fn velocities_m_per_s(&self) -> &[f64] {
        &self.velocities_m_per_s
    }

    /// Returns the system.
    pub fn system(&self) -> &System {
        &self.system
    }

    /// Returns the system for mutation, for example to run an integrator.
    pub fn system_mut(&mut self) -> &mut System {
        &mut self.system
    }

    /// Returns the flat positions for mutation.
    pub fn positions_mut(&mut self) -> &mut [f64] {
        &mut self.positions_m
    }

    /// Returns the flat velocities for mutation.
    pub fn velocities_mut(&mut self) -> &mut [f64] {
        &mut self.velocities_m_per_s
    }

    /// Returns the instantaneous kinetic temperature in kelvin.
    pub fn temperature_k(&self) -> Result<f64, EngineError> {
        self.system.temperature_k(&self.velocities_m_per_s)
    }

    /// Returns the potential energy in joules.
    pub fn potential_energy_j(&mut self) -> Result<f64, EngineError> {
        self.system.energy_j(&self.positions_m)
    }
}

/// Builds a periodic SPC-like water box.
///
/// `molecule_count` molecules sit on a simple cubic lattice in a cubic box of
/// side `box_length_m`. Velocities are drawn from a Maxwell-Boltzmann
/// distribution at `temperature_k`. The `seed` fixes the molecule orientations
/// and the velocities, so the case is reproducible.
///
/// The box must be large enough for a stable cutoff. The cutoff is 0.35 of the
/// box side and the skin is 0.10, so the list reach stays inside half the box.
pub fn spc_water_box(
    molecule_count: usize,
    box_length_m: f64,
    temperature_k: f64,
    seed: u64,
) -> Result<WaterBox, EngineError> {
    if molecule_count == 0 {
        return Err(EngineError::InvalidMoleculeCount { molecule_count });
    }
    if !box_length_m.is_finite() || box_length_m <= 0.0 {
        return Err(EngineError::InvalidPeriodicBox {
            index: 0,
            length_m: box_length_m,
        });
    }
    if !temperature_k.is_finite() || temperature_k < 0.0 {
        return Err(EngineError::InvalidTemperature { temperature_k });
    }

    let periodic_box = PeriodicBox::new([box_length_m, box_length_m, box_length_m])?;
    let cutoff_m = 0.35 * box_length_m;
    let switch_on_m = 0.9 * cutoff_m;
    let skin_m = 0.1 * box_length_m;
    let cutoff = Cutoff::new(cutoff_m, switch_on_m)?;

    let atom_count = 3 * molecule_count;
    let mut positions_m = vec![0.0; 3 * atom_count];
    let mut vdw_params = vec![VdwParams::new(0.0, 0.0, 0.0); atom_count];
    let mut charges_c = vec![0.0; atom_count];
    let mut masses_kg = vec![HYDROGEN_MASS_KG; atom_count];
    let mut bonds = BondStretchTerm::new();
    let mut triples = Vec::with_capacity(molecule_count);
    let mut angle_params = Vec::with_capacity(molecule_count);

    let mut rng = Rng::new(seed);
    let side = (molecule_count as f64).cbrt().ceil() as usize;
    let spacing_m = box_length_m / side as f64;
    let cos_angle = SPC_ANGLE_RAD.cos();
    let sin_angle = SPC_ANGLE_RAD.sin();

    for molecule in 0..molecule_count {
        let ix = molecule % side;
        let iy = (molecule / side) % side;
        let iz = (molecule / (side * side)) % side;
        let origin_m = [
            (ix as f64 + 0.5) * spacing_m,
            (iy as f64 + 0.5) * spacing_m,
            (iz as f64 + 0.5) * spacing_m,
        ];
        let phi = std::f64::consts::TAU * rng.next_f64();
        let (sin_phi, cos_phi) = phi.sin_cos();

        let oxygen = 3 * molecule;
        let hydrogen_one = oxygen + 1;
        let hydrogen_two = oxygen + 2;

        let local_m = [
            [0.0, 0.0, 0.0],
            [SPC_BOND_LENGTH_M, 0.0, 0.0],
            [
                SPC_BOND_LENGTH_M * cos_angle,
                SPC_BOND_LENGTH_M * sin_angle,
                0.0,
            ],
        ];
        for (site, local) in local_m.iter().enumerate() {
            let atom = oxygen + site;
            let base = 3 * atom;
            positions_m[base] = origin_m[0] + local[0] * cos_phi - local[1] * sin_phi;
            positions_m[base + 1] = origin_m[1] + local[0] * sin_phi + local[1] * cos_phi;
            positions_m[base + 2] = origin_m[2] + local[2];
        }

        masses_kg[oxygen] = OXYGEN_MASS_KG;
        vdw_params[oxygen] = VdwParams::new(SPC_VDW_A_J, SPC_VDW_B_PER_M, SPC_VDW_C_J_M6);
        charges_c[oxygen] = SPC_OXYGEN_CHARGE_C;
        charges_c[hydrogen_one] = SPC_HYDROGEN_CHARGE_C;
        charges_c[hydrogen_two] = SPC_HYDROGEN_CHARGE_C;

        bonds.add_bond(
            oxygen as u32,
            hydrogen_one as u32,
            SPC_BOND_K_N_PER_M,
            SPC_BOND_LENGTH_M,
        )?;
        bonds.add_bond(
            oxygen as u32,
            hydrogen_two as u32,
            SPC_BOND_K_N_PER_M,
            SPC_BOND_LENGTH_M,
        )?;
        triples.push([hydrogen_one as u32, oxygen as u32, hydrogen_two as u32]);
        angle_params.push(AngleBendParams::new(SPC_ANGLE_K_J_PER_RAD2, SPC_ANGLE_RAD));
    }

    let angles = AngleBendTerm::from_triples(&triples, &angle_params)?;
    let van_der_waals = VanDerWaalsTerm::from_params(&vdw_params, cutoff, periodic_box)?;
    let electrostatic = ElectrostaticTerm::from_charges(&charges_c, cutoff, periodic_box)?;

    let mut system = System::with_masses_kg(atom_count, &masses_kg)?;
    system.set_bond_stretch(bonds)?;
    system.set_angle_bend(angles)?;
    system.set_van_der_waals(van_der_waals)?;
    system.set_electrostatic(electrostatic)?;
    system.set_skin_m(skin_m)?;
    system.set_periodic_box(periodic_box)?;

    for molecule in 0..molecule_count {
        let oxygen = 3 * molecule;
        system.add_exclusion(oxygen as u32, (oxygen + 1) as u32)?;
        system.add_exclusion(oxygen as u32, (oxygen + 2) as u32)?;
        system.add_exclusion((oxygen + 1) as u32, (oxygen + 2) as u32)?;
    }

    let mut velocities_m_per_s = vec![0.0; 3 * atom_count];
    if temperature_k > 0.0 {
        for (atom, &mass_kg) in masses_kg.iter().enumerate() {
            let sigma_m_per_s = (BOLTZMANN_J_PER_K * temperature_k / mass_kg).sqrt();
            let base = 3 * atom;
            velocities_m_per_s[base] = sigma_m_per_s * rng.next_normal();
            velocities_m_per_s[base + 1] = sigma_m_per_s * rng.next_normal();
            velocities_m_per_s[base + 2] = sigma_m_per_s * rng.next_normal();
        }
    }

    Ok(WaterBox {
        system,
        positions_m,
        velocities_m_per_s,
        periodic_box,
        molecule_count,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::integrator::VelocityVerlet;

    fn box_length_m(molecule_count: usize) -> f64 {
        let volume_m3 = molecule_count as f64 * 33.4e-30;
        volume_m3.cbrt()
    }

    #[test]
    fn an_empty_or_invalid_box_is_rejected() {
        assert!(matches!(
            spc_water_box(0, 1.0e-9, 300.0, 1),
            Err(EngineError::InvalidMoleculeCount { .. })
        ));
        assert!(matches!(
            spc_water_box(10, -1.0, 300.0, 1),
            Err(EngineError::InvalidPeriodicBox { .. })
        ));
        assert!(matches!(
            spc_water_box(10, 1.0e-9, -1.0, 1),
            Err(EngineError::InvalidTemperature { .. })
        ));
    }

    #[test]
    fn the_water_box_has_three_atoms_per_molecule_and_is_periodic() {
        let water = spc_water_box(200, box_length_m(200), 300.0, 0x5A17E2).expect("valid box");
        assert_eq!(water.molecule_count(), 200);
        assert_eq!(water.atom_count(), 600);
        assert!(water.periodic_box().is_periodic());
        assert_eq!(water.positions_m().len(), 1800);
        assert_eq!(water.velocities_m_per_s().len(), 1800);
        assert_eq!(water.system().exclusion_count(), 3 * 200);
    }

    #[test]
    fn the_water_box_energy_is_finite_and_nonzero() {
        let mut water = spc_water_box(200, box_length_m(200), 300.0, 0x5A17E2).expect("valid box");
        let energy_j = water.potential_energy_j().expect("valid energy");
        let temperature_k = water.temperature_k().expect("valid temperature");
        println!(
            "water box: {energy_j:e} J potential, {temperature_k} K kinetic, {} neighbor pairs",
            water.system().neighbor_pair_count()
        );
        assert!(energy_j.is_finite());
        assert!(energy_j.abs() > 1.0e-20);
        assert!(water.system().neighbor_pair_count() > 0);
        assert!((temperature_k - 300.0).abs() / 300.0 < 0.1);
    }

    #[test]
    fn a_short_nve_run_keeps_the_energy_near_constant() {
        let molecule_count = 200;
        let water = spc_water_box(molecule_count, box_length_m(molecule_count), 300.0, 0x4E4E)
            .expect("valid box");
        let mut system = water.system().clone();
        let mut positions_m = water.positions_m().to_vec();
        let mut velocities_m_per_s = water.velocities_m_per_s().to_vec();
        let integrator = VelocityVerlet::new(1.0e-17).expect("valid step");

        let initial_energy_j = system
            .total_energy_j(&positions_m, &velocities_m_per_s)
            .expect("valid");
        let initial_potential_j = system.energy_j(&positions_m).expect("valid");
        let steps = 400;
        let mut max_drift_j = 0.0_f64;
        for _ in 0..steps {
            integrator
                .step(&mut system, &mut positions_m, &mut velocities_m_per_s)
                .expect("valid step");
            let energy_j = system
                .total_energy_j(&positions_m, &velocities_m_per_s)
                .expect("valid");
            max_drift_j = max_drift_j.max((energy_j - initial_energy_j).abs());
        }
        let relative_drift = max_drift_j / initial_energy_j.abs();
        println!(
            "water NVE: potential {initial_potential_j:e} J, total {initial_energy_j:e} J, max drift {max_drift_j:e} J, relative {relative_drift:e} over {steps} steps"
        );
        assert!(relative_drift < 1.0e-3, "relative drift {relative_drift:e}");
    }
}
