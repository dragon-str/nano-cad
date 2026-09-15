//! L1 atomistic engine: pair lists, forces, integrators.
//!
//! The engine reads positions from `nanocad-model` and computes interaction
//! terms. Positions are in SI metres, energies in joules, and forces in
//! newtons.
#![forbid(unsafe_code)]

mod angle_bend;
mod berendsen;
mod bond_stretch;
mod electrostatic;
mod error;
mod geometry;
mod integrator;
mod langevin;
mod minimize;
mod neighbor;
mod nonbonded;
mod out_of_plane;
mod pair_kernel;
mod reference;
mod rng;
mod system;
mod torsion;
mod van_der_waals;

#[cfg(test)]
mod test_support;

pub use angle_bend::{AngleBendParams, AngleBendTerm};
pub use berendsen::BerendsenThermostat;
pub use bond_stretch::{BondInfo, BondStretchParams, BondStretchTerm};
pub use electrostatic::{ElectrostaticTerm, COULOMB_CONSTANT_N_M2_PER_C2};
pub use error::EngineError;
pub use integrator::VelocityVerlet;
pub use langevin::LangevinThermostat;
pub use minimize::{minimize, minimize_with, MinimizeMethod, MinimizeOptions, MinimizeResult};
pub use neighbor::VerletList;
pub use nonbonded::{Cutoff, PeriodicBox};
pub use out_of_plane::{OutOfPlaneParams, OutOfPlaneTerm};
pub use reference::{
    spc_water_box, WaterBox, ELEMENTARY_CHARGE_C, HYDROGEN_MASS_KG, OXYGEN_MASS_KG,
    SPC_ANGLE_K_J_PER_RAD2, SPC_ANGLE_RAD, SPC_BOND_K_N_PER_M, SPC_BOND_LENGTH_M,
    SPC_HYDROGEN_CHARGE_C, SPC_OXYGEN_CHARGE_C, SPC_VDW_A_J, SPC_VDW_B_PER_M, SPC_VDW_C_J_M6,
};
pub use rng::Rng;
pub use system::{System, BOLTZMANN_J_PER_K};
pub use torsion::{TorsionParams, TorsionTerm};
pub use van_der_waals::{VanDerWaalsTerm, VdwParams};

/// Returns the crate version string.
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_is_not_empty() {
        assert!(!version().is_empty());
    }
}
