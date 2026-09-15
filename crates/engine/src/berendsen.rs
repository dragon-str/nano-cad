//! A Berendsen weak-coupling thermostat over a [`System`] and its velocities.
//!
//! One step is a velocity-Verlet step followed by a velocity rescale. The
//! scale factor is `lambda = sqrt(1 + (dt / tau) (T_target / T - 1))`.
//!
//! Berendsen is a relaxation thermostat. It drives the temperature to the
//! target, but it does not sample the canonical ensemble correctly: the
//! kinetic-energy fluctuations are too small. Do not use it to collect
//! equilibrium statistics. Use it to prepare a state or to remove drift.

use crate::error::EngineError;
use crate::integrator::VelocityVerlet;
use crate::system::System;

/// A Berendsen weak-coupling thermostat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BerendsenThermostat {
    integrator: VelocityVerlet,
    temperature_k: f64,
    tau_s: f64,
}

impl BerendsenThermostat {
    /// Builds a thermostat.
    ///
    /// The time step is in seconds, the temperature in kelvin, and the coupling
    /// time constant in seconds. The time step and the temperature must be
    /// positive and finite. The time constant must be positive and finite.
    pub fn new(dt_s: f64, temperature_k: f64, tau_s: f64) -> Result<Self, EngineError> {
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(EngineError::InvalidTemperature { temperature_k });
        }
        if !tau_s.is_finite() || tau_s <= 0.0 {
            return Err(EngineError::InvalidCouplingTime { tau_s });
        }
        let integrator = VelocityVerlet::new(dt_s)?;
        Ok(Self {
            integrator,
            temperature_k,
            tau_s,
        })
    }

    /// Returns the time step in seconds.
    pub fn dt_s(&self) -> f64 {
        self.integrator.dt_s()
    }

    /// Returns the target temperature in kelvin.
    pub fn temperature_k(&self) -> f64 {
        self.temperature_k
    }

    /// Returns the coupling time constant in seconds.
    pub fn tau_s(&self) -> f64 {
        self.tau_s
    }

    /// Advances the system by one step.
    ///
    /// Updates `positions_m` and `velocities_m_per_s` in place. A mass must be
    /// set on the system. An invalid input returns an error.
    pub fn step(
        &mut self,
        system: &mut System,
        positions_m: &mut [f64],
        velocities_m_per_s: &mut [f64],
    ) -> Result<(), EngineError> {
        self.integrator
            .step(system, positions_m, velocities_m_per_s)?;

        let current_k = system.temperature_k(velocities_m_per_s)?;
        if current_k <= 0.0 {
            return Ok(());
        }
        let lambda = (1.0 + self.dt_s() / self.tau_s * (self.temperature_k / current_k - 1.0))
            .max(0.0)
            .sqrt();
        for velocity in velocities_m_per_s.iter_mut() {
            *velocity *= lambda;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::CARBON_MASS_KG;

    fn free_system(atom_count: usize) -> System {
        let masses = vec![CARBON_MASS_KG; atom_count];
        System::with_masses_kg(atom_count, &masses).expect("valid masses")
    }

    #[test]
    fn invalid_settings_are_rejected() {
        assert!(matches!(
            BerendsenThermostat::new(1.0e-15, 0.0, 1.0e-13),
            Err(EngineError::InvalidTemperature { .. })
        ));
        assert!(matches!(
            BerendsenThermostat::new(1.0e-15, 300.0, 0.0),
            Err(EngineError::InvalidCouplingTime { .. })
        ));
        assert!(matches!(
            BerendsenThermostat::new(0.0, 300.0, 1.0e-13),
            Err(EngineError::NonPositiveTimestep { .. })
        ));
    }

    #[test]
    fn the_temperature_relaxes_to_the_target() {
        let atom_count = 40;
        let target_k = 300.0;
        let mut system = free_system(atom_count);
        let mut positions_m = vec![0.0; 3 * atom_count];
        let mut velocities_m_per_s = vec![0.0; 3 * atom_count];
        let cold_k = 10.0;
        let sigma_m_per_s = (crate::system::BOLTZMANN_J_PER_K * cold_k / CARBON_MASS_KG).sqrt();
        let mut rng = crate::rng::Rng::new(0xC01D);
        for velocity in velocities_m_per_s.iter_mut() {
            *velocity = sigma_m_per_s * rng.next_normal();
        }

        let dt_s = 5.0e-15;
        let tau_s = 2.0e-13;
        let mut thermostat =
            BerendsenThermostat::new(dt_s, target_k, tau_s).expect("valid thermostat");
        let start_k = system.temperature_k(&velocities_m_per_s).expect("valid");
        for _ in 0..4000 {
            thermostat
                .step(&mut system, &mut positions_m, &mut velocities_m_per_s)
                .expect("valid step");
        }
        let final_k = system.temperature_k(&velocities_m_per_s).expect("valid");
        let relative_error = (final_k - target_k).abs() / target_k;
        println!("Berendsen start {start_k} K, final {final_k} K, target {target_k} K");
        assert!(start_k < 30.0);
        assert!(relative_error < 0.02, "relative error {relative_error}");
    }
}
