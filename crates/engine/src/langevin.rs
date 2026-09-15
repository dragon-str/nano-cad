//! A Langevin thermostat over a [`System`] and its velocities.
//!
//! One step is a velocity-Verlet step followed by an exact
//! Ornstein-Uhlenbeck update of every velocity. The update is
//! `v <- c1 v + c2 sqrt(k_B T / m) xi`, with `c1 = exp(-gamma dt)`,
//! `c2 = sqrt(1 - c1^2)`, and `xi` a standard normal value. The friction is
//! `gamma`. The noise is seeded, so a run is reproducible.
//!
//! The stationary velocity distribution is Maxwell-Boltzmann at the target
//! temperature.

use crate::error::EngineError;
use crate::integrator::VelocityVerlet;
use crate::rng::Rng;
use crate::system::{System, BOLTZMANN_J_PER_K};

/// A Langevin thermostat.
#[derive(Clone, Debug, PartialEq)]
pub struct LangevinThermostat {
    integrator: VelocityVerlet,
    temperature_k: f64,
    friction_per_s: f64,
    rng: Rng,
}

impl LangevinThermostat {
    /// Builds a thermostat.
    ///
    /// The time step is in seconds, the temperature in kelvin, and the friction
    /// in reciprocal seconds. The time step and the temperature must be
    /// positive and finite. The friction must be non-negative and finite. The
    /// seed fixes the noise.
    pub fn new(
        dt_s: f64,
        temperature_k: f64,
        friction_per_s: f64,
        seed: u64,
    ) -> Result<Self, EngineError> {
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(EngineError::InvalidTemperature { temperature_k });
        }
        if !friction_per_s.is_finite() || friction_per_s < 0.0 {
            return Err(EngineError::InvalidFriction { friction_per_s });
        }
        let integrator = VelocityVerlet::new(dt_s)?;
        Ok(Self {
            integrator,
            temperature_k,
            friction_per_s,
            rng: Rng::new(seed),
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

    /// Returns the friction in reciprocal seconds.
    pub fn friction_per_s(&self) -> f64 {
        self.friction_per_s
    }

    /// Returns the current random-generator state, for reproducibility checks.
    pub fn rng_state(&self) -> u64 {
        self.rng.state()
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

        let masses_kg = system.masses_kg();
        let dt_s = self.dt_s();
        let decay = (-self.friction_per_s * dt_s).exp();
        let noise_scale = (1.0 - decay * decay).sqrt();
        let temperature_k = self.temperature_k;
        for (atom, &mass_kg) in masses_kg.iter().enumerate() {
            let sigma_m_per_s = noise_scale * (BOLTZMANN_J_PER_K * temperature_k / mass_kg).sqrt();
            let base = atom * 3;
            for axis in 0..3 {
                let index = base + axis;
                let noise_m_per_s = sigma_m_per_s * self.rng.next_normal();
                velocities_m_per_s[index] = decay * velocities_m_per_s[index] + noise_m_per_s;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::{Rng as TestRng, CARBON_MASS_KG};

    fn free_system(atom_count: usize) -> System {
        let masses = vec![CARBON_MASS_KG; atom_count];
        System::with_masses_kg(atom_count, &masses).expect("valid masses")
    }

    fn velocity_spread_m_per_s(temperature_k: f64, mass_kg: f64) -> f64 {
        (BOLTZMANN_J_PER_K * temperature_k / mass_kg).sqrt()
    }

    #[test]
    fn invalid_settings_are_rejected() {
        assert!(matches!(
            LangevinThermostat::new(1.0e-15, 0.0, 1.0e12, 1),
            Err(EngineError::InvalidTemperature { .. })
        ));
        assert!(matches!(
            LangevinThermostat::new(1.0e-15, 300.0, -1.0, 1),
            Err(EngineError::InvalidFriction { .. })
        ));
        assert!(matches!(
            LangevinThermostat::new(0.0, 300.0, 1.0e12, 1),
            Err(EngineError::NonPositiveTimestep { .. })
        ));
    }

    #[test]
    fn the_same_seed_gives_the_same_trajectory() {
        let mut first = free_system(4);
        let mut second = free_system(4);
        let mut first_positions = vec![0.0; 12];
        let mut second_positions = vec![0.0; 12];
        let mut first_velocities = vec![1.0; 12];
        let mut second_velocities = vec![1.0; 12];
        let mut a = LangevinThermostat::new(1.0e-15, 300.0, 1.0e13, 42).expect("valid");
        let mut b = LangevinThermostat::new(1.0e-15, 300.0, 1.0e13, 42).expect("valid");
        for _ in 0..100 {
            a.step(&mut first, &mut first_positions, &mut first_velocities)
                .expect("valid step");
            b.step(&mut second, &mut second_positions, &mut second_velocities)
                .expect("valid step");
        }
        assert_eq!(first_velocities, second_velocities);
    }

    #[test]
    fn the_long_run_temperature_matches_the_target() {
        let atom_count = 40;
        let target_k = 300.0;
        let mut system = free_system(atom_count);
        let mut positions_m = vec![0.0; 3 * atom_count];
        let mut velocities_m_per_s = vec![0.0; 3 * atom_count];
        let mut thermostat =
            LangevinThermostat::new(5.0e-15, target_k, 5.0e12, 0xA11CE).expect("valid");

        let equilibration = 2000;
        for _ in 0..equilibration {
            thermostat
                .step(&mut system, &mut positions_m, &mut velocities_m_per_s)
                .expect("valid step");
        }

        let samples = 20_000;
        let mut temperature_sum_k = 0.0;
        for _ in 0..samples {
            thermostat
                .step(&mut system, &mut positions_m, &mut velocities_m_per_s)
                .expect("valid step");
            temperature_sum_k += system
                .temperature_k(&velocities_m_per_s)
                .expect("valid temperature");
        }
        let mean_k = temperature_sum_k / samples as f64;
        let relative_error = (mean_k - target_k).abs() / target_k;
        println!("Langevin mean temperature {mean_k} K, target {target_k} K");
        assert!(relative_error < 0.03, "relative error {relative_error}");
    }

    #[test]
    fn the_velocity_variance_is_boltzmann_like() {
        let atom_count = 40;
        let target_k = 300.0;
        let mut system = free_system(atom_count);
        let mut positions_m = vec![0.0; 3 * atom_count];
        let mut velocities_m_per_s = vec![0.0; 3 * atom_count];
        let mut thermostat =
            LangevinThermostat::new(5.0e-15, target_k, 5.0e12, 0xBEE5).expect("valid");

        for _ in 0..2000 {
            thermostat
                .step(&mut system, &mut positions_m, &mut velocities_m_per_s)
                .expect("valid step");
        }

        let samples = 20_000;
        let mut sum = 0.0;
        let mut sum_sq = 0.0;
        for _ in 0..samples {
            thermostat
                .step(&mut system, &mut positions_m, &mut velocities_m_per_s)
                .expect("valid step");
            for atom in 0..atom_count {
                let vx = velocities_m_per_s[3 * atom];
                sum += vx;
                sum_sq += vx * vx;
            }
        }
        let count = (samples * atom_count) as f64;
        let mean = sum / count;
        let variance = sum_sq / count - mean * mean;
        let expected_m2_per_s2 = BOLTZMANN_J_PER_K * target_k / CARBON_MASS_KG;
        let relative_error = (variance - expected_m2_per_s2).abs() / expected_m2_per_s2;
        println!(
            "Langevin vx variance {variance:e} m^2/s^2, expected {expected_m2_per_s2:e} m^2/s^2"
        );
        assert!(relative_error < 0.05, "relative error {relative_error}");
    }

    #[test]
    fn a_zero_friction_thermostat_conserves_energy() {
        let atom_count = 8;
        let mut system = free_system(atom_count);
        let mut rng = TestRng::new(0x5EED);
        let mut positions_m = Vec::with_capacity(3 * atom_count);
        let mut velocities_m_per_s = Vec::with_capacity(3 * atom_count);
        let sigma_m_per_s = velocity_spread_m_per_s(300.0, CARBON_MASS_KG);
        for _ in 0..3 * atom_count {
            positions_m.push(0.0);
            velocities_m_per_s.push(sigma_m_per_s * (rng.next_f64() - 0.5));
        }
        let initial_energy_j = system
            .total_energy_j(&positions_m, &velocities_m_per_s)
            .expect("valid");
        let mut thermostat = LangevinThermostat::new(5.0e-15, 300.0, 0.0, 7).expect("valid");
        for _ in 0..1000 {
            thermostat
                .step(&mut system, &mut positions_m, &mut velocities_m_per_s)
                .expect("valid step");
        }
        let final_energy_j = system
            .total_energy_j(&positions_m, &velocities_m_per_s)
            .expect("valid");
        let relative_drift = (final_energy_j - initial_energy_j).abs() / initial_energy_j.abs();
        assert!(relative_drift < 1.0e-9, "relative drift {relative_drift}");
    }
}
