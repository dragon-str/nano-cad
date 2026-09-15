use crate::error::EngineError;
use crate::system::System;

/// The velocity Verlet integrator in SI units.
///
/// Positions are in metres, velocities in metres per second, and the time step
/// in seconds. The masses live on the [`System`], in kilograms.
///
/// One step does a half velocity kick, a position drift, a new force
/// evaluation, and a second half kick. This is symplectic, so the energy
/// oscillates but does not drift on average.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VelocityVerlet {
    dt_s: f64,
}

impl VelocityVerlet {
    /// Builds an integrator for a time step in seconds.
    ///
    /// The time step must be finite and positive. An invalid value returns an
    /// error and never panics.
    pub fn new(dt_s: f64) -> Result<Self, EngineError> {
        if !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(EngineError::NonPositiveTimestep { dt_s });
        }
        Ok(Self { dt_s })
    }

    /// Returns the time step in seconds.
    pub fn dt_s(&self) -> f64 {
        self.dt_s
    }

    /// Advances the system by one step.
    ///
    /// Updates `positions_m` and `velocities_m_per_s` in place. A mass must be
    /// set on the system. Both buffers must hold three values per atom and be
    /// finite. This returns an error instead of panicking.
    pub fn step(
        &self,
        system: &mut System,
        positions_m: &mut [f64],
        velocities_m_per_s: &mut [f64],
    ) -> Result<(), EngineError> {
        let atom_count = system.atom_count();
        let expected = 3 * atom_count;
        if positions_m.len() != expected {
            return Err(EngineError::BufferSizeMismatch {
                len: positions_m.len(),
                expected,
            });
        }
        if velocities_m_per_s.len() != expected {
            return Err(EngineError::BufferSizeMismatch {
                len: velocities_m_per_s.len(),
                expected,
            });
        }
        if !system.has_masses() {
            return Err(EngineError::MassesNotSet);
        }
        for (index, value) in positions_m.iter().enumerate() {
            if !value.is_finite() {
                return Err(EngineError::NonFinitePosition { index });
            }
        }
        for (index, value) in velocities_m_per_s.iter().enumerate() {
            if !value.is_finite() {
                return Err(EngineError::NonFiniteVelocity { index });
            }
        }

        let masses_kg = system.masses_kg().to_vec();
        let dt_s = self.dt_s;
        let half_dt_s = 0.5 * dt_s;

        let forces_n = system.forces_n(positions_m)?;
        for (atom, &mass_kg) in masses_kg.iter().enumerate() {
            let base = atom * 3;
            for axis in 0..3 {
                let index = base + axis;
                let half_velocity =
                    velocities_m_per_s[index] + half_dt_s * forces_n[index] / mass_kg;
                positions_m[index] += half_velocity * dt_s;
                velocities_m_per_s[index] = half_velocity;
            }
        }

        let forces_n = system.forces_n(positions_m)?;
        for (atom, &mass_kg) in masses_kg.iter().enumerate() {
            let base = atom * 3;
            for axis in 0..3 {
                let index = base + axis;
                velocities_m_per_s[index] += half_dt_s * forces_n[index] / mass_kg;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::small_molecule_positions_m;

    #[test]
    fn a_non_positive_or_non_finite_time_step_is_rejected() {
        assert!(matches!(
            VelocityVerlet::new(0.0),
            Err(EngineError::NonPositiveTimestep { .. })
        ));
        assert!(matches!(
            VelocityVerlet::new(-1.0),
            Err(EngineError::NonPositiveTimestep { .. })
        ));
        assert!(matches!(
            VelocityVerlet::new(f64::NAN),
            Err(EngineError::NonPositiveTimestep { .. })
        ));
    }

    #[test]
    fn a_step_without_masses_is_rejected() {
        let mut system = System::new(2);
        let mut positions_m = vec![0.0; 6];
        let mut velocities_m_per_s = vec![0.0; 6];
        let integrator = VelocityVerlet::new(1.0e-15).expect("valid step");
        assert!(matches!(
            integrator.step(&mut system, &mut positions_m, &mut velocities_m_per_s),
            Err(EngineError::MassesNotSet)
        ));
    }

    #[test]
    fn a_buffer_of_the_wrong_size_is_rejected() {
        let positions_m = small_molecule_positions_m();
        let masses = vec![crate::test_support::CARBON_MASS_KG; 6];
        let mut system = crate::test_support::bonded_system(&positions_m, &masses);
        let integrator = VelocityVerlet::new(1.0e-15).expect("valid step");
        let mut short_positions_m = vec![0.0; 3];
        let mut short_velocities_m_per_s = vec![0.0; 3];
        assert!(matches!(
            integrator.step(
                &mut system,
                &mut short_positions_m,
                &mut short_velocities_m_per_s
            ),
            Err(EngineError::BufferSizeMismatch { .. })
        ));
    }

    #[test]
    fn a_free_atom_moves_at_constant_velocity() {
        let mut system = System::with_masses_kg(1, &[1.0e-26]).expect("valid mass");
        let integrator = VelocityVerlet::new(1.0e-12).expect("valid step");
        let mut positions_m = vec![0.0, 0.0, 0.0];
        let mut velocities_m_per_s = vec![1.0, 0.0, 0.0];
        for _ in 0..10 {
            integrator
                .step(&mut system, &mut positions_m, &mut velocities_m_per_s)
                .expect("valid step");
        }
        assert!((positions_m[0] - 10.0e-12).abs() < 1.0e-24);
        assert!((velocities_m_per_s[0] - 1.0).abs() < 1.0e-15);
    }

    #[test]
    fn the_nve_total_energy_drift_stays_below_the_bound() {
        let positions_m = small_molecule_positions_m();
        let masses = vec![crate::test_support::CARBON_MASS_KG; 6];
        let mut system = crate::test_support::bonded_system(&positions_m, &masses);
        let mut positions_m = positions_m;
        let mut velocities_m_per_s = vec![
            120.0, -80.0, 40.0, -60.0, 100.0, -30.0, 90.0, 50.0, -70.0, -110.0, 20.0, 60.0, 70.0,
            -40.0, 30.0, -50.0, 80.0, -20.0,
        ];
        let integrator = VelocityVerlet::new(5.0e-17).expect("valid step");

        let initial_energy_j = system
            .total_energy_j(&positions_m, &velocities_m_per_s)
            .expect("valid");
        let mut max_drift_j = 0.0_f64;
        let steps = 4_000;
        for _ in 0..steps {
            integrator
                .step(&mut system, &mut positions_m, &mut velocities_m_per_s)
                .expect("valid step");
            let energy_j = system
                .total_energy_j(&positions_m, &velocities_m_per_s)
                .expect("valid");
            let drift_j = (energy_j - initial_energy_j).abs();
            if drift_j > max_drift_j {
                max_drift_j = drift_j;
            }
        }
        let relative_drift = max_drift_j / initial_energy_j.abs();
        println!(
            "NVE drift: initial {initial_energy_j:e} J, max absolute {max_drift_j:e} J, relative {relative_drift:e} over {steps} steps"
        );

        assert!(initial_energy_j.abs() > 1.0e-24);
        assert!(relative_drift < 1.0e-4, "relative drift {relative_drift:e}");
    }
}
