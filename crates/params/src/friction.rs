//! Friction extraction from a driven sliding contact.
//!
//! This implements step 3 of the extraction pipeline in `PARAMETERS.md`: slide
//! one part over another at a set load and speed, then record the mean and the
//! spread. The lateral force is the force in the drive spring. Its magnitude
//! gives the friction force, and its spread gives the uncertainty.
//!
//! The result is an estimate from a simulation. It is not a certified material
//! constant.

use nanocad_engine::{EngineError, System};

use crate::error::ParamError;
use crate::provenance::{extraction_provenance, Method, Provenance};
use crate::quantity::Quantity;
use crate::stats::{mean, sample_std_dev};

/// The driven-contact setup for a friction extraction.
///
/// The slider atom is pulled along x by a spring to a stage that moves at
/// `stage_velocity_m_per_s`. A constant `normal_load_n` pushes it down onto the
/// substrate. A drag `drag_n_s_per_m` acts on its velocity. Only the slider
/// moves. Every other atom is fixed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrictionConfig {
    /// The atom index of the sliding body.
    pub slider_atom: u32,
    /// The normal load, in newtons. It must be positive.
    pub normal_load_n: f64,
    /// The viscous drag on the slider, in newton seconds per metre.
    pub drag_n_s_per_m: f64,
    /// The lateral stage speed, in metres per second.
    pub stage_velocity_m_per_s: f64,
    /// The drive-spring stiffness, in newtons per metre.
    pub drive_stiffness_n_per_m: f64,
    /// The integration time step, in seconds.
    pub dt_s: f64,
    /// The number of steps before the measurement window.
    pub equilibration_steps: usize,
    /// The number of measured steps. It must be at least two.
    pub measurement_steps: usize,
}

impl Default for FrictionConfig {
    fn default() -> Self {
        Self {
            slider_atom: 0,
            normal_load_n: 1.0e-11,
            drag_n_s_per_m: 2.0e-13,
            stage_velocity_m_per_s: 1.0,
            drive_stiffness_n_per_m: 1.0e-1,
            dt_s: 1.0e-14,
            equilibration_steps: 10_000,
            measurement_steps: 50_000,
        }
    }
}

/// The friction summary of a sliding run.
#[derive(Clone, Debug, PartialEq)]
pub struct FrictionResult {
    /// The friction coefficient, dimensionless, with the spread as uncertainty.
    pub friction_coefficient: Quantity,
    /// The mean magnitude of the lateral force, in newtons.
    pub friction_force_mean_n: f64,
    /// The sample standard deviation of the lateral-force magnitude, in newtons.
    pub friction_force_spread_n: f64,
    /// The normal load, in newtons.
    pub normal_load_n: f64,
    /// The number of lateral-force samples.
    pub samples: usize,
    /// The provenance of the extraction.
    pub provenance: Provenance,
}

/// Summarizes signed lateral-force samples into a friction result.
///
/// The magnitudes of the samples give the friction force. The mean magnitude
/// divided by the load gives the coefficient. The sample standard deviation of
/// the magnitudes is the spread, and the spread divided by the load is the
/// uncertainty on the coefficient.
pub fn summarize_friction(
    lateral_forces_n: &[f64],
    normal_load_n: f64,
) -> Result<FrictionResult, ParamError> {
    let moments = friction_moments(lateral_forces_n, normal_load_n)?;
    let notes = format!(
        "summary of {} lateral-force samples at a normal load of {:e} N; \
         the mean magnitude is the friction force",
        moments.samples, moments.normal_load_n
    );
    moments.into_result(extraction_provenance(Method::Fit, notes))
}

/// Runs a driven sliding contact and extracts a friction coefficient.
///
/// The function integrates the slider atom with a velocity-Verlet step that
/// adds three external forces to the engine force: the drive spring, the normal
/// load, and the drag. Every other atom is held fixed. It discards the
/// equilibration window, then measures the drive-spring force for the
/// measurement window. See [`summarize_friction`] for the summary.
///
/// The engine has no public integrator that accepts an external force, so this
/// module carries its own small step. It uses the engine for the interaction
/// force only.
pub fn extract_friction(
    system: &mut System,
    positions_m: &[f64],
    config: &FrictionConfig,
) -> Result<FrictionResult, ParamError> {
    validate_config(system, positions_m, config)?;
    let lateral_forces_n = simulate_sliding(system, positions_m, config)?;
    let moments = friction_moments(&lateral_forces_n, config.normal_load_n)?;
    let notes = format!(
        "driven sliding contact over {} measured steps of {:e} s at {:e} m/s; \
         drag {:e} N*s/m; {} lateral-force samples",
        config.measurement_steps,
        config.dt_s,
        config.stage_velocity_m_per_s,
        config.drag_n_s_per_m,
        moments.samples
    );
    moments.into_result(extraction_provenance(Method::Md, notes))
}

struct FrictionMoments {
    friction_coefficient: f64,
    spread_coefficient: f64,
    mean_force_n: f64,
    spread_force_n: f64,
    normal_load_n: f64,
    samples: usize,
}

impl FrictionMoments {
    fn into_result(self, provenance: Provenance) -> Result<FrictionResult, ParamError> {
        let friction_coefficient = Quantity::derived(self.friction_coefficient, "1")?
            .with_uncertainty_si(self.spread_coefficient);
        Ok(FrictionResult {
            friction_coefficient,
            friction_force_mean_n: self.mean_force_n,
            friction_force_spread_n: self.spread_force_n,
            normal_load_n: self.normal_load_n,
            samples: self.samples,
            provenance,
        })
    }
}

fn friction_moments(
    lateral_forces_n: &[f64],
    normal_load_n: f64,
) -> Result<FrictionMoments, ParamError> {
    if !normal_load_n.is_finite() || normal_load_n <= 0.0 {
        return Err(ParamError::InvalidExtraction {
            reason: "the normal load must be finite and positive".to_string(),
        });
    }
    if lateral_forces_n.len() < 2 {
        return Err(ParamError::InvalidExtraction {
            reason: "at least two lateral-force samples are needed".to_string(),
        });
    }
    if !lateral_forces_n.iter().all(|force| force.is_finite()) {
        return Err(ParamError::InvalidExtraction {
            reason: "a lateral-force sample is not finite".to_string(),
        });
    }
    let magnitudes_n: Vec<f64> = lateral_forces_n.iter().map(|force| force.abs()).collect();
    let mean_force_n = mean(&magnitudes_n).ok_or_else(|| ParamError::InvalidExtraction {
        reason: "the lateral-force samples are empty".to_string(),
    })?;
    let spread_force_n =
        sample_std_dev(&magnitudes_n).ok_or_else(|| ParamError::InvalidExtraction {
            reason: "the lateral-force spread needs at least two samples".to_string(),
        })?;
    Ok(FrictionMoments {
        friction_coefficient: mean_force_n / normal_load_n,
        spread_coefficient: spread_force_n / normal_load_n,
        mean_force_n,
        spread_force_n,
        normal_load_n,
        samples: lateral_forces_n.len(),
    })
}

fn simulate_sliding(
    system: &mut System,
    base_positions_m: &[f64],
    config: &FrictionConfig,
) -> Result<Vec<f64>, ParamError> {
    let slider = config.slider_atom as usize;
    let mass_kg = system.masses_kg()[slider];
    let mut positions_m = base_positions_m.to_vec();
    let mut velocities_m_per_s = vec![0.0; positions_m.len()];
    let start_x_m = positions_m[3 * slider];
    let half_dt_s = 0.5 * config.dt_s;
    let total_steps = config.equilibration_steps + config.measurement_steps;
    let mut time_s = 0.0;
    let mut samples_n = Vec::with_capacity(config.measurement_steps);

    let mut force_n = slider_force(
        system,
        &positions_m,
        &velocities_m_per_s,
        config,
        time_s,
        start_x_m,
    )?;
    for step in 0..total_steps {
        for (axis, force) in force_n.iter().enumerate() {
            let index = 3 * slider + axis;
            let half_velocity = velocities_m_per_s[index] + half_dt_s * force / mass_kg;
            positions_m[index] += half_velocity * config.dt_s;
            velocities_m_per_s[index] = half_velocity;
        }
        time_s += config.dt_s;
        force_n = slider_force(
            system,
            &positions_m,
            &velocities_m_per_s,
            config,
            time_s,
            start_x_m,
        )?;
        for (axis, force) in force_n.iter().enumerate() {
            let index = 3 * slider + axis;
            velocities_m_per_s[index] += half_dt_s * force / mass_kg;
        }
        if step >= config.equilibration_steps {
            let stage_x_m = start_x_m + config.stage_velocity_m_per_s * time_s;
            let spring_force_n =
                config.drive_stiffness_n_per_m * (stage_x_m - positions_m[3 * slider]);
            samples_n.push(spring_force_n);
        }
    }
    if !samples_n.iter().all(|force| force.is_finite()) {
        return Err(ParamError::InvalidExtraction {
            reason: "the sliding contact became unstable".to_string(),
        });
    }
    Ok(samples_n)
}

fn slider_force(
    system: &mut System,
    positions_m: &[f64],
    velocities_m_per_s: &[f64],
    config: &FrictionConfig,
    time_s: f64,
    start_x_m: f64,
) -> Result<[f64; 3], ParamError> {
    let slider = config.slider_atom as usize;
    let engine_forces_n = system.forces_n(positions_m)?;
    let base = 3 * slider;
    let stage_x_m = start_x_m + config.stage_velocity_m_per_s * time_s;
    let drag = config.drag_n_s_per_m;
    let spring_n = config.drive_stiffness_n_per_m * (stage_x_m - positions_m[base]);
    let mut force_n = [0.0; 3];
    force_n[0] = engine_forces_n[base] + spring_n - drag * velocities_m_per_s[base];
    force_n[1] = engine_forces_n[base + 1] - drag * velocities_m_per_s[base + 1];
    force_n[2] =
        engine_forces_n[base + 2] - config.normal_load_n - drag * velocities_m_per_s[base + 2];
    Ok(force_n)
}

fn validate_config(
    system: &System,
    positions_m: &[f64],
    config: &FrictionConfig,
) -> Result<(), ParamError> {
    if positions_m.len() != 3 * system.atom_count() {
        return Err(EngineError::BufferSizeMismatch {
            len: positions_m.len(),
            expected: 3 * system.atom_count(),
        }
        .into());
    }
    if config.slider_atom as usize >= system.atom_count() {
        return Err(ParamError::InvalidExtraction {
            reason: "the slider atom is outside the atom count".to_string(),
        });
    }
    if !system.has_masses() {
        return Err(EngineError::MassesNotSet.into());
    }
    if !config.normal_load_n.is_finite() || config.normal_load_n <= 0.0 {
        return Err(ParamError::InvalidExtraction {
            reason: "the normal load must be finite and positive".to_string(),
        });
    }
    if !config.drag_n_s_per_m.is_finite() || config.drag_n_s_per_m <= 0.0 {
        return Err(ParamError::InvalidExtraction {
            reason: "the drag must be finite and positive".to_string(),
        });
    }
    if !config.stage_velocity_m_per_s.is_finite() || config.stage_velocity_m_per_s <= 0.0 {
        return Err(ParamError::InvalidExtraction {
            reason: "the stage velocity must be finite and positive".to_string(),
        });
    }
    if !config.drive_stiffness_n_per_m.is_finite() || config.drive_stiffness_n_per_m <= 0.0 {
        return Err(ParamError::InvalidExtraction {
            reason: "the drive stiffness must be finite and positive".to_string(),
        });
    }
    if !config.dt_s.is_finite() || config.dt_s <= 0.0 {
        return Err(ParamError::InvalidExtraction {
            reason: "the time step must be finite and positive".to_string(),
        });
    }
    if config.measurement_steps < 2 {
        return Err(ParamError::InvalidExtraction {
            reason: "at least two measurement steps are needed".to_string(),
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::{sliding_system, SLIDER_ATOM};
    use std::f64::consts::PI;

    #[test]
    fn a_sinusoidal_force_profile_gives_the_analytic_friction() {
        let amplitude_n = 3.0e-11;
        let period_m = 2.0e-10;
        let samples: Vec<f64> = (0..2000)
            .map(|index| {
                let x_m = period_m * index as f64 / 2000.0;
                amplitude_n * (2.0 * PI * x_m / period_m).sin()
            })
            .collect();
        let load_n = 1.0e-10;
        let result = summarize_friction(&samples, load_n).expect("valid samples");
        let expected_mean_n = 2.0 * amplitude_n / PI;
        let expected_spread_n = amplitude_n * (0.5 - 4.0 / (PI * PI)).sqrt();
        println!(
            "friction profile: mean {:e} N, spread {:e} N, coefficient {:e}",
            result.friction_force_mean_n,
            result.friction_force_spread_n,
            result.friction_coefficient.value_si()
        );
        assert!((result.friction_force_mean_n - expected_mean_n).abs() / expected_mean_n < 1.0e-3);
        assert!(
            (result.friction_force_spread_n - expected_spread_n).abs() / expected_spread_n < 1.0e-2
        );
        let expected_coefficient = expected_mean_n / load_n;
        assert!(
            (result.friction_coefficient.value_si() - expected_coefficient).abs()
                / expected_coefficient
                < 1.0e-3
        );
    }

    #[test]
    fn the_sliding_contact_repeats_and_reports_a_spread() {
        let config = FrictionConfig {
            slider_atom: SLIDER_ATOM,
            equilibration_steps: 10_000,
            measurement_steps: 30_000,
            ..FrictionConfig::default()
        };
        let (mut system, positions_m) = sliding_system();
        let first = extract_friction(&mut system, &positions_m, &config).expect("valid run");
        let (mut system, positions_m) = sliding_system();
        let second = extract_friction(&mut system, &positions_m, &config).expect("valid run");
        println!(
            "sliding: coefficient {:e}, mean {:e} N, spread {:e} N, samples {}",
            first.friction_coefficient.value_si(),
            first.friction_force_mean_n,
            first.friction_force_spread_n,
            first.samples
        );
        assert_eq!(first.samples, second.samples);
        assert_eq!(first.normal_load_n, second.normal_load_n);
        assert_eq!(first.friction_force_mean_n, second.friction_force_mean_n);
        assert_eq!(
            first.friction_force_spread_n,
            second.friction_force_spread_n
        );
        assert_eq!(
            first.friction_coefficient.value_si(),
            second.friction_coefficient.value_si()
        );
        assert!(first.friction_force_mean_n > 0.0);
        assert!(first.friction_force_mean_n.is_finite());
        assert!(first.friction_force_spread_n > 0.0);
        assert!(first.friction_coefficient.value_si() > 0.0);
        let expected_coefficient = first.friction_force_mean_n / config.normal_load_n;
        assert!(
            (first.friction_coefficient.value_si() - expected_coefficient).abs()
                <= 1.0e-15 * expected_coefficient.abs()
        );
    }

    #[test]
    fn an_invalid_friction_setup_is_rejected() {
        let (mut system, positions_m) = sliding_system();
        let zero_load = FrictionConfig {
            slider_atom: SLIDER_ATOM,
            normal_load_n: 0.0,
            ..FrictionConfig::default()
        };
        assert!(matches!(
            extract_friction(&mut system, &positions_m, &zero_load),
            Err(ParamError::InvalidExtraction { .. })
        ));
        let bad_atom = FrictionConfig {
            slider_atom: 999,
            ..FrictionConfig::default()
        };
        assert!(matches!(
            extract_friction(&mut system, &positions_m, &bad_atom),
            Err(ParamError::InvalidExtraction { .. })
        ));
    }

    #[test]
    fn a_single_sample_cannot_give_a_spread() {
        assert!(matches!(
            summarize_friction(&[1.0e-11], 1.0e-10),
            Err(ParamError::InvalidExtraction { .. })
        ));
    }
}
