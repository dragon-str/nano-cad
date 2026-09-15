//! Thermal-property extraction under NVT.
//!
//! This implements step 5 of the extraction pipeline in `PARAMETERS.md`: run
//! NVT dynamics and extract the specific heat and the thermal conductivity.
//!
//! # Specific heat
//!
//! The estimator is the kinetic-energy fluctuation formula
//!
//! ```text
//! C_v = ( <E_k^2> - <E_k>^2 ) / ( k_B T^2 )
//! ```
//!
//! divided by the total mass. In the canonical ensemble the kinetic energy is a
//! sum of `3 N` quadratic velocity terms, so this estimator recovers the
//! translational heat capacity `(3/2) N k_B` for any system. It does not
//! include the configurational (potential-energy) contribution, which also
//! needs `Var(U)`. The result is therefore a kinetic estimate of the specific
//! heat. For a monatomic ideal gas it is the full classical value. This is an
//! estimate from classical equipartition, not a measured material constant.
//!
//! # Thermal conductivity
//!
//! The estimator is a non-equilibrium (direct) method. Two Langevin baths hold
//! the two ends of the sample at different temperatures. In the steady state
//! the power that the hot bath adds equals the heat current along the sample.
//! The conductivity is that current divided by the measured temperature
//! gradient and the stated cross-sectional area. A full Green-Kubo heat-flux
//! operator needs the per-term virial, which the engine does not expose, so
//! this direct method replaces it.
//!
//! Both results are estimates from a simulation. They are not certified
//! material constants. See [`extract_thermal`] for the model limits.

use nanocad_engine::{
    BerendsenThermostat, EngineError, LangevinThermostat, Rng, System, VelocityVerlet,
    BOLTZMANN_J_PER_K,
};

use crate::error::ParamError;
use crate::provenance::{extraction_provenance, Method, Provenance};
use crate::quantity::Quantity;
use crate::stats::{mean, sample_std_dev};

/// The fewest samples any thermal estimator accepts.
///
/// The block uncertainty needs several blocks. This is the floor, not a
/// recommendation.
const MIN_THERMAL_SAMPLES: usize = 32;

/// The number of contiguous blocks used for the uncertainty of a thermal
/// estimate. Block averaging absorbs the correlation between nearby samples.
const THERMAL_BLOCK_COUNT: usize = 16;

/// The NVT setup and the two-bath setup for a thermal extraction.
///
/// The specific-heat run uses a Berendsen pre-equilibration, then a Langevin
/// NVT production. A Berendsen thermostat does not sample the canonical
/// ensemble, so it prepares the state only. The conductivity run uses a custom
/// two-bath Langevin loop.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalConfig {
    /// The integration time step, in seconds.
    pub dt_s: f64,
    /// The NVT target temperature for the specific heat, in kelvin.
    pub temperature_k: f64,
    /// The Langevin friction for the NVT run, in reciprocal seconds.
    pub friction_per_s: f64,
    /// The Berendsen coupling time for the pre-equilibration, in seconds.
    pub berendsen_tau_s: f64,
    /// The number of Berendsen pre-equilibration steps.
    pub equilibration_steps: usize,
    /// The number of Langevin production steps for the specific heat.
    pub production_steps: usize,
    /// The sample stride for the kinetic energy. It must be at least one.
    pub sample_interval: usize,
    /// The hot-bath target temperature, in kelvin.
    pub hot_temperature_k: f64,
    /// The cold-bath target temperature, in kelvin.
    pub cold_temperature_k: f64,
    /// The Langevin friction of both conductivity baths, in reciprocal
    /// seconds.
    pub bath_friction_per_s: f64,
    /// The number of steps before the conductivity measurement window.
    pub conductivity_equilibration_steps: usize,
    /// The number of conductivity production steps.
    pub conductivity_production_steps: usize,
    /// The sample stride for the heat current. It must be at least one.
    pub conductivity_sample_interval: usize,
    /// The number of atoms in each bath slab. The two slabs sit at the two
    /// ends of the atom-index range.
    pub slab_atoms: usize,
    /// The stated cross-sectional area for the conductivity. It is a model
    /// choice, not a measured value.
    pub cross_section_area_m2: f64,
    /// The seed for the thermostat noise and the initial velocities.
    pub seed: u64,
}

impl Default for ThermalConfig {
    fn default() -> Self {
        Self {
            dt_s: 5.0e-15,
            temperature_k: 300.0,
            friction_per_s: 5.0e12,
            berendsen_tau_s: 2.0e-13,
            equilibration_steps: 4_000,
            production_steps: 80_000,
            sample_interval: 20,
            hot_temperature_k: 400.0,
            cold_temperature_k: 200.0,
            bath_friction_per_s: 1.0e13,
            conductivity_equilibration_steps: 10_000,
            conductivity_production_steps: 40_000,
            conductivity_sample_interval: 1,
            slab_atoms: 2,
            cross_section_area_m2: 1.0e-20,
            seed: 0x7E57_1CE5,
        }
    }
}

/// The specific heat that an NVT run produced.
#[derive(Clone, Debug, PartialEq)]
pub struct SpecificHeatResult {
    /// The specific heat `C_v = Var(E_k) / (k_B T^2) / mass`, with the block
    /// standard error as uncertainty.
    pub specific_heat_j_kg_k: Quantity,
    /// The measured mean kinetic temperature, in kelvin.
    pub mean_temperature_k: f64,
    /// The measured mean kinetic energy, in joules.
    pub mean_kinetic_energy_j: f64,
    /// The number of kinetic-energy samples.
    pub samples: usize,
    /// The total mass, in kilograms.
    pub mass_kg: f64,
    /// The atom count.
    pub atoms: usize,
    /// The provenance of the extraction.
    pub provenance: Provenance,
}

/// The thermal conductivity that a two-bath run produced.
#[derive(Clone, Debug, PartialEq)]
pub struct ConductivityResult {
    /// The thermal conductivity, with the block standard error of the heat
    /// current propagated into it.
    pub thermal_conductivity_w_m_k: Quantity,
    /// The measured temperature gradient along the sample, in kelvin per metre.
    pub temperature_gradient_k_per_m: f64,
    /// The measured heat current, in watts.
    pub heat_current_w: f64,
    /// The measured hot-slab temperature, in kelvin.
    pub hot_temperature_k: f64,
    /// The measured cold-slab temperature, in kelvin.
    pub cold_temperature_k: f64,
    /// The distance between the two slab centres, in metres.
    pub slab_distance_m: f64,
    /// The number of heat-current samples.
    pub samples: usize,
    /// The provenance of the extraction.
    pub provenance: Provenance,
}

/// The two thermal properties of one part, from one NVT run.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalResult {
    /// The specific heat.
    pub specific_heat_j_kg_k: Quantity,
    /// The thermal conductivity.
    pub thermal_conductivity_w_m_k: Quantity,
    /// The measured mean kinetic temperature of the NVT run, in kelvin.
    pub mean_temperature_k: f64,
    /// The measured temperature gradient of the two-bath run, in kelvin per
    /// metre.
    pub temperature_gradient_k_per_m: f64,
    /// The measured heat current of the two-bath run, in watts.
    pub heat_current_w: f64,
    /// The number of kinetic-energy samples.
    pub kinetic_energy_samples: usize,
    /// The number of heat-current samples.
    pub heat_current_samples: usize,
    /// The provenance of the extraction.
    pub provenance: Provenance,
}

/// Extracts the specific heat from the kinetic-energy fluctuation under NVT.
///
/// The function initializes the velocities at the target temperature from the
/// input geometry, runs a Berendsen pre-equilibration, then samples the kinetic
/// energy during a Langevin NVT production. The estimator is
/// `Var(E_k) / (k_B T^2)`, divided by the total mass. The uncertainty is the
/// block standard error over [`THERMAL_BLOCK_COUNT`] contiguous blocks.
///
/// The run is deterministic for a fixed seed. The result is an estimate from
/// classical equipartition. It holds the kinetic contribution only.
pub fn extract_specific_heat(
    system: &mut System,
    positions_m: &[f64],
    config: &ThermalConfig,
) -> Result<SpecificHeatResult, ParamError> {
    validate_specific_heat(system, positions_m, config)?;
    let atom_count = system.atom_count();
    let total_mass_kg: f64 = system.masses_kg().iter().sum();

    let mut positions_m = positions_m.to_vec();
    let mut velocities_m_per_s = thermal_velocities_m_per_s(system, config);

    let mut berendsen =
        BerendsenThermostat::new(config.dt_s, config.temperature_k, config.berendsen_tau_s)?;
    for _ in 0..config.equilibration_steps {
        berendsen.step(&mut *system, &mut positions_m, &mut velocities_m_per_s)?;
    }

    let mut langevin = LangevinThermostat::new(
        config.dt_s,
        config.temperature_k,
        config.friction_per_s,
        config.seed,
    )?;
    let mut kinetic_energies_j =
        Vec::with_capacity(config.production_steps / config.sample_interval + 1);
    for step in 0..config.production_steps {
        langevin.step(&mut *system, &mut positions_m, &mut velocities_m_per_s)?;
        if step.is_multiple_of(config.sample_interval) {
            kinetic_energies_j.push(system.kinetic_energy_j(&velocities_m_per_s)?);
        }
    }

    summarize_specific_heat(&kinetic_energies_j, atom_count, total_mass_kg)
}

/// Runs a two-bath non-equilibrium run and extracts a thermal conductivity.
///
/// Two Langevin baths hold the first `slab_atoms` and the last `slab_atoms`
/// atoms at `hot_temperature_k` and `cold_temperature_k`. The middle atoms are
/// ballistic. In the steady state the power the hot bath adds equals the heat
/// current. The conductivity is `q * distance / (area * delta_T)`, with the
/// measured slab temperatures and the measured slab-centre distance.
///
/// The result is an estimate from a simulation. For a one-dimensional sample
/// the value is an effective coefficient for the stated geometry and the
/// stated bath coupling, not a bulk material constant. A harmonic chain shows
/// anomalous (length-dependent) transport, so the value does not converge to a
/// thermodynamic limit.
pub fn extract_thermal_conductivity(
    system: &mut System,
    positions_m: &[f64],
    config: &ThermalConfig,
) -> Result<ConductivityResult, ParamError> {
    validate_conductivity(system, positions_m, config)?;
    let atom_count = system.atom_count();
    let hot_count = config.slab_atoms;
    let cold_start = atom_count - config.slab_atoms;

    let hot_center_m = slab_center_m(positions_m, 0, hot_count);
    let cold_center_m = slab_center_m(positions_m, cold_start, atom_count);
    let slab_distance_m = (cold_center_m - hot_center_m).abs();

    let mut positions = positions_m.to_vec();
    let mut velocities_m_per_s = thermal_velocities_m_per_s(system, config);
    let masses_kg = system.masses_kg().to_vec();
    let integrator = VelocityVerlet::new(config.dt_s)?;
    let decay = (-config.bath_friction_per_s * config.dt_s).exp();
    let noise_scale = (1.0 - decay * decay).sqrt();
    let bath = BathStep { decay, noise_scale };
    let mut rng = Rng::new(config.seed ^ 0x9E37_79B9_7F4A_7C15);

    let total_steps = config
        .conductivity_equilibration_steps
        .saturating_add(config.conductivity_production_steps);
    let mut heat_currents_w = Vec::with_capacity(
        config.conductivity_production_steps / config.conductivity_sample_interval + 1,
    );
    let mut hot_temperature_sum_k = 0.0;
    let mut cold_temperature_sum_k = 0.0;
    let mut temperature_samples = 0usize;

    for step in 0..total_steps {
        integrator.step(&mut *system, &mut positions, &mut velocities_m_per_s)?;
        let hot_delta_j = bath.apply(
            &mut velocities_m_per_s,
            &masses_kg,
            0,
            hot_count,
            config.hot_temperature_k,
            &mut rng,
        );
        let cold_delta_j = bath.apply(
            &mut velocities_m_per_s,
            &masses_kg,
            cold_start,
            atom_count,
            config.cold_temperature_k,
            &mut rng,
        );
        if step < config.conductivity_equilibration_steps {
            continue;
        }
        let in_window = step - config.conductivity_equilibration_steps;
        if !in_window.is_multiple_of(config.conductivity_sample_interval) {
            continue;
        }
        let heat_current_w = (hot_delta_j - cold_delta_j) / (2.0 * config.dt_s);
        heat_currents_w.push(heat_current_w);

        let hot_kinetic_j = slab_kinetic_energy_j(&velocities_m_per_s, &masses_kg, 0, hot_count);
        let cold_kinetic_j =
            slab_kinetic_energy_j(&velocities_m_per_s, &masses_kg, cold_start, atom_count);
        hot_temperature_sum_k += slab_temperature_k(hot_kinetic_j, hot_count);
        cold_temperature_sum_k += slab_temperature_k(cold_kinetic_j, config.slab_atoms);
        temperature_samples += 1;
    }

    summarize_conductivity(
        &heat_currents_w,
        hot_temperature_sum_k,
        cold_temperature_sum_k,
        temperature_samples,
        slab_distance_m,
        config,
    )
}

/// Runs both thermal extractions and combines them into one result.
///
/// See [`extract_specific_heat`] and [`extract_thermal_conductivity`] for the
/// two models and their limits. The combined provenance names both models.
pub fn extract_thermal(
    system: &mut System,
    positions_m: &[f64],
    config: &ThermalConfig,
) -> Result<ThermalResult, ParamError> {
    let specific = extract_specific_heat(&mut *system, positions_m, config)?;
    let conductivity = extract_thermal_conductivity(&mut *system, positions_m, config)?;
    let notes = format!(
        "NVT kinetic-energy fluctuation for the specific heat over {} samples at {:.3} K; \
         two-bath non-equilibrium estimate for the conductivity over {} samples with area {:e} m^2; \
         both are estimates, not material constants",
        specific.samples, specific.mean_temperature_k, conductivity.samples, config.cross_section_area_m2
    );
    Ok(ThermalResult {
        specific_heat_j_kg_k: specific.specific_heat_j_kg_k,
        thermal_conductivity_w_m_k: conductivity.thermal_conductivity_w_m_k,
        mean_temperature_k: specific.mean_temperature_k,
        temperature_gradient_k_per_m: conductivity.temperature_gradient_k_per_m,
        heat_current_w: conductivity.heat_current_w,
        kinetic_energy_samples: specific.samples,
        heat_current_samples: conductivity.samples,
        provenance: extraction_provenance(Method::Md, notes),
    })
}

/// Builds a result from kinetic-energy samples. It isolates the estimator from
/// the integration run for focused tests.
fn summarize_specific_heat(
    kinetic_energies_j: &[f64],
    atom_count: usize,
    total_mass_kg: f64,
) -> Result<SpecificHeatResult, ParamError> {
    if kinetic_energies_j.len() < MIN_THERMAL_SAMPLES {
        return Err(ParamError::ThermalSamplesTooFew {
            required: MIN_THERMAL_SAMPLES,
            measured: kinetic_energies_j.len(),
        });
    }
    if !kinetic_energies_j
        .iter()
        .all(|energy_j| energy_j.is_finite())
    {
        return Err(ParamError::NonFiniteThermalObservable {
            observable: "kinetic energy",
        });
    }
    if atom_count == 0 || !total_mass_kg.is_finite() || total_mass_kg <= 0.0 {
        return Err(ParamError::InvalidExtraction {
            reason: "the total mass must be finite and positive".to_string(),
        });
    }

    let mean_kinetic_energy_j =
        mean(kinetic_energies_j).ok_or(ParamError::ThermalSamplesTooFew {
            required: MIN_THERMAL_SAMPLES,
            measured: 0,
        })?;
    let mean_temperature_k =
        2.0 * mean_kinetic_energy_j / (3.0 * atom_count as f64 * BOLTZMANN_J_PER_K);
    if !mean_temperature_k.is_finite() || mean_temperature_k <= 0.0 {
        return Err(ParamError::ThermalGradientNotEstablished {
            hot_k: mean_temperature_k,
            cold_k: 0.0,
        });
    }

    let variance = sample_variance(kinetic_energies_j).ok_or(ParamError::ThermalSamplesTooFew {
        required: MIN_THERMAL_SAMPLES,
        measured: kinetic_energies_j.len(),
    })?;
    let specific_heat_si =
        variance / (BOLTZMANN_J_PER_K * mean_temperature_k * mean_temperature_k * total_mass_kg);
    let specific_heat_uncertainty_si =
        block_specific_heat_standard_error(kinetic_energies_j, atom_count, total_mass_kg)
            .unwrap_or(0.0);

    let notes = format!(
        "kinetic-energy fluctuation over {} samples at {:.3} K for {} atoms; \
         C_v = Var(E_k) / (k_B T^2) / mass; classical equipartition estimate, \
         kinetic contribution only",
        kinetic_energies_j.len(),
        mean_temperature_k,
        atom_count
    );
    Ok(SpecificHeatResult {
        specific_heat_j_kg_k: Quantity::derived(specific_heat_si, "J/(kg*K)")?
            .with_uncertainty_si(specific_heat_uncertainty_si),
        mean_temperature_k,
        mean_kinetic_energy_j,
        samples: kinetic_energies_j.len(),
        mass_kg: total_mass_kg,
        atoms: atom_count,
        provenance: extraction_provenance(Method::Md, notes),
    })
}

fn summarize_conductivity(
    heat_currents_w: &[f64],
    hot_temperature_sum_k: f64,
    cold_temperature_sum_k: f64,
    temperature_samples: usize,
    slab_distance_m: f64,
    config: &ThermalConfig,
) -> Result<ConductivityResult, ParamError> {
    if heat_currents_w.len() < MIN_THERMAL_SAMPLES {
        return Err(ParamError::ThermalSamplesTooFew {
            required: MIN_THERMAL_SAMPLES,
            measured: heat_currents_w.len(),
        });
    }
    if !heat_currents_w
        .iter()
        .all(|current_w| current_w.is_finite())
    {
        return Err(ParamError::NonFiniteThermalObservable {
            observable: "heat current",
        });
    }
    if temperature_samples == 0 {
        return Err(ParamError::ThermalSamplesTooFew {
            required: 1,
            measured: 0,
        });
    }
    let hot_temperature_k = hot_temperature_sum_k / temperature_samples as f64;
    let cold_temperature_k = cold_temperature_sum_k / temperature_samples as f64;
    let delta_temperature_k = hot_temperature_k - cold_temperature_k;
    if !delta_temperature_k.is_finite() || delta_temperature_k <= 0.0 {
        return Err(ParamError::ThermalGradientNotEstablished {
            hot_k: hot_temperature_k,
            cold_k: cold_temperature_k,
        });
    }
    let area_m2 = config.cross_section_area_m2;
    if !area_m2.is_finite()
        || area_m2 <= 0.0
        || !slab_distance_m.is_finite()
        || slab_distance_m <= 0.0
    {
        return Err(ParamError::InvalidExtraction {
            reason: "the conductivity needs a positive area and slab distance".to_string(),
        });
    }

    let heat_current_w = mean(heat_currents_w).ok_or(ParamError::ThermalSamplesTooFew {
        required: MIN_THERMAL_SAMPLES,
        measured: heat_currents_w.len(),
    })?;
    let temperature_gradient_k_per_m = delta_temperature_k / slab_distance_m;
    let conductivity_si = heat_current_w / (area_m2 * temperature_gradient_k_per_m);
    let heat_current_uncertainty_w =
        block_standard_error(heat_currents_w, THERMAL_BLOCK_COUNT).unwrap_or(0.0);
    let conductivity_uncertainty_si =
        heat_current_uncertainty_w / (area_m2 * temperature_gradient_k_per_m);

    let notes = format!(
        "two-bath non-equilibrium run over {} samples; hot {:.3} K, cold {:.3} K, \
         gradient {:e} K/m, heat current {:e} W, area {:e} m^2; effective estimate \
         for the stated geometry, not a bulk constant",
        heat_currents_w.len(),
        hot_temperature_k,
        cold_temperature_k,
        temperature_gradient_k_per_m,
        heat_current_w,
        area_m2
    );
    Ok(ConductivityResult {
        thermal_conductivity_w_m_k: Quantity::derived(conductivity_si, "W/(m*K)")?
            .with_uncertainty_si(conductivity_uncertainty_si),
        temperature_gradient_k_per_m,
        heat_current_w,
        hot_temperature_k,
        cold_temperature_k,
        slab_distance_m,
        samples: heat_currents_w.len(),
        provenance: extraction_provenance(Method::Md, notes),
    })
}

/// The two precomputed coefficients of one Ornstein-Uhlenbeck Langevin step.
#[derive(Clone, Copy, Debug)]
struct BathStep {
    /// The velocity decay `exp(-gamma dt)`.
    decay: f64,
    /// The noise scale `sqrt(1 - decay^2)`.
    noise_scale: f64,
}

impl BathStep {
    /// Applies one update to a slab of atoms.
    ///
    /// Returns the change in the slab kinetic energy, in joules. The update
    /// matches [`LangevinThermostat`]:
    /// `v <- decay v + sqrt(1 - decay^2) sqrt(k_B T / m) xi`.
    fn apply(
        self,
        velocities_m_per_s: &mut [f64],
        masses_kg: &[f64],
        start: usize,
        end: usize,
        temperature_k: f64,
        rng: &mut Rng,
    ) -> f64 {
        let before_j = slab_kinetic_energy_j(velocities_m_per_s, masses_kg, start, end);
        for (atom, &mass_kg) in masses_kg.iter().enumerate().take(end).skip(start) {
            let sigma_m_per_s =
                self.noise_scale * (BOLTZMANN_J_PER_K * temperature_k / mass_kg).sqrt();
            let base = atom * 3;
            for axis in 0..3 {
                let index = base + axis;
                velocities_m_per_s[index] =
                    self.decay * velocities_m_per_s[index] + sigma_m_per_s * rng.next_normal();
            }
        }
        let after_j = slab_kinetic_energy_j(velocities_m_per_s, masses_kg, start, end);
        after_j - before_j
    }
}

fn thermal_velocities_m_per_s(system: &System, config: &ThermalConfig) -> Vec<f64> {
    let mut rng = Rng::new(config.seed ^ 0x51AB_5EED_0BAD_F00D);
    let mut velocities_m_per_s = Vec::with_capacity(3 * system.atom_count());
    for &mass_kg in system.masses_kg() {
        let sigma_m_per_s = (BOLTZMANN_J_PER_K * config.temperature_k / mass_kg).sqrt();
        for _ in 0..3 {
            velocities_m_per_s.push(sigma_m_per_s * rng.next_normal());
        }
    }
    velocities_m_per_s
}

fn slab_center_m(positions_m: &[f64], start: usize, end: usize) -> f64 {
    let mut sum_m = 0.0;
    for atom in start..end {
        sum_m += positions_m[3 * atom];
    }
    sum_m / (end - start) as f64
}

fn slab_kinetic_energy_j(
    velocities_m_per_s: &[f64],
    masses_kg: &[f64],
    start: usize,
    end: usize,
) -> f64 {
    let mut kinetic_j = 0.0;
    for (atom, &mass_kg) in masses_kg.iter().enumerate().take(end).skip(start) {
        let base = atom * 3;
        let vx = velocities_m_per_s[base];
        let vy = velocities_m_per_s[base + 1];
        let vz = velocities_m_per_s[base + 2];
        kinetic_j += 0.5 * mass_kg * (vx * vx + vy * vy + vz * vz);
    }
    kinetic_j
}

fn slab_temperature_k(slab_kinetic_j: f64, slab_atoms: usize) -> f64 {
    2.0 * slab_kinetic_j / (3.0 * slab_atoms as f64 * BOLTZMANN_J_PER_K)
}

/// The population variance with the `n - 1` denominator, or [`None`] below two
/// values.
fn sample_variance(values: &[f64]) -> Option<f64> {
    sample_std_dev(values).map(|spread| spread * spread)
}

/// The standard error of the mean from `block_count` contiguous blocks.
fn block_standard_error(values: &[f64], block_count: usize) -> Option<f64> {
    if block_count < 2 || values.len() < 2 * block_count {
        return None;
    }
    let block_len = values.len() / block_count;
    let mut block_means = Vec::with_capacity(block_count);
    for block in 0..block_count {
        let start = block * block_len;
        let end = if block == block_count - 1 {
            values.len()
        } else {
            start + block_len
        };
        block_means.push(mean(&values[start..end])?);
    }
    let spread = sample_std_dev(&block_means)?;
    Some(spread / (block_count as f64).sqrt())
}

/// The block standard error of the specific heat itself.
fn block_specific_heat_standard_error(
    kinetic_energies_j: &[f64],
    atom_count: usize,
    total_mass_kg: f64,
) -> Option<f64> {
    let block_count = THERMAL_BLOCK_COUNT;
    if kinetic_energies_j.len() < 2 * block_count {
        return None;
    }
    let block_len = kinetic_energies_j.len() / block_count;
    let mut block_specific_heats = Vec::with_capacity(block_count);
    for block in 0..block_count {
        let start = block * block_len;
        let end = if block == block_count - 1 {
            kinetic_energies_j.len()
        } else {
            start + block_len
        };
        let values = &kinetic_energies_j[start..end];
        let mean_j = mean(values)?;
        let mean_temperature_k = 2.0 * mean_j / (3.0 * atom_count as f64 * BOLTZMANN_J_PER_K);
        if mean_temperature_k <= 0.0 {
            return None;
        }
        let variance = sample_variance(values)?;
        block_specific_heats.push(
            variance
                / (BOLTZMANN_J_PER_K * mean_temperature_k * mean_temperature_k * total_mass_kg),
        );
    }
    let spread = sample_std_dev(&block_specific_heats)?;
    Some(spread / (block_count as f64).sqrt())
}

fn validate_specific_heat(
    system: &System,
    positions_m: &[f64],
    config: &ThermalConfig,
) -> Result<(), ParamError> {
    if positions_m.len() != 3 * system.atom_count() {
        return Err(EngineError::BufferSizeMismatch {
            len: positions_m.len(),
            expected: 3 * system.atom_count(),
        }
        .into());
    }
    if !system.has_masses() {
        return Err(EngineError::MassesNotSet.into());
    }
    if system.atom_count() == 0 {
        return Err(ParamError::InvalidExtraction {
            reason: "a thermal extraction needs at least one atom".to_string(),
        });
    }
    if !config.dt_s.is_finite() || config.dt_s <= 0.0 {
        return Err(ParamError::InvalidExtraction {
            reason: "the time step must be finite and positive".to_string(),
        });
    }
    if !config.temperature_k.is_finite() || config.temperature_k <= 0.0 {
        return Err(ParamError::InvalidExtraction {
            reason: "the temperature must be finite and positive".to_string(),
        });
    }
    if !config.friction_per_s.is_finite() || config.friction_per_s <= 0.0 {
        return Err(ParamError::InvalidExtraction {
            reason: "the Langevin friction must be finite and positive".to_string(),
        });
    }
    if !config.berendsen_tau_s.is_finite() || config.berendsen_tau_s <= 0.0 {
        return Err(ParamError::InvalidExtraction {
            reason: "the Berendsen coupling time must be finite and positive".to_string(),
        });
    }
    if config.sample_interval == 0 {
        return Err(ParamError::InvalidExtraction {
            reason: "the sample interval must be at least one".to_string(),
        });
    }
    let sample_count = config.production_steps / config.sample_interval + 1;
    if sample_count < MIN_THERMAL_SAMPLES {
        return Err(ParamError::ThermalSamplesTooFew {
            required: MIN_THERMAL_SAMPLES,
            measured: sample_count,
        });
    }
    Ok(())
}

fn validate_conductivity(
    system: &System,
    positions_m: &[f64],
    config: &ThermalConfig,
) -> Result<(), ParamError> {
    if positions_m.len() != 3 * system.atom_count() {
        return Err(EngineError::BufferSizeMismatch {
            len: positions_m.len(),
            expected: 3 * system.atom_count(),
        }
        .into());
    }
    if !system.has_masses() {
        return Err(EngineError::MassesNotSet.into());
    }
    let atom_count = system.atom_count();
    if config.slab_atoms == 0 || 2 * config.slab_atoms >= atom_count {
        return Err(ParamError::ThermalSlabsTooLarge {
            requested: 2 * config.slab_atoms + 1,
            available: atom_count,
        });
    }
    if !config.dt_s.is_finite() || config.dt_s <= 0.0 {
        return Err(ParamError::InvalidExtraction {
            reason: "the time step must be finite and positive".to_string(),
        });
    }
    if !config.hot_temperature_k.is_finite()
        || !config.cold_temperature_k.is_finite()
        || config.hot_temperature_k <= config.cold_temperature_k
        || config.cold_temperature_k <= 0.0
    {
        return Err(ParamError::ThermalGradientNotEstablished {
            hot_k: config.hot_temperature_k,
            cold_k: config.cold_temperature_k,
        });
    }
    if !config.bath_friction_per_s.is_finite() || config.bath_friction_per_s <= 0.0 {
        return Err(ParamError::InvalidExtraction {
            reason: "the bath friction must be finite and positive".to_string(),
        });
    }
    if !config.cross_section_area_m2.is_finite() || config.cross_section_area_m2 <= 0.0 {
        return Err(ParamError::InvalidExtraction {
            reason: "the cross-sectional area must be finite and positive".to_string(),
        });
    }
    if config.conductivity_sample_interval == 0 {
        return Err(ParamError::InvalidExtraction {
            reason: "the conductivity sample interval must be at least one".to_string(),
        });
    }
    let sample_count =
        config.conductivity_production_steps / config.conductivity_sample_interval + 1;
    if sample_count < MIN_THERMAL_SAMPLES {
        return Err(ParamError::ThermalSamplesTooFew {
            required: MIN_THERMAL_SAMPLES,
            measured: sample_count,
        });
    }
    let hot_count = config.slab_atoms;
    let cold_start = atom_count - config.slab_atoms;
    let distance_m = (slab_center_m(positions_m, cold_start, atom_count)
        - slab_center_m(positions_m, 0, hot_count))
    .abs();
    if !distance_m.is_finite() || distance_m <= 0.0 {
        return Err(ParamError::ThermalGradientNotEstablished {
            hot_k: config.hot_temperature_k,
            cold_k: config.cold_temperature_k,
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::{uniform_chain, CARBON_MASS_KG};

    const IDEAL_ATOMS: usize = 40;
    const IDEAL_MASS_KG: f64 = CARBON_MASS_KG;

    fn ideal_gas() -> System {
        System::with_masses_kg(IDEAL_ATOMS, &[IDEAL_MASS_KG; IDEAL_ATOMS]).expect("valid masses")
    }

    fn ideal_positions_m() -> Vec<f64> {
        vec![0.0; 3 * IDEAL_ATOMS]
    }

    fn ideal_gas_config() -> ThermalConfig {
        ThermalConfig {
            production_steps: 120_000,
            sample_interval: 25,
            ..ThermalConfig::default()
        }
    }

    #[test]
    fn a_monatomic_ideal_gas_recovers_the_equipartition_specific_heat() {
        let mut system = ideal_gas();
        let result = extract_specific_heat(&mut system, &ideal_positions_m(), &ideal_gas_config())
            .expect("valid run");
        let expected_j_kg_k = 1.5 * BOLTZMANN_J_PER_K / IDEAL_MASS_KG;
        let relative_error =
            (result.specific_heat_j_kg_k.value_si() - expected_j_kg_k).abs() / expected_j_kg_k;
        println!(
            "specific heat: extracted {:e} J/(kg*K), expected {:e} J/(kg*K), \
             uncertainty {:e}, relative error {:e}, T {:.3} K, samples {}",
            result.specific_heat_j_kg_k.value_si(),
            expected_j_kg_k,
            result.specific_heat_j_kg_k.uncertainty_si().unwrap_or(0.0),
            relative_error,
            result.mean_temperature_k,
            result.samples
        );
        assert!(relative_error < 0.05, "relative error {relative_error:e}");
        assert!(result.specific_heat_j_kg_k.uncertainty_si().unwrap_or(0.0) > 0.0);
        assert_eq!(result.provenance.method, Method::Md);
    }

    #[test]
    fn the_specific_heat_run_repeats_for_a_fixed_seed() {
        let config = ideal_gas_config();
        let positions_m = ideal_positions_m();
        let mut first_system = ideal_gas();
        let first =
            extract_specific_heat(&mut first_system, &positions_m, &config).expect("valid run");
        let mut second_system = ideal_gas();
        let second =
            extract_specific_heat(&mut second_system, &positions_m, &config).expect("valid run");
        assert_eq!(first.specific_heat_j_kg_k, second.specific_heat_j_kg_k);
        assert_eq!(first.mean_kinetic_energy_j, second.mean_kinetic_energy_j);
        assert_eq!(first.samples, second.samples);
    }

    #[test]
    fn a_chain_gives_a_two_bath_conductivity_that_repeats_and_is_positive() {
        let config = ThermalConfig::default();
        let (mut system, positions_m) = uniform_chain();
        let first =
            extract_thermal_conductivity(&mut system, &positions_m, &config).expect("valid run");
        let (mut system, positions_m) = uniform_chain();
        let second =
            extract_thermal_conductivity(&mut system, &positions_m, &config).expect("valid run");
        println!(
            "conductivity: {:e} W/(m*K), uncertainty {:e}, heat current {:e} W, \
             gradient {:e} K/m, distance {:e} m, samples {}",
            first.thermal_conductivity_w_m_k.value_si(),
            first
                .thermal_conductivity_w_m_k
                .uncertainty_si()
                .unwrap_or(0.0),
            first.heat_current_w,
            first.temperature_gradient_k_per_m,
            first.slab_distance_m,
            first.samples
        );
        assert_eq!(
            first.thermal_conductivity_w_m_k,
            second.thermal_conductivity_w_m_k
        );
        assert_eq!(first.heat_current_w, second.heat_current_w);
        assert_eq!(
            first.temperature_gradient_k_per_m,
            second.temperature_gradient_k_per_m
        );
        assert!(first.thermal_conductivity_w_m_k.value_si().is_finite());
        assert!(first.thermal_conductivity_w_m_k.value_si() >= 0.0);
        assert!(first.temperature_gradient_k_per_m > 0.0);
        assert!(first.heat_current_w > 0.0);
    }

    #[test]
    fn the_conductivity_stays_within_a_factor_across_seeds() {
        let (mut system, positions_m) = uniform_chain();
        let first =
            extract_thermal_conductivity(&mut system, &positions_m, &ThermalConfig::default())
                .expect("valid run");
        let (mut system, positions_m) = uniform_chain();
        let second = extract_thermal_conductivity(
            &mut system,
            &positions_m,
            &ThermalConfig {
                seed: 0x0BAD_F00D_1234_5678,
                ..ThermalConfig::default()
            },
        )
        .expect("valid run");
        let ratio = first.thermal_conductivity_w_m_k.value_si()
            / second.thermal_conductivity_w_m_k.value_si();
        println!(
            "conductivity cross-seed: {:e} and {:e} W/(m*K), ratio {ratio:.4}",
            first.thermal_conductivity_w_m_k.value_si(),
            second.thermal_conductivity_w_m_k.value_si()
        );
        assert!(ratio > 0.0 && ratio.is_finite());
        assert!(
            (0.25..=4.0).contains(&ratio),
            "the two seeded runs differ by more than a factor of four: {ratio}"
        );
    }

    #[test]
    fn both_properties_fill_the_thermal_fields_of_a_record() {
        let config = ThermalConfig {
            production_steps: 40_000,
            conductivity_production_steps: 20_000,
            ..ThermalConfig::default()
        };
        let (mut system, positions_m) = uniform_chain();
        let result = extract_thermal(&mut system, &positions_m, &config).expect("valid run");
        println!(
            "thermal record: C_v {:e} J/(kg*K) +/- {:e}, k {:e} W/(m*K) +/- {:e}",
            result.specific_heat_j_kg_k.value_si(),
            result.specific_heat_j_kg_k.uncertainty_si().unwrap_or(0.0),
            result.thermal_conductivity_w_m_k.value_si(),
            result
                .thermal_conductivity_w_m_k
                .uncertainty_si()
                .unwrap_or(0.0)
        );
        assert_eq!(result.specific_heat_j_kg_k.unit(), "J/(kg*K)");
        assert_eq!(result.thermal_conductivity_w_m_k.unit(), "W/(m*K)");
        assert!(result.specific_heat_j_kg_k.uncertainty_si().is_some());
        assert!(result.thermal_conductivity_w_m_k.uncertainty_si().is_some());
        assert_eq!(result.provenance.method, Method::Md);
        assert_eq!(
            result.provenance.validation,
            crate::provenance::Validation::Unverified
        );
    }

    #[test]
    fn an_impossible_thermal_setup_is_rejected() {
        let mut system = ideal_gas();
        let zero_friction = ThermalConfig {
            friction_per_s: 0.0,
            ..ideal_gas_config()
        };
        assert!(matches!(
            extract_specific_heat(&mut system, &ideal_positions_m(), &zero_friction),
            Err(ParamError::InvalidExtraction { .. })
        ));

        let (mut system, positions_m) = uniform_chain();
        let tiny_bath = ThermalConfig {
            slab_atoms: 20,
            ..ThermalConfig::default()
        };
        assert!(matches!(
            extract_thermal_conductivity(&mut system, &positions_m, &tiny_bath),
            Err(ParamError::ThermalSlabsTooLarge { .. })
        ));

        let no_gradient = ThermalConfig {
            hot_temperature_k: 200.0,
            cold_temperature_k: 300.0,
            ..ThermalConfig::default()
        };
        assert!(matches!(
            extract_thermal_conductivity(&mut system, &positions_m, &no_gradient),
            Err(ParamError::ThermalGradientNotEstablished { .. })
        ));
    }
}
