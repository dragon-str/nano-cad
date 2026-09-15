//! Stiffness extraction from a controlled strain sweep.
//!
//! This implements step 2 of the extraction pipeline in `PARAMETERS.md`:
//! apply small strains, record the stress response, and fit the modulus. The
//! fit reports its uncertainty, which is the standard error of the slope.
//!
//! The strain axis is the x axis, applied about the sample centroid. The lateral
//! dimensions do not relax, so the result is a longitudinal stiffness at the
//! stated cross-sectional area. The area and the reference length are inputs.

use nanocad_engine::{EngineError, System};

use crate::error::ParamError;
use crate::provenance::{extraction_provenance, Method, Provenance};
use crate::quantity::Quantity;
use crate::stats::{fit_line, Rng};
use crate::strain::stress_from_strain;

/// The strain sweep and sample geometry for a stiffness extraction.
///
/// The sample volume is `reference_length_m * cross_section_area_m2`. The area
/// is a stated model choice. Change it and the modulus changes with it, because
/// a smaller area carries the same force at a higher stress.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StiffnessConfig {
    /// The lowest applied strain. It may be zero or negative.
    pub strain_min: f64,
    /// The highest applied strain. It must exceed `strain_min`.
    pub strain_max: f64,
    /// The number of strain samples. It must be at least three.
    pub samples: usize,
    /// The reference sample length along the strain axis, in metres.
    pub reference_length_m: f64,
    /// The cross-sectional area normal to the strain axis, in square metres.
    pub cross_section_area_m2: f64,
    /// The central-difference step on the strain, dimensionless.
    pub derivative_step: f64,
    /// The half-width of deterministic synthetic stress noise, in pascals.
    ///
    /// It is zero for a noiseless run. A positive value models a noisy
    /// measurement and makes the fit uncertainty meaningful.
    pub noise_amplitude_pa: f64,
    /// The seed for the synthetic noise.
    pub seed: u64,
}

impl Default for StiffnessConfig {
    fn default() -> Self {
        Self {
            strain_min: -1.0e-3,
            strain_max: 1.0e-3,
            samples: 11,
            reference_length_m: 3.0e-9,
            cross_section_area_m2: 1.0e-20,
            derivative_step: 1.0e-5,
            noise_amplitude_pa: 0.0,
            seed: 0x5EED_1234,
        }
    }
}

/// The elastic modulus that a strain sweep produced, with its fit quality.
#[derive(Clone, Debug, PartialEq)]
pub struct StiffnessResult {
    /// The fitted elastic modulus, with the slope standard error attached.
    pub elastic_modulus_pa: Quantity,
    /// The fitted stress intercept in pascals. A relaxed sample gives near zero.
    pub fit_intercept_pa: f64,
    /// The lowest and highest strain that were sampled.
    pub strain_range: (f64, f64),
    /// The number of strain samples that were fitted.
    pub samples: usize,
    /// The coefficient of determination of the linear fit.
    pub r_squared: f64,
    /// The provenance of the extraction.
    pub provenance: Provenance,
}

/// Extracts an elastic modulus from a strain sweep on an engine system.
///
/// The function applies a uniaxial x strain to `positions_m`, evaluates the
/// engine energy at each strain, and takes the stress as the central difference
/// of the energy with respect to the strain, divided by the sample volume. It
/// then fits `stress = intercept + modulus * strain` and reports the modulus
/// with the standard error of the slope.
///
/// The result is an estimate. It is not a certified material constant. The
/// stated area fixes the scale, and the model holds the lateral dimensions
/// fixed, so the value is a constrained longitudinal modulus for this geometry.
pub fn extract_stiffness(
    system: &mut System,
    positions_m: &[f64],
    config: &StiffnessConfig,
) -> Result<StiffnessResult, ParamError> {
    validate_config(system, positions_m, config)?;
    let volume_m3 = config.reference_length_m * config.cross_section_area_m2;
    let mut rng = Rng::new(config.seed);
    let mut strains = Vec::with_capacity(config.samples);
    let mut stresses_pa = Vec::with_capacity(config.samples);
    for index in 0..config.samples {
        let fraction = index as f64 / (config.samples - 1) as f64;
        let strain = config.strain_min + fraction * (config.strain_max - config.strain_min);
        let stress_pa = stress_from_strain(
            system,
            positions_m,
            strain,
            config.derivative_step,
            volume_m3,
        )?;
        let noise_pa = if config.noise_amplitude_pa > 0.0 {
            rng.uniform_symmetric(config.noise_amplitude_pa)
        } else {
            0.0
        };
        strains.push(strain);
        stresses_pa.push(stress_pa + noise_pa);
    }
    let fit = fit_line(&strains, &stresses_pa).ok_or_else(|| ParamError::InvalidExtraction {
        reason: "the strain sweep did not give a fit".to_string(),
    })?;
    let elastic_modulus_pa =
        Quantity::derived(fit.slope, "Pa")?.with_uncertainty_si(fit.slope_standard_error);
    let notes = format!(
        "strain sweep over [{:e}, {:e}] with {} samples and area {:e} m^2; \
         modulus is the least-squares slope of stress versus strain",
        config.strain_min, config.strain_max, config.samples, config.cross_section_area_m2
    );
    Ok(StiffnessResult {
        elastic_modulus_pa,
        fit_intercept_pa: fit.intercept,
        strain_range: (config.strain_min, config.strain_max),
        samples: config.samples,
        r_squared: fit.r_squared,
        provenance: extraction_provenance(Method::Fit, notes),
    })
}

fn validate_config(
    system: &System,
    positions_m: &[f64],
    config: &StiffnessConfig,
) -> Result<(), ParamError> {
    if positions_m.len() != 3 * system.atom_count() {
        return Err(EngineError::BufferSizeMismatch {
            len: positions_m.len(),
            expected: 3 * system.atom_count(),
        }
        .into());
    }
    if !config.strain_min.is_finite() || !config.strain_max.is_finite() {
        return Err(ParamError::InvalidExtraction {
            reason: "the strain range must be finite".to_string(),
        });
    }
    if config.strain_max <= config.strain_min {
        return Err(ParamError::InvalidExtraction {
            reason: "strain_max must exceed strain_min".to_string(),
        });
    }
    if config.samples < 3 {
        return Err(ParamError::InvalidExtraction {
            reason: "at least three strain samples are needed for a fit".to_string(),
        });
    }
    if !config.reference_length_m.is_finite() || config.reference_length_m <= 0.0 {
        return Err(ParamError::InvalidExtraction {
            reason: "the reference length must be finite and positive".to_string(),
        });
    }
    if !config.cross_section_area_m2.is_finite() || config.cross_section_area_m2 <= 0.0 {
        return Err(ParamError::InvalidExtraction {
            reason: "the cross-sectional area must be finite and positive".to_string(),
        });
    }
    if !config.derivative_step.is_finite() || config.derivative_step <= 0.0 {
        return Err(ParamError::InvalidExtraction {
            reason: "the derivative step must be finite and positive".to_string(),
        });
    }
    if !config.noise_amplitude_pa.is_finite() || config.noise_amplitude_pa < 0.0 {
        return Err(ParamError::InvalidExtraction {
            reason: "the noise half-width must be finite and nonnegative".to_string(),
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::{true_modulus_pa, uniform_chain, AREA_M2, BOND_COUNT, R0_M};

    fn config() -> StiffnessConfig {
        StiffnessConfig {
            reference_length_m: BOND_COUNT as f64 * R0_M,
            cross_section_area_m2: AREA_M2,
            derivative_step: 1.0e-5,
            ..StiffnessConfig::default()
        }
    }

    #[test]
    fn a_harmonic_chain_recovers_its_known_modulus() {
        let (mut system, positions_m) = uniform_chain();
        let result = extract_stiffness(&mut system, &positions_m, &config()).expect("valid run");
        let expected_pa = true_modulus_pa();
        let relative_error =
            (result.elastic_modulus_pa.value_si() - expected_pa).abs() / expected_pa;
        println!(
            "stiffness: extracted {:e} Pa, expected {:e} Pa, relative error {:e}",
            result.elastic_modulus_pa.value_si(),
            expected_pa,
            relative_error
        );
        assert!(relative_error < 1.0e-6, "relative error {relative_error:e}");
        assert!(result.r_squared > 0.999_999);
    }

    #[test]
    fn the_fit_uncertainty_covers_the_known_modulus_under_noise() {
        let (mut system, positions_m) = uniform_chain();
        let expected_pa = true_modulus_pa();
        let noisy = StiffnessConfig {
            strain_min: -1.0e-2,
            strain_max: 1.0e-2,
            samples: 41,
            noise_amplitude_pa: 2.0e-3 * expected_pa,
            seed: 0xABCD_EF01,
            ..config()
        };
        let result = extract_stiffness(&mut system, &positions_m, &noisy).expect("valid run");
        let uncertainty_pa = result
            .elastic_modulus_pa
            .uncertainty_si()
            .expect("the fit has an uncertainty");
        let deviation_pa = (result.elastic_modulus_pa.value_si() - expected_pa).abs();
        println!(
            "stiffness with noise: extracted {:e} Pa, uncertainty {:e} Pa, deviation {:e} Pa",
            result.elastic_modulus_pa.value_si(),
            uncertainty_pa,
            deviation_pa
        );
        assert!(uncertainty_pa > 0.0);
        assert!(
            deviation_pa <= 3.0 * uncertainty_pa,
            "the known modulus lies outside three standard errors"
        );
    }

    #[test]
    fn an_invalid_sweep_is_rejected() {
        let (mut system, positions_m) = uniform_chain();
        let bad = StiffnessConfig {
            strain_max: -1.0,
            ..config()
        };
        assert!(matches!(
            extract_stiffness(&mut system, &positions_m, &bad),
            Err(ParamError::InvalidExtraction { .. })
        ));
        let too_few = StiffnessConfig {
            samples: 2,
            ..config()
        };
        assert!(matches!(
            extract_stiffness(&mut system, &positions_m, &too_few),
            Err(ParamError::InvalidExtraction { .. })
        ));
    }
}
