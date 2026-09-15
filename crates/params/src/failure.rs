//! Failure-stress extraction from a load sweep.
//!
//! This implements step 4 of the extraction pipeline in `PARAMETERS.md`:
//! increase the load until a failure criterion triggers, then record the
//! stress. The criterion is a per-bond strain limit. A bond fails when its
//! local strain reaches its stated rupture strain.
//!
//! The result is an estimate from a simulation. It is not a certified material
//! constant.

use nanocad_engine::{EngineError, System};

use crate::error::ParamError;
use crate::provenance::{extraction_provenance, Method, Provenance};
use crate::quantity::Quantity;
use crate::strain::{stress_from_strain, uniaxial_strain_positions};

/// One bond limit for a failure criterion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BondLimit {
    /// The first atom index.
    pub u: u32,
    /// The second atom index.
    pub v: u32,
    /// The equilibrium bond length, in metres.
    pub r0_m: f64,
    /// The rupture strain. A local strain at or above this value fails the bond.
    pub rupture_strain: f64,
}

/// The load sweep and sample geometry for a failure extraction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FailureConfig {
    /// The strain increment on the first sweep, dimensionless.
    pub strain_step: f64,
    /// The highest applied strain. No failure below it is an error.
    pub max_strain: f64,
    /// The reference sample length along the strain axis, in metres.
    pub reference_length_m: f64,
    /// The cross-sectional area normal to the strain axis, in square metres.
    pub cross_section_area_m2: f64,
    /// The central-difference step on the strain for the stress.
    pub derivative_step: f64,
}

impl Default for FailureConfig {
    fn default() -> Self {
        Self {
            strain_step: 1.0e-3,
            max_strain: 0.2,
            reference_length_m: 3.0e-9,
            cross_section_area_m2: 1.0e-20,
            derivative_step: 1.0e-5,
        }
    }
}

/// The failure stress that a load sweep produced.
#[derive(Clone, Debug, PartialEq)]
pub struct FailureResult {
    /// The failure stress, with the sweep resolution as uncertainty.
    pub failure_stress_pa: Quantity,
    /// The applied strain at failure.
    pub failure_strain: f64,
    /// The index into the bond list of the first bond that failed.
    pub critical_bond: usize,
    /// The local strain of that bond at failure.
    pub critical_bond_strain: f64,
    /// The provenance of the extraction.
    pub provenance: Provenance,
}

/// Builds one bond limit per bond in the system bond-stretch term.
///
/// Every bond gets the same rupture strain. The endpoints and the equilibrium
/// length come from the system. An empty bond term or a non-positive rupture
/// strain returns an error. The caller can edit the returned limits, for
/// example to set the rupture strain of one bond.
pub fn bond_limits_from_system(
    system: &System,
    rupture_strain: f64,
) -> Result<Vec<BondLimit>, ParamError> {
    if !rupture_strain.is_finite() || rupture_strain <= 0.0 {
        return Err(ParamError::InvalidExtraction {
            reason: "the rupture strain must be finite and positive".to_string(),
        });
    }
    if system.bond_count() == 0 {
        return Err(ParamError::InvalidExtraction {
            reason: "the system has no bonds to fail".to_string(),
        });
    }
    let mut bonds = Vec::with_capacity(system.bond_count());
    for bond in 0..system.bond_count() {
        let info = system
            .bond(bond)
            .ok_or_else(|| ParamError::InvalidExtraction {
                reason: "a system bond index is out of range".to_string(),
            })?;
        bonds.push(BondLimit {
            u: info.u,
            v: info.v,
            r0_m: info.r0_m,
            rupture_strain,
        });
    }
    Ok(bonds)
}

/// Extracts a failure stress with the bond limits taken from the system.
///
/// Every bond in the system bond-stretch term gets the same `rupture_strain`.
/// See [`extract_failure_stress`] for the sweep and the result.
pub fn extract_failure_stress_from_system(
    system: &mut System,
    positions_m: &[f64],
    config: &FailureConfig,
    rupture_strain: f64,
) -> Result<FailureResult, ParamError> {
    let bonds = bond_limits_from_system(system, rupture_strain)?;
    extract_failure_stress(system, positions_m, config, &bonds)
}

/// Extracts a failure stress by loading a sample until a bond limit triggers.
///
/// The function applies a growing uniaxial x strain. At each step it measures
/// the local strain of every bond in `bonds`. The first step where a local
/// strain reaches its rupture strain brackets the failure. A bisection then
/// refines the failure strain. The stress there comes from the central
/// difference of the engine energy with respect to the strain, divided by the
/// sample volume.
///
/// The stated uncertainty is half the stress change across one sweep step. It
/// is the strain resolution of the search, not a statistical fit error.
pub fn extract_failure_stress(
    system: &mut System,
    positions_m: &[f64],
    config: &FailureConfig,
    bonds: &[BondLimit],
) -> Result<FailureResult, ParamError> {
    validate_config(system, positions_m, config, bonds)?;
    let volume_m3 = config.reference_length_m * config.cross_section_area_m2;

    let step_count = (config.max_strain / config.strain_step).ceil() as usize;
    let mut previous_strain = 0.0;
    let mut failure_strain = None;
    for step in 1..=step_count {
        let strain = (step as f64 * config.strain_step).min(config.max_strain);
        if exceeds_any_limit(positions_m, strain, bonds)? {
            failure_strain = Some(strain);
            break;
        }
        previous_strain = strain;
        if strain >= config.max_strain {
            break;
        }
    }
    let failure_strain = failure_strain.ok_or_else(|| ParamError::InvalidExtraction {
        reason: format!(
            "no bond failed below the maximum strain {:e}",
            config.max_strain
        ),
    })?;

    let refined_strain =
        refine_failure_strain(positions_m, previous_strain, failure_strain, bonds)?;
    let critical_bond = critical_bond_index(positions_m, refined_strain, bonds)?;
    let critical_bond_strain =
        bond_local_strain(positions_m, refined_strain, bonds[critical_bond])?;

    let failure_stress_pa = stress_from_strain(
        system,
        positions_m,
        refined_strain,
        config.derivative_step,
        volume_m3,
    )?;
    let bracket_stress_pa = stress_from_strain(
        system,
        positions_m,
        previous_strain,
        config.derivative_step,
        volume_m3,
    )?;
    let uncertainty_pa = 0.5 * (failure_stress_pa - bracket_stress_pa).abs();

    let notes = format!(
        "strain sweep to failure with step {:e}; bond {} failed at strain {:e}; \
         stress from the energy derivative at area {:e} m^2",
        config.strain_step, critical_bond, refined_strain, config.cross_section_area_m2
    );
    Ok(FailureResult {
        failure_stress_pa: Quantity::derived(failure_stress_pa, "Pa")?
            .with_uncertainty_si(uncertainty_pa),
        failure_strain: refined_strain,
        critical_bond,
        critical_bond_strain,
        provenance: extraction_provenance(Method::Fit, notes),
    })
}

fn exceeds_any_limit(
    base_positions_m: &[f64],
    strain: f64,
    bonds: &[BondLimit],
) -> Result<bool, ParamError> {
    for bond in bonds {
        if bond_local_strain(base_positions_m, strain, *bond)? >= bond.rupture_strain {
            return Ok(true);
        }
    }
    Ok(false)
}

fn critical_bond_index(
    base_positions_m: &[f64],
    strain: f64,
    bonds: &[BondLimit],
) -> Result<usize, ParamError> {
    let mut best_index = 0;
    let mut best_margin = f64::NEG_INFINITY;
    for (index, bond) in bonds.iter().enumerate() {
        let margin = bond_local_strain(base_positions_m, strain, *bond)? - bond.rupture_strain;
        if margin > best_margin {
            best_margin = margin;
            best_index = index;
        }
    }
    Ok(best_index)
}

fn refine_failure_strain(
    base_positions_m: &[f64],
    mut low: f64,
    mut high: f64,
    bonds: &[BondLimit],
) -> Result<f64, ParamError> {
    for _ in 0..80 {
        let middle = 0.5 * (low + high);
        if exceeds_any_limit(base_positions_m, middle, bonds)? {
            high = middle;
        } else {
            low = middle;
        }
    }
    Ok(high)
}

fn bond_local_strain(
    base_positions_m: &[f64],
    strain: f64,
    bond: BondLimit,
) -> Result<f64, ParamError> {
    if bond.u == bond.v {
        return Err(ParamError::InvalidExtraction {
            reason: "a failure bond joins an atom to itself".to_string(),
        });
    }
    if !bond.r0_m.is_finite() || bond.r0_m <= 0.0 {
        return Err(ParamError::InvalidExtraction {
            reason: "a failure bond has a non-positive equilibrium length".to_string(),
        });
    }
    let strained_m = uniaxial_strain_positions(base_positions_m, strain);
    let base_u = bond.u as usize * 3;
    let base_v = bond.v as usize * 3;
    let dx = strained_m[base_u] - strained_m[base_v];
    let dy = strained_m[base_u + 1] - strained_m[base_v + 1];
    let dz = strained_m[base_u + 2] - strained_m[base_v + 2];
    let distance_m = (dx * dx + dy * dy + dz * dz).sqrt();
    Ok((distance_m - bond.r0_m) / bond.r0_m)
}

fn validate_config(
    system: &System,
    positions_m: &[f64],
    config: &FailureConfig,
    bonds: &[BondLimit],
) -> Result<(), ParamError> {
    if positions_m.len() != 3 * system.atom_count() {
        return Err(EngineError::BufferSizeMismatch {
            len: positions_m.len(),
            expected: 3 * system.atom_count(),
        }
        .into());
    }
    if bonds.is_empty() {
        return Err(ParamError::InvalidExtraction {
            reason: "a failure criterion needs at least one bond".to_string(),
        });
    }
    for bond in bonds {
        if bond.u as usize >= system.atom_count() || bond.v as usize >= system.atom_count() {
            return Err(ParamError::InvalidExtraction {
                reason: "a failure bond index is outside the atom count".to_string(),
            });
        }
    }
    if !config.strain_step.is_finite() || config.strain_step <= 0.0 {
        return Err(ParamError::InvalidExtraction {
            reason: "the strain step must be finite and positive".to_string(),
        });
    }
    if !config.max_strain.is_finite() || config.max_strain <= 0.0 {
        return Err(ParamError::InvalidExtraction {
            reason: "the maximum strain must be finite and positive".to_string(),
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
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::{true_modulus_pa, uniform_chain, AREA_M2, BOND_COUNT, R0_M};

    const RUPTURE_STRAIN: f64 = 0.05;

    fn bonds() -> Vec<BondLimit> {
        (0..BOND_COUNT)
            .map(|atom| BondLimit {
                u: atom as u32,
                v: atom as u32 + 1,
                r0_m: R0_M,
                rupture_strain: RUPTURE_STRAIN,
            })
            .collect()
    }

    fn config() -> FailureConfig {
        FailureConfig {
            strain_step: 1.0e-3,
            max_strain: 0.2,
            reference_length_m: BOND_COUNT as f64 * R0_M,
            cross_section_area_m2: AREA_M2,
            derivative_step: 1.0e-5,
        }
    }

    #[test]
    fn a_known_rupture_strain_gives_the_known_failure_stress() {
        let (mut system, positions_m) = uniform_chain();
        let result =
            extract_failure_stress(&mut system, &positions_m, &config(), &bonds()).expect("fails");
        let expected_pa = true_modulus_pa() * RUPTURE_STRAIN;
        let uncertainty_pa = result
            .failure_stress_pa
            .uncertainty_si()
            .expect("the sweep has an uncertainty");
        let deviation_pa = (result.failure_stress_pa.value_si() - expected_pa).abs();
        println!(
            "failure: extracted {:e} Pa, expected {:e} Pa, uncertainty {:e} Pa, strain {:e}",
            result.failure_stress_pa.value_si(),
            expected_pa,
            uncertainty_pa,
            result.failure_strain
        );
        assert!((result.failure_strain - RUPTURE_STRAIN).abs() < 1.0e-9);
        assert!(
            deviation_pa <= uncertainty_pa,
            "the known failure stress lies outside the stated uncertainty"
        );
    }

    #[test]
    fn no_failure_below_the_maximum_strain_is_an_error() {
        let (mut system, positions_m) = uniform_chain();
        let short = FailureConfig {
            max_strain: 1.0e-2,
            ..config()
        };
        assert!(matches!(
            extract_failure_stress(&mut system, &positions_m, &short, &bonds()),
            Err(ParamError::InvalidExtraction { .. })
        ));
    }

    #[test]
    fn the_system_bond_limits_match_the_hand_built_limits() {
        let (system, _positions_m) = uniform_chain();
        let limits = bond_limits_from_system(&system, RUPTURE_STRAIN).expect("valid limits");
        assert_eq!(limits, bonds());
    }

    #[test]
    fn a_system_driven_failure_run_matches_the_explicit_bond_path() {
        let (mut system, positions_m) = uniform_chain();
        let explicit =
            extract_failure_stress(&mut system, &positions_m, &config(), &bonds()).expect("fails");
        let (mut system, positions_m) = uniform_chain();
        let derived = extract_failure_stress_from_system(
            &mut system,
            &positions_m,
            &config(),
            RUPTURE_STRAIN,
        )
        .expect("fails");
        assert_eq!(derived.critical_bond, explicit.critical_bond);
        assert_eq!(derived.failure_strain, explicit.failure_strain);
        assert_eq!(
            derived.failure_stress_pa.value_si(),
            explicit.failure_stress_pa.value_si()
        );
    }

    #[test]
    fn a_system_without_bonds_or_a_bad_rupture_strain_is_rejected() {
        let (mut system, positions_m) = uniform_chain();
        assert!(matches!(
            bond_limits_from_system(&system, 0.0),
            Err(ParamError::InvalidExtraction { .. })
        ));
        let empty = System::with_masses_kg(2, &[1.0e-26, 1.0e-26]).expect("valid");
        assert!(matches!(
            bond_limits_from_system(&empty, RUPTURE_STRAIN),
            Err(ParamError::InvalidExtraction { .. })
        ));
        assert!(matches!(
            extract_failure_stress_from_system(&mut system, &positions_m, &config(), -0.1),
            Err(ParamError::InvalidExtraction { .. })
        ));
    }

    #[test]
    fn a_heterogeneous_criterion_reports_the_weakest_bond() {
        let (mut system, positions_m) = uniform_chain();
        let mut mixed = bonds();
        mixed[7].rupture_strain = 0.01;
        mixed[3].rupture_strain = 0.02;
        let result =
            extract_failure_stress(&mut system, &positions_m, &config(), &mixed).expect("fails");
        assert_eq!(result.critical_bond, 7);
        assert!((result.failure_strain - 0.01).abs() < 1.0e-9);
    }
}
