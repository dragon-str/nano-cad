//! Shared strain helpers for stiffness and failure extraction.
//!
//! The strain is a uniaxial stretch of the x axis about the sample centroid.
//! The y and z coordinates are held fixed, so the extraction measures a
//! longitudinal response at constant lateral dimensions. This is a stated
//! model, not a full elastic tensor.

use nanocad_engine::{EngineError, System};

/// Returns a copy of the positions with a uniaxial strain on the x axis.
///
/// The strain is applied about the mean x coordinate, so the sample does not
/// translate. The lateral coordinates do not relax.
pub(crate) fn uniaxial_strain_positions(base_positions_m: &[f64], strain: f64) -> Vec<f64> {
    let atom_count = base_positions_m.len() / 3;
    let mut positions_m = base_positions_m.to_vec();
    if atom_count == 0 {
        return positions_m;
    }
    let x_mean_m: f64 = (0..atom_count)
        .map(|atom| base_positions_m[3 * atom])
        .sum::<f64>()
        / atom_count as f64;
    for atom in 0..atom_count {
        let x_m = base_positions_m[3 * atom];
        positions_m[3 * atom] = x_mean_m + (x_m - x_mean_m) * (1.0 + strain);
    }
    positions_m
}

/// Returns the axial stress in pascals at a strain.
///
/// The stress follows `stress = (1 / V) dU / d(strain)`, with `V` the sample
/// volume. A central difference of the engine energy gives the derivative. The
/// volume is `cross_section_area * reference_length`.
pub(crate) fn stress_from_strain(
    system: &mut System,
    base_positions_m: &[f64],
    strain: f64,
    derivative_step: f64,
    volume_m3: f64,
) -> Result<f64, EngineError> {
    let forward_m = uniaxial_strain_positions(base_positions_m, strain + derivative_step);
    let backward_m = uniaxial_strain_positions(base_positions_m, strain - derivative_step);
    let forward_j = system.energy_j(&forward_m)?;
    let backward_j = system.energy_j(&backward_m)?;
    Ok((forward_j - backward_j) / (2.0 * derivative_step) / volume_m3)
}
