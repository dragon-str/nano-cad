use crate::error::EngineError;
use crate::system::System;

const ARMIJO_C1: f64 = 1.0e-4;
const CURVATURE_C2: f64 = 0.4;
const MAX_LINE_SEARCH_ITERATIONS: usize = 80;
const MAX_ALPHA_M: f64 = 1.0e-8;

/// Options for the conjugate-gradient minimizer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MinimizeOptions {
    /// Maximum number of outer iterations.
    pub max_iterations: usize,
    /// Convergence threshold on the infinity norm of the gradient, in newtons.
    pub gradient_tolerance_n: f64,
    /// First trial line-search step along a unit direction, in metres.
    pub initial_step_m: f64,
}

impl Default for MinimizeOptions {
    fn default() -> Self {
        Self {
            max_iterations: 200,
            gradient_tolerance_n: 1.0e-14,
            initial_step_m: 1.0e-13,
        }
    }
}

/// The outcome of a minimization run.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MinimizeResult {
    /// The final potential energy in joules.
    pub energy_j: f64,
    /// The final gradient infinity norm in newtons.
    pub gradient_norm_n: f64,
    /// The number of outer iterations that ran.
    pub iterations: usize,
    /// Whether the gradient norm reached the tolerance.
    pub converged: bool,
}

/// Minimizes the total potential of a system in place.
///
/// The method is nonlinear conjugate gradient with a Polak-Ribiere+ update and
/// a strong-Wolfe line search. `positions_m` is the start geometry and receives
/// the result. The direction is normalized before the search, so the line-search
/// step is a length in metres.
///
/// The function returns the final energy, the final gradient norm, the
/// iteration count, and whether the gradient norm reached the tolerance. An
/// invalid option or a line-search failure returns an error. It never panics.
pub fn minimize(
    system: &mut System,
    positions_m: &mut [f64],
    options: &MinimizeOptions,
) -> Result<MinimizeResult, EngineError> {
    validate_options(options)?;
    if positions_m.len() != 3 * system.atom_count() {
        return Err(EngineError::BufferSizeMismatch {
            len: positions_m.len(),
            expected: 3 * system.atom_count(),
        });
    }

    let coordinate_count = positions_m.len();
    let (mut energy_j, mut gradient) = system.energy_and_gradient_j(positions_m)?;
    let mut gradient_norm_n = infinity_norm(&gradient);
    if gradient_norm_n <= options.gradient_tolerance_n {
        return Ok(MinimizeResult {
            energy_j,
            gradient_norm_n,
            iterations: 0,
            converged: true,
        });
    }

    let mut direction = vec![0.0; coordinate_count];
    let mut previous_gradient = vec![0.0; coordinate_count];
    let mut iterations = 0;
    let mut converged = false;

    while iterations < options.max_iterations {
        iterations += 1;
        let restart = iterations == 1 || (iterations - 1) % coordinate_count == 0;
        if restart {
            direction.copy_from_slice(&gradient);
            for value in &mut direction {
                *value = -*value;
            }
        } else {
            let beta = polak_ribiere(&gradient, &previous_gradient);
            for coordinate in 0..coordinate_count {
                direction[coordinate] = -gradient[coordinate] + beta * direction[coordinate];
            }
        }

        let mut slope = dot(&gradient, &direction);
        if !slope.is_finite() || slope >= 0.0 {
            direction.copy_from_slice(&gradient);
            for value in &mut direction {
                *value = -*value;
            }
            slope = dot(&gradient, &direction);
        }

        let norm = l2_norm(&direction);
        if norm <= 0.0 {
            break;
        }
        let search_direction: Vec<f64> = direction.iter().map(|value| value / norm).collect();
        let slope_unit = slope / norm;

        let outcome = line_search(
            system,
            positions_m,
            &search_direction,
            energy_j,
            slope_unit,
            options.initial_step_m,
            iterations,
        )?;

        previous_gradient.copy_from_slice(&gradient);
        for coordinate in 0..coordinate_count {
            positions_m[coordinate] += outcome.alpha_m * search_direction[coordinate];
        }
        energy_j = outcome.energy_j;
        gradient = outcome.gradient;
        gradient_norm_n = infinity_norm(&gradient);
        if gradient_norm_n <= options.gradient_tolerance_n {
            converged = true;
            break;
        }
    }

    Ok(MinimizeResult {
        energy_j,
        gradient_norm_n,
        iterations,
        converged,
    })
}

struct LineSearchOutcome {
    alpha_m: f64,
    energy_j: f64,
    gradient: Vec<f64>,
}

#[allow(clippy::too_many_arguments)]
fn line_search(
    system: &mut System,
    positions_m: &[f64],
    direction: &[f64],
    energy0_j: f64,
    slope0: f64,
    initial_step_m: f64,
    iteration: usize,
) -> Result<LineSearchOutcome, EngineError> {
    if !slope0.is_finite() || slope0 >= 0.0 {
        return Err(EngineError::LineSearchFailed { iteration });
    }
    let mut trial = vec![0.0; positions_m.len()];
    let mut alpha_previous_m = 0.0;
    let mut energy_previous_j = energy0_j;
    let mut alpha_m = initial_step_m;

    for _ in 0..MAX_LINE_SEARCH_ITERATIONS {
        let (energy_j, gradient) = evaluate(system, positions_m, direction, alpha_m, &mut trial)?;
        if energy_j > energy0_j + ARMIJO_C1 * alpha_m * slope0
            || (alpha_previous_m > 0.0 && energy_j >= energy_previous_j)
        {
            return zoom(
                system,
                positions_m,
                direction,
                energy0_j,
                slope0,
                alpha_previous_m,
                alpha_m,
                &mut trial,
                iteration,
            );
        }
        let slope = dot(&gradient, direction);
        if slope.abs() <= -CURVATURE_C2 * slope0 {
            return Ok(LineSearchOutcome {
                alpha_m,
                energy_j,
                gradient,
            });
        }
        if slope >= 0.0 {
            return zoom(
                system,
                positions_m,
                direction,
                energy0_j,
                slope0,
                alpha_m,
                alpha_previous_m.max(alpha_m * 0.5),
                &mut trial,
                iteration,
            );
        }
        alpha_previous_m = alpha_m;
        energy_previous_j = energy_j;
        alpha_m = (alpha_m * 2.0).min(MAX_ALPHA_M);
    }
    Err(EngineError::LineSearchFailed { iteration })
}

#[allow(clippy::too_many_arguments)]
fn zoom(
    system: &mut System,
    positions_m: &[f64],
    direction: &[f64],
    energy0_j: f64,
    slope0: f64,
    mut low_m: f64,
    mut high_m: f64,
    trial: &mut [f64],
    iteration: usize,
) -> Result<LineSearchOutcome, EngineError> {
    if low_m > high_m {
        std::mem::swap(&mut low_m, &mut high_m);
    }
    let mut energy_low_j = energy0_j;
    for _ in 0..MAX_LINE_SEARCH_ITERATIONS {
        let alpha_m = 0.5 * (low_m + high_m);
        let (energy_j, gradient) = evaluate(system, positions_m, direction, alpha_m, trial)?;
        if energy_j > energy0_j + ARMIJO_C1 * alpha_m * slope0 || energy_j >= energy_low_j {
            high_m = alpha_m;
        } else {
            let slope = dot(&gradient, direction);
            if slope.abs() <= -CURVATURE_C2 * slope0 {
                return Ok(LineSearchOutcome {
                    alpha_m,
                    energy_j,
                    gradient,
                });
            }
            if slope * (high_m - low_m) >= 0.0 {
                high_m = low_m;
            }
            low_m = alpha_m;
            energy_low_j = energy_j;
        }
    }
    Err(EngineError::LineSearchFailed { iteration })
}

fn evaluate(
    system: &mut System,
    positions_m: &[f64],
    direction: &[f64],
    alpha_m: f64,
    trial: &mut [f64],
) -> Result<(f64, Vec<f64>), EngineError> {
    for coordinate in 0..positions_m.len() {
        trial[coordinate] = positions_m[coordinate] + alpha_m * direction[coordinate];
    }
    system.energy_and_gradient_j(trial)
}

fn validate_options(options: &MinimizeOptions) -> Result<(), EngineError> {
    if options.max_iterations == 0 {
        return Err(EngineError::InvalidMaxIterations);
    }
    if !options.gradient_tolerance_n.is_finite() || options.gradient_tolerance_n <= 0.0 {
        return Err(EngineError::InvalidGradientTolerance {
            tolerance_n: options.gradient_tolerance_n,
        });
    }
    if !options.initial_step_m.is_finite() || options.initial_step_m <= 0.0 {
        return Err(EngineError::InvalidInitialStep {
            step_m: options.initial_step_m,
        });
    }
    Ok(())
}

fn polak_ribiere(gradient: &[f64], previous_gradient: &[f64]) -> f64 {
    let mut numerator = 0.0;
    let mut denominator = 0.0;
    for coordinate in 0..gradient.len() {
        numerator += gradient[coordinate] * (gradient[coordinate] - previous_gradient[coordinate]);
        denominator += previous_gradient[coordinate] * previous_gradient[coordinate];
    }
    if denominator > 0.0 {
        (numerator / denominator).max(0.0)
    } else {
        0.0
    }
}

fn l2_norm(values: &[f64]) -> f64 {
    let mut norm_sq = 0.0;
    for value in values {
        norm_sq += value * value;
    }
    norm_sq.sqrt()
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    let mut sum = 0.0;
    for index in 0..a.len() {
        sum += a[index] * b[index];
    }
    sum
}

fn infinity_norm(values: &[f64]) -> f64 {
    let mut norm = 0.0;
    for value in values {
        if value.abs() > norm {
            norm = value.abs();
        }
    }
    norm
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::{full_system, small_molecule_positions_m, CARBON_MASS_KG};

    #[test]
    fn invalid_options_are_rejected() {
        let mut system = System::new(0);
        let mut positions_m: Vec<f64> = Vec::new();
        let no_iterations = MinimizeOptions {
            max_iterations: 0,
            ..MinimizeOptions::default()
        };
        assert!(matches!(
            minimize(&mut system, &mut positions_m, &no_iterations),
            Err(EngineError::InvalidMaxIterations)
        ));
        let bad_tolerance = MinimizeOptions {
            gradient_tolerance_n: 0.0,
            ..MinimizeOptions::default()
        };
        assert!(matches!(
            minimize(&mut system, &mut positions_m, &bad_tolerance),
            Err(EngineError::InvalidGradientTolerance { .. })
        ));
        let bad_step = MinimizeOptions {
            initial_step_m: -1.0,
            ..MinimizeOptions::default()
        };
        assert!(matches!(
            minimize(&mut system, &mut positions_m, &bad_step),
            Err(EngineError::InvalidInitialStep { .. })
        ));
    }

    #[test]
    fn a_wrong_position_buffer_is_rejected() {
        let mut system =
            System::with_masses_kg(2, &[CARBON_MASS_KG, CARBON_MASS_KG]).expect("valid");
        let mut positions_m = vec![0.0; 3];
        assert!(matches!(
            minimize(&mut system, &mut positions_m, &MinimizeOptions::default()),
            Err(EngineError::BufferSizeMismatch { .. })
        ));
    }

    #[test]
    fn a_diatomic_reaches_its_equilibrium_length() {
        let r0_m = 1.5e-10;
        let mut system =
            System::with_masses_kg(2, &[CARBON_MASS_KG, CARBON_MASS_KG]).expect("valid");
        let mut bonds = crate::bond_stretch::BondStretchTerm::new();
        bonds.add_bond(0, 1, 300.0, r0_m).expect("valid bond");
        system.set_bond_stretch(bonds).expect("fits");
        let mut positions_m = vec![-0.6 * r0_m, 0.0, 0.0, 0.6 * r0_m, 0.0, 0.0];

        let result = minimize(
            &mut system,
            &mut positions_m,
            &MinimizeOptions {
                gradient_tolerance_n: 1.0e-16,
                ..MinimizeOptions::default()
            },
        )
        .expect("valid minimization");

        let dx = positions_m[3] - positions_m[0];
        let dy = positions_m[4] - positions_m[1];
        let dz = positions_m[5] - positions_m[2];
        let distance_m = (dx * dx + dy * dy + dz * dz).sqrt();
        assert!(result.converged, "gradient norm {}", result.gradient_norm_n);
        assert!((distance_m - r0_m).abs() / r0_m < 1.0e-6);
    }

    #[test]
    fn a_perturbed_molecule_reaches_a_lower_energy_with_a_small_gradient() {
        let positions_m = small_molecule_positions_m();
        let masses = vec![CARBON_MASS_KG; 6];
        let mut system = full_system(&positions_m, &masses);
        let start_energy_j = system.energy_j(&positions_m).expect("valid");

        let mut perturbed_m = positions_m.clone();
        let mut rng = crate::test_support::Rng::new(0x2B1A5);
        for value in &mut perturbed_m {
            *value += rng.symmetric(2.0e-12);
        }
        let perturbed_energy_j = system.energy_j(&perturbed_m).expect("valid");
        assert!(perturbed_energy_j > start_energy_j);

        let options = MinimizeOptions {
            max_iterations: 400,
            gradient_tolerance_n: 1.0e-14,
            initial_step_m: 1.0e-13,
        };
        let result = minimize(&mut system, &mut perturbed_m, &options).expect("valid minimization");
        println!(
            "minimizer: start {start_energy_j:e} J, perturbed {perturbed_energy_j:e} J, final {:e} J, gradient {:e} N, iterations {}",
            result.energy_j, result.gradient_norm_n, result.iterations
        );

        assert!(result.energy_j < perturbed_energy_j);
        assert!(result.gradient_norm_n <= options.gradient_tolerance_n);
        assert!(result.converged);
    }
}
