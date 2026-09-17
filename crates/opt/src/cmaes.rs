//! A compact CMA-ES minimizer for a small parameter vector.
//!
//! The covariance matrix adaptation evolution strategy samples a population
//! around a mean, keeps the best half, and adapts the step size and the
//! covariance from the survivors. It handles a rugged, non-differentiable
//! cost, which is what a discrete generator score is.
//!
//! The search is deterministic for a given seed, so a result reproduces.
//! Every candidate is clamped to its bounds, and a rejected candidate is
//! dropped. A parameter that must be a whole number is rounded.

use std::f64::consts::TAU;

/// The number of generations without a new best cost before the search stops.
const STALL_GENERATIONS: usize = 8;

/// The range and the type of one parameter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ParameterBounds {
    /// The inclusive lower bound.
    pub min: f64,
    /// The inclusive upper bound.
    pub max: f64,
    /// Whether the value must be a whole number.
    pub integer: bool,
}

impl ParameterBounds {
    /// Builds a parameter range.
    pub const fn new(min: f64, max: f64, integer: bool) -> Self {
        Self { min, max, integer }
    }

    /// Rounds an integer parameter and clamps it to the range.
    pub fn clamp(&self, value: f64) -> f64 {
        let rounded = if self.integer { value.round() } else { value };
        rounded.max(self.min).min(self.max)
    }

    /// Returns the width of the range.
    pub fn span(&self) -> f64 {
        self.max - self.min
    }
}

/// A cost function over a parameter vector. A lower cost is better.
///
/// A return of `None` marks a rejected candidate, for example a parameter set
/// the generator refuses. A blanket implementation covers a closure.
pub trait Objective {
    /// Returns the cost of one parameter vector, or `None` when it is invalid.
    fn evaluate(&self, parameters: &[f64]) -> Option<f64>;
}

impl<F> Objective for F
where
    F: Fn(&[f64]) -> Option<f64>,
{
    fn evaluate(&self, parameters: &[f64]) -> Option<f64> {
        self(parameters)
    }
}

/// The settings of one search.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CmaEsOptions {
    /// The number of candidates in one generation.
    pub population: usize,
    /// The maximum number of generations.
    pub generations: usize,
    /// The initial step size, as a fraction of the mean parameter range.
    pub sigma_fraction: f64,
    /// The random seed. The same seed gives the same result.
    pub seed: u64,
}

impl Default for CmaEsOptions {
    fn default() -> Self {
        Self {
            population: 8,
            generations: 12,
            sigma_fraction: 0.25,
            seed: 1,
        }
    }
}

/// The result of one search.
#[derive(Clone, Debug, PartialEq)]
pub struct CmaEsReport {
    /// The best parameter vector found.
    pub best_parameters: Vec<f64>,
    /// The cost of the best vector. It is infinite when nothing was valid.
    pub best_cost: f64,
    /// The number of generations run.
    pub generations: usize,
    /// The number of candidates evaluated.
    pub evaluations: usize,
    /// The number of candidates the objective rejected.
    pub rejected: usize,
    /// The best cost after each generation.
    pub history: Vec<f64>,
    /// True when the best cost stopped improving for `STALL_GENERATIONS`
    /// generations before the end.
    pub converged: bool,
}

/// Runs CMA-ES from an initial vector.
///
/// The objective returns `None` for an invalid candidate. The result is the
/// best valid candidate. When no candidate is valid, the best cost is
/// infinite and the best vector is the initial one.
pub fn minimize<O: Objective>(
    initial: &[f64],
    bounds: &[ParameterBounds],
    options: &CmaEsOptions,
    objective: &O,
) -> CmaEsReport {
    let dimension = initial.len();
    if dimension == 0 || bounds.len() != dimension {
        return CmaEsReport {
            best_parameters: initial.to_vec(),
            best_cost: f64::INFINITY,
            generations: 0,
            evaluations: 0,
            rejected: 0,
            history: Vec::new(),
            converged: false,
        };
    }

    let population = options.population.max(4);
    let mu = population / 2;
    let weights = weights_of(mu);
    let mu_eff = 1.0 / weights.iter().map(|w| w * w).sum::<f64>();
    let strategy = Strategy::new(dimension, mu_eff);

    let mean_span: f64 = bounds.iter().map(|b| b.span()).sum::<f64>() / dimension as f64;
    let mut sigma = (options.sigma_fraction * mean_span).max(f64::EPSILON);
    let mut mean = initial.to_vec();
    let mut covariance = identity(dimension);
    let mut path_sigma = vec![0.0; dimension];
    let mut path_covariance = vec![0.0; dimension];

    let mut rng = Rng::new(options.seed);
    let mut best_parameters = initial.to_vec();
    let mut best_cost = f64::INFINITY;
    let mut evaluations = 0usize;
    let mut rejected = 0usize;
    let mut history = Vec::new();
    let mut converged = false;
    let mut stall = 0usize;
    let mut generation = 0usize;

    while generation < options.generations {
        let (eigenvalues, eigenvectors) = jacobi_eigen(&covariance, dimension);
        let deviations: Vec<f64> = eigenvalues
            .iter()
            .map(|value| value.max(1.0e-14).sqrt())
            .collect();

        let mut candidates: Vec<(f64, Vec<f64>, Vec<f64>)> = Vec::new();
        for _ in 0..population {
            let z: Vec<f64> = (0..dimension).map(|_| rng.normal()).collect();
            let scaled: Vec<f64> = (0..dimension)
                .map(|axis| deviations[axis] * z[axis])
                .collect();
            let rotated = apply(&eigenvectors, &scaled, dimension);
            let candidate: Vec<f64> = (0..dimension)
                .map(|index| bounds[index].clamp(mean[index] + sigma * rotated[index]))
                .collect();
            let step: Vec<f64> = (0..dimension)
                .map(|index| (candidate[index] - mean[index]) / sigma)
                .collect();
            evaluations += 1;
            match objective.evaluate(&candidate) {
                Some(cost) if cost.is_finite() => candidates.push((cost, candidate, step)),
                _ => rejected += 1,
            }
        }

        if candidates.is_empty() {
            break;
        }
        candidates.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        if candidates[0].0 < best_cost {
            best_cost = candidates[0].0;
            best_parameters = candidates[0].1.clone();
            stall = 0;
        } else {
            stall += 1;
            if stall >= STALL_GENERATIONS {
                converged = true;
            }
        }
        history.push(best_cost);
        generation += 1;

        let selected = mu.min(candidates.len());
        let mut new_mean = vec![0.0; dimension];
        let mut weight_sum = 0.0;
        for index in 0..selected {
            let weight = weights[index];
            weight_sum += weight;
            for (axis, value) in new_mean.iter_mut().enumerate() {
                *value += weight * candidates[index].1[axis];
            }
        }
        for value in new_mean.iter_mut() {
            *value /= weight_sum;
        }
        let step_mean: Vec<f64> = (0..dimension)
            .map(|axis| (new_mean[axis] - mean[axis]) / sigma)
            .collect();

        let inverse_sqrt_step = {
            let rotated = apply_transpose(&eigenvectors, &step_mean, dimension);
            let solved: Vec<f64> = (0..dimension)
                .map(|axis| rotated[axis] / deviations[axis].max(1.0e-14))
                .collect();
            apply(&eigenvectors, &solved, dimension)
        };
        for axis in 0..dimension {
            path_sigma[axis] = (1.0 - strategy.c_sigma) * path_sigma[axis]
                + strategy.sigma_scale * inverse_sqrt_step[axis];
        }
        let path_norm = norm(&path_sigma);
        let normalizer = (1.0 - (1.0 - strategy.c_sigma).powi(2 * generation as i32)).sqrt();
        let h_sigma = if normalizer > 0.0
            && path_norm / normalizer < (1.4 + 2.0 / (dimension as f64 + 1.0)) * strategy.chi_n
        {
            1.0
        } else {
            0.0
        };
        for axis in 0..dimension {
            path_covariance[axis] = (1.0 - strategy.c_covariance) * path_covariance[axis]
                + h_sigma * strategy.covariance_scale * step_mean[axis];
        }

        let decay = 1.0 - strategy.c_rank_one - strategy.c_rank_mu;
        let mut next = vec![0.0; dimension * dimension];
        for row in 0..dimension {
            for column in 0..dimension {
                let index = row * dimension + column;
                next[index] = decay * covariance[index]
                    + strategy.c_rank_one
                        * (path_covariance[row] * path_covariance[column]
                            + (1.0 - h_sigma)
                                * strategy.c_covariance
                                * (2.0 - strategy.c_covariance)
                                * covariance[index]);
            }
        }
        for index in 0..selected {
            let weight = weights[index];
            for row in 0..dimension {
                for column in 0..dimension {
                    next[row * dimension + column] +=
                        weight * candidates[index].2[row] * candidates[index].2[column];
                }
            }
        }
        covariance = next;
        sigma *= (strategy.c_sigma / strategy.d_sigma * (path_norm / strategy.chi_n - 1.0)).exp();
        sigma = sigma.clamp(1.0e-12, mean_span);
        mean = new_mean;
        if converged {
            break;
        }
    }

    CmaEsReport {
        best_parameters,
        best_cost,
        generations: generation,
        evaluations,
        rejected,
        history,
        converged,
    }
}

/// The derived CMA-ES scalars.
struct Strategy {
    c_sigma: f64,
    d_sigma: f64,
    c_covariance: f64,
    c_rank_one: f64,
    c_rank_mu: f64,
    sigma_scale: f64,
    covariance_scale: f64,
    chi_n: f64,
}

impl Strategy {
    fn new(dimension: usize, mu_eff: f64) -> Self {
        let n = dimension as f64;
        let c_sigma = (mu_eff + 2.0) / (n + mu_eff + 5.0);
        let d_sigma = 1.0 + 2.0 * (((mu_eff - 1.0) / (n + 1.0)).sqrt() - 1.0).max(0.0) + c_sigma;
        let c_covariance = (4.0 + mu_eff / n) / (n + 4.0 + 2.0 * mu_eff / n);
        let c_rank_one = 2.0 / ((n + 1.3).powi(2) + mu_eff);
        let c_rank_mu = (1.0 - c_rank_one)
            .min(2.0 * (mu_eff - 2.0 + 1.0 / mu_eff) / ((n + 2.0).powi(2) + mu_eff));
        let sigma_scale = (c_sigma * (2.0 - c_sigma) * mu_eff).sqrt();
        let covariance_scale = (c_covariance * (2.0 - c_covariance) * mu_eff).sqrt();
        let chi_n = n.sqrt() * (1.0 - 1.0 / (4.0 * n) + 1.0 / (21.0 * n * n));
        Self {
            c_sigma,
            d_sigma,
            c_covariance,
            c_rank_one,
            c_rank_mu,
            sigma_scale,
            covariance_scale,
            chi_n,
        }
    }
}

/// Returns the positive, normalized selection weights.
fn weights_of(mu: usize) -> Vec<f64> {
    let raw: Vec<f64> = (0..mu)
        .map(|index| ((mu as f64 + 0.5).ln() - ((index + 1) as f64).ln()).max(0.0))
        .collect();
    let total: f64 = raw.iter().sum();
    if total <= 0.0 {
        return vec![1.0 / mu as f64; mu];
    }
    raw.iter().map(|value| value / total).collect()
}

/// Returns an identity matrix in row-major order.
fn identity(dimension: usize) -> Vec<f64> {
    let mut matrix = vec![0.0; dimension * dimension];
    for index in 0..dimension {
        matrix[index * dimension + index] = 1.0;
    }
    matrix
}

/// Applies a row-major matrix to a vector.
fn apply(matrix: &[f64], vector: &[f64], dimension: usize) -> Vec<f64> {
    (0..dimension)
        .map(|row| {
            (0..dimension)
                .map(|column| matrix[row * dimension + column] * vector[column])
                .sum()
        })
        .collect()
}

/// Applies the transpose of a row-major matrix to a vector.
fn apply_transpose(matrix: &[f64], vector: &[f64], dimension: usize) -> Vec<f64> {
    (0..dimension)
        .map(|column| {
            (0..dimension)
                .map(|row| matrix[row * dimension + column] * vector[row])
                .sum()
        })
        .collect()
}

/// Returns the Euclidean norm.
fn norm(vector: &[f64]) -> f64 {
    vector.iter().map(|value| value * value).sum::<f64>().sqrt()
}

/// Returns the eigenvalues and the eigenvectors of a symmetric matrix.
///
/// The input is row-major and is not modified. The eigenvalues come out in
/// ascending order, and eigenvector `k` is column `k` of the second result.
fn jacobi_eigen(matrix: &[f64], dimension: usize) -> (Vec<f64>, Vec<f64>) {
    let mut a = matrix.to_vec();
    let mut vectors = identity(dimension);
    for _ in 0..100 {
        let mut off = 0.0;
        for row in 0..dimension {
            for column in (row + 1)..dimension {
                off += a[row * dimension + column] * a[row * dimension + column];
            }
        }
        if off.sqrt() < 1.0e-14 {
            break;
        }
        for p in 0..dimension {
            for q in (p + 1)..dimension {
                let apq = a[p * dimension + q];
                if apq.abs() < 1.0e-30 {
                    continue;
                }
                let app = a[p * dimension + p];
                let aqq = a[q * dimension + q];
                let theta = (aqq - app) / (2.0 * apq);
                let sign = if theta >= 0.0 { 1.0 } else { -1.0 };
                let t = sign / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..dimension {
                    let akp = a[k * dimension + p];
                    let akq = a[k * dimension + q];
                    a[k * dimension + p] = c * akp - s * akq;
                    a[k * dimension + q] = s * akp + c * akq;
                }
                for k in 0..dimension {
                    let apk = a[p * dimension + k];
                    let aqk = a[q * dimension + k];
                    a[p * dimension + k] = c * apk - s * aqk;
                    a[q * dimension + k] = s * apk + c * aqk;
                }
                for k in 0..dimension {
                    let vkp = vectors[k * dimension + p];
                    let vkq = vectors[k * dimension + q];
                    vectors[k * dimension + p] = c * vkp - s * vkq;
                    vectors[k * dimension + q] = s * vkp + c * vkq;
                }
            }
        }
    }
    let mut order: Vec<usize> = (0..dimension).collect();
    let diagonal: Vec<f64> = (0..dimension).map(|i| a[i * dimension + i]).collect();
    order.sort_by(|a, b| {
        diagonal[*a]
            .partial_cmp(&diagonal[*b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let eigenvalues: Vec<f64> = order.iter().map(|i| diagonal[*i]).collect();
    let mut sorted = vec![0.0; dimension * dimension];
    for (target, source) in order.iter().enumerate() {
        for row in 0..dimension {
            sorted[row * dimension + target] = vectors[row * dimension + *source];
        }
    }
    (eigenvalues, sorted)
}

/// A deterministic SplitMix64 generator with a normal sampler.
struct Rng {
    state: u64,
}

impl Rng {
    fn new(seed: u64) -> Self {
        Self {
            state: seed ^ 0x9E37_79B9_7F4A_7C15,
        }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / 9_007_199_254_740_992.0
    }

    fn normal(&mut self) -> f64 {
        let u1 = self.uniform().max(1.0e-12);
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (TAU * u2).cos()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_bowl_is_minimized() {
        let bounds = vec![ParameterBounds::new(-5.0, 5.0, false); 3];
        let target = [0.3, -1.0, 2.0];
        let objective = |x: &[f64]| {
            Some(
                x.iter()
                    .zip(target.iter())
                    .map(|(value, want)| (value - want) * (value - want))
                    .sum::<f64>(),
            )
        };
        let options = CmaEsOptions {
            population: 12,
            generations: 120,
            ..Default::default()
        };
        let report = minimize(&[0.0, 0.0, 0.0], &bounds, &options, &objective);
        assert!(report.best_cost < 1.0e-3, "cost is {}", report.best_cost);
        for (value, want) in report.best_parameters.iter().zip(target.iter()) {
            assert!((value - want).abs() < 0.05, "{value} is not {want}");
        }
        assert!(report.history.first().copied().unwrap_or(0.0) > report.best_cost);
    }

    #[test]
    fn an_integer_parameter_stays_whole() {
        let bounds = vec![ParameterBounds::new(0.0, 10.0, true)];
        let objective = |x: &[f64]| Some((x[0] - 3.4) * (x[0] - 3.4));
        let options = CmaEsOptions {
            population: 8,
            generations: 30,
            ..Default::default()
        };
        let report = minimize(&[5.0], &bounds, &options, &objective);
        assert_eq!(report.best_parameters[0].fract(), 0.0);
        assert!(
            (report.best_parameters[0] - 3.0).abs() < 0.5,
            "the integer result is {}",
            report.best_parameters[0]
        );
    }

    #[test]
    fn a_rejected_candidate_is_skipped() {
        let bounds = vec![ParameterBounds::new(0.0, 1.0, false)];
        let objective = |x: &[f64]| {
            if x[0] > 0.5 {
                None
            } else {
                Some((x[0] - 0.2) * (x[0] - 0.2))
            }
        };
        let options = CmaEsOptions {
            population: 8,
            generations: 20,
            ..Default::default()
        };
        let report = minimize(&[0.3], &bounds, &options, &objective);
        assert!(report.rejected > 0, "nothing was rejected");
        assert!(report.best_parameters[0] <= 0.5);
        assert!(report.best_cost < 0.05, "cost is {}", report.best_cost);
    }

    #[test]
    fn an_empty_problem_is_refused() {
        let report = minimize(&[], &[], &CmaEsOptions::default(), &|_: &[f64]| Some(0.0));
        assert_eq!(report.generations, 0);
        assert!(report.best_cost.is_infinite());
    }
}
