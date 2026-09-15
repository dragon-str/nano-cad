//! Statistics and deterministic noise helpers for parameter extraction.
//!
//! These helpers are internal. They have no physical meaning on their own.

/// A small deterministic xorshift random generator.
///
/// The generator produces synthetic measurement noise. It is seeded, so a run
/// repeats exactly. A zero seed is forced odd, because the xorshift state must
/// never be zero.
pub(crate) struct Rng(u64);

impl Rng {
    /// Builds a generator from a seed.
    pub(crate) const fn new(seed: u64) -> Self {
        Self(seed | 1)
    }

    /// Returns the next value in `[0, 1)`.
    pub(crate) fn next_f64(&mut self) -> f64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        (x >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Returns a value in `[-amplitude, amplitude)`.
    pub(crate) fn uniform_symmetric(&mut self, amplitude: f64) -> f64 {
        (self.next_f64() - 0.5) * 2.0 * amplitude
    }
}

/// Returns the arithmetic mean, or [`None`] for an empty slice.
pub(crate) fn mean(values: &[f64]) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    Some(values.iter().sum::<f64>() / values.len() as f64)
}

/// Returns the sample standard deviation, or [`None`] below two values.
pub(crate) fn sample_std_dev(values: &[f64]) -> Option<f64> {
    if values.len() < 2 {
        return None;
    }
    let average = mean(values)?;
    let sum_squares: f64 = values
        .iter()
        .map(|value| (value - average) * (value - average))
        .sum();
    Some((sum_squares / (values.len() as f64 - 1.0)).sqrt())
}

/// A least-squares straight-line fit `y = intercept + slope * x`.
pub(crate) struct LineFit {
    /// The fitted slope.
    pub(crate) slope: f64,
    /// The fitted intercept.
    pub(crate) intercept: f64,
    /// The standard error of the slope.
    pub(crate) slope_standard_error: f64,
    /// The coefficient of determination.
    pub(crate) r_squared: f64,
}

/// Fits a straight line by least squares.
///
/// It needs at least three points to leave one residual degree of freedom, so
/// that the slope standard error is defined. A degenerate x spread returns
/// [`None`].
pub(crate) fn fit_line(x: &[f64], y: &[f64]) -> Option<LineFit> {
    if x.len() < 3 || y.len() != x.len() {
        return None;
    }
    let x_mean = mean(x)?;
    let y_mean = mean(y)?;
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    let mut syy = 0.0;
    for (x_value, y_value) in x.iter().zip(y.iter()) {
        let dx = x_value - x_mean;
        let dy = y_value - y_mean;
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }
    if sxx <= 0.0 {
        return None;
    }
    let slope = sxy / sxx;
    let intercept = y_mean - slope * x_mean;
    let sum_squared_residuals = (syy - slope * sxy).max(0.0);
    let residual_variance = sum_squared_residuals / (x.len() as f64 - 2.0);
    let slope_standard_error = (residual_variance / sxx).sqrt();
    let r_squared = if syy > 0.0 {
        (1.0 - sum_squared_residuals / syy).clamp(0.0, 1.0)
    } else {
        1.0
    };
    Some(LineFit {
        slope,
        intercept,
        slope_standard_error,
        r_squared,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mean_and_spread_match_a_known_sample() {
        let values = [1.0, 2.0, 3.0, 4.0];
        assert_eq!(mean(&values), Some(2.5));
        assert_eq!(mean(&[]), None);
        let spread = sample_std_dev(&values).expect("four values");
        assert!((spread - 1.290_994_448_735_813).abs() < 1.0e-12);
        assert_eq!(sample_std_dev(&[1.0]), None);
    }

    #[test]
    fn a_perfect_line_fits_with_zero_residual() {
        let x = [-1.0, 0.0, 1.0, 2.0];
        let y = [1.0, 3.0, 5.0, 7.0];
        let fit = fit_line(&x, &y).expect("nondegenerate");
        assert!((fit.slope - 2.0).abs() < 1.0e-12);
        assert!((fit.intercept - 3.0).abs() < 1.0e-12);
        assert!(fit.slope_standard_error.abs() < 1.0e-12);
        assert!((fit.r_squared - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn a_degenerate_spread_returns_none() {
        assert!(fit_line(&[1.0, 1.0, 1.0], &[1.0, 2.0, 3.0]).is_none());
        assert!(fit_line(&[1.0, 2.0], &[1.0, 2.0]).is_none());
    }

    #[test]
    fn the_seeded_generator_repeats() {
        let mut first = Rng::new(0x1234);
        let mut second = Rng::new(0x1234);
        for _ in 0..8 {
            assert_eq!(first.next_f64(), second.next_f64());
        }
    }
}
