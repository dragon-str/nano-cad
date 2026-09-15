//! A small deterministic random generator for the thermostats.
//!
//! The generator is a SplitMix64 core with a Box-Muller normal transform. It
//! has no dependency and no global state. The same seed gives the same stream
//! on every platform.

/// A seeded SplitMix64 generator with a normal-distribution transform.
#[derive(Clone, Debug, PartialEq)]
pub struct Rng {
    state: u64,
    spare_normal: Option<f64>,
}

impl Rng {
    /// Builds a generator from a seed.
    pub const fn new(seed: u64) -> Self {
        Self {
            state: seed,
            spare_normal: None,
        }
    }

    /// Returns the next 64-bit value.
    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Returns the next value in `[0, 1)`.
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Returns the next standard normal value.
    ///
    /// The transform is Box-Muller. The second value is held for the next
    /// call, so the stream stays deterministic.
    pub fn next_normal(&mut self) -> f64 {
        if let Some(spare) = self.spare_normal.take() {
            return spare;
        }
        let mut u1 = self.next_f64();
        while u1 <= f64::MIN_POSITIVE {
            u1 = self.next_f64();
        }
        let u2 = self.next_f64();
        let magnitude = (-2.0 * u1.ln()).sqrt();
        let angle = std::f64::consts::TAU * u2;
        self.spare_normal = Some(magnitude * angle.sin());
        magnitude * angle.cos()
    }

    /// Returns the seed state.
    pub const fn state(&self) -> u64 {
        self.state
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_same_seed_gives_the_same_stream() {
        let mut a = Rng::new(0x1234_5678_9ABC_DEF0);
        let mut b = Rng::new(0x1234_5678_9ABC_DEF0);
        for _ in 0..1000 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn different_seeds_give_different_streams() {
        let mut a = Rng::new(1);
        let mut b = Rng::new(2);
        assert_ne!(a.next_u64(), b.next_u64());
    }

    #[test]
    fn the_uniform_stream_stays_in_range_with_a_sane_mean() {
        let mut rng = Rng::new(0xC0FF_EE00_1234_5678);
        let samples = 100_000;
        let mut sum = 0.0;
        for _ in 0..samples {
            let value = rng.next_f64();
            assert!((0.0..1.0).contains(&value));
            sum += value;
        }
        let mean = sum / samples as f64;
        assert!((mean - 0.5).abs() < 0.01, "mean {mean}");
    }

    #[test]
    fn the_normal_stream_has_zero_mean_and_unit_variance() {
        let mut rng = Rng::new(0xDEAD_BEEF_0BAD_F00D);
        let samples = 200_000;
        let mut sum = 0.0;
        let mut sum_sq = 0.0;
        for _ in 0..samples {
            let value = rng.next_normal();
            sum += value;
            sum_sq += value * value;
        }
        let mean = sum / samples as f64;
        let variance = sum_sq / samples as f64 - mean * mean;
        assert!(mean.abs() < 0.01, "mean {mean}");
        assert!((variance - 1.0).abs() < 0.02, "variance {variance}");
    }
}
