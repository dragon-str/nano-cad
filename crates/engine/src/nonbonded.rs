use crate::error::EngineError;

/// Builds a flat list of every unordered atom pair.
///
/// Entry `2k` is smaller than entry `2k + 1`. A pair list of this form is a
/// valid input to the pair-list path of the non-bonded terms.
pub(crate) fn all_pairs(atom_count: usize) -> Vec<u32> {
    let mut pairs = Vec::with_capacity(atom_count.saturating_mul(atom_count.saturating_sub(1)));
    for i in 0..atom_count {
        for j in (i + 1)..atom_count {
            pairs.push(i as u32);
            pairs.push(j as u32);
        }
    }
    pairs
}

/// An orthorhombic periodic box with a length per axis.
///
/// A length of zero means that the axis is not periodic. A pair displacement is
/// reduced by the minimum-image convention on each periodic axis.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PeriodicBox {
    lengths_m: [f64; 3],
}

impl Default for PeriodicBox {
    fn default() -> Self {
        Self::non_periodic()
    }
}

impl PeriodicBox {
    /// Builds a box that is not periodic on any axis.
    pub const fn non_periodic() -> Self {
        Self {
            lengths_m: [0.0; 3],
        }
    }

    /// Builds a box from three axis lengths in SI metres.
    ///
    /// Every length must be finite and non-negative. A zero length disables
    /// periodicity on that axis.
    pub fn new(lengths_m: [f64; 3]) -> Result<Self, EngineError> {
        for (index, length_m) in lengths_m.iter().enumerate() {
            if !length_m.is_finite() || *length_m < 0.0 {
                return Err(EngineError::InvalidPeriodicBox {
                    index,
                    length_m: *length_m,
                });
            }
        }
        Ok(Self { lengths_m })
    }

    /// Returns the three axis lengths in SI metres.
    pub const fn lengths_m(&self) -> [f64; 3] {
        self.lengths_m
    }

    /// Reports whether any axis is periodic.
    pub fn is_periodic(&self) -> bool {
        self.lengths_m.iter().any(|length_m| *length_m > 0.0)
    }

    /// Reduces a pair displacement by the minimum-image convention.
    ///
    /// A non-periodic axis is unchanged.
    pub fn minimum_image(&self, delta_m: [f64; 3]) -> [f64; 3] {
        let mut out_m = delta_m;
        for (axis, value_m) in out_m.iter_mut().enumerate() {
            let length_m = self.lengths_m[axis];
            if length_m > 0.0 {
                *value_m -= length_m * (*value_m / length_m).round();
            }
        }
        out_m
    }
}

/// A cutoff with a smooth switching region for a non-bonded potential.
///
/// The potential is multiplied by a polynomial switching function `S(r)`. The
/// function is one below `switch_on_m`, falls smoothly to zero at `cutoff_m`,
/// and stays zero beyond the cutoff.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cutoff {
    cutoff_m: f64,
    switch_on_m: f64,
}

impl Cutoff {
    /// Builds a cutoff from a cutoff distance and a switch-on distance.
    ///
    /// The cutoff must be positive and finite. The switch-on distance must be
    /// finite, non-negative, and below the cutoff.
    pub fn new(cutoff_m: f64, switch_on_m: f64) -> Result<Self, EngineError> {
        if !cutoff_m.is_finite() || cutoff_m <= 0.0 {
            return Err(EngineError::NonPositiveCutoff { cutoff_m });
        }
        if !switch_on_m.is_finite() || switch_on_m < 0.0 || switch_on_m >= cutoff_m {
            return Err(EngineError::InvalidSwitchingRange {
                cutoff_m,
                switch_on_m,
            });
        }
        Ok(Self {
            cutoff_m,
            switch_on_m,
        })
    }

    /// Returns the cutoff distance in SI metres.
    pub const fn cutoff_m(&self) -> f64 {
        self.cutoff_m
    }

    /// Returns the switch-on distance in SI metres.
    pub const fn switch_on_m(&self) -> f64 {
        self.switch_on_m
    }

    /// Returns the switching value `S(r)` and its derivative `dS/dr`.
    ///
    /// The polynomial form is continuous with a continuous first derivative at
    /// both ends of the switching region.
    pub fn switch_value(&self, r_m: f64) -> (f64, f64) {
        if r_m <= self.switch_on_m {
            return (1.0, 0.0);
        }
        if r_m >= self.cutoff_m {
            return (0.0, 0.0);
        }
        let r_on = self.switch_on_m;
        let r_cut = self.cutoff_m;
        let denom = r_cut * r_cut - r_on * r_on;
        let denom_cubed = denom * denom * denom;
        let d = r_cut * r_cut - r_m * r_m;
        let e = r_cut * r_cut + 2.0 * r_m * r_m - 3.0 * r_on * r_on;
        let value = d * d * e / denom_cubed;
        let derivative = 12.0 * r_m * d * (r_on * r_on - r_m * r_m) / denom_cubed;
        (value, derivative)
    }

    /// Checks that the cutoff does not exceed half of any periodic box length.
    ///
    /// The minimum-image convention is only valid when the cutoff is at most
    /// half the shortest periodic length.
    pub fn validate_box(&self, periodic_box: &PeriodicBox) -> Result<(), EngineError> {
        for length_m in periodic_box.lengths_m() {
            if length_m > 0.0 && self.cutoff_m > 0.5 * length_m {
                return Err(EngineError::CutoffExceedsHalfBox {
                    cutoff_m: self.cutoff_m,
                    length_m,
                });
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_negative_periodic_length_is_rejected() {
        assert!(matches!(
            PeriodicBox::new([1.0, -1.0, 0.0]),
            Err(EngineError::InvalidPeriodicBox { .. })
        ));
    }

    #[test]
    fn an_invalid_switching_range_is_rejected() {
        assert!(matches!(
            Cutoff::new(1.0e-9, 1.0e-9),
            Err(EngineError::InvalidSwitchingRange { .. })
        ));
        assert!(matches!(
            Cutoff::new(1.0e-9, 2.0e-9),
            Err(EngineError::InvalidSwitchingRange { .. })
        ));
        assert!(matches!(
            Cutoff::new(1.0e-9, -1.0),
            Err(EngineError::InvalidSwitchingRange { .. })
        ));
    }

    #[test]
    fn the_minimum_image_wraps_to_the_near_image() {
        let periodic_box = PeriodicBox::new([1.0e-9, 0.0, 0.0]).expect("valid box");
        let delta_m = periodic_box.minimum_image([0.9e-9, 0.0, 0.0]);
        assert!((delta_m[0] + 0.1e-9).abs() < 1.0e-24);
        assert_eq!(delta_m[1], 0.0);
    }

    #[test]
    fn the_switching_function_is_one_in_the_core_and_zero_at_the_cutoff() {
        let cutoff = Cutoff::new(1.0e-9, 0.8e-9).expect("valid cutoff");
        assert_eq!(cutoff.switch_value(0.5e-9), (1.0, 0.0));
        assert_eq!(cutoff.switch_value(1.0e-9), (0.0, 0.0));
        let (value, derivative) = cutoff.switch_value(0.9e-9);
        assert!(value > 0.0 && value < 1.0);
        assert!(derivative != 0.0);
    }

    #[test]
    fn a_cutoff_longer_than_half_the_box_is_rejected() {
        let cutoff = Cutoff::new(0.6e-9, 0.5e-9).expect("valid cutoff");
        let periodic_box = PeriodicBox::new([1.0e-9, 0.0, 0.0]).expect("valid box");
        assert!(matches!(
            cutoff.validate_box(&periodic_box),
            Err(EngineError::CutoffExceedsHalfBox { .. })
        ));
    }
}
