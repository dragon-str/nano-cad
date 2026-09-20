//! The loading-sector estimate: the time for one pocket to fill from a gas,
//! and the part of a rotor turn that this time costs.
//!
//! The model is a single binding site in a gas. The arrival rate is the gas
//! kinetic flux onto the pocket mouth, `n * v / 4`, multiplied by the mouth
//! area and a sticking probability. The site reaches the target occupancy when
//! that many molecules have arrived. The model states no diffusion in a liquid
//! and no crowding. Treat every result as an estimate. See ADR-0081.

use std::f64::consts::PI;

use nanocad_parts::guest;

/// The gas constant, in joules per mole kelvin.
pub const GAS_CONSTANT_J_PER_MOL_K: f64 = 8.314_462_618;

/// The Boltzmann constant, in joules per kelvin.
pub const BOLTZMANN_J_PER_K: f64 = 1.380_649e-23;

/// The gas pressures, temperature and geometry of a loading estimate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LoadingTarget {
    /// The oxygen partial pressure, in pascals.
    pub o2_pressure_pa: f64,
    /// The carbon dioxide partial pressure, in pascals.
    pub co2_pressure_pa: f64,
    /// The gas temperature, in kelvin.
    pub temperature_k: f64,
    /// The radius of the pocket mouth, in metres.
    pub opening_radius_m: f64,
    /// The fraction of the arriving molecules that stay on the site.
    pub sticking_probability: f64,
    /// The rotor angular rate, in radians per second.
    pub angular_rate_rad_per_s: f64,
    /// The site occupancy that counts as loaded, from zero to one.
    pub target_occupancy: f64,
}

impl Default for LoadingTarget {
    fn default() -> Self {
        Self {
            o2_pressure_pa: 1.33e4,
            co2_pressure_pa: 5.33e3,
            temperature_k: 310.15,
            opening_radius_m: 6.0e-10,
            sticking_probability: 1.0,
            angular_rate_rad_per_s: 1.0e6,
            target_occupancy: 0.63,
        }
    }
}

/// The loading estimate for one gas.
#[derive(Clone, Debug, PartialEq)]
pub struct GasLoadingReport {
    /// The gas id in the guest table.
    pub gas: String,
    /// The number density, in molecules per cubic metre.
    pub number_density_per_m3: f64,
    /// The mean molecular speed, in metres per second.
    pub mean_speed_m_per_s: f64,
    /// The capture rate of one site, in cubic metres per second.
    pub capture_rate_m3_per_s: f64,
    /// The time to reach the target occupancy, in seconds.
    pub saturation_time_s: f64,
    /// The rotor angle that the saturation time spans, in radians.
    pub loading_angle_rad: f64,
}

/// The loading estimate for the two gases.
#[derive(Clone, Debug, PartialEq)]
pub struct LoadingReport {
    /// The mouth area, in square metres.
    pub opening_area_m2: f64,
    /// The oxygen estimate.
    pub o2: GasLoadingReport,
    /// The carbon dioxide estimate.
    pub co2: GasLoadingReport,
    /// The wider of the two loading angles, in radians.
    pub loading_angle_rad: f64,
    /// The rotor angular rate that the angles use, in radians per second.
    pub angular_rate_rad_per_s: f64,
}

/// Estimates the loading of one pocket from the oxygen and carbon dioxide
/// partial pressures.
pub fn loading_report(target: &LoadingTarget) -> LoadingReport {
    let area = PI * target.opening_radius_m * target.opening_radius_m;
    let gas = |id: &str, pressure_pa: f64| -> GasLoadingReport {
        let molar_mass = guest(id).map_or(f64::NAN, |g| g.molar_mass_kg_per_mol);
        let number_density = pressure_pa / (BOLTZMANN_J_PER_K * target.temperature_k);
        let mean_speed =
            (8.0 * GAS_CONSTANT_J_PER_MOL_K * target.temperature_k / (PI * molar_mass)).sqrt();
        let capture_rate = target.sticking_probability * area * mean_speed / 4.0;
        let saturation_time = target.target_occupancy / (capture_rate * number_density);
        GasLoadingReport {
            gas: id.to_string(),
            number_density_per_m3: number_density,
            mean_speed_m_per_s: mean_speed,
            capture_rate_m3_per_s: capture_rate,
            saturation_time_s: saturation_time,
            loading_angle_rad: target.angular_rate_rad_per_s * saturation_time,
        }
    };
    let o2 = gas("o2", target.o2_pressure_pa);
    let co2 = gas("co2", target.co2_pressure_pa);
    let loading_angle_rad = o2.loading_angle_rad.max(co2.loading_angle_rad);
    LoadingReport {
        opening_area_m2: area,
        o2,
        co2,
        loading_angle_rad,
        angular_rate_rad_per_s: target.angular_rate_rad_per_s,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_oxygen_saturation_time_is_short() {
        let report = loading_report(&LoadingTarget::default());
        assert!(report.o2.saturation_time_s > 1.0e-10);
        assert!(report.o2.saturation_time_s < 1.0e-7);
    }

    #[test]
    fn the_loading_angle_is_a_small_part_of_a_turn() {
        let report = loading_report(&LoadingTarget::default());
        assert!(report.loading_angle_rad > 0.0);
        assert!(report.loading_angle_rad < PI / 2.0);
    }

    #[test]
    fn the_lower_partial_pressure_fills_more_slowly() {
        let report = loading_report(&LoadingTarget::default());
        assert!(report.co2.saturation_time_s > report.o2.saturation_time_s);
    }

    #[test]
    fn a_faster_rotor_needs_a_wider_angle() {
        let slow = loading_report(&LoadingTarget {
            angular_rate_rad_per_s: 1.0e6,
            ..LoadingTarget::default()
        });
        let fast = loading_report(&LoadingTarget {
            angular_rate_rad_per_s: 2.0e6,
            ..LoadingTarget::default()
        });
        assert!((fast.loading_angle_rad - 2.0 * slow.loading_angle_rad).abs() < 1.0e-12);
    }

    #[test]
    fn the_mouth_area_sets_the_capture_rate() {
        let small = loading_report(&LoadingTarget::default());
        let large = loading_report(&LoadingTarget {
            opening_radius_m: 1.2e-9,
            ..LoadingTarget::default()
        });
        let factor = (1.2e-9 / 6.0e-10_f64).powi(2);
        assert!(
            (large.o2.capture_rate_m3_per_s / small.o2.capture_rate_m3_per_s - factor).abs()
                < 1.0e-9
        );
    }
}
