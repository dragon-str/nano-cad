//! The delivery model: how much of a bound guest a passive site releases over a
//! stated concentration swing.
//!
//! A site that binds a gas follows the Hill curve
//! `theta(P) = P^n / (P50^n + P^n)`. The delivered fraction of one cycle is
//! the occupancy at the loading pressure minus the occupancy at the unloading
//! pressure. The loading pressure is high (the lung) and the unloading pressure
//! is low (the tissue) for oxygen, and the reverse for carbon dioxide.
//!
//! The model states no kinetics, no hysteresis and no competing ligand. An
//! occupancy is an equilibrium value. Treat every result as an estimate. See
//! ADR-0082.

use crate::loading::BOLTZMANN_J_PER_K;

/// The standard pressure, in pascals. One bar.
pub const STANDARD_PRESSURE_PA: f64 = 100_000.0;

/// The pressure that a binding energy gives a half-full site, in pascals.
///
/// The binding energy is the free energy of one guest at the standard pressure,
/// relative to the bound state. A more negative energy gives a lower `P50` and
/// a stronger binding.
pub fn p50_from_binding_energy(delta_g_j: f64, temperature_k: f64) -> f64 {
    STANDARD_PRESSURE_PA * (delta_g_j / (BOLTZMANN_J_PER_K * temperature_k)).exp()
}

/// The binding energy that gives a site the stated `P50`, in joules.
pub fn binding_energy_from_p50(p50_pa: f64, temperature_k: f64) -> f64 {
    BOLTZMANN_J_PER_K * temperature_k * (p50_pa / STANDARD_PRESSURE_PA).ln()
}

/// The `P50` that maximises the delivered fraction over a swing.
///
/// The optimum is the geometric mean of the two pressures, for every Hill
/// coefficient. The result is in pascals.
pub fn optimal_p50_pa(load_pressure_pa: f64, unload_pressure_pa: f64) -> f64 {
    (load_pressure_pa * unload_pressure_pa).sqrt()
}

/// The site occupancy at one pressure, from zero to one.
pub fn occupancy(pressure_pa: f64, p50_pa: f64, hill_coefficient: f64) -> f64 {
    let load = pressure_pa.powf(hill_coefficient);
    let half = p50_pa.powf(hill_coefficient);
    load / (half + load)
}

/// The pressures, the temperature and the site parameters of one delivery.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DeliveryTarget {
    /// The loading partial pressure, in pascals.
    pub load_pressure_pa: f64,
    /// The unloading partial pressure, in pascals.
    pub unload_pressure_pa: f64,
    /// The body temperature, in kelvin.
    pub temperature_k: f64,
    /// The half-fill pressure of the site, in pascals.
    pub p50_pa: f64,
    /// The Hill coefficient. One is a non-cooperative site.
    pub hill_coefficient: f64,
    /// The number of guest molecules that one pocket holds when it is full.
    pub pocket_capacity: f64,
}

impl Default for DeliveryTarget {
    /// A non-cooperative oxygen site tuned to the lung and the tissue.
    fn default() -> Self {
        let load_pressure_pa = 13.3e3;
        let unload_pressure_pa = 5.3e3;
        Self {
            load_pressure_pa,
            unload_pressure_pa,
            temperature_k: 310.15,
            p50_pa: optimal_p50_pa(load_pressure_pa, unload_pressure_pa),
            hill_coefficient: 1.0,
            pocket_capacity: 10.0,
        }
    }
}

/// The delivery result for one site and one swing.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DeliveryReport {
    /// The loading partial pressure, in pascals.
    pub load_pressure_pa: f64,
    /// The unloading partial pressure, in pascals.
    pub unload_pressure_pa: f64,
    /// The half-fill pressure, in pascals.
    pub p50_pa: f64,
    /// The binding energy of the site, in joules.
    pub binding_energy_j: f64,
    /// The site occupancy at the loading pressure.
    pub occupancy_loaded: f64,
    /// The site occupancy at the unloading pressure.
    pub occupancy_unloaded: f64,
    /// The fraction of the pocket capacity that leaves in one cycle.
    pub delivered_fraction: f64,
    /// The molecules that leave one pocket in one cycle.
    pub delivered_molecules_per_cycle: f64,
}

impl DeliveryReport {
    /// The pressure ratio of the swing, loading over unloading.
    pub fn pressure_ratio(&self) -> f64 {
        self.load_pressure_pa / self.unload_pressure_pa
    }
}

/// Measures the delivered fraction of one site over one swing.
pub fn delivery_report(target: &DeliveryTarget) -> DeliveryReport {
    let occupancy_loaded = occupancy(
        target.load_pressure_pa,
        target.p50_pa,
        target.hill_coefficient,
    );
    let occupancy_unloaded = occupancy(
        target.unload_pressure_pa,
        target.p50_pa,
        target.hill_coefficient,
    );
    let delivered_fraction = occupancy_loaded - occupancy_unloaded;
    DeliveryReport {
        load_pressure_pa: target.load_pressure_pa,
        unload_pressure_pa: target.unload_pressure_pa,
        p50_pa: target.p50_pa,
        binding_energy_j: binding_energy_from_p50(target.p50_pa, target.temperature_k),
        occupancy_loaded,
        occupancy_unloaded,
        delivered_fraction,
        delivered_molecules_per_cycle: delivered_fraction * target.pocket_capacity,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The oxygen swing of the body, in pascals.
    const LUNG_O2_PA: f64 = 13.3e3;
    const TISSUE_O2_PA: f64 = 5.3e3;

    #[test]
    fn the_optimum_is_the_geometric_mean() {
        let p50 = optimal_p50_pa(LUNG_O2_PA, TISSUE_O2_PA);
        assert!((p50 - (LUNG_O2_PA * TISSUE_O2_PA).sqrt()).abs() < 1.0e-9);
        let moved = optimal_p50_pa(LUNG_O2_PA, TISSUE_O2_PA);
        assert!((p50 - moved).abs() < 1.0e-12);
    }

    #[test]
    fn a_non_cooperative_oxygen_site_delivers_the_langmuir_bound() {
        let report = delivery_report(&DeliveryTarget::default());
        assert!(
            (report.delivered_fraction - 0.2260).abs() < 1.0e-3,
            "the delivered fraction is {}",
            report.delivered_fraction
        );
        assert!(report.occupancy_loaded < 0.65);
        assert!(report.occupancy_unloaded > 0.35);
        assert!((report.pressure_ratio() - 2.509).abs() < 1.0e-2);
    }

    #[test]
    fn a_hill_site_delivers_more_than_a_flat_site() {
        let flat = delivery_report(&DeliveryTarget {
            load_pressure_pa: LUNG_O2_PA,
            unload_pressure_pa: 3.0e3,
            ..DeliveryTarget::default()
        });
        let cooperative = delivery_report(&DeliveryTarget {
            load_pressure_pa: LUNG_O2_PA,
            unload_pressure_pa: 3.0e3,
            hill_coefficient: 2.8,
            p50_pa: optimal_p50_pa(LUNG_O2_PA, 3.0e3),
            ..DeliveryTarget::default()
        });
        assert!(
            cooperative.delivered_fraction > flat.delivered_fraction,
            "cooperativity must sharpen the swing: {} against {}",
            cooperative.delivered_fraction,
            flat.delivered_fraction
        );
    }

    #[test]
    fn the_carbon_dioxide_swing_delivers_almost_nothing() {
        let oxygen = delivery_report(&DeliveryTarget {
            load_pressure_pa: LUNG_O2_PA,
            unload_pressure_pa: TISSUE_O2_PA,
            ..DeliveryTarget::default()
        });
        let carbon_dioxide = delivery_report(&DeliveryTarget {
            load_pressure_pa: 6.1e3,
            unload_pressure_pa: 5.3e3,
            ..DeliveryTarget::default()
        });
        assert!(
            carbon_dioxide.delivered_fraction < 0.05,
            "the carbon dioxide swing is small: {}",
            carbon_dioxide.delivered_fraction
        );
        assert!(oxygen.delivered_fraction > 5.0 * carbon_dioxide.delivered_fraction);
    }

    #[test]
    fn a_stronger_bond_lowers_the_half_fill_pressure() {
        let weak = p50_from_binding_energy(-5.0e-21, 310.15);
        let strong = p50_from_binding_energy(-1.5e-20, 310.15);
        assert!(
            strong < weak,
            "a stronger bond must bind at a lower pressure"
        );
    }

    #[test]
    fn the_binding_energy_and_the_half_fill_pressure_round_trip() {
        let p50 = 8.0e3;
        let energy = binding_energy_from_p50(p50, 310.15);
        let back = p50_from_binding_energy(energy, 310.15);
        assert!((back - p50).abs() / p50 < 1.0e-12);
        assert!(energy < 0.0, "binding at the body pressures is favorable");
    }

    #[test]
    fn no_swing_delivers_nothing() {
        let report = delivery_report(&DeliveryTarget {
            load_pressure_pa: 5.3e3,
            unload_pressure_pa: 5.3e3,
            ..DeliveryTarget::default()
        });
        assert!(report.delivered_fraction.abs() < 1.0e-15);
    }

    #[test]
    fn the_pocket_capacity_scales_the_delivered_molecules() {
        let small = delivery_report(&DeliveryTarget {
            pocket_capacity: 10.0,
            ..DeliveryTarget::default()
        });
        let large = delivery_report(&DeliveryTarget {
            pocket_capacity: 20.0,
            ..DeliveryTarget::default()
        });
        assert!(
            (large.delivered_molecules_per_cycle - 2.0 * small.delivered_molecules_per_cycle).abs()
                < 1.0e-12
        );
    }
}
