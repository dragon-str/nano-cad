//! Competitive occupancy of one binding site.
//!
//! One site holds at most one guest. Each candidate guest reaches the site
//! from its own partial pressure. For a site in liquid water, the partial
//! pressure of water is its vapour pressure: a molecule in the liquid and a
//! molecule in the saturated vapour have the same chemical potential, so the
//! gas-phase reference is exact for the exchange equilibrium.
//!
//! The binding free energy of guest `i` against one bar is
//! `dG_i = dH_i + T * S_i`, where `dH_i` is the binding enthalpy (negative
//! when the guest binds) and `S_i` is the standard gas entropy that the guest
//! loses when it binds. The model assumes the guest keeps no entropy in the
//! site, so it states an upper bound on the entropy cost.
//!
//! The occupancy is the multi-component Langmuir form
//! `theta_i = x_i / (1 + sum_j x_j)` with `x_i = (p_i / P0) exp(-dG_i / kT)`.
//!
//! The model states no kinetics, no crowding and no cooperativity. Treat
//! every result as an estimate. See ADR-0086.

use crate::delivery::STANDARD_PRESSURE_PA;
use crate::loading::BOLTZMANN_J_PER_K;

/// The Avogadro constant, per mole.
pub const AVOGADRO_PER_MOL: f64 = 6.022_140_76e23;

/// One guest that competes for the site.
#[derive(Debug, Clone, PartialEq)]
pub struct Competitor {
    /// The guest id, for example `o2` or `water`.
    pub id: String,
    /// The partial pressure of the guest at the site, in pascals.
    pub pressure_pa: f64,
    /// The binding enthalpy of one guest, in joules. Negative binds.
    pub binding_enthalpy_j: f64,
    /// The standard gas entropy that one guest loses on binding, in J/K.
    pub entropy_loss_j_per_k: f64,
}

impl Competitor {
    /// A competitor from molar quantities: kJ/mol and J/(mol K).
    pub fn from_molar(
        id: &str,
        pressure_pa: f64,
        binding_enthalpy_kj_per_mol: f64,
        entropy_loss_j_per_mol_k: f64,
    ) -> Self {
        Self {
            id: id.to_string(),
            pressure_pa,
            binding_enthalpy_j: binding_enthalpy_kj_per_mol * 1.0e3 / AVOGADRO_PER_MOL,
            entropy_loss_j_per_k: entropy_loss_j_per_mol_k / AVOGADRO_PER_MOL,
        }
    }

    /// The binding free energy against one bar at `temperature_k`, in joules.
    pub fn binding_free_energy_j(&self, temperature_k: f64) -> f64 {
        self.binding_enthalpy_j + temperature_k * self.entropy_loss_j_per_k
    }

    /// The Boltzmann weight of the bound state against the empty site.
    #[cfg(test)]
    fn weight(&self, temperature_k: f64) -> f64 {
        let kt = BOLTZMANN_J_PER_K * temperature_k;
        (self.pressure_pa / STANDARD_PRESSURE_PA)
            * (-self.binding_free_energy_j(temperature_k) / kt).exp()
    }
}

/// The occupancy of one competitor.
#[derive(Debug, Clone, PartialEq)]
pub struct CompetitorOccupancy {
    /// The guest id.
    pub id: String,
    /// The binding free energy against one bar, in joules.
    pub binding_free_energy_j: f64,
    /// The fraction of time the guest holds the site.
    pub occupancy: f64,
}

/// The competitive occupancy of one site.
#[derive(Debug, Clone, PartialEq)]
pub struct CompetitionReport {
    /// The temperature, in kelvin.
    pub temperature_k: f64,
    /// The fraction of time the site is empty.
    pub empty: f64,
    /// One row for each competitor, in input order.
    pub guests: Vec<CompetitorOccupancy>,
}

impl CompetitionReport {
    /// The occupancy ratio of guest `a` over guest `b`.
    ///
    /// Returns `None` when an id is absent or when `b` has zero occupancy.
    pub fn ratio(&self, a: &str, b: &str) -> Option<f64> {
        let find = |id: &str| self.guests.iter().find(|g| g.id == id);
        let (a, b) = (find(a)?, find(b)?);
        (b.occupancy > 0.0).then(|| a.occupancy / b.occupancy)
    }
}

/// The competitive occupancy of one site at `temperature_k`.
///
/// The weights are summed in log space, so a very strong binder does not
/// overflow.
pub fn competition_report(competitors: &[Competitor], temperature_k: f64) -> CompetitionReport {
    let kt = BOLTZMANN_J_PER_K * temperature_k;
    let logs: Vec<f64> = competitors
        .iter()
        .map(|c| {
            (c.pressure_pa / STANDARD_PRESSURE_PA).ln()
                - c.binding_free_energy_j(temperature_k) / kt
        })
        .collect();
    // log of the partition sum 1 + sum exp(l_i), with the empty state as l=0.
    let top = logs.iter().copied().fold(0.0_f64, f64::max);
    let sum = (-top).exp() + logs.iter().map(|l| (l - top).exp()).sum::<f64>();
    let log_z = top + sum.ln();
    let guests = competitors
        .iter()
        .zip(&logs)
        .map(|(c, l)| CompetitorOccupancy {
            id: c.id.clone(),
            binding_free_energy_j: c.binding_free_energy_j(temperature_k),
            occupancy: (l - log_z).exp(),
        })
        .collect();
    CompetitionReport {
        temperature_k,
        empty: (-log_z).exp(),
        guests,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const BODY_K: f64 = 310.15;

    #[test]
    fn one_guest_matches_the_langmuir_form() {
        let c = Competitor::from_molar("o2", 13_300.0, -70.0, 205.2);
        let r = competition_report(std::slice::from_ref(&c), BODY_K);
        let x = c.weight(BODY_K);
        assert!((r.guests[0].occupancy - x / (1.0 + x)).abs() < 1e-12);
        assert!((r.empty + r.guests[0].occupancy - 1.0).abs() < 1e-12);
    }

    #[test]
    fn the_occupancies_sum_to_one() {
        let list = [
            Competitor::from_molar("o2", 13_300.0, -96.0, 205.2),
            Competitor::from_molar("water", 6_280.0, -10.0, 188.8),
            Competitor::from_molar("co2", 5_330.0, 0.0, 213.8),
        ];
        let r = competition_report(&list, BODY_K);
        let total: f64 = r.empty + r.guests.iter().map(|g| g.occupancy).sum::<f64>();
        assert!((total - 1.0).abs() < 1e-12);
    }

    #[test]
    fn a_very_strong_binder_does_not_overflow() {
        let list = [
            Competitor::from_molar("strong", 1.0e5, -2000.0, 0.0),
            Competitor::from_molar("weak", 1.0e5, 0.0, 0.0),
        ];
        let r = competition_report(&list, BODY_K);
        assert!(r.guests[0].occupancy.is_finite());
        assert!((r.guests[0].occupancy - 1.0).abs() < 1e-12);
        assert!(r.guests[1].occupancy >= 0.0);
    }

    #[test]
    fn the_ratio_follows_the_free_energy_gap() {
        // Equal pressures and a gap of one kT give a ratio of e.
        let kt_kj = BOLTZMANN_J_PER_K * BODY_K * AVOGADRO_PER_MOL / 1.0e3;
        let list = [
            Competitor::from_molar("a", 1.0e3, -kt_kj, 0.0),
            Competitor::from_molar("b", 1.0e3, 0.0, 0.0),
        ];
        let r = competition_report(&list, BODY_K);
        let ratio = r.ratio("a", "b").unwrap();
        assert!((ratio - std::f64::consts::E).abs() < 1e-9);
        assert!(r.ratio("a", "missing").is_none());
    }

    #[test]
    fn the_metal_site_prefers_oxygen_over_water() {
        // ADR-0086 kill test: QM enthalpies on the Fe(NH3)4 proxy.
        let list = [
            Competitor::from_molar("o2", 13_300.0, -95.6, 205.2),
            Competitor::from_molar("water", 6_280.0, -10.3, 188.8),
        ];
        let r = competition_report(&list, BODY_K);
        assert!(r.ratio("o2", "water").unwrap() > 1.0e8);
    }
}
