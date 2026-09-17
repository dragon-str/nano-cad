//! The planetary gear search: a CMA-ES over the tooth geometry.

use nanocad_meter::{Clearance, MovingAtoms, SlipBarrier, SlipBarrierTarget};
use nanocad_parts::{ParameterSet, PartGenerator, PlanetaryGenerator};

use crate::cmaes::{Objective, ParameterBounds};

/// The sun rate that the viewer and the scene use, in radians per second.
pub const SUN_RAD_PER_S: f64 = -0.34;

/// The temperature of the thermal energy, in kelvin.
pub const TEMPERATURE_K: f64 = 300.0;

/// The clearance shortfall that costs one penalty unit, in metres.
const CLEARANCE_SCALE_M: f64 = 1.0e-10;

/// The cost of one clearance shortfall of `CLEARANCE_SCALE_M`, in `k T`.
const CLEARANCE_PENALTY_KT: f64 = 100.0;

/// The cost of one gear design.
#[derive(Clone, Debug, PartialEq)]
pub struct PlanetaryScore {
    /// The slip-barrier energy over one sun tooth pitch, in joules.
    pub barrier_j: f64,
    /// The barrier in units of the thermal energy `k T`.
    pub barrier_over_kt: f64,
    /// The least atom-to-atom distance over one revolution, in metres.
    pub clearance_m: f64,
    /// The number of atoms in the set.
    pub atom_count: usize,
    /// The parameter vector that produced the score.
    pub parameters: Vec<f64>,
    /// True when the clearance meets its target.
    pub feasible: bool,
    /// The clearance shortfall below the target, in metres. It is zero when
    /// the design is feasible.
    pub clearance_shortfall_m: f64,
}

/// A CMA-ES search over the tooth geometry of the planetary set.
///
/// The module and the layer count stay at their defaults, so every candidate
/// has a comparable size. The objective is the slip barrier in units of the
/// thermal energy, plus a penalty for a clearance shortfall, so the search can
/// move through an infeasible region toward a feasible design. A parameter set
/// that the generator refuses is rejected.
#[derive(Clone, Debug, PartialEq)]
pub struct PlanetarySearch {
    /// The searched parameter names, in vector order.
    pub names: Vec<&'static str>,
    /// The range of each searched parameter.
    pub bounds: Vec<ParameterBounds>,
    /// The number of time samples in the slip sweep.
    pub slip_steps: usize,
    /// The temperature of the thermal energy, in kelvin.
    pub temperature_k: f64,
}

impl Default for PlanetarySearch {
    fn default() -> Self {
        Self {
            names: vec!["sun_teeth", "planet_teeth", "planet_count"],
            bounds: vec![
                ParameterBounds::new(9.0, 16.0, true),
                ParameterBounds::new(7.0, 12.0, true),
                ParameterBounds::new(2.0, 4.0, true),
            ],
            slip_steps: 60,
            temperature_k: TEMPERATURE_K,
        }
    }
}

impl PlanetarySearch {
    /// Builds a search with the default geometry ranges.
    pub fn new() -> Self {
        Self::default()
    }

    /// Returns the starting vector, at every generator default.
    ///
    /// The default design is feasible, so the search starts from a valid
    /// point. A name that has no spec falls back to the middle of its range.
    pub fn initial(&self) -> Vec<f64> {
        let specs = PlanetaryGenerator.parameters();
        self.names
            .iter()
            .zip(self.bounds.iter())
            .map(|(name, bounds)| {
                let default = specs
                    .iter()
                    .find(|spec| spec.name == *name)
                    .map(|spec| spec.default)
                    .unwrap_or(0.5 * (bounds.min + bounds.max));
                bounds.clamp(default)
            })
            .collect()
    }

    /// Turns a parameter vector into a generator parameter set.
    pub fn parameters_of(&self, vector: &[f64]) -> ParameterSet {
        let mut parameters = ParameterSet::new();
        for (name, value) in self.names.iter().zip(vector.iter()) {
            parameters.set(*name, *value);
        }
        parameters
    }

    /// Scores one design. Returns `None` when it is invalid or infeasible.
    pub fn score(&self, vector: &[f64]) -> Option<PlanetaryScore> {
        let set = PlanetaryGenerator.build(&self.parameters_of(vector)).ok()?;
        let moving = MovingAtoms::planetary(&set, &set.design, SUN_RAD_PER_S);
        let clearance = Clearance::default();
        let report = clearance.measure(&moving);
        let feasible = report.passes(&clearance.target);
        let clearance_shortfall_m = if feasible {
            0.0
        } else {
            (clearance.target.minimum_m - report.minimum_m).max(0.0)
        };
        let slip = SlipBarrier::new(SlipBarrierTarget {
            steps: self.slip_steps,
            temperature_k: self.temperature_k,
            ..Default::default()
        });
        let barrier = slip.measure(&moving, set.design.sun_teeth());
        Some(PlanetaryScore {
            barrier_j: barrier.barrier_j,
            barrier_over_kt: barrier.barrier_over_kt(self.temperature_k),
            clearance_m: report.minimum_m,
            atom_count: set.part.topology.atom_count(),
            parameters: vector.to_vec(),
            feasible,
            clearance_shortfall_m,
        })
    }
}

impl Objective for PlanetarySearch {
    fn evaluate(&self, parameters: &[f64]) -> Option<f64> {
        self.score(parameters).map(|score| {
            let penalty = score.clearance_shortfall_m / CLEARANCE_SCALE_M * CLEARANCE_PENALTY_KT;
            score.barrier_over_kt + penalty
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cmaes::{minimize, CmaEsOptions};

    #[test]
    fn the_starting_vector_is_inside_every_range() {
        let search = PlanetarySearch::default();
        let initial = search.initial();
        assert_eq!(initial.len(), search.bounds.len());
        for (value, bounds) in initial.iter().zip(search.bounds.iter()) {
            assert!(*value >= bounds.min && *value <= bounds.max);
            if bounds.integer {
                assert_eq!(value.fract(), 0.0);
            }
        }
    }

    #[test]
    fn a_vector_becomes_a_parameter_set() {
        let search = PlanetarySearch::default();
        let parameters = search.parameters_of(&[12.0, 9.0, 3.0]);
        assert_eq!(parameters.get("sun_teeth"), Some(12.0));
        assert_eq!(parameters.get("planet_teeth"), Some(9.0));
        assert_eq!(parameters.get("planet_count"), Some(3.0));
        assert_eq!(parameters.len(), 3);
    }

    #[test]
    fn an_impossible_design_is_rejected() {
        let search = PlanetarySearch::default();
        assert!(search.score(&[0.0, 0.0, 0.0]).is_none());
        assert!(search.score(&[-1.0, -1.0, -1.0]).is_none());
    }

    /// This test generates real gear sets, so it runs in release only.
    /// `cargo run --release -p nanocad-opt --example optimize_planetary` is
    /// the gate.
    #[test]
    #[ignore = "generates a full gear set; run the example in release"]
    fn a_short_search_runs_and_respects_the_ranges() {
        let search = PlanetarySearch::default();
        let initial = search.initial();
        let options = CmaEsOptions {
            population: 4,
            generations: 2,
            seed: 7,
            ..Default::default()
        };
        let report = minimize(&initial, &search.bounds, &options, &search);
        assert_eq!(report.evaluations, 8);
        for (value, bounds) in report.best_parameters.iter().zip(search.bounds.iter()) {
            assert!(*value >= bounds.min && *value <= bounds.max);
        }
    }
}
