//! Parameter search for generated parts.
//!
//! The optimize stage searches the generator parameters for a better score.
//! The search is a deterministic CMA-ES over a small parameter vector. Every
//! candidate goes through the real generator and the real metrics, so a
//! result is a measured score, not a model output.
#![forbid(unsafe_code)]

pub mod cmaes;
pub mod planetary;

pub use cmaes::{minimize, CmaEsOptions, CmaEsReport, Objective, ParameterBounds};
pub use planetary::{PlanetaryScore, PlanetarySearch, SUN_RAD_PER_S, TEMPERATURE_K};
