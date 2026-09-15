//! Planetary device assembly and URDF export.
//!
//! This exposes the M6-05 and M6-07 paths to the Python and agent surface: a
//! planetary gear set becomes an L2 rigid-body device, and the device exports
//! to URDF. The existing `assemble` tool only groups parts into a document, so
//! the device assembly needed its own binding.

use std::collections::HashMap;

use nanocad_jigs as jigs;
use nanocad_parts as parts;
use pyo3::prelude::*;

use super::util;

/// An assembled planetary gearbox, with the parameter library behind it.
#[pyclass(name = "PlanetaryAssembly", skip_from_py_object)]
#[derive(Clone)]
pub struct PlanetaryAssembly {
    inner: jigs::PlanetaryAssembly,
}

#[pymethods]
impl PlanetaryAssembly {
    /// The analytic sun-to-carrier gear ratio, with the ring fixed.
    #[getter]
    fn gear_ratio(&self) -> f64 {
        self.inner.analytic_gear_ratio()
    }

    #[getter]
    fn body_count(&self) -> usize {
        self.inner.system.body_count()
    }

    #[getter]
    fn planet_count(&self) -> usize {
        self.inner.design.planet_count()
    }

    #[getter]
    fn sun_teeth(&self) -> usize {
        self.inner.design.sun_teeth()
    }

    #[getter]
    fn planet_teeth(&self) -> usize {
        self.inner.design.planet_teeth()
    }

    #[getter]
    fn ring_teeth(&self) -> usize {
        self.inner.design.ring_teeth()
    }

    /// The sun joint angle, in radians.
    fn sun_angle_rad(&self) -> PyResult<f64> {
        self.inner.sun_angle_rad().map_err(util::err)
    }

    /// The carrier joint angle, in radians.
    fn carrier_angle_rad(&self) -> PyResult<f64> {
        self.inner.carrier_angle_rad().map_err(util::err)
    }

    /// The ring joint angle, in radians.
    fn ring_angle_rad(&self) -> PyResult<f64> {
        self.inner.ring_angle_rad().map_err(util::err)
    }

    /// The angle of one planet joint, in radians.
    fn planet_angle_rad(&self, planet: usize) -> PyResult<f64> {
        self.inner.planet_angle_rad(planet).map_err(util::err)
    }

    /// Sets every body to a rotation that satisfies the gear constraints.
    fn set_consistent_rotation(&mut self, sun_angular_velocity_rad_per_s: f64) -> PyResult<()> {
        self.inner
            .set_consistent_rotation(sun_angular_velocity_rad_per_s)
            .map_err(util::err)
    }

    /// Advances the device by one step with `iterations` constraint passes.
    fn step(&mut self, dt_s: f64, iterations: usize) -> PyResult<()> {
        self.inner.system.step(dt_s, iterations).map_err(util::err)
    }

    /// The largest gear-constraint error, in metres.
    fn max_gear_constraint_error_m(&self) -> PyResult<f64> {
        self.inner
            .system
            .max_gear_constraint_error_m()
            .map_err(util::err)
    }

    /// The largest revolute anchor error, in metres.
    fn max_revolute_anchor_error_m(&self) -> PyResult<f64> {
        self.inner
            .system
            .max_revolute_anchor_error_m()
            .map_err(util::err)
    }

    /// The largest revolute axis error, in radians.
    fn max_revolute_axis_error_rad(&self) -> PyResult<f64> {
        self.inner
            .system
            .max_revolute_axis_error_rad()
            .map_err(util::err)
    }

    /// Exports the assembly as URDF text.
    #[pyo3(signature = (name=None))]
    fn export_urdf(&self, name: Option<&str>) -> PyResult<String> {
        match name {
            Some(name) => jigs::export_urdf_named(&self.inner, name),
            None => jigs::export_urdf(&self.inner),
        }
        .map_err(util::err)
    }
}

/// Assembles a planetary gearbox from a parameter map.
///
/// The keys are the `planetary` generator parameters. Every length is SI.
#[pyfunction]
#[pyo3(signature = (specs=None))]
pub fn assemble_planetary(specs: Option<HashMap<String, f64>>) -> PyResult<PlanetaryAssembly> {
    let mut set = parts::ParameterSet::new();
    for (name, value) in specs.unwrap_or_default() {
        set.set(name, value);
    }
    Ok(PlanetaryAssembly {
        inner: jigs::assemble_planetary(&set).map_err(util::err)?,
    })
}
