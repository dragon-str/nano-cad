//! Conservative jigs: anchors and springs.

use nanocad_jigs::{AnchorJig as AnchorJigInner, Jig, SpringJig as SpringJigInner};
use nanocad_model::Topology as ModelTopology;
use pyo3::prelude::*;

use super::model::Topology;
use super::util;

/// Holds selected atoms at fixed positions with harmonic springs.
#[pyclass(name = "AnchorJig", skip_from_py_object)]
pub struct AnchorJig {
    topology: ModelTopology,
    inner: AnchorJigInner,
}

#[pymethods]
impl AnchorJig {
    #[new]
    fn new(topology: PyRef<'_, Topology>, k_n_per_m: f64) -> PyResult<Self> {
        Ok(Self {
            topology: topology.inner.clone(),
            inner: AnchorJigInner::new(&topology.inner, k_n_per_m).map_err(util::err)?,
        })
    }

    #[getter]
    fn k_n_per_m(&self) -> f64 {
        self.inner.k_n_per_m()
    }

    #[getter]
    fn held_count(&self) -> usize {
        self.inner.held_count()
    }

    /// Holds one atom at a position in SI metres.
    fn hold_atom(&mut self, atom: u32, position_m: Vec<f64>) -> PyResult<()> {
        self.inner
            .hold_atom(&self.topology, atom, util::vec3(position_m)?)
            .map_err(util::err)
    }

    /// Holds one atom at its current position.
    fn hold_atom_current(&mut self, atom: u32) -> PyResult<()> {
        self.inner
            .hold_atom_current(&self.topology, atom)
            .map_err(util::err)
    }

    /// Returns the index of a held atom, or `None`.
    fn held_atom(&self, index: usize) -> Option<u32> {
        self.inner.held_atom(index)
    }

    /// Returns the hold position of a held atom, or `None`.
    fn hold_position_m(&self, index: usize) -> Option<[f64; 3]> {
        self.inner.hold_position_m(index)
    }

    /// Returns the potential energy in joules.
    fn energy(&self, positions_m: Vec<f64>) -> PyResult<f64> {
        self.inner.energy_j(&positions_m).map_err(util::err)
    }

    /// Returns the force in newtons.
    fn forces(&self, positions_m: Vec<f64>) -> PyResult<Vec<f64>> {
        self.inner.forces_n(&positions_m).map_err(util::err)
    }

    /// Returns the gradient in newtons.
    fn gradient(&self, positions_m: Vec<f64>) -> PyResult<Vec<f64>> {
        self.inner.gradient_j_per_m(&positions_m).map_err(util::err)
    }
}

/// Links atoms with harmonic springs, with a fixed rest length.
#[pyclass(name = "SpringJig", skip_from_py_object)]
pub struct SpringJig {
    topology: ModelTopology,
    inner: SpringJigInner,
}

#[pymethods]
impl SpringJig {
    #[new]
    fn new(topology: PyRef<'_, Topology>) -> Self {
        Self {
            topology: topology.inner.clone(),
            inner: SpringJigInner::new(&topology.inner),
        }
    }

    #[getter]
    fn spring_count(&self) -> usize {
        self.inner.spring_count()
    }

    /// Adds a spring between two atoms.
    fn add_spring(&mut self, u: u32, v: u32, k_n_per_m: f64, rest_length_m: f64) -> PyResult<()> {
        self.inner
            .add_spring(&self.topology, u, v, k_n_per_m, rest_length_m)
            .map_err(util::err)
    }

    /// Adds a spring between one atom and a fixed point.
    fn add_spring_to_point(
        &mut self,
        u: u32,
        point_m: Vec<f64>,
        k_n_per_m: f64,
        rest_length_m: f64,
    ) -> PyResult<()> {
        self.inner
            .add_spring_to_point(
                &self.topology,
                u,
                util::vec3(point_m)?,
                k_n_per_m,
                rest_length_m,
            )
            .map_err(util::err)
    }

    /// Returns the stiffness of a spring, or `None`.
    fn k_n_per_m(&self, spring: usize) -> Option<f64> {
        self.inner.k_n_per_m(spring)
    }

    /// Returns the rest length of a spring, or `None`.
    fn rest_length_m(&self, spring: usize) -> Option<f64> {
        self.inner.rest_length_m(spring)
    }

    /// Returns the potential energy in joules.
    fn energy(&self, positions_m: Vec<f64>) -> PyResult<f64> {
        self.inner.energy_j(&positions_m).map_err(util::err)
    }

    /// Returns the force in newtons.
    fn forces(&self, positions_m: Vec<f64>) -> PyResult<Vec<f64>> {
        self.inner.forces_n(&positions_m).map_err(util::err)
    }

    /// Returns the gradient in newtons.
    fn gradient(&self, positions_m: Vec<f64>) -> PyResult<Vec<f64>> {
        self.inner.gradient_j_per_m(&positions_m).map_err(util::err)
    }
}
