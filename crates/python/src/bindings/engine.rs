//! The molecular-mechanics engine: system assembly, minimization, and
//! integration.

use nanocad_engine as engine;
use pyo3::prelude::*;

use super::util;

fn non_periodic() -> engine::PeriodicBox {
    engine::PeriodicBox::non_periodic()
}

/// A force-field system over a fixed atom count.
#[pyclass(name = "System", skip_from_py_object)]
#[derive(Clone)]
pub struct System {
    pub(crate) inner: engine::System,
}

#[pymethods]
impl System {
    #[new]
    #[pyo3(signature = (atom_count, masses_kg=None))]
    fn new(atom_count: usize, masses_kg: Option<Vec<f64>>) -> PyResult<Self> {
        let inner = match masses_kg {
            Some(masses_kg) => {
                engine::System::with_masses_kg(atom_count, &masses_kg).map_err(util::err)?
            }
            None => engine::System::new(atom_count),
        };
        Ok(Self { inner })
    }

    #[getter]
    fn atom_count(&self) -> usize {
        self.inner.atom_count()
    }

    #[getter]
    fn has_masses(&self) -> bool {
        self.inner.has_masses()
    }

    #[getter]
    fn masses_kg(&self) -> Vec<f64> {
        self.inner.masses_kg().to_vec()
    }

    #[getter]
    fn exclusion_count(&self) -> usize {
        self.inner.exclusion_count()
    }

    #[getter]
    fn neighbor_pair_count(&self) -> usize {
        self.inner.neighbor_pair_count()
    }

    /// Sets the masses, one per atom.
    fn set_masses_kg(&mut self, masses_kg: Vec<f64>) -> PyResult<()> {
        self.inner.set_masses_kg(&masses_kg).map_err(util::err)
    }

    /// Sets the neighbour-list skin, in metres.
    fn set_skin_m(&mut self, skin_m: f64) -> PyResult<()> {
        self.inner.set_skin_m(skin_m).map_err(util::err)
    }

    /// Rebuilds the neighbour list from a flat position buffer.
    fn rebuild_neighbors(&mut self, positions_m: Vec<f64>) -> PyResult<()> {
        self.inner
            .rebuild_neighbors(&positions_m)
            .map_err(util::err)
    }

    /// Adds a non-bonded exclusion pair.
    fn add_exclusion(&mut self, i: u32, j: u32) -> PyResult<()> {
        self.inner.add_exclusion(i, j).map_err(util::err)
    }

    /// Sets a harmonic bond-stretch term. Each entry is
    /// `(u, v, k_n_per_m, r0_m)`.
    fn set_bond_stretch(&mut self, bonds: Vec<(u32, u32, f64, f64)>) -> PyResult<()> {
        let mut term = engine::BondStretchTerm::new();
        for (u, v, k_n_per_m, r0_m) in bonds {
            term.add_bond(u, v, k_n_per_m, r0_m).map_err(util::err)?;
        }
        self.inner.set_bond_stretch(term).map_err(util::err)
    }

    /// Sets a harmonic angle-bend term. Each entry is
    /// `(i, j, k, k_j_per_rad2, theta0_rad)`.
    fn set_angle_bend(&mut self, angles: Vec<(u32, u32, u32, f64, f64)>) -> PyResult<()> {
        let mut term = engine::AngleBendTerm::new();
        for (i, j, k, force_constant, theta0_rad) in angles {
            term.add_angle(i, j, k, force_constant, theta0_rad)
                .map_err(util::err)?;
        }
        self.inner.set_angle_bend(term).map_err(util::err)
    }

    /// Sets a torsion term. Each entry is
    /// `(i, j, k, l, v1_j, v2_j, v3_j)`.
    fn set_torsion(&mut self, torsions: Vec<(u32, u32, u32, u32, f64, f64, f64)>) -> PyResult<()> {
        let mut term = engine::TorsionTerm::new();
        for (i, j, k, l, v1_j, v2_j, v3_j) in torsions {
            term.add_torsion(i, j, k, l, v1_j, v2_j, v3_j)
                .map_err(util::err)?;
        }
        self.inner.set_torsion(term).map_err(util::err)
    }

    /// Sets an out-of-plane improper term. Each entry is
    /// `(i, j, k, l, k_j_per_rad2, chi0_rad)`.
    fn set_out_of_plane(&mut self, impropers: Vec<(u32, u32, u32, u32, f64, f64)>) -> PyResult<()> {
        let mut term = engine::OutOfPlaneTerm::new();
        for (i, j, k, l, k_j_per_rad2, chi0_rad) in impropers {
            term.add_improper(i, j, k, l, k_j_per_rad2, chi0_rad)
                .map_err(util::err)?;
        }
        self.inner.set_out_of_plane(term).map_err(util::err)
    }

    /// Sets a Buckingham van der Waals term. `params` holds
    /// `(a_j, b_per_m, c_j_m6)` per atom.
    fn set_van_der_waals(
        &mut self,
        cutoff_m: f64,
        switch_on_m: f64,
        params: Vec<(f64, f64, f64)>,
    ) -> PyResult<()> {
        let cutoff = engine::Cutoff::new(cutoff_m, switch_on_m).map_err(util::err)?;
        let mut term = engine::VanDerWaalsTerm::new(cutoff, non_periodic()).map_err(util::err)?;
        for (a_j, b_per_m, c_j_m6) in params {
            term.add_atom(a_j, b_per_m, c_j_m6).map_err(util::err)?;
        }
        self.inner.set_van_der_waals(term).map_err(util::err)
    }

    /// Sets an electrostatic term from one charge per atom, in coulombs.
    fn set_electrostatic(
        &mut self,
        cutoff_m: f64,
        switch_on_m: f64,
        charges_c: Vec<f64>,
    ) -> PyResult<()> {
        let cutoff = engine::Cutoff::new(cutoff_m, switch_on_m).map_err(util::err)?;
        let term = engine::ElectrostaticTerm::from_charges(&charges_c, cutoff, non_periodic())
            .map_err(util::err)?;
        self.inner.set_electrostatic(term).map_err(util::err)
    }

    /// Sets an orthorhombic periodic box. A zero length axis is not periodic.
    fn set_periodic_box(&mut self, lengths_m: Vec<f64>) -> PyResult<()> {
        let box_m = engine::PeriodicBox::new(util::vec3(lengths_m)?).map_err(util::err)?;
        self.inner.set_periodic_box(box_m).map_err(util::err)
    }

    /// Returns the potential energy in joules.
    fn energy(&mut self, positions_m: Vec<f64>) -> PyResult<f64> {
        self.inner.energy_j(&positions_m).map_err(util::err)
    }

    /// Returns the gradient in newtons.
    fn gradient(&mut self, positions_m: Vec<f64>) -> PyResult<Vec<f64>> {
        self.inner.gradient_j_per_m(&positions_m).map_err(util::err)
    }

    /// Returns the force in newtons.
    fn forces(&mut self, positions_m: Vec<f64>) -> PyResult<Vec<f64>> {
        self.inner.forces_n(&positions_m).map_err(util::err)
    }

    /// Returns `(energy_j, gradient)`.
    fn energy_and_gradient(&mut self, positions_m: Vec<f64>) -> PyResult<(f64, Vec<f64>)> {
        self.inner
            .energy_and_gradient_j(&positions_m)
            .map_err(util::err)
    }

    /// Returns the kinetic energy in joules.
    fn kinetic_energy(&self, velocities_m_per_s: Vec<f64>) -> PyResult<f64> {
        self.inner
            .kinetic_energy_j(&velocities_m_per_s)
            .map_err(util::err)
    }

    /// Returns the temperature in kelvin.
    fn temperature(&self, velocities_m_per_s: Vec<f64>) -> PyResult<f64> {
        self.inner
            .temperature_k(&velocities_m_per_s)
            .map_err(util::err)
    }

    /// Returns `kinetic + potential`.
    fn total_energy(
        &mut self,
        positions_m: Vec<f64>,
        velocities_m_per_s: Vec<f64>,
    ) -> PyResult<f64> {
        self.inner
            .total_energy_j(&positions_m, &velocities_m_per_s)
            .map_err(util::err)
    }
}

/// The outcome of a minimization run.
#[pyclass(name = "MinimizeResult", skip_from_py_object)]
#[derive(Clone)]
pub struct MinimizeResult {
    /// The final potential energy in joules.
    #[pyo3(get)]
    energy_j: f64,
    /// The final gradient infinity norm in newtons.
    #[pyo3(get)]
    gradient_norm_n: f64,
    /// The number of outer iterations.
    #[pyo3(get)]
    iterations: usize,
    /// Whether the gradient norm reached the tolerance.
    #[pyo3(get)]
    converged: bool,
}

fn run_minimize(
    system: &mut engine::System,
    mut positions_m: Vec<f64>,
    options: &engine::MinimizeOptions,
) -> PyResult<(Vec<f64>, MinimizeResult)> {
    let result = engine::minimize(system, &mut positions_m, options).map_err(util::err)?;
    Ok((
        positions_m,
        MinimizeResult {
            energy_j: result.energy_j,
            gradient_norm_n: result.gradient_norm_n,
            iterations: result.iterations,
            converged: result.converged,
        },
    ))
}

/// Minimizes the potential energy. Returns `(positions_m, result)`.
#[pyfunction]
#[pyo3(signature = (system, positions_m, max_iterations=200, gradient_tolerance_n=1.0e-14, initial_step_m=1.0e-13))]
pub fn minimize(
    mut system: PyRefMut<'_, System>,
    positions_m: Vec<f64>,
    max_iterations: usize,
    gradient_tolerance_n: f64,
    initial_step_m: f64,
) -> PyResult<(Vec<f64>, MinimizeResult)> {
    let options = engine::MinimizeOptions {
        max_iterations,
        gradient_tolerance_n,
        initial_step_m,
    };
    run_minimize(&mut system.inner, positions_m, &options)
}

/// The alias of [`minimize`] used by the agent tool surface.
#[pyfunction]
#[pyo3(signature = (system, positions_m, max_iterations=200, gradient_tolerance_n=1.0e-14, initial_step_m=1.0e-13))]
pub fn relax(
    system: PyRefMut<'_, System>,
    positions_m: Vec<f64>,
    max_iterations: usize,
    gradient_tolerance_n: f64,
    initial_step_m: f64,
) -> PyResult<(Vec<f64>, MinimizeResult)> {
    minimize(
        system,
        positions_m,
        max_iterations,
        gradient_tolerance_n,
        initial_step_m,
    )
}

/// Velocity Verlet integration.
#[pyclass(name = "VelocityVerlet", skip_from_py_object)]
#[derive(Clone)]
pub struct VelocityVerlet {
    inner: engine::VelocityVerlet,
}

#[pymethods]
impl VelocityVerlet {
    #[new]
    fn new(dt_s: f64) -> PyResult<Self> {
        Ok(Self {
            inner: engine::VelocityVerlet::new(dt_s).map_err(util::err)?,
        })
    }

    #[getter]
    fn dt_s(&self) -> f64 {
        self.inner.dt_s()
    }

    /// Advances the system by one step, in place.
    fn step(
        &self,
        mut system: PyRefMut<'_, System>,
        mut positions_m: Vec<f64>,
        mut velocities_m_per_s: Vec<f64>,
    ) -> PyResult<(Vec<f64>, Vec<f64>)> {
        self.inner
            .step(&mut system.inner, &mut positions_m, &mut velocities_m_per_s)
            .map_err(util::err)?;
        Ok((positions_m, velocities_m_per_s))
    }
}

fn run_integrate(
    system: &mut engine::System,
    mut positions_m: Vec<f64>,
    mut velocities_m_per_s: Vec<f64>,
    dt_s: f64,
    steps: usize,
) -> PyResult<(Vec<f64>, Vec<f64>)> {
    let integrator = engine::VelocityVerlet::new(dt_s).map_err(util::err)?;
    for _ in 0..steps {
        integrator
            .step(system, &mut positions_m, &mut velocities_m_per_s)
            .map_err(util::err)?;
    }
    Ok((positions_m, velocities_m_per_s))
}

/// Integrates for `steps` Velocity Verlet steps. Returns
/// `(positions_m, velocities_m_per_s)`.
#[pyfunction]
pub fn integrate(
    mut system: PyRefMut<'_, System>,
    positions_m: Vec<f64>,
    velocities_m_per_s: Vec<f64>,
    dt_s: f64,
    steps: usize,
) -> PyResult<(Vec<f64>, Vec<f64>)> {
    run_integrate(
        &mut system.inner,
        positions_m,
        velocities_m_per_s,
        dt_s,
        steps,
    )
}

/// The alias of [`integrate`] used by the agent tool surface.
#[pyfunction]
pub fn run_md(
    system: PyRefMut<'_, System>,
    positions_m: Vec<f64>,
    velocities_m_per_s: Vec<f64>,
    dt_s: f64,
    steps: usize,
) -> PyResult<(Vec<f64>, Vec<f64>)> {
    integrate(system, positions_m, velocities_m_per_s, dt_s, steps)
}

/// The Coulomb constant, in SI units.
#[pyfunction]
pub fn coulomb_constant() -> f64 {
    engine::COULOMB_CONSTANT_N_M2_PER_C2
}
