//! The parameter store: quantities, provenance, records, and libraries.

use nanocad_params as params;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use super::engine;
use super::util;

fn quantity(inner: params::Quantity) -> ParamQuantity {
    ParamQuantity { inner }
}

/// A physical value in SI base units with a presentation unit and an optional
/// uncertainty. An unknown uncertainty is `None`, never zero.
#[pyclass(name = "ParamQuantity", skip_from_py_object)]
#[derive(Clone)]
pub struct ParamQuantity {
    pub(crate) inner: params::Quantity,
}

#[pymethods]
impl ParamQuantity {
    /// Builds a quantity. A base unit symbol converts the value to SI. A
    /// derived symbol such as `"Pa"` takes the value already in SI.
    #[new]
    fn new(value: f64, unit: &str) -> PyResult<Self> {
        let inner = match util::unit_from_symbol(unit) {
            Ok(base) => params::Quantity::new(value, base),
            Err(_) => params::Quantity::derived(value, unit).map_err(util::err)?,
        };
        Ok(Self { inner })
    }

    /// Builds a quantity from a value already in SI base units.
    #[staticmethod]
    fn from_si(value_si: f64, unit: &str) -> PyResult<Self> {
        if let Ok(base) = util::unit_from_symbol(unit) {
            return Ok(Self {
                inner: params::Quantity::from_si(value_si, base),
            });
        }
        Ok(Self {
            inner: params::Quantity::derived(value_si, unit).map_err(util::err)?,
        })
    }

    /// Returns a copy with an uncertainty in the presentation unit.
    fn with_uncertainty(&self, uncertainty: f64) -> Self {
        Self {
            inner: self.inner.clone().with_uncertainty(uncertainty),
        }
    }

    /// Returns a copy with an uncertainty already in SI base units.
    fn with_uncertainty_si(&self, uncertainty_si: f64) -> Self {
        Self {
            inner: self.inner.clone().with_uncertainty_si(uncertainty_si),
        }
    }

    /// The value in the presentation unit.
    #[getter]
    fn value(&self) -> f64 {
        self.inner.value()
    }

    /// The value in SI base units.
    #[getter]
    fn value_si(&self) -> f64 {
        self.inner.value_si()
    }

    /// The uncertainty in the presentation unit, or `None`.
    #[getter]
    fn uncertainty(&self) -> Option<f64> {
        self.inner.uncertainty()
    }

    /// The uncertainty in SI base units, or `None`.
    #[getter]
    fn uncertainty_si(&self) -> Option<f64> {
        self.inner.uncertainty_si()
    }

    #[getter]
    fn unit(&self) -> String {
        self.inner.unit().to_owned()
    }

    /// The seven dimension exponents, or `None` for an unknown unit.
    #[getter]
    fn dimension(&self) -> Option<Vec<i8>> {
        self.inner.dimension().map(|dimension| {
            vec![
                dimension.length,
                dimension.mass,
                dimension.time,
                dimension.current,
                dimension.temperature,
                dimension.amount,
                dimension.luminous_intensity,
            ]
        })
    }

    fn __repr__(&self) -> String {
        match self.inner.uncertainty() {
            Some(uncertainty) => format!(
                "ParamQuantity({}, unit={}, uncertainty={})",
                self.inner.value(),
                self.inner.unit(),
                uncertainty
            ),
            None => format!(
                "ParamQuantity({}, unit={})",
                self.inner.value(),
                self.inner.unit()
            ),
        }
    }
}

/// How a value or a record was produced.
#[pyclass(name = "Provenance", skip_from_py_object)]
#[derive(Clone)]
pub struct Provenance {
    pub(crate) inner: params::Provenance,
}

#[pymethods]
impl Provenance {
    #[new]
    #[pyo3(signature = (source, method, code_version, timestamp, uncertainty=None, validation="unverified", force_field=None, notes=""))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        source: String,
        method: &str,
        code_version: String,
        timestamp: String,
        uncertainty: Option<PyRef<'_, ParamQuantity>>,
        validation: &str,
        force_field: Option<String>,
        notes: &str,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: params::Provenance {
                source,
                method: util::method_from_str(method)?,
                code_version,
                force_field,
                timestamp,
                uncertainty: uncertainty.map(|quantity| quantity.inner.clone()),
                validation: util::validation_from_str(validation)?,
                notes: notes.to_owned(),
            },
        })
    }

    #[getter]
    fn source(&self) -> String {
        self.inner.source.clone()
    }

    #[getter]
    fn method(&self) -> String {
        util::method_name(self.inner.method).to_owned()
    }

    #[getter]
    fn code_version(&self) -> String {
        self.inner.code_version.clone()
    }

    #[getter]
    fn force_field(&self) -> Option<String> {
        self.inner.force_field.clone()
    }

    #[getter]
    fn timestamp(&self) -> String {
        self.inner.timestamp.clone()
    }

    #[getter]
    fn uncertainty(&self) -> Option<ParamQuantity> {
        self.inner.uncertainty.clone().map(quantity)
    }

    #[getter]
    fn validation(&self) -> String {
        util::validation_name(self.inner.validation).to_owned()
    }

    #[getter]
    fn notes(&self) -> String {
        self.inner.notes.clone()
    }
}

/// The measured or computed properties of one part.
#[pyclass(name = "PartRecord", skip_from_py_object)]
#[derive(Clone)]
pub struct PartRecord {
    pub(crate) inner: params::PartRecord,
}

fn take_quantity(quantity: &PyRef<'_, ParamQuantity>) -> params::Quantity {
    quantity.inner.clone()
}

#[pymethods]
impl PartRecord {
    #[new]
    #[pyo3(signature = (part_id, geometry_ref, material, atoms, mass_kg, inertia_kg_m2, elastic_modulus_pa, shear_modulus_pa, poisson_ratio, failure_stress_pa, friction_coefficient, thermal_conductivity_w_m_k, specific_heat_j_kg_k, provenance, method=None, validation=None))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        part_id: String,
        geometry_ref: String,
        material: String,
        atoms: u64,
        mass_kg: PyRef<'_, ParamQuantity>,
        inertia_kg_m2: Vec<PyRef<'_, ParamQuantity>>,
        elastic_modulus_pa: PyRef<'_, ParamQuantity>,
        shear_modulus_pa: PyRef<'_, ParamQuantity>,
        poisson_ratio: PyRef<'_, ParamQuantity>,
        failure_stress_pa: PyRef<'_, ParamQuantity>,
        friction_coefficient: PyRef<'_, ParamQuantity>,
        thermal_conductivity_w_m_k: PyRef<'_, ParamQuantity>,
        specific_heat_j_kg_k: PyRef<'_, ParamQuantity>,
        provenance: PyRef<'_, Provenance>,
        method: Option<&str>,
        validation: Option<&str>,
    ) -> PyResult<Self> {
        if inertia_kg_m2.len() != 3 {
            return Err(PyValueError::new_err(
                "inertia_kg_m2 must hold three quantities",
            ));
        }
        let method = match method {
            Some(name) => util::method_from_str(name)?,
            None => provenance.inner.method,
        };
        let validation = match validation {
            Some(name) => util::validation_from_str(name)?,
            None => provenance.inner.validation,
        };
        Ok(Self {
            inner: params::PartRecord {
                schema: params::SCHEMA.to_owned(),
                version: params::VERSION,
                part_id,
                geometry_ref,
                material,
                atoms,
                mass_kg: take_quantity(&mass_kg),
                inertia_kg_m2: [
                    take_quantity(&inertia_kg_m2[0]),
                    take_quantity(&inertia_kg_m2[1]),
                    take_quantity(&inertia_kg_m2[2]),
                ],
                elastic_modulus_pa: take_quantity(&elastic_modulus_pa),
                shear_modulus_pa: take_quantity(&shear_modulus_pa),
                poisson_ratio: take_quantity(&poisson_ratio),
                failure_stress_pa: take_quantity(&failure_stress_pa),
                friction_coefficient: take_quantity(&friction_coefficient),
                thermal_conductivity_w_m_k: take_quantity(&thermal_conductivity_w_m_k),
                specific_heat_j_kg_k: take_quantity(&specific_heat_j_kg_k),
                method,
                validation,
                provenance: provenance.inner.clone(),
            },
        })
    }

    #[getter]
    fn schema(&self) -> String {
        self.inner.schema.clone()
    }

    #[getter]
    fn version(&self) -> u32 {
        self.inner.version
    }

    #[getter]
    fn part_id(&self) -> String {
        self.inner.part_id.clone()
    }

    #[getter]
    fn geometry_ref(&self) -> String {
        self.inner.geometry_ref.clone()
    }

    #[getter]
    fn material(&self) -> String {
        self.inner.material.clone()
    }

    #[getter]
    fn atoms(&self) -> u64 {
        self.inner.atoms
    }

    #[getter]
    fn mass_kg(&self) -> ParamQuantity {
        quantity(self.inner.mass_kg.clone())
    }

    #[getter]
    fn inertia_kg_m2(&self) -> Vec<ParamQuantity> {
        self.inner
            .inertia_kg_m2
            .iter()
            .cloned()
            .map(quantity)
            .collect()
    }

    #[getter]
    fn elastic_modulus_pa(&self) -> ParamQuantity {
        quantity(self.inner.elastic_modulus_pa.clone())
    }

    #[getter]
    fn shear_modulus_pa(&self) -> ParamQuantity {
        quantity(self.inner.shear_modulus_pa.clone())
    }

    #[getter]
    fn poisson_ratio(&self) -> ParamQuantity {
        quantity(self.inner.poisson_ratio.clone())
    }

    #[getter]
    fn failure_stress_pa(&self) -> ParamQuantity {
        quantity(self.inner.failure_stress_pa.clone())
    }

    #[getter]
    fn friction_coefficient(&self) -> ParamQuantity {
        quantity(self.inner.friction_coefficient.clone())
    }

    #[getter]
    fn thermal_conductivity_w_m_k(&self) -> ParamQuantity {
        quantity(self.inner.thermal_conductivity_w_m_k.clone())
    }

    #[getter]
    fn specific_heat_j_kg_k(&self) -> ParamQuantity {
        quantity(self.inner.specific_heat_j_kg_k.clone())
    }

    #[getter]
    fn method(&self) -> String {
        util::method_name(self.inner.method).to_owned()
    }

    #[getter]
    fn validation(&self) -> String {
        util::validation_name(self.inner.validation).to_owned()
    }

    #[getter]
    fn provenance(&self) -> Provenance {
        Provenance {
            inner: self.inner.provenance.clone(),
        }
    }

    /// Serializes the record to compact JSON.
    fn to_json(&self) -> PyResult<String> {
        self.inner.to_json().map_err(util::err)
    }

    /// Serializes the record to indented JSON.
    fn to_json_pretty(&self) -> PyResult<String> {
        self.inner.to_json_pretty().map_err(util::err)
    }

    /// Serializes the record to a JSON byte buffer.
    fn to_bytes(&self) -> PyResult<Vec<u8>> {
        self.inner.to_bytes().map_err(util::err)
    }

    /// Deserializes a record from JSON and checks its header.
    #[staticmethod]
    fn from_json(json: &str) -> PyResult<Self> {
        Ok(Self {
            inner: params::PartRecord::from_json(json).map_err(util::err)?,
        })
    }

    /// Deserializes a record from a JSON byte buffer and checks its header.
    #[staticmethod]
    fn from_bytes(bytes: &[u8]) -> PyResult<Self> {
        Ok(Self {
            inner: params::PartRecord::from_bytes(bytes).map_err(util::err)?,
        })
    }

    /// Rejects a foreign schema name or an unsupported version.
    fn check_header(&self) -> PyResult<()> {
        self.inner.check_header().map_err(util::err)
    }
}

/// A versioned store of many records, keyed by part id and version.
#[pyclass(name = "ParameterLibrary", skip_from_py_object)]
#[derive(Clone)]
pub struct ParameterLibrary {
    inner: params::ParameterLibrary,
}

#[pymethods]
impl ParameterLibrary {
    #[new]
    fn new() -> Self {
        Self {
            inner: params::ParameterLibrary::new(),
        }
    }

    #[getter]
    fn schema(&self) -> String {
        self.inner.schema().to_owned()
    }

    #[getter(version)]
    fn library_version(&self) -> u32 {
        self.inner.version()
    }

    #[getter]
    fn revision(&self) -> u64 {
        self.inner.revision()
    }

    #[getter]
    fn len(&self) -> usize {
        self.inner.len()
    }

    #[getter]
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// Returns every part id, sorted.
    fn part_ids(&self) -> Vec<String> {
        self.inner.part_ids().map(str::to_owned).collect()
    }

    /// Inserts one record at a version.
    fn insert(&mut self, record: PyRef<'_, PartRecord>, version: u32) -> PyResult<()> {
        self.inner
            .insert(record.inner.clone(), version)
            .map_err(util::err)
    }

    /// Returns the latest revision of a part, or `None`.
    fn get(&self, part_id: &str) -> Option<PartRecord> {
        self.inner.get(part_id).map(|record| PartRecord {
            inner: record.clone(),
        })
    }

    /// Returns one version of a part, or `None`.
    fn get_version(&self, part_id: &str, version: u32) -> Option<PartRecord> {
        self.inner
            .get_version(part_id, version)
            .map(|record| PartRecord {
                inner: record.clone(),
            })
    }

    /// Returns the latest version number of a part, or `None`.
    fn latest_version(&self, part_id: &str) -> Option<u32> {
        self.inner.latest_version(part_id)
    }

    /// Returns every stored record, in key order.
    fn list(&self) -> Vec<PartRecord> {
        self.inner
            .list()
            .into_iter()
            .map(|record| PartRecord {
                inner: record.clone(),
            })
            .collect()
    }

    /// Removes every version of a part. Returns whether it was present.
    fn remove(&mut self, part_id: &str) -> bool {
        self.inner.remove(part_id)
    }

    /// Removes one version of a part and returns it, or `None`.
    fn remove_version(&mut self, part_id: &str, version: u32) -> Option<PartRecord> {
        self.inner
            .remove_version(part_id, version)
            .map(|record| PartRecord { inner: record })
    }

    /// Checks a stored record against an expected geometry hash.
    fn verify_geometry(
        &self,
        part_id: &str,
        version: u32,
        expected_geometry_ref: &str,
    ) -> PyResult<()> {
        self.inner
            .verify_geometry(part_id, version, expected_geometry_ref)
            .map_err(util::err)
    }

    /// Rejects a foreign library schema or version.
    fn check_header(&self) -> PyResult<()> {
        self.inner.check_header().map_err(util::err)
    }

    fn to_json(&self) -> PyResult<String> {
        self.inner.to_json().map_err(util::err)
    }

    fn to_json_pretty(&self) -> PyResult<String> {
        self.inner.to_json_pretty().map_err(util::err)
    }

    fn to_bytes(&self) -> PyResult<Vec<u8>> {
        self.inner.to_bytes().map_err(util::err)
    }

    #[staticmethod]
    fn from_json(json: &str) -> PyResult<Self> {
        Ok(Self {
            inner: params::ParameterLibrary::from_json(json).map_err(util::err)?,
        })
    }

    #[staticmethod]
    fn from_bytes(bytes: &[u8]) -> PyResult<Self> {
        Ok(Self {
            inner: params::ParameterLibrary::from_bytes(bytes).map_err(util::err)?,
        })
    }

    /// Saves the library to a file.
    fn save(&self, path: &str) -> PyResult<()> {
        self.inner.save_to_file(path).map_err(util::err)
    }

    /// Loads a library from a file.
    #[staticmethod]
    fn load(path: &str) -> PyResult<Self> {
        Ok(Self {
            inner: params::ParameterLibrary::load_from_file(path).map_err(util::err)?,
        })
    }
}

/// Checks one record and returns every problem as a readable string.
#[pyfunction]
pub fn check_record(record: PyRef<'_, PartRecord>) -> Vec<String> {
    params::check(&record.inner)
        .iter()
        .map(|problem| problem.to_string())
        .collect()
}

/// Reports whether a record is consistent.
#[pyfunction]
pub fn is_consistent(record: PyRef<'_, PartRecord>) -> bool {
    params::is_consistent(&record.inner)
}

fn provenance(inner: params::Provenance) -> Provenance {
    Provenance { inner }
}

/// The elastic modulus that a strain sweep produced.
#[pyclass(name = "StiffnessResult", skip_from_py_object)]
#[derive(Clone)]
pub struct StiffnessResult {
    inner: params::StiffnessResult,
}

#[pymethods]
impl StiffnessResult {
    #[getter]
    fn elastic_modulus_pa(&self) -> ParamQuantity {
        quantity(self.inner.elastic_modulus_pa.clone())
    }

    #[getter]
    fn fit_intercept_pa(&self) -> f64 {
        self.inner.fit_intercept_pa
    }

    #[getter]
    fn strain_min(&self) -> f64 {
        self.inner.strain_range.0
    }

    #[getter]
    fn strain_max(&self) -> f64 {
        self.inner.strain_range.1
    }

    #[getter]
    fn samples(&self) -> usize {
        self.inner.samples
    }

    #[getter]
    fn r_squared(&self) -> f64 {
        self.inner.r_squared
    }

    #[getter]
    fn provenance(&self) -> Provenance {
        provenance(self.inner.provenance.clone())
    }
}

/// Extracts an elastic modulus from a controlled strain sweep.
#[pyfunction]
#[pyo3(signature = (system, positions_m, strain_min=-1.0e-3, strain_max=1.0e-3, samples=11, reference_length_m=3.0e-9, cross_section_area_m2=1.0e-20, derivative_step=1.0e-5, noise_amplitude_pa=0.0, seed=0x5EED_1234))]
#[allow(clippy::too_many_arguments)]
pub fn extract_stiffness(
    mut system: PyRefMut<'_, engine::System>,
    positions_m: Vec<f64>,
    strain_min: f64,
    strain_max: f64,
    samples: usize,
    reference_length_m: f64,
    cross_section_area_m2: f64,
    derivative_step: f64,
    noise_amplitude_pa: f64,
    seed: u64,
) -> PyResult<StiffnessResult> {
    let config = params::StiffnessConfig {
        strain_min,
        strain_max,
        samples,
        reference_length_m,
        cross_section_area_m2,
        derivative_step,
        noise_amplitude_pa,
        seed,
    };
    let inner =
        params::extract_stiffness(&mut system.inner, &positions_m, &config).map_err(util::err)?;
    Ok(StiffnessResult { inner })
}

/// The friction coefficient that a driven sliding run produced.
#[pyclass(name = "FrictionResult", skip_from_py_object)]
#[derive(Clone)]
pub struct FrictionResult {
    inner: params::FrictionResult,
}

#[pymethods]
impl FrictionResult {
    #[getter]
    fn friction_coefficient(&self) -> ParamQuantity {
        quantity(self.inner.friction_coefficient.clone())
    }

    #[getter]
    fn friction_force_mean_n(&self) -> f64 {
        self.inner.friction_force_mean_n
    }

    #[getter]
    fn friction_force_spread_n(&self) -> f64 {
        self.inner.friction_force_spread_n
    }

    #[getter]
    fn normal_load_n(&self) -> f64 {
        self.inner.normal_load_n
    }

    #[getter]
    fn samples(&self) -> usize {
        self.inner.samples
    }

    #[getter]
    fn provenance(&self) -> Provenance {
        provenance(self.inner.provenance.clone())
    }
}

/// Extracts a friction coefficient from a driven sliding contact.
#[pyfunction]
#[pyo3(signature = (system, positions_m, slider_atom=0, normal_load_n=1.0e-11, drag_n_s_per_m=2.0e-13, stage_velocity_m_per_s=1.0, drive_stiffness_n_per_m=1.0e-1, dt_s=1.0e-14, equilibration_steps=10_000, measurement_steps=50_000))]
#[allow(clippy::too_many_arguments)]
pub fn extract_friction(
    mut system: PyRefMut<'_, engine::System>,
    positions_m: Vec<f64>,
    slider_atom: u32,
    normal_load_n: f64,
    drag_n_s_per_m: f64,
    stage_velocity_m_per_s: f64,
    drive_stiffness_n_per_m: f64,
    dt_s: f64,
    equilibration_steps: usize,
    measurement_steps: usize,
) -> PyResult<FrictionResult> {
    let config = params::FrictionConfig {
        slider_atom,
        normal_load_n,
        drag_n_s_per_m,
        stage_velocity_m_per_s,
        drive_stiffness_n_per_m,
        dt_s,
        equilibration_steps,
        measurement_steps,
    };
    let inner =
        params::extract_friction(&mut system.inner, &positions_m, &config).map_err(util::err)?;
    Ok(FrictionResult { inner })
}

/// Summarizes signed lateral-force samples into a friction result.
#[pyfunction]
pub fn summarize_friction(
    lateral_forces_n: Vec<f64>,
    normal_load_n: f64,
) -> PyResult<FrictionResult> {
    let inner = params::summarize_friction(&lateral_forces_n, normal_load_n).map_err(util::err)?;
    Ok(FrictionResult { inner })
}

/// The failure stress that a load sweep produced.
#[pyclass(name = "FailureResult", skip_from_py_object)]
#[derive(Clone)]
pub struct FailureResult {
    inner: params::FailureResult,
}

#[pymethods]
impl FailureResult {
    #[getter]
    fn failure_stress_pa(&self) -> ParamQuantity {
        quantity(self.inner.failure_stress_pa.clone())
    }

    #[getter]
    fn failure_strain(&self) -> f64 {
        self.inner.failure_strain
    }

    #[getter]
    fn critical_bond(&self) -> usize {
        self.inner.critical_bond
    }

    #[getter]
    fn critical_bond_strain(&self) -> f64 {
        self.inner.critical_bond_strain
    }

    #[getter]
    fn provenance(&self) -> Provenance {
        provenance(self.inner.provenance.clone())
    }
}

fn failure_config(
    strain_step: f64,
    max_strain: f64,
    reference_length_m: f64,
    cross_section_area_m2: f64,
    derivative_step: f64,
) -> params::FailureConfig {
    params::FailureConfig {
        strain_step,
        max_strain,
        reference_length_m,
        cross_section_area_m2,
        derivative_step,
    }
}

/// Extracts a failure stress with explicit per-bond limits.
#[pyfunction]
#[pyo3(signature = (system, positions_m, bonds, strain_step=1.0e-3, max_strain=0.2, reference_length_m=3.0e-9, cross_section_area_m2=1.0e-20, derivative_step=1.0e-5))]
#[allow(clippy::too_many_arguments)]
pub fn extract_failure_stress(
    mut system: PyRefMut<'_, engine::System>,
    positions_m: Vec<f64>,
    bonds: Vec<(u32, u32, f64, f64)>,
    strain_step: f64,
    max_strain: f64,
    reference_length_m: f64,
    cross_section_area_m2: f64,
    derivative_step: f64,
) -> PyResult<FailureResult> {
    let limits: Vec<params::BondLimit> = bonds
        .into_iter()
        .map(|(u, v, r0_m, rupture_strain)| params::BondLimit {
            u,
            v,
            r0_m,
            rupture_strain,
        })
        .collect();
    let config = failure_config(
        strain_step,
        max_strain,
        reference_length_m,
        cross_section_area_m2,
        derivative_step,
    );
    let inner = params::extract_failure_stress(&mut system.inner, &positions_m, &config, &limits)
        .map_err(util::err)?;
    Ok(FailureResult { inner })
}

/// Extracts a failure stress with the bond limits taken from the system.
#[pyfunction]
#[pyo3(signature = (system, positions_m, rupture_strain, strain_step=1.0e-3, max_strain=0.2, reference_length_m=3.0e-9, cross_section_area_m2=1.0e-20, derivative_step=1.0e-5))]
#[allow(clippy::too_many_arguments)]
pub fn extract_failure_stress_from_system(
    mut system: PyRefMut<'_, engine::System>,
    positions_m: Vec<f64>,
    rupture_strain: f64,
    strain_step: f64,
    max_strain: f64,
    reference_length_m: f64,
    cross_section_area_m2: f64,
    derivative_step: f64,
) -> PyResult<FailureResult> {
    let config = failure_config(
        strain_step,
        max_strain,
        reference_length_m,
        cross_section_area_m2,
        derivative_step,
    );
    let inner = params::extract_failure_stress_from_system(
        &mut system.inner,
        &positions_m,
        &config,
        rupture_strain,
    )
    .map_err(util::err)?;
    Ok(FailureResult { inner })
}

/// The specific heat that an NVT run produced.
#[pyclass(name = "SpecificHeatResult", skip_from_py_object)]
#[derive(Clone)]
pub struct SpecificHeatResult {
    inner: params::SpecificHeatResult,
}

#[pymethods]
impl SpecificHeatResult {
    #[getter]
    fn specific_heat_j_kg_k(&self) -> ParamQuantity {
        quantity(self.inner.specific_heat_j_kg_k.clone())
    }

    #[getter]
    fn mean_temperature_k(&self) -> f64 {
        self.inner.mean_temperature_k
    }

    #[getter]
    fn mean_kinetic_energy_j(&self) -> f64 {
        self.inner.mean_kinetic_energy_j
    }

    #[getter]
    fn samples(&self) -> usize {
        self.inner.samples
    }

    #[getter]
    fn mass_kg(&self) -> f64 {
        self.inner.mass_kg
    }

    #[getter]
    fn atoms(&self) -> usize {
        self.inner.atoms
    }

    #[getter]
    fn provenance(&self) -> Provenance {
        provenance(self.inner.provenance.clone())
    }
}

/// The thermal conductivity that a two-bath run produced.
#[pyclass(name = "ConductivityResult", skip_from_py_object)]
#[derive(Clone)]
pub struct ConductivityResult {
    inner: params::ConductivityResult,
}

#[pymethods]
impl ConductivityResult {
    #[getter]
    fn thermal_conductivity_w_m_k(&self) -> ParamQuantity {
        quantity(self.inner.thermal_conductivity_w_m_k.clone())
    }

    #[getter]
    fn temperature_gradient_k_per_m(&self) -> f64 {
        self.inner.temperature_gradient_k_per_m
    }

    #[getter]
    fn heat_current_w(&self) -> f64 {
        self.inner.heat_current_w
    }

    #[getter]
    fn hot_temperature_k(&self) -> f64 {
        self.inner.hot_temperature_k
    }

    #[getter]
    fn cold_temperature_k(&self) -> f64 {
        self.inner.cold_temperature_k
    }

    #[getter]
    fn slab_distance_m(&self) -> f64 {
        self.inner.slab_distance_m
    }

    #[getter]
    fn samples(&self) -> usize {
        self.inner.samples
    }

    #[getter]
    fn provenance(&self) -> Provenance {
        provenance(self.inner.provenance.clone())
    }
}

/// The two thermal properties of one part, from one run.
#[pyclass(name = "ThermalResult", skip_from_py_object)]
#[derive(Clone)]
pub struct ThermalResult {
    inner: params::ThermalResult,
}

#[pymethods]
impl ThermalResult {
    #[getter]
    fn specific_heat_j_kg_k(&self) -> ParamQuantity {
        quantity(self.inner.specific_heat_j_kg_k.clone())
    }

    #[getter]
    fn thermal_conductivity_w_m_k(&self) -> ParamQuantity {
        quantity(self.inner.thermal_conductivity_w_m_k.clone())
    }

    #[getter]
    fn mean_temperature_k(&self) -> f64 {
        self.inner.mean_temperature_k
    }

    #[getter]
    fn temperature_gradient_k_per_m(&self) -> f64 {
        self.inner.temperature_gradient_k_per_m
    }

    #[getter]
    fn heat_current_w(&self) -> f64 {
        self.inner.heat_current_w
    }

    #[getter]
    fn kinetic_energy_samples(&self) -> usize {
        self.inner.kinetic_energy_samples
    }

    #[getter]
    fn heat_current_samples(&self) -> usize {
        self.inner.heat_current_samples
    }

    #[getter]
    fn provenance(&self) -> Provenance {
        provenance(self.inner.provenance.clone())
    }
}

#[allow(clippy::too_many_arguments)]
fn thermal_config(
    dt_s: f64,
    temperature_k: f64,
    friction_per_s: f64,
    berendsen_tau_s: f64,
    equilibration_steps: usize,
    production_steps: usize,
    sample_interval: usize,
    hot_temperature_k: f64,
    cold_temperature_k: f64,
    bath_friction_per_s: f64,
    conductivity_equilibration_steps: usize,
    conductivity_production_steps: usize,
    conductivity_sample_interval: usize,
    slab_atoms: usize,
    cross_section_area_m2: f64,
    seed: u64,
) -> params::ThermalConfig {
    params::ThermalConfig {
        dt_s,
        temperature_k,
        friction_per_s,
        berendsen_tau_s,
        equilibration_steps,
        production_steps,
        sample_interval,
        hot_temperature_k,
        cold_temperature_k,
        bath_friction_per_s,
        conductivity_equilibration_steps,
        conductivity_production_steps,
        conductivity_sample_interval,
        slab_atoms,
        cross_section_area_m2,
        seed,
    }
}

/// Extracts the specific heat from the kinetic-energy fluctuation under NVT.
#[pyfunction]
#[pyo3(signature = (system, positions_m, dt_s=5.0e-15, temperature_k=300.0, friction_per_s=5.0e12, berendsen_tau_s=2.0e-13, equilibration_steps=4_000, production_steps=80_000, sample_interval=20, hot_temperature_k=400.0, cold_temperature_k=200.0, bath_friction_per_s=1.0e13, conductivity_equilibration_steps=10_000, conductivity_production_steps=40_000, conductivity_sample_interval=1, slab_atoms=2, cross_section_area_m2=1.0e-20, seed=0x7E57_1CE5))]
#[allow(clippy::too_many_arguments)]
pub fn extract_specific_heat(
    mut system: PyRefMut<'_, engine::System>,
    positions_m: Vec<f64>,
    dt_s: f64,
    temperature_k: f64,
    friction_per_s: f64,
    berendsen_tau_s: f64,
    equilibration_steps: usize,
    production_steps: usize,
    sample_interval: usize,
    hot_temperature_k: f64,
    cold_temperature_k: f64,
    bath_friction_per_s: f64,
    conductivity_equilibration_steps: usize,
    conductivity_production_steps: usize,
    conductivity_sample_interval: usize,
    slab_atoms: usize,
    cross_section_area_m2: f64,
    seed: u64,
) -> PyResult<SpecificHeatResult> {
    let config = thermal_config(
        dt_s,
        temperature_k,
        friction_per_s,
        berendsen_tau_s,
        equilibration_steps,
        production_steps,
        sample_interval,
        hot_temperature_k,
        cold_temperature_k,
        bath_friction_per_s,
        conductivity_equilibration_steps,
        conductivity_production_steps,
        conductivity_sample_interval,
        slab_atoms,
        cross_section_area_m2,
        seed,
    );
    let inner = params::extract_specific_heat(&mut system.inner, &positions_m, &config)
        .map_err(util::err)?;
    Ok(SpecificHeatResult { inner })
}

/// Extracts a thermal conductivity from a two-bath non-equilibrium run.
#[pyfunction]
#[pyo3(signature = (system, positions_m, dt_s=5.0e-15, temperature_k=300.0, friction_per_s=5.0e12, berendsen_tau_s=2.0e-13, equilibration_steps=4_000, production_steps=80_000, sample_interval=20, hot_temperature_k=400.0, cold_temperature_k=200.0, bath_friction_per_s=1.0e13, conductivity_equilibration_steps=10_000, conductivity_production_steps=40_000, conductivity_sample_interval=1, slab_atoms=2, cross_section_area_m2=1.0e-20, seed=0x7E57_1CE5))]
#[allow(clippy::too_many_arguments)]
pub fn extract_thermal_conductivity(
    mut system: PyRefMut<'_, engine::System>,
    positions_m: Vec<f64>,
    dt_s: f64,
    temperature_k: f64,
    friction_per_s: f64,
    berendsen_tau_s: f64,
    equilibration_steps: usize,
    production_steps: usize,
    sample_interval: usize,
    hot_temperature_k: f64,
    cold_temperature_k: f64,
    bath_friction_per_s: f64,
    conductivity_equilibration_steps: usize,
    conductivity_production_steps: usize,
    conductivity_sample_interval: usize,
    slab_atoms: usize,
    cross_section_area_m2: f64,
    seed: u64,
) -> PyResult<ConductivityResult> {
    let config = thermal_config(
        dt_s,
        temperature_k,
        friction_per_s,
        berendsen_tau_s,
        equilibration_steps,
        production_steps,
        sample_interval,
        hot_temperature_k,
        cold_temperature_k,
        bath_friction_per_s,
        conductivity_equilibration_steps,
        conductivity_production_steps,
        conductivity_sample_interval,
        slab_atoms,
        cross_section_area_m2,
        seed,
    );
    let inner = params::extract_thermal_conductivity(&mut system.inner, &positions_m, &config)
        .map_err(util::err)?;
    Ok(ConductivityResult { inner })
}

/// Runs both thermal extractions and combines them into one result.
#[pyfunction]
#[pyo3(signature = (system, positions_m, dt_s=5.0e-15, temperature_k=300.0, friction_per_s=5.0e12, berendsen_tau_s=2.0e-13, equilibration_steps=4_000, production_steps=80_000, sample_interval=20, hot_temperature_k=400.0, cold_temperature_k=200.0, bath_friction_per_s=1.0e13, conductivity_equilibration_steps=10_000, conductivity_production_steps=40_000, conductivity_sample_interval=1, slab_atoms=2, cross_section_area_m2=1.0e-20, seed=0x7E57_1CE5))]
#[allow(clippy::too_many_arguments)]
pub fn extract_thermal(
    mut system: PyRefMut<'_, engine::System>,
    positions_m: Vec<f64>,
    dt_s: f64,
    temperature_k: f64,
    friction_per_s: f64,
    berendsen_tau_s: f64,
    equilibration_steps: usize,
    production_steps: usize,
    sample_interval: usize,
    hot_temperature_k: f64,
    cold_temperature_k: f64,
    bath_friction_per_s: f64,
    conductivity_equilibration_steps: usize,
    conductivity_production_steps: usize,
    conductivity_sample_interval: usize,
    slab_atoms: usize,
    cross_section_area_m2: f64,
    seed: u64,
) -> PyResult<ThermalResult> {
    let config = thermal_config(
        dt_s,
        temperature_k,
        friction_per_s,
        berendsen_tau_s,
        equilibration_steps,
        production_steps,
        sample_interval,
        hot_temperature_k,
        cold_temperature_k,
        bath_friction_per_s,
        conductivity_equilibration_steps,
        conductivity_production_steps,
        conductivity_sample_interval,
        slab_atoms,
        cross_section_area_m2,
        seed,
    );
    let inner =
        params::extract_thermal(&mut system.inner, &positions_m, &config).map_err(util::err)?;
    Ok(ThermalResult { inner })
}
