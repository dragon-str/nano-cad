//! The parameter store: quantities, provenance, records, and libraries.

use nanocad_params as params;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

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
