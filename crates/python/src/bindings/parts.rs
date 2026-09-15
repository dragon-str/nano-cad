//! Parametric part generators and part validation.

use nanocad_parts as parts;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::HashMap;

use super::model::{Document, Part};
use super::util;

/// A default reference bond length for validation, in metres. It is the
/// diamondoid carbon-carbon bond.
const DEFAULT_REFERENCE_M: f64 = 1.544e-10;

fn generators() -> [&'static dyn parts::PartGenerator; 6] {
    [
        &parts::SpurGearGenerator,
        &parts::GearProfileGenerator,
        &parts::PlanetaryGenerator,
        &parts::NanotubeGenerator,
        &parts::DiamondGenerator,
        &parts::GraphiteGenerator,
    ]
}

fn generator(id: &str) -> Option<&'static dyn parts::PartGenerator> {
    generators()
        .into_iter()
        .find(|generator| generator.id() == id)
}

fn parameter_set(specs: &HashMap<String, f64>) -> parts::ParameterSet {
    let mut set = parts::ParameterSet::new();
    for (name, value) in specs {
        set.set(name.clone(), *value);
    }
    set
}

/// The specification of one generator parameter.
#[pyclass(name = "ParameterInfo", skip_from_py_object)]
#[derive(Clone)]
pub struct ParameterInfo {
    name: String,
    unit: Option<String>,
    default_value: f64,
    min: f64,
    max: f64,
    integer: bool,
    description: String,
}

#[pymethods]
impl ParameterInfo {
    #[getter]
    fn name(&self) -> String {
        self.name.clone()
    }

    /// The presentation unit symbol, or `None` when the parameter is
    /// dimensionless.
    #[getter]
    fn unit(&self) -> Option<String> {
        self.unit.clone()
    }

    /// The default value, in SI base units or dimensionless.
    #[getter]
    fn default_value(&self) -> f64 {
        self.default_value
    }

    #[getter]
    fn min(&self) -> f64 {
        self.min
    }

    #[getter]
    fn max(&self) -> f64 {
        self.max
    }

    /// Whether the value must be a whole number.
    #[getter]
    fn integer(&self) -> bool {
        self.integer
    }

    #[getter]
    fn description(&self) -> String {
        self.description.clone()
    }

    fn __repr__(&self) -> String {
        match &self.unit {
            Some(unit) => format!("ParameterInfo({}, unit={unit})", self.name),
            None => format!("ParameterInfo({})", self.name),
        }
    }
}

fn parameter_info(spec: &parts::ParameterSpec) -> ParameterInfo {
    ParameterInfo {
        name: spec.name.to_owned(),
        unit: spec.unit.map(|unit| unit.symbol().to_owned()),
        default_value: spec.default,
        min: spec.min,
        max: spec.max,
        integer: spec.integer,
        description: spec.description.to_owned(),
    }
}

/// A registered part generator and its parameters.
#[pyclass(name = "GeneratorInfo", skip_from_py_object)]
#[derive(Clone)]
pub struct GeneratorInfo {
    id: String,
    name: String,
    parameters: Vec<ParameterInfo>,
}

#[pymethods]
impl GeneratorInfo {
    #[getter]
    fn id(&self) -> String {
        self.id.clone()
    }

    #[getter]
    fn name(&self) -> String {
        self.name.clone()
    }

    #[getter]
    fn parameters(&self) -> Vec<ParameterInfo> {
        self.parameters.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "GeneratorInfo({}, {} parameters)",
            self.id,
            self.parameters.len()
        )
    }
}

/// Returns every registered generator with its parameters.
#[pyfunction]
pub fn list_generators() -> Vec<GeneratorInfo> {
    generators()
        .into_iter()
        .map(|generator| GeneratorInfo {
            id: generator.id().to_owned(),
            name: generator.name().to_owned(),
            parameters: generator.parameters().iter().map(parameter_info).collect(),
        })
        .collect()
}

/// Returns the parameters of one generator, or `None` for an unknown id.
#[pyfunction]
pub fn generator_parameters(id: &str) -> Option<Vec<ParameterInfo>> {
    generator(id).map(|generator| generator.parameters().iter().map(parameter_info).collect())
}

/// Generates a part from a generator id and a parameter map.
#[pyfunction]
pub fn generate_part(name: &str, specs: HashMap<String, f64>) -> PyResult<Part> {
    let generator = generator(name).ok_or_else(|| {
        pyo3::exceptions::PyValueError::new_err(format!("unknown generator {name:?}"))
    })?;
    let inner = generator
        .generate(&parameter_set(&specs))
        .map_err(util::err)?;
    Ok(Part { inner })
}

/// Generates a gear from the agent-facing gear registry.
#[pyfunction]
pub fn generate_gear(name: &str, specs: HashMap<String, f64>) -> PyResult<Part> {
    let inner = parts::generate_gear(name, &parameter_set(&specs)).map_err(util::err)?;
    Ok(Part { inner })
}

/// A resolved planetary design, with the derived geometry.
#[pyclass(name = "PlanetaryGeometry", skip_from_py_object)]
#[derive(Clone)]
pub struct PlanetaryGeometry {
    inner: parts::PlanetaryDesign,
}

#[pymethods]
impl PlanetaryGeometry {
    #[getter]
    fn module_m(&self) -> f64 {
        self.inner.module_m()
    }

    #[getter]
    fn sun_teeth(&self) -> usize {
        self.inner.sun_teeth()
    }

    #[getter]
    fn planet_teeth(&self) -> usize {
        self.inner.planet_teeth()
    }

    #[getter]
    fn planet_count(&self) -> usize {
        self.inner.planet_count()
    }

    #[getter]
    fn pressure_angle_rad(&self) -> f64 {
        self.inner.pressure_angle_rad()
    }

    #[getter]
    fn ring_teeth(&self) -> usize {
        self.inner.ring_teeth()
    }

    /// Whether the assembly constraint holds.
    #[getter]
    fn constraint_holds(&self) -> bool {
        self.inner.constraint_holds()
    }

    #[getter]
    fn sun_pitch_radius_m(&self) -> f64 {
        self.inner.sun_pitch_radius_m()
    }

    #[getter]
    fn planet_pitch_radius_m(&self) -> f64 {
        self.inner.planet_pitch_radius_m()
    }

    #[getter]
    fn ring_pitch_radius_m(&self) -> f64 {
        self.inner.ring_pitch_radius_m()
    }

    #[getter]
    fn carrier_radius_m(&self) -> f64 {
        self.inner.carrier_radius_m()
    }

    #[getter]
    fn planet_outer_radius_m(&self) -> f64 {
        self.inner.planet_outer_radius_m()
    }

    /// The carrier-to-sun speed ratio.
    #[getter]
    fn gear_ratio(&self) -> f64 {
        self.inner.gear_ratio()
    }

    #[getter]
    fn minimum_planet_spacing_m(&self) -> f64 {
        self.inner.minimum_planet_spacing_m()
    }

    #[getter]
    fn planet_clearance_m(&self) -> f64 {
        self.inner.planet_clearance_m()
    }

    /// Returns the angular position of one planet, in radians.
    fn planet_angle_rad(&self, index: usize) -> PyResult<f64> {
        if index >= self.inner.planet_count() {
            return Err(pyo3::exceptions::PyIndexError::new_err(format!(
                "planet index {index} is out of range"
            )));
        }
        Ok(self.inner.planet_angle_rad(index))
    }

    /// Returns the plane position of one planet, in metres.
    fn planet_center_m(&self, index: usize) -> PyResult<(f64, f64)> {
        if index >= self.inner.planet_count() {
            return Err(pyo3::exceptions::PyIndexError::new_err(format!(
                "planet index {index} is out of range"
            )));
        }
        let center = self.inner.planet_center_m(index);
        Ok((center[0], center[1]))
    }
}

/// Resolves a planetary design from its parameters.
///
/// Every value is in SI: `module_m` in metres, `pressure_angle_rad` in radians.
#[pyfunction]
#[pyo3(signature = (module_m, sun_teeth, planet_teeth, planet_count, pressure_angle_rad))]
pub fn planetary_geometry(
    module_m: f64,
    sun_teeth: usize,
    planet_teeth: usize,
    planet_count: usize,
    pressure_angle_rad: f64,
) -> PyResult<PlanetaryGeometry> {
    Ok(PlanetaryGeometry {
        inner: parts::PlanetaryDesign::new(
            module_m,
            sun_teeth,
            planet_teeth,
            planet_count,
            pressure_angle_rad,
        )
        .map_err(util::err)?,
    })
}

/// Validates a part. Returns every violation as a readable string.
#[pyfunction]
#[pyo3(signature = (part, reference_m=None, strain_tolerance_relative=parts::DEFAULT_STRAIN_TOLERANCE_RELATIVE))]
pub fn validate_part(
    part: PyRef<'_, Part>,
    reference_m: Option<f64>,
    strain_tolerance_relative: f64,
) -> Vec<String> {
    let report = parts::validate_part(
        &part.inner,
        reference_m.unwrap_or(DEFAULT_REFERENCE_M),
        strain_tolerance_relative,
    );
    report
        .violations()
        .iter()
        .map(|violation| format!("{violation:?}"))
        .collect()
}

/// Returns the valence of an element by atomic number, or `None`.
#[pyfunction]
pub fn element_valence(atomic_number: u8) -> Option<u32> {
    nanocad_model::Element::from_atomic_number(atomic_number).and_then(parts::element_valence)
}

/// Returns the full parameter schema of a generator as a Python dict.
#[pyfunction]
pub fn generator_schema(py: Python<'_>, id: &str) -> PyResult<Option<Py<PyDict>>> {
    let Some(generator) = generator(id) else {
        return Ok(None);
    };
    let dict = PyDict::new(py);
    dict.set_item("id", generator.id())?;
    dict.set_item("name", generator.name())?;
    let parameters = PyDict::new(py);
    for spec in generator.parameters() {
        parameters.set_item(spec.name, parameter_info(spec))?;
    }
    dict.set_item("parameters", parameters)?;
    Ok(Some(dict.unbind()))
}

/// Alias of [`generate_part`] used by the agent tool surface.
#[pyfunction]
#[pyo3(signature = (name, specs=None))]
pub fn create_part(name: &str, specs: Option<HashMap<String, f64>>) -> PyResult<Part> {
    generate_part(name, specs.unwrap_or_default())
}

/// Generates a nanotube part with the default parameters.
#[pyfunction]
#[pyo3(signature = (specs=None))]
pub fn generate_nanotube(specs: Option<HashMap<String, f64>>) -> PyResult<Part> {
    generate_part("nanotube", specs.unwrap_or_default())
}

/// Generates a lattice part. `kind` is `"diamond"` or `"graphite"`.
#[pyfunction]
#[pyo3(signature = (kind, specs=None))]
pub fn generate_lattice(kind: &str, specs: Option<HashMap<String, f64>>) -> PyResult<Part> {
    generate_part(kind, specs.unwrap_or_default())
}

/// Assembles parts into a new document named `"assembly"`.
#[pyfunction]
pub fn assemble(parts: Vec<PyRef<'_, Part>>) -> Document {
    let mut inner = nanocad_model::Document::new("assembly");
    for part in parts {
        inner.add_part(part.inner.clone());
    }
    Document { inner }
}

/// Measures the geometry of a part. Returns counts, the bounding box, the
/// centroid, and the total charge. No property is estimated.
#[pyfunction]
pub fn measure(py: Python<'_>, part: PyRef<'_, Part>) -> PyResult<Py<PyDict>> {
    let topology = &part.inner.topology;
    let positions_m = topology.positions_m();
    let atom_count = topology.atom_count();

    let mut min_m = [f64::INFINITY; 3];
    let mut max_m = [f64::NEG_INFINITY; 3];
    let mut centroid_m = [0.0_f64; 3];
    for atom in 0..atom_count {
        for (axis, value) in centroid_m.iter_mut().enumerate() {
            let coordinate_m = positions_m[3 * atom + axis];
            min_m[axis] = min_m[axis].min(coordinate_m);
            max_m[axis] = max_m[axis].max(coordinate_m);
            *value += coordinate_m;
        }
    }
    let size_m = if atom_count == 0 {
        min_m = [0.0; 3];
        max_m = [0.0; 3];
        [0.0; 3]
    } else {
        for value in &mut centroid_m {
            *value /= atom_count as f64;
        }
        [
            max_m[0] - min_m[0],
            max_m[1] - min_m[1],
            max_m[2] - min_m[2],
        ]
    };
    let total_charge_c: f64 = topology.charges_c().iter().sum();

    let dict = PyDict::new(py);
    dict.set_item("name", part.inner.name.clone())?;
    dict.set_item("atom_count", atom_count)?;
    dict.set_item("bond_count", part.inner.bond_count())?;
    dict.set_item("bounding_box_min_m", min_m)?;
    dict.set_item("bounding_box_max_m", max_m)?;
    dict.set_item("bounding_box_size_m", size_m)?;
    dict.set_item("centroid_m", centroid_m)?;
    dict.set_item("total_charge_c", total_charge_c)?;
    Ok(dict.unbind())
}
