//! SI units, `Quantity`, and conversions.

use nanocad_units as units;
use pyo3::prelude::*;

use super::util;

/// A physical value with a presentation unit. The stored value is SI.
#[pyclass(name = "Quantity", skip_from_py_object)]
#[derive(Clone)]
pub struct Quantity {
    pub(crate) inner: units::Quantity,
}

#[pymethods]
impl Quantity {
    #[new]
    fn new(value: f64, unit: &str) -> PyResult<Self> {
        Ok(Self {
            inner: units::Quantity::new(value, util::unit_from_symbol(unit)?),
        })
    }

    /// The value in the presentation unit.
    #[getter]
    fn value(&self) -> f64 {
        self.inner.value_si() / self.inner.unit().factor_to_si()
    }

    /// The value in SI base units.
    #[getter]
    fn value_si(&self) -> f64 {
        self.inner.value_si()
    }

    /// The presentation unit symbol.
    #[getter]
    fn unit(&self) -> String {
        self.inner.unit().symbol().to_owned()
    }

    /// Returns a quantity with the same physical value and a new unit.
    fn to(&self, unit: &str) -> PyResult<Self> {
        Ok(Self {
            inner: self
                .inner
                .convert_to(util::unit_from_symbol(unit)?)
                .map_err(util::err)?,
        })
    }

    /// Reads the value expressed in another unit.
    fn value_in(&self, unit: &str) -> PyResult<f64> {
        self.inner
            .value_in(util::unit_from_symbol(unit)?)
            .map_err(util::err)
    }

    fn __repr__(&self) -> String {
        format!("{}", self.inner)
    }
}

/// Converts a value from one unit to another. Both must share a dimension.
#[pyfunction]
pub fn convert_units(value: f64, from_unit: &str, to_unit: &str) -> PyResult<f64> {
    let quantity = units::Quantity::new(value, util::unit_from_symbol(from_unit)?);
    quantity
        .value_in(util::unit_from_symbol(to_unit)?)
        .map_err(util::err)
}

/// Returns a value in SI base units.
#[pyfunction]
pub fn to_si(value: f64, unit: &str) -> PyResult<f64> {
    Ok(units::Quantity::new(value, util::unit_from_symbol(unit)?).value_si())
}

/// Returns a quantity from a value already in SI base units.
#[pyfunction]
pub fn from_si(value_si: f64, unit: &str) -> PyResult<Quantity> {
    Ok(Quantity {
        inner: units::Quantity::from_si(value_si, util::unit_from_symbol(unit)?),
    })
}

/// Returns every accepted unit symbol.
#[pyfunction]
pub fn unit_symbols() -> Vec<String> {
    util::UNIT_SYMBOLS
        .iter()
        .map(|symbol| (*symbol).to_owned())
        .collect()
}

/// Returns the seven SI dimension exponents of a unit.
#[pyfunction]
pub fn unit_dimension(unit: &str) -> PyResult<Vec<i8>> {
    let dimension = util::unit_from_symbol(unit)?.dimension();
    Ok(vec![
        dimension.length,
        dimension.mass,
        dimension.time,
        dimension.current,
        dimension.temperature,
        dimension.amount,
        dimension.luminous_intensity,
    ])
}
