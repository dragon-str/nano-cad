//! Document model: atoms, bonds, topology, parts.

use nanocad_model as model;
use pyo3::prelude::*;
use std::collections::BTreeMap;

use super::util;

fn element_from_atomic_number(atomic_number: u8) -> PyResult<model::Element> {
    model::Element::from_atomic_number(atomic_number).ok_or_else(|| {
        pyo3::exceptions::PyValueError::new_err(format!(
            "atomic number {atomic_number} is out of range"
        ))
    })
}

/// One atom as a value.
#[pyclass(name = "Atom", skip_from_py_object)]
#[derive(Clone)]
pub struct Atom {
    pub(crate) inner: model::Atom,
}

#[pymethods]
impl Atom {
    #[new]
    fn new(element: u8, position_m: Vec<f64>, charge_c: f64, atom_type: String) -> PyResult<Self> {
        Ok(Self {
            inner: model::Atom::new(
                element_from_atomic_number(element)?,
                util::vec3(position_m)?,
                charge_c,
                atom_type,
            ),
        })
    }

    /// The atomic number.
    #[getter]
    fn element(&self) -> u8 {
        self.inner.element.atomic_number()
    }

    /// The element symbol, or `None` when the element is a dummy.
    #[getter]
    fn symbol(&self) -> Option<String> {
        self.inner.element.symbol().map(str::to_owned)
    }

    /// The position in SI metres.
    #[getter]
    fn position_m(&self) -> (f64, f64, f64) {
        let p = self.inner.position_m;
        (p[0], p[1], p[2])
    }

    /// The partial charge in SI coulombs.
    #[getter]
    fn charge_c(&self) -> f64 {
        self.inner.charge_c
    }

    /// The force-field atom type name.
    #[getter]
    fn atom_type(&self) -> String {
        self.inner.atom_type.clone()
    }
}

/// Atoms and bonds as struct-of-arrays, addressed by index.
#[pyclass(name = "Topology", skip_from_py_object)]
#[derive(Clone)]
pub struct Topology {
    pub(crate) inner: model::Topology,
}

#[pymethods]
impl Topology {
    #[new]
    fn new() -> Self {
        Self {
            inner: model::Topology::new(),
        }
    }

    #[getter]
    fn atom_count(&self) -> usize {
        self.inner.atom_count()
    }

    #[getter]
    fn bond_count(&self) -> usize {
        self.inner.bond_count()
    }

    /// Appends an atom and returns its index.
    fn add_atom(
        &mut self,
        element: u8,
        position_m: Vec<f64>,
        charge_c: f64,
        atom_type: String,
    ) -> PyResult<usize> {
        Ok(self.inner.add_atom(model::Atom::new(
            element_from_atomic_number(element)?,
            util::vec3(position_m)?,
            charge_c,
            atom_type,
        )) as usize)
    }

    /// Appends a bond and returns its index.
    fn add_bond(&mut self, u: u32, v: u32, order: u8, bond_type: &str) -> PyResult<usize> {
        Ok(self
            .inner
            .add_bond(model::Bond::new(
                u,
                v,
                order,
                util::bond_type_from_str(bond_type)?,
            ))
            .map_err(util::err)? as usize)
    }

    /// Returns the atom at an index, or `None`.
    fn atom(&self, index: usize) -> Option<Atom> {
        self.inner.atom(index).map(|inner| Atom { inner })
    }

    /// Returns the position at an index in SI metres, or `None`.
    fn position_m(&self, index: usize) -> Option<(f64, f64, f64)> {
        self.inner
            .position_m(index)
            .map(|position_m| (position_m[0], position_m[1], position_m[2]))
    }

    /// Returns the charge at an index in SI coulombs, or `None`.
    fn charge_c(&self, index: usize) -> Option<f64> {
        self.inner.charge_c(index)
    }

    /// Returns the atomic number at an index, or `None`.
    fn element(&self, index: usize) -> Option<u8> {
        self.inner
            .element(index)
            .map(|element| element.atomic_number())
    }

    /// Returns the atom type at an index, or `None`.
    fn atom_type(&self, index: usize) -> Option<String> {
        self.inner.atom_type(index).map(str::to_owned)
    }

    /// Returns the bond at an index as `(u, v, order, bond_type)`, or `None`.
    fn bond(&self, index: usize) -> Option<(u32, u32, u8, String)> {
        self.inner.bond(index).map(|bond| {
            (
                bond.u,
                bond.v,
                bond.order,
                util::bond_type_name(bond.bond_type).to_owned(),
            )
        })
    }

    /// Returns every atom in index order.
    fn atoms(&self) -> Vec<Atom> {
        self.inner.atoms().map(|inner| Atom { inner }).collect()
    }

    /// Returns every bond in index order.
    fn bonds(&self) -> Vec<(u32, u32, u8, String)> {
        self.inner
            .bonds()
            .map(|bond| {
                (
                    bond.u,
                    bond.v,
                    bond.order,
                    util::bond_type_name(bond.bond_type).to_owned(),
                )
            })
            .collect()
    }

    /// Returns the flat position buffer.
    fn positions_m(&self) -> Vec<f64> {
        self.inner.positions_m().to_vec()
    }

    /// Returns the flat charge buffer.
    fn charges_c(&self) -> Vec<f64> {
        self.inner.charges_c().to_vec()
    }

    /// Returns the atomic numbers of every atom.
    fn elements(&self) -> Vec<u8> {
        self.inner
            .elements()
            .iter()
            .map(|element| element.atomic_number())
            .collect()
    }
}

/// A named collection of atoms and bonds, plus a material.
#[pyclass(name = "Part", skip_from_py_object)]
#[derive(Clone)]
pub struct Part {
    pub(crate) inner: model::Part,
}

#[pymethods]
impl Part {
    #[new]
    #[pyo3(signature = (name, topology=None))]
    fn new(name: String, topology: Option<PyRef<'_, Topology>>) -> Self {
        let inner = match topology {
            Some(topology) => topology.inner.clone(),
            None => model::Topology::new(),
        };
        Self {
            inner: model::Part::new(name, inner),
        }
    }

    #[getter]
    fn name(&self) -> String {
        self.inner.name.clone()
    }

    #[setter]
    fn set_name(&mut self, name: String) {
        self.inner.name = name;
    }

    #[getter]
    fn material(&self) -> String {
        self.inner.material.clone()
    }

    #[setter]
    fn set_material(&mut self, material: String) {
        self.inner.material = material;
    }

    #[getter]
    fn topology(&self) -> Topology {
        Topology {
            inner: self.inner.topology.clone(),
        }
    }

    #[getter]
    fn metadata(&self) -> BTreeMap<String, String> {
        self.inner.metadata.clone()
    }

    #[getter]
    fn atom_count(&self) -> usize {
        self.inner.atom_count()
    }

    #[getter]
    fn bond_count(&self) -> usize {
        self.inner.bond_count()
    }

    /// Inserts or replaces one metadata entry.
    fn set_metadata(&mut self, key: String, value: String) {
        self.inner.metadata.insert(key, value);
    }

    /// Replaces every atom position from a flat buffer, three values per atom.
    ///
    /// The buffer length must be three times the atom count.
    fn set_positions_m(&mut self, positions_m: Vec<f64>) -> PyResult<()> {
        let atom_count = self.inner.atom_count();
        if positions_m.len() != atom_count * 3 {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "expected {} positions, got {}",
                atom_count * 3,
                positions_m.len()
            )));
        }
        for atom in 0..atom_count {
            self.inner
                .topology
                .set_position_m(
                    atom,
                    [
                        positions_m[3 * atom],
                        positions_m[3 * atom + 1],
                        positions_m[3 * atom + 2],
                    ],
                )
                .map_err(util::err)?;
        }
        Ok(())
    }

    /// Encodes the part to the private little-endian buffer.
    fn to_bytes(&self) -> PyResult<Vec<u8>> {
        self.inner.to_bytes().map_err(util::err)
    }

    /// Decodes a part from the private little-endian buffer.
    #[staticmethod]
    fn from_bytes(bytes: &[u8]) -> PyResult<Self> {
        Ok(Self {
            inner: model::Part::from_bytes(bytes).map_err(util::err)?,
        })
    }
}

/// A versioned, serializable design graph.
#[pyclass(name = "Document", skip_from_py_object)]
#[derive(Clone)]
pub struct Document {
    pub(crate) inner: model::Document,
}

#[pymethods]
impl Document {
    #[new]
    fn new(name: String) -> Self {
        Self {
            inner: model::Document::new(name),
        }
    }

    #[getter]
    fn name(&self) -> String {
        self.inner.name.clone()
    }

    #[setter]
    fn set_name(&mut self, name: String) {
        self.inner.name = name;
    }

    #[getter]
    fn schema_version(&self) -> u32 {
        self.inner.schema_version
    }

    #[getter]
    fn metadata(&self) -> BTreeMap<String, String> {
        self.inner.metadata.clone()
    }

    #[getter]
    fn part_count(&self) -> usize {
        self.inner.part_count()
    }

    /// Appends a part and returns its index.
    fn add_part(&mut self, part: PyRef<'_, Part>) -> usize {
        self.inner.add_part(part.inner.clone())
    }

    /// Returns a part by index, or `None`.
    fn part(&self, index: usize) -> Option<Part> {
        self.inner.part(index).map(|inner| Part {
            inner: inner.clone(),
        })
    }

    /// Returns the first part with a name, or `None`.
    fn part_by_name(&self, name: &str) -> Option<Part> {
        self.inner.part_by_name(name).map(|inner| Part {
            inner: inner.clone(),
        })
    }

    /// Returns every part in index order.
    fn parts(&self) -> Vec<Part> {
        self.inner
            .parts
            .iter()
            .map(|inner| Part {
                inner: inner.clone(),
            })
            .collect()
    }

    /// Returns the names of every part.
    fn part_names(&self) -> Vec<String> {
        self.inner
            .parts
            .iter()
            .map(|part| part.name.clone())
            .collect()
    }

    /// Inserts or replaces one metadata entry.
    fn set_metadata(&mut self, key: String, value: String) {
        self.inner.metadata.insert(key, value);
    }

    /// Encodes the document to the private little-endian buffer.
    fn to_bytes(&self) -> PyResult<Vec<u8>> {
        self.inner.to_bytes().map_err(util::err)
    }

    /// Decodes a document from the private little-endian buffer.
    #[staticmethod]
    fn from_bytes(bytes: &[u8]) -> PyResult<Self> {
        Ok(Self {
            inner: model::Document::from_bytes(bytes).map_err(util::err)?,
        })
    }
}

/// Returns the document schema version.
#[pyfunction]
pub fn document_schema_version() -> u32 {
    model::SCHEMA_VERSION
}
