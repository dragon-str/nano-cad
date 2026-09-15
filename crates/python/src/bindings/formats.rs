//! File formats: NCZ, MMP, and XYZ, as files and as bytes.

use nanocad_format as format;
use pyo3::prelude::*;

use super::model::Document;
use super::util;

fn document(inner: nanocad_model::Document) -> Document {
    Document { inner }
}

fn extension(path: &str) -> String {
    std::path::Path::new(path)
        .extension()
        .and_then(|ext| ext.to_str())
        .map(str::to_ascii_lowercase)
        .unwrap_or_default()
}

/// Reads an NCZ archive from a file path.
#[pyfunction]
pub fn read_ncz(path: &str) -> PyResult<Document> {
    let bytes = std::fs::read(path).map_err(util::err)?;
    Ok(document(format::from_ncz_bytes(&bytes).map_err(util::err)?))
}

/// Writes a document to an NCZ archive at a file path.
#[pyfunction]
pub fn write_ncz(path: &str, document: PyRef<'_, Document>) -> PyResult<()> {
    let bytes = format::to_ncz_bytes(&document.inner).map_err(util::err)?;
    std::fs::write(path, bytes).map_err(util::err)
}

/// Encodes a document to NCZ bytes.
#[pyfunction]
pub fn ncz_to_bytes(document: PyRef<'_, Document>) -> PyResult<Vec<u8>> {
    format::to_ncz_bytes(&document.inner).map_err(util::err)
}

/// Decodes a document from NCZ bytes.
#[pyfunction]
pub fn ncz_from_bytes(bytes: &[u8]) -> PyResult<Document> {
    Ok(document(format::from_ncz_bytes(bytes).map_err(util::err)?))
}

/// Imports an XYZ text into a one-part document.
#[pyfunction]
pub fn import_xyz(text: &str) -> PyResult<Document> {
    Ok(document(format::import_xyz(text).map_err(util::err)?))
}

/// Exports a document to XYZ text.
#[pyfunction]
pub fn export_xyz(document: PyRef<'_, Document>) -> PyResult<String> {
    format::export_xyz(&document.inner).map_err(util::err)
}

/// Imports an MMP text into a document.
#[pyfunction]
pub fn import_mmp(text: &str) -> PyResult<Document> {
    Ok(document(format::import_mmp(text).map_err(util::err)?))
}

/// Exports a document to MMP text.
#[pyfunction]
pub fn export_mmp(document: PyRef<'_, Document>) -> PyResult<String> {
    format::export_mmp(&document.inner).map_err(util::err)
}

/// Saves a document. The format follows the file extension.
#[pyfunction]
pub fn save(path: &str, document: PyRef<'_, Document>) -> PyResult<()> {
    match extension(path).as_str() {
        "ncz" => write_ncz(path, document),
        "xyz" => std::fs::write(path, export_xyz(document)?.into_bytes()).map_err(util::err),
        "mmp" => std::fs::write(path, export_mmp(document)?.into_bytes()).map_err(util::err),
        other => Err(pyo3::exceptions::PyValueError::new_err(format!(
            "cannot save format {other:?}; expected ncz, xyz, or mmp"
        ))),
    }
}

/// Loads a document. The format follows the file extension.
#[pyfunction]
pub fn load(path: &str) -> PyResult<Document> {
    let bytes = std::fs::read(path).map_err(util::err)?;
    match extension(path).as_str() {
        "ncz" => ncz_from_bytes(&bytes),
        "xyz" => {
            let text = String::from_utf8(bytes).map_err(util::err)?;
            import_xyz(&text)
        }
        "mmp" => {
            let text = String::from_utf8(bytes).map_err(util::err)?;
            import_mmp(&text)
        }
        other => Err(pyo3::exceptions::PyValueError::new_err(format!(
            "cannot load format {other:?}; expected ncz, xyz, or mmp"
        ))),
    }
}

/// The alias of [`save`] used by the agent tool surface.
#[pyfunction]
pub fn export(path: &str, document: PyRef<'_, Document>) -> PyResult<()> {
    save(path, document)
}

/// Loads a document and returns its part names.
#[pyfunction]
pub fn list_parts(path: &str) -> PyResult<Vec<String>> {
    Ok(load(path)?
        .inner
        .parts
        .iter()
        .map(|part| part.name.clone())
        .collect())
}
