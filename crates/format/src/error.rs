use thiserror::Error;

use nanocad_model::ModelError;

/// Errors from reading or writing a file format.
///
/// Library code returns these instead of panicking.
#[derive(Debug, Error)]
pub enum FormatError {
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
    #[error("zip error: {0}")]
    Zip(#[from] zip::result::ZipError),
    #[error("json error: {0}")]
    Json(#[from] serde_json::Error),
    #[error("model error: {0}")]
    Model(#[from] ModelError),
    #[error("ncz entry {0} is missing from the archive")]
    MissingEntry(String),
    #[error("ncz header is not valid utf-8")]
    InvalidUtf8,
    #[error("ncz format tag {found:?} is not supported; expected {expected:?}")]
    UnsupportedFormat { found: String, expected: String },
    #[error("ncz schema version {found} is not supported; expected {expected}")]
    UnsupportedSchemaVersion { found: u32, expected: u32 },
    #[error("ncz array {name} has {found} values; expected {expected}")]
    BadArrayLength {
        name: String,
        found: usize,
        expected: usize,
    },
    #[error("element with atomic number {0} is not supported")]
    UnknownElement(u8),
    #[error("element symbol {0:?} is not known")]
    UnknownElementSymbol(String),
    #[error("an mmp line is malformed at line {line}: {message}")]
    MmpParse { line: usize, message: String },
    #[error("an mmp file names no atoms or is empty")]
    MmpEmpty,
    #[error("an xyz line is malformed at line {line}: {message}")]
    XyzParse { line: usize, message: String },
    #[error("an xyz count line does not match the number of atom lines")]
    XyzCountMismatch,
    #[error("a count or length does not fit the encoded width")]
    IntegerOverflow,
}
