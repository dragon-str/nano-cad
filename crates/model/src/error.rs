use thiserror::Error;

/// Errors from model validation and from document decoding.
///
/// Library code returns these instead of panicking.
#[derive(Debug, Error, PartialEq, Eq)]
pub enum ModelError {
    #[error("atom index {index} is out of bounds for {len} atoms")]
    AtomIndexOutOfBounds { index: usize, len: usize },
    #[error("bond index {index} is out of bounds for {len} bonds")]
    BondIndexOutOfBounds { index: usize, len: usize },
    #[error("bond endpoint {endpoint} is out of bounds for {atom_count} atoms")]
    BondEndpointOutOfBounds { endpoint: u32, atom_count: usize },
    #[error("bond order {0} is not supported")]
    InvalidBondOrder(u8),
    #[error("bond type tag {0} is not supported")]
    InvalidBondType(u8),
    #[error("element with atomic number {0} is not supported")]
    UnknownElement(u8),
    #[error("payload is not valid UTF-8")]
    InvalidUtf8,
    #[error("unexpected end of buffer at byte {offset}")]
    UnexpectedEof { offset: usize },
    #[error("bad magic bytes {0:?}")]
    BadMagic([u8; 4]),
    #[error("unsupported schema version {found}; expected {expected}")]
    UnsupportedSchemaVersion { found: u32, expected: u32 },
    #[error("buffer has {remaining} trailing bytes after the document")]
    TrailingBytes { remaining: usize },
    #[error("a length or count does not fit the encoded width")]
    IntegerOverflow,
}
