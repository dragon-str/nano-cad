//! File formats: NCZ, MMP, XYZ, PDB, and later URDF.
//!
//! This crate converts between the public file formats and the in-memory
//! `nanocad-model` document. It keeps SI units internally and converts at the
//! file boundary.
//!
//! - [`ncz`] is the native, versioned archive format.
//! - [`mmp`] reads and writes NanoEngineer molecular parts.
//! - [`xyz`] reads and writes the plain-text atom list.
//! - [`pdb`] reads and writes Protein Data Bank coordinate files.
#![forbid(unsafe_code)]

mod error;
mod periodic;

pub mod mmp;
pub mod ncz;
pub mod pdb;
pub mod xyz;

pub use error::FormatError;
pub use mmp::{export_mmp, import_mmp};
pub use ncz::{from_ncz_bytes, read_ncz, to_ncz_bytes, write_ncz, NCZ_FORMAT, NCZ_SCHEMA_VERSION};
pub use pdb::{export_pdb, import_pdb};
pub use xyz::{export_xyz, import_xyz};

/// Returns the crate version string.
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_is_not_empty() {
        assert!(!version().is_empty());
    }
}
