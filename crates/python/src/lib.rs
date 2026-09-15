//! PyO3 bindings to the nanocad Python package.
//!
//! The `python` feature builds the `nanocad._core` extension module through
//! maturin. Without that feature the crate stays a plain library, so
//! `cargo test` never needs a Python interpreter. This follows ADR-0012.
//!
//! # Unsafe policy
//!
//! Every other crate in the workspace uses `#![forbid(unsafe_code)]`. PyO3's
//! macro expansion uses `unsafe` internally, so this crate cannot forbid it
//! when the `python` feature is on. The default build still forbids it. There
//! is no hand-written `unsafe` in this crate.
#![cfg_attr(not(feature = "python"), forbid(unsafe_code))]
#![cfg_attr(feature = "python", allow(unsafe_code))]

/// Returns the crate version string.
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[cfg(feature = "python")]
mod bindings;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_is_not_empty() {
        assert!(!version().is_empty());
    }
}
