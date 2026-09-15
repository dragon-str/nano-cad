use thiserror::Error;

/// Errors from parameter-record decoding and construction.
///
/// Library code returns these instead of panicking.
#[derive(Debug, Error)]
pub enum ParamError {
    /// A unit symbol is not in the unit registry.
    #[error("unknown unit symbol {symbol:?}")]
    UnknownUnit { symbol: String },
    /// The record names a schema that this crate does not implement.
    #[error("unsupported schema {found:?}; expected {expected:?}")]
    UnsupportedSchema { found: String, expected: String },
    /// The record names a schema version that this crate does not implement.
    #[error("unsupported version {found}; expected {expected}")]
    UnsupportedVersion { found: u32, expected: u32 },
    /// The JSON payload is malformed or misses a required field.
    #[error("invalid parameter record JSON: {0}")]
    Json(#[from] serde_json::Error),
    /// The store already holds this part id at this version.
    #[error("duplicate part record {part_id:?} at version {version}")]
    DuplicateRecord { part_id: String, version: u32 },
    /// A stored geometry hash does not match the expected hash.
    #[error(
        "geometry reference mismatch for {part_id:?} version {version}: expected {expected:?}, found {found:?}"
    )]
    GeometryMismatch {
        part_id: String,
        version: u32,
        expected: String,
        found: String,
    },
    /// The library holds no record for this part id and version.
    #[error("unknown part {part_id:?} at version {version}")]
    UnknownPart { part_id: String, version: u32 },
    /// The library names a schema that this crate does not implement.
    #[error("unsupported library schema {found:?}; expected {expected:?}")]
    UnsupportedLibrarySchema { found: String, expected: String },
    /// The library names a schema version that this crate does not implement.
    #[error("unsupported library version {found}; expected {expected}")]
    UnsupportedLibraryVersion { found: u32, expected: u32 },
    /// A loaded library is internally inconsistent, for example a key does not
    /// match the record part id.
    #[error("inconsistent parameter library: {reason}")]
    LibraryIntegrity { reason: String },
    /// A file could not be read or written.
    #[error("parameter library i/o error: {0}")]
    Io(#[from] std::io::Error),
    /// The engine failed during an extraction run.
    #[error("engine error during parameter extraction: {0}")]
    Engine(#[from] nanocad_engine::EngineError),
    /// An extraction configuration or its inputs are not usable.
    #[error("invalid parameter extraction: {reason}")]
    InvalidExtraction { reason: String },
    /// A thermal extraction produced fewer samples than the estimator needs.
    #[error("thermal extraction needs at least {required} samples, got {measured}")]
    ThermalSamplesTooFew { required: usize, measured: usize },
    /// The two thermal bath slabs do not fit inside the atom count.
    #[error("thermal convection needs {requested} slab atoms but the system has {available}")]
    ThermalSlabsTooLarge { requested: usize, available: usize },
    /// A thermal observable is not finite, so the estimate is undefined.
    #[error("thermal observable {observable} is not finite")]
    NonFiniteThermalObservable { observable: &'static str },
    /// The two baths did not hold a positive temperature difference, so the
    /// conductivity is undefined.
    #[error("the thermal gradient did not establish: hot {hot_k} K, cold {cold_k} K")]
    ThermalGradientNotEstablished { hot_k: f64, cold_k: f64 },
}
