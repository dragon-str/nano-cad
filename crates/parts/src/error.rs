use nanocad_model::ModelError;
use thiserror::Error;

/// Errors from part generation and parameter validation.
///
/// Generator code returns these instead of panicking. Invalid size input is a
/// value error, never a panic.
#[derive(Debug, Error, PartialEq)]
pub enum PartError {
    #[error("unknown parameter '{0}'")]
    UnknownParameter(String),
    #[error("missing parameter '{0}'")]
    MissingParameter(String),
    #[error("parameter '{name}' value {value} is outside [{min}, {max}]")]
    ParameterOutOfRange {
        name: String,
        value: f64,
        min: f64,
        max: f64,
    },
    #[error("parameter '{name}' must be an integer, got {value}")]
    NonIntegerParameter { name: String, value: f64 },
    #[error("unknown generator id '{0}'")]
    UnknownGenerator(String),
    #[error("invalid geometry: {0}")]
    InvalidGeometry(String),
    #[error("tooth count {teeth} is below the undercut limit {min_teeth}")]
    ToothCountBelowUndercut { teeth: usize, min_teeth: usize },
    #[error("bore radius {bore_m} m is not smaller than the gear root radius {root_m} m")]
    BoreTooLarge { bore_m: f64, root_m: f64 },
    #[error(
        "planetary constraint violated: ring teeth {ring_teeth} != sun teeth {sun_teeth} \
         + 2 * planet teeth {planet_teeth}"
    )]
    PlanetaryConstraintViolated {
        sun_teeth: usize,
        planet_teeth: usize,
        ring_teeth: usize,
    },
    #[error(
        "planet count {planet_count} does not divide the sun-plus-ring tooth total \
         {total_teeth}, so the planets cannot space equally"
    )]
    PlanetSpacingNotPossible {
        planet_count: usize,
        total_teeth: usize,
    },
    #[error(
        "planets overlap: adjacent centre distance {center_m} m is below the clearance \
         {clearance_m} m"
    )]
    PlanetsOverlap { center_m: f64, clearance_m: f64 },
    #[error(transparent)]
    Model(#[from] ModelError),
}
