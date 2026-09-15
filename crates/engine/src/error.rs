use thiserror::Error;

/// Errors from engine setup and from force-term evaluation.
///
/// Library code returns these instead of panicking.
#[derive(Debug, Error, PartialEq)]
pub enum EngineError {
    #[error("cutoff must be a positive, finite length, got {cutoff_m} m")]
    NonPositiveCutoff { cutoff_m: f64 },
    #[error("a water box must hold at least one molecule, got {molecule_count}")]
    InvalidMoleculeCount { molecule_count: usize },
    #[error("skin must be a non-negative, finite length, got {skin_m} m")]
    InvalidSkin { skin_m: f64 },
    #[error("switch-on distance {switch_on_m} m must be at least zero and below the cutoff {cutoff_m} m")]
    InvalidSwitchingRange { cutoff_m: f64, switch_on_m: f64 },
    #[error("periodic box length {index} must be a non-negative, finite length, got {length_m} m")]
    InvalidPeriodicBox { index: usize, length_m: f64 },
    #[error("cutoff {cutoff_m} m exceeds half the periodic box length {length_m} m")]
    CutoffExceedsHalfBox { cutoff_m: f64, length_m: f64 },
    #[error("position buffer length {len} is not three times the atom count")]
    PositionBufferSizeMismatch { len: usize, atom_count: usize },
    #[error("position {index} is not finite")]
    NonFinitePosition { index: usize },
    #[error("bond {bond} has an endpoint outside the atom count")]
    BondEndpointOutOfBounds { bond: usize, atom_count: usize },
    #[error("bond parameters for bond {bond} are invalid")]
    InvalidBondParameters { bond: usize },
    #[error("bond parameter count {provided} does not match {expected} bonds")]
    BondParameterCountMismatch { provided: usize, expected: usize },
    #[error("bond {bond} has coincident atoms, so its direction is undefined")]
    CoincidentBondAtoms { bond: usize },
    #[error("angle parameters for angle {angle} are invalid")]
    InvalidAngleParameters { angle: usize },
    #[error("angle parameter count {provided} does not match {expected} angles")]
    AngleParameterCountMismatch { provided: usize, expected: usize },
    #[error("angle {angle} is degenerate, so its bend angle is undefined")]
    AngleGeometryUndefined { angle: usize },
    #[error("torsion parameters for torsion {torsion} are invalid")]
    InvalidTorsionParameters { torsion: usize },
    #[error("torsion parameter count {provided} does not match {expected} torsions")]
    TorsionParameterCountMismatch { provided: usize, expected: usize },
    #[error("torsion {torsion} is degenerate, so its dihedral angle is undefined")]
    TorsionGeometryUndefined { torsion: usize },
    #[error("out-of-plane parameters for improper {improper} are invalid")]
    InvalidOutOfPlaneParameters { improper: usize },
    #[error("out-of-plane parameter count {provided} does not match {expected} impropers")]
    OutOfPlaneParameterCountMismatch { provided: usize, expected: usize },
    #[error("out-of-plane improper {improper} is degenerate, so its angle is undefined")]
    OutOfPlaneGeometryUndefined { improper: usize },
    #[error("van der Waals parameters for atom {atom} are invalid")]
    InvalidVdwParameters { atom: usize },
    #[error("van der Waals parameter count {provided} does not match {expected} atoms")]
    VdwParameterCountMismatch { provided: usize, expected: usize },
    #[error("charge for atom {atom} is not finite")]
    InvalidCharge { atom: usize },
    #[error("charge count {provided} does not match {expected} atoms")]
    ChargeCountMismatch { provided: usize, expected: usize },
    #[error("atoms {i} and {j} are coincident, so their non-bonded energy is undefined")]
    CoincidentNonbondedAtoms { i: u32, j: u32 },
    #[error("atom count {atom_count} exceeds the u32 pair-list index width")]
    TooManyAtoms { atom_count: usize },
    #[error("the grid cell count for these positions is too large")]
    GridTooLarge,
    #[error("term {term} covers {provided} atoms, but the system has {expected}")]
    TermAtomCountMismatch {
        term: &'static str,
        provided: usize,
        expected: usize,
    },
    #[error("pair list length {len} is not even")]
    PairListSizeMismatch { len: usize },
    #[error("pair index {index} at entry {entry} is outside the atom count {atom_count}")]
    PairIndexOutOfBounds {
        entry: usize,
        index: u32,
        atom_count: usize,
    },
    #[error("buffer length {len} does not match the expected length {expected}")]
    BufferSizeMismatch { len: usize, expected: usize },
    #[error("velocity {index} is not finite")]
    NonFiniteVelocity { index: usize },
    #[error("external force {index} is not finite")]
    NonFiniteExternalForce { index: usize },
    #[error("no masses are set on the system")]
    MassesNotSet,
    #[error("mass count {provided} does not match {expected} atoms")]
    MassCountMismatch { provided: usize, expected: usize },
    #[error("mass for atom {atom} is not finite and positive")]
    InvalidMass { atom: usize },
    #[error("thread count must be at least one, got {thread_count}")]
    InvalidThreadCount { thread_count: usize },
    #[error("a force worker thread panicked")]
    WorkerPanicked,
    #[error("temperature must be a positive, finite value, got {temperature_k} K")]
    InvalidTemperature { temperature_k: f64 },
    #[error("friction must be a non-negative, finite value, got {friction_per_s} / s")]
    InvalidFriction { friction_per_s: f64 },
    #[error("coupling time constant must be a positive, finite duration, got {tau_s} s")]
    InvalidCouplingTime { tau_s: f64 },
    #[error("time step must be a positive, finite duration, got {dt_s} s")]
    NonPositiveTimestep { dt_s: f64 },
    #[error("gradient tolerance must be finite and positive, got {tolerance_n} N")]
    InvalidGradientTolerance { tolerance_n: f64 },
    #[error("the minimizer iteration cap must be at least one")]
    InvalidMaxIterations,
    #[error("the initial line-search step must be finite and positive, got {step_m} m")]
    InvalidInitialStep { step_m: f64 },
    #[error("the line search did not find a lower point at iteration {iteration}")]
    LineSearchFailed { iteration: usize },
}
