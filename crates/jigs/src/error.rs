use nanocad_units::Unit;
use thiserror::Error;

/// Errors from jig setup and from jig evaluation.
///
/// Library code returns these instead of panicking.
#[derive(Debug, Error, PartialEq)]
pub enum JigError {
    #[error("position buffer length {len} is not three times the atom count {atom_count}")]
    PositionBufferSizeMismatch { len: usize, atom_count: usize },
    #[error("position {index} is not finite")]
    NonFinitePosition { index: usize },
    #[error("atom index {atom} is outside the atom count {atom_count}")]
    AtomIndexOutOfBounds { atom: usize, atom_count: usize },
    #[error("atom {atom} has a hold position that is not finite")]
    NonFiniteHoldPosition { atom: usize },
    #[error("spring {spring} has a fixed point that is not finite")]
    NonFiniteFixedPoint { spring: usize },
    #[error("stiffness {k_n_per_m} N/m must be finite and non-negative")]
    InvalidStiffness { k_n_per_m: f64 },
    #[error("rest length {rest_length_m} m must be finite and positive")]
    InvalidRestLength { rest_length_m: f64 },
    #[error("spring {spring} connects atom {atom} to itself")]
    SelfConnectedSpring { spring: usize, atom: usize },
    #[error("spring {spring} has coincident atoms, so its direction is undefined")]
    CoincidentSpringAtoms { spring: usize },
    #[error("motor axis {axis:?} must be a finite, non-zero vector")]
    InvalidAxis { axis: [f64; 3] },
    #[error("motor angular velocity {angular_velocity_rad_per_s} rad/s must be finite")]
    InvalidAngularVelocity { angular_velocity_rad_per_s: f64 },
    #[error("motor angle {angle_rad} rad must be finite")]
    InvalidAngle { angle_rad: f64 },
    #[error("motor torque {torque_n_m} N*m must be finite")]
    InvalidTorque { torque_n_m: f64 },
    #[error("motor pivot {pivot_m:?} m must be finite")]
    InvalidPivot { pivot_m: [f64; 3] },
    #[error("motor drives no atoms")]
    EmptyMotor,
    #[error("motor geometry is degenerate: the driven atoms have no spread about the axis")]
    DegenerateMotorGeometry,
    #[error("unit mismatch: {from} cannot be read as {to}")]
    UnitMismatch { from: Unit, to: Unit },
}
