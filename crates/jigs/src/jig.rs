use crate::error::JigError;

/// The kind of a jig.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum JigKind {
    Anchor,
    Motor,
    Spring,
    Thermostat,
    RigidBody,
    Joint,
}

/// A jig over a [`nanocad_model::Topology`].
///
/// A jig contributes an energy and a force to a simulation. Positions are in
/// SI metres, energies in joules, and forces in newtons. The force buffer uses
/// the same flat layout as the position buffer: entry `3 * atom + axis`.
///
/// # Conservative and non-conservative jigs
///
/// [`Jig::energy_j`] and [`Jig::gradient_j_per_m`] describe the conservative
/// part of the jig. For a potential jig, [`Jig::forces_n`] is the negative of
/// the gradient, and the finite-difference gradient check applies.
///
/// A motor is a torque source, not a potential. Its `energy_j` and
/// `gradient_j_per_m` are zero. It overrides [`Jig::forces_n`] to return the
/// torque-distributed force. The report must distinguish the two cases, so
/// [`Jig::kind`] names the jig and a caller can test it.
///
/// # Buffer sizes
///
/// Every jig stores the atom count of the topology it was built from. A
/// position or force buffer must hold exactly `3 * atom_count` values, or the
/// call returns [`JigError::PositionBufferSizeMismatch`].
pub trait Jig {
    /// Returns the kind of this jig.
    fn kind(&self) -> JigKind;

    /// Returns the number of atoms the jig was built from.
    fn atom_count(&self) -> usize;

    /// Returns the potential energy in joules for a flat position buffer.
    ///
    /// A non-conservative jig returns zero.
    fn energy_j(&self, positions_m: &[f64]) -> Result<f64, JigError>;

    /// Returns the potential gradient in newtons for a flat position buffer.
    ///
    /// A non-conservative jig returns a zero buffer.
    fn gradient_j_per_m(&self, positions_m: &[f64]) -> Result<Vec<f64>, JigError>;

    /// Returns the applied force in newtons for a flat position buffer.
    ///
    /// The default is the negative of [`Jig::gradient_j_per_m`]. A
    /// non-conservative jig overrides this method.
    fn forces_n(&self, positions_m: &[f64]) -> Result<Vec<f64>, JigError> {
        let mut forces = self.gradient_j_per_m(positions_m)?;
        for force in &mut forces {
            *force = -*force;
        }
        Ok(forces)
    }
}

/// Checks that a flat position buffer matches an atom count and is finite.
pub(crate) fn validate_positions(positions_m: &[f64], atom_count: usize) -> Result<(), JigError> {
    if positions_m.len() != 3 * atom_count {
        return Err(JigError::PositionBufferSizeMismatch {
            len: positions_m.len(),
            atom_count,
        });
    }
    for (index, value) in positions_m.iter().enumerate() {
        if !value.is_finite() {
            return Err(JigError::NonFinitePosition { index });
        }
    }
    Ok(())
}

/// Reads the position of an atom from a validated flat buffer.
pub(crate) fn atom_position_m(positions_m: &[f64], atom: u32) -> [f64; 3] {
    let base = atom as usize * 3;
    [
        positions_m[base],
        positions_m[base + 1],
        positions_m[base + 2],
    ]
}

/// Returns the Euclidean distance between two points in metres.
pub(crate) fn distance_m(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}
