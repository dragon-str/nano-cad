use nanocad_model::Topology;

use crate::error::JigError;
use crate::jig::{atom_position_m, validate_positions, Jig, JigKind};

/// A motor jig: it drives selected atoms about an axis.
///
/// The motor is a torque source, not a potential. Its state is the rotation
/// angle `angle_rad`, the angular velocity `angular_velocity_rad_per_s`, the
/// unit `axis`, the `pivot_m` point on the axis, and the constant source
/// `torque_n_m`. [`Jig::energy_j`] and [`Jig::gradient_j_per_m`] are zero,
/// because a constant torque is path dependent. [`Jig::forces_n`] returns the
/// torque-distributed force instead. The report must distinguish the two
/// cases; [`Jig::kind`] names the jig as [`JigKind::Motor`].
///
/// The motor spreads its torque over the driven atoms. It removes the axial
/// part of each atom offset first, then applies the minimum-norm force set
/// that gives zero net force and the exact torque about the axis. This needs
/// at least two driven atoms with different in-plane offsets, or the call
/// returns [`JigError::DegenerateMotorGeometry`].
#[derive(Clone, Debug, PartialEq)]
pub struct MotorJig {
    atom_count: usize,
    axis: [f64; 3],
    pivot_m: [f64; 3],
    angular_velocity_rad_per_s: f64,
    torque_n_m: f64,
    angle_rad: f64,
    atom_indices: Vec<u32>,
}

impl MotorJig {
    /// Creates a motor over a topology.
    ///
    /// `atom_indices` selects the driven atoms. The axis must be a finite,
    /// non-zero vector. It is normalized to unit length. The pivot must be
    /// finite. The selection must not be empty.
    pub fn new(
        topology: &Topology,
        atom_indices: &[u32],
        axis: [f64; 3],
        pivot_m: [f64; 3],
    ) -> Result<Self, JigError> {
        if atom_indices.is_empty() {
            return Err(JigError::EmptyMotor);
        }
        let atom_count = topology.atom_count();
        for atom in atom_indices {
            if *atom as usize >= atom_count {
                return Err(JigError::AtomIndexOutOfBounds {
                    atom: *atom as usize,
                    atom_count,
                });
            }
        }
        let axis = normalize_axis(axis)?;
        if pivot_m.iter().any(|value| !value.is_finite()) {
            return Err(JigError::InvalidPivot { pivot_m });
        }
        Ok(Self {
            atom_count,
            axis,
            pivot_m,
            angular_velocity_rad_per_s: 0.0,
            torque_n_m: 0.0,
            angle_rad: 0.0,
            atom_indices: atom_indices.to_vec(),
        })
    }

    /// Returns the unit rotation axis.
    pub fn axis(&self) -> [f64; 3] {
        self.axis
    }

    /// Returns the pivot point in SI metres.
    pub fn pivot_m(&self) -> [f64; 3] {
        self.pivot_m
    }

    /// Returns the rotation angle in radians.
    pub fn angle_rad(&self) -> f64 {
        self.angle_rad
    }

    /// Returns the angular velocity in radians per second.
    pub fn angular_velocity_rad_per_s(&self) -> f64 {
        self.angular_velocity_rad_per_s
    }

    /// Returns the constant source torque about the axis in newton metres.
    pub fn torque_n_m(&self) -> f64 {
        self.torque_n_m
    }

    /// Returns the torque as a vector in newton metres.
    pub fn torque_vector_n_m(&self) -> [f64; 3] {
        [
            self.torque_n_m * self.axis[0],
            self.torque_n_m * self.axis[1],
            self.torque_n_m * self.axis[2],
        ]
    }

    /// Returns the driven atom indices.
    pub fn atom_indices(&self) -> &[u32] {
        &self.atom_indices
    }

    /// Replaces the rotation angle in radians. The value must be finite.
    pub fn set_angle_rad(&mut self, angle_rad: f64) -> Result<(), JigError> {
        if !angle_rad.is_finite() {
            return Err(JigError::InvalidAngle { angle_rad });
        }
        self.angle_rad = angle_rad;
        Ok(())
    }

    /// Replaces the angular velocity in radians per second. The value must be
    /// finite.
    pub fn set_angular_velocity_rad_per_s(
        &mut self,
        angular_velocity_rad_per_s: f64,
    ) -> Result<(), JigError> {
        if !angular_velocity_rad_per_s.is_finite() {
            return Err(JigError::InvalidAngularVelocity {
                angular_velocity_rad_per_s,
            });
        }
        self.angular_velocity_rad_per_s = angular_velocity_rad_per_s;
        Ok(())
    }

    /// Replaces the constant source torque in newton metres. The value must be
    /// finite.
    pub fn set_torque_n_m(&mut self, torque_n_m: f64) -> Result<(), JigError> {
        if !torque_n_m.is_finite() {
            return Err(JigError::InvalidTorque { torque_n_m });
        }
        self.torque_n_m = torque_n_m;
        Ok(())
    }

    /// Advances the angle by `angular_velocity_rad_per_s * dt_s`.
    ///
    /// `dt_s` must be a finite, non-negative time. The method does not
    /// integrate the atoms. It only advances the motor state.
    pub fn advance_angle_rad(&mut self, dt_s: f64) -> Result<(), JigError> {
        if !dt_s.is_finite() || dt_s < 0.0 {
            return Err(JigError::InvalidAngle { angle_rad: dt_s });
        }
        self.angle_rad += self.angular_velocity_rad_per_s * dt_s;
        Ok(())
    }

    /// Returns the torque-distributed force in newtons for a flat position
    /// buffer.
    pub fn torque_forces_n(&self, positions_m: &[f64]) -> Result<Vec<f64>, JigError> {
        validate_positions(positions_m, self.atom_count)?;
        let mut offsets = Vec::with_capacity(self.atom_indices.len());
        for atom in &self.atom_indices {
            let position_m = atom_position_m(positions_m, *atom);
            let d = [
                position_m[0] - self.pivot_m[0],
                position_m[1] - self.pivot_m[1],
                position_m[2] - self.pivot_m[2],
            ];
            let axial = dot(d, self.axis);
            offsets.push([
                d[0] - axial * self.axis[0],
                d[1] - axial * self.axis[1],
                d[2] - axial * self.axis[2],
            ]);
        }
        let count = offsets.len();
        let mut mean = [0.0; 3];
        for offset in &offsets {
            mean[0] += offset[0];
            mean[1] += offset[1];
            mean[2] += offset[2];
        }
        let inv_count = 1.0 / count as f64;
        mean[0] *= inv_count;
        mean[1] *= inv_count;
        mean[2] *= inv_count;

        let mut response_sq = 0.0;
        for offset in &mut offsets {
            offset[0] -= mean[0];
            offset[1] -= mean[1];
            offset[2] -= mean[2];
            response_sq += dot(*offset, *offset);
        }
        if response_sq <= 0.0 {
            return Err(JigError::DegenerateMotorGeometry);
        }

        let mut forces = vec![0.0; 3 * self.atom_count];
        let scale = self.torque_n_m / response_sq;
        for (slot, atom) in self.atom_indices.iter().enumerate() {
            let force = cross(self.axis, offsets[slot]);
            let base = *atom as usize * 3;
            forces[base] = scale * force[0];
            forces[base + 1] = scale * force[1];
            forces[base + 2] = scale * force[2];
        }
        Ok(forces)
    }
}

impl Jig for MotorJig {
    fn kind(&self) -> JigKind {
        JigKind::Motor
    }

    fn atom_count(&self) -> usize {
        self.atom_count
    }

    fn energy_j(&self, positions_m: &[f64]) -> Result<f64, JigError> {
        validate_positions(positions_m, self.atom_count)?;
        Ok(0.0)
    }

    fn gradient_j_per_m(&self, positions_m: &[f64]) -> Result<Vec<f64>, JigError> {
        validate_positions(positions_m, self.atom_count)?;
        Ok(vec![0.0; 3 * self.atom_count])
    }

    fn forces_n(&self, positions_m: &[f64]) -> Result<Vec<f64>, JigError> {
        self.torque_forces_n(positions_m)
    }
}

fn normalize_axis(axis: [f64; 3]) -> Result<[f64; 3], JigError> {
    if axis.iter().any(|value| !value.is_finite()) {
        return Err(JigError::InvalidAxis { axis });
    }
    let norm = dot(axis, axis).sqrt();
    if norm <= 0.0 {
        return Err(JigError::InvalidAxis { axis });
    }
    Ok([axis[0] / norm, axis[1] / norm, axis[2] / norm])
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use nanocad_model::{Atom, Element};

    fn rotor_topology(atom_count: usize, radius_m: f64) -> Topology {
        let mut topology = Topology::new();
        for atom in 0..atom_count {
            let angle = 2.0 * std::f64::consts::PI * atom as f64 / atom_count as f64;
            let position_m = [radius_m * angle.cos(), radius_m * angle.sin(), 0.0];
            topology.add_atom(Atom::new(Element::CARBON, position_m, 0.0, "C3"));
        }
        topology
    }

    fn rotor_motor(atom_count: usize, radius_m: f64, torque_n_m: f64) -> (Topology, MotorJig) {
        let topology = rotor_topology(atom_count, radius_m);
        let indices: Vec<u32> = (0..atom_count as u32).collect();
        let mut motor = MotorJig::new(&topology, &indices, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0])
            .expect("valid motor");
        motor.set_torque_n_m(torque_n_m).expect("valid torque");
        (topology, motor)
    }

    #[test]
    fn a_motor_is_a_torque_source_with_zero_potential() {
        let (topology, motor) = rotor_motor(3, 1.0e-9, 1.0e-21);
        assert_eq!(motor.kind(), JigKind::Motor);
        let positions_m = topology.positions_m();
        assert_eq!(motor.energy_j(positions_m).expect("valid"), 0.0);
        assert!(motor
            .gradient_j_per_m(positions_m)
            .expect("valid")
            .iter()
            .all(|value| *value == 0.0));
    }

    #[test]
    fn the_axis_is_normalized() {
        let topology = rotor_topology(3, 1.0e-9);
        let indices: Vec<u32> = (0..3).collect();
        let motor = MotorJig::new(&topology, &indices, [0.0, 0.0, 7.0], [0.0, 0.0, 0.0])
            .expect("valid motor");
        assert_eq!(motor.axis(), [0.0, 0.0, 1.0]);
    }

    #[test]
    fn an_invalid_axis_is_rejected() {
        let topology = rotor_topology(3, 1.0e-9);
        let indices: Vec<u32> = (0..3).collect();
        assert!(matches!(
            MotorJig::new(&topology, &indices, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
            Err(JigError::InvalidAxis { .. })
        ));
        assert!(matches!(
            MotorJig::new(&topology, &indices, [f64::NAN, 0.0, 1.0], [0.0, 0.0, 0.0]),
            Err(JigError::InvalidAxis { .. })
        ));
        assert!(matches!(
            MotorJig::new(
                &topology,
                &indices,
                [0.0, 0.0, 1.0],
                [f64::INFINITY, 0.0, 0.0]
            ),
            Err(JigError::InvalidPivot { .. })
        ));
        assert!(matches!(
            MotorJig::new(&topology, &[], [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]),
            Err(JigError::EmptyMotor)
        ));
        assert!(matches!(
            MotorJig::new(&topology, &[9], [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]),
            Err(JigError::AtomIndexOutOfBounds { .. })
        ));
    }

    #[test]
    fn the_torque_forces_have_zero_net_force_and_the_exact_axial_torque() {
        let radius_m = 1.0e-9;
        let torque_n_m = 4.0e-21;
        let (topology, motor) = rotor_motor(3, radius_m, torque_n_m);
        let positions_m = topology.positions_m();
        let forces = motor.forces_n(positions_m).expect("valid");

        let mut net_force = [0.0; 3];
        let mut net_torque = [0.0; 3];
        for (atom, indices) in (0..3).zip(motor.atom_indices().iter()) {
            let base = *indices as usize * 3;
            let force = [forces[base], forces[base + 1], forces[base + 2]];
            net_force[0] += force[0];
            net_force[1] += force[1];
            net_force[2] += force[2];
            let position_m = topology.position_m(atom).expect("atom exists");
            let moment = cross(position_m, force);
            net_torque[0] += moment[0];
            net_torque[1] += moment[1];
            net_torque[2] += moment[2];
        }
        for component in net_force {
            assert!(component.abs() < 1.0e-30, "net force {net_force:?}");
        }
        assert!(
            (net_torque[2] - torque_n_m).abs() < 1.0e-36,
            "{net_torque:?}"
        );
        assert!(net_torque[0].abs() < 1.0e-36);
        assert!(net_torque[1].abs() < 1.0e-36);
    }

    #[test]
    fn the_torque_forces_are_tangential_to_the_axis() {
        let radius_m = 1.0e-9;
        let (topology, motor) = rotor_motor(4, radius_m, 2.0e-21);
        let positions_m = topology.positions_m();
        let forces = motor.forces_n(positions_m).expect("valid");
        for atom in 0..4 {
            let base = atom * 3;
            let force = [forces[base], forces[base + 1], forces[base + 2]];
            let position_m = topology.position_m(atom).expect("atom exists");
            let radial = dot(position_m, force);
            assert!(radial.abs() < 1.0e-30);
            assert!(force[2] == 0.0);
        }
    }

    #[test]
    fn a_degenerate_motor_geometry_is_rejected() {
        let mut topology = Topology::new();
        topology.add_atom(Atom::new(Element::CARBON, [0.0, 0.0, 0.0], 0.0, "C3"));
        topology.add_atom(Atom::new(Element::CARBON, [0.0, 0.0, 1.0e-10], 0.0, "C3"));
        let motor = MotorJig::new(&topology, &[0, 1], [0.0, 0.0, 1.0], [0.0, 0.0, 0.0])
            .expect("valid motor setup");
        let positions_m = topology.positions_m();
        assert!(matches!(
            motor.forces_n(positions_m),
            Err(JigError::DegenerateMotorGeometry)
        ));
    }

    #[test]
    fn invalid_state_values_are_rejected() {
        let (_, mut motor) = rotor_motor(3, 1.0e-9, 1.0e-21);
        assert!(matches!(
            motor.set_angle_rad(f64::NAN),
            Err(JigError::InvalidAngle { .. })
        ));
        assert!(matches!(
            motor.set_angular_velocity_rad_per_s(f64::INFINITY),
            Err(JigError::InvalidAngularVelocity { .. })
        ));
        assert!(matches!(
            motor.set_torque_n_m(f64::NAN),
            Err(JigError::InvalidTorque { .. })
        ));
        assert!(matches!(
            motor.advance_angle_rad(-1.0),
            Err(JigError::InvalidAngle { .. })
        ));
    }

    #[test]
    fn advancing_the_angle_uses_the_angular_velocity() {
        let (_, mut motor) = rotor_motor(3, 1.0e-9, 1.0e-21);
        motor
            .set_angular_velocity_rad_per_s(2.0)
            .expect("valid velocity");
        motor.advance_angle_rad(0.5).expect("valid time step");
        assert_eq!(motor.angle_rad(), 1.0);
        motor.set_angle_rad(0.0).expect("valid angle");
        assert_eq!(motor.angle_rad(), 0.0);
    }

    #[test]
    fn the_torque_vector_points_along_the_axis() {
        let (_, motor) = rotor_motor(3, 1.0e-9, 3.0e-21);
        assert_eq!(motor.torque_vector_n_m(), [0.0, 0.0, 3.0e-21]);
    }

    #[test]
    fn a_position_buffer_of_the_wrong_size_is_rejected() {
        let (_, motor) = rotor_motor(3, 1.0e-9, 1.0e-21);
        assert!(matches!(
            motor.energy_j(&[0.0, 0.0, 0.0]),
            Err(JigError::PositionBufferSizeMismatch { .. })
        ));
    }
}
