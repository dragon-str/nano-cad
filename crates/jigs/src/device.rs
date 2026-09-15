//! The L2 device layer: rigid bodies and joints.
//!
//! A [`RigidBody`] has a mass, an inertia tensor, a center of mass, a pose
//! (position and quaternion orientation), a linear velocity, and an angular
//! velocity. All values are SI: kilograms, metres, seconds, radians.
//!
//! # Integrator
//!
//! The free step is semi-implicit (symplectic) Euler. The world angular
//! momentum is the primary angular state, and it obeys `dL/dt = tau` directly.
//! A body under zero torque therefore keeps its world angular momentum to
//! floating-point precision. The angular velocity is derived from
//! `omega = I_world^-1 L` after each orientation update.
//!
//! # Constraint method
//!
//! A joint constrains the relative motion of two bodies. The joints here use a
//! **position-level projection**. After integration, the solver runs a number
//! of Gauss-Seidel passes. Each pass removes the constraint error in the
//! mass-weighted least-squares sense. This method is simple and directly
//! drives the position error to zero. It does not remove the constraint error
//! velocity, so a constrained free body can carry a small residual velocity
//! between steps. Callers that need this must add a velocity-level pass.
//!
//! # Joints
//!
//! - [`RevoluteJoint`]: one rotation about a common axis, anchor points
//!   coincident. One degree of freedom.
//! - [`PrismaticJoint`]: translation along one axis only, relative orientation
//!   fixed. One degree of freedom.
//! - [`GearCoupling`]: two revolute joints obey
//!   `r_a * theta_a + r_b * theta_b = 0`.
//!
//! [`RigidBodySystem`] owns the bodies and the joints and runs the step loop.

use nanocad_units::{Quantity, Unit};
use thiserror::Error;

/// Errors from device setup and from the device step.
#[derive(Debug, Error, PartialEq)]
pub enum DeviceError {
    #[error("mass {mass_kg} kg must be finite and positive")]
    InvalidMass { mass_kg: f64 },
    #[error("inertia tensor {inertia_kg_m2:?} kg*m^2 must be finite")]
    NonFiniteInertia { inertia_kg_m2: [[f64; 3]; 3] },
    #[error("inertia tensor {inertia_kg_m2:?} kg*m^2 must be symmetric")]
    NonSymmetricInertia { inertia_kg_m2: [[f64; 3]; 3] },
    #[error("inertia tensor {inertia_kg_m2:?} kg*m^2 must be positive definite")]
    SingularInertia { inertia_kg_m2: [[f64; 3]; 3] },
    #[error("center of mass {center_of_mass_m:?} m must be finite")]
    NonFiniteCenterOfMass { center_of_mass_m: [f64; 3] },
    #[error("position {position_m:?} m must be finite")]
    NonFinitePosition { position_m: [f64; 3] },
    #[error("linear velocity {linear_velocity_m_per_s:?} m/s must be finite")]
    NonFiniteLinearVelocity { linear_velocity_m_per_s: [f64; 3] },
    #[error("angular velocity {angular_velocity_rad_per_s:?} rad/s must be finite")]
    NonFiniteAngularVelocity {
        angular_velocity_rad_per_s: [f64; 3],
    },
    #[error("orientation norm {norm} must be finite and positive")]
    InvalidOrientation { norm: f64 },
    #[error("time step {dt_s} s must be finite and non-negative")]
    InvalidTimeStep { dt_s: f64 },
    #[error("body index {body} is outside the body count {body_count}")]
    BodyIndexOutOfBounds { body: usize, body_count: usize },
    #[error("joint index {joint} is outside the joint count {joint_count}")]
    JointIndexOutOfBounds { joint: usize, joint_count: usize },
    #[error("a joint cannot connect body {body} to itself")]
    SelfJoint { body: usize },
    #[error("axis {axis:?} must be a finite, non-zero vector")]
    InvalidAxis { axis: [f64; 3] },
    #[error("anchor {anchor_m:?} m must be finite")]
    NonFiniteAnchor { anchor_m: [f64; 3] },
    #[error("gear radius {radius_m} m must be finite and positive")]
    InvalidRadius { radius_m: f64 },
    #[error("a gear coupling cannot reference joint {joint} twice")]
    DuplicateGearJoint { joint: usize },
    #[error("a gear constraint needs at least two joint terms, got {terms}")]
    TooFewGearTerms { terms: usize },
    #[error("a gear constraint cannot reference joint {joint} twice")]
    DuplicateGearConstraintJoint { joint: usize },
    #[error("gear constraint coefficient {coefficient_m} m for joint {joint} must be finite and non-zero")]
    InvalidGearCoefficient { joint: usize, coefficient_m: f64 },
    #[error("unit mismatch: {from} cannot be read as {to}")]
    UnitMismatch { from: Unit, to: Unit },
}

/// A rotation represented as a unit quaternion `w + x i + y j + z k`.
///
/// The quaternion maps a vector from the body frame to the world frame. Every
/// method treats the stored value as a rotation and normalizes it first, so a
/// caller may hold an unnormalized value between operations.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Quat {
    pub w: f64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Quat {
    /// The identity rotation.
    pub const IDENTITY: Quat = Quat {
        w: 1.0,
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };

    /// Builds a quaternion from its four components.
    pub const fn new(w: f64, x: f64, y: f64, z: f64) -> Self {
        Self { w, x, y, z }
    }

    /// Builds the identity rotation.
    pub const fn identity() -> Self {
        Self::IDENTITY
    }

    /// Builds a rotation of `angle_rad` about `axis`.
    ///
    /// The axis is normalized. A zero axis gives the identity rotation.
    pub fn from_axis_angle(axis: [f64; 3], angle_rad: f64) -> Self {
        let norm = norm3(axis);
        if !norm.is_finite() || norm <= 1.0e-300 {
            return Self::IDENTITY;
        }
        let (sin_half, cos_half) = (0.5 * angle_rad).sin_cos();
        let scale = sin_half / norm;
        Self::new(cos_half, scale * axis[0], scale * axis[1], scale * axis[2])
    }

    /// Returns the Euclidean norm of the quaternion.
    pub fn norm(&self) -> f64 {
        (self.w * self.w + self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// Returns true when every component is finite.
    pub fn is_finite(&self) -> bool {
        self.w.is_finite() && self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }

    /// Returns the rotation with unit norm.
    ///
    /// A zero or non-finite quaternion maps to the identity.
    pub fn normalized(&self) -> Self {
        let norm = self.norm();
        if !norm.is_finite() || norm <= 1.0e-300 {
            return Self::IDENTITY;
        }
        let inv = 1.0 / norm;
        Self::new(self.w * inv, self.x * inv, self.y * inv, self.z * inv)
    }

    /// Returns the conjugate, which is the inverse for a unit quaternion.
    pub fn conjugate(&self) -> Self {
        Self::new(self.w, -self.x, -self.y, -self.z)
    }

    /// Returns the quaternion scaled by `factor`.
    pub fn scale(&self, factor: f64) -> Self {
        Self::new(
            self.w * factor,
            self.x * factor,
            self.y * factor,
            self.z * factor,
        )
    }

    /// Rotates a world-frame vector into the body frame.
    pub fn inverse_rotate(&self, v: [f64; 3]) -> [f64; 3] {
        self.conjugate().rotate(v)
    }

    /// Rotates a body-frame vector into the world frame.
    pub fn rotate(&self, v: [f64; 3]) -> [f64; 3] {
        mat3_mul_vec(self.to_rotation_matrix(), v)
    }

    /// Returns the equivalent rotation matrix.
    pub fn to_rotation_matrix(&self) -> [[f64; 3]; 3] {
        let q = self.normalized();
        let (w, x, y, z) = (q.w, q.x, q.y, q.z);
        [
            [
                1.0 - 2.0 * (y * y + z * z),
                2.0 * (x * y - w * z),
                2.0 * (x * z + w * y),
            ],
            [
                2.0 * (x * y + w * z),
                1.0 - 2.0 * (x * x + z * z),
                2.0 * (y * z - w * x),
            ],
            [
                2.0 * (x * z - w * y),
                2.0 * (y * z + w * x),
                1.0 - 2.0 * (x * x + y * y),
            ],
        ]
    }

    /// Returns the shortest axis-angle form `(unit_axis, angle_rad)`.
    ///
    /// The angle lies in `(-pi, pi]`. A zero rotation returns the z axis.
    pub fn to_axis_angle(&self) -> ([f64; 3], f64) {
        let q = self.normalized();
        let vector_norm = (q.x * q.x + q.y * q.y + q.z * q.z).sqrt();
        if vector_norm <= 1.0e-300 {
            return ([0.0, 0.0, 1.0], 0.0);
        }
        let mut axis = [q.x / vector_norm, q.y / vector_norm, q.z / vector_norm];
        let mut angle = 2.0 * vector_norm.atan2(q.w);
        if angle > std::f64::consts::PI {
            angle -= 2.0 * std::f64::consts::PI;
            axis = scale3(axis, -1.0);
        }
        (axis, angle)
    }
}

impl std::ops::Mul for Quat {
    type Output = Quat;

    fn mul(self, rhs: Quat) -> Quat {
        Quat::new(
            self.w * rhs.w - self.x * rhs.x - self.y * rhs.y - self.z * rhs.z,
            self.w * rhs.x + self.x * rhs.w + self.y * rhs.z - self.z * rhs.y,
            self.w * rhs.y - self.x * rhs.z + self.y * rhs.w + self.z * rhs.x,
            self.w * rhs.z + self.x * rhs.y - self.y * rhs.x + self.z * rhs.w,
        )
    }
}

impl std::ops::Add for Quat {
    type Output = Quat;

    fn add(self, rhs: Quat) -> Quat {
        Quat::new(
            self.w + rhs.w,
            self.x + rhs.x,
            self.y + rhs.y,
            self.z + rhs.z,
        )
    }
}

/// A rigid body in the L2 device layer.
///
/// `position_m` is the world position of the center of mass.
/// `center_of_mass_m` is the center of mass in the body frame, relative to the
/// body's geometric origin. A local point maps to the world by
/// `position_m + R * (local - center_of_mass_m)`.
///
/// The inertia tensor is about the center of mass, in the body frame. It must
/// be symmetric and positive definite.
///
/// A fixed body does not move. Its inverse mass and inverse inertia are zero,
/// and [`RigidBody::integrate`] returns without a change.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidBody {
    mass_kg: f64,
    inertia_kg_m2: [[f64; 3]; 3],
    inverse_inertia_kg_m2: [[f64; 3]; 3],
    center_of_mass_m: [f64; 3],
    position_m: [f64; 3],
    orientation: Quat,
    linear_velocity_m_per_s: [f64; 3],
    angular_momentum_kg_m2_per_s: [f64; 3],
    force_n: [f64; 3],
    torque_n_m: [f64; 3],
    fixed: bool,
}

impl RigidBody {
    /// Creates a body.
    ///
    /// The mass must be finite and positive. The inertia tensor must be finite,
    /// symmetric, and positive definite. Every vector must be finite. The
    /// orientation is normalized.
    pub fn new(
        mass_kg: f64,
        inertia_kg_m2: [[f64; 3]; 3],
        center_of_mass_m: [f64; 3],
        position_m: [f64; 3],
        orientation: Quat,
    ) -> Result<Self, DeviceError> {
        if !mass_kg.is_finite() || mass_kg <= 0.0 {
            return Err(DeviceError::InvalidMass { mass_kg });
        }
        if inertia_kg_m2
            .iter()
            .flatten()
            .any(|value| !value.is_finite())
        {
            return Err(DeviceError::NonFiniteInertia { inertia_kg_m2 });
        }
        if !is_symmetric(inertia_kg_m2) {
            return Err(DeviceError::NonSymmetricInertia { inertia_kg_m2 });
        }
        let inverse_inertia_kg_m2 = mat3_inverse(inertia_kg_m2)
            .filter(|_| is_positive_definite(inertia_kg_m2))
            .ok_or(DeviceError::SingularInertia { inertia_kg_m2 })?;
        if center_of_mass_m.iter().any(|value| !value.is_finite()) {
            return Err(DeviceError::NonFiniteCenterOfMass { center_of_mass_m });
        }
        if position_m.iter().any(|value| !value.is_finite()) {
            return Err(DeviceError::NonFinitePosition { position_m });
        }
        if !orientation.is_finite() || orientation.norm() <= 1.0e-300 {
            return Err(DeviceError::InvalidOrientation {
                norm: orientation.norm(),
            });
        }
        Ok(Self {
            mass_kg,
            inertia_kg_m2,
            inverse_inertia_kg_m2,
            center_of_mass_m,
            position_m,
            orientation: orientation.normalized(),
            linear_velocity_m_per_s: [0.0; 3],
            angular_momentum_kg_m2_per_s: [0.0; 3],
            force_n: [0.0; 3],
            torque_n_m: [0.0; 3],
            fixed: false,
        })
    }

    /// Builds a body with a diagonal inertia tensor `(ixx, iyy, izz)`.
    pub fn with_diagonal_inertia(
        mass_kg: f64,
        inertia_diagonal_kg_m2: [f64; 3],
        position_m: [f64; 3],
        orientation: Quat,
    ) -> Result<Self, DeviceError> {
        let [ixx, iyy, izz] = inertia_diagonal_kg_m2;
        let inertia_kg_m2 = [[ixx, 0.0, 0.0], [0.0, iyy, 0.0], [0.0, 0.0, izz]];
        Self::new(mass_kg, inertia_kg_m2, [0.0; 3], position_m, orientation)
    }

    /// Returns the mass in kilograms.
    pub fn mass_kg(&self) -> f64 {
        self.mass_kg
    }

    /// Returns the mass as a quantity.
    pub fn mass_quantity(&self) -> Quantity {
        Quantity::from_si(self.mass_kg, Unit::Kilogram)
    }

    /// Returns the inverse mass in inverse kilograms. A fixed body returns 0.
    pub fn inverse_mass_kg(&self) -> f64 {
        if self.fixed {
            0.0
        } else {
            1.0 / self.mass_kg
        }
    }

    /// Returns true when the body is fixed.
    pub fn is_fixed(&self) -> bool {
        self.fixed
    }

    /// Fixes or releases the body.
    ///
    /// A fixed body ignores its accumulated force and torque.
    pub fn set_fixed(&mut self, fixed: bool) {
        self.fixed = fixed;
        if fixed {
            self.linear_velocity_m_per_s = [0.0; 3];
            self.angular_momentum_kg_m2_per_s = [0.0; 3];
        }
    }

    /// Returns the body-frame inertia tensor about the center of mass.
    pub fn inertia_kg_m2(&self) -> [[f64; 3]; 3] {
        self.inertia_kg_m2
    }

    /// Returns the inverse body-frame inertia tensor.
    pub fn inverse_inertia_kg_m2(&self) -> [[f64; 3]; 3] {
        self.inverse_inertia_kg_m2
    }

    /// Returns the world-frame inertia tensor `R I R^T`.
    pub fn inertia_world_kg_m2(&self) -> [[f64; 3]; 3] {
        let r = self.orientation.to_rotation_matrix();
        mat3_mul(mat3_mul(r, self.inertia_kg_m2), mat3_transpose(r))
    }

    /// Returns the inverse world-frame inertia tensor. A fixed body returns 0.
    pub fn inverse_inertia_world(&self) -> [[f64; 3]; 3] {
        if self.fixed {
            return [[0.0; 3]; 3];
        }
        let r = self.orientation.to_rotation_matrix();
        mat3_mul(mat3_mul(r, self.inverse_inertia_kg_m2), mat3_transpose(r))
    }

    /// Returns the center of mass in the body frame.
    pub fn center_of_mass_m(&self) -> [f64; 3] {
        self.center_of_mass_m
    }

    /// Returns the world position of the center of mass.
    pub fn position_m(&self) -> [f64; 3] {
        self.position_m
    }

    /// Returns the world position of the center of mass as a quantity per axis.
    pub fn position_quantity_m(&self) -> [Quantity; 3] {
        [
            Quantity::from_si(self.position_m[0], Unit::Metre),
            Quantity::from_si(self.position_m[1], Unit::Metre),
            Quantity::from_si(self.position_m[2], Unit::Metre),
        ]
    }

    /// Returns the orientation quaternion.
    pub fn orientation(&self) -> Quat {
        self.orientation
    }

    /// Returns the world linear velocity in metres per second.
    pub fn linear_velocity_m_per_s(&self) -> [f64; 3] {
        self.linear_velocity_m_per_s
    }

    /// Returns the world angular velocity `I_world^-1 L` in radians per second.
    pub fn angular_velocity_rad_per_s(&self) -> [f64; 3] {
        mat3_mul_vec(
            self.inverse_inertia_world(),
            self.angular_momentum_kg_m2_per_s,
        )
    }

    /// Returns the world angular momentum in kilogram metre squared per second.
    pub fn angular_momentum_kg_m2_per_s(&self) -> [f64; 3] {
        self.angular_momentum_kg_m2_per_s
    }

    /// Returns the accumulated force in newtons.
    pub fn force_n(&self) -> [f64; 3] {
        self.force_n
    }

    /// Returns the accumulated torque about the center of mass in newton metres.
    pub fn torque_n_m(&self) -> [f64; 3] {
        self.torque_n_m
    }

    /// Replaces the world position of the center of mass.
    pub fn set_position_m(&mut self, position_m: [f64; 3]) -> Result<(), DeviceError> {
        if position_m.iter().any(|value| !value.is_finite()) {
            return Err(DeviceError::NonFinitePosition { position_m });
        }
        self.position_m = position_m;
        Ok(())
    }

    /// Replaces the orientation. The value is normalized.
    pub fn set_orientation(&mut self, orientation: Quat) -> Result<(), DeviceError> {
        if !orientation.is_finite() || orientation.norm() <= 1.0e-300 {
            return Err(DeviceError::InvalidOrientation {
                norm: orientation.norm(),
            });
        }
        self.orientation = orientation.normalized();
        Ok(())
    }

    /// Replaces the world linear velocity in metres per second.
    pub fn set_linear_velocity_m_per_s(
        &mut self,
        linear_velocity_m_per_s: [f64; 3],
    ) -> Result<(), DeviceError> {
        if linear_velocity_m_per_s
            .iter()
            .any(|value| !value.is_finite())
        {
            return Err(DeviceError::NonFiniteLinearVelocity {
                linear_velocity_m_per_s,
            });
        }
        self.linear_velocity_m_per_s = linear_velocity_m_per_s;
        Ok(())
    }

    /// Replaces the world angular velocity in radians per second.
    ///
    /// The method stores the matching angular momentum `I_world * omega`.
    pub fn set_angular_velocity_rad_per_s(
        &mut self,
        angular_velocity_rad_per_s: [f64; 3],
    ) -> Result<(), DeviceError> {
        if angular_velocity_rad_per_s
            .iter()
            .any(|value| !value.is_finite())
        {
            return Err(DeviceError::NonFiniteAngularVelocity {
                angular_velocity_rad_per_s,
            });
        }
        self.angular_momentum_kg_m2_per_s =
            mat3_mul_vec(self.inertia_world_kg_m2(), angular_velocity_rad_per_s);
        Ok(())
    }

    /// Maps a local body-frame point to the world frame.
    pub fn point_world_m(&self, local_body_m: [f64; 3]) -> [f64; 3] {
        let offset = sub3(local_body_m, self.center_of_mass_m);
        add3(self.position_m, self.orientation.rotate(offset))
    }

    /// Adds a force through the center of mass, in newtons.
    pub fn add_force_n(&mut self, force_n: [f64; 3]) {
        self.force_n = add3(self.force_n, force_n);
    }

    /// Adds a pure torque about the center of mass, in newton metres.
    pub fn add_torque_n_m(&mut self, torque_n_m: [f64; 3]) {
        self.torque_n_m = add3(self.torque_n_m, torque_n_m);
    }

    /// Adds a force at a world point. The offset produces a torque.
    pub fn add_force_at_point_n(&mut self, force_n: [f64; 3], point_world_m: [f64; 3]) {
        self.force_n = add3(self.force_n, force_n);
        let arm = sub3(point_world_m, self.position_m);
        self.torque_n_m = add3(self.torque_n_m, cross(arm, force_n));
    }

    /// Clears the accumulated force and torque.
    pub fn clear_forces(&mut self) {
        self.force_n = [0.0; 3];
        self.torque_n_m = [0.0; 3];
    }

    /// Returns the linear momentum `m * v` in kilogram metres per second.
    pub fn momentum_kg_m_per_s(&self) -> [f64; 3] {
        scale3(self.linear_velocity_m_per_s, self.mass_kg)
    }

    /// Returns the kinetic energy in joules.
    pub fn kinetic_energy_j(&self) -> f64 {
        let linear = self.mass_kg * dot(self.linear_velocity_m_per_s, self.linear_velocity_m_per_s);
        let angular = dot(
            self.angular_velocity_rad_per_s(),
            self.angular_momentum_kg_m2_per_s,
        );
        0.5 * linear + 0.5 * angular
    }

    /// Rotates the body by `angle_rad` about a world `point_m`.
    ///
    /// The orientation turns, and the center of mass orbits the point.
    pub fn rotate_about_point(&mut self, axis: [f64; 3], angle_rad: f64, point_m: [f64; 3]) {
        let delta = Quat::from_axis_angle(axis, angle_rad);
        let arm = sub3(self.position_m, point_m);
        self.position_m = add3(point_m, delta.rotate(arm));
        self.orientation = (delta * self.orientation).normalized();
    }

    /// Applies a world rotation given by an axis-angle vector.
    pub(crate) fn apply_rotation_delta(&mut self, rotation_delta_rad: [f64; 3]) {
        let angle = norm3(rotation_delta_rad);
        if angle > 1.0e-15 {
            let axis = scale3(rotation_delta_rad, 1.0 / angle);
            let delta = Quat::from_axis_angle(axis, angle);
            self.orientation = (delta * self.orientation).normalized();
        }
    }

    /// Advances the body by `dt_s` with the accumulated force and torque.
    ///
    /// The method uses semi-implicit Euler. The position level constraint
    /// projection runs later, in [`RigidBodySystem::step`].
    pub fn integrate(&mut self, dt_s: f64) -> Result<(), DeviceError> {
        if !dt_s.is_finite() || dt_s < 0.0 {
            return Err(DeviceError::InvalidTimeStep { dt_s });
        }
        if self.fixed {
            return Ok(());
        }
        self.linear_velocity_m_per_s = add3(
            self.linear_velocity_m_per_s,
            scale3(self.force_n, self.inverse_mass_kg() * dt_s),
        );
        self.position_m = add3(self.position_m, scale3(self.linear_velocity_m_per_s, dt_s));

        self.angular_momentum_kg_m2_per_s = add3(
            self.angular_momentum_kg_m2_per_s,
            scale3(self.torque_n_m, dt_s),
        );
        let omega = self.angular_velocity_rad_per_s();
        let omega_quat = Quat::new(0.0, omega[0], omega[1], omega[2]);
        let delta = (omega_quat * self.orientation).scale(0.5 * dt_s);
        self.orientation = (self.orientation + delta).normalized();
        Ok(())
    }
}

/// A revolute joint between two bodies.
///
/// The joint holds the two anchor points together and keeps the two axes
/// parallel. The relative rotation about the common axis is free, so the joint
/// has one degree of freedom.
///
/// Each body stores its anchor and axis in its own frame. Build the joint from
/// local frames, or use [`RevoluteJoint::from_world`] to project one world
/// anchor and axis onto both frames.
#[derive(Clone, Debug, PartialEq)]
pub struct RevoluteJoint {
    body_a: usize,
    body_b: usize,
    anchor_a_body_m: [f64; 3],
    anchor_b_body_m: [f64; 3],
    axis_a_body: [f64; 3],
    axis_b_body: [f64; 3],
}

impl RevoluteJoint {
    /// Creates a revolute joint from local-frame data.
    ///
    /// The two bodies must differ and must exist in `bodies`. The axes must be
    /// finite and non-zero. They are normalized. The anchors must be finite.
    pub fn new(
        bodies: &[RigidBody],
        body_a: usize,
        body_b: usize,
        anchor_a_body_m: [f64; 3],
        anchor_b_body_m: [f64; 3],
        axis_a_body: [f64; 3],
        axis_b_body: [f64; 3],
    ) -> Result<Self, DeviceError> {
        let joint = Self {
            body_a,
            body_b,
            anchor_a_body_m: check_anchor(anchor_a_body_m)?,
            anchor_b_body_m: check_anchor(anchor_b_body_m)?,
            axis_a_body: check_axis(axis_a_body)?,
            axis_b_body: check_axis(axis_b_body)?,
        };
        joint.check_bodies(bodies.len())?;
        Ok(joint)
    }

    /// Creates a revolute joint that starts satisfied.
    ///
    /// The world `anchor_m` and `axis` are projected onto each body frame. The
    /// joint error is zero at construction.
    pub fn from_world(
        bodies: &[RigidBody],
        body_a: usize,
        body_b: usize,
        anchor_m: [f64; 3],
        axis: [f64; 3],
    ) -> Result<Self, DeviceError> {
        let a = body_at(bodies, body_a)?;
        let b = body_at(bodies, body_b)?;
        let anchor_a_body_m = a.orientation.inverse_rotate(sub3(anchor_m, a.position_m));
        let anchor_b_body_m = b.orientation.inverse_rotate(sub3(anchor_m, b.position_m));
        let axis = check_axis(axis)?;
        let axis_a_body = a.orientation.inverse_rotate(axis);
        let axis_b_body = b.orientation.inverse_rotate(axis);
        Self::new(
            bodies,
            body_a,
            body_b,
            add3(anchor_a_body_m, a.center_of_mass_m),
            add3(anchor_b_body_m, b.center_of_mass_m),
            axis_a_body,
            axis_b_body,
        )
    }

    fn check_bodies(&self, body_count: usize) -> Result<(), DeviceError> {
        if self.body_a == self.body_b {
            return Err(DeviceError::SelfJoint { body: self.body_a });
        }
        if self.body_a >= body_count || self.body_b >= body_count {
            return Err(DeviceError::BodyIndexOutOfBounds {
                body: self.body_a.max(self.body_b),
                body_count,
            });
        }
        Ok(())
    }

    /// Returns the first body index.
    pub fn body_a(&self) -> usize {
        self.body_a
    }

    /// Returns the second body index.
    pub fn body_b(&self) -> usize {
        self.body_b
    }

    /// Returns the unit axis of the first body in its own frame.
    pub fn axis_a_body(&self) -> [f64; 3] {
        self.axis_a_body
    }

    /// Returns the unit axis of the second body in its own frame.
    pub fn axis_b_body(&self) -> [f64; 3] {
        self.axis_b_body
    }

    /// Returns the world anchor point of the first body.
    pub fn anchor_a_world_m(&self, bodies: &[RigidBody]) -> Result<[f64; 3], DeviceError> {
        let a = body_at(bodies, self.body_a)?;
        Ok(a.point_world_m(self.anchor_a_body_m))
    }

    /// Returns the world anchor point of the second body.
    pub fn anchor_b_world_m(&self, bodies: &[RigidBody]) -> Result<[f64; 3], DeviceError> {
        let b = body_at(bodies, self.body_b)?;
        Ok(b.point_world_m(self.anchor_b_body_m))
    }

    /// Returns the world axis of the first body.
    pub fn axis_a_world(&self, bodies: &[RigidBody]) -> Result<[f64; 3], DeviceError> {
        let a = body_at(bodies, self.body_a)?;
        Ok(a.orientation.rotate(self.axis_a_body))
    }

    /// Returns the world axis of the second body.
    pub fn axis_b_world(&self, bodies: &[RigidBody]) -> Result<[f64; 3], DeviceError> {
        let b = body_at(bodies, self.body_b)?;
        Ok(b.orientation.rotate(self.axis_b_body))
    }

    /// Returns the distance between the two anchor points in metres.
    pub fn anchor_error_m(&self, bodies: &[RigidBody]) -> Result<f64, DeviceError> {
        let a = self.anchor_a_world_m(bodies)?;
        let b = self.anchor_b_world_m(bodies)?;
        Ok(norm3(sub3(b, a)))
    }

    /// Returns the angle between the two axes in radians.
    pub fn axis_error_rad(&self, bodies: &[RigidBody]) -> Result<f64, DeviceError> {
        let a = self.axis_a_world(bodies)?;
        let b = self.axis_b_world(bodies)?;
        let sin_theta = norm3(cross(b, a));
        let cos_theta = dot(a, b);
        Ok(sin_theta.atan2(cos_theta))
    }

    /// Returns the relative rotation angle of body B about the common axis.
    ///
    /// The sign is positive when body B turns counterclockwise about the axis
    /// of body A. The value wraps at `+/- pi`.
    pub fn angle_rad(&self, bodies: &[RigidBody]) -> Result<f64, DeviceError> {
        let a = body_at(bodies, self.body_a)?;
        let b = body_at(bodies, self.body_b)?;
        let relative = a.orientation.conjugate() * b.orientation;
        let sin_half = relative.x * self.axis_a_body[0]
            + relative.y * self.axis_a_body[1]
            + relative.z * self.axis_a_body[2];
        Ok(2.0 * sin_half.atan2(relative.w))
    }

    /// Projects the joint error to zero in the mass-weighted least-squares sense.
    ///
    /// The anchor point uses three point rows. The axis uses the shortest
    /// rotation that aligns body B to body A. `relaxation` in `(0, 1]` scales
    /// the correction. One call is one Gauss-Seidel pass.
    pub fn project(&self, bodies: &mut [RigidBody], relaxation: f64) -> Result<(), DeviceError> {
        let local_a = sub3(
            self.anchor_a_body_m,
            body_at(bodies, self.body_a)?.center_of_mass_m,
        );
        let local_b = sub3(
            self.anchor_b_body_m,
            body_at(bodies, self.body_b)?.center_of_mass_m,
        );
        let axes = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        project_point_rows(
            bodies,
            self.body_a,
            local_a,
            self.body_b,
            local_b,
            &axes,
            relaxation,
        )?;
        self.align_axes(bodies, relaxation)
    }

    fn align_axes(&self, bodies: &mut [RigidBody], relaxation: f64) -> Result<(), DeviceError> {
        let (a, b) = two_bodies_mut(bodies, self.body_a, self.body_b)?;
        let axis_a = a.orientation.rotate(self.axis_a_body);
        let axis_b = b.orientation.rotate(self.axis_b_body);
        let cross_ab = cross(axis_b, axis_a);
        let sin_theta = norm3(cross_ab);
        let theta = sin_theta.atan2(dot(axis_a, axis_b));
        if theta <= 1.0e-12 || sin_theta <= 1.0e-15 {
            return Ok(());
        }
        let direction = scale3(cross_ab, 1.0 / sin_theta);
        let inverse_inertia_a = a.inverse_inertia_world();
        let inverse_inertia_b = b.inverse_inertia_world();
        let weight_a = dot(direction, mat3_mul_vec(inverse_inertia_a, direction));
        let weight_b = dot(direction, mat3_mul_vec(inverse_inertia_b, direction));
        let total = weight_a + weight_b;
        if total <= 1.0e-30 {
            return Ok(());
        }
        let angle_a = theta * weight_a / total * relaxation;
        let angle_b = theta * weight_b / total * relaxation;
        a.apply_rotation_delta(scale3(direction, angle_a));
        b.apply_rotation_delta(scale3(direction, angle_b));
        Ok(())
    }
}

/// A prismatic joint between two bodies.
///
/// The joint allows translation along one common axis. It holds the two anchor
/// points together in the plane perpendicular to the axis, and it fixes the
/// relative orientation. The joint has one degree of freedom.
#[derive(Clone, Debug, PartialEq)]
pub struct PrismaticJoint {
    body_a: usize,
    body_b: usize,
    anchor_a_body_m: [f64; 3],
    anchor_b_body_m: [f64; 3],
    axis_a_body: [f64; 3],
    axis_b_body: [f64; 3],
}

impl PrismaticJoint {
    /// Creates a prismatic joint from local-frame data.
    ///
    /// The two bodies must differ and must exist in `bodies`. The axes must be
    /// finite and non-zero. They are normalized. The anchors must be finite.
    pub fn new(
        bodies: &[RigidBody],
        body_a: usize,
        body_b: usize,
        anchor_a_body_m: [f64; 3],
        anchor_b_body_m: [f64; 3],
        axis_a_body: [f64; 3],
        axis_b_body: [f64; 3],
    ) -> Result<Self, DeviceError> {
        let joint = Self {
            body_a,
            body_b,
            anchor_a_body_m: check_anchor(anchor_a_body_m)?,
            anchor_b_body_m: check_anchor(anchor_b_body_m)?,
            axis_a_body: check_axis(axis_a_body)?,
            axis_b_body: check_axis(axis_b_body)?,
        };
        joint.check_bodies(bodies.len())?;
        Ok(joint)
    }

    /// Creates a prismatic joint that starts satisfied.
    pub fn from_world(
        bodies: &[RigidBody],
        body_a: usize,
        body_b: usize,
        anchor_m: [f64; 3],
        axis: [f64; 3],
    ) -> Result<Self, DeviceError> {
        let a = body_at(bodies, body_a)?;
        let b = body_at(bodies, body_b)?;
        let anchor_a_body_m = a.orientation.inverse_rotate(sub3(anchor_m, a.position_m));
        let anchor_b_body_m = b.orientation.inverse_rotate(sub3(anchor_m, b.position_m));
        let axis = check_axis(axis)?;
        let axis_a_body = a.orientation.inverse_rotate(axis);
        let axis_b_body = b.orientation.inverse_rotate(axis);
        Self::new(
            bodies,
            body_a,
            body_b,
            add3(anchor_a_body_m, a.center_of_mass_m),
            add3(anchor_b_body_m, b.center_of_mass_m),
            axis_a_body,
            axis_b_body,
        )
    }

    fn check_bodies(&self, body_count: usize) -> Result<(), DeviceError> {
        if self.body_a == self.body_b {
            return Err(DeviceError::SelfJoint { body: self.body_a });
        }
        if self.body_a >= body_count || self.body_b >= body_count {
            return Err(DeviceError::BodyIndexOutOfBounds {
                body: self.body_a.max(self.body_b),
                body_count,
            });
        }
        Ok(())
    }

    /// Returns the first body index.
    pub fn body_a(&self) -> usize {
        self.body_a
    }

    /// Returns the second body index.
    pub fn body_b(&self) -> usize {
        self.body_b
    }

    /// Returns the unit axis of the first body in its own frame.
    pub fn axis_a_body(&self) -> [f64; 3] {
        self.axis_a_body
    }

    /// Returns the world axis of the first body.
    pub fn axis_a_world(&self, bodies: &[RigidBody]) -> Result<[f64; 3], DeviceError> {
        let a = body_at(bodies, self.body_a)?;
        Ok(a.orientation.rotate(self.axis_a_body))
    }

    /// Returns the signed translation along the axis in metres.
    pub fn translation_m(&self, bodies: &[RigidBody]) -> Result<f64, DeviceError> {
        let a = body_at(bodies, self.body_a)?;
        let b = body_at(bodies, self.body_b)?;
        let anchor_a = a.point_world_m(self.anchor_a_body_m);
        let anchor_b = b.point_world_m(self.anchor_b_body_m);
        let axis_a_world = a.orientation.rotate(self.axis_a_body);
        Ok(dot(sub3(anchor_b, anchor_a), axis_a_world))
    }

    /// Returns the distance between the anchors perpendicular to the axis.
    pub fn position_error_m(&self, bodies: &[RigidBody]) -> Result<f64, DeviceError> {
        let a = body_at(bodies, self.body_a)?;
        let b = body_at(bodies, self.body_b)?;
        let anchor_a = a.point_world_m(self.anchor_a_body_m);
        let anchor_b = b.point_world_m(self.anchor_b_body_m);
        let axis_a_world = a.orientation.rotate(self.axis_a_body);
        let delta = sub3(anchor_b, anchor_a);
        Ok(norm3(sub3(
            delta,
            scale3(axis_a_world, dot(delta, axis_a_world)),
        )))
    }

    /// Returns the relative orientation error in radians.
    pub fn orientation_error_rad(&self, bodies: &[RigidBody]) -> Result<f64, DeviceError> {
        let a = body_at(bodies, self.body_a)?;
        let b = body_at(bodies, self.body_b)?;
        let relative = a.orientation.conjugate() * b.orientation;
        let (_, angle) = relative.to_axis_angle();
        Ok(angle)
    }

    /// Projects the joint error to zero in the mass-weighted least-squares sense.
    pub fn project(&self, bodies: &mut [RigidBody], relaxation: f64) -> Result<(), DeviceError> {
        let local_a = sub3(
            self.anchor_a_body_m,
            body_at(bodies, self.body_a)?.center_of_mass_m,
        );
        let local_b = sub3(
            self.anchor_b_body_m,
            body_at(bodies, self.body_b)?.center_of_mass_m,
        );
        let axis = self.axis_a_world(bodies)?;
        let (u, v) = perpendicular_basis(axis);
        project_point_rows(
            bodies,
            self.body_a,
            local_a,
            self.body_b,
            local_b,
            &[u, v],
            relaxation,
        )?;
        self.lock_orientation(bodies, relaxation)
    }

    fn lock_orientation(
        &self,
        bodies: &mut [RigidBody],
        relaxation: f64,
    ) -> Result<(), DeviceError> {
        let (a, b) = two_bodies_mut(bodies, self.body_a, self.body_b)?;
        let relative = a.orientation.conjugate() * b.orientation;
        let world_delta = a.orientation * relative.conjugate() * a.orientation.conjugate();
        let (direction, angle) = world_delta.to_axis_angle();
        if angle <= 1.0e-12 {
            return Ok(());
        }
        let inverse_inertia_a = a.inverse_inertia_world();
        let inverse_inertia_b = b.inverse_inertia_world();
        let weight_a = dot(direction, mat3_mul_vec(inverse_inertia_a, direction));
        let weight_b = dot(direction, mat3_mul_vec(inverse_inertia_b, direction));
        let total = weight_a + weight_b;
        if total <= 1.0e-30 {
            return Ok(());
        }
        let fraction_a = weight_a / total;
        a.apply_rotation_delta(scale3(direction, -fraction_a * angle * relaxation));
        b.apply_rotation_delta(scale3(direction, (1.0 - fraction_a) * angle * relaxation));
        Ok(())
    }
}

/// A gear coupling between two revolute joints.
///
/// The coupling constrains the two joint angles with
/// `radius_a * theta_a + radius_b * theta_b = 0`. For two external gears of
/// pitch radius `r`, this gives the gear ratio `theta_b / theta_a = -r_a / r_b`.
///
/// The coupling is a position-level constraint on the angle. Each projection
/// turns the second body of each revolute joint about its own axis, through the
/// joint anchor. That keeps the revolute constraints satisfied.
#[derive(Clone, Debug, PartialEq)]
pub struct GearCoupling {
    joint_a: usize,
    joint_b: usize,
    radius_a_m: f64,
    radius_b_m: f64,
}

impl GearCoupling {
    /// Creates a gear coupling between two revolute joints.
    ///
    /// The radii must be finite and positive. The two joint indices must
    /// differ.
    pub fn new(
        joint_a: usize,
        joint_b: usize,
        radius_a_m: f64,
        radius_b_m: f64,
    ) -> Result<Self, DeviceError> {
        for radius_m in [radius_a_m, radius_b_m] {
            if !radius_m.is_finite() || radius_m <= 0.0 {
                return Err(DeviceError::InvalidRadius { radius_m });
            }
        }
        if joint_a == joint_b {
            return Err(DeviceError::DuplicateGearJoint { joint: joint_a });
        }
        Ok(Self {
            joint_a,
            joint_b,
            radius_a_m,
            radius_b_m,
        })
    }

    /// Returns the index of the first revolute joint.
    pub fn joint_a(&self) -> usize {
        self.joint_a
    }

    /// Returns the index of the second revolute joint.
    pub fn joint_b(&self) -> usize {
        self.joint_b
    }

    /// Returns the first pitch radius in metres.
    pub fn radius_a_m(&self) -> f64 {
        self.radius_a_m
    }

    /// Returns the second pitch radius in metres.
    pub fn radius_b_m(&self) -> f64 {
        self.radius_b_m
    }

    /// Returns the angle of the first joint in radians.
    pub fn angle_a_rad(
        &self,
        revolute: &[RevoluteJoint],
        bodies: &[RigidBody],
    ) -> Result<f64, DeviceError> {
        self.revolute_at(revolute, self.joint_a)?.angle_rad(bodies)
    }

    /// Returns the angle of the second joint in radians.
    pub fn angle_b_rad(
        &self,
        revolute: &[RevoluteJoint],
        bodies: &[RigidBody],
    ) -> Result<f64, DeviceError> {
        self.revolute_at(revolute, self.joint_b)?.angle_rad(bodies)
    }

    /// Returns the constraint error `r_a * theta_a + r_b * theta_b` in metres.
    pub fn error_m(
        &self,
        revolute: &[RevoluteJoint],
        bodies: &[RigidBody],
    ) -> Result<f64, DeviceError> {
        let theta_a = self.angle_a_rad(revolute, bodies)?;
        let theta_b = self.angle_b_rad(revolute, bodies)?;
        Ok(self.radius_a_m * theta_a + self.radius_b_m * theta_b)
    }

    /// Projects the coupling error to zero.
    pub fn project(
        &self,
        revolute: &[RevoluteJoint],
        bodies: &mut [RigidBody],
        relaxation: f64,
    ) -> Result<(), DeviceError> {
        let error = self.error_m(revolute, bodies)?;
        let denominator = self.radius_a_m * self.radius_a_m + self.radius_b_m * self.radius_b_m;
        let delta_a = -error * self.radius_a_m / denominator * relaxation;
        let delta_b = -error * self.radius_b_m / denominator * relaxation;
        self.rotate_joint_body(revolute, bodies, self.joint_a, delta_a)?;
        self.rotate_joint_body(revolute, bodies, self.joint_b, delta_b)?;
        Ok(())
    }

    fn rotate_joint_body(
        &self,
        revolute: &[RevoluteJoint],
        bodies: &mut [RigidBody],
        joint: usize,
        delta_rad: f64,
    ) -> Result<(), DeviceError> {
        apply_joint_angle_delta(revolute, bodies, joint, delta_rad)
    }

    fn revolute_at<'a>(
        &self,
        revolute: &'a [RevoluteJoint],
        joint: usize,
    ) -> Result<&'a RevoluteJoint, DeviceError> {
        revolute
            .get(joint)
            .ok_or(DeviceError::JointIndexOutOfBounds {
                joint,
                joint_count: revolute.len(),
            })
    }
}

/// One joint term of a [`GearConstraint`]: a revolute joint and its signed
/// coefficient in metres.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GearTerm {
    joint: usize,
    coefficient_m: f64,
}

impl GearTerm {
    /// Builds one signed term. The coefficient must be finite and non-zero.
    pub fn new(joint: usize, coefficient_m: f64) -> Result<Self, DeviceError> {
        if !coefficient_m.is_finite() || coefficient_m == 0.0 {
            return Err(DeviceError::InvalidGearCoefficient {
                joint,
                coefficient_m,
            });
        }
        Ok(Self {
            joint,
            coefficient_m,
        })
    }

    /// Returns the revolute joint index.
    pub fn joint(&self) -> usize {
        self.joint
    }

    /// Returns the signed coefficient in metres.
    pub fn coefficient_m(&self) -> f64 {
        self.coefficient_m
    }
}

/// A general linear gear constraint over revolute joints.
///
/// The constraint is `sum_i c_i * theta_i = 0`, with `theta_i` the angle of
/// revolute joint `i` and `c_i` its signed coefficient. A two-term constraint
/// with opposite signs is the same relation as [`GearCoupling`]. More terms
/// express a planetary gear mesh, where the carrier angle must be included.
///
/// The projection is the same mass-independent least-squares correction used by
/// [`GearCoupling`]. It drives the constraint error to zero. A fixed body never
/// moves: when the second body of a joint is fixed, the method turns the first
/// body by the opposite angle instead.
#[derive(Clone, Debug, PartialEq)]
pub struct GearConstraint {
    terms: Vec<GearTerm>,
}

impl GearConstraint {
    /// Builds a constraint from `(joint, coefficient_m)` pairs.
    ///
    /// The constraint needs at least two terms. A joint must not appear twice,
    /// so each coefficient is unique. A zero or non-finite coefficient is
    /// rejected.
    pub fn new(terms: impl IntoIterator<Item = (usize, f64)>) -> Result<Self, DeviceError> {
        let mut collected: Vec<GearTerm> = Vec::new();
        for (joint, coefficient_m) in terms {
            if collected.iter().any(|term| term.joint == joint) {
                return Err(DeviceError::DuplicateGearConstraintJoint { joint });
            }
            collected.push(GearTerm::new(joint, coefficient_m)?);
        }
        if collected.len() < 2 {
            return Err(DeviceError::TooFewGearTerms {
                terms: collected.len(),
            });
        }
        Ok(Self { terms: collected })
    }

    /// Returns the terms in construction order.
    pub fn terms(&self) -> &[GearTerm] {
        &self.terms
    }

    /// Returns the constraint error `sum_i c_i * theta_i` in metres.
    pub fn error_m(
        &self,
        revolute: &[RevoluteJoint],
        bodies: &[RigidBody],
    ) -> Result<f64, DeviceError> {
        let mut error_m = 0.0;
        for term in &self.terms {
            let joint = revolute_at(revolute, term.joint)?;
            error_m += term.coefficient_m * joint.angle_rad(bodies)?;
        }
        Ok(error_m)
    }

    /// Projects the constraint error to zero in the least-squares sense.
    pub fn project(
        &self,
        revolute: &[RevoluteJoint],
        bodies: &mut [RigidBody],
        relaxation: f64,
    ) -> Result<(), DeviceError> {
        let error_m = self.error_m(revolute, bodies)?;
        let denominator: f64 = self
            .terms
            .iter()
            .map(|term| term.coefficient_m * term.coefficient_m)
            .sum();
        if denominator <= 0.0 {
            return Ok(());
        }
        for term in &self.terms {
            let delta_rad = -term.coefficient_m * error_m / denominator * relaxation;
            apply_joint_angle_delta(revolute, bodies, term.joint, delta_rad)?;
        }
        Ok(())
    }
}

/// Turns the angle of one revolute joint by `delta_rad`.
///
/// A positive delta increases the joint angle. The method turns the second
/// body by `+delta_rad`. When that body is fixed, it turns the first body by
/// `-delta_rad` instead, so a fixed body never moves.
fn apply_joint_angle_delta(
    revolute: &[RevoluteJoint],
    bodies: &mut [RigidBody],
    joint: usize,
    delta_rad: f64,
) -> Result<(), DeviceError> {
    let joint = revolute_at(revolute, joint)?;
    let body_a = joint.body_a();
    let body_b = joint.body_b();
    let body_count = bodies.len();
    let body_b_fixed = bodies
        .get(body_b)
        .ok_or(DeviceError::BodyIndexOutOfBounds {
            body: body_b,
            body_count,
        })?
        .is_fixed();
    let body_a_fixed = bodies
        .get(body_a)
        .ok_or(DeviceError::BodyIndexOutOfBounds {
            body: body_a,
            body_count,
        })?
        .is_fixed();
    let axis = joint.axis_a_world(bodies)?;
    let anchor = joint.anchor_a_world_m(bodies)?;
    if !body_b_fixed {
        bodies[body_b].rotate_about_point(axis, delta_rad, anchor);
    } else if !body_a_fixed {
        bodies[body_a].rotate_about_point(axis, -delta_rad, anchor);
    }
    Ok(())
}

fn revolute_at(revolute: &[RevoluteJoint], joint: usize) -> Result<&RevoluteJoint, DeviceError> {
    revolute
        .get(joint)
        .ok_or(DeviceError::JointIndexOutOfBounds {
            joint,
            joint_count: revolute.len(),
        })
}

/// A device: rigid bodies, joints, and the step loop.
///
/// The system owns the bodies. A joint references bodies by index. A gear
/// coupling references revolute joints by index.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct RigidBodySystem {
    bodies: Vec<RigidBody>,
    revolute_joints: Vec<RevoluteJoint>,
    prismatic_joints: Vec<PrismaticJoint>,
    gear_couplings: Vec<GearCoupling>,
    gear_constraints: Vec<GearConstraint>,
}

impl RigidBodySystem {
    /// Creates an empty system.
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a body and returns its index.
    pub fn add_body(&mut self, body: RigidBody) -> usize {
        self.bodies.push(body);
        self.bodies.len() - 1
    }

    /// Returns the number of bodies.
    pub fn body_count(&self) -> usize {
        self.bodies.len()
    }

    /// Returns the bodies.
    pub fn bodies(&self) -> &[RigidBody] {
        &self.bodies
    }

    /// Returns the bodies for mutation.
    pub fn bodies_mut(&mut self) -> &mut [RigidBody] {
        &mut self.bodies
    }

    /// Returns a body by index.
    pub fn body(&self, index: usize) -> Option<&RigidBody> {
        self.bodies.get(index)
    }

    /// Returns a body by index for mutation.
    pub fn body_mut(&mut self, index: usize) -> Option<&mut RigidBody> {
        self.bodies.get_mut(index)
    }

    /// Adds a revolute joint and returns its index.
    pub fn add_revolute_joint(&mut self, joint: RevoluteJoint) -> usize {
        self.revolute_joints.push(joint);
        self.revolute_joints.len() - 1
    }

    /// Adds a prismatic joint and returns its index.
    pub fn add_prismatic_joint(&mut self, joint: PrismaticJoint) -> usize {
        self.prismatic_joints.push(joint);
        self.prismatic_joints.len() - 1
    }

    /// Adds a gear coupling and returns its index.
    ///
    /// The joint indices must name revolute joints of this system.
    pub fn add_gear_coupling(&mut self, coupling: GearCoupling) -> Result<usize, DeviceError> {
        for joint in [coupling.joint_a(), coupling.joint_b()] {
            if joint >= self.revolute_joints.len() {
                return Err(DeviceError::JointIndexOutOfBounds {
                    joint,
                    joint_count: self.revolute_joints.len(),
                });
            }
        }
        self.gear_couplings.push(coupling);
        Ok(self.gear_couplings.len() - 1)
    }

    /// Adds a general gear constraint and returns its index.
    ///
    /// Every joint index in the constraint must name a revolute joint of this
    /// system.
    pub fn add_gear_constraint(
        &mut self,
        constraint: GearConstraint,
    ) -> Result<usize, DeviceError> {
        for term in constraint.terms() {
            if term.joint() >= self.revolute_joints.len() {
                return Err(DeviceError::JointIndexOutOfBounds {
                    joint: term.joint(),
                    joint_count: self.revolute_joints.len(),
                });
            }
        }
        self.gear_constraints.push(constraint);
        Ok(self.gear_constraints.len() - 1)
    }

    /// Returns the general gear constraints.
    pub fn gear_constraints(&self) -> &[GearConstraint] {
        &self.gear_constraints
    }

    /// Returns the revolute joints.
    pub fn revolute_joints(&self) -> &[RevoluteJoint] {
        &self.revolute_joints
    }

    /// Returns the prismatic joints.
    pub fn prismatic_joints(&self) -> &[PrismaticJoint] {
        &self.prismatic_joints
    }

    /// Returns the gear couplings.
    pub fn gear_couplings(&self) -> &[GearCoupling] {
        &self.gear_couplings
    }

    /// Advances every body and then projects the constraints.
    ///
    /// `iterations` is the number of Gauss-Seidel passes. A larger value gives
    /// a smaller constraint error. A value of 4 to 50 is typical.
    pub fn step(&mut self, dt_s: f64, iterations: usize) -> Result<(), DeviceError> {
        for body in &mut self.bodies {
            body.integrate(dt_s)?;
        }
        for _ in 0..iterations {
            self.project()?;
        }
        Ok(())
    }

    fn project(&mut self) -> Result<(), DeviceError> {
        let relaxation = 1.0;
        for joint in &self.revolute_joints {
            joint.project(&mut self.bodies, relaxation)?;
        }
        for joint in &self.prismatic_joints {
            joint.project(&mut self.bodies, relaxation)?;
        }
        for coupling in &self.gear_couplings {
            coupling.project(&self.revolute_joints, &mut self.bodies, relaxation)?;
        }
        for constraint in &self.gear_constraints {
            constraint.project(&self.revolute_joints, &mut self.bodies, relaxation)?;
        }
        Ok(())
    }

    /// Returns the largest revolute anchor error in metres.
    pub fn max_revolute_anchor_error_m(&self) -> Result<f64, DeviceError> {
        let mut max_error_m = 0.0_f64;
        for joint in &self.revolute_joints {
            max_error_m = max_error_m.max(joint.anchor_error_m(&self.bodies)?);
        }
        Ok(max_error_m)
    }

    /// Returns the largest revolute axis error in radians.
    pub fn max_revolute_axis_error_rad(&self) -> Result<f64, DeviceError> {
        let mut max_error_rad = 0.0_f64;
        for joint in &self.revolute_joints {
            max_error_rad = max_error_rad.max(joint.axis_error_rad(&self.bodies)?);
        }
        Ok(max_error_rad)
    }

    /// Returns the largest prismatic perpendicular position error in metres.
    pub fn max_prismatic_position_error_m(&self) -> Result<f64, DeviceError> {
        let mut max_error_m = 0.0_f64;
        for joint in &self.prismatic_joints {
            max_error_m = max_error_m.max(joint.position_error_m(&self.bodies)?);
        }
        Ok(max_error_m)
    }

    /// Returns the largest prismatic orientation error in radians.
    pub fn max_prismatic_orientation_error_rad(&self) -> Result<f64, DeviceError> {
        let mut max_error_rad = 0.0_f64;
        for joint in &self.prismatic_joints {
            max_error_rad = max_error_rad.max(joint.orientation_error_rad(&self.bodies)?);
        }
        Ok(max_error_rad)
    }

    /// Returns the largest gear coupling error in metres.
    pub fn max_gear_error_m(&self) -> Result<f64, DeviceError> {
        let mut max_error_m = 0.0_f64;
        for coupling in &self.gear_couplings {
            max_error_m = max_error_m.max(coupling.error_m(&self.revolute_joints, &self.bodies)?);
        }
        Ok(max_error_m)
    }

    /// Returns the largest general gear constraint error in metres.
    pub fn max_gear_constraint_error_m(&self) -> Result<f64, DeviceError> {
        let mut max_error_m = 0.0_f64;
        for constraint in &self.gear_constraints {
            max_error_m = max_error_m.max(constraint.error_m(&self.revolute_joints, &self.bodies)?);
        }
        Ok(max_error_m)
    }
}

fn body_at(bodies: &[RigidBody], index: usize) -> Result<&RigidBody, DeviceError> {
    bodies.get(index).ok_or(DeviceError::BodyIndexOutOfBounds {
        body: index,
        body_count: bodies.len(),
    })
}

fn two_bodies_mut(
    bodies: &mut [RigidBody],
    index_a: usize,
    index_b: usize,
) -> Result<(&mut RigidBody, &mut RigidBody), DeviceError> {
    let body_count = bodies.len();
    if index_a >= body_count || index_b >= body_count {
        return Err(DeviceError::BodyIndexOutOfBounds {
            body: index_a.max(index_b),
            body_count,
        });
    }
    if index_a < index_b {
        let (left, right) = bodies.split_at_mut(index_b);
        Ok((&mut left[index_a], &mut right[0]))
    } else {
        let (left, right) = bodies.split_at_mut(index_a);
        Ok((&mut right[0], &mut left[index_b]))
    }
}

fn check_axis(axis: [f64; 3]) -> Result<[f64; 3], DeviceError> {
    if axis.iter().any(|value| !value.is_finite()) {
        return Err(DeviceError::InvalidAxis { axis });
    }
    let norm = norm3(axis);
    if norm <= 1.0e-300 {
        return Err(DeviceError::InvalidAxis { axis });
    }
    Ok(scale3(axis, 1.0 / norm))
}

fn check_anchor(anchor_m: [f64; 3]) -> Result<[f64; 3], DeviceError> {
    if anchor_m.iter().any(|value| !value.is_finite()) {
        return Err(DeviceError::NonFiniteAnchor { anchor_m });
    }
    Ok(anchor_m)
}

fn project_point_rows(
    bodies: &mut [RigidBody],
    body_a: usize,
    local_a_com_m: [f64; 3],
    body_b: usize,
    local_b_com_m: [f64; 3],
    directions: &[[f64; 3]],
    relaxation: f64,
) -> Result<(), DeviceError> {
    let (a, b) = two_bodies_mut(bodies, body_a, body_b)?;
    for direction in directions {
        let r_a = a.orientation.rotate(local_a_com_m);
        let r_b = b.orientation.rotate(local_b_com_m);
        let p_a = add3(a.position_m, r_a);
        let p_b = add3(b.position_m, r_b);
        let error = dot(*direction, sub3(p_b, p_a));
        if error == 0.0 {
            continue;
        }
        let stiffness = point_k_matrix(a, b, r_a, r_b);
        let denominator = dot(*direction, mat3_mul_vec(stiffness, *direction));
        if denominator <= 1.0e-30 {
            continue;
        }
        let impulse = scale3(*direction, error / denominator);
        a.position_m = add3(
            a.position_m,
            scale3(impulse, a.inverse_mass_kg() * relaxation),
        );
        let rotation_a = mat3_mul_vec(a.inverse_inertia_world(), cross(r_a, impulse));
        a.apply_rotation_delta(scale3(rotation_a, relaxation));
        b.position_m = sub3(
            b.position_m,
            scale3(impulse, b.inverse_mass_kg() * relaxation),
        );
        let rotation_b = mat3_mul_vec(b.inverse_inertia_world(), cross(r_b, impulse));
        b.apply_rotation_delta(scale3(rotation_b, -relaxation));
    }
    Ok(())
}

fn point_k_matrix(a: &RigidBody, b: &RigidBody, r_a: [f64; 3], r_b: [f64; 3]) -> [[f64; 3]; 3] {
    let total_inverse_mass = a.inverse_mass_kg() + b.inverse_mass_kg();
    let mut k = scale_mat3(IDENTITY3, total_inverse_mass);
    let term_a = mat3_mul(
        mat3_mul(skew(r_a), a.inverse_inertia_world()),
        mat3_transpose(skew(r_a)),
    );
    let term_b = mat3_mul(
        mat3_mul(skew(r_b), b.inverse_inertia_world()),
        mat3_transpose(skew(r_b)),
    );
    k = add_mat3(k, add_mat3(term_a, term_b));
    k
}

fn perpendicular_basis(axis: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let helper = if axis[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let u = {
        let raw = cross(axis, helper);
        scale3(raw, 1.0 / norm3(raw))
    };
    let v = cross(axis, u);
    (u, v)
}

const IDENTITY3: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scale3(a: [f64; 3], factor: f64) -> [f64; 3] {
    [a[0] * factor, a[1] * factor, a[2] * factor]
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

fn norm3(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn skew(v: [f64; 3]) -> [[f64; 3]; 3] {
    [[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]]
}

fn scale_mat3(m: [[f64; 3]; 3], factor: f64) -> [[f64; 3]; 3] {
    [
        scale3(m[0], factor),
        scale3(m[1], factor),
        scale3(m[2], factor),
    ]
}

fn add_mat3(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [add3(a[0], b[0]), add3(a[1], b[1]), add3(a[2], b[2])]
}

fn mat3_transpose(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]],
    ]
}

fn mat3_mul_vec(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [dot(m[0], v), dot(m[1], v), dot(m[2], v)]
}

fn mat3_mul(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let bt = mat3_transpose(b);
    [
        [dot(a[0], bt[0]), dot(a[0], bt[1]), dot(a[0], bt[2])],
        [dot(a[1], bt[0]), dot(a[1], bt[1]), dot(a[1], bt[2])],
        [dot(a[2], bt[0]), dot(a[2], bt[1]), dot(a[2], bt[2])],
    ]
}

fn mat3_inverse(m: [[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
    let a = m[0][0];
    let b = m[0][1];
    let c = m[0][2];
    let d = m[1][0];
    let e = m[1][1];
    let f = m[1][2];
    let g = m[2][0];
    let h = m[2][1];
    let i = m[2][2];
    let determinant = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
    if !determinant.is_finite() || determinant.abs() < 1.0e-300 {
        return None;
    }
    let inverse_determinant = 1.0 / determinant;
    Some([
        [
            (e * i - f * h) * inverse_determinant,
            (c * h - b * i) * inverse_determinant,
            (b * f - c * e) * inverse_determinant,
        ],
        [
            (f * g - d * i) * inverse_determinant,
            (a * i - c * g) * inverse_determinant,
            (c * d - a * f) * inverse_determinant,
        ],
        [
            (d * h - e * g) * inverse_determinant,
            (b * g - a * h) * inverse_determinant,
            (a * e - b * d) * inverse_determinant,
        ],
    ])
}

fn is_symmetric(m: [[f64; 3]; 3]) -> bool {
    for (row, values) in m.iter().enumerate() {
        for (column, value) in values.iter().enumerate() {
            let other = m[column][row];
            let scale = value.abs().max(other.abs()).max(1.0);
            if (value - other).abs() > 1.0e-12 * scale {
                return false;
            }
        }
    }
    true
}

fn is_positive_definite(m: [[f64; 3]; 3]) -> bool {
    let leading_1 = m[0][0];
    let leading_2 = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    let determinant = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    leading_1 > 0.0 && leading_2 > 0.0 && determinant > 0.0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn diagonal_body(mass_kg: f64, inertia_kg_m2: [f64; 3]) -> RigidBody {
        RigidBody::with_diagonal_inertia(mass_kg, inertia_kg_m2, [0.0; 3], Quat::IDENTITY)
            .expect("valid body")
    }

    fn fixed_ground(position_m: [f64; 3]) -> RigidBody {
        let mut ground = diagonal_body(1000.0, [1000.0, 1000.0, 1000.0]);
        ground.set_position_m(position_m).expect("valid position");
        ground.set_fixed(true);
        ground
    }

    #[test]
    fn quaternion_rotation_round_trips() {
        let rotation = Quat::from_axis_angle([0.0, 0.0, 5.0], 0.7);
        let v = [1.0, 2.0, 3.0];
        let world = rotation.rotate(v);
        let back = rotation.inverse_rotate(world);
        for axis in 0..3 {
            assert!((back[axis] - v[axis]).abs() < 1.0e-12);
        }
        assert!((rotation.norm() - 1.0).abs() < 1.0e-12);
        let (axis, angle) = rotation.to_axis_angle();
        assert!((axis[2] - 1.0).abs() < 1.0e-12);
        assert!((angle - 0.7).abs() < 1.0e-12);
    }

    #[test]
    fn free_translation_conserves_momentum() {
        let mut body = diagonal_body(2.0, [1.0, 1.0, 1.0]);
        body.set_linear_velocity_m_per_s([1.0, 2.0, 3.0])
            .expect("valid velocity");
        let momentum_before = body.momentum_kg_m_per_s();
        let energy_before = body.kinetic_energy_j();
        let dt_s = 1.0e-3;
        for _ in 0..10_000 {
            body.integrate(dt_s).expect("valid step");
        }
        let momentum_after = body.momentum_kg_m_per_s();
        for axis in 0..3 {
            assert!((momentum_after[axis] - momentum_before[axis]).abs() < 1.0e-9);
        }
        assert!((body.kinetic_energy_j() - energy_before).abs() < 1.0e-9);
        let expected_x = 1.0 * 10_000.0 * dt_s;
        assert!((body.position_m()[0] - expected_x).abs() < 1.0e-9);
    }

    #[test]
    fn free_rotation_conserves_angular_momentum() {
        let mut body = diagonal_body(1.0, [1.0, 2.0, 3.0]);
        body.set_angular_velocity_rad_per_s([1.0, 0.5, 0.2])
            .expect("valid velocity");
        let momentum_before = body.angular_momentum_kg_m2_per_s();
        let dt_s = 1.0e-4;
        for _ in 0..10_000 {
            body.integrate(dt_s).expect("valid step");
        }
        let momentum_after = body.angular_momentum_kg_m2_per_s();
        for axis in 0..3 {
            let scale = momentum_before[axis].abs().max(1.0);
            assert!(
                (momentum_after[axis] - momentum_before[axis]).abs() < 1.0e-10 * scale,
                "axis {axis}: before {} after {}",
                momentum_before[axis],
                momentum_after[axis]
            );
        }
        assert!((body.orientation().norm() - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn free_rotation_is_not_constant_for_a_triaxial_body() {
        let mut body = diagonal_body(1.0, [1.0, 2.0, 3.0]);
        body.set_angular_velocity_rad_per_s([1.0, 0.5, 0.2])
            .expect("valid velocity");
        let velocity_before = body.angular_velocity_rad_per_s();
        let dt_s = 1.0e-4;
        for _ in 0..10_000 {
            body.integrate(dt_s).expect("valid step");
        }
        let velocity_after = body.angular_velocity_rad_per_s();
        let changed =
            (0..3).any(|axis| (velocity_after[axis] - velocity_before[axis]).abs() > 1.0e-3);
        assert!(changed, "the angular velocity must precess");
    }

    #[test]
    fn a_fixed_body_does_not_move() {
        let mut body = diagonal_body(1.0, [1.0, 1.0, 1.0]);
        body.set_position_m([1.0, 2.0, 3.0]).expect("valid");
        body.set_fixed(true);
        body.add_force_n([10.0, 20.0, 30.0]);
        body.add_torque_n_m([1.0, 2.0, 3.0]);
        body.integrate(0.1).expect("valid step");
        assert_eq!(body.position_m(), [1.0, 2.0, 3.0]);
        assert_eq!(body.orientation(), Quat::IDENTITY);
        assert_eq!(body.inverse_mass_kg(), 0.0);
    }

    #[test]
    fn the_center_of_mass_offset_maps_a_local_point() {
        let body = RigidBody::new(
            1.0,
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [0.5, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            Quat::from_axis_angle([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_2),
        )
        .expect("valid body");
        let world = body.point_world_m([0.5, 0.0, 0.0]);
        assert!((world[0] - 1.0).abs() < 1.0e-12);
        assert!((world[1] - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn invalid_bodies_are_rejected() {
        assert!(matches!(
            RigidBody::new(
                0.0,
                [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                [0.0; 3],
                [0.0; 3],
                Quat::IDENTITY
            ),
            Err(DeviceError::InvalidMass { .. })
        ));
        assert!(matches!(
            RigidBody::new(
                1.0,
                [[1.0, 0.5, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                [0.0; 3],
                [0.0; 3],
                Quat::IDENTITY
            ),
            Err(DeviceError::NonSymmetricInertia { .. })
        ));
        assert!(matches!(
            RigidBody::new(
                1.0,
                [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]],
                [0.0; 3],
                [0.0; 3],
                Quat::IDENTITY
            ),
            Err(DeviceError::SingularInertia { .. })
        ));
        assert!(matches!(
            RigidBody::with_diagonal_inertia(
                1.0,
                [1.0, 1.0, 1.0],
                [f64::NAN, 0.0, 0.0],
                Quat::IDENTITY
            ),
            Err(DeviceError::NonFinitePosition { .. })
        ));
        assert!(matches!(
            diagonal_body(1.0, [1.0, 1.0, 1.0]).set_angular_velocity_rad_per_s([
                f64::NAN,
                0.0,
                0.0
            ]),
            Err(DeviceError::NonFiniteAngularVelocity { .. })
        ));
    }

    #[test]
    fn a_revolute_joint_starts_satisfied_and_rejects_bad_setup() {
        let ground = fixed_ground([0.0, 0.0, 0.0]);
        let rotor = diagonal_body(1.0, [1.0, 1.0, 1.0]);
        let bodies = vec![ground, rotor];
        let joint = RevoluteJoint::from_world(&bodies, 0, 1, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0])
            .expect("valid joint");
        assert!(joint.anchor_error_m(&bodies).expect("valid") < 1.0e-15);
        assert!(joint.axis_error_rad(&bodies).expect("valid") < 1.0e-15);
        assert!(matches!(
            RevoluteJoint::from_world(&bodies, 0, 0, [0.0; 3], [0.0, 0.0, 1.0]),
            Err(DeviceError::SelfJoint { .. })
        ));
        assert!(matches!(
            RevoluteJoint::from_world(&bodies, 0, 1, [0.0; 3], [0.0, 0.0, 0.0]),
            Err(DeviceError::InvalidAxis { .. })
        ));
        assert!(matches!(
            RevoluteJoint::from_world(&bodies, 0, 5, [0.0; 3], [0.0, 0.0, 1.0]),
            Err(DeviceError::BodyIndexOutOfBounds { .. })
        ));
    }

    #[test]
    fn a_revolute_anchor_error_stays_small_over_many_steps() {
        let ground = fixed_ground([0.0, 0.0, 0.0]);
        let mut rotor = diagonal_body(1.0, [0.5, 0.5, 0.5]);
        rotor
            .set_position_m([0.1, 0.0, 0.0])
            .expect("valid position");
        rotor
            .set_linear_velocity_m_per_s([1.0, 2.0, 0.5])
            .expect("valid velocity");
        rotor
            .set_angular_velocity_rad_per_s([0.3, -0.2, 2.0])
            .expect("valid velocity");
        let joint = RevoluteJoint::from_world(
            &[ground.clone(), rotor.clone()],
            0,
            1,
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        )
        .expect("valid joint");
        let mut bodies = vec![ground, rotor];
        let dt_s = 1.0e-3;
        let mut max_anchor_error_m = 0.0_f64;
        let mut max_axis_error_rad = 0.0_f64;
        for _ in 0..2_000 {
            for body in &mut bodies {
                body.integrate(dt_s).expect("valid step");
            }
            for _ in 0..32 {
                joint.project(&mut bodies, 1.0).expect("valid projection");
            }
            max_anchor_error_m =
                max_anchor_error_m.max(joint.anchor_error_m(&bodies).expect("valid"));
            max_axis_error_rad =
                max_axis_error_rad.max(joint.axis_error_rad(&bodies).expect("valid"));
        }
        assert!(
            max_anchor_error_m < 1.0e-9,
            "maximum anchor error {max_anchor_error_m} m"
        );
        assert!(
            max_axis_error_rad < 1.0e-9,
            "maximum axis error {max_axis_error_rad} rad"
        );
    }

    #[test]
    fn a_revolute_joint_reports_the_relative_angle() {
        let ground = fixed_ground([0.0, 0.0, 0.0]);
        let mut rotor = diagonal_body(1.0, [1.0, 1.0, 1.0]);
        rotor
            .set_orientation(Quat::from_axis_angle([0.0, 0.0, 1.0], 0.6))
            .expect("valid orientation");
        let bodies = vec![ground, rotor];
        let joint = RevoluteJoint::from_world(&bodies, 0, 1, [0.0; 3], [0.0, 0.0, 1.0])
            .expect("valid joint");
        assert!((joint.angle_rad(&bodies).expect("valid") - 0.6).abs() < 1.0e-12);
    }

    #[test]
    fn a_prismatic_joint_constrains_the_perpendicular_motion() {
        let ground = fixed_ground([0.0, 0.0, 0.0]);
        let mut slider = diagonal_body(1.0, [1.0, 1.0, 1.0]);
        slider
            .set_linear_velocity_m_per_s([0.5, 0.3, 0.2])
            .expect("valid velocity");
        slider
            .set_angular_velocity_rad_per_s([0.0, 0.4, 0.0])
            .expect("valid velocity");
        let joint = PrismaticJoint::from_world(
            &[ground.clone(), slider.clone()],
            0,
            1,
            [0.0; 3],
            [1.0, 0.0, 0.0],
        )
        .expect("valid joint");
        let mut bodies = vec![ground, slider];
        let dt_s = 1.0e-3;
        let mut max_position_error_m = 0.0_f64;
        let mut max_orientation_error_rad = 0.0_f64;
        for _ in 0..2_000 {
            for body in &mut bodies {
                body.integrate(dt_s).expect("valid step");
            }
            for _ in 0..32 {
                joint.project(&mut bodies, 1.0).expect("valid projection");
            }
            max_position_error_m =
                max_position_error_m.max(joint.position_error_m(&bodies).expect("valid"));
            max_orientation_error_rad =
                max_orientation_error_rad.max(joint.orientation_error_rad(&bodies).expect("valid"));
        }
        assert!(
            max_position_error_m < 1.0e-9,
            "maximum position error {max_position_error_m} m"
        );
        assert!(
            max_orientation_error_rad < 1.0e-9,
            "maximum orientation error {max_orientation_error_rad} rad"
        );
        let translation_m = joint.translation_m(&bodies).expect("valid");
        assert!((translation_m - 0.5 * 2.0).abs() < 1.0e-6);
    }

    #[test]
    fn a_gear_coupling_holds_the_ratio_over_many_steps() {
        let ground_a = fixed_ground([0.0, 0.0, 0.0]);
        let ground_b = fixed_ground([0.2, 0.0, 0.0]);
        let mut gear_a = diagonal_body(1.0, [0.005, 0.005, 0.005]);
        let gear_b = diagonal_body(1.0, [0.005, 0.005, 0.005]);
        gear_a
            .set_angular_velocity_rad_per_s([0.0, 0.0, 2.0])
            .expect("valid velocity");
        let setup = vec![
            ground_a.clone(),
            gear_a.clone(),
            ground_b.clone(),
            gear_b.clone(),
        ];
        let joint_a = RevoluteJoint::from_world(&setup, 0, 1, [0.0; 3], [0.0, 0.0, 1.0])
            .expect("valid joint");
        let joint_b = RevoluteJoint::from_world(&setup, 2, 3, [0.2, 0.0, 0.0], [0.0, 0.0, 1.0])
            .expect("valid joint");
        let revolute = vec![joint_a.clone(), joint_b.clone()];
        let radius_a_m = 0.05;
        let radius_b_m = 0.10;
        let coupling = GearCoupling::new(0, 1, radius_a_m, radius_b_m).expect("valid coupling");
        let mut bodies = vec![ground_a, gear_a, ground_b, gear_b];
        let dt_s = 1.0e-3;
        let mut max_error_m = 0.0_f64;
        for _ in 0..2_000 {
            for body in &mut bodies {
                body.integrate(dt_s).expect("valid step");
            }
            for _ in 0..32 {
                for joint in &revolute {
                    joint.project(&mut bodies, 1.0).expect("valid projection");
                }
                coupling
                    .project(&revolute, &mut bodies, 1.0)
                    .expect("valid projection");
            }
            max_error_m = max_error_m.max(coupling.error_m(&revolute, &bodies).expect("valid"));
        }
        assert!(
            max_error_m < 1.0e-9,
            "maximum gear constraint error {max_error_m} m"
        );
        let theta_a = coupling.angle_a_rad(&revolute, &bodies).expect("valid");
        let theta_b = coupling.angle_b_rad(&revolute, &bodies).expect("valid");
        let expected_b = -(radius_a_m / radius_b_m) * theta_a;
        assert!(
            (theta_b - expected_b).abs() < 1.0e-6,
            "theta_a {theta_a} theta_b {theta_b} expected {expected_b}"
        );
    }

    #[test]
    fn a_gear_coupling_rejects_bad_setup() {
        assert!(matches!(
            GearCoupling::new(0, 1, -1.0, 1.0),
            Err(DeviceError::InvalidRadius { .. })
        ));
        assert!(matches!(
            GearCoupling::new(0, 1, 1.0, f64::NAN),
            Err(DeviceError::InvalidRadius { .. })
        ));
        assert!(matches!(
            GearCoupling::new(2, 2, 0.05, 0.1),
            Err(DeviceError::DuplicateGearJoint { .. })
        ));
    }

    #[test]
    fn a_system_runs_the_whole_step_loop() {
        let ground = fixed_ground([0.0, 0.0, 0.0]);
        let mut rotor = diagonal_body(1.0, [1.0, 1.0, 1.0]);
        rotor
            .set_linear_velocity_m_per_s([0.2, -0.1, 0.0])
            .expect("valid velocity");
        let mut system = RigidBodySystem::new();
        assert_eq!(system.add_body(ground), 0);
        assert_eq!(system.add_body(rotor), 1);
        let joint = RevoluteJoint::from_world(system.bodies(), 0, 1, [0.0; 3], [0.0, 0.0, 1.0])
            .expect("valid joint");
        system.add_revolute_joint(joint);
        for _ in 0..1_000 {
            system.step(1.0e-3, 32).expect("valid step");
        }
        assert!(system.max_revolute_anchor_error_m().expect("valid") < 1.0e-9);
        assert!(system.max_revolute_axis_error_rad().expect("valid") < 1.0e-9);
        assert_eq!(system.body_count(), 2);
    }

    #[test]
    fn adding_a_gear_coupling_checks_the_joint_index() {
        let mut system = RigidBodySystem::new();
        system.add_body(fixed_ground([0.0; 3]));
        system.add_body(diagonal_body(1.0, [1.0, 1.0, 1.0]));
        let coupling = GearCoupling::new(0, 1, 0.05, 0.1).expect("valid coupling");
        assert!(matches!(
            system.add_gear_coupling(coupling),
            Err(DeviceError::JointIndexOutOfBounds { .. })
        ));
    }

    #[test]
    fn a_three_term_gear_constraint_holds_over_many_steps() {
        let ground = fixed_ground([0.0, 0.0, 0.0]);
        let mut wheel_a = diagonal_body(1.0, [0.005, 0.005, 0.005]);
        let wheel_b = diagonal_body(1.0, [0.005, 0.005, 0.005]);
        wheel_a
            .set_angular_velocity_rad_per_s([0.0, 0.0, 2.0])
            .expect("valid velocity");
        let setup = vec![ground.clone(), wheel_a.clone(), wheel_b.clone()];
        let joint_a = RevoluteJoint::from_world(&setup, 0, 1, [0.0; 3], [0.0, 0.0, 1.0])
            .expect("valid joint");
        let joint_b = RevoluteJoint::from_world(&setup, 0, 2, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0])
            .expect("valid joint");
        let revolute = vec![joint_a, joint_b];
        let constraint = GearConstraint::new([(0, 0.05), (1, 0.025)]).expect("valid constraint");
        let mut bodies = vec![ground, wheel_a, wheel_b];
        let dt_s = 1.0e-3;
        let mut max_error_m = 0.0_f64;
        for _ in 0..2_000 {
            for body in &mut bodies {
                body.integrate(dt_s).expect("valid step");
            }
            for _ in 0..32 {
                for joint in &revolute {
                    joint.project(&mut bodies, 1.0).expect("valid projection");
                }
                constraint
                    .project(&revolute, &mut bodies, 1.0)
                    .expect("valid projection");
            }
            max_error_m = max_error_m.max(constraint.error_m(&revolute, &bodies).expect("valid"));
        }
        assert!(
            max_error_m < 1.0e-12,
            "maximum gear constraint error {max_error_m} m"
        );
        let theta_a = revolute[0].angle_rad(&bodies).expect("valid");
        let theta_b = revolute[1].angle_rad(&bodies).expect("valid");
        let expected_b = -(0.05 / 0.025) * theta_a;
        assert!((theta_b - expected_b).abs() < 1.0e-9);
    }

    #[test]
    fn a_gear_constraint_leaves_a_fixed_second_body_alone() {
        let ground = fixed_ground([0.0, 0.0, 0.0]);
        let mut driver = diagonal_body(1.0, [0.005, 0.005, 0.005]);
        driver
            .set_orientation(Quat::from_axis_angle([0.0, 0.0, 1.0], 0.1))
            .expect("valid orientation");
        let other_ground = fixed_ground([0.3, 0.0, 0.0]);
        let mut anchored = diagonal_body(1.0, [0.005, 0.005, 0.005]);
        anchored.set_fixed(true);
        let setup = vec![
            ground.clone(),
            driver.clone(),
            other_ground.clone(),
            anchored.clone(),
        ];
        let joint_driver = RevoluteJoint::from_world(&setup, 0, 1, [0.0; 3], [0.0, 0.0, 1.0])
            .expect("valid joint");
        let joint_fixed = RevoluteJoint::from_world(&setup, 2, 3, [0.3, 0.0, 0.0], [0.0, 0.0, 1.0])
            .expect("valid joint");
        let revolute = vec![joint_driver, joint_fixed];
        let constraint = GearConstraint::new([(0, 0.05), (1, 0.05)]).expect("valid constraint");
        let mut bodies = vec![ground, driver, other_ground, anchored];
        let fixed_before = bodies[3].clone();
        for _ in 0..10 {
            constraint
                .project(&revolute, &mut bodies, 1.0)
                .expect("valid projection");
        }
        assert_eq!(bodies[3].position_m(), fixed_before.position_m());
        assert_eq!(bodies[3].orientation(), fixed_before.orientation());
        assert!(revolute[0].angle_rad(&bodies).expect("valid").abs() > 0.0);
    }

    #[test]
    fn a_gear_constraint_rejects_bad_terms() {
        assert!(matches!(
            GearConstraint::new([(0, 1.0)]),
            Err(DeviceError::TooFewGearTerms { terms: 1 })
        ));
        assert!(matches!(
            GearConstraint::new([(0, 1.0), (0, 2.0)]),
            Err(DeviceError::DuplicateGearConstraintJoint { joint: 0 })
        ));
        assert!(matches!(
            GearConstraint::new([(0, 1.0), (1, 0.0)]),
            Err(DeviceError::InvalidGearCoefficient { joint: 1, .. })
        ));
        assert!(matches!(
            GearConstraint::new([(0, 1.0), (1, f64::NAN)]),
            Err(DeviceError::InvalidGearCoefficient { joint: 1, .. })
        ));
    }
}
