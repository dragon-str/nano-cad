//! The L2 device layer: rigid bodies and joints.

use nanocad_jigs as jigs;
use pyo3::prelude::*;

use super::util;

fn quat(values: Vec<f64>) -> PyResult<jigs::Quat> {
    let [w, x, y, z] = util::vec4(values)?;
    Ok(jigs::Quat::new(w, x, y, z))
}

fn mat3_rows(matrix: [[f64; 3]; 3]) -> Vec<Vec<f64>> {
    matrix.iter().map(|row| row.to_vec()).collect()
}

/// A unit quaternion rotation.
#[pyclass(name = "Quat", skip_from_py_object)]
#[derive(Clone)]
pub struct Quat {
    inner: jigs::Quat,
}

#[pymethods]
impl Quat {
    #[new]
    fn new(w: f64, x: f64, y: f64, z: f64) -> Self {
        Self {
            inner: jigs::Quat::new(w, x, y, z),
        }
    }

    /// The identity rotation.
    #[staticmethod]
    fn identity() -> Self {
        Self {
            inner: jigs::Quat::identity(),
        }
    }

    /// A rotation of `angle_rad` about a normalised axis.
    #[staticmethod]
    fn from_axis_angle(axis: Vec<f64>, angle_rad: f64) -> PyResult<Self> {
        Ok(Self {
            inner: jigs::Quat::from_axis_angle(util::vec3(axis)?, angle_rad),
        })
    }

    #[getter]
    fn w(&self) -> f64 {
        self.inner.w
    }

    #[getter]
    fn x(&self) -> f64 {
        self.inner.x
    }

    #[getter]
    fn y(&self) -> f64 {
        self.inner.y
    }

    #[getter]
    fn z(&self) -> f64 {
        self.inner.z
    }

    #[getter]
    fn norm(&self) -> f64 {
        self.inner.norm()
    }

    /// Returns a unit-length copy.
    fn normalized(&self) -> Self {
        Self {
            inner: self.inner.normalized(),
        }
    }

    /// Returns the inverse rotation.
    fn conjugate(&self) -> Self {
        Self {
            inner: self.inner.conjugate(),
        }
    }

    /// Rotates a vector.
    fn rotate(&self, vector: Vec<f64>) -> PyResult<(f64, f64, f64)> {
        let out = self.inner.rotate(util::vec3(vector)?);
        Ok((out[0], out[1], out[2]))
    }

    /// Rotates a vector by the inverse rotation.
    fn inverse_rotate(&self, vector: Vec<f64>) -> PyResult<(f64, f64, f64)> {
        let out = self.inner.inverse_rotate(util::vec3(vector)?);
        Ok((out[0], out[1], out[2]))
    }

    /// Returns `(axis, angle_rad)`.
    fn to_axis_angle(&self) -> ([f64; 3], f64) {
        self.inner.to_axis_angle()
    }

    /// Returns the 3x3 rotation matrix.
    fn to_rotation_matrix(&self) -> Vec<Vec<f64>> {
        mat3_rows(self.inner.to_rotation_matrix())
    }
}

/// A rigid body with mass, inertia, and state.
#[pyclass(name = "RigidBody", skip_from_py_object)]
#[derive(Clone)]
pub struct RigidBody {
    inner: jigs::RigidBody,
}

#[pymethods]
impl RigidBody {
    #[new]
    #[pyo3(signature = (mass_kg, inertia_kg_m2, position_m, center_of_mass_m=None, orientation=None))]
    fn new(
        mass_kg: f64,
        inertia_kg_m2: Vec<Vec<f64>>,
        position_m: Vec<f64>,
        center_of_mass_m: Option<Vec<f64>>,
        orientation: Option<Vec<f64>>,
    ) -> PyResult<Self> {
        let center_of_mass_m = match center_of_mass_m {
            Some(values) => util::vec3(values)?,
            None => [0.0; 3],
        };
        let orientation = match orientation {
            Some(values) => quat(values)?,
            None => jigs::Quat::identity(),
        };
        Ok(Self {
            inner: jigs::RigidBody::new(
                mass_kg,
                util::mat3(inertia_kg_m2)?,
                center_of_mass_m,
                util::vec3(position_m)?,
                orientation,
            )
            .map_err(util::err)?,
        })
    }

    /// Builds a body with a diagonal inertia tensor.
    #[staticmethod]
    #[pyo3(signature = (mass_kg, inertia_diagonal_kg_m2, position_m, orientation=None))]
    fn diagonal(
        mass_kg: f64,
        inertia_diagonal_kg_m2: Vec<f64>,
        position_m: Vec<f64>,
        orientation: Option<Vec<f64>>,
    ) -> PyResult<Self> {
        let orientation = match orientation {
            Some(values) => quat(values)?,
            None => jigs::Quat::identity(),
        };
        Ok(Self {
            inner: jigs::RigidBody::with_diagonal_inertia(
                mass_kg,
                util::vec3(inertia_diagonal_kg_m2)?,
                util::vec3(position_m)?,
                orientation,
            )
            .map_err(util::err)?,
        })
    }

    #[getter]
    fn mass_kg(&self) -> f64 {
        self.inner.mass_kg()
    }

    #[getter]
    fn is_fixed(&self) -> bool {
        self.inner.is_fixed()
    }

    #[setter]
    fn set_is_fixed(&mut self, fixed: bool) {
        self.inner.set_fixed(fixed);
    }

    #[getter]
    fn center_of_mass_m(&self) -> [f64; 3] {
        self.inner.center_of_mass_m()
    }

    #[getter]
    fn position_m(&self) -> [f64; 3] {
        self.inner.position_m()
    }

    #[setter]
    fn set_position_m(&mut self, position_m: Vec<f64>) -> PyResult<()> {
        self.inner
            .set_position_m(util::vec3(position_m)?)
            .map_err(util::err)
    }

    #[getter]
    fn orientation(&self) -> Quat {
        Quat {
            inner: self.inner.orientation(),
        }
    }

    #[getter]
    fn linear_velocity_m_per_s(&self) -> [f64; 3] {
        self.inner.linear_velocity_m_per_s()
    }

    #[setter]
    fn set_linear_velocity_m_per_s(&mut self, velocity: Vec<f64>) -> PyResult<()> {
        self.inner
            .set_linear_velocity_m_per_s(util::vec3(velocity)?)
            .map_err(util::err)
    }

    #[getter]
    fn angular_velocity_rad_per_s(&self) -> [f64; 3] {
        self.inner.angular_velocity_rad_per_s()
    }

    #[setter]
    fn set_angular_velocity_rad_per_s(&mut self, velocity: Vec<f64>) -> PyResult<()> {
        self.inner
            .set_angular_velocity_rad_per_s(util::vec3(velocity)?)
            .map_err(util::err)
    }

    /// Sets the orientation from an axis and an angle, in radians.
    fn set_orientation_axis_angle(&mut self, axis: Vec<f64>, angle_rad: f64) -> PyResult<()> {
        self.inner
            .set_orientation(jigs::Quat::from_axis_angle(util::vec3(axis)?, angle_rad))
            .map_err(util::err)
    }

    #[getter]
    fn inertia_kg_m2(&self) -> Vec<Vec<f64>> {
        mat3_rows(self.inner.inertia_kg_m2())
    }

    #[getter]
    fn inertia_world_kg_m2(&self) -> Vec<Vec<f64>> {
        mat3_rows(self.inner.inertia_world_kg_m2())
    }

    #[getter]
    fn force_n(&self) -> [f64; 3] {
        self.inner.force_n()
    }

    #[getter]
    fn torque_n_m(&self) -> [f64; 3] {
        self.inner.torque_n_m()
    }

    /// Adds a force in newtons.
    fn add_force_n(&mut self, force_n: Vec<f64>) -> PyResult<()> {
        self.inner.add_force_n(util::vec3(force_n)?);
        Ok(())
    }

    /// Adds a torque in newton metres.
    fn add_torque_n_m(&mut self, torque_n_m: Vec<f64>) -> PyResult<()> {
        self.inner.add_torque_n_m(util::vec3(torque_n_m)?);
        Ok(())
    }

    /// Adds a force at a world point.
    fn add_force_at_point_n(&mut self, force_n: Vec<f64>, point_m: Vec<f64>) -> PyResult<()> {
        self.inner
            .add_force_at_point_n(util::vec3(force_n)?, util::vec3(point_m)?);
        Ok(())
    }

    /// Clears the accumulated force and torque.
    fn clear_forces(&mut self) {
        self.inner.clear_forces();
    }

    /// Returns the momentum in kilogram metres per second.
    fn momentum(&self) -> [f64; 3] {
        self.inner.momentum_kg_m_per_s()
    }

    /// Returns the kinetic energy in joules.
    fn kinetic_energy(&self) -> f64 {
        self.inner.kinetic_energy_j()
    }

    /// Integrates the body by one step.
    fn integrate(&mut self, dt_s: f64) -> PyResult<()> {
        self.inner.integrate(dt_s).map_err(util::err)
    }
}

/// A system of rigid bodies with joints.
#[pyclass(name = "RigidBodySystem", skip_from_py_object)]
#[derive(Clone)]
pub struct RigidBodySystem {
    inner: jigs::RigidBodySystem,
}

#[pymethods]
impl RigidBodySystem {
    #[new]
    fn new() -> Self {
        Self {
            inner: jigs::RigidBodySystem::new(),
        }
    }

    #[getter]
    fn body_count(&self) -> usize {
        self.inner.body_count()
    }

    #[getter]
    fn revolute_joint_count(&self) -> usize {
        self.inner.revolute_joints().len()
    }

    #[getter]
    fn prismatic_joint_count(&self) -> usize {
        self.inner.prismatic_joints().len()
    }

    #[getter]
    fn gear_coupling_count(&self) -> usize {
        self.inner.gear_couplings().len()
    }

    /// Adds a body and returns its index.
    fn add_body(&mut self, body: PyRef<'_, RigidBody>) -> usize {
        self.inner.add_body(body.inner.clone())
    }

    /// Returns a body by index, or `None`.
    fn body(&self, index: usize) -> Option<RigidBody> {
        self.inner.body(index).map(|inner| RigidBody {
            inner: inner.clone(),
        })
    }

    /// Sets a body position by index.
    fn set_body_position_m(&mut self, index: usize, position_m: Vec<f64>) -> PyResult<()> {
        let position_m = util::vec3(position_m)?;
        let body = self
            .inner
            .body_mut(index)
            .ok_or_else(|| pyo3::exceptions::PyIndexError::new_err("body index out of range"))?;
        body.set_position_m(position_m).map_err(util::err)
    }

    /// Adds a revolute joint that starts satisfied.
    fn add_revolute_joint(
        &mut self,
        body_a: usize,
        body_b: usize,
        anchor_m: Vec<f64>,
        axis: Vec<f64>,
    ) -> PyResult<usize> {
        let joint = jigs::RevoluteJoint::from_world(
            self.inner.bodies(),
            body_a,
            body_b,
            util::vec3(anchor_m)?,
            util::vec3(axis)?,
        )
        .map_err(util::err)?;
        Ok(self.inner.add_revolute_joint(joint))
    }

    /// Adds a prismatic joint that starts satisfied.
    fn add_prismatic_joint(
        &mut self,
        body_a: usize,
        body_b: usize,
        anchor_m: Vec<f64>,
        axis: Vec<f64>,
    ) -> PyResult<usize> {
        let joint = jigs::PrismaticJoint::from_world(
            self.inner.bodies(),
            body_a,
            body_b,
            util::vec3(anchor_m)?,
            util::vec3(axis)?,
        )
        .map_err(util::err)?;
        Ok(self.inner.add_prismatic_joint(joint))
    }

    /// Adds a gear coupling between two revolute joints.
    fn add_gear_coupling(
        &mut self,
        joint_a: usize,
        joint_b: usize,
        radius_a_m: f64,
        radius_b_m: f64,
    ) -> PyResult<usize> {
        let coupling =
            jigs::GearCoupling::new(joint_a, joint_b, radius_a_m, radius_b_m).map_err(util::err)?;
        self.inner.add_gear_coupling(coupling).map_err(util::err)
    }

    /// Advances the system by one step with `iterations` constraint passes.
    fn step(&mut self, dt_s: f64, iterations: usize) -> PyResult<()> {
        self.inner.step(dt_s, iterations).map_err(util::err)
    }

    /// Returns the angle of a revolute joint, in radians.
    fn revolute_angle_rad(&self, joint: usize) -> PyResult<f64> {
        let joint_ref =
            self.inner.revolute_joints().get(joint).ok_or_else(|| {
                pyo3::exceptions::PyIndexError::new_err("joint index out of range")
            })?;
        joint_ref.angle_rad(self.inner.bodies()).map_err(util::err)
    }

    /// Returns the translation of a prismatic joint, in metres.
    fn prismatic_translation_m(&self, joint: usize) -> PyResult<f64> {
        let joint_ref =
            self.inner.prismatic_joints().get(joint).ok_or_else(|| {
                pyo3::exceptions::PyIndexError::new_err("joint index out of range")
            })?;
        joint_ref
            .translation_m(self.inner.bodies())
            .map_err(util::err)
    }

    /// Returns the largest revolute anchor error in metres.
    fn max_revolute_anchor_error_m(&self) -> PyResult<f64> {
        self.inner.max_revolute_anchor_error_m().map_err(util::err)
    }

    /// Returns the largest revolute axis error in radians.
    fn max_revolute_axis_error_rad(&self) -> PyResult<f64> {
        self.inner.max_revolute_axis_error_rad().map_err(util::err)
    }

    /// Returns the largest prismatic position error in metres.
    fn max_prismatic_position_error_m(&self) -> PyResult<f64> {
        self.inner
            .max_prismatic_position_error_m()
            .map_err(util::err)
    }

    /// Returns the largest prismatic orientation error in radians.
    fn max_prismatic_orientation_error_rad(&self) -> PyResult<f64> {
        self.inner
            .max_prismatic_orientation_error_rad()
            .map_err(util::err)
    }

    /// Returns the largest gear constraint error in metres.
    fn max_gear_error_m(&self) -> PyResult<f64> {
        self.inner.max_gear_error_m().map_err(util::err)
    }

    /// Returns every body position as a flat buffer.
    fn positions_m(&self) -> Vec<f64> {
        self.inner
            .bodies()
            .iter()
            .flat_map(|body| body.position_m())
            .collect()
    }
}

/// Adds a device constraint by kind: `"revolute"` or `"prismatic"`.
#[pyfunction]
#[pyo3(signature = (device, kind, body_a, body_b, anchor_m, axis))]
pub fn add_constraint(
    mut device: PyRefMut<'_, RigidBodySystem>,
    kind: &str,
    body_a: usize,
    body_b: usize,
    anchor_m: Vec<f64>,
    axis: Vec<f64>,
) -> PyResult<usize> {
    let anchor_m = util::vec3(anchor_m)?;
    let axis = util::vec3(axis)?;
    match kind {
        "revolute" => {
            let joint = jigs::RevoluteJoint::from_world(
                device.inner.bodies(),
                body_a,
                body_b,
                anchor_m,
                axis,
            )
            .map_err(util::err)?;
            Ok(device.inner.add_revolute_joint(joint))
        }
        "prismatic" => {
            let joint = jigs::PrismaticJoint::from_world(
                device.inner.bodies(),
                body_a,
                body_b,
                anchor_m,
                axis,
            )
            .map_err(util::err)?;
            Ok(device.inner.add_prismatic_joint(joint))
        }
        other => Err(pyo3::exceptions::PyValueError::new_err(format!(
            "unknown constraint kind {other:?}; expected revolute or prismatic"
        ))),
    }
}
