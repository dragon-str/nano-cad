//! A rotor test machine: a bead-and-spring ring driven by a [`MotorJig`] under
//! a Langevin thermostat.
//!
//! The test machine is the M3-07 rotor. The ring atoms are point masses in SI
//! kilograms. Ring bonds and radial springs hold the ring shape. The motor
//! supplies a torque about the axis with the existing [`MotorJig`] force
//! distribution. A rotor-scale Langevin coupling supplies the 300 K
//! fluctuations.
//!
//! The engine crate owns the cell list and the production thermostat. This
//! module holds the small rotor-scale coupling that the device test needs, so
//! the test does not depend on unfinished engine work. See the module docs of
//! `device` for the L2 layer.
//!
//! The thermostat is an exact Ornstein-Uhlenbeck update on each velocity
//! component. It is stable for any `gamma * dt`, so the test machine runs many
//! frames without a blow-up.

use nanocad_model::{Atom, Element, Topology};

use crate::device::DeviceError;
use crate::error::JigError;
use crate::jig::Jig;
use crate::motor::MotorJig;
use crate::spring::SpringJig;
use thiserror::Error;

/// The Boltzmann constant in joules per kelvin. This is the exact 2019 SI
/// value. It is defined here because `nanocad-units` holds units, not physics
/// constants.
pub const BOLTZMANN_J_PER_K: f64 = 1.380_649e-23;

/// The mass of one carbon-12 atom in kilograms, `12 u`.
pub const CARBON_ATOM_MASS_KG: f64 = 1.992_646_879_92e-26;

/// Errors from rotor-machine setup and from one rotor step.
#[derive(Debug, Error, PartialEq)]
pub enum RotorError {
    #[error("rotor atom count {atom_count} must be at least three")]
    InvalidAtomCount { atom_count: usize },
    #[error("rotor radius {radius_m} m must be finite and positive")]
    InvalidRadius { radius_m: f64 },
    #[error("rotor stiffness {k_n_per_m} N/m must be finite and non-negative")]
    InvalidStiffness { k_n_per_m: f64 },
    #[error("mass {mass_kg} kg must be finite and positive")]
    InvalidMass { mass_kg: f64 },
    #[error("inertia {inertia_kg_m2} kg*m^2 must be finite and positive")]
    InvalidInertia { inertia_kg_m2: f64 },
    #[error("thermostat temperature {temperature_k} K must be finite and non-negative")]
    InvalidTemperature { temperature_k: f64 },
    #[error("thermostat rate {gamma_per_s} 1/s must be finite and non-negative")]
    InvalidGamma { gamma_per_s: f64 },
    #[error("target angular velocity {angular_velocity_rad_per_s} rad/s must be finite")]
    InvalidTargetAngularVelocity { angular_velocity_rad_per_s: f64 },
    #[error("servo gain {gain_n_m_s_per_rad} N*m*s/rad must be finite and positive")]
    InvalidServoGain { gain_n_m_s_per_rad: f64 },
    #[error("maximum motor torque {max_torque_n_m} N*m must be finite and positive")]
    InvalidMaxTorque { max_torque_n_m: f64 },
    #[error("time step {dt_s} s must be finite and positive")]
    InvalidTimeStep { dt_s: f64 },
    #[error("rotor state is not finite")]
    NonFiniteState,
    #[error(transparent)]
    Jig(#[from] JigError),
    #[error(transparent)]
    Device(#[from] DeviceError),
}

/// A Langevin thermostat over one degree of freedom.
///
/// The update is the exact Ornstein-Uhlenbeck solution. For a velocity it
/// gives `v' = c1 v + c2 xi`, with `c1 = exp(-gamma dt)` and
/// `c2 = sqrt((1 - c1^2) k_B T / m)`. The stationary variance is `k_B T / m`.
/// The rate must be non-negative and the temperature must be non-negative.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LangevinCoupling {
    gamma_per_s: f64,
    temperature_k: f64,
}

impl LangevinCoupling {
    /// Builds a Langevin coupling from a rate and a temperature.
    pub fn new(gamma_per_s: f64, temperature_k: f64) -> Result<Self, RotorError> {
        if !gamma_per_s.is_finite() || gamma_per_s < 0.0 {
            return Err(RotorError::InvalidGamma { gamma_per_s });
        }
        if !temperature_k.is_finite() || temperature_k < 0.0 {
            return Err(RotorError::InvalidTemperature { temperature_k });
        }
        Ok(Self {
            gamma_per_s,
            temperature_k,
        })
    }

    /// Returns the friction rate in inverse seconds.
    pub fn gamma_per_s(&self) -> f64 {
        self.gamma_per_s
    }

    /// Returns the thermostat temperature in kelvin.
    pub fn temperature_k(&self) -> f64 {
        self.temperature_k
    }

    /// Returns the Ornstein-Uhlenbeck decay and noise factors for `dt_s`.
    fn factors(&self, dt_s: f64) -> (f64, f64) {
        let decay = (-self.gamma_per_s * dt_s).exp();
        let noise_weight = (1.0 - decay * decay).max(0.0);
        (decay, noise_weight)
    }

    /// Applies one velocity update to a velocity component.
    ///
    /// `normal_sample` must be a standard normal sample. The mass is in
    /// kilograms. The result is the new velocity component.
    pub fn step_velocity_m_per_s(
        &self,
        velocity_m_per_s: f64,
        mass_kg: f64,
        dt_s: f64,
        normal_sample: f64,
    ) -> Result<f64, RotorError> {
        if !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(RotorError::InvalidTimeStep { dt_s });
        }
        if !mass_kg.is_finite() || mass_kg <= 0.0 {
            return Err(RotorError::InvalidMass { mass_kg });
        }
        if !normal_sample.is_finite() {
            return Err(RotorError::NonFiniteState);
        }
        let (decay, noise_weight) = self.factors(dt_s);
        let noise = (noise_weight * BOLTZMANN_J_PER_K * self.temperature_k / mass_kg).sqrt();
        Ok(decay * velocity_m_per_s + noise * normal_sample)
    }

    /// Applies one update to an angular momentum.
    ///
    /// The inertia is in kilogram metre squared. The result is the new angular
    /// momentum in kilogram metre squared per second.
    pub fn step_angular_momentum_kg_m2_per_s(
        &self,
        angular_momentum_kg_m2_per_s: f64,
        inertia_kg_m2: f64,
        dt_s: f64,
        normal_sample: f64,
    ) -> Result<f64, RotorError> {
        if !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(RotorError::InvalidTimeStep { dt_s });
        }
        if !inertia_kg_m2.is_finite() || inertia_kg_m2 <= 0.0 {
            return Err(RotorError::InvalidInertia { inertia_kg_m2 });
        }
        if !normal_sample.is_finite() {
            return Err(RotorError::NonFiniteState);
        }
        let (decay, noise_weight) = self.factors(dt_s);
        let noise = (noise_weight * BOLTZMANN_J_PER_K * self.temperature_k * inertia_kg_m2).sqrt();
        Ok(decay * angular_momentum_kg_m2_per_s + noise * normal_sample)
    }
}

/// The configuration of a [`RotorMachine`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RotorConfig {
    /// The number of ring atoms.
    pub atom_count: usize,
    /// The ring radius in metres.
    pub radius_m: f64,
    /// The stiffness of each ring bond in newtons per metre.
    pub bond_stiffness_n_per_m: f64,
    /// The stiffness of each radial spring in newtons per metre.
    pub radial_stiffness_n_per_m: f64,
    /// The motor target angular velocity in radians per second.
    pub target_angular_velocity_rad_per_s: f64,
    /// The proportional servo gain in newton metre seconds per radian.
    pub servo_gain_n_m_s_per_rad: f64,
    /// The largest motor torque in newton metres.
    pub max_torque_n_m: f64,
    /// The rotor-scale thermostat.
    pub thermostat: LangevinCoupling,
}

impl RotorConfig {
    /// Builds a default configuration for a ring of `atom_count` atoms.
    ///
    /// The target velocity is zero until the caller sets it. The thermostat is
    /// at rest.
    pub fn new(atom_count: usize, radius_m: f64) -> Result<Self, RotorError> {
        Ok(Self {
            atom_count,
            radius_m,
            bond_stiffness_n_per_m: 1.0e-2,
            radial_stiffness_n_per_m: 1.0e-2,
            target_angular_velocity_rad_per_s: 0.0,
            servo_gain_n_m_s_per_rad: 1.0e-33,
            max_torque_n_m: 1.0e-22,
            thermostat: LangevinCoupling::new(0.0, 0.0)?,
        })
    }
}

/// The measured result of one rotor step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RotorStep {
    /// The motor torque used for this step in newton metres.
    pub motor_torque_n_m: f64,
    /// The norm of the net motor force in newtons. It is a constraint residual.
    pub motor_force_residual_n: f64,
    /// The axial angular velocity in radians per second.
    pub angular_velocity_rad_per_s: f64,
    /// The total kinetic energy in joules.
    pub kinetic_energy_j: f64,
    /// The spring potential energy in joules.
    pub potential_energy_j: f64,
    /// The thermal kinetic temperature in kelvin, relative to the rigid
    /// rotation.
    pub temperature_k: f64,
    /// The largest ring-bond strain, `|d - d0| / d0`.
    pub max_bond_strain: f64,
}

/// A rotor: a bead-and-spring ring driven by a motor under a thermostat.
///
/// The machine integrates the atoms with a velocity-Verlet half-step and an
/// exact Ornstein-Uhlenbeck velocity update. The motor torque comes from the
/// proportional servo `clamp(gain * (target - omega), +/- max_torque)`. The
/// force field comes from the existing [`MotorJig`], so the net force and the
/// net axial torque are the jig's own numbers.
#[derive(Clone, Debug)]
pub struct RotorMachine {
    config: RotorConfig,
    spring: SpringJig,
    motor: MotorJig,
    positions_m: Vec<f64>,
    velocities_m_per_s: Vec<f64>,
    bond_pairs: Vec<(u32, u32)>,
    bond_rest_length_m: Vec<f64>,
    radial_rest_length_m: f64,
    mass_kg: f64,
    rng: Xorshift,
}

impl RotorMachine {
    /// Builds a rotor with the given configuration.
    pub fn new(config: RotorConfig, seed: u64) -> Result<Self, RotorError> {
        if config.atom_count < 3 {
            return Err(RotorError::InvalidAtomCount {
                atom_count: config.atom_count,
            });
        }
        if !config.radius_m.is_finite() || config.radius_m <= 0.0 {
            return Err(RotorError::InvalidRadius {
                radius_m: config.radius_m,
            });
        }
        for stiffness in [
            config.bond_stiffness_n_per_m,
            config.radial_stiffness_n_per_m,
        ] {
            if !stiffness.is_finite() || stiffness < 0.0 {
                return Err(RotorError::InvalidStiffness {
                    k_n_per_m: stiffness,
                });
            }
        }
        if !config.target_angular_velocity_rad_per_s.is_finite() {
            return Err(RotorError::InvalidTargetAngularVelocity {
                angular_velocity_rad_per_s: config.target_angular_velocity_rad_per_s,
            });
        }
        if !config.servo_gain_n_m_s_per_rad.is_finite() || config.servo_gain_n_m_s_per_rad <= 0.0 {
            return Err(RotorError::InvalidServoGain {
                gain_n_m_s_per_rad: config.servo_gain_n_m_s_per_rad,
            });
        }
        if !config.max_torque_n_m.is_finite() || config.max_torque_n_m <= 0.0 {
            return Err(RotorError::InvalidMaxTorque {
                max_torque_n_m: config.max_torque_n_m,
            });
        }

        let mass_kg = CARBON_ATOM_MASS_KG;
        let count = config.atom_count;
        let mut topology = Topology::new();
        let mut positions_m = vec![0.0; 3 * count];
        for atom in 0..count {
            let angle_rad = 2.0 * std::f64::consts::PI * atom as f64 / count as f64;
            let position_m = [
                config.radius_m * angle_rad.cos(),
                config.radius_m * angle_rad.sin(),
                0.0,
            ];
            topology.add_atom(Atom::new(Element::CARBON, position_m, 0.0, "C3"));
            let base = 3 * atom;
            positions_m[base] = position_m[0];
            positions_m[base + 1] = position_m[1];
            positions_m[base + 2] = position_m[2];
        }

        let mut spring = SpringJig::new(&topology);
        let mut bond_pairs = Vec::with_capacity(count);
        let mut bond_rest_length_m = Vec::with_capacity(count);
        for atom in 0..count {
            let next = (atom + 1) % count;
            let rest_length_m = 2.0 * config.radius_m * (std::f64::consts::PI / count as f64).sin();
            spring.add_spring(
                &topology,
                atom as u32,
                next as u32,
                config.bond_stiffness_n_per_m,
                rest_length_m,
            )?;
            bond_pairs.push((atom as u32, next as u32));
            bond_rest_length_m.push(rest_length_m);
            spring.add_spring_to_point(
                &topology,
                atom as u32,
                [0.0, 0.0, 0.0],
                config.radial_stiffness_n_per_m,
                config.radius_m,
            )?;
        }

        let atom_indices: Vec<u32> = (0..count as u32).collect();
        let motor = MotorJig::new(&topology, &atom_indices, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0])?;

        Ok(Self {
            config,
            spring,
            motor,
            positions_m,
            velocities_m_per_s: vec![0.0; 3 * count],
            bond_pairs,
            bond_rest_length_m,
            radial_rest_length_m: config.radius_m,
            mass_kg,
            rng: Xorshift::new(seed),
        })
    }

    /// Returns the configuration.
    pub fn config(&self) -> &RotorConfig {
        &self.config
    }

    /// Returns the flat position buffer `[x, y, z, ...]` in metres.
    pub fn positions_m(&self) -> &[f64] {
        &self.positions_m
    }

    /// Returns the flat velocity buffer `[vx, vy, vz, ...]` in metres per second.
    pub fn velocities_m_per_s(&self) -> &[f64] {
        &self.velocities_m_per_s
    }

    /// Returns the number of ring atoms.
    pub fn atom_count(&self) -> usize {
        self.config.atom_count
    }

    /// Returns the rotor moment of inertia about the z axis in kilogram metre
    /// squared, from the current atom positions.
    pub fn inertia_kg_m2(&self) -> f64 {
        let mut inertia_kg_m2 = 0.0;
        for atom in 0..self.config.atom_count {
            let base = 3 * atom;
            let x_m = self.positions_m[base];
            let y_m = self.positions_m[base + 1];
            inertia_kg_m2 += self.mass_kg * (x_m * x_m + y_m * y_m);
        }
        inertia_kg_m2
    }

    /// Returns the axial angular momentum `sum m (x vy - y vx)`.
    pub fn angular_momentum_kg_m2_per_s(&self) -> f64 {
        let mut angular_momentum_kg_m2_per_s = 0.0;
        for atom in 0..self.config.atom_count {
            let base = 3 * atom;
            let x_m = self.positions_m[base];
            let y_m = self.positions_m[base + 1];
            let vx_m_per_s = self.velocities_m_per_s[base];
            let vy_m_per_s = self.velocities_m_per_s[base + 1];
            angular_momentum_kg_m2_per_s += self.mass_kg * (x_m * vy_m_per_s - y_m * vx_m_per_s);
        }
        angular_momentum_kg_m2_per_s
    }

    /// Returns the axial angular velocity in radians per second.
    pub fn angular_velocity_rad_per_s(&self) -> f64 {
        let inertia_kg_m2 = self.inertia_kg_m2();
        if inertia_kg_m2 <= 0.0 {
            return 0.0;
        }
        self.angular_momentum_kg_m2_per_s() / inertia_kg_m2
    }

    /// Returns the total kinetic energy in joules.
    pub fn kinetic_energy_j(&self) -> f64 {
        let mut kinetic_energy_j = 0.0;
        for value in &self.velocities_m_per_s {
            kinetic_energy_j += 0.5 * self.mass_kg * value * value;
        }
        kinetic_energy_j
    }

    /// Returns the thermal kinetic temperature in kelvin.
    ///
    /// The method removes the center-of-mass velocity and the rigid rotation.
    /// The remaining fluctuation has `3 N - 3` degrees of freedom, because the
    /// center of mass is fixed. This is the standard fluctuation temperature.
    pub fn temperature_k(&self) -> Result<f64, RotorError> {
        let count = self.config.atom_count;
        let mut com_velocity_m_per_s = [0.0; 3];
        for atom in 0..count {
            let base = 3 * atom;
            for (axis, value) in com_velocity_m_per_s.iter_mut().enumerate() {
                *value += self.velocities_m_per_s[base + axis];
            }
        }
        for value in com_velocity_m_per_s.iter_mut() {
            *value /= count as f64;
        }
        let omega_rad_per_s = self.angular_velocity_rad_per_s();
        let mut kinetic_fluctuation_j = 0.0;
        for atom in 0..count {
            let base = 3 * atom;
            let x_m = self.positions_m[base];
            let y_m = self.positions_m[base + 1];
            let relative_m_per_s = [
                self.velocities_m_per_s[base] - com_velocity_m_per_s[0] + omega_rad_per_s * y_m,
                self.velocities_m_per_s[base + 1] - com_velocity_m_per_s[1] - omega_rad_per_s * x_m,
                self.velocities_m_per_s[base + 2] - com_velocity_m_per_s[2],
            ];
            for value in relative_m_per_s {
                kinetic_fluctuation_j += 0.5 * self.mass_kg * value * value;
            }
        }
        let degrees_of_freedom = (3 * count - 3).max(1) as f64;
        Ok(2.0 * kinetic_fluctuation_j / (degrees_of_freedom * BOLTZMANN_J_PER_K))
    }

    /// Returns the largest ring-bond strain, `|d - d0| / d0`.
    pub fn max_bond_strain(&self) -> f64 {
        let mut max_strain = 0.0_f64;
        for (index, (u, v)) in self.bond_pairs.iter().enumerate() {
            let a = self.atom_position_m(*u);
            let b = self.atom_position_m(*v);
            let dx = a[0] - b[0];
            let dy = a[1] - b[1];
            let dz = a[2] - b[2];
            let distance_m = (dx * dx + dy * dy + dz * dz).sqrt();
            let rest_length_m = self.bond_rest_length_m[index];
            let strain = (distance_m - rest_length_m).abs() / rest_length_m;
            max_strain = max_strain.max(strain);
        }
        max_strain
    }

    /// Returns the largest radial deviation `|r - r0| / r0`.
    pub fn max_radial_strain(&self) -> f64 {
        let mut max_strain = 0.0_f64;
        for atom in 0..self.config.atom_count {
            let position_m = self.atom_position_m(atom as u32);
            let radius_m = (position_m[0] * position_m[0] + position_m[1] * position_m[1]).sqrt();
            let strain = (radius_m - self.radial_rest_length_m).abs() / self.radial_rest_length_m;
            max_strain = max_strain.max(strain);
        }
        max_strain
    }

    fn atom_position_m(&self, atom: u32) -> [f64; 3] {
        let base = atom as usize * 3;
        [
            self.positions_m[base],
            self.positions_m[base + 1],
            self.positions_m[base + 2],
        ]
    }

    fn atom_velocity_m_per_s(&self, atom: u32) -> [f64; 3] {
        let base = atom as usize * 3;
        [
            self.velocities_m_per_s[base],
            self.velocities_m_per_s[base + 1],
            self.velocities_m_per_s[base + 2],
        ]
    }

    fn set_atom_velocity_m_per_s(&mut self, atom: u32, velocity_m_per_s: [f64; 3]) {
        let base = atom as usize * 3;
        self.velocities_m_per_s[base] = velocity_m_per_s[0];
        self.velocities_m_per_s[base + 1] = velocity_m_per_s[1];
        self.velocities_m_per_s[base + 2] = velocity_m_per_s[2];
    }

    /// Sets a rigid-body rotation about the z axis for every atom.
    ///
    /// This is a helper for the test setup, not an equilibrium operation.
    pub fn set_rigid_rotation(&mut self, angular_velocity_rad_per_s: f64) {
        for atom in 0..self.config.atom_count {
            let position_m = self.atom_position_m(atom as u32);
            let velocity_m_per_s = [
                -angular_velocity_rad_per_s * position_m[1],
                angular_velocity_rad_per_s * position_m[0],
                0.0,
            ];
            self.set_atom_velocity_m_per_s(atom as u32, velocity_m_per_s);
        }
    }

    /// Returns the potential energy of the springs in joules.
    pub fn potential_energy_j(&self) -> Result<f64, RotorError> {
        Ok(self.spring.energy_j(&self.positions_m)?)
    }

    /// Returns the net motor force in newtons. It must be near zero.
    pub fn motor_force_residual_n(&self) -> Result<f64, RotorError> {
        let forces = self.motor.forces_n(&self.positions_m)?;
        let mut net = [0.0; 3];
        for atom in 0..self.config.atom_count {
            let base = 3 * atom;
            net[0] += forces[base];
            net[1] += forces[base + 1];
            net[2] += forces[base + 2];
        }
        Ok((net[0] * net[0] + net[1] * net[1] + net[2] * net[2]).sqrt())
    }

    /// Advances the rotor by one frame `dt_s`.
    ///
    /// The motor torque is the proportional servo value for this frame. The
    /// method returns the measured step report. It never panics.
    pub fn step(&mut self, dt_s: f64) -> Result<RotorStep, RotorError> {
        if !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(RotorError::InvalidTimeStep { dt_s });
        }
        let torque_n_m = self.servo_torque_n_m();
        self.motor.set_torque_n_m(torque_n_m)?;

        let force_n = self.total_force_n()?;
        let half_dt = 0.5 * dt_s;
        for atom in 0..self.config.atom_count {
            let velocity_m_per_s = self.atom_velocity_m_per_s(atom as u32);
            let mut updated = [0.0; 3];
            for axis in 0..3 {
                let base = 3 * atom + axis;
                updated[axis] = velocity_m_per_s[axis] + force_n[base] / self.mass_kg * half_dt;
            }
            self.set_atom_velocity_m_per_s(atom as u32, updated);
        }

        let thermostat = self.config.thermostat;
        for atom in 0..self.config.atom_count {
            let velocity_m_per_s = self.atom_velocity_m_per_s(atom as u32);
            let mut updated = [0.0; 3];
            for axis in 0..3 {
                let sample = self.rng.normal();
                updated[axis] = thermostat.step_velocity_m_per_s(
                    velocity_m_per_s[axis],
                    self.mass_kg,
                    dt_s,
                    sample,
                )?;
            }
            self.set_atom_velocity_m_per_s(atom as u32, updated);
        }

        for atom in 0..self.config.atom_count {
            let velocity_m_per_s = self.atom_velocity_m_per_s(atom as u32);
            let base = 3 * atom;
            self.positions_m[base] += velocity_m_per_s[0] * dt_s;
            self.positions_m[base + 1] += velocity_m_per_s[1] * dt_s;
            self.positions_m[base + 2] += velocity_m_per_s[2] * dt_s;
        }

        let force_n = self.total_force_n()?;
        for atom in 0..self.config.atom_count {
            let velocity_m_per_s = self.atom_velocity_m_per_s(atom as u32);
            let mut updated = [0.0; 3];
            for axis in 0..3 {
                let base = 3 * atom + axis;
                updated[axis] = velocity_m_per_s[axis] + force_n[base] / self.mass_kg * half_dt;
            }
            self.set_atom_velocity_m_per_s(atom as u32, updated);
        }

        if self.positions_m.iter().any(|value| !value.is_finite())
            || self
                .velocities_m_per_s
                .iter()
                .any(|value| !value.is_finite())
        {
            return Err(RotorError::NonFiniteState);
        }

        Ok(RotorStep {
            motor_torque_n_m: torque_n_m,
            motor_force_residual_n: self.motor_force_residual_n()?,
            angular_velocity_rad_per_s: self.angular_velocity_rad_per_s(),
            kinetic_energy_j: self.kinetic_energy_j(),
            potential_energy_j: self.potential_energy_j()?,
            temperature_k: self.temperature_k()?,
            max_bond_strain: self.max_bond_strain(),
        })
    }

    fn servo_torque_n_m(&self) -> f64 {
        let error_rad_per_s =
            self.config.target_angular_velocity_rad_per_s - self.angular_velocity_rad_per_s();
        let raw_n_m = self.config.servo_gain_n_m_s_per_rad * error_rad_per_s;
        raw_n_m.clamp(-self.config.max_torque_n_m, self.config.max_torque_n_m)
    }

    fn total_force_n(&self) -> Result<Vec<f64>, RotorError> {
        let mut force_n = self.spring.forces_n(&self.positions_m)?;
        let motor_force_n = self.motor.forces_n(&self.positions_m)?;
        for (total, motor) in force_n.iter_mut().zip(motor_force_n.iter()) {
            *total += *motor;
        }
        Ok(force_n)
    }
}

/// A small deterministic xorshift random source with a normal sampler.
#[derive(Clone, Debug)]
struct Xorshift {
    state: u64,
    spare: Option<f64>,
}

impl Xorshift {
    fn new(seed: u64) -> Self {
        Self {
            state: seed.max(1),
            spare: None,
        }
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }

    fn next_unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    fn normal(&mut self) -> f64 {
        if let Some(spare) = self.spare.take() {
            return spare;
        }
        let u1 = self.next_unit().max(f64::MIN_POSITIVE);
        let u2 = self.next_unit();
        let radius = (-2.0 * u1.ln()).sqrt();
        let angle = 2.0 * std::f64::consts::PI * u2;
        self.spare = Some(radius * angle.sin());
        radius * angle.cos()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn machine_at_300_k() -> RotorMachine {
        let atom_count = 8;
        let radius_m = 1.0e-9;
        let inertia_kg_m2 = atom_count as f64 * CARBON_ATOM_MASS_KG * radius_m * radius_m;
        let target_angular_velocity_rad_per_s = (BOLTZMANN_J_PER_K * 300.0 / inertia_kg_m2).sqrt();
        let mut config = RotorConfig::new(atom_count, radius_m).expect("valid config");
        config.bond_stiffness_n_per_m = 400.0;
        config.radial_stiffness_n_per_m = 400.0;
        config.thermostat = LangevinCoupling::new(1.0e12, 300.0).expect("valid thermostat");
        config.target_angular_velocity_rad_per_s = target_angular_velocity_rad_per_s;
        config.servo_gain_n_m_s_per_rad = inertia_kg_m2 * config.thermostat.gamma_per_s();
        config.max_torque_n_m = 4.0
            * inertia_kg_m2
            * config.thermostat.gamma_per_s()
            * target_angular_velocity_rad_per_s;
        RotorMachine::new(config, 0x5eed).expect("valid machine")
    }

    #[test]
    fn invalid_configuration_is_rejected() {
        assert!(RotorConfig::new(2, 1.0e-9).is_ok());
        assert!(RotorMachine::new(RotorConfig::new(2, 1.0e-9).unwrap(), 1).is_err());
        let mut config = RotorConfig::new(8, 1.0e-9).unwrap();
        config.radius_m = -1.0;
        assert!(matches!(
            RotorMachine::new(config, 1),
            Err(RotorError::InvalidRadius { .. })
        ));
    }

    #[test]
    fn the_langevin_coupling_has_the_equipartition_variance() {
        let coupling = LangevinCoupling::new(5.0e12, 300.0).expect("valid coupling");
        let mass_kg = CARBON_ATOM_MASS_KG;
        let dt_s = 1.0e-13;
        let mut rng = Xorshift::new(0xabc);
        let mut sum_sq: f64 = 0.0;
        let samples = 200_000;
        let mut velocity_m_per_s = 0.0;
        for _ in 0..samples {
            velocity_m_per_s = coupling
                .step_velocity_m_per_s(velocity_m_per_s, mass_kg, dt_s, rng.normal())
                .expect("valid step");
            if sum_sq.is_finite() {
                sum_sq += velocity_m_per_s * velocity_m_per_s;
            }
        }
        let measured_m2_per_s2 = sum_sq / samples as f64;
        let expected_m2_per_s2 = BOLTZMANN_J_PER_K * 300.0 / mass_kg;
        let relative = (measured_m2_per_s2 - expected_m2_per_s2).abs() / expected_m2_per_s2;
        assert!(relative < 0.02, "relative variance error {relative}");
    }

    #[test]
    fn a_motor_drives_the_rotor_and_it_stays_stable_at_300_k() {
        let mut machine = machine_at_300_k();
        machine.set_rigid_rotation(machine.config().target_angular_velocity_rad_per_s);
        let dt_s = 5.0e-16;
        let frames = 200_000;
        let mut max_speed_m_per_s = 0.0_f64;
        let mut max_kinetic_j = 0.0_f64;
        let mut max_bond_strain = 0.0_f64;
        let mut max_force_residual_n = 0.0_f64;
        let mut temperature_sum_k = 0.0_f64;
        let mut temperature_samples = 0usize;
        for frame in 0..frames {
            let step = machine.step(dt_s).expect("stable step");
            assert!(step.kinetic_energy_j.is_finite());
            assert!(step.temperature_k.is_finite());
            max_bond_strain = max_bond_strain.max(step.max_bond_strain);
            max_force_residual_n = max_force_residual_n.max(step.motor_force_residual_n);
            if frame > frames / 2 {
                max_speed_m_per_s = max_speed_m_per_s.max(
                    machine
                        .velocities_m_per_s()
                        .iter()
                        .map(|value| value.abs())
                        .fold(0.0, f64::max),
                );
                max_kinetic_j = max_kinetic_j.max(step.kinetic_energy_j);
                temperature_sum_k += step.temperature_k;
                temperature_samples += 1;
            }
        }
        let mean_temperature_k = temperature_sum_k / temperature_samples as f64;
        assert!(
            max_bond_strain < 5.0e-2,
            "maximum bond strain {max_bond_strain}"
        );
        assert!(
            max_force_residual_n < 1.0e-24,
            "maximum motor force residual {max_force_residual_n} N"
        );
        assert!(
            (mean_temperature_k - 300.0).abs() / 300.0 < 0.15,
            "mean temperature {mean_temperature_k} K"
        );
        assert!(max_speed_m_per_s.is_finite());
        assert!(max_kinetic_j.is_finite());
        println!(
            "rotor: frames {frames} mean_T {mean_temperature_k:.2} K max_bond_strain {max_bond_strain:.3e} max_force_residual {max_force_residual_n:.3e} N max_kinetic {max_kinetic_j:.3e} J max_speed {max_speed_m_per_s:.3e} m/s"
        );
    }
}
