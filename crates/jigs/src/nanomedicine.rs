//! NM-01: one nanomedicine slice. A machine is coupled to a coarse host.
//!
//! The host is a coarse continuum-like environment. It supplies a viscous
//! Stokes drag, an optional flow field, an optional tether spring, and an
//! optional constant drive force. The machine is an L2 [`RigidBody`] held in a
//! [`RigidBodySystem`].
//!
//! # Host force model
//!
//! The drag coefficient of a sphere of radius `r` in a medium of dynamic
//! viscosity `eta` is the Stokes value
//!
//! ```text
//! b = 6 * pi * eta * r          [N*s/m]
//! ```
//!
//! The drag force on the body is relative to the flow,
//! `F_drag = b * (v_flow - v_body)`. The tether is a linear spring to a fixed
//! anchor, `F_tether = -k * (x_body - x_anchor)`. The drive is a constant force.
//! The total host force is the sum. The force acts through the center of mass,
//! so the host supplies no torque.
//!
//! # Thermal noise
//!
//! [`HostEnvironment::with_temperature_kbt_j`] sets the thermal energy
//! `k_B * T`. The step then adds a Brownian force. For the discrete velocity
//! update `v' = (1 - b dt / m) v + xi`, the fluctuation-dissipation relation
//! fixes the per-axis noise variance
//! `Var(xi) = (k_B T / m) (1 - (1 - b dt / m)^2)`. The stationary velocity
//! variance is then `k_B T / m`, the equipartition value. The noise needs a
//! positive drag coefficient. A [`NanoMachine`] carries a seeded random source,
//! so a thermal run is reproducible. A zero-temperature run is deterministic.
//!
//! # Analytic check
//!
//! With no tether and a constant drive `F_d`, the body relaxes to the terminal
//! velocity `v_terminal = F_d / b`. The relaxation time constant is
//! `tau = m / b`. A test compares the simulated velocity to `v_terminal`.
//!
//! # Integrator and energy accounting
//!
//! The slice uses the semi-implicit Euler step of [`RigidBodySystem`]. The work
//! done by the host force over one step uses the trapezoidal rule on the force
//! and the displacement, `W = 0.5 * (F_before + F_after) . dx`. The energy
//! balance is `final_kinetic - initial_kinetic - work`.

use thiserror::Error;

use crate::device::{DeviceError, Quat, RigidBody, RigidBodySystem};
use crate::rotor::BOLTZMANN_J_PER_K;

/// Errors from nanomedicine host setup and from a slice run.
#[derive(Debug, Error, PartialEq)]
pub enum NanoMedicineError {
    #[error("host radius {radius_m} m must be finite and positive")]
    InvalidRadius { radius_m: f64 },
    #[error("host viscosity {viscosity_pa_s} Pa*s must be finite and non-negative")]
    InvalidViscosity { viscosity_pa_s: f64 },
    #[error("mass {mass_kg} kg must be finite and positive")]
    InvalidMass { mass_kg: f64 },
    #[error("inertia {inertia_kg_m2:?} kg*m^2 must be finite and positive")]
    InvalidInertia { inertia_kg_m2: [f64; 3] },
    #[error("position {position_m:?} m must be finite")]
    InvalidPosition { position_m: [f64; 3] },
    #[error("velocity {velocity_m_per_s:?} m/s must be finite")]
    InvalidVelocity { velocity_m_per_s: [f64; 3] },
    #[error("flow velocity {flow_velocity_m_per_s:?} m/s must be finite")]
    InvalidFlowVelocity { flow_velocity_m_per_s: [f64; 3] },
    #[error("drive force {drive_force_n:?} N must be finite")]
    InvalidDriveForce { drive_force_n: [f64; 3] },
    #[error("tether anchor {anchor_m:?} m must be finite")]
    InvalidTetherAnchor { anchor_m: [f64; 3] },
    #[error("tether stiffness {stiffness_n_per_m} N/m must be finite and non-negative")]
    InvalidTetherStiffness { stiffness_n_per_m: f64 },
    #[error("time step {dt_s} s must be finite and positive")]
    InvalidTimeStep { dt_s: f64 },
    #[error("thermal energy k_B*T {kbt_j} J must be finite and non-negative")]
    InvalidTemperature { kbt_j: f64 },
    #[error("body index {body} is outside the body count {body_count}")]
    BodyIndexOutOfBounds { body: usize, body_count: usize },
    #[error("the host has no drag, so a terminal velocity does not exist")]
    NoDrag,
    #[error("the slice state is not finite")]
    NonFiniteState,
    #[error(transparent)]
    Device(#[from] DeviceError),
}

/// A coarse host environment around a nanoscale machine.
///
/// The environment is a viscous medium. It can also carry a uniform flow, a
/// tether spring, and a constant drive force. The struct is cheap to copy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HostEnvironment {
    viscosity_pa_s: f64,
    radius_m: f64,
    flow_velocity_m_per_s: [f64; 3],
    drive_force_n: [f64; 3],
    tether_anchor_m: Option<[f64; 3]>,
    tether_stiffness_n_per_m: f64,
    kbt_j: f64,
}

impl HostEnvironment {
    /// Builds a Stokes host with the given dynamic viscosity and body radius.
    ///
    /// The viscosity must be finite and non-negative. The radius must be finite
    /// and positive. The flow, the drive, and the tether start at zero.
    pub fn stokes(viscosity_pa_s: f64, radius_m: f64) -> Result<Self, NanoMedicineError> {
        if !viscosity_pa_s.is_finite() || viscosity_pa_s < 0.0 {
            return Err(NanoMedicineError::InvalidViscosity { viscosity_pa_s });
        }
        if !radius_m.is_finite() || radius_m <= 0.0 {
            return Err(NanoMedicineError::InvalidRadius { radius_m });
        }
        Ok(Self {
            viscosity_pa_s,
            radius_m,
            flow_velocity_m_per_s: [0.0; 3],
            drive_force_n: [0.0; 3],
            tether_anchor_m: None,
            tether_stiffness_n_per_m: 0.0,
            kbt_j: 0.0,
        })
    }

    /// Sets the thermal energy `k_B * T` of the medium, in joules.
    ///
    /// The value must be finite and non-negative. A value of zero, the default,
    /// makes the host deterministic. A positive value turns on the Brownian
    /// force. The force amplitude obeys the fluctuation-dissipation relation
    /// for the Stokes drag, so the equilibrium velocity variance is `k_B T / m`
    /// on each axis. Thermal noise needs a positive drag coefficient: without
    /// drag the medium cannot exchange energy with the body.
    pub fn with_temperature_kbt_j(mut self, kbt_j: f64) -> Result<Self, NanoMedicineError> {
        if !kbt_j.is_finite() || kbt_j < 0.0 {
            return Err(NanoMedicineError::InvalidTemperature { kbt_j });
        }
        self.kbt_j = kbt_j;
        Ok(self)
    }

    /// Returns the thermal energy `k_B * T` in joules.
    pub fn kbt_j(&self) -> f64 {
        self.kbt_j
    }

    /// Returns the medium temperature in kelvin, `k_B T / k_B`.
    pub fn temperature_k(&self) -> f64 {
        self.kbt_j / BOLTZMANN_J_PER_K
    }

    /// Sets the uniform flow velocity of the medium in metres per second.
    pub fn with_flow_velocity(
        mut self,
        flow_velocity_m_per_s: [f64; 3],
    ) -> Result<Self, NanoMedicineError> {
        if flow_velocity_m_per_s.iter().any(|value| !value.is_finite()) {
            return Err(NanoMedicineError::InvalidFlowVelocity {
                flow_velocity_m_per_s,
            });
        }
        self.flow_velocity_m_per_s = flow_velocity_m_per_s;
        Ok(self)
    }

    /// Sets a constant drive force on the machine in newtons.
    pub fn with_drive_force(mut self, drive_force_n: [f64; 3]) -> Result<Self, NanoMedicineError> {
        if drive_force_n.iter().any(|value| !value.is_finite()) {
            return Err(NanoMedicineError::InvalidDriveForce { drive_force_n });
        }
        self.drive_force_n = drive_force_n;
        Ok(self)
    }

    /// Adds a linear tether spring to a fixed anchor.
    ///
    /// The stiffness must be finite and non-negative. A stiffness of zero is
    /// the same as no tether.
    pub fn with_tether(
        mut self,
        anchor_m: [f64; 3],
        stiffness_n_per_m: f64,
    ) -> Result<Self, NanoMedicineError> {
        if anchor_m.iter().any(|value| !value.is_finite()) {
            return Err(NanoMedicineError::InvalidTetherAnchor { anchor_m });
        }
        if !stiffness_n_per_m.is_finite() || stiffness_n_per_m < 0.0 {
            return Err(NanoMedicineError::InvalidTetherStiffness { stiffness_n_per_m });
        }
        self.tether_anchor_m = Some(anchor_m);
        self.tether_stiffness_n_per_m = stiffness_n_per_m;
        Ok(self)
    }

    /// Returns the dynamic viscosity in pascal seconds.
    pub fn viscosity_pa_s(&self) -> f64 {
        self.viscosity_pa_s
    }

    /// Returns the body radius in metres.
    pub fn radius_m(&self) -> f64 {
        self.radius_m
    }

    /// Returns the uniform flow velocity in metres per second.
    pub fn flow_velocity_m_per_s(&self) -> [f64; 3] {
        self.flow_velocity_m_per_s
    }

    /// Returns the constant drive force in newtons.
    pub fn drive_force_n(&self) -> [f64; 3] {
        self.drive_force_n
    }

    /// Returns the tether anchor in metres, or `None` when there is no tether.
    pub fn tether_anchor_m(&self) -> Option<[f64; 3]> {
        self.tether_anchor_m
    }

    /// Returns the tether stiffness in newtons per metre.
    pub fn tether_stiffness_n_per_m(&self) -> f64 {
        self.tether_stiffness_n_per_m
    }

    /// Returns the Stokes drag coefficient `6 * pi * eta * r` in N*s/m.
    pub fn drag_coefficient_n_s_per_m(&self) -> f64 {
        6.0 * std::f64::consts::PI * self.viscosity_pa_s * self.radius_m
    }

    /// Returns the drag force relative to the flow, in newtons.
    pub fn drag_force_n(&self, velocity_m_per_s: [f64; 3]) -> [f64; 3] {
        let coefficient = self.drag_coefficient_n_s_per_m();
        scale3(
            sub3(self.flow_velocity_m_per_s, velocity_m_per_s),
            coefficient,
        )
    }

    /// Returns the tether force in newtons, or zero when there is no tether.
    pub fn tether_force_n(&self, position_m: [f64; 3]) -> [f64; 3] {
        match self.tether_anchor_m {
            Some(anchor_m) => scale3(sub3(anchor_m, position_m), self.tether_stiffness_n_per_m),
            None => [0.0; 3],
        }
    }

    /// Returns the total host force on the machine in newtons.
    ///
    /// The force is the sum of the drag, the tether, and the drive.
    pub fn force_on(&self, velocity_m_per_s: [f64; 3], position_m: [f64; 3]) -> [f64; 3] {
        add3(
            add3(
                self.drag_force_n(velocity_m_per_s),
                self.tether_force_n(position_m),
            ),
            self.drive_force_n,
        )
    }

    /// Returns the terminal velocity `F_drive / b`, in metres per second.
    ///
    /// The terminal velocity exists only when the drag coefficient is positive.
    /// The method ignores the tether: a tether has no terminal velocity.
    pub fn terminal_velocity_m_per_s(&self) -> Result<[f64; 3], NanoMedicineError> {
        let coefficient = self.drag_coefficient_n_s_per_m();
        if coefficient <= 0.0 {
            return Err(NanoMedicineError::NoDrag);
        }
        Ok(scale3(self.drive_force_n, 1.0 / coefficient))
    }

    /// Returns the drag relaxation time `m / b` in seconds.
    pub fn relaxation_time_s(&self, mass_kg: f64) -> Result<f64, NanoMedicineError> {
        let coefficient = self.drag_coefficient_n_s_per_m();
        if coefficient <= 0.0 {
            return Err(NanoMedicineError::NoDrag);
        }
        if !mass_kg.is_finite() || mass_kg <= 0.0 {
            return Err(NanoMedicineError::InvalidMass { mass_kg });
        }
        Ok(mass_kg / coefficient)
    }
}

/// The measured result of one coupled slice.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SliceResult {
    /// The final world position of the center of mass in metres.
    pub final_position_m: [f64; 3],
    /// The final world linear velocity in metres per second.
    pub final_velocity_m_per_s: [f64; 3],
    /// The final drag force in newtons.
    pub drag_force_n: [f64; 3],
    /// The final tether force in newtons.
    pub tether_force_n: [f64; 3],
    /// The constant drive force in newtons.
    pub drive_force_n: [f64; 3],
    /// The work done by the total host force in joules.
    pub work_j: f64,
    /// The kinetic energy at the start of the slice in joules.
    pub initial_kinetic_energy_j: f64,
    /// The kinetic energy at the end of the slice in joules.
    pub kinetic_energy_j: f64,
    /// The energy balance `kinetic_energy_j - initial_kinetic_energy_j - work_j`
    /// in joules. It is zero when the work accounts for every energy change.
    pub energy_balance_j: f64,
    /// The number of steps run.
    pub steps: u64,
    /// The time step in seconds.
    pub dt_s: f64,
}

/// A nanomedicine machine: one rigid body and its coarse host.
///
/// The machine owns a [`RigidBodySystem`] with a single body and the
/// [`HostEnvironment`] that acts on it. Use [`NanoMachine::run`] to run the
/// coupled slice.
#[derive(Clone, Debug)]
pub struct NanoMachine {
    system: RigidBodySystem,
    body: usize,
    environment: HostEnvironment,
    rng: Xorshift,
}

/// The fixed seed of a machine that a caller does not seed explicitly.
const DEFAULT_RANDOM_SEED: u64 = 0x9E37_79B9_7F4A_7C15;

impl NanoMachine {
    /// Builds a machine from a host, a mass, a diagonal inertia, and a pose.
    ///
    /// The mass must be finite and positive. The inertia diagonal must be finite
    /// and positive. The position must be finite. The orientation must be finite
    /// and non-zero. Every failure returns a typed error, never a panic.
    pub fn new(
        environment: HostEnvironment,
        mass_kg: f64,
        inertia_diagonal_kg_m2: [f64; 3],
        position_m: [f64; 3],
        orientation: Quat,
    ) -> Result<Self, NanoMedicineError> {
        if !mass_kg.is_finite() || mass_kg <= 0.0 {
            return Err(NanoMedicineError::InvalidMass { mass_kg });
        }
        if inertia_diagonal_kg_m2
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(NanoMedicineError::InvalidInertia {
                inertia_kg_m2: inertia_diagonal_kg_m2,
            });
        }
        if position_m.iter().any(|value| !value.is_finite()) {
            return Err(NanoMedicineError::InvalidPosition { position_m });
        }
        let body = RigidBody::with_diagonal_inertia(
            mass_kg,
            inertia_diagonal_kg_m2,
            position_m,
            orientation,
        )?;
        let mut system = RigidBodySystem::new();
        let body = system.add_body(body);
        Ok(Self {
            system,
            body,
            environment,
            rng: Xorshift::new(DEFAULT_RANDOM_SEED),
        })
    }

    /// Returns a copy of the machine with a new Brownian random seed.
    ///
    /// Two machines with the same seed and the same host produce the same
    /// thermal trajectory. A run with zero temperature does not read the seed.
    pub fn with_seed(mut self, seed: u64) -> Self {
        self.rng = Xorshift::new(seed);
        self
    }

    /// Returns the device system that holds the body.
    pub fn system(&self) -> &RigidBodySystem {
        &self.system
    }

    /// Returns the device system for mutation.
    pub fn system_mut(&mut self) -> &mut RigidBodySystem {
        &mut self.system
    }

    /// Returns the body index of the machine in its system.
    pub fn body_index(&self) -> usize {
        self.body
    }

    /// Returns the machine body.
    ///
    /// The body is always present, so this method uses the system accessor and
    /// returns a reference. A missing body is a construction bug, not user
    /// input.
    pub fn body(&self) -> &RigidBody {
        &self.system.bodies()[self.body]
    }

    /// Returns the machine body for mutation.
    pub fn body_mut(&mut self) -> &mut RigidBody {
        &mut self.system.bodies_mut()[self.body]
    }

    /// Returns the coarse host.
    pub fn environment(&self) -> &HostEnvironment {
        &self.environment
    }

    /// Replaces the coarse host.
    pub fn set_environment(&mut self, environment: HostEnvironment) {
        self.environment = environment;
    }

    /// Sets the body linear velocity in metres per second.
    pub fn set_linear_velocity_m_per_s(
        &mut self,
        velocity_m_per_s: [f64; 3],
    ) -> Result<(), NanoMedicineError> {
        if velocity_m_per_s.iter().any(|value| !value.is_finite()) {
            return Err(NanoMedicineError::InvalidVelocity { velocity_m_per_s });
        }
        self.body_mut()
            .set_linear_velocity_m_per_s(velocity_m_per_s)?;
        Ok(())
    }

    /// Sets the body world position of the center of mass in metres.
    pub fn set_position_m(&mut self, position_m: [f64; 3]) -> Result<(), NanoMedicineError> {
        if position_m.iter().any(|value| !value.is_finite()) {
            return Err(NanoMedicineError::InvalidPosition { position_m });
        }
        self.body_mut().set_position_m(position_m)?;
        Ok(())
    }

    /// Runs the coupled slice for `steps` steps of `dt_s` and returns the result.
    ///
    /// When the host has a positive `k_B T`, the step adds a Brownian force.
    /// The force uses the machine random source, so two machines with the same
    /// seed produce the same trajectory. A zero-temperature run is
    /// deterministic and ignores the seed.
    pub fn run(&mut self, dt_s: f64, steps: u64) -> Result<SliceResult, NanoMedicineError> {
        run_slice_with_rng(
            &mut self.system,
            self.body,
            &self.environment,
            dt_s,
            steps,
            &mut self.rng,
        )
    }
}

/// Steps a coupled body and host for `steps` steps and returns the result.
///
/// The function applies the host force to the body at every step, integrates
/// with the semi-implicit Euler step of [`RigidBodySystem`], and accumulates the
/// work of the host force. It returns a typed error and never panics.
pub fn run_slice(
    system: &mut RigidBodySystem,
    body: usize,
    environment: &HostEnvironment,
    dt_s: f64,
    steps: u64,
) -> Result<SliceResult, NanoMedicineError> {
    let mut rng = Xorshift::new(DEFAULT_RANDOM_SEED);
    run_slice_with_rng(system, body, environment, dt_s, steps, &mut rng)
}

/// Steps a coupled body and host with an explicit Brownian random source.
///
/// The deterministic [`run_slice`] uses a fixed seed. This function lets a
/// [`NanoMachine`] carry the seed across calls.
fn run_slice_with_rng(
    system: &mut RigidBodySystem,
    body: usize,
    environment: &HostEnvironment,
    dt_s: f64,
    steps: u64,
    rng: &mut Xorshift,
) -> Result<SliceResult, NanoMedicineError> {
    if !dt_s.is_finite() || dt_s <= 0.0 {
        return Err(NanoMedicineError::InvalidTimeStep { dt_s });
    }
    let initial = read_state(system, body)?;
    let mut position_m = initial.position_m;
    let mut velocity_m_per_s = initial.velocity_m_per_s;
    let mut work_j = 0.0;
    let mass_kg = system
        .body(body)
        .ok_or(NanoMedicineError::BodyIndexOutOfBounds {
            body,
            body_count: system.body_count(),
        })?
        .mass_kg();
    let thermal = environment.kbt_j() > 0.0 && environment.drag_coefficient_n_s_per_m() > 0.0;

    for _ in 0..steps {
        let host_force_n = environment.force_on(velocity_m_per_s, position_m);
        let mut applied_force_n = host_force_n;
        if thermal {
            applied_force_n = add3(
                applied_force_n,
                thermal_force_n(environment, mass_kg, dt_s, rng),
            );
        }
        {
            let body_count = system.body_count();
            let rigid_body = system
                .body_mut(body)
                .ok_or(NanoMedicineError::BodyIndexOutOfBounds { body, body_count })?;
            rigid_body.clear_forces();
            rigid_body.add_force_n(applied_force_n);
        }
        system.step(dt_s, 0)?;
        let after = read_state(system, body)?;
        let force_after_n = environment.force_on(after.velocity_m_per_s, after.position_m);
        let displacement_m = sub3(after.position_m, position_m);
        work_j += 0.5 * dot(add3(host_force_n, force_after_n), displacement_m);
        position_m = after.position_m;
        velocity_m_per_s = after.velocity_m_per_s;
    }

    let final_state = read_state(system, body)?;
    if !final_state.position_m.iter().all(|value| value.is_finite())
        || !final_state
            .velocity_m_per_s
            .iter()
            .all(|value| value.is_finite())
        || !work_j.is_finite()
    {
        return Err(NanoMedicineError::NonFiniteState);
    }

    let kinetic_energy_j = final_state.kinetic_energy_j;
    Ok(SliceResult {
        final_position_m: final_state.position_m,
        final_velocity_m_per_s: final_state.velocity_m_per_s,
        drag_force_n: environment.drag_force_n(final_state.velocity_m_per_s),
        tether_force_n: environment.tether_force_n(final_state.position_m),
        drive_force_n: environment.drive_force_n(),
        work_j,
        initial_kinetic_energy_j: initial.kinetic_energy_j,
        kinetic_energy_j,
        energy_balance_j: kinetic_energy_j - initial.kinetic_energy_j - work_j,
        steps,
        dt_s,
    })
}

struct BodyState {
    position_m: [f64; 3],
    velocity_m_per_s: [f64; 3],
    kinetic_energy_j: f64,
}

fn read_state(system: &RigidBodySystem, body: usize) -> Result<BodyState, NanoMedicineError> {
    let rigid_body = system
        .body(body)
        .ok_or(NanoMedicineError::BodyIndexOutOfBounds {
            body,
            body_count: system.body_count(),
        })?;
    Ok(BodyState {
        position_m: rigid_body.position_m(),
        velocity_m_per_s: rigid_body.linear_velocity_m_per_s(),
        kinetic_energy_j: rigid_body.kinetic_energy_j(),
    })
}

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

/// Returns the Brownian force for one step, in newtons.
///
/// For zero flow and zero applied force the discrete velocity update is
/// `v' = (1 - b dt / m) v + xi`, where `b` is the Stokes drag coefficient. The
/// fluctuation-dissipation relation fixes the per-axis noise variance
/// `Var(xi) = (k_B T / m) (1 - (1 - b dt / m)^2)`, so the stationary variance
/// is `k_B T / m`. The force that gives `xi` over `dt_s` is
/// `m * sqrt(Var(xi)) / dt_s` times a standard normal vector.
fn thermal_force_n(
    environment: &HostEnvironment,
    mass_kg: f64,
    dt_s: f64,
    rng: &mut Xorshift,
) -> [f64; 3] {
    let coefficient = environment.drag_coefficient_n_s_per_m();
    let decay = 1.0 - coefficient * dt_s / mass_kg;
    let variance_m2_per_s2 = (environment.kbt_j() / mass_kg) * (1.0 - decay * decay).max(0.0);
    let scale_n = mass_kg * variance_m2_per_s2.sqrt() / dt_s;
    [
        rng.normal() * scale_n,
        rng.normal() * scale_n,
        rng.normal() * scale_n,
    ]
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

    const WATER_VISCOSITY_PA_S: f64 = 1.0e-3;
    const PROBE_RADIUS_M: f64 = 50.0e-9;
    const PROBE_MASS_KG: f64 = 1.0e-18;
    const DRIVE_FORCE_N: f64 = 1.0e-15;

    fn probe_environment() -> HostEnvironment {
        HostEnvironment::stokes(WATER_VISCOSITY_PA_S, PROBE_RADIUS_M)
            .expect("valid host")
            .with_drive_force([DRIVE_FORCE_N, 0.0, 0.0])
            .expect("valid drive")
    }

    fn probe_machine() -> NanoMachine {
        NanoMachine::new(
            probe_environment(),
            PROBE_MASS_KG,
            [1.0e-36; 3],
            [0.0; 3],
            Quat::IDENTITY,
        )
        .expect("valid machine")
    }

    fn assert_close(actual: f64, expected: f64, tolerance: f64, label: &str) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "{label}: actual {actual} expected {expected} tolerance {tolerance}"
        );
    }

    #[test]
    fn stokes_drag_reaches_the_analytic_terminal_velocity() {
        let mut machine = probe_machine();
        let environment = *machine.environment();
        let coefficient = environment.drag_coefficient_n_s_per_m();
        let expected_m_per_s = DRIVE_FORCE_N / coefficient;
        let analytic = environment
            .terminal_velocity_m_per_s()
            .expect("drag is present");
        assert_close(analytic[0], expected_m_per_s, 1.0e-30, "terminal accessor");

        let tau_s = environment
            .relaxation_time_s(PROBE_MASS_KG)
            .expect("drag is present");
        let dt_s = tau_s / 1000.0;
        let result = machine.run(dt_s, 20_000).expect("valid run");

        let relative =
            (result.final_velocity_m_per_s[0] - expected_m_per_s).abs() / expected_m_per_s;
        assert!(
            relative < 1.0e-6,
            "terminal velocity relative error {relative}, simulated {} m/s analytic {expected_m_per_s} m/s",
            result.final_velocity_m_per_s[0]
        );
        assert_close(
            result.drag_force_n[0],
            -DRIVE_FORCE_N,
            1.0e-7 * DRIVE_FORCE_N,
            "drag balances the drive",
        );
        println!(
            "stokes: simulated {:.6e} m/s analytic {:.6e} m/s relative_error {:.3e} tau {:.3e} s",
            result.final_velocity_m_per_s[0], expected_m_per_s, relative, tau_s
        );
    }

    #[test]
    fn the_host_force_work_closes_the_energy_balance() {
        let mut machine = probe_machine();
        let tau_s = machine
            .environment()
            .relaxation_time_s(PROBE_MASS_KG)
            .expect("drag is present");
        let dt_s = tau_s / 10_000.0;
        let result = machine.run(dt_s, 50_000).expect("valid run");

        let scale = result.kinetic_energy_j.abs()
            + result.work_j.abs()
            + result.initial_kinetic_energy_j.abs()
            + f64::MIN_POSITIVE;
        let relative = result.energy_balance_j.abs() / scale;
        assert!(
            relative < 1.0e-3,
            "energy balance relative residual {relative}, balance {} J work {} J kinetic {} J",
            result.energy_balance_j,
            result.work_j,
            result.kinetic_energy_j
        );
        println!(
            "energy: balance {:.3e} J work {:.3e} J kinetic {:.3e} J relative_residual {:.3e}",
            result.energy_balance_j, result.work_j, result.kinetic_energy_j, relative
        );
    }

    #[test]
    fn a_zero_force_body_conserves_momentum() {
        let environment =
            HostEnvironment::stokes(0.0, PROBE_RADIUS_M).expect("valid zero-viscosity host");
        let mut machine = NanoMachine::new(
            environment,
            PROBE_MASS_KG,
            [1.0e-36; 3],
            [0.0; 3],
            Quat::IDENTITY,
        )
        .expect("valid machine");
        machine
            .set_linear_velocity_m_per_s([1.0, -2.0, 0.5])
            .expect("valid velocity");
        let momentum_before = machine.body().momentum_kg_m_per_s();
        let result = machine.run(1.0e-12, 1_000).expect("valid run");
        for axis in 0..3 {
            assert_close(
                result.final_velocity_m_per_s[axis],
                [1.0, -2.0, 0.5][axis],
                1.0e-12,
                "conserve velocity",
            );
            let momentum_after = PROBE_MASS_KG * result.final_velocity_m_per_s[axis];
            assert_close(
                momentum_after,
                momentum_before[axis],
                1.0e-30,
                "conserve momentum",
            );
        }
        assert_close(result.work_j, 0.0, 1.0e-30, "no work without force");
    }

    #[test]
    fn the_flow_field_sets_the_relaxation_target() {
        let flow_m_per_s = [2.0e-6, 0.0, 0.0];
        let environment = HostEnvironment::stokes(WATER_VISCOSITY_PA_S, PROBE_RADIUS_M)
            .expect("valid host")
            .with_flow_velocity(flow_m_per_s)
            .expect("valid flow");
        let mut machine = NanoMachine::new(
            environment,
            PROBE_MASS_KG,
            [1.0e-36; 3],
            [0.0; 3],
            Quat::IDENTITY,
        )
        .expect("valid machine");
        let tau_s = environment
            .relaxation_time_s(PROBE_MASS_KG)
            .expect("drag is present");
        let result = machine.run(tau_s / 1000.0, 20_000).expect("valid run");
        let relative = (result.final_velocity_m_per_s[0] - flow_m_per_s[0]).abs() / flow_m_per_s[0];
        assert!(
            relative < 1.0e-6,
            "flow relaxation relative error {relative}"
        );
    }

    #[test]
    fn the_tether_force_matches_the_finite_difference_gradient() {
        let stiffness_n_per_m = 2.0e-6;
        let environment = HostEnvironment::stokes(0.0, PROBE_RADIUS_M)
            .expect("valid host")
            .with_tether([0.0; 3], stiffness_n_per_m)
            .expect("valid tether");
        let potential_j = |position_m: [f64; 3]| -> f64 {
            0.5 * stiffness_n_per_m
                * (position_m[0] * position_m[0]
                    + position_m[1] * position_m[1]
                    + position_m[2] * position_m[2])
        };
        let position_m = [1.0e-8, -3.0e-8, 2.0e-8];
        let force_n = environment.tether_force_n(position_m);
        let step_m = 1.0e-14;
        for axis in 0..3 {
            let mut plus = position_m;
            let mut minus = position_m;
            plus[axis] += step_m;
            minus[axis] -= step_m;
            let numeric = -(potential_j(plus) - potential_j(minus)) / (2.0 * step_m);
            assert_close(force_n[axis], numeric, 1.0e-18, "tether gradient axis");
        }
    }

    #[test]
    fn thermal_equilibrium_velocity_variance_matches_equipartition() {
        let temperature_k = 300.0;
        let kbt_j = BOLTZMANN_J_PER_K * temperature_k;
        let environment = HostEnvironment::stokes(WATER_VISCOSITY_PA_S, PROBE_RADIUS_M)
            .expect("valid host")
            .with_temperature_kbt_j(kbt_j)
            .expect("valid temperature");
        let coefficient = environment.drag_coefficient_n_s_per_m();
        let tau_s = PROBE_MASS_KG / coefficient;
        let dt_s = tau_s / 100.0;
        let mut machine = NanoMachine::new(
            environment,
            PROBE_MASS_KG,
            [1.0e-36; 3],
            [0.0; 3],
            Quat::IDENTITY,
        )
        .expect("valid machine")
        .with_seed(0x5EED_5EED);

        machine.run(dt_s, 2_000).expect("warmup");
        let samples = 200_000u64;
        let mut sum = 0.0;
        let mut sum_sq = 0.0;
        for _ in 0..samples {
            machine.run(dt_s, 1).expect("step");
            let velocity_x_m_per_s = machine.body().linear_velocity_m_per_s()[0];
            sum += velocity_x_m_per_s;
            sum_sq += velocity_x_m_per_s * velocity_x_m_per_s;
        }
        let count = samples as f64;
        let mean_m_per_s = sum / count;
        let variance_m2_per_s2 = sum_sq / count - mean_m_per_s * mean_m_per_s;
        let expected_m2_per_s2 = kbt_j / PROBE_MASS_KG;
        let relative = (variance_m2_per_s2 - expected_m2_per_s2).abs() / expected_m2_per_s2;
        assert!(
            relative < 0.12,
            "velocity variance {variance_m2_per_s2:e} m^2/s^2 differs from k_B T / m \
             {expected_m2_per_s2:e} m^2/s^2 by relative {relative}"
        );
        assert!(
            (machine.environment().temperature_k() - temperature_k).abs() / temperature_k < 1.0e-12
        );
        println!(
            "brownian: variance {variance_m2_per_s2:.6e} m^2/s^2 k_B T / m {expected_m2_per_s2:.6e} \
             relative_error {relative:.3e}"
        );
    }

    #[test]
    fn zero_temperature_is_deterministic_and_matches_terminal_velocity() {
        let run_once = || {
            let mut machine = probe_machine().with_seed(7);
            let tau_s = machine
                .environment()
                .relaxation_time_s(PROBE_MASS_KG)
                .expect("drag is present");
            machine.run(tau_s / 1000.0, 20_000).expect("valid run")
        };
        let first = run_once();
        let second = run_once();
        assert_eq!(first.final_velocity_m_per_s, second.final_velocity_m_per_s);
        assert_eq!(first.final_position_m, second.final_position_m);

        let expected_m_per_s = DRIVE_FORCE_N / probe_environment().drag_coefficient_n_s_per_m();
        let relative =
            (first.final_velocity_m_per_s[0] - expected_m_per_s).abs() / expected_m_per_s;
        assert!(
            relative < 1.0e-6,
            "terminal velocity relative error {relative}, simulated {} m/s analytic \
             {expected_m_per_s} m/s",
            first.final_velocity_m_per_s[0]
        );
    }

    #[test]
    fn a_seeded_thermal_run_is_reproducible() {
        let kbt_j = BOLTZMANN_J_PER_K * 300.0;
        let build = |seed: u64| {
            let environment = HostEnvironment::stokes(WATER_VISCOSITY_PA_S, PROBE_RADIUS_M)
                .expect("valid host")
                .with_temperature_kbt_j(kbt_j)
                .expect("valid temperature");
            let mut machine = NanoMachine::new(
                environment,
                PROBE_MASS_KG,
                [1.0e-36; 3],
                [0.0; 3],
                Quat::IDENTITY,
            )
            .expect("valid machine")
            .with_seed(seed);
            machine.run(1.0e-13, 100).expect("valid run")
        };
        let first = build(11);
        let second = build(11);
        let third = build(12);
        assert_eq!(first.final_velocity_m_per_s, second.final_velocity_m_per_s);
        assert_ne!(first.final_velocity_m_per_s, third.final_velocity_m_per_s);
    }

    #[test]
    fn invalid_inputs_return_a_typed_error() {
        assert!(matches!(
            HostEnvironment::stokes(-1.0e-3, PROBE_RADIUS_M),
            Err(NanoMedicineError::InvalidViscosity { .. })
        ));
        assert!(matches!(
            HostEnvironment::stokes(f64::NAN, PROBE_RADIUS_M),
            Err(NanoMedicineError::InvalidViscosity { .. })
        ));
        assert!(matches!(
            HostEnvironment::stokes(WATER_VISCOSITY_PA_S, -PROBE_RADIUS_M),
            Err(NanoMedicineError::InvalidRadius { .. })
        ));
        assert!(matches!(
            HostEnvironment::stokes(WATER_VISCOSITY_PA_S, 0.0),
            Err(NanoMedicineError::InvalidRadius { .. })
        ));
        assert!(matches!(
            probe_environment().with_drive_force([f64::INFINITY, 0.0, 0.0]),
            Err(NanoMedicineError::InvalidDriveForce { .. })
        ));
        assert!(matches!(
            probe_environment().with_tether([0.0; 3], -1.0),
            Err(NanoMedicineError::InvalidTetherStiffness { .. })
        ));
        assert!(matches!(
            probe_environment().with_temperature_kbt_j(-1.0e-21),
            Err(NanoMedicineError::InvalidTemperature { .. })
        ));
        assert!(matches!(
            probe_environment().with_temperature_kbt_j(f64::NAN),
            Err(NanoMedicineError::InvalidTemperature { .. })
        ));
        assert!(matches!(
            NanoMachine::new(
                probe_environment(),
                0.0,
                [1.0e-36; 3],
                [0.0; 3],
                Quat::IDENTITY
            ),
            Err(NanoMedicineError::InvalidMass { .. })
        ));
        assert!(matches!(
            NanoMachine::new(
                probe_environment(),
                PROBE_MASS_KG,
                [f64::NAN, 1.0e-36, 1.0e-36],
                [0.0; 3],
                Quat::IDENTITY
            ),
            Err(NanoMedicineError::InvalidInertia { .. })
        ));
        assert!(matches!(
            NanoMachine::new(
                probe_environment(),
                PROBE_MASS_KG,
                [1.0e-36; 3],
                [f64::NAN, 0.0, 0.0],
                Quat::IDENTITY
            ),
            Err(NanoMedicineError::InvalidPosition { .. })
        ));
        let mut machine = probe_machine();
        assert!(matches!(
            machine.run(0.0, 10),
            Err(NanoMedicineError::InvalidTimeStep { .. })
        ));
        assert!(matches!(
            machine.run(f64::NAN, 10),
            Err(NanoMedicineError::InvalidTimeStep { .. })
        ));
        assert!(matches!(
            HostEnvironment::stokes(0.0, PROBE_RADIUS_M)
                .expect("valid")
                .terminal_velocity_m_per_s(),
            Err(NanoMedicineError::NoDrag)
        ));
    }
}
