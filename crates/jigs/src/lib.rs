//! Jigs: anchors, motors, springs, thermostats.
//!
//! A jig is a constraint or a driver on a simulation. It reads positions from
//! `nanocad-model` and contributes an energy and a force. Positions are in SI
//! metres, energies in joules, and forces in newtons. See the [`Jig`] trait
//! for the contract.
//!
//! The anchor and spring jigs are conservative: they define a potential and a
//! gradient. The motor jig is a torque source, not a potential; its energy and
//! gradient are zero, and it supplies torque-distributed forces instead.
//!
//! The `device` module holds the L2 layer: rigid bodies and revolute,
//! prismatic, and gear joints. See [`RigidBody`] and [`RigidBodySystem`].
#![forbid(unsafe_code)]

mod anchor;
mod assembly;
mod device;
mod error;
mod jig;
mod motor;
mod rotor;
mod scene;
mod spring;
mod urdf;

pub use anchor::AnchorJig;
pub use assembly::{
    assemble_from_records, assemble_planetary, planetary_records, AssemblyError, BodyRole,
    PlanetaryAssembly, ASSEMBLY_RECORD_VERSION, CARRIER_PART_ID, PLANET_PART_ID_PREFIX,
    RING_PART_ID, SUN_PART_ID,
};
pub use device::{
    DeviceError, GearConstraint, GearCoupling, GearTerm, PrismaticJoint, Quat, RevoluteJoint,
    RigidBody, RigidBodySystem,
};
pub use error::JigError;
pub use jig::{Jig, JigKind};
pub use motor::MotorJig;
pub use rotor::{
    LangevinCoupling, RotorConfig, RotorError, RotorMachine, RotorStep, BOLTZMANN_J_PER_K,
    CARBON_ATOM_MASS_KG,
};
pub use scene::{
    build_scene, scene_from_json, scene_to_json, scene_to_json_pretty, write_scene_json,
    AtomisticLayer, CoarseBody, CoarseCircle, CoarseCylinder, CoarseLayer, DeviceLayer, Scene,
    SceneAtom, SceneBody, SceneDesign, SceneError, SceneGearConstraint, SceneGearCoupling,
    SceneGearTerm, SceneJoint, SCENE_SCHEMA, SCENE_VERSION,
};
pub use spring::SpringJig;
pub use urdf::{
    export_urdf, export_urdf_named, export_urdf_with_limits, format_urdf, parse_urdf, urdf_robot,
    UrdfError, UrdfJoint, UrdfJointKind, UrdfLimits, UrdfLink, UrdfRobot, CARRIER_LINK,
    GROUND_LINK, PLANET_LINK_PREFIX, RING_LINK, SUN_LINK, URDF_DEFAULT_EFFORT_N_M,
    URDF_DEFAULT_VELOCITY_RAD_PER_S,
};

/// Returns the crate version string.
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_is_not_empty() {
        assert!(!version().is_empty());
    }
}
