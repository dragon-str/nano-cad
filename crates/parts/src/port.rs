//! Ports and port mating.
//!
//! A [`Port`] is a connection point on a part. A [`PortFrame`] is the planar
//! transform that brings a plug port onto a socket port. The crate holds SI
//! positions and unit-suffixed names, following ADR-0003.

use std::f64::consts::PI;

/// The motion a port allows.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Dof {
    /// The joint holds both parts rigidly together.
    Fixed,
    /// The joint allows one rotation about the port axis.
    Revolute,
    /// The joint allows one translation along the port axis.
    Prismatic,
    /// The joint couples one rotation to one translation along the axis.
    Screw,
}

/// A connection point on a part. The axis points out of the part.
#[derive(Clone, Debug, PartialEq)]
pub struct Port {
    /// The port name. It is unique inside one part.
    pub name: String,
    /// The port position in the part frame, in metres.
    pub origin_m: [f64; 3],
    /// The outward port direction. It must be a unit vector.
    pub axis_m: [f64; 3],
    /// The motion the port allows.
    pub dof: Dof,
    /// The standoff between the port origin and the mating face, in metres.
    pub gap_m: f64,
}

impl Port {
    /// A port with the given name, at the origin, along +z, fixed, with no gap.
    pub fn named(name: &str) -> Self {
        Self {
            name: name.to_string(),
            origin_m: [0.0, 0.0, 0.0],
            axis_m: [0.0, 0.0, 1.0],
            dof: Dof::Fixed,
            gap_m: 0.0,
        }
    }
}

/// The transform that brings a plug onto a socket.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PortFrame {
    /// The rotation about the part z axis, in radians.
    pub rotation_rad: f64,
    /// The translation to add to the plug origin, in metres.
    pub offset_m: [f64; 3],
}

/// Returns the transform that mates a plug to a socket.
///
/// The caller must pass unit axes. The function does not normalise them.
pub fn connect(plug: &Port, socket: &Port) -> PortFrame {
    mate_offset_m(plug, socket, 0.0)
}

/// Returns the transform that mates a plug to a socket with an extra separation.
///
/// The caller must pass unit axes. The function does not normalise them. The
/// separation adds to the socket gap and the plug gap along the socket axis.
pub fn mate_offset_m(plug: &Port, socket: &Port, separation_m: f64) -> PortFrame {
    let reach_m = socket.gap_m + plug.gap_m + separation_m;
    let point_m = [
        socket.origin_m[0] + socket.axis_m[0] * reach_m,
        socket.origin_m[1] + socket.axis_m[1] * reach_m,
        socket.origin_m[2] + socket.axis_m[2] * reach_m,
    ];
    let offset_m = [
        point_m[0] - plug.origin_m[0],
        point_m[1] - plug.origin_m[1],
        point_m[2] - plug.origin_m[2],
    ];
    let rotation_rad = mating_rotation_rad(&plug.axis_m, &socket.axis_m);
    PortFrame {
        rotation_rad,
        offset_m,
    }
}

/// Returns the rotation about z that makes the plug axis anti-parallel to the
/// socket axis.
///
/// A rotation about z cannot align a vertical axis. When either axis has zero
/// length in the xy plane, this function returns 0.0, and the caller must
/// handle a tilt about x or y. The function stays planar, as the contract
/// requires.
fn mating_rotation_rad(plug_axis_m: &[f64; 3], socket_axis_m: &[f64; 3]) -> f64 {
    let plug_xy_m = (plug_axis_m[0] * plug_axis_m[0] + plug_axis_m[1] * plug_axis_m[1]).sqrt();
    let socket_xy_m =
        (socket_axis_m[0] * socket_axis_m[0] + socket_axis_m[1] * socket_axis_m[1]).sqrt();
    if plug_xy_m == 0.0 || socket_xy_m == 0.0 {
        return 0.0;
    }
    let yaw_plug_rad = plug_axis_m[1].atan2(plug_axis_m[0]);
    let yaw_socket_rad = socket_axis_m[1].atan2(socket_axis_m[0]);
    let yaw_target_rad = wrap_to_pi(yaw_socket_rad + PI);
    wrap_to_pi(yaw_target_rad - yaw_plug_rad)
}

/// Maps an angle in radians into the range `(-PI, PI]`.
fn wrap_to_pi(angle_rad: f64) -> f64 {
    let two_pi = 2.0 * PI;
    let mut wrapped_rad = angle_rad.rem_euclid(two_pi);
    if wrapped_rad > PI {
        wrapped_rad -= two_pi;
    }
    wrapped_rad
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_plug_meets_a_socket_at_the_gap() {
        let mut socket = Port::named("socket");
        socket.gap_m = 2e-10;
        let plug = Port::named("plug");
        let frame = connect(&plug, &socket);
        assert_eq!(frame.offset_m, [0.0, 0.0, 2e-10]);
    }

    #[test]
    fn the_gaps_add() {
        let mut socket = Port::named("socket");
        socket.gap_m = 1e-10;
        let mut plug = Port::named("plug");
        plug.gap_m = 2e-10;
        let frame = connect(&plug, &socket);
        assert_eq!(frame.offset_m, [0.0, 0.0, 3e-10]);
    }

    #[test]
    fn an_extra_separation_moves_along_the_axis() {
        let socket = Port::named("socket");
        let plug = Port::named("plug");
        let frame = mate_offset_m(&plug, &socket, 5e-10);
        assert_eq!(frame.offset_m, [0.0, 0.0, 5e-10]);
    }

    #[test]
    fn an_x_axis_socket_flips_the_plug_yaw() {
        let socket = Port {
            name: "socket".to_string(),
            origin_m: [0.0, 0.0, 0.0],
            axis_m: [1.0, 0.0, 0.0],
            dof: Dof::Fixed,
            gap_m: 0.0,
        };
        let plug = Port {
            name: "plug".to_string(),
            origin_m: [0.0, 0.0, 0.0],
            axis_m: [1.0, 0.0, 0.0],
            dof: Dof::Fixed,
            gap_m: 0.0,
        };
        let frame = connect(&plug, &socket);
        assert!((frame.rotation_rad - PI).abs() < 1e-9);
    }

    #[test]
    fn a_vertical_axis_has_no_rotation() {
        let socket = Port::named("socket");
        let plug = Port::named("plug");
        let frame = connect(&plug, &socket);
        assert_eq!(frame.rotation_rad, 0.0);
    }

    #[test]
    fn the_offset_subtracts_the_plug_origin() {
        let mut socket = Port::named("socket");
        socket.gap_m = 3e-10;
        let mut plug = Port::named("plug");
        plug.origin_m = [1e-9, 0.0, 0.0];
        let frame = connect(&plug, &socket);
        assert_eq!(frame.offset_m, [-1e-9, 0.0, 3e-10]);
    }
}
