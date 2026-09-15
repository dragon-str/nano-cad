//! Export an assembled gearbox to URDF, and parse it back.
//!
//! The exporter writes one `<link>` per device body and one `<joint>` per
//! revolute joint. The mass and the inertia come from the [`PartRecord`] of
//! the body. The link and joint frames coincide with the world frame at the
//! zero configuration, so each joint origin is the child placement relative
//! to its parent.
//!
//! The parser is a small XML reader. It reads the URDF subset that the
//! exporter writes. It is not a general XML parser.
//!
//! [`PartRecord`]: nanocad_params::PartRecord

use thiserror::Error;

use crate::assembly::{
    BodyRole, PlanetaryAssembly, CARRIER_PART_ID, PLANET_PART_ID_PREFIX, RING_PART_ID, SUN_PART_ID,
};

/// The effort limit that the exporter writes when the caller gives none.
///
/// This value is an estimate. It is not a measured gear strength and it is not
/// a claim about the real device.
pub const URDF_DEFAULT_EFFORT_N_M: f64 = 1.0e-6;
/// The velocity limit that the exporter writes when the caller gives none.
///
/// This value is an estimate. It is not a measured gear speed.
pub const URDF_DEFAULT_VELOCITY_RAD_PER_S: f64 = 1.0e13;

/// Errors from URDF export and parse.
#[derive(Debug, Error, PartialEq)]
pub enum UrdfError {
    #[error(transparent)]
    Device(#[from] crate::device::DeviceError),
    #[error("assembly error: {0}")]
    Assembly(String),
    #[error("the assembly has no part record for body {0}")]
    MissingRecord(String),
    #[error("the record for {part_id:?} is invalid: {reason}")]
    InvalidRecord { part_id: String, reason: String },
    #[error("cannot parse {field:?} as a number in {value:?}")]
    InvalidNumber { field: String, value: String },
    #[error("malformed XML: {0}")]
    Malformed(String),
    #[error("the root element is {actual:?}, not \"robot\"")]
    NotRobot { actual: String },
    #[error("unsupported joint type {0:?}")]
    UnsupportedJointType(String),
    #[error("missing attribute {attr:?} on element {element:?}")]
    MissingAttribute { element: String, attr: String },
}

/// The joint limits that the exporter writes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UrdfLimits {
    /// The effort limit in newton metres.
    pub effort_n_m: f64,
    /// The velocity limit in radians per second.
    pub velocity_rad_per_s: f64,
    /// The lower joint limit in radians.
    pub lower_rad: f64,
    /// The upper joint limit in radians.
    pub upper_rad: f64,
}

impl Default for UrdfLimits {
    fn default() -> Self {
        Self {
            effort_n_m: URDF_DEFAULT_EFFORT_N_M,
            velocity_rad_per_s: URDF_DEFAULT_VELOCITY_RAD_PER_S,
            lower_rad: -std::f64::consts::PI,
            upper_rad: std::f64::consts::PI,
        }
    }
}

/// The kind of a URDF joint.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum UrdfJointKind {
    /// A bounded revolute joint.
    Revolute,
    /// A bounded prismatic joint.
    Prismatic,
    /// A fixed joint.
    Fixed,
    /// An unbounded revolute joint.
    Continuous,
}

impl UrdfJointKind {
    /// Returns the URDF type string.
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Revolute => "revolute",
            Self::Prismatic => "prismatic",
            Self::Fixed => "fixed",
            Self::Continuous => "continuous",
        }
    }

    fn parse(value: &str) -> Result<Self, UrdfError> {
        match value {
            "revolute" => Ok(Self::Revolute),
            "prismatic" => Ok(Self::Prismatic),
            "fixed" => Ok(Self::Fixed),
            "continuous" => Ok(Self::Continuous),
            other => Err(UrdfError::UnsupportedJointType(other.to_owned())),
        }
    }
}

/// One URDF link.
#[derive(Clone, Debug, PartialEq)]
pub struct UrdfLink {
    /// The link name.
    pub name: String,
    /// The mass in kilograms. It is zero for a link with no inertia.
    pub mass_kg: f64,
    /// The three principal moments of inertia, in kilograms metre squared.
    pub inertia_kg_m2: [f64; 3],
}

/// One URDF joint.
#[derive(Clone, Debug, PartialEq)]
pub struct UrdfJoint {
    /// The joint name.
    pub name: String,
    /// The joint kind.
    pub kind: UrdfJointKind,
    /// The parent link name.
    pub parent: String,
    /// The child link name.
    pub child: String,
    /// The joint origin in the parent frame, in metres.
    pub origin_m: [f64; 3],
    /// The joint axis in the child frame.
    pub axis: [f64; 3],
    /// The joint limits.
    pub limits: UrdfLimits,
}

/// A parsed URDF robot.
#[derive(Clone, Debug, PartialEq)]
pub struct UrdfRobot {
    /// The robot name.
    pub name: String,
    /// The links.
    pub links: Vec<UrdfLink>,
    /// The joints.
    pub joints: Vec<UrdfJoint>,
}

/// The link name of the fixed world body.
pub const GROUND_LINK: &str = "ground";
/// The link name of the sun gear.
pub const SUN_LINK: &str = "sun";
/// The link name prefix of the planet gears.
pub const PLANET_LINK_PREFIX: &str = "planet_";
/// The link name of the ring gear.
pub const RING_LINK: &str = "ring";
/// The link name of the carrier.
pub const CARRIER_LINK: &str = "carrier";

fn link_name(role: BodyRole) -> String {
    match role {
        BodyRole::Ground => GROUND_LINK.to_owned(),
        BodyRole::Sun => SUN_LINK.to_owned(),
        BodyRole::Planet(index) => format!("{PLANET_LINK_PREFIX}{index}"),
        BodyRole::Ring => RING_LINK.to_owned(),
        BodyRole::Carrier => CARRIER_LINK.to_owned(),
    }
}

fn record_part_id(role: BodyRole) -> Option<String> {
    match role {
        BodyRole::Ground => None,
        BodyRole::Sun => Some(SUN_PART_ID.to_owned()),
        BodyRole::Planet(index) => Some(format!("{PLANET_PART_ID_PREFIX}{index}")),
        BodyRole::Ring => Some(RING_PART_ID.to_owned()),
        BodyRole::Carrier => Some(CARRIER_PART_ID.to_owned()),
    }
}

fn joint_name(role: BodyRole) -> String {
    match role {
        BodyRole::Planet(index) => format!("planet_joint_{index}"),
        BodyRole::Sun => "sun_joint".to_owned(),
        BodyRole::Ring => "ring_joint".to_owned(),
        BodyRole::Carrier => "carrier_joint".to_owned(),
        BodyRole::Ground => "ground_joint".to_owned(),
    }
}

/// Builds the URDF model from an assembly. It does not format it.
pub fn urdf_robot(
    assembly: &PlanetaryAssembly,
    name: &str,
    limits: UrdfLimits,
) -> Result<UrdfRobot, UrdfError> {
    if !limits.effort_n_m.is_finite()
        || !limits.velocity_rad_per_s.is_finite()
        || !limits.lower_rad.is_finite()
        || !limits.upper_rad.is_finite()
    {
        return Err(UrdfError::InvalidRecord {
            part_id: name.to_owned(),
            reason: "joint limits must be finite".to_owned(),
        });
    }
    if limits.lower_rad > limits.upper_rad {
        return Err(UrdfError::InvalidRecord {
            part_id: name.to_owned(),
            reason: "lower limit is above the upper limit".to_owned(),
        });
    }

    let mut links = Vec::with_capacity(assembly.roles.len());
    for role in &assembly.roles {
        let part_id = record_part_id(*role);
        let (mass_kg, inertia_kg_m2) = match part_id {
            None => (0.0, [0.0, 0.0, 0.0]),
            Some(part_id) => {
                let record = assembly
                    .library
                    .get(&part_id)
                    .ok_or_else(|| UrdfError::MissingRecord(part_id.clone()))?;
                let mass_kg = record.mass_kg.value_si();
                let inertia = [
                    record.inertia_kg_m2[0].value_si(),
                    record.inertia_kg_m2[1].value_si(),
                    record.inertia_kg_m2[2].value_si(),
                ];
                if !mass_kg.is_finite() || mass_kg <= 0.0 {
                    return Err(UrdfError::InvalidRecord {
                        part_id,
                        reason: "mass must be finite and positive".to_owned(),
                    });
                }
                for value in inertia {
                    if !value.is_finite() || value <= 0.0 {
                        return Err(UrdfError::InvalidRecord {
                            part_id,
                            reason: "inertia must be finite and positive".to_owned(),
                        });
                    }
                }
                (mass_kg, inertia)
            }
        };
        links.push(UrdfLink {
            name: link_name(*role),
            mass_kg,
            inertia_kg_m2,
        });
    }

    let mut joints = Vec::new();
    for role in &assembly.roles {
        if *role == BodyRole::Ground {
            continue;
        }
        let (parent, _child) = match role {
            BodyRole::Sun => (assembly.ground_body, assembly.sun_body),
            BodyRole::Carrier => (assembly.ground_body, assembly.carrier_body),
            BodyRole::Ring => (assembly.ground_body, assembly.ring_body),
            BodyRole::Planet(index) => (assembly.carrier_body, assembly.planet_bodies[*index]),
            BodyRole::Ground => continue,
        };
        let joint = assembly
            .system
            .revolute_joints()
            .get(revolute_index(assembly, *role))
            .ok_or_else(|| UrdfError::Assembly("joint is not available".to_owned()))?;
        let anchor_m = joint.anchor_a_world_m(assembly.system.bodies())?;
        let parent_position_m = assembly
            .system
            .body(parent)
            .ok_or_else(|| UrdfError::Assembly("parent body is not available".to_owned()))?
            .position_m();
        let origin_m = [
            anchor_m[0] - parent_position_m[0],
            anchor_m[1] - parent_position_m[1],
            anchor_m[2] - parent_position_m[2],
        ];
        joints.push(UrdfJoint {
            name: joint_name(*role),
            kind: UrdfJointKind::Revolute,
            parent: link_name(role_of(assembly, parent)),
            child: link_name(*role),
            origin_m,
            axis: [0.0, 0.0, 1.0],
            limits,
        });
    }

    Ok(UrdfRobot {
        name: name.to_owned(),
        links,
        joints,
    })
}

fn revolute_index(assembly: &PlanetaryAssembly, role: BodyRole) -> usize {
    match role {
        BodyRole::Sun => assembly.sun_joint,
        BodyRole::Carrier => assembly.carrier_joint,
        BodyRole::Ring => assembly.ring_joint,
        BodyRole::Planet(index) => assembly.planet_joints[index],
        BodyRole::Ground => 0,
    }
}

fn role_of(assembly: &PlanetaryAssembly, body: usize) -> BodyRole {
    assembly
        .roles
        .get(body)
        .copied()
        .unwrap_or(BodyRole::Ground)
}

/// Formats a URDF model as XML text.
pub fn format_urdf(robot: &UrdfRobot) -> String {
    let mut out = String::new();
    out.push_str("<?xml version=\"1.0\"?>\n");
    out.push_str(&format!("<robot name=\"{}\">\n", escape(&robot.name)));
    for link in &robot.links {
        if link.mass_kg > 0.0 {
            out.push_str(&format!("  <link name=\"{}\">\n", escape(&link.name)));
            out.push_str("    <inertial>\n");
            out.push_str("      <origin xyz=\"0 0 0\" rpy=\"0 0 0\"/>\n");
            out.push_str(&format!(
                "      <mass value=\"{}\"/>\n",
                format_number(link.mass_kg)
            ));
            out.push_str(&format!(
                "      <inertia ixx=\"{}\" ixy=\"0\" ixz=\"0\" iyy=\"{}\" iyz=\"0\" izz=\"{}\"/>\n",
                format_number(link.inertia_kg_m2[0]),
                format_number(link.inertia_kg_m2[1]),
                format_number(link.inertia_kg_m2[2]),
            ));
            out.push_str("    </inertial>\n");
            out.push_str("  </link>\n");
        } else {
            out.push_str(&format!("  <link name=\"{}\"/>\n", escape(&link.name)));
        }
    }
    for joint in &robot.joints {
        out.push_str(&format!(
            "  <joint name=\"{}\" type=\"{}\">\n",
            escape(&joint.name),
            joint.kind.as_str()
        ));
        out.push_str(&format!(
            "    <parent link=\"{}\"/>\n",
            escape(&joint.parent)
        ));
        out.push_str(&format!("    <child link=\"{}\"/>\n", escape(&joint.child)));
        out.push_str(&format!(
            "    <origin xyz=\"{} {} {}\" rpy=\"0 0 0\"/>\n",
            format_number(joint.origin_m[0]),
            format_number(joint.origin_m[1]),
            format_number(joint.origin_m[2]),
        ));
        out.push_str(&format!(
            "    <axis xyz=\"{} {} {}\"/>\n",
            format_number(joint.axis[0]),
            format_number(joint.axis[1]),
            format_number(joint.axis[2]),
        ));
        if joint.kind != UrdfJointKind::Fixed {
            out.push_str(&format!(
                "    <limit lower=\"{}\" upper=\"{}\" effort=\"{}\" velocity=\"{}\"/>\n",
                format_number(joint.limits.lower_rad),
                format_number(joint.limits.upper_rad),
                format_number(joint.limits.effort_n_m),
                format_number(joint.limits.velocity_rad_per_s),
            ));
        }
        out.push_str("  </joint>\n");
    }
    out.push_str("</robot>\n");
    out
}

/// Exports the assembly as a URDF string with the default robot name.
pub fn export_urdf(assembly: &PlanetaryAssembly) -> Result<String, UrdfError> {
    export_urdf_named(assembly, "planetary_gearbox")
}

/// Exports the assembly as a URDF string with a chosen robot name.
pub fn export_urdf_named(assembly: &PlanetaryAssembly, name: &str) -> Result<String, UrdfError> {
    export_urdf_with_limits(assembly, name, UrdfLimits::default())
}

/// Exports the assembly as a URDF string with chosen joint limits.
pub fn export_urdf_with_limits(
    assembly: &PlanetaryAssembly,
    name: &str,
    limits: UrdfLimits,
) -> Result<String, UrdfError> {
    let robot = urdf_robot(assembly, name, limits)?;
    Ok(format_urdf(&robot))
}

fn escape(value: &str) -> String {
    value
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

fn format_number(value: f64) -> String {
    if value == 0.0 {
        return "0".to_owned();
    }
    let mut text = format!("{value:.17e}");
    if let Some(dot) = text.find('.') {
        let exponent = text.find('e').unwrap_or(text.len());
        let mut end = exponent;
        while end > dot + 1 && text.as_bytes()[end - 1] == b'0' {
            end -= 1;
        }
        if end == dot + 1 {
            end = dot;
        }
        text.replace_range(end..exponent, "");
    }
    text
}

#[derive(Clone, Debug)]
struct Element {
    name: String,
    attrs: Vec<(String, String)>,
    children: Vec<Element>,
}

impl Element {
    fn attr(&self, name: &str) -> Option<&str> {
        self.attrs
            .iter()
            .find(|(key, _)| key == name)
            .map(|(_, value)| value.as_str())
    }

    fn require_attr(&self, name: &str) -> Result<&str, UrdfError> {
        self.attr(name).ok_or_else(|| UrdfError::MissingAttribute {
            element: self.name.clone(),
            attr: name.to_owned(),
        })
    }

    fn child(&self, name: &str) -> Option<&Element> {
        self.children.iter().find(|child| child.name == name)
    }
}

fn skip_trivia(xml: &[u8], pos: &mut usize) -> Result<(), UrdfError> {
    loop {
        while matches!(xml.get(*pos), Some(b' ' | b'\t' | b'\n' | b'\r')) {
            *pos += 1;
        }
        if xml.get(*pos) == Some(&b'<') && xml.get(*pos + 1) == Some(&b'?') {
            while *pos < xml.len() && xml.get(*pos) != Some(&b'>') {
                *pos += 1;
            }
            *pos += 1;
            continue;
        }
        if xml.get(*pos) == Some(&b'<') && xml.get(*pos + 1) == Some(&b'!') {
            while *pos + 1 < xml.len() && !(xml[*pos] == b'-' && xml[*pos + 1] == b'>') {
                *pos += 1;
            }
            *pos += 2;
            continue;
        }
        return Ok(());
    }
}

fn parse_name(xml: &[u8], pos: &mut usize) -> Result<String, UrdfError> {
    let start = *pos;
    while let Some(byte) = xml.get(*pos) {
        let ok = byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'-' | b'.' | b':');
        if !ok {
            break;
        }
        *pos += 1;
    }
    if *pos == start {
        return Err(UrdfError::Malformed("expected a name".to_owned()));
    }
    std::str::from_utf8(&xml[start..*pos])
        .map(str::to_owned)
        .map_err(|_| UrdfError::Malformed("the name is not valid UTF-8".to_owned()))
}

fn parse_quoted(xml: &[u8], pos: &mut usize) -> Result<String, UrdfError> {
    let quote = *xml
        .get(*pos)
        .ok_or_else(|| UrdfError::Malformed("expected an attribute value".to_owned()))?;
    if quote != b'"' && quote != b'\'' {
        return Err(UrdfError::Malformed(
            "attribute values must be quoted".to_owned(),
        ));
    }
    *pos += 1;
    let start = *pos;
    while xml.get(*pos).is_some() && xml[*pos] != quote {
        *pos += 1;
    }
    let value = std::str::from_utf8(&xml[start..*pos])
        .map_err(|_| UrdfError::Malformed("attribute value is not valid UTF-8".to_owned()))?
        .to_owned();
    *pos += 1;
    Ok(value)
}

fn parse_element(xml: &[u8], pos: &mut usize) -> Result<Element, UrdfError> {
    skip_trivia(xml, pos)?;
    if xml.get(*pos) != Some(&b'<') {
        return Err(UrdfError::Malformed("expected an element".to_owned()));
    }
    *pos += 1;
    let name = parse_name(xml, pos)?;
    let mut attrs = Vec::new();
    loop {
        while matches!(xml.get(*pos), Some(b' ' | b'\t' | b'\n' | b'\r')) {
            *pos += 1;
        }
        match xml.get(*pos) {
            Some(b'/') => {
                *pos += 1;
                if xml.get(*pos) != Some(&b'>') {
                    return Err(UrdfError::Malformed("expected \">\"".to_owned()));
                }
                *pos += 1;
                return Ok(Element {
                    name,
                    attrs,
                    children: Vec::new(),
                });
            }
            Some(b'>') => {
                *pos += 1;
                break;
            }
            Some(_) => {
                let key = parse_name(xml, pos)?;
                while matches!(xml.get(*pos), Some(b' ' | b'\t' | b'\n' | b'\r')) {
                    *pos += 1;
                }
                if xml.get(*pos) != Some(&b'=') {
                    return Err(UrdfError::Malformed("expected \"=\"".to_owned()));
                }
                *pos += 1;
                while matches!(xml.get(*pos), Some(b' ' | b'\t' | b'\n' | b'\r')) {
                    *pos += 1;
                }
                let value = parse_quoted(xml, pos)?;
                attrs.push((key, value));
            }
            None => return Err(UrdfError::Malformed("unterminated start tag".to_owned())),
        }
    }

    let mut children = Vec::new();
    loop {
        skip_trivia(xml, pos)?;
        match xml.get(*pos) {
            Some(b'<') if xml.get(*pos + 1) == Some(&b'/') => {
                *pos += 2;
                let close = parse_name(xml, pos)?;
                skip_trivia(xml, pos)?;
                if xml.get(*pos) != Some(&b'>') {
                    return Err(UrdfError::Malformed("expected \">\"".to_owned()));
                }
                *pos += 1;
                if close != name {
                    return Err(UrdfError::Malformed(format!(
                        "closing tag {close:?} does not match {name:?}"
                    )));
                }
                return Ok(Element {
                    name,
                    attrs,
                    children,
                });
            }
            Some(b'<') => children.push(parse_element(xml, pos)?),
            Some(_) => *pos += 1,
            None => {
                return Err(UrdfError::Malformed(format!(
                    "unterminated element {name:?}"
                )))
            }
        }
    }
}

fn parse_f64(text: &str, field: &str) -> Result<f64, UrdfError> {
    text.trim()
        .parse::<f64>()
        .map_err(|_| UrdfError::InvalidNumber {
            field: field.to_owned(),
            value: text.to_owned(),
        })
}

fn parse_vector(text: &str, field: &str) -> Result<[f64; 3], UrdfError> {
    let mut values = [0.0; 3];
    let mut count = 0usize;
    for part in text.split_whitespace() {
        if count >= 3 {
            return Err(UrdfError::InvalidNumber {
                field: field.to_owned(),
                value: text.to_owned(),
            });
        }
        values[count] = parse_f64(part, field)?;
        count += 1;
    }
    if count != 3 {
        return Err(UrdfError::InvalidNumber {
            field: field.to_owned(),
            value: text.to_owned(),
        });
    }
    Ok(values)
}

/// Parses a URDF string that the exporter wrote.
pub fn parse_urdf(xml: &str) -> Result<UrdfRobot, UrdfError> {
    let bytes = xml.as_bytes();
    let mut pos = 0usize;
    let root = parse_element(bytes, &mut pos)?;
    if root.name != "robot" {
        return Err(UrdfError::NotRobot { actual: root.name });
    }
    let name = root
        .attr("name")
        .ok_or_else(|| UrdfError::MissingAttribute {
            element: "robot".to_owned(),
            attr: "name".to_owned(),
        })?
        .to_owned();

    let mut links = Vec::new();
    let mut joints = Vec::new();
    for child in &root.children {
        match child.name.as_str() {
            "link" => {
                let name = child.require_attr("name")?.to_owned();
                if let Some(inertial) = child.child("inertial") {
                    let mass = inertial.child("mass").ok_or_else(|| {
                        UrdfError::Malformed(format!("link {name:?} has no mass"))
                    })?;
                    let mass_kg = parse_f64(mass.require_attr("value")?, "mass")?;
                    let inertia = inertial.child("inertia").ok_or_else(|| {
                        UrdfError::Malformed(format!("link {name:?} has no inertia"))
                    })?;
                    let inertia_kg_m2 = [
                        parse_f64(inertia.require_attr("ixx")?, "ixx")?,
                        parse_f64(inertia.require_attr("iyy")?, "iyy")?,
                        parse_f64(inertia.require_attr("izz")?, "izz")?,
                    ];
                    links.push(UrdfLink {
                        name,
                        mass_kg,
                        inertia_kg_m2,
                    });
                } else {
                    links.push(UrdfLink {
                        name,
                        mass_kg: 0.0,
                        inertia_kg_m2: [0.0, 0.0, 0.0],
                    });
                }
            }
            "joint" => {
                let name = child.require_attr("name")?.to_owned();
                let kind = UrdfJointKind::parse(child.require_attr("type")?)?;
                let parent = child
                    .child("parent")
                    .ok_or_else(|| UrdfError::Malformed(format!("joint {name:?} has no parent")))?
                    .require_attr("link")?
                    .to_owned();
                let child_link = child
                    .child("child")
                    .ok_or_else(|| UrdfError::Malformed(format!("joint {name:?} has no child")))?
                    .require_attr("link")?
                    .to_owned();
                let origin_m = match child.child("origin") {
                    Some(origin) => parse_vector(origin.attr("xyz").unwrap_or("0 0 0"), "xyz")?,
                    None => [0.0, 0.0, 0.0],
                };
                let axis = match child.child("axis") {
                    Some(axis) => parse_vector(axis.attr("xyz").unwrap_or("0 0 1"), "xyz")?,
                    None => [0.0, 0.0, 1.0],
                };
                let limits = match child.child("limit") {
                    Some(limit) => UrdfLimits {
                        lower_rad: parse_f64(
                            limit.attr("lower").unwrap_or("-3.141592653589793"),
                            "lower",
                        )?,
                        upper_rad: parse_f64(
                            limit.attr("upper").unwrap_or("3.141592653589793"),
                            "upper",
                        )?,
                        effort_n_m: parse_f64(limit.attr("effort").unwrap_or("0"), "effort")?,
                        velocity_rad_per_s: parse_f64(
                            limit.attr("velocity").unwrap_or("0"),
                            "velocity",
                        )?,
                    },
                    None => UrdfLimits {
                        effort_n_m: 0.0,
                        velocity_rad_per_s: 0.0,
                        lower_rad: 0.0,
                        upper_rad: 0.0,
                    },
                };
                joints.push(UrdfJoint {
                    name,
                    kind,
                    parent,
                    child: child_link,
                    origin_m,
                    axis,
                    limits,
                });
            }
            _ => {}
        }
    }

    Ok(UrdfRobot {
        name,
        links,
        joints,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::assemble_planetary;
    use nanocad_parts::ParameterSet;

    fn assembly() -> PlanetaryAssembly {
        assemble_planetary(&ParameterSet::new().with("planet_count", 3.0)).expect("assemble")
    }

    #[test]
    fn the_export_has_one_link_per_body_and_one_joint_per_revolute() {
        let assembly = assembly();
        let robot = urdf_robot(&assembly, "gearbox", UrdfLimits::default()).expect("robot");
        assert_eq!(robot.name, "gearbox");
        assert_eq!(robot.links.len(), 7);
        assert_eq!(robot.joints.len(), 6);
        let names: Vec<&str> = robot.links.iter().map(|link| link.name.as_str()).collect();
        assert_eq!(
            names,
            vec!["ground", "sun", "planet_0", "planet_1", "planet_2", "ring", "carrier"]
        );
        let sun = robot
            .links
            .iter()
            .find(|link| link.name == "sun")
            .expect("sun");
        assert!(sun.mass_kg > 0.0);
        assert!(sun.inertia_kg_m2[2] > 0.0);
    }

    #[test]
    fn the_export_round_trips_through_the_parser() {
        let assembly = assembly();
        let robot = urdf_robot(&assembly, "gearbox", UrdfLimits::default()).expect("robot");
        let xml = format_urdf(&robot);
        let parsed = parse_urdf(&xml).expect("parse");
        assert_eq!(parsed, robot);
    }

    #[test]
    fn a_planet_joint_origin_is_the_planet_centre_in_the_carrier_frame() {
        let assembly = assembly();
        let robot = urdf_robot(&assembly, "gearbox", UrdfLimits::default()).expect("robot");
        let joint = robot
            .joints
            .iter()
            .find(|joint| joint.name == "planet_joint_0")
            .expect("planet joint");
        assert_eq!(joint.parent, "carrier");
        assert_eq!(joint.child, "planet_0");
        let expected = assembly.design.planet_center_m(0);
        assert!((joint.origin_m[0] - expected[0]).abs() < 1.0e-18);
        assert!((joint.origin_m[1] - expected[1]).abs() < 1.0e-18);
        assert_eq!(joint.origin_m[2], 0.0);
    }

    #[test]
    fn the_parser_rejects_a_non_robot_root() {
        assert!(matches!(
            parse_urdf("<world name=\"x\"/>"),
            Err(UrdfError::NotRobot { .. })
        ));
    }

    #[test]
    fn the_parser_rejects_an_unknown_joint_type() {
        let xml = "<robot name=\"x\"><joint name=\"j\" type=\"worm\">\
                   <parent link=\"a\"/><child link=\"b\"/></joint></robot>";
        assert!(matches!(
            parse_urdf(xml),
            Err(UrdfError::UnsupportedJointType(_))
        ));
    }

    #[test]
    fn export_rejects_bad_limits() {
        let assembly = assembly();
        let limits = UrdfLimits {
            lower_rad: 1.0,
            upper_rad: -1.0,
            ..UrdfLimits::default()
        };
        assert!(matches!(
            export_urdf_with_limits(&assembly, "x", limits),
            Err(UrdfError::InvalidRecord { .. })
        ));
    }
}
