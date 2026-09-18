//! Functional groups on a diamond surface.
//!
//! A pocket that selects one molecule and not another needs more than a
//! cavity. Its wall needs chemical character: a hydroxyl group is a hydrogen
//! bond donor and acceptor, an amino group is a base, a methyl group is
//! apolar. This module builds such groups on the free tetrahedral directions
//! of a diamond surface.
//!
//! Every group attaches through one bond to a surface carbon. The bond
//! lengths come from [`nanocad_model::bond_length_m`], so the group and the
//! lattice share one source of radii. The remaining bonds of each heavy atom
//! take the tetrahedral angle, and each one carries a hydrogen, so every
//! atom reaches its usual valence.

use nanocad_model::{bond_length_m, Atom, Bond, BondType, Element, Topology};

use crate::error::PartError;

/// A group that caps one free tetrahedral direction of a surface carbon.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FunctionalGroup {
    /// A hydroxyl group, `-O-H`. A hydrogen bond donor and acceptor.
    Hydroxyl,
    /// An amino group, `-N-H2`. A hydrogen bond donor and a base.
    Amino,
    /// A methyl group, `-C-H3`. Apolar.
    Methyl,
    /// A fluorine atom, `-F`. A weak hydrogen bond acceptor.
    Fluoro,
    /// A chlorine atom, `-Cl`. Apolar and polarisable.
    Chloro,
    /// A thiol group, `-S-H`. A weak hydrogen bond donor.
    Thiol,
}

/// The tetrahedral bond angle, `acos(-1/3)`, in radians.
const TETRAHEDRAL_ANGLE_RAD: f64 = 1.910_633_149_330_218;

impl FunctionalGroup {
    /// Returns a short label for the group.
    pub fn label(self) -> &'static str {
        match self {
            Self::Hydroxyl => "hydroxyl",
            Self::Amino => "amino",
            Self::Methyl => "methyl",
            Self::Fluoro => "fluoro",
            Self::Chloro => "chloro",
            Self::Thiol => "thiol",
        }
    }

    /// Returns the heavy element of the group.
    pub fn heavy_element(self) -> Element {
        match self {
            Self::Hydroxyl => Element::OXYGEN,
            Self::Amino => Element::NITROGEN,
            Self::Methyl => Element::CARBON,
            Self::Fluoro => Element::FLUORINE,
            Self::Chloro => Element::CHLORINE,
            Self::Thiol => Element::SULFUR,
        }
    }

    /// Returns the number of hydrogens on the heavy atom.
    ///
    /// The count fills the heavy atom to its usual valence: one bond to the
    /// host, and one hydrogen for each remaining bond.
    pub fn hydrogen_count(self) -> usize {
        let valence: u8 = match self {
            Self::Hydroxyl | Self::Thiol => 2,
            Self::Amino => 3,
            Self::Methyl => 4,
            Self::Fluoro | Self::Chloro => 1,
        };
        valence.saturating_sub(1) as usize
    }

    /// Returns the atom type string of the heavy atom.
    pub fn heavy_type(self) -> &'static str {
        match self {
            Self::Hydroxyl => "O_hydroxyl",
            Self::Amino => "N_amino",
            Self::Methyl => "C_methyl",
            Self::Fluoro => "F_",
            Self::Chloro => "Cl_",
            Self::Thiol => "S_thiol",
        }
    }
}

/// Attaches one group to the host atom along `direction_m`.
///
/// The direction is a unit vector from the host towards the group. The
/// returned indices are the new atoms, the heavy atom first. A halide group
/// returns one atom and no hydrogen.
pub(crate) fn add_group(
    topology: &mut Topology,
    host: u32,
    host_position_m: [f64; 3],
    direction_m: [f64; 3],
    group: FunctionalGroup,
) -> Result<Vec<u32>, PartError> {
    let element = group.heavy_element();
    let host_element = topology
        .element(host as usize)
        .ok_or(PartError::InvalidGeometry(
            "the host atom has no element".to_owned(),
        ))?;
    let contact_m = bond_length_m(host_element, element).ok_or(PartError::InvalidGeometry(
        format!("no bond length for {host_element} and {element}"),
    ))?;
    let direction = normalize(direction_m).ok_or(PartError::InvalidGeometry(
        "the group direction has no length".to_owned(),
    ))?;
    let position_m = along(host_position_m, direction, contact_m);
    let heavy = topology.add_atom(Atom::new(element, position_m, 0.0, group.heavy_type()));
    topology.add_bond(Bond::new(host, heavy, 1, BondType::Single))?;
    let mut added = vec![heavy];

    let hydrogen_count = group.hydrogen_count();
    if hydrogen_count == 0 {
        return Ok(added);
    }
    let hydrogen_m = bond_length_m(element, Element::HYDROGEN).ok_or(
        PartError::InvalidGeometry(format!("no bond length for {element} and hydrogen")),
    )?;
    for branch in branch_directions(direction, hydrogen_count) {
        let hydrogen_m = along(position_m, branch, hydrogen_m);
        let hydrogen = topology.add_atom(Atom::new(Element::HYDROGEN, hydrogen_m, 0.0, "H"));
        topology.add_bond(Bond::new(heavy, hydrogen, 1, BondType::Single))?;
        added.push(hydrogen);
    }
    Ok(added)
}

/// Returns unit vectors for the bonds that remain on a heavy atom.
///
/// `axis` points from the host to the heavy atom. `count` bonds point away
/// from the host at the tetrahedral angle, spread evenly in azimuth.
fn branch_directions(axis: [f64; 3], count: usize) -> Vec<[f64; 3]> {
    if count == 0 {
        return Vec::new();
    }
    let (u, v) = perpendicular_frame(axis);
    let cos_angle = TETRAHEDRAL_ANGLE_RAD.cos();
    let sin_angle = TETRAHEDRAL_ANGLE_RAD.sin();
    let mut directions = Vec::with_capacity(count);
    for slot in 0..count {
        let azimuth = std::f64::consts::TAU * slot as f64 / count as f64;
        let (sin_azimuth, cos_azimuth) = azimuth.sin_cos();
        directions.push([
            sin_angle * (cos_azimuth * u[0] + sin_azimuth * v[0]) + cos_angle * axis[0],
            sin_angle * (cos_azimuth * u[1] + sin_azimuth * v[1]) + cos_angle * axis[1],
            sin_angle * (cos_azimuth * u[2] + sin_azimuth * v[2]) + cos_angle * axis[2],
        ]);
    }
    directions
}

/// Returns two unit vectors, both perpendicular to `axis` and to each other.
fn perpendicular_frame(axis: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let seed = if axis[0].abs() <= axis[1].abs() && axis[0].abs() <= axis[2].abs() {
        [1.0, 0.0, 0.0]
    } else if axis[1].abs() <= axis[2].abs() {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    let u = normalize([
        axis[1] * seed[2] - axis[2] * seed[1],
        axis[2] * seed[0] - axis[0] * seed[2],
        axis[0] * seed[1] - axis[1] * seed[0],
    ])
    .unwrap_or([1.0, 0.0, 0.0]);
    let v = [
        axis[1] * u[2] - axis[2] * u[1],
        axis[2] * u[0] - axis[0] * u[2],
        axis[0] * u[1] - axis[1] * u[0],
    ];
    (u, v)
}

fn normalize(vector: [f64; 3]) -> Option<[f64; 3]> {
    let length = (vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]).sqrt();
    if length < 1.0e-30 {
        return None;
    }
    Some([vector[0] / length, vector[1] / length, vector[2] / length])
}

fn along(origin_m: [f64; 3], direction: [f64; 3], distance_m: f64) -> [f64; 3] {
    [
        origin_m[0] + direction[0] * distance_m,
        origin_m[1] + direction[1] * distance_m,
        origin_m[2] + direction[2] * distance_m,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    const HOST_POSITION_M: [f64; 3] = [0.0, 0.0, 0.0];
    const AXIS: [f64; 3] = [1.0, 0.0, 0.0];

    fn attach(group: FunctionalGroup) -> (Topology, Vec<u32>, u32) {
        let mut topology = Topology::new();
        let host = topology.add_atom(Atom::new(Element::CARBON, HOST_POSITION_M, 0.0, "C"));
        let added = add_group(&mut topology, host, HOST_POSITION_M, AXIS, group)
            .expect("the group attaches");
        (topology, added, host)
    }

    fn distance_m(topology: &Topology, u: u32, v: u32) -> f64 {
        let a = topology.position_m(u as usize).expect("u has a position");
        let b = topology.position_m(v as usize).expect("v has a position");
        let dx = a[0] - b[0];
        let dy = a[1] - b[1];
        let dz = a[2] - b[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    #[test]
    fn every_group_reaches_its_usual_valence() {
        for group in [
            FunctionalGroup::Hydroxyl,
            FunctionalGroup::Amino,
            FunctionalGroup::Methyl,
            FunctionalGroup::Fluoro,
            FunctionalGroup::Chloro,
            FunctionalGroup::Thiol,
        ] {
            let (topology, added, host) = attach(group);
            let heavy = added[0];
            let expected = 1 + group.hydrogen_count();
            let bonded = topology
                .bonds()
                .filter(|bond| bond.u == heavy || bond.v == heavy)
                .count();
            assert_eq!(
                bonded,
                expected,
                "{} has the wrong bond count",
                group.label()
            );
            assert!(
                topology
                    .bonds()
                    .any(|bond| (bond.u == host && bond.v == heavy)
                        || (bond.v == host && bond.u == heavy)),
                "{} is not bonded to the host",
                group.label()
            );
        }
    }

    #[test]
    fn the_contact_bond_uses_the_table_length() {
        let (topology, added, host) = attach(FunctionalGroup::Hydroxyl);
        let expected_m = bond_length_m(Element::CARBON, Element::OXYGEN).expect("both are known");
        let actual_m = distance_m(&topology, host, added[0]);
        assert!(
            (actual_m - expected_m).abs() < 1.0e-18,
            "contact {actual_m}"
        );
    }

    #[test]
    fn a_hydroxyl_group_points_along_its_direction() {
        let (topology, added, _host) = attach(FunctionalGroup::Hydroxyl);
        let position_m = topology
            .position_m(added[0] as usize)
            .expect("it has a position");
        assert!(position_m[0] > 0.0, "the oxygen is behind the host");
        assert!(position_m[1].abs() < 1.0e-20 && position_m[2].abs() < 1.0e-20);
    }

    #[test]
    fn the_hydrogens_sit_at_the_tetrahedral_angle() {
        let (topology, added, _host) = attach(FunctionalGroup::Methyl);
        let heavy = added[0];
        let heavy_position_m = topology
            .position_m(heavy as usize)
            .expect("it has a position");
        let axis = normalize([
            heavy_position_m[0] - HOST_POSITION_M[0],
            heavy_position_m[1] - HOST_POSITION_M[1],
            heavy_position_m[2] - HOST_POSITION_M[2],
        ])
        .expect("the axis has a length");
        for &hydrogen in &added[1..] {
            let hydrogen_m = topology
                .position_m(hydrogen as usize)
                .expect("it has a position");
            let bond = normalize([
                hydrogen_m[0] - heavy_position_m[0],
                hydrogen_m[1] - heavy_position_m[1],
                hydrogen_m[2] - heavy_position_m[2],
            ])
            .expect("the bond has a length");
            let cosine = axis[0] * bond[0] + axis[1] * bond[1] + axis[2] * bond[2];
            let angle = cosine.clamp(-1.0, 1.0).acos();
            assert!(
                (angle - TETRAHEDRAL_ANGLE_RAD).abs() < 1.0e-9,
                "angle {angle} rad is not tetrahedral"
            );
        }
    }

    #[test]
    fn a_halide_group_has_no_hydrogen() {
        for group in [FunctionalGroup::Fluoro, FunctionalGroup::Chloro] {
            let (_topology, added, _host) = attach(group);
            assert_eq!(added.len(), 1, "{} added a hydrogen", group.label());
            assert_eq!(group.hydrogen_count(), 0);
        }
    }

    #[test]
    fn the_hydrogens_are_not_on_top_of_the_heavy_atom() {
        for group in [
            FunctionalGroup::Hydroxyl,
            FunctionalGroup::Amino,
            FunctionalGroup::Methyl,
            FunctionalGroup::Thiol,
        ] {
            let (topology, added, _host) = attach(group);
            let heavy = added[0];
            for &hydrogen in &added[1..] {
                let separation_m = distance_m(&topology, heavy, hydrogen);
                assert!(
                    separation_m > 0.9e-10,
                    "{} has a hydrogen at {separation_m} m",
                    group.label()
                );
            }
        }
    }

    #[test]
    fn a_zero_direction_is_refused() {
        let mut topology = Topology::new();
        let host = topology.add_atom(Atom::new(Element::CARBON, HOST_POSITION_M, 0.0, "C"));
        let result = add_group(
            &mut topology,
            host,
            HOST_POSITION_M,
            [0.0, 0.0, 0.0],
            FunctionalGroup::Hydroxyl,
        );
        assert!(result.is_err());
    }

    #[test]
    fn the_labels_are_distinct() {
        let labels = [
            FunctionalGroup::Hydroxyl.label(),
            FunctionalGroup::Amino.label(),
            FunctionalGroup::Methyl.label(),
            FunctionalGroup::Fluoro.label(),
            FunctionalGroup::Chloro.label(),
            FunctionalGroup::Thiol.label(),
        ];
        for (index, label) in labels.iter().enumerate() {
            assert!(!label.is_empty());
            for other in &labels[index + 1..] {
                assert_ne!(label, other);
            }
        }
    }
}
