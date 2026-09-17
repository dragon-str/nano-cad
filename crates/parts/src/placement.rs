//! Interactive part placement by port mating.
//!
//! A [`PlacedPart`] holds a part and the port frame that mates one of its ports
//! to a socket port. Applying the frame rotates the part about the z axis
//! through the origin, then translates it. The crate holds SI positions and
//! unit-suffixed names, following ADR-0003.

use nanocad_model::{Document, Part};

use crate::port::{mate_offset_m, Port, PortFrame};

/// A part together with the frame that places it on a socket.
#[derive(Clone, Debug, PartialEq)]
pub struct PlacedPart {
    /// The placed part name.
    pub name: String,
    /// The part geometry, in its own frame.
    pub part: Part,
    /// The frame that brings the plug port onto the socket port.
    pub frame: PortFrame,
}

/// Places a part by mating `plug` to `socket` with an extra separation.
///
/// The frame comes from [`mate_offset_m`]. The separation is in metres.
pub fn place(part: Part, plug: &Port, socket: &Port, separation_m: f64) -> PlacedPart {
    let name = part.name.clone();
    let frame = mate_offset_m(plug, socket, separation_m);
    PlacedPart { name, part, frame }
}

impl PlacedPart {
    /// Returns the part geometry moved into the socket frame.
    ///
    /// The returned part keeps the name, the material, and the metadata.
    pub fn transformed_part(&self) -> Part {
        let mut part = self.part.clone();
        let count = part.topology.atom_count();
        let mut moved_m = Vec::with_capacity(count);
        for index in 0..count {
            if let Some(position_m) = part.topology.position_m(index) {
                moved_m.push(self.transformed_position_m(position_m));
            }
        }
        let (positions, _tail) = part.topology.positions_m_mut().as_chunks_mut::<3>();
        for (chunk, position_m) in positions.iter_mut().zip(moved_m) {
            *chunk = position_m;
        }
        part
    }

    /// Returns one point moved into the socket frame.
    ///
    /// The point rotates about the z axis through the origin by the frame
    /// rotation, then takes the frame offset.
    pub fn transformed_position_m(&self, point_m: [f64; 3]) -> [f64; 3] {
        let [x_m, y_m, z_m] = point_m;
        let (sin_rad, cos_rad) = self.frame.rotation_rad.sin_cos();
        let rotated_m = [
            x_m * cos_rad - y_m * sin_rad,
            x_m * sin_rad + y_m * cos_rad,
            z_m,
        ];
        [
            rotated_m[0] + self.frame.offset_m[0],
            rotated_m[1] + self.frame.offset_m[1],
            rotated_m[2] + self.frame.offset_m[2],
        ]
    }
}

/// Builds a document that holds every placed part, in order.
///
/// The document is named `"assembly"`. Each part is stored in its socket frame.
pub fn assembly_document(placed: &[PlacedPart]) -> Document {
    let mut document = Document::new("assembly");
    for placed_part in placed {
        document.add_part(placed_part.transformed_part());
    }
    document
}

#[cfg(test)]
mod tests {
    use super::*;
    use nanocad_model::{Atom, Element, Topology};

    fn one_atom_part(name: &str, material: &str, position_m: [f64; 3]) -> Part {
        let mut topology = Topology::new();
        topology.add_atom(Atom::new(Element::CARBON, position_m, 0.0, "C"));
        Part::new(name, topology).with_material(material)
    }

    #[test]
    fn a_placed_part_moves_to_the_socket() {
        let part = one_atom_part("peg", "diamondoid", [1e-10, 0.0, 0.0]);
        let plug = Port::named("plug");
        let socket = Port::named("socket");
        let placed = place(part, &plug, &socket, 2e-10);
        assert_eq!(placed.name, "peg");
        assert_eq!(placed.frame.offset_m, [0.0, 0.0, 2e-10]);
        let moved_m = placed.transformed_position_m([1e-10, 0.0, 0.0]);
        assert_eq!(moved_m, [1e-10, 0.0, 2e-10]);
    }

    #[test]
    fn a_rotation_about_z_moves_the_plug() {
        let part = one_atom_part("peg", "diamondoid", [1e-10, 0.0, 0.0]);
        let plug = Port {
            name: "plug".to_string(),
            origin_m: [0.0, 0.0, 0.0],
            axis_m: [1.0, 0.0, 0.0],
            dof: crate::port::Dof::Fixed,
            gap_m: 0.0,
        };
        let socket = Port {
            name: "socket".to_string(),
            origin_m: [0.0, 0.0, 0.0],
            axis_m: [1.0, 0.0, 0.0],
            dof: crate::port::Dof::Fixed,
            gap_m: 0.0,
        };
        let placed = place(part, &plug, &socket, 0.0);
        assert!((placed.frame.rotation_rad.abs() - std::f64::consts::PI).abs() < 1e-9);
        let moved_m = placed.transformed_position_m([1e-10, 0.0, 0.0]);
        assert!((moved_m[0] + 1e-10).abs() < 1e-21, "x {}", moved_m[0]);
        assert!(moved_m[1].abs() < 1e-25, "y {}", moved_m[1]);
        assert_eq!(moved_m[2], 0.0);
    }

    #[test]
    fn the_transformed_part_keeps_its_name_and_material() {
        let part = one_atom_part("washer", "diamondoid", [1e-10, 0.0, 0.0]);
        let plug = Port::named("plug");
        let socket = Port::named("socket");
        let placed = place(part, &plug, &socket, 3e-10);
        let transformed = placed.transformed_part();
        assert_eq!(transformed.name, "washer");
        assert_eq!(transformed.material, "diamondoid");
        assert_eq!(transformed.atom_count(), 1);
        assert_eq!(
            transformed.topology.position_m(0),
            Some([1e-10, 0.0, 3e-10])
        );
    }

    #[test]
    fn an_assembly_holds_every_placed_part() {
        let first = one_atom_part("a", "diamondoid", [1e-10, 0.0, 0.0]);
        let second = one_atom_part("b", "diamondoid", [0.0, 1e-10, 0.0]);
        let plug = Port::named("plug");
        let socket = Port::named("socket");
        let placed = [
            place(first, &plug, &socket, 1e-10),
            place(second, &plug, &socket, 2e-10),
        ];
        let document = assembly_document(&placed);
        assert_eq!(document.name, "assembly");
        assert_eq!(document.part_count(), 2);
        assert_eq!(document.part(0).map(|part| part.name.as_str()), Some("a"));
        assert_eq!(document.part(1).map(|part| part.name.as_str()), Some("b"));
    }
}
