//! The binding pocket.
//!
//! A binding pocket is a cylindrical well cut into a diamond block. The well
//! wall carries one kind of functional group, so the pocket selects by shape
//! and by chemistry together.
//!
//! The wall decoration is geometric. A free tetrahedral direction decorates
//! when a short step along it lands inside the cavity, because that direction
//! points at the space the guest occupies. Every other free direction takes a
//! hydrogen, so every carbon keeps a valence of four.
//!
//! The pocket is a model. The guest geometry and the guest charges come from
//! `guest_data`, and the nonbonded parameters are UFF. No result from this
//! module is validated against experiment. See ADR-0063 and ADR-0064.
//!
//! All lengths are SI metres and all angles are radians.

use nanocad_model::{pair_params, vdw_distance_m, Atom, Element, Part, Topology};
use nanocad_units::Unit;

use crate::diamond_solid;
use crate::error::PartError;
use crate::generator::PartGenerator;
use crate::group::FunctionalGroup;
use crate::guest::Guest;
use crate::lattice_fill::fill_solid;
use crate::parameter::{ParameterSet, ParameterSpec};
use crate::port::Port;
use crate::shape::{Box3, Cylinder, Difference, Solid};

/// The longest bond in diamond, in metres. A probe shorter than this stays
/// inside the cell that owns the direction.
const PROBE_STEP_M: f64 = 0.6 * 1.544e-10;

/// The parameter specs of [`BindingPocketGenerator`], in a stable order.
static BINDING_POCKET_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new(
        "pocket_radius_m",
        Some(Unit::Metre),
        1.0e-9,
        1.5e-10,
        4.0e-9,
        false,
        "radius of the cylindrical well in metres",
    ),
    ParameterSpec::new(
        "pocket_depth_m",
        Some(Unit::Metre),
        1.2e-9,
        1.5e-10,
        8.0e-9,
        false,
        "depth of the well from the top face in metres",
    ),
    ParameterSpec::new(
        "wall_m",
        Some(Unit::Metre),
        0.8e-9,
        1.5e-10,
        4.0e-9,
        false,
        "diamond thickness around the well in metres",
    ),
    ParameterSpec::new(
        "thickness_m",
        Some(Unit::Metre),
        2.0e-9,
        4.0e-10,
        1.0e-8,
        false,
        "height of the block along z in metres",
    ),
    ParameterSpec::new(
        "wall_group",
        None,
        0.0,
        0.0,
        5.0,
        true,
        "index of the functional group on the wall, from FunctionalGroup::index",
    ),
];

/// Builds a diamond block with a decorated cylindrical well.
#[derive(Clone, Copy, Debug, Default)]
pub struct BindingPocketGenerator;

impl BindingPocketGenerator {
    /// Returns the mouth of the well as a port.
    ///
    /// The port sits on the top face, on the well axis, and it points along
    /// `+z` and out of the pocket, so a part that enters along the axis meets
    /// the well.
    pub fn opening_port(&self) -> Port {
        let thickness_m = self.default_m("thickness_m");
        let mut port = Port::named("opening");
        port.origin_m = [0.0, 0.0, 0.5 * thickness_m];
        port
    }

    /// Returns the wall group that the default parameters select.
    pub fn wall_group(&self) -> FunctionalGroup {
        FunctionalGroup::from_index(self.default_m("wall_group") as usize)
            .unwrap_or(FunctionalGroup::Hydroxyl)
    }

    /// Returns the centre of the well in the part frame, in metres.
    pub fn cavity_center_m(&self) -> [f64; 3] {
        let thickness_m = self.default_m("thickness_m");
        let depth_m = self.default_m("pocket_depth_m");
        [0.0, 0.0, 0.5 * thickness_m - 0.5 * depth_m]
    }

    /// Returns the clearance a guest has inside the default well, in metres.
    ///
    /// The clearance is the well radius minus the guest extent. A negative
    /// clearance means the guest does not fit in the well.
    pub fn guest_clearance_m(&self, guest: &Guest) -> f64 {
        self.default_m("pocket_radius_m") - guest.extent_m()
    }

    fn default_m(&self, name: &str) -> f64 {
        self.resolve(&ParameterSet::new())
            .ok()
            .and_then(|resolved| resolved.get(name))
            .unwrap_or(0.0)
    }
}

/// Returns the smallest well radius that holds `guest` with `clearance_m` to
/// spare.
pub fn pocket_radius_for(guest: &Guest, clearance_m: f64) -> f64 {
    guest.extent_m() + clearance_m
}

/// Returns the contact distance between two carbon walls, in metres.
///
/// This is the UFF nonbonded minimum of carbon, so it is the gap at which two
/// diamond walls touch. A guest needs at least this much room on each side.
pub fn wall_contact_distance_m() -> Option<f64> {
    vdw_distance_m(Element::CARBON)
}

/// Returns the well depth of the wall-guest pair, in joules.
///
/// The guest contributes its heaviest element, so the value states the
/// strongest wall attraction the guest can make. The mixing rule is the UFF
/// rule, so the value is `sqrt(D_wall D_guest) / 4`.
pub fn wall_well_depth_j(guest: &Guest) -> Option<f64> {
    let heaviest = guest
        .atoms
        .iter()
        .map(|atom| atom.element)
        .max_by_key(|element| element.atomic_number())?;
    Some(pair_params(Element::CARBON, heaviest)?.0)
}

impl PartGenerator for BindingPocketGenerator {
    fn id(&self) -> &'static str {
        "binding_pocket"
    }

    fn name(&self) -> &'static str {
        "binding pocket"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        BINDING_POCKET_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        let resolved = self.resolve(parameters)?;
        let pocket_radius_m = resolved.require("pocket_radius_m")?;
        let pocket_depth_m = resolved.require("pocket_depth_m")?;
        let wall_m = resolved.require("wall_m")?;
        let thickness_m = resolved.require("thickness_m")?;
        let group_index = resolved.require("wall_group")? as usize;

        let group = FunctionalGroup::from_index(group_index).ok_or_else(|| {
            PartError::InvalidGeometry(format!(
                "the wall group index {group_index} names no functional group"
            ))
        })?;
        if pocket_depth_m > thickness_m {
            return Err(PartError::InvalidGeometry(format!(
                "the well depth {pocket_depth_m} m is deeper than the block \
                 thickness {thickness_m} m, so the well would cut through the floor"
            )));
        }

        let half_m = pocket_radius_m + wall_m;
        let block = Box3::from_center_m([0.0, 0.0, 0.0], [2.0 * half_m, 2.0 * half_m, thickness_m]);
        let cavity = Cylinder {
            center_m: [0.0, 0.0, 0.5 * thickness_m - 0.5 * pocket_depth_m],
            radius_m: pocket_radius_m,
            height_m: pocket_depth_m,
        };
        let solid = Difference {
            outer: Box::new(block),
            inner: Box::new(cavity),
        };

        let sites_m = fill_solid(&solid);
        let mut topology = Topology::new();
        let mut carbons = Vec::with_capacity(sites_m.len());
        for site_m in sites_m {
            carbons.push(topology.add_atom(Atom::new(Element::CARBON, site_m, 0.0, "C")));
        }
        if carbons.is_empty() {
            return Err(PartError::InvalidGeometry(format!(
                "the pocket cut no atom from the diamond lattice: the well \
                 radius {pocket_radius_m} m with the wall {wall_m} m fills the block"
            )));
        }

        diamond_solid::bond_atoms(&mut topology, &carbons)?;
        let free = diamond_solid::free_directions(&topology, &carbons);

        let mut plan = diamond_solid::CapPlan::new();
        for direction in &free {
            let Some(host_m) = topology.position_m(direction.host as usize) else {
                continue;
            };
            let probe_m = [
                host_m[0] + PROBE_STEP_M * direction.direction_m[0],
                host_m[1] + PROBE_STEP_M * direction.direction_m[1],
                host_m[2] + PROBE_STEP_M * direction.direction_m[2],
            ];
            if cavity.contains_m(probe_m) {
                plan.insert((direction.host, direction.slot), group);
            }
        }
        if plan.is_empty() {
            return Err(PartError::InvalidGeometry(
                "no free direction points into the well, so the well wall \
                 carries no group: make the well wider or deeper"
                    .to_string(),
            ));
        }

        diamond_solid::cap_free(&mut topology, &free, &plan)?;
        charge_wall(&mut topology, group);
        Ok(Part::new("binding_pocket", topology).with_material("diamond with functional groups"))
    }
}

/// Give the wall groups their stated model charges.
///
/// Each group is modelled as that group on a methyl carbon, so one charge set
/// serves every group of that kind. The group heavy atom and its hydrogens
/// take the charges of the model, and the lattice carbon that carries the
/// group takes the host charge. Every other atom keeps zero charge.
///
/// The function returns the number of atoms that it charged. The charges are a
/// model and they are not validated against experiment. See ADR-0063.
pub fn charge_wall(topology: &mut Topology, group: FunctionalGroup) -> usize {
    let host_charge_c = crate::group_data::GROUP_CHARGES[group.index()].host_charge_c;
    let heavy_charge_c = crate::group_data::GROUP_CHARGES[group.index()].heavy_charge_c;
    let hydrogen_charge_c = crate::group_data::GROUP_CHARGES[group.index()].hydrogen_charge_c;
    let heavy_type = group.heavy_type();

    let mut heavy_atoms = Vec::new();
    for index in 0..topology.atom_count() {
        if topology.atom_type(index) == Some(heavy_type) {
            heavy_atoms.push(index);
        }
    }

    let mut updates: Vec<(usize, f64)> = Vec::new();
    for &heavy in &heavy_atoms {
        updates.push((heavy, heavy_charge_c));
        for bond in topology.bonds() {
            let other = if bond.u as usize == heavy {
                bond.v as usize
            } else if bond.v as usize == heavy {
                bond.u as usize
            } else {
                continue;
            };
            match topology.element(other) {
                Some(Element::HYDROGEN) => updates.push((other, hydrogen_charge_c)),
                Some(Element::CARBON) if topology.atom_type(other) == Some("C") => {
                    updates.push((other, host_charge_c));
                }
                _ => {}
            }
        }
    }

    for (index, charge_c) in &updates {
        let _ = topology.set_charge_c(*index, *charge_c);
    }
    updates.len()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::guest::{guest, guests};
    use crate::port::Dof;

    fn origin_m() -> [f64; 3] {
        let thickness_m = BindingPocketGenerator.default_m("thickness_m");
        [0.0, 0.0, 0.5 * thickness_m]
    }

    fn radial_m(position_m: [f64; 3]) -> f64 {
        (position_m[0] * position_m[0] + position_m[1] * position_m[1]).sqrt()
    }

    #[test]
    fn the_pocket_has_a_cavity() {
        let part = BindingPocketGenerator
            .generate_with_defaults()
            .expect("valid");
        let radius_m = BindingPocketGenerator.default_m("pocket_radius_m");
        let center_m = BindingPocketGenerator.cavity_center_m();
        let depth_m = BindingPocketGenerator.default_m("pocket_depth_m");
        let inside = (0..part.atom_count())
            .filter(|index| part.topology.element(*index) == Some(Element::CARBON))
            .filter_map(|index| part.topology.position_m(index))
            .filter(|position_m| {
                (position_m[2] - center_m[2]).abs() <= 0.5 * depth_m
                    && radial_m(*position_m) < radius_m - 1.0e-11
            })
            .count();
        assert_eq!(inside, 0, "the well must hold no lattice site");
    }

    #[test]
    fn the_well_reaches_the_top_face() {
        let part = BindingPocketGenerator
            .generate_with_defaults()
            .expect("valid");
        let radius_m = BindingPocketGenerator.default_m("pocket_radius_m");
        let top_m = origin_m()[2];
        let open = (0..part.atom_count())
            .filter(|index| part.topology.element(*index) == Some(Element::CARBON))
            .filter_map(|index| part.topology.position_m(index))
            .filter(|position_m| {
                position_m[2] > top_m - 1.0e-10 && radial_m(*position_m) < radius_m
            })
            .count();
        assert_eq!(open, 0, "the well must open at the top face");
    }

    #[test]
    fn a_wider_well_leaves_a_wider_void() {
        let mut parameters = ParameterSet::new();
        parameters.set("pocket_radius_m", 1.6e-9);
        let wide = BindingPocketGenerator.generate(&parameters).expect("valid");
        let center_m = BindingPocketGenerator.cavity_center_m();
        let depth_m = BindingPocketGenerator.default_m("pocket_depth_m");
        let inside = (0..wide.atom_count())
            .filter(|index| wide.topology.element(*index) == Some(Element::CARBON))
            .filter_map(|index| wide.topology.position_m(index))
            .filter(|position_m| {
                (position_m[2] - center_m[2]).abs() <= 0.5 * depth_m
                    && radial_m(*position_m) < 1.6e-9 - 1.0e-11
            })
            .count();
        assert_eq!(inside, 0, "the wide well must hold no lattice site");
        assert!(
            wide.atom_count()
                > BindingPocketGenerator
                    .generate_with_defaults()
                    .expect("valid")
                    .atom_count(),
            "a wider well needs a wider wall, so it holds more atoms"
        );
    }

    #[test]
    fn the_wall_carries_the_chosen_group() {
        let mut parameters = ParameterSet::new();
        parameters.set("wall_group", FunctionalGroup::Fluoro.index() as f64);
        let part = BindingPocketGenerator.generate(&parameters).expect("valid");
        let fluorine = (0..part.atom_count())
            .filter(|index| part.topology.element(*index) == Some(Element::FLUORINE))
            .count();
        assert!(fluorine > 0, "the wall must carry the chosen group");
        assert_eq!(
            (0..part.atom_count())
                .filter(|index| part.topology.element(*index) == Some(Element::OXYGEN))
                .count(),
            0,
            "a fluorine wall holds no oxygen"
        );
    }

    #[test]
    fn an_unknown_wall_group_is_refused() {
        let mut parameters = ParameterSet::new();
        parameters.set("wall_group", 9.0);
        assert!(BindingPocketGenerator.generate(&parameters).is_err());
    }

    #[test]
    fn a_well_through_the_floor_is_refused() {
        let mut parameters = ParameterSet::new();
        parameters.set("pocket_depth_m", 4.0e-9);
        parameters.set("thickness_m", 2.0e-9);
        assert!(BindingPocketGenerator.generate(&parameters).is_err());
    }

    #[test]
    fn the_opening_port_points_out_of_the_pocket() {
        let port = BindingPocketGenerator.opening_port();
        assert_eq!(port.name, "opening");
        assert_eq!(port.axis_m, [0.0, 0.0, 1.0]);
        assert_eq!(port.dof, Dof::Fixed);
        assert!((port.origin_m[2] - origin_m()[2]).abs() < 1.0e-12);
    }

    #[test]
    fn every_group_index_round_trips() {
        for group in crate::group::FUNCTIONAL_GROUPS {
            assert_eq!(FunctionalGroup::from_index(group.index()), Some(group));
        }
        assert_eq!(FunctionalGroup::from_index(6), None);
    }

    #[test]
    fn the_wall_contact_distance_is_the_carbon_distance() {
        let contact_m = wall_contact_distance_m().expect("carbon is in the table");
        assert!((contact_m - 3.851e-10).abs() < 1.0e-13);
    }

    #[test]
    fn a_small_guest_fits_and_a_large_guest_does_not() {
        let methanol = guest("methanol").expect("known guest");
        let cyclohexane = guest("cyclohexane").expect("known guest");
        assert!(BindingPocketGenerator.guest_clearance_m(methanol) > 0.0);
        assert!(
            BindingPocketGenerator.guest_clearance_m(cyclohexane)
                < BindingPocketGenerator.guest_clearance_m(methanol),
            "a larger guest has less clearance"
        );
    }

    #[test]
    fn the_isomer_pair_needs_different_wells() {
        let ethanol = guest("ethanol").expect("known guest");
        let ether = guest("dimethyl_ether").expect("known guest");
        let clearance_m = 1.0e-10;
        let for_ethanol = pocket_radius_for(ethanol, clearance_m);
        let for_ether = pocket_radius_for(ether, clearance_m);
        assert!(
            (for_ethanol - for_ether).abs() > 1.0e-13,
            "the two isomers have different extents, so one well radius cannot \
             suit both"
        );
    }

    #[test]
    fn every_guest_states_a_wall_well_depth() {
        for guest in guests() {
            assert!(wall_well_depth_j(guest).is_some());
        }
    }

    #[test]
    fn a_guest_scale_well_builds_and_decorates() {
        let mut parameters = ParameterSet::new();
        parameters.set("pocket_radius_m", 2.3e-10);
        let part = BindingPocketGenerator.generate(&parameters).expect("valid");
        let oxygen = (0..part.atom_count())
            .filter(|index| part.topology.element(*index) == Some(Element::OXYGEN))
            .count();
        assert!(oxygen > 0, "a guest-scale well still carries its groups");
        assert!(part.atom_count() < 3000);
    }

    #[test]
    fn a_guest_scale_well_separates_the_size_classes() {
        let radius_m = 2.3e-10;
        let methanol = guest("methanol").expect("known guest");
        let ethanol = guest("ethanol").expect("known guest");
        let ether = guest("dimethyl_ether").expect("known guest");
        let benzene = guest("benzene").expect("known guest");
        let cyclohexane = guest("cyclohexane").expect("known guest");
        for small in [methanol, ethanol, ether] {
            assert!(small.extent_m() < radius_m, "{} fits", small.id);
        }
        for large in [benzene, cyclohexane] {
            assert!(large.extent_m() > radius_m, "{} does not fit", large.id);
        }
        assert_eq!(guests().len(), 5, "the five guests cover both size classes");
    }

    #[test]
    fn the_well_radius_resolves_in_lattice_steps() {
        let mut narrow = ParameterSet::new();
        narrow.set("pocket_radius_m", 2.2e-10);
        let mut wide = ParameterSet::new();
        wide.set("pocket_radius_m", 2.3e-10);
        let narrow = BindingPocketGenerator.generate(&narrow).expect("valid");
        let wide = BindingPocketGenerator.generate(&wide).expect("valid");
        assert_eq!(
            narrow.atom_count(),
            wide.atom_count(),
            "the diamond lattice is discrete, so a tenth of an Angstrom can              change nothing"
        );
    }

    #[test]
    fn every_carbon_keeps_a_valid_valence() {
        let part = BindingPocketGenerator
            .generate_with_defaults()
            .expect("valid");
        let mut degree = vec![0_u32; part.atom_count()];
        for bond in part.topology.bonds() {
            degree[bond.u as usize] += 1;
            degree[bond.v as usize] += 1;
        }
        for (index, count) in degree.iter().enumerate() {
            let element = part
                .topology
                .element(index)
                .expect("every atom has an element");
            let want = nanocad_model::valence(element).expect("the table holds every element here");
            assert_eq!(
                *count,
                want as u32,
                "atom {index} of element {} has degree {}",
                element.symbol().unwrap_or("?"),
                count
            );
        }
    }
}
