//! The sorting rotor disk.
//!
//! The sorting rotor is the device in Freitas, *Nanomedicine* Volume I,
//! Section 3.4.2, after Drexler. A disk carries a row of binding pockets along
//! its rim. The disk turns about its axis, and each pocket carries its bound
//! molecule from the outside solution to an inner chamber.
//!
//! This generator builds the disk and its pockets. It does not build a binding
//! site, and it does not model a molecule. A pocket is a cylindrical void that
//! opens to the rim.
//!
//! The disk is cut from the diamond cubic lattice, so every bond is the diamond
//! first-shell length and every carbon keeps a valence of four. All lengths are
//! SI metres and all angles are radians.

use std::f64::consts::PI;

use nanocad_model::{Atom, Element, Part, Topology};
use nanocad_units::Unit;

use crate::diamond_solid;
use crate::error::PartError;
use crate::generator::PartGenerator;
use crate::group::FunctionalGroup;
use crate::lattice_fill::fill_solid;
use crate::parameter::{ParameterSet, ParameterSpec};
use crate::pocket::charge_wall;
use crate::port::{Dof, Port};
use crate::shape::{Box3, Cylinder, Difference, RadialCylinder, Solid, Union};

/// The longest bond in diamond, in metres. A probe shorter than this stays
/// inside the cell that owns the direction.
const PROBE_STEP_M: f64 = 0.6 * 1.544e-10;

/// The parameter specs of [`SortingRotorGenerator`], in a stable order.
static SORTING_ROTOR_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new(
        "disc_radius_m",
        Some(Unit::Metre),
        7.0e-9,
        2.0e-9,
        4.0e-8,
        false,
        "outside radius of the rotor disk in metres",
    ),
    ParameterSpec::new(
        "thickness_m",
        Some(Unit::Metre),
        2.0e-9,
        7.0e-10,
        2.0e-8,
        false,
        "axial thickness of the disk in metres",
    ),
    ParameterSpec::new(
        "pocket_count",
        None,
        12.0,
        3.0,
        24.0,
        true,
        "number of binding pockets on the rim",
    ),
    ParameterSpec::new(
        "pocket_radius_m",
        Some(Unit::Metre),
        1.0e-9,
        5.0e-10,
        5.0e-9,
        false,
        "radius of each pocket void in metres",
    ),
    ParameterSpec::new(
        "pocket_orbit_m",
        Some(Unit::Metre),
        6.5e-9,
        1.0e-9,
        4.0e-8,
        false,
        "radius of the pocket centre circle in metres",
    ),
    ParameterSpec::new(
        "bore_radius_m",
        Some(Unit::Metre),
        1.5e-9,
        0.0,
        2.0e-8,
        false,
        "radius of the central bore in metres, which holds the drive shaft",
    ),
    ParameterSpec::new(
        "ejection_bore_radius_m",
        Some(Unit::Metre),
        6.0e-10,
        0.0,
        5.0e-9,
        false,
        "radius of each radial ejection bore in metres, or zero for no bore",
    ),
    ParameterSpec::new(
        "wall_group",
        None,
        0.0,
        0.0,
        5.0,
        true,
        "functional group on the pocket walls, as an index into the group list",
    ),
    ParameterSpec::new(
        "keyway_width_m",
        Some(Unit::Metre),
        3.0e-10,
        0.0,
        2.0e-9,
        false,
        "width of the keyway in the central bore in metres, or zero for no keyway",
    ),
    ParameterSpec::new(
        "keyway_depth_m",
        Some(Unit::Metre),
        3.0e-10,
        1.0e-11,
        1.0e-9,
        false,
        "radial depth of the keyway past the bore wall in metres",
    ),
];

/// Generates a sorting rotor disk.
///
/// The disk is the solid between the central bore and the outside radius, less
/// one cylindrical pocket for every slot. A pocket whose centre circle plus its
/// radius reaches past the outside radius opens to the rim, which is the
/// binding face the solution sees.
///
/// Each pocket also has a radial ejection bore, unless the bore radius is zero.
/// The bore runs from the central bore out to the inner wall of the pocket, so a
/// rod inside it can thrust a bound guest out of the pocket.
///
/// The disk turns about `z`, so the axis port lies on that axis.
#[derive(Clone, Copy, Debug, Default)]
pub struct SortingRotorGenerator;

impl SortingRotorGenerator {
    /// Returns the axis port.
    ///
    /// The port sits at the origin, its axis is `+z`, and it allows one
    /// rotation about that axis. The gap is zero.
    pub fn axis_port(&self) -> Port {
        let mut port = Port::named("axis");
        port.dof = Dof::Revolute;
        port
    }

    /// Returns the port of one pocket.
    ///
    /// The port sits at the pocket centre, its axis is `+z`, and it holds the
    /// pocket rigidly. The angle of pocket `index` is `2 pi index / count`.
    pub fn pocket_port(&self, index: usize, count: usize, pocket_orbit_m: f64) -> Port {
        let angle_rad = 2.0 * PI * index as f64 / count.max(1) as f64;
        let (sin_rad, cos_rad) = angle_rad.sin_cos();
        let mut port = Port::named(&format!("pocket_{index}"));
        port.origin_m = [pocket_orbit_m * cos_rad, pocket_orbit_m * sin_rad, 0.0];
        port
    }

    /// Returns the functional group that the pocket walls carry.
    ///
    /// The parameter holds an index into
    /// [`crate::group::FUNCTIONAL_GROUPS`], because a parameter set holds
    /// numbers. An index outside the list falls back to the first group here,
    /// because [`PartGenerator::generate`] refuses it.
    pub fn wall_group(&self) -> FunctionalGroup {
        let index = self.default_m("wall_group").round().max(0.0) as usize;
        FunctionalGroup::from_index(index).unwrap_or(FunctionalGroup::Hydroxyl)
    }

    /// Returns the centre of one pocket, in metres.
    pub fn pocket_center_m(index: usize, count: usize, pocket_orbit_m: f64) -> [f64; 3] {
        let angle_rad = 2.0 * PI * index as f64 / count.max(1) as f64;
        let (sin_rad, cos_rad) = angle_rad.sin_cos();
        [pocket_orbit_m * cos_rad, pocket_orbit_m * sin_rad, 0.0]
    }

    /// Returns the radial ejection bore of one pocket, in the rotor frame.
    ///
    /// The bore runs from the central bore out to the inner wall of the pocket,
    /// so a rod inside it can thrust a guest out of the pocket. Each end is
    /// extended by the bore radius, so the bore opens into both the central
    /// bore and the pocket at every lattice layer.
    pub fn ejection_bore(
        index: usize,
        count: usize,
        bore_radius_m: f64,
        pocket_orbit_m: f64,
        pocket_radius_m: f64,
        ejection_bore_radius_m: f64,
    ) -> RadialCylinder {
        let angle_rad = 2.0 * PI * index as f64 / count.max(1) as f64;
        let inner_radius_m = bore_radius_m;
        let outer_radius_m = pocket_orbit_m - pocket_radius_m;
        RadialCylinder {
            azimuth_rad: angle_rad,
            center_radius_m: 0.5 * (inner_radius_m + outer_radius_m),
            radius_m: ejection_bore_radius_m,
            half_length_m: 0.5 * (outer_radius_m - inner_radius_m) + ejection_bore_radius_m,
            center_z_m: 0.0,
        }
    }
}

impl SortingRotorGenerator {
    /// Resolves one declared default, so a helper can state a geometry.
    fn default_m(&self, name: &str) -> f64 {
        self.resolve(&ParameterSet::new())
            .ok()
            .and_then(|resolved| resolved.get(name))
            .unwrap_or(0.0)
    }
}

impl PartGenerator for SortingRotorGenerator {
    fn id(&self) -> &'static str {
        "sorting_rotor"
    }

    fn name(&self) -> &'static str {
        "sorting rotor"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        SORTING_ROTOR_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        let resolved = self.resolve(parameters)?;
        let disc_radius_m = resolved.require("disc_radius_m")?;
        let thickness_m = resolved.require("thickness_m")?;
        let pocket_count = resolved.require("pocket_count")? as usize;
        let pocket_radius_m = resolved.require("pocket_radius_m")?;
        let pocket_orbit_m = resolved.require("pocket_orbit_m")?;
        let bore_radius_m = resolved.require("bore_radius_m")?;
        let ejection_bore_radius_m = resolved.require("ejection_bore_radius_m")?;
        let keyway_width_m = resolved.require("keyway_width_m")?;
        let keyway_depth_m = resolved.require("keyway_depth_m")?;
        let wall_index = resolved.require("wall_group")?.round();
        let Some(group) = FunctionalGroup::from_index(wall_index.max(0.0) as usize) else {
            return Err(PartError::InvalidGeometry(format!(
                "the wall group index {wall_index} names no functional group"
            )));
        };

        if pocket_orbit_m < bore_radius_m + pocket_radius_m {
            return Err(PartError::InvalidGeometry(format!(
                "the pocket orbit {pocket_orbit_m} m is inside the bore {bore_radius_m} m \
                 plus the pocket radius {pocket_radius_m} m, so a pocket would breach the bore"
            )));
        }
        if pocket_orbit_m - pocket_radius_m >= disc_radius_m {
            return Err(PartError::InvalidGeometry(format!(
                "the pocket orbit {pocket_orbit_m} m minus the pocket radius {pocket_radius_m} m \
                 is at or outside the disc radius {disc_radius_m} m, so the pockets would not bite"
            )));
        }
        if ejection_bore_radius_m > 0.0 {
            if ejection_bore_radius_m >= pocket_radius_m {
                return Err(PartError::InvalidGeometry(format!(
                    "the ejection bore radius {ejection_bore_radius_m} m is not smaller than the \
                     pocket radius {pocket_radius_m} m, so a bore would eat the whole pocket"
                )));
            }
            if pocket_orbit_m - pocket_radius_m <= bore_radius_m {
                return Err(PartError::InvalidGeometry(format!(
                    "the inner pocket wall {} m is not outside the central bore {bore_radius_m} m, \
                     so an ejection bore has no length",
                    pocket_orbit_m - pocket_radius_m
                )));
            }
        }

        let center_m = [0.0, 0.0, 0.0];
        let disc = Cylinder {
            center_m,
            radius_m: disc_radius_m,
            height_m: thickness_m,
        };
        let mut void: Box<dyn Solid> = Box::new(Cylinder {
            center_m,
            radius_m: bore_radius_m,
            height_m: thickness_m,
        });
        if keyway_width_m > 0.0 {
            if keyway_depth_m <= 0.0 {
                return Err(PartError::InvalidGeometry(
                    "the keyway has a width but no depth".to_string(),
                ));
            }
            if bore_radius_m + keyway_depth_m >= disc_radius_m {
                return Err(PartError::InvalidGeometry(format!(
                    "the keyway reaches {} m and leaves no disk at the bore",
                    bore_radius_m + keyway_depth_m
                )));
            }
            let keyway = Box3 {
                min_m: [
                    bore_radius_m - keyway_depth_m,
                    -0.5 * keyway_width_m,
                    -0.5 * thickness_m,
                ],
                max_m: [
                    bore_radius_m + keyway_depth_m,
                    0.5 * keyway_width_m,
                    0.5 * thickness_m,
                ],
            };
            void = Box::new(Union {
                a: void,
                b: Box::new(keyway),
            });
        }
        let mut pockets = Vec::with_capacity(pocket_count);
        for index in 0..pocket_count {
            let pocket = Cylinder {
                center_m: Self::pocket_center_m(index, pocket_count, pocket_orbit_m),
                radius_m: pocket_radius_m,
                height_m: thickness_m,
            };
            void = Box::new(Union {
                a: void,
                b: Box::new(pocket),
            });
            pockets.push(pocket);
        }
        if ejection_bore_radius_m > 0.0 {
            for index in 0..pocket_count {
                let bore = Self::ejection_bore(
                    index,
                    pocket_count,
                    bore_radius_m,
                    pocket_orbit_m,
                    pocket_radius_m,
                    ejection_bore_radius_m,
                );
                void = Box::new(Union {
                    a: void,
                    b: Box::new(bore),
                });
            }
        }

        let solid = Difference {
            outer: Box::new(disc),
            inner: void,
        };
        fill_part(
            &format!("sorting-rotor-{pocket_count}"),
            &solid,
            &pockets,
            group,
            "the disc is thinner than one lattice layer",
        )
    }
}

/// Fills a solid region with the diamond lattice, decorates the pocket walls,
/// and caps every remaining free direction with hydrogen.
///
/// A free tetrahedral direction takes `group` when a short step along it lands
/// inside one of `pockets`, because that direction faces the space that a guest
/// occupies. Every other free direction takes a hydrogen.
///
/// The region must have a finite bounding box. The function returns an error
/// when the region cut no lattice site, or when no free direction faces a
/// pocket, so an undecorated rotor cannot reach a caller.
fn fill_part(
    name: &str,
    solid: &dyn Solid,
    pockets: &[Cylinder],
    group: FunctionalGroup,
    empty_note: &str,
) -> Result<Part, PartError> {
    let sites_m = fill_solid(solid);
    let mut topology = Topology::new();
    let mut carbons = Vec::with_capacity(sites_m.len());
    for site_m in sites_m {
        carbons.push(topology.add_atom(Atom::new(Element::CARBON, site_m, 0.0, "C")));
    }
    if carbons.is_empty() {
        return Err(PartError::InvalidGeometry(format!(
            "{name} cut no atom from the diamond lattice: {empty_note}"
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
        if pockets.iter().any(|pocket| pocket.contains_m(probe_m)) {
            plan.insert((direction.host, direction.slot), group);
        }
    }
    if plan.is_empty() {
        return Err(PartError::InvalidGeometry(
            "no free direction points into a pocket, so the pocket walls carry \
             no group: make the pockets wider or the orbit smaller"
                .to_string(),
        ));
    }

    diamond_solid::cap_free(&mut topology, &free, &plan)?;
    charge_wall(&mut topology, group);
    Ok(Part::new(name, topology).with_material("diamond with functional groups"))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Counts the atoms whose force-field type starts with the group's prefix.
    fn group_atoms(part: &Part, group: FunctionalGroup) -> usize {
        let prefix = group.heavy_type().trim_end_matches('_');
        (0..part.atom_count())
            .filter(|&index| {
                part.topology
                    .atom_type(index)
                    .is_some_and(|name| name.starts_with(prefix))
            })
            .count()
    }

    /// Generates the rotor with one parameter changed.
    fn rotor_with(name: &str, value: f64) -> Part {
        let mut parameters = ParameterSet::new();
        parameters.set(name, value);
        SortingRotorGenerator
            .generate(&parameters)
            .expect("the rotor builds")
    }

    fn radial_xy_m(position_m: [f64; 3]) -> f64 {
        (position_m[0] * position_m[0] + position_m[1] * position_m[1]).sqrt()
    }

    #[test]
    fn the_rotor_has_a_void_for_every_pocket() {
        let pocket_count = 12;
        let pocket_radius_m = 1.0e-9;
        let pocket_orbit_m = 6.5e-9;
        let part = SortingRotorGenerator
            .generate_with_defaults()
            .expect("rotor");
        for index in 0..pocket_count {
            let center_m =
                SortingRotorGenerator::pocket_center_m(index, pocket_count, pocket_orbit_m);
            let half_m = 0.5 * pocket_radius_m;
            let inside = part.topology.atoms().filter(|atom| {
                let dx = atom.position_m[0] - center_m[0];
                let dy = atom.position_m[1] - center_m[1];
                (dx * dx + dy * dy).sqrt() < half_m
            });
            assert_eq!(
                inside.count(),
                0,
                "pocket {index} holds a carbon at its centre"
            );
        }
    }

    #[test]
    fn the_rotor_respects_its_radius() {
        let disc_radius_m = 7.0e-9;
        let part = SortingRotorGenerator
            .generate_with_defaults()
            .expect("rotor");
        let farthest_m = part
            .topology
            .atoms()
            .map(|atom| radial_xy_m(atom.position_m))
            .fold(0.0_f64, f64::max);
        assert!(
            farthest_m <= disc_radius_m + 1.0e-10,
            "an atom sits at {farthest_m:e} m, past the disc radius {disc_radius_m:e} m"
        );
    }

    #[test]
    fn the_bore_is_open() {
        let bore_radius_m = 1.5e-9;
        let part = SortingRotorGenerator
            .generate_with_defaults()
            .expect("rotor");
        let inside = part
            .topology
            .atoms()
            .filter(|atom| radial_xy_m(atom.position_m) < 0.5 * bore_radius_m);
        assert_eq!(inside.count(), 0, "the bore holds a carbon");
    }

    #[test]
    fn every_ejection_bore_is_open() {
        let pocket_count = 12;
        let part = SortingRotorGenerator
            .generate_with_defaults()
            .expect("rotor");
        let mut occupied = 0;
        for index in 0..pocket_count {
            let bore = SortingRotorGenerator::ejection_bore(
                index,
                pocket_count,
                1.5e-9,
                6.5e-9,
                1.0e-9,
                5.0e-10,
            );
            occupied += part
                .topology
                .atoms()
                .filter(|atom| atom.element == Element::CARBON)
                .filter(|atom| bore.contains_m(atom.position_m))
                .count();
        }
        assert_eq!(occupied, 0, "an ejection bore holds a carbon");
    }

    #[test]
    fn a_zero_ejection_bore_removes_nothing() {
        let bored = SortingRotorGenerator
            .generate_with_defaults()
            .expect("rotor");
        let solid = rotor_with("ejection_bore_radius_m", 0.0);
        assert!(
            bored.atom_count() < solid.atom_count(),
            "the bores removed no atom: {} against {}",
            bored.atom_count(),
            solid.atom_count()
        );
    }

    #[test]
    fn the_keyway_is_void_and_removes_material() {
        let keyed = SortingRotorGenerator
            .generate_with_defaults()
            .expect("rotor");
        let plain = rotor_with("keyway_width_m", 0.0);
        assert!(
            keyed.atom_count() < plain.atom_count(),
            "the keyway removed no atom: {} against {}",
            keyed.atom_count(),
            plain.atom_count()
        );
        let bore_radius_m = 1.5e-9;
        let keyway_depth_m = 3.0e-10;
        let keyway_width_m = 3.0e-10;
        let occupied = keyed
            .topology
            .atoms()
            .filter(|atom| atom.element == Element::CARBON)
            .filter(|atom| {
                atom.position_m[0] > bore_radius_m - 0.5 * keyway_depth_m
                    && atom.position_m[0] < bore_radius_m + keyway_depth_m
                    && atom.position_m[1].abs() < 0.25 * keyway_width_m
            })
            .count();
        assert_eq!(occupied, 0, "the keyway holds a carbon");
    }

    #[test]
    fn more_pockets_change_the_rotor() {
        let parameters = ParameterSet::new()
            .with("pocket_count", 6.0)
            .with("pocket_radius_m", 1.3e-9)
            .with("pocket_orbit_m", 6.2e-9);
        let six = SortingRotorGenerator.generate(&parameters).expect("six");
        let parameters = ParameterSet::new()
            .with("pocket_count", 12.0)
            .with("pocket_radius_m", 1.3e-9)
            .with("pocket_orbit_m", 6.2e-9);
        let twelve = SortingRotorGenerator.generate(&parameters).expect("twelve");
        assert_ne!(six.atom_count(), twelve.atom_count());
    }

    #[test]
    fn a_pocket_that_breaches_the_bore_is_refused() {
        let parameters = ParameterSet::new()
            .with("pocket_orbit_m", 2.0e-9)
            .with("pocket_radius_m", 1.0e-9)
            .with("bore_radius_m", 1.5e-9);
        assert!(SortingRotorGenerator.generate(&parameters).is_err());
    }

    #[test]
    fn the_axis_port_is_revolute() {
        let port = SortingRotorGenerator.axis_port();
        assert_eq!(port.name, "axis");
        assert_eq!(port.dof, Dof::Revolute);
        assert_eq!(port.origin_m, [0.0, 0.0, 0.0]);
        assert_eq!(port.axis_m, [0.0, 0.0, 1.0]);
        assert_eq!(port.gap_m, 0.0);
    }

    #[test]
    fn a_pocket_port_sits_on_the_pocket_circle() {
        let port = SortingRotorGenerator.pocket_port(3, 12, 6.5e-9);
        assert_eq!(port.name, "pocket_3");
        assert_eq!(port.dof, Dof::Fixed);
        assert!((radial_xy_m(port.origin_m) - 6.5e-9).abs() < 1.0e-18);
        let quarter_turn = SortingRotorGenerator.pocket_port(3, 12, 6.5e-9);
        assert!(quarter_turn.origin_m[0].abs() < 1.0e-15);
        assert!((quarter_turn.origin_m[1] - 6.5e-9).abs() < 1.0e-18);
    }

    #[test]
    fn every_carbon_keeps_a_valid_valence() {
        let part = SortingRotorGenerator
            .generate_with_defaults()
            .expect("rotor");
        let mut degree = vec![0usize; part.atom_count()];
        for bond in part.topology.bonds() {
            degree[bond.u as usize] += 1;
            degree[bond.v as usize] += 1;
        }
        for (index, &count) in degree.iter().enumerate() {
            if part.topology.element(index) == Some(Element::CARBON) {
                assert!(
                    (1..=4).contains(&count),
                    "carbon {index} keeps a degree of {count}, outside 1 to 4"
                );
            }
        }
    }

    use crate::FUNCTIONAL_GROUPS;

    #[test]
    fn the_pocket_walls_carry_the_chosen_group() {
        let rotor = SortingRotorGenerator.generate_with_defaults().unwrap();
        assert_eq!(
            SortingRotorGenerator.wall_group(),
            FunctionalGroup::Hydroxyl
        );
        assert!(group_atoms(&rotor, FunctionalGroup::Hydroxyl) > 0);
        assert!(
            group_atoms(&rotor, FunctionalGroup::Hydroxyl) < rotor.atom_count() / 4,
            "only the pocket walls are decorated"
        );
    }

    #[test]
    fn a_different_wall_group_changes_the_wall() {
        let hydroxyl = rotor_with("wall_group", FunctionalGroup::Hydroxyl.index() as f64);
        let fluoro = rotor_with("wall_group", FunctionalGroup::Fluoro.index() as f64);
        assert!(group_atoms(&fluoro, FunctionalGroup::Fluoro) > 0);
        assert_eq!(group_atoms(&fluoro, FunctionalGroup::Hydroxyl), 0);
        assert_ne!(hydroxyl.atom_count(), fluoro.atom_count());
    }

    #[test]
    fn every_wall_group_builds_a_rotor() {
        for group in FUNCTIONAL_GROUPS {
            let rotor = rotor_with("wall_group", group.index() as f64);
            assert!(
                group_atoms(&rotor, group) > 0,
                "{} decorates",
                group.label()
            );
        }
    }

    #[test]
    fn an_unknown_wall_group_is_refused() {
        let mut parameters = ParameterSet::new();
        parameters.set("wall_group", FUNCTIONAL_GROUPS.len() as f64);
        let error = SortingRotorGenerator
            .generate(&parameters)
            .expect_err("the index names no group");
        assert!(
            matches!(
                error,
                PartError::InvalidGeometry(_) | PartError::ParameterOutOfRange { .. }
            ),
            "the refusal names the group index"
        );
    }

    #[test]
    fn the_pocket_wall_carries_the_model_charges() {
        let rotor = SortingRotorGenerator.generate_with_defaults().unwrap();
        let mut charged = 0;
        let mut neutral = 0;
        for index in 0..rotor.atom_count() {
            let Some(charge_c) = rotor.topology.charge_c(index) else {
                continue;
            };
            if charge_c == 0.0 {
                neutral += 1;
            } else {
                charged += 1;
            }
        }
        assert!(charged > 0, "the wall groups carry a charge");
        assert!(charged < neutral, "the lattice carbons stay neutral");
    }
}
