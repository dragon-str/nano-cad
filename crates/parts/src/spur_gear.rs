use std::f64::consts::PI;

use nanocad_model::{Atom, Bond, BondType, Element, Part, Topology};
use nanocad_units::Unit;

use crate::error::PartError;
use crate::gear_profile::GearProfile;
use crate::generator::PartGenerator;
use crate::geometry::resample_closed_polyline;
use crate::parameter::{ParameterSet, ParameterSpec};

/// The diamond carbon-carbon bond length, in metres.
///
/// Source: the diamond cubic lattice constant `a = 3.567e-10 m` and the
/// first-shell distance `a sqrt(3) / 4 = 1.544e-10 m`.
const DIAMONDOID_BOND_M: f64 = 1.544e-10;

/// The default pressure angle, 20 degrees, in radians.
const DEFAULT_PRESSURE_ANGLE_RAD: f64 = 20.0 * PI / 180.0;

static SPUR_GEAR_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new(
        "module_m",
        Some(Unit::Metre),
        0.5e-9,
        0.05e-9,
        5.0e-9,
        false,
        "gear module in metres",
    ),
    ParameterSpec::new(
        "tooth_count",
        None,
        20.0,
        3.0,
        120.0,
        true,
        "number of teeth",
    ),
    ParameterSpec::new(
        "pressure_angle_rad",
        None,
        DEFAULT_PRESSURE_ANGLE_RAD,
        10.0 * PI / 180.0,
        25.0 * PI / 180.0,
        false,
        "pressure angle in radians",
    ),
    ParameterSpec::new(
        "face_width_m",
        Some(Unit::Metre),
        0.6e-9,
        0.1e-9,
        5.0e-9,
        false,
        "tooth face width in metres",
    ),
    ParameterSpec::new(
        "bore_radius_m",
        Some(Unit::Metre),
        0.5e-9,
        0.05e-9,
        5.0e-9,
        false,
        "central bore radius in metres",
    ),
    ParameterSpec::new(
        "hub_radius_m",
        Some(Unit::Metre),
        1.2e-9,
        0.05e-9,
        5.0e-9,
        false,
        "hub outside radius in metres",
    ),
    ParameterSpec::new(
        "hub_width_m",
        Some(Unit::Metre),
        1.0e-9,
        0.1e-9,
        5.0e-9,
        false,
        "hub axial width in metres",
    ),
];

/// Builds an atomistic spur gear: teeth, face width, bore, and hub.
///
/// The gear is a diamondoid analog. A rim of carbon atoms follows the involute
/// tooth outline at the diamond bond spacing. Radial spokes join the rim to a
/// cylindrical hub that surrounds the empty bore. The part is skeletal, in the
/// style of the NanoEngineer gear moieties, not a filled solid.
#[derive(Clone, Copy, Debug, Default)]
pub struct SpurGearGenerator;

impl PartGenerator for SpurGearGenerator {
    fn id(&self) -> &'static str {
        "spur_gear"
    }

    fn name(&self) -> &'static str {
        "Atomistic spur gear"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        SPUR_GEAR_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        let resolved = self.resolve(parameters)?;
        let module_m = resolved.require("module_m")?;
        let teeth = resolved.require("tooth_count")? as usize;
        let pressure_angle_rad = resolved.require("pressure_angle_rad")?;
        let face_width_m = resolved.require("face_width_m")?;
        let bore_radius_m = resolved.require("bore_radius_m")?;
        let hub_radius_m = resolved.require("hub_radius_m")?;
        let hub_width_m = resolved.require("hub_width_m")?;

        let profile = GearProfile::new(module_m, teeth, pressure_angle_rad)?;
        profile.check_undercut()?;
        let root_m = profile.root_radius_m();
        let outer_m = profile.outer_radius_m();
        if bore_radius_m >= root_m {
            return Err(PartError::BoreTooLarge {
                bore_m: bore_radius_m,
                root_m,
            });
        }
        if hub_radius_m <= bore_radius_m || hub_radius_m >= root_m {
            return Err(PartError::InvalidGeometry(format!(
                "hub radius {hub_radius_m} m must be between the bore radius \
                 {bore_radius_m} m and the root radius {root_m} m"
            )));
        }

        let bond_m = DIAMONDOID_BOND_M;
        let pitch_angle_rad = 2.0 * PI / teeth as f64;

        let outline = profile.outline_points(6, 3)?;
        let mut rim_points = resample_closed_polyline(&outline, bond_m);
        for tooth in 0..teeth {
            let angle_rad = tooth as f64 * pitch_angle_rad;
            snap_nearest_angle(&mut rim_points, angle_rad, outer_m);
        }

        let body_layers = layer_positions(face_width_m, bond_m);
        let hub_layers = layer_positions(hub_width_m, bond_m);

        let mut topology = Topology::new();
        let mut rim_index: Vec<Vec<u32>> = Vec::with_capacity(body_layers.len());
        for z_m in &body_layers {
            let mut layer = Vec::with_capacity(rim_points.len());
            for point in &rim_points {
                layer.push(add_carbon(&mut topology, point[0], point[1], *z_m));
            }
            rim_index.push(layer);
        }

        let hub_samples = (((2.0 * PI * hub_radius_m / bond_m).round()) as usize).max(8);
        let mut hub_index: Vec<Vec<u32>> = Vec::with_capacity(hub_layers.len());
        for z_m in &hub_layers {
            let mut layer = Vec::with_capacity(hub_samples);
            for sample in 0..hub_samples {
                let angle_rad = 2.0 * PI * sample as f64 / hub_samples as f64;
                layer.push(add_carbon(
                    &mut topology,
                    hub_radius_m * angle_rad.cos(),
                    hub_radius_m * angle_rad.sin(),
                    *z_m,
                ));
            }
            hub_index.push(layer);
        }

        let spoke_count = (teeth / 2).clamp(4, 12);
        let spoke_steps = (((root_m - hub_radius_m) / bond_m).round() as usize).max(1);
        let mut spoke_index: Vec<Vec<Vec<u32>>> = Vec::with_capacity(body_layers.len());
        for z_m in &body_layers {
            let mut layer = Vec::with_capacity(spoke_count);
            for spoke in 0..spoke_count {
                let angle_rad = 2.0 * PI * spoke as f64 / spoke_count as f64;
                let mut chain = Vec::with_capacity(spoke_steps + 1);
                for step in 0..=spoke_steps {
                    let t = step as f64 / spoke_steps as f64;
                    let radius_m = hub_radius_m + (root_m - hub_radius_m) * t;
                    chain.push(add_carbon(
                        &mut topology,
                        radius_m * angle_rad.cos(),
                        radius_m * angle_rad.sin(),
                        *z_m,
                    ));
                }
                layer.push(chain);
            }
            spoke_index.push(layer);
        }

        let mut bonds: Vec<(u32, u32)> = Vec::new();
        for (layer_number, layer) in rim_index.iter().enumerate() {
            bonds.extend(ring_bonds(layer));
            if layer_number > 0 {
                bonds.extend(column_bonds(&rim_index[layer_number - 1], layer));
            }
        }
        for (layer_number, layer) in hub_index.iter().enumerate() {
            bonds.extend(ring_bonds(layer));
            if layer_number > 0 {
                bonds.extend(column_bonds(&hub_index[layer_number - 1], layer));
            }
        }
        for (layer_number, layer) in spoke_index.iter().enumerate() {
            for chain in layer {
                bonds.extend(chain.windows(2).map(|pair| (pair[0], pair[1])));
            }
            if layer_number > 0 {
                let previous = &spoke_index[layer_number - 1];
                for (spoke, chain) in layer.iter().enumerate() {
                    for (step, atom) in chain.iter().enumerate() {
                        bonds.push((previous[spoke][step], *atom));
                    }
                }
            }
        }
        for (u, v) in bonds {
            if u != v {
                topology.add_bond(Bond::new(u, v, 1, BondType::Single))?;
            }
        }

        let name = format!("spur-gear-z{teeth}");
        Ok(Part::new(name, topology).with_material("diamondoid"))
    }
}

fn add_carbon(topology: &mut Topology, x_m: f64, y_m: f64, z_m: f64) -> u32 {
    topology.add_atom(Atom::new(Element::CARBON, [x_m, y_m, z_m], 0.0, "C"))
}

/// Returns evenly spaced axial positions across a total width.
fn layer_positions(total_m: f64, spacing_m: f64) -> Vec<f64> {
    let count = ((total_m / spacing_m).round() as usize).max(1) + 1;
    if count <= 1 {
        return vec![0.0];
    }
    (0..count)
        .map(|index| {
            let t = index as f64 / (count - 1) as f64;
            -total_m / 2.0 + total_m * t
        })
        .collect()
}

/// Bonds each ring member to the next, with a closing edge.
fn ring_bonds(layer: &[u32]) -> Vec<(u32, u32)> {
    let count = layer.len();
    (0..count)
        .map(|index| (layer[index], layer[(index + 1) % count]))
        .collect()
}

/// Bonds matching members of two consecutive layers.
fn column_bonds(previous: &[u32], current: &[u32]) -> Vec<(u32, u32)> {
    previous
        .iter()
        .zip(current.iter())
        .map(|(a, b)| (*a, *b))
        .collect()
}

/// Moves the rim sample nearest to a tooth centre onto the tip circle.
fn snap_nearest_angle(points_m: &mut [[f64; 2]], angle_rad: f64, radius_m: f64) {
    if points_m.is_empty() {
        return;
    }
    let mut best_index = 0usize;
    let mut best_delta_rad = f64::MAX;
    for (index, point) in points_m.iter().enumerate() {
        let point_angle_rad = point[1].atan2(point[0]);
        let mut delta_rad = (point_angle_rad - angle_rad).abs();
        if delta_rad > PI {
            delta_rad = 2.0 * PI - delta_rad;
        }
        if delta_rad < best_delta_rad {
            best_delta_rad = delta_rad;
            best_index = index;
        }
    }
    points_m[best_index] = [radius_m * angle_rad.cos(), radius_m * angle_rad.sin()];
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_parameters() -> ParameterSet {
        ParameterSet::new()
    }

    fn design_profile() -> GearProfile {
        GearProfile::new(0.5e-9, 20, DEFAULT_PRESSURE_ANGLE_RAD).expect("valid profile")
    }

    fn radial_m(position_m: [f64; 3]) -> f64 {
        (position_m[0] * position_m[0] + position_m[1] * position_m[1]).sqrt()
    }

    #[test]
    fn atom_count_is_stable_and_nonzero() {
        let first = SpurGearGenerator
            .generate(&default_parameters())
            .expect("gear");
        let second = SpurGearGenerator
            .generate(&default_parameters())
            .expect("gear");
        assert_eq!(first.atom_count(), second.atom_count());
        assert!(first.atom_count() > 100, "atoms {}", first.atom_count());
        assert!(first.bond_count() > 0, "bonds {}", first.bond_count());
        assert_eq!(first.name, "spur-gear-z20");
        assert_eq!(first.material, "diamondoid");
    }

    #[test]
    fn the_outer_radius_matches_the_design() {
        let part = SpurGearGenerator
            .generate(&default_parameters())
            .expect("gear");
        let design = design_profile();
        let max_m = part
            .topology
            .atoms()
            .map(|atom| radial_m(atom.position_m))
            .fold(0.0_f64, f64::max);
        assert!(
            (max_m - design.outer_radius_m()).abs() <= design.outer_radius_m() * 1.0e-12,
            "outer radius {max_m:e} m, design {:e} m",
            design.outer_radius_m()
        );
    }

    #[test]
    fn the_bore_is_empty_of_atoms() {
        let part = SpurGearGenerator
            .generate(&default_parameters())
            .expect("gear");
        let bore_m = 0.5e-9;
        for atom in part.topology.atoms() {
            let radius_m = radial_m(atom.position_m);
            assert!(
                radius_m > bore_m,
                "atom at radius {radius_m:e} m is inside bore {bore_m:e} m"
            );
        }
    }

    #[test]
    fn a_bore_larger_than_the_gear_is_rejected() {
        let bad = ParameterSet::new().with("bore_radius_m", 4.5e-9);
        assert!(matches!(
            SpurGearGenerator.generate(&bad),
            Err(PartError::BoreTooLarge { .. })
        ));
    }

    #[test]
    fn a_tooth_count_below_undercut_is_rejected() {
        let bad = ParameterSet::new().with("tooth_count", 5.0);
        assert!(matches!(
            SpurGearGenerator.generate(&bad),
            Err(PartError::ToothCountBelowUndercut { .. })
        ));
    }

    #[test]
    fn a_hub_inside_the_bore_is_rejected() {
        let bad = ParameterSet::new().with("hub_radius_m", 0.1e-9);
        assert!(SpurGearGenerator.generate(&bad).is_err());
    }

    #[test]
    fn invalid_inputs_are_errors_not_panics() {
        let fraction = ParameterSet::new().with("tooth_count", 19.5);
        assert!(SpurGearGenerator.generate(&fraction).is_err());
        let zero = ParameterSet::new().with("face_width_m", 0.0);
        assert!(SpurGearGenerator.generate(&zero).is_err());
    }

    #[test]
    fn a_generator_reports_its_identity() {
        assert_eq!(SpurGearGenerator.id(), "spur_gear");
        assert_eq!(SpurGearGenerator.name(), "Atomistic spur gear");
        assert_eq!(SpurGearGenerator.parameters().len(), 7);
    }
}
