use std::f64::consts::PI;
use std::ops::Range;

use nanocad_model::{Atom, Bond, BondType, Element, Part, Topology};
use nanocad_units::Unit;

use crate::error::PartError;
use crate::gear_profile::GearProfile;
use crate::generator::PartGenerator;
use crate::parameter::{ParameterSet, ParameterSpec};

/// The external full-depth addendum coefficient `h_a*`. The addendum is
/// `h_a* m`. Source: ISO 21771 standard full-depth tooth form.
const ADDENDUM_COEFFICIENT: f64 = 1.0;

/// The default pressure angle, 20 degrees, in radians.
const DEFAULT_PRESSURE_ANGLE_RAD: f64 = 20.0 * PI / 180.0;

/// The carrier bond spacing heuristic, as a fraction of the module.
const CARRIER_SPACING_MODULE_FRACTION: f64 = 1.0;

/// Reports whether the coaxial planetary constraint `N_ring = N_sun +
/// 2 N_planet` holds.
pub fn planetary_constraint_holds(
    sun_teeth: usize,
    planet_teeth: usize,
    ring_teeth: usize,
) -> bool {
    ring_teeth == sun_teeth + 2 * planet_teeth
}

/// The derived geometry of an equally spaced planetary gear set.
///
/// The ring tooth count is derived from the coaxial constraint, so the set is
/// always self-consistent. [`PlanetaryDesign::new`] then checks the
/// equal-spacing assembly condition and the planet-to-planet clearance. All
/// lengths are SI metres.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanetaryDesign {
    module_m: f64,
    sun_teeth: usize,
    planet_teeth: usize,
    planet_count: usize,
    pressure_angle_rad: f64,
}

impl PlanetaryDesign {
    /// Builds a design and rejects a set that cannot assemble.
    pub fn new(
        module_m: f64,
        sun_teeth: usize,
        planet_teeth: usize,
        planet_count: usize,
        pressure_angle_rad: f64,
    ) -> Result<Self, PartError> {
        if !module_m.is_finite() || module_m <= 0.0 {
            return Err(PartError::InvalidGeometry(format!(
                "module {module_m} m must be finite and positive"
            )));
        }
        if sun_teeth < 3 || planet_teeth < 3 {
            return Err(PartError::InvalidGeometry(format!(
                "sun teeth {sun_teeth} and planet teeth {planet_teeth} must each be at least 3"
            )));
        }
        if planet_count == 0 {
            return Err(PartError::InvalidGeometry(
                "planet count must be at least one".to_owned(),
            ));
        }
        if !pressure_angle_rad.is_finite()
            || pressure_angle_rad <= 0.0
            || pressure_angle_rad >= std::f64::consts::FRAC_PI_2
        {
            return Err(PartError::InvalidGeometry(format!(
                "pressure angle {pressure_angle_rad} rad must be in (0, pi/2)"
            )));
        }
        let design = Self {
            module_m,
            sun_teeth,
            planet_teeth,
            planet_count,
            pressure_angle_rad,
        };
        GearProfile::new(module_m, sun_teeth, pressure_angle_rad)?.check_undercut()?;
        GearProfile::new(module_m, planet_teeth, pressure_angle_rad)?.check_undercut()?;
        design.check_assembly()?;
        Ok(design)
    }

    /// The module, in metres.
    pub fn module_m(&self) -> f64 {
        self.module_m
    }

    /// The sun tooth count.
    pub fn sun_teeth(&self) -> usize {
        self.sun_teeth
    }

    /// The planet tooth count.
    pub fn planet_teeth(&self) -> usize {
        self.planet_teeth
    }

    /// The planet count.
    pub fn planet_count(&self) -> usize {
        self.planet_count
    }

    /// The pressure angle, in radians.
    pub fn pressure_angle_rad(&self) -> f64 {
        self.pressure_angle_rad
    }

    /// The ring tooth count from the coaxial constraint.
    pub fn ring_teeth(&self) -> usize {
        self.sun_teeth + 2 * self.planet_teeth
    }

    /// Reports whether the coaxial constraint holds for this design. It always
    /// holds, because [`PlanetaryDesign::ring_teeth`] derives the ring count.
    pub fn constraint_holds(&self) -> bool {
        planetary_constraint_holds(self.sun_teeth, self.planet_teeth, self.ring_teeth())
    }

    /// The sun pitch radius `m N_sun / 2`, in metres.
    pub fn sun_pitch_radius_m(&self) -> f64 {
        self.module_m * self.sun_teeth as f64 / 2.0
    }

    /// The planet pitch radius `m N_planet / 2`, in metres.
    pub fn planet_pitch_radius_m(&self) -> f64 {
        self.module_m * self.planet_teeth as f64 / 2.0
    }

    /// The ring pitch radius `m N_ring / 2`, in metres.
    pub fn ring_pitch_radius_m(&self) -> f64 {
        self.module_m * self.ring_teeth() as f64 / 2.0
    }

    /// The carrier (pin) radius `r_sun + r_planet`, in metres.
    pub fn carrier_radius_m(&self) -> f64 {
        self.sun_pitch_radius_m() + self.planet_pitch_radius_m()
    }

    /// The planet outside radius `m (N_planet / 2 + h_a*)`, in metres.
    pub fn planet_outer_radius_m(&self) -> f64 {
        self.planet_pitch_radius_m() + ADDENDUM_COEFFICIENT * self.module_m
    }

    /// The gear ratio `(N_sun + N_ring) / N_sun` with the ring held fixed.
    pub fn gear_ratio(&self) -> f64 {
        (self.sun_teeth + self.ring_teeth()) as f64 / self.sun_teeth as f64
    }

    /// The angular position of a planet on the carrier, in radians.
    pub fn planet_angle_rad(&self, index: usize) -> f64 {
        2.0 * PI * index as f64 / self.planet_count as f64
    }

    /// The centre of a planet on the carrier, in metres.
    pub fn planet_center_m(&self, index: usize) -> [f64; 2] {
        let angle_rad = self.planet_angle_rad(index);
        [
            self.carrier_radius_m() * angle_rad.cos(),
            self.carrier_radius_m() * angle_rad.sin(),
        ]
    }

    /// The smallest centre-to-centre distance between two planets, in metres.
    ///
    /// For equally spaced planets it is the chord between neighbours,
    /// `2 r_c sin(pi / n)`.
    pub fn minimum_planet_spacing_m(&self) -> f64 {
        if self.planet_count < 2 {
            return f64::INFINITY;
        }
        2.0 * self.carrier_radius_m() * (PI / self.planet_count as f64).sin()
    }

    /// The planet-to-planet clearance radius `2 r_planet_outer`, in metres.
    pub fn planet_clearance_m(&self) -> f64 {
        2.0 * self.planet_outer_radius_m()
    }

    /// Checks the equal-spacing and clearance rules, and returns every failure.
    pub fn check_assembly(&self) -> Result<(), PartError> {
        let total_teeth = self.sun_teeth + self.ring_teeth();
        if !total_teeth.is_multiple_of(self.planet_count) {
            return Err(PartError::PlanetSpacingNotPossible {
                planet_count: self.planet_count,
                total_teeth,
            });
        }
        let spacing_m = self.minimum_planet_spacing_m();
        let clearance_m = self.planet_clearance_m();
        if spacing_m < clearance_m {
            return Err(PartError::PlanetsOverlap {
                center_m: spacing_m,
                clearance_m,
            });
        }
        Ok(())
    }
}

static PLANETARY_PARAMETERS: &[ParameterSpec] = &[
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
        "sun_teeth",
        None,
        24.0,
        3.0,
        120.0,
        true,
        "sun gear tooth count",
    ),
    ParameterSpec::new(
        "planet_teeth",
        None,
        18.0,
        3.0,
        120.0,
        true,
        "planet gear tooth count",
    ),
    ParameterSpec::new(
        "planet_count",
        None,
        3.0,
        1.0,
        12.0,
        true,
        "number of equally spaced planets",
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
        "samples_per_flank",
        None,
        4.0,
        1.0,
        64.0,
        true,
        "samples along each involute flank",
    ),
    ParameterSpec::new(
        "samples_per_arc",
        None,
        2.0,
        1.0,
        64.0,
        true,
        "samples along each tip and root arc",
    ),
    ParameterSpec::new(
        "carrier_offset_m",
        Some(Unit::Metre),
        0.8e-9,
        0.0,
        5.0e-9,
        false,
        "axial offset of the carrier plate in metres",
    ),
];

/// A generated planetary set with the atom ranges of each body.
#[derive(Clone, Debug, PartialEq)]
pub struct PlanetarySet {
    /// The generated part: sun, planets, ring, and carrier.
    pub part: Part,
    /// The derived design.
    pub design: PlanetaryDesign,
    /// The atom range of the sun.
    pub sun_atoms: Range<usize>,
    /// The atom ranges of the planets, in carrier order.
    pub planet_atoms: Vec<Range<usize>>,
    /// The atom range of the ring.
    pub ring_atoms: Range<usize>,
    /// The atom range of the carrier.
    pub carrier_atoms: Range<usize>,
}

/// Builds a skeletal planetary gear set: sun, planets, ring, and carrier.
///
/// The sun and planets are external involute outlines. The ring is the same
/// external outline reflected through its pitch circle, which turns it into an
/// internal gear. The carrier holds one pin at each planet centre and sits at
/// an axial offset, as a real carrier plate does. The part is skeletal, in the
/// style of the other generators, not a filled solid.
#[derive(Clone, Copy, Debug, Default)]
pub struct PlanetaryGenerator;

impl PlanetaryGenerator {
    /// Resolves the inputs and builds the set with its body ranges.
    pub fn build(&self, parameters: &ParameterSet) -> Result<PlanetarySet, PartError> {
        let resolved = self.resolve(parameters)?;
        let module_m = resolved.require("module_m")?;
        let sun_teeth = resolved.require("sun_teeth")? as usize;
        let planet_teeth = resolved.require("planet_teeth")? as usize;
        let planet_count = resolved.require("planet_count")? as usize;
        let pressure_angle_rad = resolved.require("pressure_angle_rad")?;
        let samples_per_flank = resolved.require("samples_per_flank")? as usize;
        let samples_per_arc = resolved.require("samples_per_arc")? as usize;
        let carrier_offset_m = resolved.require("carrier_offset_m")?;

        let design = PlanetaryDesign::new(
            module_m,
            sun_teeth,
            planet_teeth,
            planet_count,
            pressure_angle_rad,
        )?;

        let mut topology = Topology::new();

        let sun_profile = GearProfile::new(module_m, sun_teeth, pressure_angle_rad)?;
        let sun_points =
            trim_closed_loop(sun_profile.outline_points(samples_per_flank, samples_per_arc)?);
        let sun_indices = add_closed_loop(&mut topology, &sun_points, 0.0, 0.0, [0.0, 0.0], "sun")?;
        let sun_atoms = index_range(&sun_indices);

        let planet_profile = GearProfile::new(module_m, planet_teeth, pressure_angle_rad)?;
        let planet_points =
            trim_closed_loop(planet_profile.outline_points(samples_per_flank, samples_per_arc)?);
        let mut planet_atoms = Vec::with_capacity(planet_count);
        for planet in 0..planet_count {
            let angle_rad = design.planet_angle_rad(planet);
            let center_m = design.planet_center_m(planet);
            let indices = add_closed_loop(
                &mut topology,
                &planet_points,
                0.0,
                angle_rad,
                center_m,
                "planet",
            )?;
            planet_atoms.push(index_range(&indices));
        }

        let ring_profile = GearProfile::new(module_m, design.ring_teeth(), pressure_angle_rad)?;
        let external_ring =
            trim_closed_loop(ring_profile.outline_points(samples_per_flank, samples_per_arc)?);
        let internal_ring = reflect_to_internal(&external_ring, design.ring_pitch_radius_m());
        let ring_indices =
            add_closed_loop(&mut topology, &internal_ring, 0.0, 0.0, [0.0, 0.0], "ring")?;
        let ring_atoms = index_range(&ring_indices);

        let carrier_start = topology.atom_count();
        add_carrier(&mut topology, &design, carrier_offset_m)?;
        let carrier_atoms = carrier_start..topology.atom_count();

        let name = format!("planetary-s{sun_teeth}-p{planet_teeth}x{planet_count}");
        let mut part = Part::new(name, topology).with_material("diamondoid");
        part.metadata
            .insert("generator".to_owned(), "planetary".to_owned());
        part.metadata
            .insert("ring_teeth".to_owned(), design.ring_teeth().to_string());
        part.metadata.insert(
            "gear_ratio".to_owned(),
            format!("{:.6}", design.gear_ratio()),
        );

        Ok(PlanetarySet {
            part,
            design,
            sun_atoms,
            planet_atoms,
            ring_atoms,
            carrier_atoms,
        })
    }
}

impl PartGenerator for PlanetaryGenerator {
    fn id(&self) -> &'static str {
        "planetary"
    }

    fn name(&self) -> &'static str {
        "Planetary gear set"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        PLANETARY_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        Ok(self.build(parameters)?.part)
    }
}

/// Trims a repeated closing point from a closed polyline.
fn trim_closed_loop(mut points_m: Vec<[f64; 2]>) -> Vec<[f64; 2]> {
    if points_m.len() >= 2 {
        let first = points_m[0];
        let last = points_m[points_m.len() - 1];
        if (first[0] - last[0]).abs() < 1.0e-15 && (first[1] - last[1]).abs() < 1.0e-15 {
            points_m.pop();
        }
    }
    points_m
}

/// Maps an external gear outline to an internal one by reflection through the
/// pitch circle. A standard external tooth at radius `r` becomes an internal
/// tooth at radius `2 r_pitch - r`, so the addendum points inward.
fn reflect_to_internal(points_m: &[[f64; 2]], pitch_radius_m: f64) -> Vec<[f64; 2]> {
    points_m
        .iter()
        .map(|point| {
            let radius_m = (point[0] * point[0] + point[1] * point[1]).sqrt();
            if radius_m <= 0.0 {
                return *point;
            }
            let scale = (2.0 * pitch_radius_m - radius_m) / radius_m;
            [point[0] * scale, point[1] * scale]
        })
        .collect()
}

/// Adds a closed loop of carbon atoms and bonds it around. Returns the indices.
fn add_closed_loop(
    topology: &mut Topology,
    points_m: &[[f64; 2]],
    z_m: f64,
    rotation_rad: f64,
    center_m: [f64; 2],
    atom_type: &str,
) -> Result<Vec<u32>, PartError> {
    let (sin_rot, cos_rot) = rotation_rad.sin_cos();
    let mut indices = Vec::with_capacity(points_m.len());
    for point in points_m {
        let x_m = point[0] * cos_rot - point[1] * sin_rot + center_m[0];
        let y_m = point[0] * sin_rot + point[1] * cos_rot + center_m[1];
        indices.push(topology.add_atom(Atom::new(
            Element::CARBON,
            [x_m, y_m, z_m],
            0.0,
            atom_type,
        )));
    }
    let count = indices.len();
    for position in 0..count {
        let u = indices[position];
        let v = indices[(position + 1) % count];
        topology.add_bond(Bond::new(u, v, 1, BondType::Single))?;
    }
    Ok(indices)
}

/// Adds a uniform closed ring of carbon atoms and bonds it around.
fn add_ring(
    topology: &mut Topology,
    center_m: [f64; 2],
    radius_m: f64,
    samples: usize,
    z_m: f64,
    atom_type: &str,
) -> Result<Vec<u32>, PartError> {
    let count = samples.max(3);
    let mut indices = Vec::with_capacity(count);
    for sample in 0..count {
        let angle_rad = 2.0 * PI * sample as f64 / count as f64;
        indices.push(topology.add_atom(Atom::new(
            Element::CARBON,
            [
                center_m[0] + radius_m * angle_rad.cos(),
                center_m[1] + radius_m * angle_rad.sin(),
                z_m,
            ],
            0.0,
            atom_type,
        )));
    }
    for position in 0..count {
        topology.add_bond(Bond::new(
            indices[position],
            indices[(position + 1) % count],
            1,
            BondType::Single,
        ))?;
    }
    Ok(indices)
}

/// Adds the carrier: a race ring at the pin radius and one pin per planet.
fn add_carrier(
    topology: &mut Topology,
    design: &PlanetaryDesign,
    offset_m: f64,
) -> Result<(), PartError> {
    let carrier_radius_m = design.carrier_radius_m();
    let circumference_m = 2.0 * PI * carrier_radius_m;
    let spacing_m = design.module_m() * CARRIER_SPACING_MODULE_FRACTION;
    let race_samples = ((circumference_m / spacing_m).round() as usize).max(8);
    let race = add_ring(
        topology,
        [0.0, 0.0],
        carrier_radius_m,
        race_samples,
        offset_m,
        "carrier",
    )?;

    let pin_radius_m = 0.4 * design.planet_pitch_radius_m();
    for planet in 0..design.planet_count() {
        let center_m = design.planet_center_m(planet);
        let pin = add_ring(topology, center_m, pin_radius_m, 8, offset_m, "carrier")?;
        connect_nearest(topology, &race, &pin)?;
    }
    Ok(())
}

/// Bonds the closest atom of two rings.
fn connect_nearest(topology: &mut Topology, a: &[u32], b: &[u32]) -> Result<(), PartError> {
    let mut best: Option<(u32, u32)> = None;
    let mut best_m = f64::MAX;
    for &u in a {
        let Some(position_u_m) = topology.position_m(u as usize) else {
            continue;
        };
        for &v in b {
            let Some(position_v_m) = topology.position_m(v as usize) else {
                continue;
            };
            let dx = position_u_m[0] - position_v_m[0];
            let dy = position_u_m[1] - position_v_m[1];
            let dz = position_u_m[2] - position_v_m[2];
            let distance_m = (dx * dx + dy * dy + dz * dz).sqrt();
            if distance_m < best_m {
                best_m = distance_m;
                best = Some((u, v));
            }
        }
    }
    if let Some((u, v)) = best {
        topology.add_bond(Bond::new(u, v, 1, BondType::Single))?;
    }
    Ok(())
}

/// Returns the contiguous index range that a list of atom indices covers.
fn index_range(indices: &[u32]) -> Range<usize> {
    match (indices.first(), indices.last()) {
        (Some(first), Some(last)) => *first as usize..*last as usize + 1,
        _ => 0..0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_parameters() -> ParameterSet {
        ParameterSet::new().with("planet_count", 3.0)
    }

    fn distance_m(a_m: [f64; 3], b_m: [f64; 3]) -> f64 {
        let dx = a_m[0] - b_m[0];
        let dy = a_m[1] - b_m[1];
        let dz = a_m[2] - b_m[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    #[test]
    fn the_coaxial_constraint_holds_for_a_generated_set() {
        let build = PlanetaryGenerator
            .build(&default_parameters())
            .expect("generate");
        let design = build.design;
        assert_eq!(design.ring_teeth(), 24 + 2 * 18);
        assert!(design.constraint_holds());
        assert!(planetary_constraint_holds(24, 18, 60));
        assert!(!planetary_constraint_holds(24, 18, 61));
        assert_eq!(
            build.part.metadata.get("ring_teeth").map(String::as_str),
            Some("60")
        );
    }

    #[test]
    fn the_gear_ratio_is_sun_plus_ring_over_sun() {
        let build = PlanetaryGenerator
            .build(&default_parameters())
            .expect("generate");
        let expected = (24.0 + 60.0) / 24.0;
        assert!((build.design.gear_ratio() - expected).abs() < 1.0e-12);
        assert_eq!(build.design.gear_ratio(), 3.5);
    }

    #[test]
    fn the_planets_do_not_overlap() {
        let build = PlanetaryGenerator
            .build(&default_parameters())
            .expect("generate");
        let design = build.design;
        assert!(design.minimum_planet_spacing_m() > design.planet_clearance_m());
        assert!(build.planet_atoms.len() >= 2);
        for first in 0..build.planet_atoms.len() {
            for second in (first + 1)..build.planet_atoms.len() {
                let range_a = &build.planet_atoms[first];
                let range_b = &build.planet_atoms[second];
                let mut minimum_m = f64::MAX;
                for a in range_a.clone() {
                    let Some(a_m) = build.part.topology.position_m(a) else {
                        continue;
                    };
                    for b in range_b.clone() {
                        let Some(b_m) = build.part.topology.position_m(b) else {
                            continue;
                        };
                        let delta_m = distance_m(a_m, b_m);
                        if delta_m < minimum_m {
                            minimum_m = delta_m;
                        }
                    }
                }
                assert!(
                    minimum_m > 1.0e-10,
                    "planets {first} and {second} are {minimum_m:e} m apart"
                );
            }
        }
    }

    #[test]
    fn the_planets_sit_on_the_carrier_circle() {
        let build = PlanetaryGenerator
            .build(&default_parameters())
            .expect("generate");
        let design = build.design;
        for planet in 0..design.planet_count() {
            let expected_m = design.planet_center_m(planet);
            let range = &build.planet_atoms[planet];
            let mean_x: f64 = range
                .clone()
                .filter_map(|index| build.part.topology.position_m(index))
                .map(|position_m| position_m[0])
                .sum::<f64>()
                / range.len() as f64;
            let mean_y: f64 = range
                .clone()
                .filter_map(|index| build.part.topology.position_m(index))
                .map(|position_m| position_m[1])
                .sum::<f64>()
                / range.len() as f64;
            assert!((mean_x - expected_m[0]).abs() < 1.0e-11);
            assert!((mean_y - expected_m[1]).abs() < 1.0e-11);
        }
    }

    #[test]
    fn the_carrier_sits_at_the_axial_offset() {
        let build = PlanetaryGenerator
            .build(&default_parameters())
            .expect("generate");
        assert!(!build.carrier_atoms.is_empty());
        for index in build.carrier_atoms.clone() {
            let z_m = build.part.topology.position_m(index).expect("carrier atom")[2];
            assert!((z_m - 0.8e-9).abs() < 1.0e-15);
        }
    }

    #[test]
    fn the_part_has_every_body_with_no_duplicate_bonds() {
        let build = PlanetaryGenerator
            .build(&default_parameters())
            .expect("generate");
        assert!(build.part.atom_count() > 100);
        assert!(build.part.bond_count() > 0);
        assert_eq!(build.part.name, "planetary-s24-p18x3");
        assert_eq!(build.part.material, "diamondoid");
        let mut seen = std::collections::BTreeSet::new();
        for bond in build.part.topology.bonds() {
            let key = if bond.u <= bond.v {
                (bond.u, bond.v)
            } else {
                (bond.v, bond.u)
            };
            assert!(seen.insert(key), "duplicate bond {key:?}");
        }
    }

    #[test]
    fn a_bad_planet_count_is_rejected() {
        let bad = ParameterSet::new().with("planet_count", 5.0);
        assert!(matches!(
            PlanetaryGenerator.build(&bad),
            Err(PartError::PlanetSpacingNotPossible { .. })
        ));
    }

    #[test]
    fn overlapping_planets_are_rejected() {
        let bad = ParameterSet::new()
            .with("sun_teeth", 18.0)
            .with("planet_teeth", 18.0)
            .with("planet_count", 6.0);
        assert!(matches!(
            PlanetaryGenerator.build(&bad),
            Err(PartError::PlanetsOverlap { .. })
        ));
    }

    #[test]
    fn invalid_inputs_are_errors_not_panics() {
        assert!(PlanetaryGenerator
            .build(&ParameterSet::new().with("module_m", 0.0))
            .is_err());
        assert!(PlanetaryGenerator
            .build(&ParameterSet::new().with("planet_count", 0.0))
            .is_err());
        assert!(PlanetaryGenerator
            .build(&ParameterSet::new().with("planet_teeth", 2.0))
            .is_err());
        assert!(PlanetaryGenerator
            .build(&ParameterSet::new().with("nope", 1.0))
            .is_err());
    }

    #[test]
    fn a_generator_reports_its_identity() {
        assert_eq!(PlanetaryGenerator.id(), "planetary");
        assert_eq!(PlanetaryGenerator.name(), "Planetary gear set");
        assert_eq!(PlanetaryGenerator.parameters().len(), 8);
    }
}
