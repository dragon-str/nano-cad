//! The selective gear pump (ADR-0086).
//!
//! Two equal solid diamond spur gears mesh inside a figure-eight housing. The
//! housing is the membrane. A tooth gap is a pocket for one small molecule.
//! On the suction side the teeth separate and each gap opens to the suction
//! port. The gap carries its guest round the housing bore to the discharge
//! side. There the teeth close the gap and push the guest out. The mesh forces
//! the release, so the binding site does not need to let go by itself.
//!
//! Each gear turns on a post of the housing. The posts join the two side
//! plates, so the housing is one connected lattice. The gears are held
//! radially by the posts and axially by the plates.
//!
//! Every clearance is a carbon-to-carbon distance. A hydrogen cap reaches
//! 1.09 Angstrom from its carbon, so two capped faces at 4.6 Angstrom leave
//! about 2.4 Angstrom between the hydrogen atoms. That is less than the
//! 2.8 Angstrom diameter of a water molecule. M27-03 measures the real gaps.
//!
//! Scope limits. The part holds no binding site yet (M27-04). It models no
//! drive, no torque, no contact force and no solvent.

use std::collections::{HashMap, HashSet, VecDeque};
use std::f64::consts::PI;
use std::ops::Range;

use nanocad_model::{Atom, Element, Part, Topology};
use nanocad_units::Unit;

use crate::diamond_solid::{self, DIAMOND_LATTICE_CONSTANT_M};
use crate::error::PartError;
use crate::gear_profile::GearProfile;
use crate::generator::PartGenerator;
use crate::lattice_fill::lattice_atoms;
use crate::parameter::{ParameterSet, ParameterSpec};
use crate::planetary::{add_gear_local, index_range, trim_closed_loop};
use crate::shape::{Box3, Cylinder, Difference, Profile, Solid, Union};

/// The default pressure angle, 30 degrees, in radians. At 30 degrees the
/// undercut limit is eight teeth, so a small gear stays sound.
const DEFAULT_PRESSURE_ANGLE_RAD: f64 = 30.0 * PI / 180.0;

/// The samples along each involute flank of the outline.
const SAMPLES_PER_FLANK: usize = 6;

/// The samples along each tip and root arc of the outline.
const SAMPLES_PER_ARC: usize = 3;

/// The iron-carbon bond of a root site, in metres. An iron-carbon sigma bond
/// in an organo-iron complex is 1.9e-10 to 2.0e-10 m.
pub const FE_C_BOND_M: f64 = 1.95e-10;

/// A lattice carbon with fewer first-shell neighbours than this is removed.
/// A carbon with one bond is a loose tag on the surface, not a part.
const MIN_CARBON_NEIGHBOURS: usize = 2;

static GEAR_PUMP_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new(
        "module_m",
        Some(Unit::Metre),
        0.6e-9,
        0.3e-9,
        3.0e-9,
        false,
        "gear module in metres; it sets the tooth gap, which is the pocket",
    ),
    ParameterSpec::new(
        "teeth",
        None,
        8.0,
        8.0,
        40.0,
        true,
        "tooth count of each gear",
    ),
    ParameterSpec::new(
        "pressure_angle_rad",
        None,
        DEFAULT_PRESSURE_ANGLE_RAD,
        20.0 * PI / 180.0,
        35.0 * PI / 180.0,
        false,
        "pressure angle in radians",
    ),
    ParameterSpec::new(
        "layers",
        None,
        8.0,
        2.0,
        64.0,
        true,
        "number of atomic (001) layers in the gear thickness",
    ),
    ParameterSpec::new(
        "backlash_m",
        Some(Unit::Metre),
        4.0e-10,
        0.0,
        2.0e-9,
        false,
        "tooth thinning of each gear at the pitch circle",
    ),
    ParameterSpec::new(
        "centre_offset_m",
        Some(Unit::Metre),
        4.0e-10,
        0.0,
        2.0e-9,
        false,
        "extra centre distance beyond the standard m z; it opens the tip-to-root gap",
    ),
    ParameterSpec::new(
        "clearance_m",
        Some(Unit::Metre),
        4.6e-10,
        2.0e-10,
        2.0e-9,
        false,
        "carbon-to-carbon clearance between a gear and the housing",
    ),
    ParameterSpec::new(
        "bore_radius_m",
        Some(Unit::Metre),
        1.0e-9,
        0.3e-9,
        5.0e-9,
        false,
        "radius of the gear bore that takes the housing post",
    ),
    ParameterSpec::new(
        "wall_m",
        Some(Unit::Metre),
        0.8e-9,
        0.3e-9,
        5.0e-9,
        false,
        "thickness of the housing wall and of each side plate",
    ),
    ParameterSpec::new(
        "metal_sites",
        None,
        1.0,
        0.0,
        1.0,
        true,
        "1 puts an iron site in the root of every tooth gap, 0 leaves the root capped",
    ),
    ParameterSpec::new(
        "port_width_m",
        Some(Unit::Metre),
        1.6e-9,
        0.3e-9,
        10.0e-9,
        false,
        "width of the suction port and of the discharge port",
    ),
];

/// The derived geometry of a gear pump. Every length is in metres.
#[derive(Clone, Debug, PartialEq)]
pub struct GearPumpDesign {
    profile: GearProfile,
    layers: usize,
    backlash_m: f64,
    centre_offset_m: f64,
    clearance_m: f64,
    bore_radius_m: f64,
    wall_m: f64,
    port_width_m: f64,
}

impl GearPumpDesign {
    /// Derives the design from resolved parameters and checks it.
    pub fn from_parameters(resolved: &ParameterSet) -> Result<Self, PartError> {
        let module_m = resolved.require("module_m")?;
        let teeth = resolved.require("teeth")?.round() as usize;
        let pressure_angle_rad = resolved.require("pressure_angle_rad")?;
        let profile = GearProfile::new(module_m, teeth, pressure_angle_rad)?;
        profile.check_undercut()?;
        let design = GearPumpDesign {
            profile,
            layers: resolved.require("layers")?.round() as usize,
            backlash_m: resolved.require("backlash_m")?,
            centre_offset_m: resolved.require("centre_offset_m")?,
            clearance_m: resolved.require("clearance_m")?,
            bore_radius_m: resolved.require("bore_radius_m")?,
            wall_m: resolved.require("wall_m")?,
            port_width_m: resolved.require("port_width_m")?,
        };
        design.check()?;
        Ok(design)
    }

    fn check(&self) -> Result<(), PartError> {
        if self.layers < 2 {
            return Err(PartError::InvalidGeometry(
                "a pump gear needs at least two layers".to_owned(),
            ));
        }
        let rim_m = self.root_radius_m() - self.bore_radius_m;
        if rim_m < DIAMOND_LATTICE_CONSTANT_M {
            return Err(PartError::BoreTooLarge {
                bore_m: self.bore_radius_m,
                root_m: self.root_radius_m(),
            });
        }
        if self.post_radius_m() < 0.5 * DIAMOND_LATTICE_CONSTANT_M {
            return Err(PartError::InvalidGeometry(format!(
                "the post radius {:.3e} m is too thin for a lattice; widen the bore",
                self.post_radius_m()
            )));
        }
        if self.tip_to_root_m() < self.clearance_m {
            return Err(PartError::InvalidGeometry(format!(
                "the tip-to-root gap {:.3e} m is below the clearance {:.3e} m; \
                 raise centre_offset_m",
                self.tip_to_root_m(),
                self.clearance_m
            )));
        }
        if self.port_width_m >= self.centre_distance_m() {
            return Err(PartError::InvalidGeometry(
                "the port is wider than the centre distance".to_owned(),
            ));
        }
        Ok(())
    }

    /// The gear profile, before the backlash thinning.
    pub fn profile(&self) -> &GearProfile {
        &self.profile
    }

    /// The tooth count of each gear.
    pub fn teeth(&self) -> usize {
        self.profile.teeth()
    }

    /// The module, in metres.
    pub fn module_m(&self) -> f64 {
        self.profile.module_m()
    }

    /// The number of atomic layers in a gear.
    pub fn layers(&self) -> usize {
        self.layers
    }

    /// The tooth thinning of each gear, in metres.
    pub fn backlash_m(&self) -> f64 {
        self.backlash_m
    }

    /// The gear-to-housing clearance, in metres.
    pub fn clearance_m(&self) -> f64 {
        self.clearance_m
    }

    /// The pitch radius of each gear, in metres.
    pub fn pitch_radius_m(&self) -> f64 {
        self.profile.pitch_radius_m()
    }

    /// The tip radius of each gear, in metres.
    pub fn outer_radius_m(&self) -> f64 {
        self.profile.outer_radius_m()
    }

    /// The root radius of each gear, in metres.
    pub fn root_radius_m(&self) -> f64 {
        self.profile.root_radius_m()
    }

    /// The gear bore radius, in metres.
    pub fn bore_radius_m(&self) -> f64 {
        self.bore_radius_m
    }

    /// The radius of each housing post, in metres.
    pub fn post_radius_m(&self) -> f64 {
        self.bore_radius_m - self.clearance_m
    }

    /// The distance between the two gear axes, `m z + offset`, in metres.
    pub fn centre_distance_m(&self) -> f64 {
        2.0 * self.pitch_radius_m() + self.centre_offset_m
    }

    /// The radial gap between one tooth tip and the root of the other gear.
    pub fn tip_to_root_m(&self) -> f64 {
        self.centre_distance_m() - self.outer_radius_m() - self.root_radius_m()
    }

    /// The centre of gear A (the driver), in metres.
    pub fn gear_a_center_m(&self) -> [f64; 2] {
        [-0.5 * self.centre_distance_m(), 0.0]
    }

    /// The centre of gear B (the follower), in metres.
    pub fn gear_b_center_m(&self) -> [f64; 2] {
        [0.5 * self.centre_distance_m(), 0.0]
    }

    /// The fixed phase of gear B, in radians.
    ///
    /// Gear A shows a tooth tip on the line of centres. Gear B must show a
    /// tooth space there. With an even tooth count the outline has a tip at
    /// `PI` too, so gear B turns by half a tooth pitch.
    pub fn gear_b_phase_rad(&self) -> f64 {
        if self.teeth() % 2 == 1 {
            0.0
        } else {
            PI / self.teeth() as f64
        }
    }

    /// The angle of gear B when gear A is at `gear_a_rad`. The gears are equal,
    /// so they turn at the same rate in opposite senses.
    pub fn gear_b_angle_rad(&self, gear_a_rad: f64) -> f64 {
        self.gear_b_phase_rad() - gear_a_rad
    }

    /// The half-thickness of a gear, from carbon plane to carbon plane.
    pub fn gear_half_thickness_m(&self) -> f64 {
        0.5 * (self.layers - 1) as f64 * diamond_solid::diamond_plane_spacing_m()
    }

    /// The radius of each lobe of the figure-eight cavity, in metres.
    pub fn cavity_radius_m(&self) -> f64 {
        self.outer_radius_m() + self.clearance_m
    }

    /// The half-height of the cavity along the gear axis, in metres.
    pub fn cavity_half_height_m(&self) -> f64 {
        self.gear_half_thickness_m() + self.clearance_m
    }

    /// The half-height of the whole housing along the gear axis, in metres.
    pub fn housing_half_height_m(&self) -> f64 {
        self.cavity_half_height_m() + self.wall_m
    }

    /// The outer radius of each housing lobe, in metres.
    pub fn housing_radius_m(&self) -> f64 {
        self.cavity_radius_m() + self.wall_m
    }

    /// The width of each port, in metres.
    pub fn port_width_m(&self) -> f64 {
        self.port_width_m
    }

    /// The area of one tooth gap in the gear plane, in square metres.
    ///
    /// The gap is the annulus between the root and tip circles, less the
    /// teeth. The estimate takes half the annulus, then adds the backlash
    /// thinning. It is an estimate for sizing, not a measured volume.
    pub fn tooth_gap_area_m2(&self) -> f64 {
        let tip_m = self.outer_radius_m();
        let root_m = self.root_radius_m();
        let annulus_m2 = PI * (tip_m * tip_m - root_m * root_m);
        let depth_m = tip_m - root_m;
        0.5 * annulus_m2 / self.teeth() as f64 + self.backlash_m * depth_m
    }
}

/// A generated gear pump with the atom ranges of each body.
#[derive(Clone, Debug, PartialEq)]
pub struct GearPumpSet {
    /// The generated part: housing, gear A and gear B.
    pub part: Part,
    /// The derived design.
    pub design: GearPumpDesign,
    /// The atom range of the housing.
    pub housing_atoms: Range<usize>,
    /// The atom range of gear A, the driver.
    pub gear_a_atoms: Range<usize>,
    /// The atom range of gear B, the follower.
    pub gear_b_atoms: Range<usize>,
    /// The iron atoms of the root sites: gear A first, then gear B, one per
    /// tooth gap in gap order. Empty when `metal_sites` is 0.
    pub metal_sites: Vec<usize>,
}

/// Builds the selective gear pump: a figure-eight housing and two meshing
/// solid diamond spur gears.
///
/// Gear A sits at `-x`, gear B at `+x`. The suction port opens at `-y` and
/// the discharge port at `+y`, both at the mesh waist. Gear A turns
/// counter-clockwise and gear B clockwise: the gaps then carry a guest from
/// the suction side round the outer bore to the discharge side.
#[derive(Clone, Copy, Debug, Default)]
pub struct GearPumpGenerator;

impl GearPumpGenerator {
    /// Resolves the inputs and builds the pump with its body ranges.
    pub fn build(&self, parameters: &ParameterSet) -> Result<GearPumpSet, PartError> {
        let resolved = self.resolve(parameters)?;
        let design = GearPumpDesign::from_parameters(&resolved)?;
        let mut topology = Topology::new();

        let housing_indices = add_housing(&mut topology, &design)?;
        let housing_atoms = index_range(&housing_indices);

        let gear_local = gear_local_sites(&design)?;
        let gear_a_indices = add_gear_local(
            &mut topology,
            &gear_local,
            0.0,
            design.gear_a_center_m(),
            design.layers(),
            "gear_a",
        )?;
        let gear_a_atoms = index_range(&gear_a_indices);
        let gear_b_indices = add_gear_local(
            &mut topology,
            &gear_local,
            design.gear_b_phase_rad(),
            design.gear_b_center_m(),
            design.layers(),
            "gear_b",
        )?;
        let gear_b_atoms = index_range(&gear_b_indices);

        let metal_sites = if resolved.require("metal_sites")? >= 0.5 {
            add_metal_sites(&mut topology, &design, &gear_a_atoms, &gear_b_atoms)?
        } else {
            Vec::new()
        };

        let name = format!("gear-pump-z{}", design.teeth());
        let mut part = Part::new(name, topology).with_material("diamond");
        let metadata = &mut part.metadata;
        metadata.insert("generator".to_owned(), "gear_pump".to_owned());
        metadata.insert("teeth".to_owned(), design.teeth().to_string());
        metadata.insert(
            "centre_distance_m".to_owned(),
            format!("{:.6e}", design.centre_distance_m()),
        );
        metadata.insert(
            "pitch_radius_m".to_owned(),
            format!("{:.6e}", design.pitch_radius_m()),
        );
        metadata.insert(
            "tooth_gap_area_m2".to_owned(),
            format!("{:.6e}", design.tooth_gap_area_m2()),
        );

        Ok(GearPumpSet {
            part,
            design,
            housing_atoms,
            gear_a_atoms,
            gear_b_atoms,
            metal_sites,
        })
    }
}

impl PartGenerator for GearPumpGenerator {
    fn id(&self) -> &'static str {
        "gear_pump"
    }

    fn name(&self) -> &'static str {
        "Selective gear pump"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        GEAR_PUMP_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        Ok(self.build(parameters)?.part)
    }
}

/// Turns one root cap of every tooth gap into an iron atom.
///
/// Gear A has a tooth tip at angle 0, so its gaps lie at `(k + 1/2) 2 PI / z`.
/// The cap nearest the root point of each gap becomes iron, moved out along
/// its cap direction to [`FE_C_BOND_M`]. Gear B is the same body turned by
/// its phase, so the same local atom becomes iron there. The iron keeps its
/// open side toward the gap: that side is the binding site (M27-04).
fn add_metal_sites(
    topology: &mut Topology,
    design: &GearPumpDesign,
    gear_a: &Range<usize>,
    gear_b: &Range<usize>,
) -> Result<Vec<usize>, PartError> {
    let center_m = design.gear_a_center_m();
    let teeth = design.teeth();
    let mut host_of: HashMap<usize, usize> = HashMap::new();
    for bond in topology.bonds() {
        let (u, v) = (bond.u as usize, bond.v as usize);
        if topology.element(u) == Some(Element::HYDROGEN) {
            host_of.insert(u, v);
        }
        if topology.element(v) == Some(Element::HYDROGEN) {
            host_of.insert(v, u);
        }
    }
    let probe_radius_m = design.root_radius_m() + 1.0e-10;
    let mut chosen = Vec::with_capacity(teeth);
    for gap in 0..teeth {
        let angle_rad = (gap as f64 + 0.5) * 2.0 * PI / teeth as f64;
        let target_m = [
            center_m[0] + probe_radius_m * angle_rad.cos(),
            center_m[1] + probe_radius_m * angle_rad.sin(),
            0.0,
        ];
        let best = gear_a
            .clone()
            .filter(|index| host_of.contains_key(index))
            .filter_map(|index| {
                let p = topology.position_m(index)?;
                let d = (p[0] - target_m[0]).powi(2)
                    + (p[1] - target_m[1]).powi(2)
                    + (p[2] - target_m[2]).powi(2);
                Some((d, index))
            })
            .min_by(|a, b| a.0.total_cmp(&b.0));
        let Some((_, index)) = best else {
            return Err(PartError::InvalidGeometry(
                "a tooth gap has no root cap for an iron site".to_owned(),
            ));
        };
        chosen.push(index - gear_a.start);
    }
    let mut sites = Vec::with_capacity(2 * teeth);
    for range in [gear_a, gear_b] {
        for &offset in &chosen {
            let cap = range.start + offset;
            let host = *host_of.get(&cap).ok_or_else(|| {
                PartError::InvalidGeometry("the gears do not share an atom order".to_owned())
            })?;
            let (Some(cap_m), Some(host_m)) = (topology.position_m(cap), topology.position_m(host))
            else {
                continue;
            };
            let length_m = ((cap_m[0] - host_m[0]).powi(2)
                + (cap_m[1] - host_m[1]).powi(2)
                + (cap_m[2] - host_m[2]).powi(2))
            .sqrt();
            let iron_m = [
                host_m[0] + FE_C_BOND_M * (cap_m[0] - host_m[0]) / length_m,
                host_m[1] + FE_C_BOND_M * (cap_m[1] - host_m[1]) / length_m,
                host_m[2] + FE_C_BOND_M * (cap_m[2] - host_m[2]) / length_m,
            ];
            topology.set_element(cap, Element::IRON)?;
            topology.set_position_m(cap, iron_m)?;
            sites.push(cap);
        }
    }
    Ok(sites)
}

/// Returns the carbon sites of one gear in its local crystal frame.
///
/// The slab starts at `z = 0`, as [`add_gear_local`] expects. The gear is the
/// thinned involute outline less the bore.
fn gear_local_sites(design: &GearPumpDesign) -> Result<Vec<[f64; 3]>, PartError> {
    let outline_m = trim_closed_loop(
        design
            .profile()
            .with_backlash(design.backlash_m())
            .outline_points(SAMPLES_PER_FLANK, SAMPLES_PER_ARC)?,
    );
    let profile = Profile {
        outline_m,
        internal: false,
        outer_radius_m: None,
    };
    let top_m = (design.layers() - 1) as f64 * diamond_solid::diamond_plane_spacing_m();
    let tip_m = design.outer_radius_m();
    let bore_m = design.bore_radius_m();
    let sites = lattice_atoms([-tip_m, -tip_m, 0.0], [tip_m, tip_m, top_m], &mut |point| {
        point[0].hypot(point[1]) > bore_m && profile.contains_m(point)
    });
    let sites = keep_well_bonded(sites);
    if sites.is_empty() {
        return Err(PartError::InvalidGeometry(
            "the gear parameters leave no lattice site".to_owned(),
        ));
    }
    Ok(sites)
}

/// Adds the housing: two side plates, a figure-eight wall with two ports, and
/// two posts. Returns the carbon and cap indices.
fn add_housing(topology: &mut Topology, design: &GearPumpDesign) -> Result<Vec<u32>, PartError> {
    let [ax, _] = design.gear_a_center_m();
    let [bx, _] = design.gear_b_center_m();
    let full_height_m = 2.0 * design.housing_half_height_m();
    let cavity_height_m = 2.0 * design.cavity_half_height_m();
    let lobe = |x_m: f64, radius_m: f64, height_m: f64| -> Box<dyn Solid> {
        Box::new(Cylinder {
            center_m: [x_m, 0.0, 0.0],
            radius_m,
            height_m,
        })
    };
    let shell: Box<dyn Solid> = Box::new(Union {
        a: lobe(ax, design.housing_radius_m(), full_height_m),
        b: lobe(bx, design.housing_radius_m(), full_height_m),
    });
    let cavity: Box<dyn Solid> = Box::new(Union {
        a: lobe(ax, design.cavity_radius_m(), cavity_height_m),
        b: lobe(bx, design.cavity_radius_m(), cavity_height_m),
    });
    // Each port is a slot through the wall at the mesh waist, the full cavity
    // height. The suction port opens at -y and the discharge port at +y.
    let reach_m = design.housing_radius_m() + DIAMOND_LATTICE_CONSTANT_M;
    let half_width_m = 0.5 * design.port_width_m();
    let half_height_m = design.cavity_half_height_m();
    let ports: Box<dyn Solid> = Box::new(Box3 {
        min_m: [-half_width_m, -reach_m, -half_height_m],
        max_m: [half_width_m, reach_m, half_height_m],
    });
    let body: Box<dyn Solid> = Box::new(Difference {
        outer: shell,
        inner: Box::new(Union {
            a: cavity,
            b: ports,
        }),
    });
    let housing = Union {
        a: body,
        b: Box::new(Union {
            a: lobe(ax, design.post_radius_m(), full_height_m),
            b: lobe(bx, design.post_radius_m(), full_height_m),
        }),
    };
    let sites = keep_well_bonded(crate::lattice_fill::fill_solid(&housing));
    if sites.is_empty() {
        return Err(PartError::InvalidGeometry(
            "the housing parameters leave no lattice site".to_owned(),
        ));
    }
    let mut indices = Vec::with_capacity(sites.len());
    for site_m in sites {
        indices.push(topology.add_atom(Atom::new(Element::CARBON, site_m, 0.0, "housing")));
    }
    let caps = diamond_solid::bond_and_cap(topology, &indices)?;
    indices.extend(caps);
    Ok(indices)
}

/// The integer key of a diamond site on the quarter-cell grid.
fn site_key(point_m: [f64; 3]) -> (i64, i64, i64) {
    let quarter_m = 0.25 * DIAMOND_LATTICE_CONSTANT_M;
    (
        (point_m[0] / quarter_m).round() as i64,
        (point_m[1] / quarter_m).round() as i64,
        (point_m[2] / quarter_m).round() as i64,
    )
}

/// Removes loosely held carbons and keeps the largest connected piece.
///
/// A diamond site has its first-shell neighbours at `(+-1, +-1, +-1)` on the
/// quarter-cell grid. A carbon with fewer than [`MIN_CARBON_NEIGHBOURS`] is
/// removed, again and again, until none is left. Then only the largest
/// connected set stays, so the part is one piece of diamond.
pub(crate) fn keep_well_bonded(sites_m: Vec<[f64; 3]>) -> Vec<[f64; 3]> {
    let mut alive: HashMap<(i64, i64, i64), [f64; 3]> = sites_m
        .iter()
        .map(|site| (site_key(*site), *site))
        .collect();
    let neighbours = |key: (i64, i64, i64), alive: &HashMap<(i64, i64, i64), [f64; 3]>| {
        let mut found = Vec::with_capacity(4);
        for dx in [-1, 1] {
            for dy in [-1, 1] {
                for dz in [-1, 1] {
                    let other = (key.0 + dx, key.1 + dy, key.2 + dz);
                    if alive.contains_key(&other) {
                        found.push(other);
                    }
                }
            }
        }
        found
    };

    let mut queue: VecDeque<(i64, i64, i64)> = alive.keys().copied().collect();
    while let Some(key) = queue.pop_front() {
        if !alive.contains_key(&key) {
            continue;
        }
        let around = neighbours(key, &alive);
        if around.len() < MIN_CARBON_NEIGHBOURS {
            alive.remove(&key);
            queue.extend(around);
        }
    }

    let mut seen: HashSet<(i64, i64, i64)> = HashSet::new();
    let mut best: Vec<(i64, i64, i64)> = Vec::new();
    for &start in alive.keys() {
        if seen.contains(&start) {
            continue;
        }
        let mut component = Vec::new();
        let mut stack = vec![start];
        seen.insert(start);
        while let Some(key) = stack.pop() {
            component.push(key);
            for other in neighbours(key, &alive) {
                if seen.insert(other) {
                    stack.push(other);
                }
            }
        }
        if component.len() > best.len() {
            best = component;
        }
    }
    let keep: HashSet<(i64, i64, i64)> = best.into_iter().collect();
    // Keep the input order, so the build is deterministic.
    sites_m
        .into_iter()
        .filter(|site| keep.contains(&site_key(*site)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_set() -> GearPumpSet {
        GearPumpGenerator
            .build(&ParameterSet::new())
            .expect("the default pump builds")
    }

    fn carbons(part: &Part, range: &Range<usize>) -> Vec<[f64; 3]> {
        range
            .clone()
            .filter(|&index| part.topology.element(index) == Some(Element::CARBON))
            .filter_map(|index| part.topology.position_m(index))
            .collect()
    }

    fn rotate_about(point: [f64; 3], center: [f64; 2], angle_rad: f64) -> [f64; 3] {
        let (sin, cos) = angle_rad.sin_cos();
        let dx = point[0] - center[0];
        let dy = point[1] - center[1];
        [
            center[0] + dx * cos - dy * sin,
            center[1] + dx * sin + dy * cos,
            point[2],
        ]
    }

    fn min_distance_m(left: &[[f64; 3]], right: &[[f64; 3]]) -> f64 {
        let mut best = f64::INFINITY;
        for a in left {
            for b in right {
                let d =
                    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt();
                best = best.min(d);
            }
        }
        best
    }

    #[test]
    fn the_default_design_is_sound() {
        let set = default_set();
        let design = &set.design;
        assert_eq!(design.teeth(), 8);
        assert!(design.tip_to_root_m() >= design.clearance_m());
        assert!(design.post_radius_m() > 0.0);
        assert!(set.housing_atoms.end == set.gear_a_atoms.start);
        assert!(set.gear_a_atoms.end == set.gear_b_atoms.start);
        assert_eq!(set.gear_b_atoms.end, set.part.atom_count());
        assert_eq!(
            set.gear_a_atoms.len(),
            set.gear_b_atoms.len(),
            "the gears are equal"
        );
    }

    #[test]
    fn the_gears_clear_each_other_through_one_tooth_pitch() {
        let set = default_set();
        let design = &set.design;
        let a = carbons(&set.part, &set.gear_a_atoms);
        let b = carbons(&set.part, &set.gear_b_atoms);
        let pitch_rad = 2.0 * PI / design.teeth() as f64;
        let steps = 12;
        let mut worst = f64::INFINITY;
        for step in 0..steps {
            let turn_rad = pitch_rad * step as f64 / steps as f64;
            let a_now: Vec<_> = a
                .iter()
                .map(|p| rotate_about(*p, design.gear_a_center_m(), turn_rad))
                .collect();
            // Gear B turns the other way at the same rate.
            let b_now: Vec<_> = b
                .iter()
                .map(|p| rotate_about(*p, design.gear_b_center_m(), -turn_rad))
                .collect();
            worst = worst.min(min_distance_m(&a_now, &b_now));
        }
        // The defaults were tuned to a 4.5 Angstrom worst gap, near the
        // 4.6 Angstrom housing clearance. Two capped faces need about
        // 3.8 Angstrom before the caps overlap at the 1.6 Angstrom threshold.
        assert!(worst > 4.2e-10, "worst gear-to-gear gap {worst:e} m");
    }

    #[test]
    fn the_gears_clear_the_housing() {
        let set = default_set();
        let housing = carbons(&set.part, &set.housing_atoms);
        for range in [&set.gear_a_atoms, &set.gear_b_atoms] {
            let gear = carbons(&set.part, range);
            let gap = min_distance_m(&housing, &gear);
            assert!(gap > 3.8e-10, "gear-to-housing gap {gap:e} m");
        }
    }

    #[test]
    fn the_ports_open_at_the_mesh_waist() {
        let set = default_set();
        let housing = carbons(&set.part, &set.housing_atoms);
        let half_width_m = 0.5 * set.design.port_width_m() - DIAMOND_LATTICE_CONSTANT_M;
        let half_height_m = set.design.cavity_half_height_m() - DIAMOND_LATTICE_CONSTANT_M;
        let blocked = housing.iter().any(|p| {
            p[0].abs() < half_width_m
                && p[2].abs() < half_height_m
                && p[1].abs() > set.design.pitch_radius_m()
        });
        assert!(!blocked, "a housing carbon sits in a port");
    }

    #[test]
    fn loose_carbons_are_removed() {
        let a = DIAMOND_LATTICE_CONSTANT_M;
        let quarter = 0.25 * a;
        let mut sites = crate::lattice_fill::fill_solid(&Box3 {
            min_m: [0.0, 0.0, 0.0],
            max_m: [a, a, a],
        });
        let before = sites.len();
        // A far island of one site must go.
        sites.push([10.0 * a, 10.0 * a, 10.0 * a]);
        let kept = keep_well_bonded(sites);
        assert!(kept.len() <= before);
        assert!(kept.iter().all(|p| p[0] < 2.0 * a + quarter));
    }

    #[test]
    fn every_tooth_gap_has_one_root_iron() {
        let set = default_set();
        let design = &set.design;
        let teeth = design.teeth();
        assert_eq!(set.metal_sites.len(), 2 * teeth);
        let pitch_rad = 2.0 * PI / teeth as f64;
        for (n, &site) in set.metal_sites.iter().enumerate() {
            assert_eq!(set.part.topology.element(site), Some(Element::IRON));
            let (center, phase) = if n < teeth {
                (design.gear_a_center_m(), 0.0)
            } else {
                (design.gear_b_center_m(), design.gear_b_phase_rad())
            };
            let p = set.part.topology.position_m(site).expect("iron");
            let radius_m = (p[0] - center[0]).hypot(p[1] - center[1]);
            assert!(radius_m > design.root_radius_m(), "{radius_m}");
            assert!(radius_m < design.root_radius_m() + 3.0e-10, "{radius_m}");
            // A gap centre sits half a pitch from a tooth tip.
            let angle = (p[1] - center[1]).atan2(p[0] - center[0]) - phase;
            let from_tip = angle.rem_euclid(pitch_rad);
            assert!(
                (from_tip - 0.5 * pitch_rad).abs() < 0.25 * pitch_rad,
                "iron {n} sits {from_tip} rad from a tip"
            );
        }
    }

    #[test]
    fn the_iron_clears_the_other_gear_through_the_mesh() {
        let set = default_set();
        let design = &set.design;
        let teeth = design.teeth();
        let irons_a: Vec<[f64; 3]> = set.metal_sites[..teeth]
            .iter()
            .filter_map(|&i| set.part.topology.position_m(i))
            .collect();
        let b = carbons(&set.part, &set.gear_b_atoms);
        let pitch_rad = 2.0 * PI / teeth as f64;
        let mut worst = f64::INFINITY;
        for step in 0..48 {
            let turn_rad = 2.0 * pitch_rad * step as f64 / 48.0;
            let a_now: Vec<_> = irons_a
                .iter()
                .map(|p| rotate_about(*p, design.gear_a_center_m(), turn_rad))
                .collect();
            let b_now: Vec<_> = b
                .iter()
                .map(|p| rotate_about(*p, design.gear_b_center_m(), -turn_rad))
                .collect();
            worst = worst.min(min_distance_m(&a_now, &b_now));
        }
        // Measured 5.46e-10 m. A bound O2 reaches only about 3.1e-10 m from
        // the iron, so the mating tooth never touches it: the mesh does not
        // release the guest (ADR-0089, ADR-0090). The test only checks that
        // the iron does not collide with the other gear.
        assert!(worst > 3.0e-10, "iron to mating carbon {worst:e} m");
    }

    #[test]
    fn the_generator_is_registered() {
        assert!(crate::registry::generator("gear_pump").is_some());
    }
}
