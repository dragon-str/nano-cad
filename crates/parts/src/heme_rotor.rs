//! The oxygen-loading rotor of the respirocyte design (ADR-0111).
//!
//! A diamond disc carries `pocket_count` radial pockets on its rim. Each
//! pocket holds one heme-like site: an idealized Fe-porphine whose plane is
//! normal to the radius, with the iron on the pocket axis. A snug bore runs
//! from the ring face out to the rim; O2 binds on the outer face of the iron.
//! Below the ring, a radial channel carries an arm: a diamond rod whose tip is
//! an NH2 group, the axial ligand. Tip in at 2.1 A from the iron loads O2, tip
//! out at 4.0 A releases it (ADR-0106).
//!
//! The arm cam (ADR-0116). Each arm carries a follower pin along z: up on even
//! pockets, down on odd ones, so the slots in the disc skins stay apart. The
//! pin passes through a radial slot in the skin and runs in a groove in the
//! housing end plate ([`ArmCamTrack`]). The groove centreline steps 1.9 A in
//! radius, so the fixed housing pushes each arm in on the loading arc and pulls
//! it out past the outlet. Both groove walls stand at H-H contact from the pin
//! ([`CAM_CONTACT_M`]), so the cam drives the arm both ways with no designed
//! play.
//!
//! Every body is returned in world coordinates. The ring is held in its slot
//! by position only: the chemical joint of the porphine to diamond is not
//! designed, so the scene holds it with a fixed joint (a placeholder).
//!
//! The site geometry is the rigid model of the quantum-chemistry scans
//! (Fe-N 2.00 A, O2 end-on at 1.85 A, 60 degrees). The relaxation keeps it.

use std::f64::consts::PI;

use nanocad_model::{Atom, Bond, BondType, Element, Part, Topology};
use nanocad_units::Unit;

use crate::diamond_solid;
use crate::error::PartError;
use crate::gear_pump::keep_well_bonded;
use crate::generator::PartGenerator;
use crate::lattice_fill::fill_solid;
use crate::parameter::{ParameterSet, ParameterSpec};
use crate::shape::{Bounds, Box3, Cylinder, Difference, Placed, RadialCylinder, Solid, Union};

const A: f64 = 1.0e-10;

static CARRIER_SITE_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new(
        "pocket_radius_m",
        Some(Unit::Metre),
        3.6e-10,
        3.2e-10,
        6.0e-10,
        false,
        "radius of the O2 pocket above the iron",
    ),
    ParameterSpec::new(
        "mouth_radius_m",
        Some(Unit::Metre),
        2.6e-10,
        2.2e-10,
        4.0e-10,
        false,
        "radius of the mouth in the top face",
    ),
    ParameterSpec::new(
        "slot_radius_m",
        Some(Unit::Metre),
        7.0e-10,
        6.0e-10,
        1.0e-9,
        false,
        "radius of the ring slot that holds the porphine",
    ),
    ParameterSpec::new(
        "guard_azimuth_rad",
        None,
        0.0,
        -3.2,
        3.2,
        false,
        "azimuth of the distal hydroxyl guard",
    ),
];

static HEME_ROTOR_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new("disc_radius_m", Some(Unit::Metre), 4.0e-9, 2.5e-9, 2.0e-8, false, "rim radius of the rotor disc"),
    ParameterSpec::new("thickness_m", Some(Unit::Metre), 2.4e-9, 1.6e-9, 1.0e-8, false, "axial thickness of the disc"),
    ParameterSpec::new("pocket_count", None, 8.0, 1.0, 64.0, true, "number of heme pockets on the rim; at 12 the ring slots of neighbours overlap (ADR-0121)"),
    ParameterSpec::new("bore_radius_m", Some(Unit::Metre), 4.0e-10, 4.0e-10, 2.0e-9, false, "radius of the bore from the ring face to the rim (carbon centres)"),
    ParameterSpec::new("bore_length_m", Some(Unit::Metre), 3.0e-10, 0.0, 2.0e-9, false, "length of the bore above the ring slot"),
    ParameterSpec::new("slot_radius_m", Some(Unit::Metre), 9.0e-10, 7.0e-10, 2.0e-9, false, "radius of the slot that holds the ring"),
    ParameterSpec::new("slot_half_m", Some(Unit::Metre), 3.6e-10, 3.0e-10, 1.0e-9, false, "half-thickness of the ring slot along the radius"),
    ParameterSpec::new("arm_radius_m", Some(Unit::Metre), 2.6e-10, 1.5e-10, 1.0e-9, false, "radius of the ligand arm lattice"),
    ParameterSpec::new("arm_clearance_m", Some(Unit::Metre), 5.5e-10, 3.0e-10, 1.0e-9, false, "carbon gap between the arm and its channel"),
    ParameterSpec::new("arm_length_m", Some(Unit::Metre), 1.4e-9, 5.0e-10, 3.0e-9, false, "length of the ligand arm"),
    ParameterSpec::new("pocket_phase_rad", None, 0.261_799_387_799_149_4, -3.2, 3.2, false, "azimuth of pocket 0; the cam profile then sets each pocket's state (ligand in and O2 bound, or ligand out and empty)"),
    ParameterSpec::new("cam_pin_offset_m", Some(Unit::Metre), 3.0e-10, 2.6e-10, 1.0e-9, false, "distance of the cam pin axis from the back end of the arm"),
    ParameterSpec::new("cam_pin_top_m", Some(Unit::Metre), 2.5e-9, 1.4e-9, 5.0e-9, false, "height of the cam pin end above the disc mid-plane"),
    ParameterSpec::new("arm_collar_gap_m", Some(Unit::Metre), 4.6e-10, 3.8e-10, 1.0e-9, false, "carbon gap between the arm and the collar that closes the arm channel from the ring slot"),
    ParameterSpec::new("arm_collar_length_m", Some(Unit::Metre), 4.0e-10, 2.0e-10, 1.5e-9, false, "length of the arm collar below the ring slot"),
    ParameterSpec::new("cam_release_margin_m", Some(Unit::Metre), 2.0e-10, 0.0, 1.0e-9, false, "extra withdrawal of the release track beyond the 4.0 A design, taken up by the lattice steps of the groove wall"),
    ParameterSpec::new("repair_donor", None, 1.0, 0.0, 1.0, true, "1 builds the iron-repair donor on the collar wall of each pocket (ADR-0108, ADR-0141)"),
    ParameterSpec::new("catalyst", None, 1.0, 0.0, 1.0, true, "1 builds the dismutase/catalase manganese centre on each pocket wall (ADR-0139, ADR-0141)"),
    ParameterSpec::new("iron_dome_m", Some(Unit::Metre), 0.0, 0.0, 1.0e-9, false, "out-of-plane displacement of the iron in an unloaded pocket (five-coordinate dome), in metres (M28-34)"),
];

/// The longest iron-to-donor distance the rotor accepts, in metres. The arm
/// collar occupies the proximal axis, so the nearest scaffold carbon sits on
/// the collar wall, further than the 3.6 A of the free carrier site (ADR-0141).
const ROTOR_DONOR_REACH_M: f64 = 9.2e-10;

/// The carbon-centre gap between the cam pin and a groove wall: the two
/// hydrogen caps at van der Waals contact (1.09 + 2.4 + 1.09 A, rounded up to
/// the running clearance of the module).
pub const CAM_CONTACT_M: f64 = 4.6e-10;
/// The azimuth where the cam pushes the arm in, in radians: inside the seal
/// the pocket crosses from the tank back to the blood (ADR-0117).
pub const CAM_INLET_RAD: f64 = 270.0 * PI / 180.0;
/// The azimuth where the cam pulls the arm out, in radians: just inside the
/// tank, after the loaded pocket has crossed the seal from the blood.
pub const CAM_OUTLET_RAD: f64 = 120.0 * PI / 180.0;
/// The half-width of each cam ramp, in radians.
pub const CAM_RAMP_HALF_RAD: f64 = 12.0 * PI / 180.0;

/// The groove of the arm cam: an annulus whose centreline radius follows the
/// cam profile, between two heights. It is a cut in the housing end plate.
#[derive(Clone, Debug)]
pub struct ArmCamTrack {
    /// The centreline radius where the arm is in (loaded), in metres.
    pub loaded_radius_m: f64,
    /// The centreline radius where the arm is out (release), in metres.
    pub released_radius_m: f64,
    /// The groove wall inside the centreline, measured from it to the wall
    /// carbon centres, in metres.
    pub inner_width_m: f64,
    /// The groove wall outside the centreline, in metres.
    pub outer_width_m: f64,
    /// The lower height of the groove, in metres.
    pub z_min_m: f64,
    /// The upper height of the groove, in metres.
    pub z_max_m: f64,
    /// The azimuth where the follower goes in, in radians.
    pub inlet_rad: f64,
    /// The azimuth where the follower goes out, in radians.
    pub outlet_rad: f64,
    /// The half-width of each ramp, in radians.
    pub ramp_half_rad: f64,
}

impl ArmCamTrack {
    /// The fraction of the stroke pushed in at azimuth `theta_rad`: 1 on the
    /// loading arc from the inlet to the outlet (rotor turning +z), 0 on the
    /// return arc, with a half-cosine ramp across each end.
    pub fn loaded_fraction(theta_rad: f64) -> f64 {
        Self::fraction_between(theta_rad, CAM_INLET_RAD, CAM_OUTLET_RAD)
    }

    /// As [`ArmCamTrack::loaded_fraction`], for any inlet and outlet.
    pub fn fraction_between(theta_rad: f64, inlet_rad: f64, outlet_rad: f64) -> f64 {
        Self::fraction_with(theta_rad, inlet_rad, outlet_rad, CAM_RAMP_HALF_RAD)
    }

    /// As [`ArmCamTrack::fraction_between`], with a ramp half-width.
    pub fn fraction_with(
        theta_rad: f64,
        inlet_rad: f64,
        outlet_rad: f64,
        ramp_half_rad: f64,
    ) -> f64 {
        let wrap = |a: f64| (a + PI).rem_euclid(2.0 * PI) - PI;
        let d_in = wrap(theta_rad - inlet_rad);
        if d_in.abs() < ramp_half_rad {
            return 0.5 - 0.5 * (PI * (d_in + ramp_half_rad) / (2.0 * ramp_half_rad)).cos();
        }
        let d_out = wrap(theta_rad - outlet_rad);
        if d_out.abs() < ramp_half_rad {
            return 0.5 + 0.5 * (PI * (d_out + ramp_half_rad) / (2.0 * ramp_half_rad)).cos();
        }
        let from_inlet = (theta_rad - inlet_rad).rem_euclid(2.0 * PI);
        let span = (outlet_rad - inlet_rad).rem_euclid(2.0 * PI);
        if from_inlet < span {
            1.0
        } else {
            0.0
        }
    }

    /// The centreline radius at azimuth `theta_rad`, in metres.
    pub fn centre_radius_m(&self, theta_rad: f64) -> f64 {
        let f = Self::fraction_with(
            theta_rad,
            self.inlet_rad,
            self.outlet_rad,
            self.ramp_half_rad,
        );
        self.released_radius_m + f * (self.loaded_radius_m - self.released_radius_m)
    }
}

impl Solid for ArmCamTrack {
    fn contains_m(&self, p: [f64; 3]) -> bool {
        if p[2] < self.z_min_m || p[2] > self.z_max_m {
            return false;
        }
        // The distance across the groove is measured normal to its
        // centreline: on a ramp the radial width grows by sqrt(1 + slope^2),
        // so a round pin fits the ramp as well as the dwell (ADR-0122).
        let r = p[0].hypot(p[1]);
        let theta = p[1].atan2(p[0]);
        let h = 1.0e-3;
        let slope = (self.centre_radius_m(theta + h) - self.centre_radius_m(theta - h))
            / (2.0 * h * r.max(1.0e-12));
        let d = (r - self.centre_radius_m(theta)) / (1.0 + slope * slope).sqrt();
        -self.inner_width_m <= d && d <= self.outer_width_m
    }

    fn bounds_m(&self) -> Bounds {
        let r = self.loaded_radius_m.max(self.released_radius_m) + self.outer_width_m;
        Bounds {
            min_m: [-r, -r, self.z_min_m],
            max_m: [r, r, self.z_max_m],
        }
    }

    fn clone_box(&self) -> Box<dyn Solid> {
        Box::new(self.clone())
    }
}

/// The ligand tip distance from the iron in the loaded state, in metres.
pub const LIGAND_IN_M: f64 = 2.1e-10;
/// The ligand tip distance from the iron in the release state, in metres.
pub const LIGAND_OUT_M: f64 = 4.0e-10;
/// The out-of-plane displacement of the iron in an unloaded pocket. The armed
/// empty site is five-coordinate, so the high-spin iron domes toward the axial
/// nitrogen (M28-34). The six-coordinate loaded site stays in the ring plane.
pub const IRON_DOME_M: f64 = 0.3e-10;
const FE_O_M: f64 = 1.85e-10;
const O_O_M: f64 = 1.23e-10;
const N_C_M: f64 = 1.47e-10;
const N_H_M: f64 = 1.01e-10;

/// A generated rotor with its sites and arms, all in world coordinates.
#[derive(Clone, Debug)]
pub struct HemeRotorSet {
    /// The diamond disc.
    pub rotor: Part,
    /// One heme-like site per pocket: Fe-porphine, plus O2 when loaded.
    pub sites: Vec<Part>,
    /// One ligand arm per pocket.
    pub arms: Vec<Part>,
    /// The iron atom index inside each site part.
    pub site_iron: Vec<usize>,
    /// The NH2 nitrogen index inside each arm part.
    pub arm_tip: Vec<usize>,
    /// For each site, the (site atom, rotor atom) pairs of the covalent
    /// attachment tethers. The site atom is a tether end carbon.
    pub site_anchor: Vec<Vec<(usize, usize)>>,
    /// For each site, its tethers (ADR-0137). The rotor atom indexes are those
    /// of `rotor`, after the replaced hydrogens were removed.
    pub site_tethers: Vec<Vec<Tether>>,
    /// For each site, the iron-repair donor (ADR-0108), when built. The rotor
    /// atom indexes are those of `rotor`, after the removed atoms were pruned.
    pub site_donor: Vec<Option<DonorSite>>,
    /// For each site, the dismutase/catalase catalyst (ADR-0139), when built.
    pub site_catalyst: Vec<Option<CatalystSite>>,
    /// The azimuth of each pocket, in radians.
    pub azimuth_rad: Vec<f64>,
    /// Whether each pocket shows the loaded state.
    pub loaded: Vec<bool>,
    /// The out-of-plane iron displacement shown in each pocket, in metres: the
    /// five-coordinate dome (M28-34) in an empty pocket, zero when loaded.
    pub site_dome_m: Vec<f64>,
    /// The radius of the iron atoms, in metres.
    pub iron_radius_m: f64,
    /// The disc rim radius and half-thickness, in metres.
    pub disc_radius_m: f64,
    /// The disc half-thickness, in metres.
    pub half_thickness_m: f64,
    /// The side of each cam pin: +1 up (+z), -1 down.
    pub pin_side: Vec<f64>,
    /// The lattice radius of the cam pin, in metres.
    pub pin_radius_m: f64,
    /// The radius of the pin axis with the arm in, in metres.
    pub pin_loaded_radius_m: f64,
    /// The radius of the pin axis on the release track, in metres. It is
    /// withdrawn by the release margin beyond the 4.0 A design.
    pub pin_released_radius_m: f64,
    /// The release margin, in metres.
    pub release_margin_m: f64,
    /// The largest |z| of a pin carbon, in metres.
    pub pin_top_m: f64,
    /// How far the pin carbons above the disc face reach inward of the pin
    /// axis (toward the rotor axis), in metres.
    pub pin_inner_m: f64,
    /// How far they reach outward of the pin axis, in metres.
    pub pin_outer_m: f64,
    /// The azimuth where the cam pushes each follower in, in radians.
    pub cam_inlet_rad: f64,
    /// The azimuth where the cam pulls each follower out, in radians.
    pub cam_outlet_rad: f64,
    /// The half-width of each cam ramp, in radians.
    pub cam_ramp_half_rad: f64,
}

impl HemeRotorSet {
    /// The two cam grooves (top, then bottom) for a housing whose chamber
    /// ends at |z| = `chamber_half_m`. The groove floor stands at
    /// [`CAM_CONTACT_M`] beyond the pin end.
    pub fn cam_tracks(&self, chamber_half_m: f64) -> [ArmCamTrack; 2] {
        let floor = self.pin_top_m + CAM_CONTACT_M;
        let top = ArmCamTrack {
            loaded_radius_m: self.pin_loaded_radius_m,
            released_radius_m: self.pin_released_radius_m,
            inner_width_m: self.pin_inner_m + CAM_CONTACT_M,
            outer_width_m: self.pin_outer_m + CAM_CONTACT_M,
            z_min_m: chamber_half_m - 1.0e-10,
            z_max_m: floor,
            inlet_rad: self.cam_inlet_rad,
            outlet_rad: self.cam_outlet_rad,
            ramp_half_rad: self.cam_ramp_half_rad,
        };
        let bottom = ArmCamTrack {
            z_min_m: -floor,
            z_max_m: -top.z_min_m,
            ..top.clone()
        };
        [top, bottom]
    }
}

/// Builds the oxygen-loading rotor with its heme sites and ligand arms.
#[derive(Clone, Copy, Debug, Default)]
pub struct HemeRotorGenerator;

/// The idealized porphine (Fe-N 2.00 A) in the ring plane, (a, b) in metres.
fn porphine_plane() -> Vec<(Element, [f64; 2])> {
    let quarter: [(Element, [f64; 2]); 9] = [
        (Element::NITROGEN, [2.00, 0.00]),
        (Element::CARBON, [2.83, 1.10]),
        (Element::CARBON, [2.83, -1.10]),
        (Element::CARBON, [4.23, 0.675]),
        (Element::CARBON, [4.23, -0.675]),
        (Element::HYDROGEN, [5.13, 1.27]),
        (Element::HYDROGEN, [5.13, -1.27]),
        (Element::CARBON, [2.43, 2.43]),
        (Element::HYDROGEN, [3.19, 3.19]),
    ];
    let mut out = Vec::new();
    for k in 0..4 {
        let (s, c) = (k as f64 * PI / 2.0).sin_cos();
        for (e, p) in quarter {
            out.push((e, [(c * p[0] - s * p[1]) * A, (s * p[0] + c * p[1]) * A]));
        }
    }
    out
}

fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

fn normalize(v: [f64; 3]) -> [f64; 3] {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if n == 0.0 {
        [0.0; 3]
    } else {
        [v[0] / n, v[1] / n, v[2] / n]
    }
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Bonds every pair closer than the covalent sum times 1.2.
fn bond_by_distance(topology: &mut Topology) -> Result<(), PartError> {
    let n = topology.atom_count();
    for i in 0..n {
        for j in (i + 1)..n {
            let (Some(p), Some(q)) = (topology.position_m(i), topology.position_m(j)) else {
                continue;
            };
            let (Some(ei), Some(ej)) = (topology.element(i), topology.element(j)) else {
                continue;
            };
            let limit = match (ei.atomic_number(), ej.atomic_number()) {
                (1, 1) => 0.0,
                (26, 7) | (7, 26) => 2.1e-10,
                (26, 8) | (8, 26) => 2.0e-10,
                (1, _) | (_, 1) => 1.25e-10,
                _ => 1.6e-10,
            };
            if dist(p, q) < limit {
                topology.add_bond(Bond::new(i as u32, j as u32, 1, BondType::Single))?;
            }
        }
    }
    Ok(())
}

impl HemeRotorGenerator {
    /// Resolves the inputs and builds every body.
    pub fn build(&self, parameters: &ParameterSet) -> Result<HemeRotorSet, PartError> {
        let r = self.resolve(parameters)?;
        let disc_m = r.require("disc_radius_m")?;
        let thick_m = r.require("thickness_m")?;
        let count = r.require("pocket_count")?.round() as usize;
        let bore_m = r.require("bore_radius_m")?;
        let bore_len_m = r.require("bore_length_m")?;
        let slot_r_m = r.require("slot_radius_m")?;
        let slot_h_m = r.require("slot_half_m")?;
        let arm_r_m = r.require("arm_radius_m")?;
        let arm_clr_m = r.require("arm_clearance_m")?;
        let arm_len_m = r.require("arm_length_m")?;
        let phase_rad = r.require("pocket_phase_rad")?;
        let release_margin_m = r.require("cam_release_margin_m")?;
        let build_donor = r.require("repair_donor")?.round() >= 0.5;
        let build_catalyst = r.require("catalyst")?.round() >= 0.5;
        let iron_dome_m = r.require("iron_dome_m")?;
        let pin_offset_m = r.require("cam_pin_offset_m")?;
        let pin_top_m = r.require("cam_pin_top_m")?;
        let half_m = 0.5 * thick_m;
        if pin_top_m < half_m + 2.0e-10 {
            return Err(PartError::InvalidGeometry(
                "the cam pin does not reach past the disc face".to_owned(),
            ));
        }
        if slot_r_m > half_m {
            return Err(PartError::InvalidGeometry(
                "the ring slot is wider than the disc is thick".to_owned(),
            ));
        }
        // The iron plane sits so that the slot and the bore end at the rim.
        let iron_r_m = disc_m - slot_h_m - bore_len_m;
        let arm_chan_r_m = arm_r_m + arm_clr_m;
        let chan_len_m = LIGAND_OUT_M + release_margin_m + arm_len_m + 1.0e-9;
        let collar_gap_m = r.require("arm_collar_gap_m")?;
        let collar_len_m = r.require("arm_collar_length_m")?;
        // The collar runs from the slot floor down `collar_len_m`.
        let collar_lo_m = iron_r_m - slot_h_m - collar_len_m;

        // The disc and its cuts.
        let mut cuts: Option<Box<dyn Solid>> = None;
        let azimuths: Vec<f64> = (0..count)
            .map(|k| (phase_rad + 2.0 * PI * k as f64 / count as f64).rem_euclid(2.0 * PI))
            .collect();
        for &az in &azimuths {
            let outer = disc_m + 1.0e-9;
            let pieces: [Box<dyn Solid>; 4] = [
                Box::new(RadialCylinder {
                    azimuth_rad: az,
                    center_radius_m: 0.5 * (iron_r_m + outer),
                    radius_m: bore_m,
                    half_length_m: 0.5 * (outer - iron_r_m),
                    center_z_m: 0.0,
                }),
                Box::new(RadialCylinder {
                    azimuth_rad: az,
                    center_radius_m: iron_r_m,
                    radius_m: slot_r_m,
                    half_length_m: slot_h_m,
                    center_z_m: 0.0,
                }),
                // The arm channel, behind the collar.
                Box::new(RadialCylinder {
                    azimuth_rad: az,
                    center_radius_m: 0.5 * (collar_lo_m + iron_r_m - chan_len_m),
                    radius_m: arm_chan_r_m,
                    half_length_m: 0.5 * (collar_lo_m - (iron_r_m - chan_len_m)),
                    center_z_m: 0.0,
                }),
                // The collar: the arm passes into the ring slot through a
                // contact-gap bore, so the pocket does not open into the
                // rotor interior (ADR-0121).
                Box::new(RadialCylinder {
                    azimuth_rad: az,
                    center_radius_m: 0.5 * (collar_lo_m + iron_r_m),
                    radius_m: arm_r_m + collar_gap_m,
                    half_length_m: 0.5 * (iron_r_m - collar_lo_m),
                    center_z_m: 0.0,
                }),
            ];
            for piece in pieces {
                cuts = Some(match cuts {
                    None => piece,
                    Some(prev) => Box::new(Union { a: prev, b: piece }),
                });
            }
        }
        let disc = Cylinder {
            center_m: [0.0, 0.0, 0.0],
            radius_m: disc_m,
            height_m: thick_m,
        };
        let solid: Box<dyn Solid> = match cuts {
            Some(inner) => Box::new(Difference {
                outer: Box::new(disc),
                inner,
            }),
            None => Box::new(disc),
        };

        // The arm is a diamond rod along the crystal [111] axis: a carbon on
        // its (111) end face has one bond straight along the axis, so the
        // C-N bond and the nitrogen lone pair face the iron. The rod is then
        // turned so [111] becomes +x.
        let a111 = [1.0 / 3f64.sqrt(); 3];
        // The pin runs along [1-10], a bond-chain direction normal to [111];
        // the frame turns [111] to +x and [1-10] to +z.
        let e3 = [1.0 / 2f64.sqrt(), -1.0 / 2f64.sqrt(), 0.0];
        let turn = [cross(e3, a111), e3];
        let turn = [a111, turn[0], turn[1]];
        let arm_solid = Union {
            a: Box::new(AxisRod {
                origin: [0.0; 3],
                axis: a111,
                radius_m: arm_r_m,
                length_m: arm_len_m,
            }),
            b: Box::new(AxisRod {
                origin: [
                    pin_offset_m * a111[0],
                    pin_offset_m * a111[1],
                    pin_offset_m * a111[2],
                ],
                axis: e3,
                radius_m: arm_r_m,
                length_m: pin_top_m,
            }),
        };
        let mut arm = Topology::new();
        let arm_c: Vec<u32> = keep_well_bonded(fill_solid(&arm_solid))
            .into_iter()
            .map(|p| arm.add_atom(Atom::new(Element::CARBON, p, 0.0, "C")))
            .collect();
        diamond_solid::bond_and_cap(&mut arm, &arm_c)?;
        for i in 0..arm.atom_count() {
            if let Some(p) = arm.position_m(i) {
                arm.set_position_m(i, apply3(&turn, p))?;
            }
        }
        // The tip: the cap most along +x among caps pointing within 10
        // degrees of +x, taken on the carbon nearest the axis.
        let mut best: Option<(f64, usize, usize)> = None;
        for bond in arm.bonds() {
            let (c_i, h_i) = match (arm.element(bond.u as usize), arm.element(bond.v as usize)) {
                (Some(Element::CARBON), Some(Element::HYDROGEN)) => {
                    (bond.u as usize, bond.v as usize)
                }
                (Some(Element::HYDROGEN), Some(Element::CARBON)) => {
                    (bond.v as usize, bond.u as usize)
                }
                _ => continue,
            };
            let (Some(c), Some(h)) = (arm.position_m(c_i), arm.position_m(h_i)) else {
                continue;
            };
            let l = dist(c, h);
            if (h[0] - c[0]) / l < 0.985 {
                continue;
            }
            let score = c[0] - 0.3 * c[1].hypot(c[2]);
            if best.is_none_or(|(b, _, _)| score > b) {
                best = Some((score, c_i, h_i));
            }
        }
        let (_, tip_c, cap) = best.ok_or_else(|| {
            PartError::InvalidGeometry("the arm has no axial tip bond".to_owned())
        })?;
        let tip_m = arm.position_m(tip_c).unwrap_or([0.0; 3]);
        let cap_m = arm.position_m(cap).unwrap_or([0.0; 3]);
        let len = dist(cap_m, tip_m);
        let d = [
            (cap_m[0] - tip_m[0]) / len,
            (cap_m[1] - tip_m[1]) / len,
            (cap_m[2] - tip_m[2]) / len,
        ];
        let n_m = [
            tip_m[0] + N_C_M * d[0],
            tip_m[1] + N_C_M * d[1],
            tip_m[2] + N_C_M * d[2],
        ];
        arm.set_element(cap, Element::NITROGEN)?;
        arm.set_atom_type(cap, "N_amine")?;
        arm.set_position_m(cap, n_m)?;
        // Two hydrogens on the nitrogen, pointing back from the iron: the
        // lone pair faces the iron along +d, so each H makes 109.5 degrees
        // with it.
        let helper = if d[2].abs() < 0.9 {
            [0.0, 0.0, 1.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        let mut w = [
            d[1] * helper[2] - d[2] * helper[1],
            d[2] * helper[0] - d[0] * helper[2],
            d[0] * helper[1] - d[1] * helper[0],
        ];
        let wn = (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt();
        w = [w[0] / wn, w[1] / wn, w[2] / wn];
        let (sa, ca) = (70.5_f64.to_radians().sin(), 70.5_f64.to_radians().cos());
        for sign in [1.0, -1.0] {
            let h = [
                n_m[0] + N_H_M * (-ca * d[0] + sign * sa * w[0]),
                n_m[1] + N_H_M * (-ca * d[1] + sign * sa * w[1]),
                n_m[2] + N_H_M * (-ca * d[2] + sign * sa * w[2]),
            ];
            let hi = arm.add_atom(Atom::new(Element::HYDROGEN, h, 0.0, "H"));
            arm.add_bond(Bond::new(cap as u32, hi, 1, BondType::Single))?;
        }
        // Shift so the nitrogen sits on the x axis at x = 0.
        let shift = [-n_m[0], -n_m[1], -n_m[2]];
        for i in 0..arm.atom_count() {
            if let Some(p) = arm.position_m(i) {
                arm.set_position_m(i, [p[0] + shift[0], p[1] + shift[1], p[2] + shift[2]])?;
            }
        }
        // The pin axis in the arm frame, and the measured pin end.
        let pin_x_m = pin_offset_m + shift[0];
        let pin_y_m = shift[1];
        let arm_pin_top_m = (0..arm.atom_count())
            .filter(|&i| arm.element(i) == Some(Element::CARBON))
            .filter_map(|i| arm.position_m(i))
            .map(|p| p[2])
            .fold(f64::NEG_INFINITY, f64::max);
        // The pin carbons above the disc face set the groove walls.
        let (mut pin_inner_m, mut pin_outer_m) = (0.0_f64, 0.0_f64);
        for i in 0..arm.atom_count() {
            if arm.element(i) != Some(Element::CARBON) {
                continue;
            }
            if let Some(p) = arm.position_m(i) {
                if p[2] > half_m {
                    pin_inner_m = pin_inner_m.max(pin_x_m - p[0]);
                    pin_outer_m = pin_outer_m.max(p[0] - pin_x_m);
                }
            }
        }
        let pin_track = |ligand_m: f64| (iron_r_m - ligand_m + pin_x_m).hypot(pin_y_m);
        let pin_loaded_r_m = pin_track(LIGAND_IN_M);
        let pin_released_r_m = pin_track(LIGAND_OUT_M + release_margin_m);
        let slot_w_m = arm_r_m + arm_clr_m;
        if pin_loaded_r_m.max(pin_released_r_m) + slot_w_m > iron_r_m - slot_h_m {
            return Err(PartError::InvalidGeometry(format!(
                "the cam pin slot reaches the ring slot: pin {pin_loaded_r_m:e} + {slot_w_m:e} > {:e} (pin x {pin_x_m:e})",
                iron_r_m - slot_h_m
            )));
        }
        // The pin slots: a radial slot through the disc skin over the stroke,
        // with the arm channel clearance, on the pin side of each pocket.
        let mut slots: Option<Box<dyn Solid>> = None;
        for (k, &az) in azimuths.iter().enumerate() {
            let side = if k % 2 == 0 { 1.0 } else { -1.0 };
            let reach = half_m + 1.0e-9;
            let (z0, z1) = if side > 0.0 {
                (0.0, reach)
            } else {
                (-reach, 0.0)
            };
            let slot: Box<dyn Solid> = Box::new(Placed {
                inner: Box::new(Box3 {
                    min_m: [pin_released_r_m - slot_w_m, pin_y_m - slot_w_m, z0],
                    max_m: [pin_loaded_r_m + slot_w_m, pin_y_m + slot_w_m, z1],
                }),
                rotation_rad: az,
                offset_m: [0.0; 3],
            });
            slots = Some(match slots {
                None => slot,
                Some(prev) => Box::new(Union { a: prev, b: slot }),
            });
        }
        let solid: Box<dyn Solid> = match slots {
            Some(inner) => Box::new(Difference {
                outer: solid,
                inner,
            }),
            None => solid,
        };
        let mut rotor = capped_lattice(fill_solid(solid.as_ref()))?;

        let mut sites = Vec::with_capacity(count);
        let mut arms = Vec::with_capacity(count);
        let mut site_iron = Vec::with_capacity(count);
        let mut site_anchor: Vec<Vec<(usize, usize)>> = Vec::with_capacity(count);
        let mut site_tethers: Vec<Vec<Tether>> = Vec::with_capacity(count);
        let mut site_donor: Vec<Option<DonorSite>> = Vec::with_capacity(count);
        let mut site_catalyst: Vec<Option<CatalystSite>> = Vec::with_capacity(count);
        let mut catalyst_hydrogens: Vec<usize> = Vec::new();
        let mut loaded = Vec::with_capacity(count);
        let mut site_dome_m = Vec::with_capacity(count);
        let mut pin_side = Vec::with_capacity(count);
        for (k, &az) in azimuths.iter().enumerate() {
            let (s, c) = az.sin_cos();
            let u = [c, s, 0.0];
            let t = [-s, c, 0.0];
            let fe = [iron_r_m * c, iron_r_m * s, 0.0];
            let is_loaded = ArmCamTrack::loaded_fraction(az) >= 0.5;
            loaded.push(is_loaded);
            // The self-sustaining groups (ADR-0141) are built on the rotor
            // before the site, so the tether search cannot claim their carbons.
            let donor = if build_donor {
                add_rotor_donor(&mut rotor, fe, u)?
            } else {
                None
            };
            let catalyst = if build_catalyst {
                match find_rotor_catalyst_carbon(&rotor, fe, u, bore_m + 1.5e-10) {
                    Some((c_i, h_i)) => {
                        rotor.set_element(c_i, Element::MANGANESE)?;
                        rotor.set_atom_type(c_i, "Mn")?;
                        catalyst_hydrogens.push(h_i);
                        let p = rotor.position_m(c_i).unwrap_or(fe);
                        Some(CatalystSite {
                            manganese: c_i,
                            fe_distance_m: dist(p, fe),
                        })
                    }
                    None => None,
                }
            } else {
                None
            };
            site_donor.push(donor);
            site_catalyst.push(catalyst);
            // An empty pocket shows the five-coordinate, domed iron (M28-34).
            let dome = if is_loaded { 0.0 } else { iron_dome_m };
            site_dome_m.push(dome);
            let (site, fe_i, tethers) = build_heme_site(fe, u, t, &rotor, is_loaded, dome)?;
            site_anchor.push(tethers.iter().map(|x| (x.site_atom, x.wall_atom)).collect());
            site_tethers.push(tethers);
            site_iron.push(fe_i);
            sites.push(
                Part::new(format!("heme_site_{k}"), site).with_material("Fe-porphine (model)"),
            );

            // The arm: +x of the arm frame points along +u; its N at the
            // ligand distance below the iron.
            // Odd pockets mirror the arm in z, so their pin points down.
            let tip_r_m = iron_r_m
                - if is_loaded {
                    LIGAND_IN_M
                } else {
                    LIGAND_OUT_M + release_margin_m
                };
            let side = if k % 2 == 0 { 1.0 } else { -1.0 };
            pin_side.push(side);
            let mut placed = arm.clone();
            for i in 0..placed.atom_count() {
                if let Some(p) = placed.position_m(i) {
                    let along = tip_r_m + p[0];
                    placed.set_position_m(
                        i,
                        [along * c - p[1] * s, along * s + p[1] * c, side * p[2]],
                    )?;
                }
            }
            arms.push(Part::new(format!("ligand_arm_{k}"), placed).with_material("diamond + NH2"));
        }

        // Each tether and the catalyst replaced a rotor hydrogen: remove those,
        // and renumber.
        let mut gone: Vec<usize> = site_tethers
            .iter()
            .flatten()
            .map(|x| x.wall_hydrogen)
            .collect();
        gone.extend(catalyst_hydrogens.iter().copied());
        let (pruned, map) = rotor.without_atoms(&gone);
        rotor = pruned;
        let remap = |i: usize| -> Result<usize, PartError> {
            map.get(i)
                .copied()
                .flatten()
                .map(|v| v as usize)
                .ok_or_else(|| PartError::InvalidGeometry("a tether anchor was removed".to_owned()))
        };
        for (pairs, tethers) in site_anchor.iter_mut().zip(site_tethers.iter_mut()) {
            for (pair, tether) in pairs.iter_mut().zip(tethers.iter_mut()) {
                tether.wall_atom = remap(tether.wall_atom)?;
                pair.1 = tether.wall_atom;
            }
        }
        for donor in site_donor.iter_mut().flatten() {
            donor.oxygen = remap(donor.oxygen)?;
            donor.hydrogen = remap(donor.hydrogen)?;
        }
        for catalyst in site_catalyst.iter_mut().flatten() {
            catalyst.manganese = remap(catalyst.manganese)?;
        }

        Ok(HemeRotorSet {
            rotor: Part::new("heme_rotor", rotor).with_material("diamond"),
            sites,
            arms,
            site_iron,
            site_anchor,
            site_tethers,
            site_donor,
            site_catalyst,
            arm_tip: vec![cap; count],
            azimuth_rad: azimuths,
            loaded,
            site_dome_m,
            iron_radius_m: iron_r_m,
            disc_radius_m: disc_m,
            half_thickness_m: half_m,
            pin_side,
            pin_radius_m: arm_r_m,
            pin_loaded_radius_m: pin_loaded_r_m,
            pin_released_radius_m: pin_released_r_m,
            release_margin_m,
            pin_top_m: arm_pin_top_m,
            pin_inner_m,
            pin_outer_m,
            cam_inlet_rad: CAM_INLET_RAD,
            cam_outlet_rad: CAM_OUTLET_RAD,
            cam_ramp_half_rad: CAM_RAMP_HALF_RAD,
        })
    }
}

/// Builds one heme-like site: an Fe-porphine with its four meso tethers to
/// `wall`, and the bound O2 when `loaded`. The ring lies in the plane spanned
/// by `t` and `cross(u, t)`; the iron-O2 axis is `u`. Returns the site, the
/// iron index, and the (site atom, wall atom) anchor pairs.
/// One covalent tether from a porphine meso carbon to a wall carbon.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Tether {
    /// The tether end carbon inside the site part, bonded to `wall_atom`.
    pub site_atom: usize,
    /// The wall carbon the tether bonds to.
    pub wall_atom: usize,
    /// The wall hydrogen the tether replaces. The generator removes it, so
    /// after assembly this index no longer names an atom of the wall.
    pub wall_hydrogen: usize,
    /// Carbons between the meso carbon and the wall carbon.
    pub bridge_atoms: usize,
    /// Valence strain of the tether (bonds, angles, meso planarity) in the
    /// harmonic MM3-scale field of `tether_energy`, in kJ/mol (ADR-0137).
    pub strain_kj_mol: f64,
    /// The soft steric penalty of the placed bridge carbons, in kJ/mol. A
    /// placement margin, not a physical energy.
    pub steric_kj_mol: f64,
}

/// One built heme-like site: its topology, its iron index, and its tethers.
type HemeSite = (Topology, usize, Vec<Tether>);

// The tether field (ADR-0137): the MM3-scale constants of the valence
// relaxation (`nanocad_meter::valence_relax`), plus an in-plane term at the sp2
// meso carbon and a soft steric wall. Lengths in angstrom, energy in kJ/mol.
const TETHER_K_BOND: f64 = 449.0 * 1.0e-20 * 6.022_140_76e23 / 1000.0; // kJ/mol/A^2
const TETHER_K_ANGLE: f64 = 0.67e-18 * 6.022_140_76e23 / 1000.0; // kJ/mol/rad^2
const TETHER_R0_CC: f64 = 1.544;
const TETHER_R0_SP2: f64 = 1.51;
const TETHER_STERIC_A: f64 = 3.0;
/// The clearance from a bridge carbon to a hydrogen, in angstrom. The bridge
/// carbons carry their own hydrogens, which are not in the field; this margin
/// keeps the chain off the wall and ring hydrogens.
const TETHER_STERIC_H_A: f64 = 2.0;
/// A tether may use this much more strain than the best one if it is shorter,
/// in kJ/mol: a shorter link holds the ring more rigidly.
const TETHER_SHORTER_MARGIN_KJ_MOL: f64 = 10.0;

fn bend(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let u = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let v = [c[0] - b[0], c[1] - b[1], c[2] - b[2]];
    let nu = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
    let nv = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if nu == 0.0 || nv == 0.0 {
        return 0.0;
    }
    ((u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (nu * nv))
        .clamp(-1.0, 1.0)
        .acos()
}

/// The fixed ends and surroundings of one candidate tether, in angstrom.
struct TetherFrame {
    meso: [f64; 3],
    ring_nb: Vec<[f64; 3]>,
    ring_normal: [f64; 3],
    wall: [f64; 3],
    wall_nb: Vec<[f64; 3]>,
    obstacles: Vec<[f64; 3]>,
    hydrogens: Vec<[f64; 3]>,
}

/// The tether energy for bridge carbons `x` (3 per atom): (valence, steric), in kJ/mol.
fn tether_terms(f: &TetherFrame, x: &[f64]) -> (f64, f64) {
    let n = x.len() / 3;
    let mut chain = Vec::with_capacity(n + 2);
    chain.push(f.meso);
    for i in 0..n {
        chain.push([x[3 * i], x[3 * i + 1], x[3 * i + 2]]);
    }
    chain.push(f.wall);
    let tetra = (-1.0_f64 / 3.0).acos();
    let mut e = 0.0;
    for i in 0..chain.len() - 1 {
        let r0 = if i == 0 { TETHER_R0_SP2 } else { TETHER_R0_CC };
        e += 0.5 * TETHER_K_BOND * (dist(chain[i], chain[i + 1]) - r0).powi(2);
    }
    for r in &f.ring_nb {
        e += 0.5 * TETHER_K_ANGLE * (bend(*r, chain[0], chain[1]) - 120.0_f64.to_radians()).powi(2);
    }
    let d1 = normalize([
        chain[1][0] - chain[0][0],
        chain[1][1] - chain[0][1],
        chain[1][2] - chain[0][2],
    ]);
    let off = (d1[0] * f.ring_normal[0] + d1[1] * f.ring_normal[1] + d1[2] * f.ring_normal[2])
        .clamp(-1.0, 1.0);
    e += 0.5 * TETHER_K_ANGLE * off.asin().powi(2);
    for i in 1..chain.len() - 1 {
        e += 0.5 * TETHER_K_ANGLE * (bend(chain[i - 1], chain[i], chain[i + 1]) - tetra).powi(2);
    }
    let last = chain.len() - 1;
    for w in &f.wall_nb {
        e += 0.5 * TETHER_K_ANGLE * (bend(*w, chain[last], chain[last - 1]) - tetra).powi(2);
    }
    let valence = e;
    let mut e = 0.0;
    for q in &chain[1..last] {
        for o in &f.obstacles {
            let d = dist(*q, *o);
            if d < TETHER_STERIC_A {
                e += 0.5 * TETHER_K_BOND * (TETHER_STERIC_A - d).powi(2);
            }
        }
        for o in &f.hydrogens {
            let d = dist(*q, *o);
            if d < TETHER_STERIC_H_A {
                e += 0.5 * TETHER_K_BOND * (TETHER_STERIC_H_A - d).powi(2);
            }
        }
    }
    (valence, e)
}

fn tether_energy(f: &TetherFrame, x: &[f64]) -> f64 {
    let (v, st) = tether_terms(f, x);
    v + st
}

/// Places `n` bridge carbons by gradient descent from several starts.
/// Returns the positions (angstrom) and the strain (kJ/mol).
fn solve_tether(f: &TetherFrame, n: usize) -> (Vec<[f64; 3]>, f64) {
    if n == 0 {
        return (Vec::new(), tether_energy(f, &[]));
    }
    let axis = normalize([
        f.wall[0] - f.meso[0],
        f.wall[1] - f.meso[1],
        f.wall[2] - f.meso[2],
    ]);
    let side = normalize(cross(axis, f.ring_normal));
    let mut best: (Vec<f64>, f64) = (Vec::new(), f64::INFINITY);
    for start in [
        f.ring_normal,
        side,
        [-f.ring_normal[0], -f.ring_normal[1], -f.ring_normal[2]],
        [-side[0], -side[1], -side[2]],
    ] {
        let mut x = Vec::with_capacity(3 * n);
        for i in 1..=n {
            let t = i as f64 / (n + 1) as f64;
            let bow = 0.8 * (std::f64::consts::PI * t).sin();
            for (k, s) in start.iter().enumerate() {
                x.push(f.meso[k] + t * (f.wall[k] - f.meso[k]) + bow * s);
            }
        }
        let mut e = tether_energy(f, &x);
        let mut step = 0.05;
        let h = 1.0e-5;
        for _ in 0..3000 {
            let mut g = vec![0.0; x.len()];
            for k in 0..x.len() {
                let keep = x[k];
                x[k] = keep + h;
                let ep = tether_energy(f, &x);
                x[k] = keep - h;
                let em = tether_energy(f, &x);
                x[k] = keep;
                g[k] = (ep - em) / (2.0 * h);
            }
            let norm = g.iter().map(|v| v * v).sum::<f64>().sqrt();
            if norm < 1.0e-6 {
                break;
            }
            let trial: Vec<f64> = x
                .iter()
                .zip(&g)
                .map(|(xi, gi)| xi - step * gi / norm)
                .collect();
            let et = tether_energy(f, &trial);
            if et < e {
                x = trial;
                e = et;
                step = (step * 1.2).min(0.2);
            } else {
                step *= 0.5;
                if step < 1.0e-7 {
                    break;
                }
            }
        }
        if e < best.1 {
            best = (x, e);
        }
    }
    let pts = (0..n)
        .map(|i| [best.0[3 * i], best.0[3 * i + 1], best.0[3 * i + 2]])
        .collect();
    (pts, best.1)
}

fn build_heme_site(
    fe: [f64; 3],
    u: [f64; 3],
    t: [f64; 3],
    wall: &Topology,
    loaded: bool,
    iron_dome_m: f64,
) -> Result<HemeSite, PartError> {
    let w = cross(u, t);
    let ring = porphine_plane();
    let mut site = Topology::new();
    // The armed empty site is five-coordinate, so the high-spin iron domes
    // toward the axial nitrogen, out of the ring plane (M28-34). The loaded
    // site is six-coordinate and stays in-plane.
    let fe_center = [
        fe[0] - iron_dome_m * u[0],
        fe[1] - iron_dome_m * u[1],
        fe[2] - iron_dome_m * u[2],
    ];
    let fe_i = site.add_atom(Atom::new(Element::IRON, fe_center, 0.0, "Fe"));
    // The porphine, less its four meso hydrogens: each meso carbon takes a
    // two-carbon tether to the wall instead.
    let mut meso_c: Vec<usize> = Vec::new();
    for (idx, (e, [a, b])) in ring.iter().enumerate() {
        if idx % 9 == 8 {
            continue;
        }
        let p = [
            fe[0] + a * t[0] + b * w[0],
            fe[1] + a * t[1] + b * w[1],
            fe[2] + a * t[2] + b * w[2],
        ];
        let ai = site.add_atom(Atom::new(*e, p, 0.0, e.symbol().unwrap_or("X")));
        if idx % 9 == 7 {
            meso_c.push(ai as usize);
        }
    }
    if loaded {
        let o1 = [
            fe_center[0] + FE_O_M * u[0],
            fe_center[1] + FE_O_M * u[1],
            fe_center[2] + FE_O_M * u[2],
        ];
        let (sn, cs) = 60.0_f64.to_radians().sin_cos();
        let o2 = [
            o1[0] + O_O_M * (cs * u[0] + sn * t[0]),
            o1[1] + O_O_M * (cs * u[1] + sn * t[1]),
            o1[2] + O_O_M * (cs * u[2] + sn * t[2]),
        ];
        site.add_atom(Atom::new(Element::OXYGEN, o1, 0.0, "O_O2"));
        site.add_atom(Atom::new(Element::OXYGEN, o2, 0.0, "O_O2"));
    }
    bond_by_distance(&mut site)?;
    // A tether at each meso carbon (ADR-0137). Each candidate anchor is a wall
    // carbon that carries a hydrogen; the tether replaces that hydrogen, so
    // the anchor keeps four bonds. For 1 to 3 bridge carbons the bridge is
    // placed by minimizing a harmonic valence strain with both ends fixed; the
    // least strained tether wins, and a shorter one wins within a margin.
    let a_per_m = 1.0e10;
    let to_a = |p: [f64; 3]| [p[0] * a_per_m, p[1] * a_per_m, p[2] * a_per_m];
    let mut wall_nbrs: Vec<Vec<usize>> = vec![Vec::new(); wall.atom_count()];
    for b in wall.bonds() {
        wall_nbrs[b.u as usize].push(b.v as usize);
        wall_nbrs[b.v as usize].push(b.u as usize);
    }
    let mut site_nbrs: Vec<Vec<usize>> = vec![Vec::new(); site.atom_count()];
    for b in site.bonds() {
        site_nbrs[b.u as usize].push(b.v as usize);
        site_nbrs[b.v as usize].push(b.u as usize);
    }
    let fe_a = to_a(fe);
    let mut used_hydrogen: Vec<usize> = Vec::new();
    let mut anchors: Vec<Tether> = Vec::new();
    for &mi in &meso_c {
        let m = site.position_m(mi).unwrap_or([0.0; 3]);
        let ma = to_a(m);
        let meso_dir = normalize([ma[0] - fe_a[0], ma[1] - fe_a[1], ma[2] - fe_a[2]]);
        let ring_nb: Vec<[f64; 3]> = site_nbrs[mi]
            .iter()
            .filter_map(|&k| site.position_m(k).map(to_a))
            .collect();
        // Heavy atoms the bridge must clear: the site so far, and the wall near the meso carbon.
        let mut site_heavy: Vec<[f64; 3]> = Vec::new();
        let mut site_h: Vec<[f64; 3]> = Vec::new();
        for k in 0..site.atom_count() {
            if k == mi || site_nbrs[mi].contains(&k) {
                continue;
            }
            if let Some(p) = site.position_m(k) {
                if site.element(k) == Some(Element::HYDROGEN) {
                    site_h.push(to_a(p));
                } else {
                    site_heavy.push(to_a(p));
                }
            }
        }
        let mut best: Option<(Tether, Vec<[f64; 3]>, [f64; 3])> = None;
        let mut options: Vec<(Tether, Vec<[f64; 3]>, [f64; 3])> = Vec::new();
        for (j, nbrs) in wall_nbrs.iter().enumerate() {
            if wall.element(j) != Some(Element::CARBON) {
                continue;
            }
            let wa = match wall.position_m(j) {
                Some(p) => to_a(p),
                None => continue,
            };
            let span = dist(ma, wa);
            if !(2.0..=6.0).contains(&span) {
                continue;
            }
            // The tether is an in-plane link: reject an anchor far off the ring
            // plane, so the proximal pocket does not pull it out of plane.
            if (wa[2] - ma[2]).abs() > 0.5 * span {
                continue;
            }
            let ahead = (wa[0] - ma[0]) * meso_dir[0]
                + (wa[1] - ma[1]) * meso_dir[1]
                + (wa[2] - ma[2]) * meso_dir[2];
            if ahead < 0.2 * span {
                continue;
            }
            let hs: Vec<usize> = nbrs
                .iter()
                .copied()
                .filter(|&k| {
                    wall.element(k) == Some(Element::HYDROGEN) && !used_hydrogen.contains(&k)
                })
                .collect();
            // The hydrogen that points most toward the meso carbon is replaced.
            let mut pick: Option<(usize, f64)> = None;
            for &hk in &hs {
                let hp = to_a(wall.position_m(hk).unwrap_or([0.0; 3]));
                let toward = normalize([ma[0] - wa[0], ma[1] - wa[1], ma[2] - wa[2]]);
                let hd = normalize([hp[0] - wa[0], hp[1] - wa[1], hp[2] - wa[2]]);
                let c = toward[0] * hd[0] + toward[1] * hd[1] + toward[2] * hd[2];
                if c > 0.0 && pick.is_none_or(|(_, best_c)| c > best_c) {
                    pick = Some((hk, c));
                }
            }
            let Some((hk, _)) = pick else { continue };
            let wall_nb: Vec<[f64; 3]> = nbrs
                .iter()
                .filter(|&&k| k != hk)
                .filter_map(|&k| wall.position_m(k).map(to_a))
                .collect();
            let mut obstacles = site_heavy.clone();
            let mut hydrogens = site_h.clone();
            for k in 0..wall.atom_count() {
                if k == j || k == hk {
                    continue;
                }
                if let Some(p) = wall.position_m(k) {
                    let pa = to_a(p);
                    if dist(pa, ma) >= 8.0 {
                        continue;
                    }
                    // Far off the ring plane, e.g. below the proximal pocket.
                    if (pa[2] - ma[2]).abs() > 1.5 {
                        continue;
                    }
                    if wall.element(k) == Some(Element::HYDROGEN) {
                        hydrogens.push(pa);
                    } else if !wall_nbrs[j].contains(&k) {
                        obstacles.push(pa);
                    }
                }
            }
            let frame = TetherFrame {
                meso: ma,
                ring_nb: ring_nb.clone(),
                ring_normal: u,
                wall: wa,
                wall_nb,
                obstacles,
                hydrogens,
            };
            for n in 1..=3usize {
                if span > 1.25 * TETHER_R0_CC * (n + 1) as f64 {
                    continue;
                }
                let (pts, _) = solve_tether(&frame, n);
                let flat: Vec<f64> = pts.iter().flat_map(|q| q.iter().copied()).collect();
                let (strain, steric) = tether_terms(&frame, &flat);
                options.push((
                    Tether {
                        site_atom: 0,
                        wall_atom: j,
                        wall_hydrogen: hk,
                        bridge_atoms: n,
                        strain_kj_mol: strain,
                        steric_kj_mol: steric,
                    },
                    pts,
                    wa,
                ));
            }
        }
        let total = |t: &Tether| t.strain_kj_mol + t.steric_kj_mol;
        let lowest = options
            .iter()
            .map(|o| total(&o.0))
            .fold(f64::INFINITY, f64::min);
        for o in options {
            if total(&o.0) > lowest + TETHER_SHORTER_MARGIN_KJ_MOL {
                continue;
            }
            let better = match &best {
                None => true,
                Some(b) => {
                    o.0.bridge_atoms < b.0.bridge_atoms
                        || (o.0.bridge_atoms == b.0.bridge_atoms && total(&o.0) < total(&b.0))
                }
            };
            if better {
                best = Some(o);
            }
        }
        let Some((mut tether, pts, wa)) = best else {
            return Err(PartError::InvalidGeometry(
                "a meso carbon has no wall carbon with a hydrogen in tether reach".to_owned(),
            ));
        };
        let wq = [wa[0] / a_per_m, wa[1] / a_per_m, wa[2] / a_per_m];
        used_hydrogen.push(tether.wall_hydrogen);
        let mut chain: Vec<usize> = Vec::new();
        for p in &pts {
            let pm = [p[0] / a_per_m, p[1] / a_per_m, p[2] / a_per_m];
            chain.push(site.add_atom(Atom::new(Element::CARBON, pm, 0.0, "C")) as usize);
        }
        let mut prev = mi;
        for &ci in &chain {
            site.add_bond(Bond::new(prev as u32, ci as u32, 1, BondType::Single))?;
            prev = ci;
        }
        for (idx, &ci) in chain.iter().enumerate() {
            let p = site.position_m(ci).unwrap_or([0.0; 3]);
            let n1 = if idx == 0 {
                m
            } else {
                site.position_m(chain[idx - 1]).unwrap_or([0.0; 3])
            };
            let n2 = if idx + 1 < chain.len() {
                site.position_m(chain[idx + 1]).unwrap_or([0.0; 3])
            } else {
                wq
            };
            let a1 = normalize([n1[0] - p[0], n1[1] - p[1], n1[2] - p[2]]);
            let a2 = normalize([n2[0] - p[0], n2[1] - p[1], n2[2] - p[2]]);
            let s = normalize([-(a1[0] + a2[0]), -(a1[1] + a2[1]), -(a1[2] + a2[2])]);
            let q = normalize(cross(a1, a2));
            for sign in [1.0, -1.0] {
                let hdir = normalize([s[0] + sign * q[0], s[1] + sign * q[1], s[2] + sign * q[2]]);
                let h = site.add_atom(Atom::new(
                    Element::HYDROGEN,
                    [
                        p[0] + 1.09e-10 * hdir[0],
                        p[1] + 1.09e-10 * hdir[1],
                        p[2] + 1.09e-10 * hdir[2],
                    ],
                    0.0,
                    "H",
                ));
                site.add_bond(Bond::new(ci as u32, h, 1, BondType::Single))?;
            }
        }
        tether.site_atom = *chain.last().unwrap_or(&mi);
        anchors.push(tether);
    }
    Ok((site, fe_i as usize, anchors))
}

/// The longest strained surface dimer the module parts accept, in metres.
/// A longer pair at a step is resolved by removing one carbon (a CH2 step)
/// instead. With 2.6 A (the default elsewhere) the relaxed rim bent by up to
/// 5.6 A; with 1.9 A the removal eroded the rotor and it collapsed; 2.1 A
/// settles in a few rounds and relaxes with RMS 1.0 A (ADR-0111).
pub const MODULE_MAX_DIMER_M: f64 = 2.1e-10;

/// Removes lattice sites until no missing site would take two caps after the
/// dimer pass, then returns the capped topology. Each round drops one host of
/// every shared site: the host with fewer carbon bonds.
pub(crate) fn capped_lattice(sites: Vec<[f64; 3]>) -> Result<Topology, PartError> {
    use std::collections::{BTreeSet, HashMap};
    let mut sites = keep_well_bonded(sites);
    for _ in 0..40 {
        let mut t = Topology::new();
        let ids: Vec<u32> = sites
            .iter()
            .map(|p| t.add_atom(Atom::new(Element::CARBON, *p, 0.0, "C")))
            .collect();
        diamond_solid::bond_and_cap_limited(&mut t, &ids, MODULE_MAX_DIMER_M)?;
        // The carbons are atoms 0..sites.len(); anything after is a cap or
        // a filled notch. Find hydrogen pairs closer than 0.9 A.
        let n_sites = sites.len();
        let mut host = HashMap::new();
        let mut degree = vec![0usize; t.atom_count()];
        for b in t.bonds() {
            let (u, v) = (b.u as usize, b.v as usize);
            degree[u] += 1;
            degree[v] += 1;
            if t.element(u) == Some(Element::HYDROGEN) {
                host.insert(u, v);
            }
            if t.element(v) == Some(Element::HYDROGEN) {
                host.insert(v, u);
            }
        }
        let hs: Vec<usize> = (0..t.atom_count())
            .filter(|&i| t.element(i) == Some(Element::HYDROGEN))
            .collect();
        let cell = |p: [f64; 3]| {
            (
                (p[0] / 3e-10).floor() as i64,
                (p[1] / 3e-10).floor() as i64,
                (p[2] / 3e-10).floor() as i64,
            )
        };
        let mut grid: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
        for &h in &hs {
            if let Some(p) = t.position_m(h) {
                grid.entry(cell(p)).or_default().push(h);
            }
        }
        let mut drop: BTreeSet<usize> = BTreeSet::new();
        for &h in &hs {
            let Some(p) = t.position_m(h) else { continue };
            let (x, y, z) = cell(p);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        for &g in grid.get(&(x + dx, y + dy, z + dz)).into_iter().flatten() {
                            if g <= h {
                                continue;
                            }
                            let Some(q) = t.position_m(g) else { continue };
                            if dist(p, q) >= 0.9e-10 {
                                continue;
                            }
                            let (a, b) = (host[&h], host[&g]);
                            if a >= n_sites
                                || b >= n_sites
                                || drop.contains(&a)
                                || drop.contains(&b)
                            {
                                continue;
                            }
                            // Drop the host with fewer bonds (ties: higher index).
                            let pick = if (degree[a], std::cmp::Reverse(a))
                                <= (degree[b], std::cmp::Reverse(b))
                            {
                                a
                            } else {
                                b
                            };
                            drop.insert(pick);
                        }
                    }
                }
            }
        }
        if drop.is_empty() {
            return Ok(t);
        }
        sites = keep_well_bonded(
            sites
                .into_iter()
                .enumerate()
                .filter(|(i, _)| !drop.contains(i))
                .map(|(_, p)| p)
                .collect(),
        );
    }
    Err(PartError::InvalidGeometry(
        "the surface still has clashing caps after 40 rounds".to_owned(),
    ))
}

/// A rod of `radius_m` from the origin along the unit `axis` for `length_m`.
#[derive(Clone, Debug)]
struct AxisRod {
    origin: [f64; 3],
    axis: [f64; 3],
    radius_m: f64,
    length_m: f64,
}

impl Solid for AxisRod {
    fn contains_m(&self, p: [f64; 3]) -> bool {
        let p = [
            p[0] - self.origin[0],
            p[1] - self.origin[1],
            p[2] - self.origin[2],
        ];
        let t = p[0] * self.axis[0] + p[1] * self.axis[1] + p[2] * self.axis[2];
        if !(0.0..=self.length_m).contains(&t) {
            return false;
        }
        let q = [
            p[0] - t * self.axis[0],
            p[1] - t * self.axis[1],
            p[2] - t * self.axis[2],
        ];
        q[0] * q[0] + q[1] * q[1] + q[2] * q[2] <= self.radius_m * self.radius_m
    }

    fn bounds_m(&self) -> Bounds {
        let r = self.radius_m;
        let e = [
            self.length_m * self.axis[0],
            self.length_m * self.axis[1],
            self.length_m * self.axis[2],
        ];
        let o = self.origin;
        Bounds {
            min_m: [
                o[0] + e[0].min(0.0) - r,
                o[1] + e[1].min(0.0) - r,
                o[2] + e[2].min(0.0) - r,
            ],
            max_m: [
                o[0] + e[0].max(0.0) + r,
                o[1] + e[1].max(0.0) + r,
                o[2] + e[2].max(0.0) + r,
            ],
        }
    }

    fn clone_box(&self) -> Box<dyn Solid> {
        Box::new(self.clone())
    }
}

fn apply3(m: &[[f64; 3]; 3], p: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * p[0] + m[0][1] * p[1] + m[0][2] * p[2],
        m[1][0] * p[0] + m[1][1] * p[1] + m[1][2] * p[2],
        m[2][0] * p[0] + m[2][1] * p[1] + m[2][2] * p[2],
    ]
}

/// The least distance from a new hydroxyl hydrogen to any other atom, in
/// metres.
pub const MIN_HYDROXYL_CLEARANCE_M: f64 = 2.0e-10;

/// The (carbon, hydrogen) cap bonds of `topology` whose hydrogen satisfies
/// `select`.
#[allow(dead_code)]
fn caps_where(topology: &Topology, select: impl Fn([f64; 3]) -> bool) -> Vec<(usize, usize)> {
    topology
        .bonds()
        .filter_map(|b| {
            let (u, v) = (b.u as usize, b.v as usize);
            match (topology.element(u), topology.element(v)) {
                (Some(Element::CARBON), Some(Element::HYDROGEN)) => Some((u, v)),
                (Some(Element::HYDROGEN), Some(Element::CARBON)) => Some((v, u)),
                _ => None,
            }
        })
        .filter(|&(_, h)| topology.position_m(h).is_some_and(&select))
        .collect()
}

/// Turns each cap hydrogen of `caps` into a hydroxyl group. The oxygen sits
/// on the old C-H direction at 1.43 A. The hydrogen sits on the C-O-H cone
/// (108.5 degrees), at the one of 24 turns about C-O farthest from every
/// neighbour. A real hydroxyl turns about C-O in picoseconds, so the start
/// only has to be free of clashes (ADR-0119). A cap whose best turn comes
/// within [`MIN_HYDROXYL_CLEARANCE_M`] of a neighbour is left as C-H.
/// Returns the number of hydroxyls made.
fn hydroxylate(topology: &mut Topology, caps: &[(usize, usize)]) -> Result<usize, PartError> {
    let mut done = 0usize;
    use std::collections::HashMap;
    let cell = |p: [f64; 3]| {
        (
            (p[0] / 4.0e-10).floor() as i64,
            (p[1] / 4.0e-10).floor() as i64,
            (p[2] / 4.0e-10).floor() as i64,
        )
    };
    let mut grid: HashMap<(i64, i64, i64), Vec<[f64; 3]>> = HashMap::new();
    for i in 0..topology.atom_count() {
        if let Some(p) = topology.position_m(i) {
            grid.entry(cell(p)).or_default().push(p);
        }
    }
    let moved =
        |grid: &mut HashMap<(i64, i64, i64), Vec<[f64; 3]>>, from: [f64; 3], to: [f64; 3]| {
            if let Some(v) = grid.get_mut(&cell(from)) {
                v.retain(|q| dist(*q, from) > 1.0e-14);
            }
            grid.entry(cell(to)).or_default().push(to);
        };
    for &(c, h) in caps {
        let (Some(cp), Some(hp)) = (topology.position_m(c), topology.position_m(h)) else {
            continue;
        };
        let d = normalize([hp[0] - cp[0], hp[1] - cp[1], hp[2] - cp[2]]);
        let o = [
            cp[0] + 1.43e-10 * d[0],
            cp[1] + 1.43e-10 * d[1],
            cp[2] + 1.43e-10 * d[2],
        ];
        let helper = if d[0].abs() < 0.9 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        let e1 = normalize(cross(d, helper));
        let e2 = cross(d, e1);
        let (sa, ca) = (71.5_f64.to_radians().sin(), 71.5_f64.to_radians().cos());
        let (ox, oy, oz) = cell(o);
        let near: Vec<[f64; 3]> = (-1..=1)
            .flat_map(|dx| {
                (-1..=1).flat_map(move |dy| (-1..=1).map(move |dz| (ox + dx, oy + dy, oz + dz)))
            })
            .filter_map(|k| grid.get(&k))
            .flatten()
            .copied()
            .filter(|q| dist(*q, o) > 1.0e-12 && dist(*q, hp) > 1.0e-14)
            .collect();
        let mut best = (f64::NEG_INFINITY, o);
        for step in 0..24 {
            let (sp, cp2) = (2.0 * PI * step as f64 / 24.0).sin_cos();
            let w = [
                cp2 * e1[0] + sp * e2[0],
                cp2 * e1[1] + sp * e2[1],
                cp2 * e1[2] + sp * e2[2],
            ];
            let hh = [
                o[0] + 0.96e-10 * (ca * d[0] + sa * w[0]),
                o[1] + 0.96e-10 * (ca * d[1] + sa * w[1]),
                o[2] + 0.96e-10 * (ca * d[2] + sa * w[2]),
            ];
            let clear = near
                .iter()
                .filter(|q| dist(**q, cp) > 1.0e-12)
                .map(|q| dist(*q, hh))
                .fold(f64::INFINITY, f64::min);
            if clear > best.0 {
                best = (clear, hh);
            }
        }
        // A crowded cap stays C-H: its hydroxyl could not turn clear, and the
        // relaxation pressed such O-H bonds below 0.9 A.
        if best.0 < MIN_HYDROXYL_CLEARANCE_M {
            continue;
        }
        topology.set_element(h, Element::OXYGEN)?;
        topology.set_atom_type(h, crate::group::FunctionalGroup::Hydroxyl.heavy_type())?;
        topology.set_position_m(h, o)?;
        moved(&mut grid, hp, o);
        let hi = topology.add_atom(Atom::new(Element::HYDROGEN, best.1, 0.0, "H"));
        topology.add_bond(Bond::new(h as u32, hi, 1, BondType::Single))?;
        grid.entry(cell(best.1)).or_default().push(best.1);
        done += 1;
    }
    Ok(done)
}

/// The geometry of the water-return plunger rotor (ADR-0122).
#[derive(Clone, Debug)]
pub struct WaterRotorSpec {
    /// The rim radius of the disc, in metres.
    pub disc_radius_m: f64,
    /// The disc thickness, in metres.
    pub thickness_m: f64,
    /// The number of plunger pockets.
    pub pocket_count: usize,
    /// The azimuth of pocket 0, in radians.
    pub phase_rad: f64,
    /// The plunger lattice radius, in metres.
    pub plunger_radius_m: f64,
    /// The carbon gap between plunger and bore (a contact-gap seal), in metres.
    pub plunger_gap_m: f64,
    /// The plunger length, in metres.
    pub plunger_length_m: f64,
    /// The plunger stroke, in metres. The pocket is this deep when withdrawn.
    pub stroke_m: f64,
    /// The plunger tip radius when flush, below the rim, in metres.
    pub flush_inset_m: f64,
    /// The cam pin offset from the plunger back end, in metres.
    pub pin_offset_m: f64,
    /// The cam pin end height above the mid-plane, in metres.
    pub pin_top_m: f64,
    /// The carbon gap between the pin and its skin slot, in metres.
    pub pin_slot_gap_m: f64,
    /// Where the cam pushes each plunger flush (blood side), in radians.
    pub inlet_rad: f64,
    /// Where the cam withdraws each plunger (tank side), in radians.
    pub outlet_rad: f64,
    /// The half-width of each cam ramp, in radians. A 0.6 nm stroke over
    /// +-20 degrees at r = 2.4 nm keeps the steepest slope near 23 degrees.
    pub ramp_half_rad: f64,
}

impl Default for WaterRotorSpec {
    fn default() -> Self {
        Self {
            disc_radius_m: 4.0e-9,
            thickness_m: 2.4e-9,
            pocket_count: 6,
            phase_rad: 15.0_f64.to_radians(),
            plunger_radius_m: 2.6e-10,
            plunger_gap_m: 4.0e-10,
            plunger_length_m: 1.9e-9,
            stroke_m: 0.6e-9,
            flush_inset_m: 1.1e-10,
            pin_offset_m: 1.0e-10,
            pin_top_m: 2.5e-9,
            pin_slot_gap_m: 4.6e-10,
            // Midway between pockets: a pin spans about +-11 degrees, so a
            // pocket 3 degrees past a ramp end touched it (ADR-0122).
            inlet_rad: 345.0_f64.to_radians(),
            outlet_rad: 165.0_f64.to_radians(),
            ramp_half_rad: 20.0_f64.to_radians(),
        }
    }
}

/// Builds the water-return rotor (ADR-0122): a positive-displacement
/// plunger rotor. Each radial bore holds a plunger on the arm-cam
/// mechanism of the O2 rotor (ADR-0116). In the tank sump the cam withdraws
/// the plunger, and the open bore fills with sump water; across the seal and
/// in the blood the cam pushes it flush, and the water leaves into the blood;
/// on the way back to the tank the pocket has no volume. The plunger runs in
/// a contact-gap bore, the seal of the arm collar (ADR-0121). Returns a set
/// with no sites; the plungers are the arms.
pub fn water_rotor(spec: &WaterRotorSpec) -> Result<HemeRotorSet, PartError> {
    let disc_m = spec.disc_radius_m;
    let half_m = 0.5 * spec.thickness_m;
    let count = spec.pocket_count;
    let bore_m = spec.plunger_radius_m + spec.plunger_gap_m;
    if bore_m >= half_m - 2.0e-10 {
        return Err(PartError::InvalidGeometry(
            "the plunger bore is wider than the disc skins allow".to_owned(),
        ));
    }
    let azimuths: Vec<f64> = (0..count)
        .map(|k| (spec.phase_rad + 2.0 * PI * k as f64 / count as f64).rem_euclid(2.0 * PI))
        .collect();
    let tip_in_m = disc_m - spec.flush_inset_m;
    let tip_out_m = tip_in_m - spec.stroke_m;

    // The plunger in its own frame: [111] to +x, the pin along [1-10] to +z.
    let a111 = [1.0 / 3f64.sqrt(); 3];
    let e3 = [1.0 / 2f64.sqrt(), -1.0 / 2f64.sqrt(), 0.0];
    let turn = [a111, cross(e3, a111), e3];
    let plunger_solid = Union {
        a: Box::new(AxisRod {
            origin: [0.0; 3],
            axis: a111,
            radius_m: spec.plunger_radius_m,
            length_m: spec.plunger_length_m,
        }),
        b: Box::new(AxisRod {
            origin: [
                spec.pin_offset_m * a111[0],
                spec.pin_offset_m * a111[1],
                spec.pin_offset_m * a111[2],
            ],
            axis: e3,
            radius_m: spec.plunger_radius_m,
            length_m: spec.pin_top_m,
        }),
    };
    let mut plunger = Topology::new();
    let ids: Vec<u32> = keep_well_bonded(fill_solid(&plunger_solid))
        .into_iter()
        .map(|p| plunger.add_atom(Atom::new(Element::CARBON, p, 0.0, "C")))
        .collect();
    diamond_solid::bond_and_cap(&mut plunger, &ids)?;
    for i in 0..plunger.atom_count() {
        if let Some(p) = plunger.position_m(i) {
            plunger.set_position_m(i, apply3(&turn, p))?;
        }
    }
    // Shift so the outermost tip carbon sits at x = 0, on the bore axis.
    let carbons: Vec<[f64; 3]> = (0..plunger.atom_count())
        .filter(|&i| plunger.element(i) == Some(Element::CARBON))
        .filter_map(|i| plunger.position_m(i))
        .collect();
    let tip_x = carbons
        .iter()
        .filter(|p| p[2].abs() < spec.plunger_radius_m + 1.0e-10)
        .map(|p| p[0])
        .fold(f64::NEG_INFINITY, f64::max);
    let axis_y = 0.0;
    for i in 0..plunger.atom_count() {
        if let Some(p) = plunger.position_m(i) {
            plunger.set_position_m(i, [p[0] - tip_x, p[1] - axis_y, p[2]])?;
        }
    }
    let pin_x = apply3(
        &turn,
        [
            spec.pin_offset_m * a111[0],
            spec.pin_offset_m * a111[1],
            spec.pin_offset_m * a111[2],
        ],
    )[0] - tip_x;
    let pin_y = 0.0;
    let back_x = -spec.plunger_length_m - 0.5e-10;
    let (mut pin_inner_m, mut pin_outer_m, mut pin_top_m) = (0.0_f64, 0.0_f64, half_m);
    for i in 0..plunger.atom_count() {
        if plunger.element(i) != Some(Element::CARBON) {
            continue;
        }
        if let Some(p) = plunger.position_m(i) {
            pin_top_m = pin_top_m.max(p[2]);
            if p[2] > half_m {
                pin_inner_m = pin_inner_m.max(pin_x - p[0]);
                pin_outer_m = pin_outer_m.max(p[0] - pin_x);
            }
        }
    }
    let pin_in_r = (tip_in_m + pin_x).hypot(pin_y);
    let pin_out_r = (tip_out_m + pin_x).hypot(pin_y);
    // The skin slot of the pin must end behind the withdrawn tip, or the
    // pocket opens into the face cavity and on to the other side (ADR-0122:
    // the first build leaked this way).
    let slot_reach_m = pin_in_r + spec.plunger_radius_m + spec.pin_slot_gap_m;
    if slot_reach_m > tip_out_m - 2.0e-10 {
        return Err(PartError::InvalidGeometry(format!(
            "the pin slot reaches the pocket: slot to {slot_reach_m:e} m, withdrawn tip at {tip_out_m:e} m"
        )));
    }

    // The disc: one bore per pocket from behind the withdrawn plunger to the
    // rim, and a pin slot through the skin on the pin side.
    let mut cuts: Option<Box<dyn Solid>> = None;
    let mut pin_side = Vec::with_capacity(count);
    for (k, &az) in azimuths.iter().enumerate() {
        let side = if k % 2 == 0 { 1.0 } else { -1.0 };
        pin_side.push(side);
        let from = tip_out_m + back_x - 3.0e-10;
        let to = disc_m + 1.0e-9;
        let slot_w = spec.plunger_radius_m + spec.pin_slot_gap_m;
        let reach = half_m + 1.0e-9;
        let (z0, z1) = if side > 0.0 {
            (0.0, reach)
        } else {
            (-reach, 0.0)
        };
        let pieces: [Box<dyn Solid>; 2] = [
            Box::new(RadialCylinder {
                azimuth_rad: az,
                center_radius_m: 0.5 * (from + to),
                radius_m: bore_m,
                half_length_m: 0.5 * (to - from),
                center_z_m: 0.0,
            }),
            Box::new(Placed {
                inner: Box::new(Box3 {
                    min_m: [pin_out_r - slot_w, pin_y - slot_w, z0],
                    max_m: [pin_in_r + slot_w, pin_y + slot_w, z1],
                }),
                rotation_rad: az,
                offset_m: [0.0; 3],
            }),
        ];
        for piece in pieces {
            cuts = Some(match cuts {
                None => piece,
                Some(prev) => Box::new(Union { a: prev, b: piece }),
            });
        }
    }
    let disc = Cylinder {
        center_m: [0.0, 0.0, 0.0],
        radius_m: disc_m,
        height_m: spec.thickness_m,
    };
    let solid: Box<dyn Solid> = match cuts {
        Some(inner) => Box::new(Difference {
            outer: Box::new(disc),
            inner,
        }),
        None => Box::new(disc),
    };
    let rotor = capped_lattice(fill_solid(solid.as_ref()))?;

    // Place each plunger by the cam: flush where the cam is in, withdrawn
    // where it is out. Odd pockets mirror in z so their pin points down.
    let mut arms = Vec::with_capacity(count);
    let mut flush = Vec::with_capacity(count);
    for (k, &az) in azimuths.iter().enumerate() {
        let f = ArmCamTrack::fraction_with(az, spec.inlet_rad, spec.outlet_rad, spec.ramp_half_rad);
        if f != 0.0 && f != 1.0 {
            return Err(PartError::InvalidGeometry(format!(
                "pocket {k} sits on a cam ramp"
            )));
        }
        flush.push(f == 1.0);
        let tip = if f == 1.0 { tip_in_m } else { tip_out_m };
        let (sn, cs) = az.sin_cos();
        let mut placed = plunger.clone();
        for i in 0..placed.atom_count() {
            if let Some(p) = placed.position_m(i) {
                let along = tip + p[0];
                placed.set_position_m(
                    i,
                    [
                        along * cs - p[1] * sn,
                        along * sn + p[1] * cs,
                        pin_side[k] * p[2],
                    ],
                )?;
            }
        }
        arms.push(Part::new(format!("plunger_{k}"), placed).with_material("diamond"));
    }
    Ok(HemeRotorSet {
        rotor: Part::new("water_rotor", rotor).with_material("diamond"),
        sites: Vec::new(),
        arms,
        site_iron: Vec::new(),
        arm_tip: Vec::new(),
        site_anchor: Vec::new(),
        site_tethers: Vec::new(),
        site_donor: Vec::new(),
        site_catalyst: Vec::new(),
        azimuth_rad: azimuths,
        loaded: flush,
        site_dome_m: Vec::new(),
        iron_radius_m: tip_in_m,
        disc_radius_m: disc_m,
        half_thickness_m: half_m,
        pin_side,
        pin_radius_m: spec.plunger_radius_m,
        pin_loaded_radius_m: pin_in_r,
        pin_released_radius_m: pin_out_r,
        pin_top_m,
        // A displacement pump tolerates play; 0.4 A of slack on each groove
        // wall clears the lattice steps (a withdrawn pin came to 1.8 A).
        pin_inner_m: pin_inner_m + 0.4e-10,
        pin_outer_m: pin_outer_m + 0.4e-10,
        release_margin_m: 0.0,
        cam_inlet_rad: spec.inlet_rad,
        cam_outlet_rad: spec.outlet_rad,
        cam_ramp_half_rad: spec.ramp_half_rad,
    })
}

/// The hull-wall section that holds the O2 rotor (ADR-0117).
///
/// The wall is a diamond slab on the plane x = 0: blood on the +x side, the
/// tank on the -x side. The rotor axis (z) lies in the wall plane, so half the
/// rim faces the blood and half the tank. Two end plates cover the rotor faces
/// and carry the arm-cam grooves. The slab closes the chamber only where the
/// rim crosses the wall: those two arcs are the seals. Everywhere else the rim
/// is open to the blood or to the tank. One drain bore through the slab,
/// beside the plates, returns water from the tank to the blood; the
/// hydrogens facing into it are hydroxyl groups (ADR-0114).
#[derive(Clone, Debug)]
pub struct HullWallSpec {
    /// The carbon gap between the rotor rim and the chamber wall, in metres.
    pub rim_gap_m: f64,
    /// The carbon gap between the rotor faces and the end plates, in metres.
    pub face_gap_m: f64,
    /// The carbon gap of the face seal: a ring on each plate at the outer
    /// edge of the rotor face, at H-H contact, so the face gap (with the cam
    /// grooves and pin slots) is closed to both blood and tank, in metres.
    pub face_seal_gap_m: f64,
    /// The radial width of the face seal ring, in metres.
    pub face_seal_width_m: f64,
    /// The lattice behind each cam groove floor, in metres.
    pub plate_backing_m: f64,
    /// How far the end plates reach past the chamber radius, in metres.
    pub plate_overhang_m: f64,
    /// The slab thickness along x, in metres. It sets the seal arcs.
    pub wall_thickness_m: f64,
    /// How far the slab reaches past the end plates along y, in metres.
    pub wall_margin_m: f64,
    /// Whether the wall has a drain bore at all. Every bore tested either
    /// carried no water (OH-lined) or let O2 through (H-capped); the rotor
    /// returns the water itself when the tank is wet (ADR-0120, ADR-0122).
    pub drain: bool,
    /// The drain bore cut radius (to the wall carbon centres), in metres.
    pub drain_radius_m: f64,
    /// The length of the narrow drain section, in metres. A long single-file
    /// bore half-empties at idle and lets an O2 sit in it; the 1.0 nm channel
    /// of ADR-0098 stays full (ADR-0119).
    pub drain_narrow_length_m: f64,
    /// The radius of the counterbore from each wall face to the narrow
    /// section, in metres.
    pub drain_counterbore_radius_m: f64,
    /// Which bore caps become hydroxyl groups: 1 every cap, n every n-th
    /// (ordered along the bore, then by angle), 0 none (an H-capped bore,
    /// the default: every OH-lined bore tested was plugged by its own
    /// turning hydrogens, while H-capped bores carry water, ADR-0120).
    pub drain_hydroxyl_every: usize,
}

impl Default for HullWallSpec {
    fn default() -> Self {
        Self {
            rim_gap_m: 4.2e-10,
            face_gap_m: 7.0e-10,
            face_seal_gap_m: 4.6e-10,
            face_seal_width_m: 0.5e-9,
            plate_backing_m: 0.5e-9,
            plate_overhang_m: 0.8e-9,
            wall_thickness_m: 2.4e-9,
            wall_margin_m: 1.6e-9,
            drain: false,
            drain_radius_m: 3.6e-10,
            drain_narrow_length_m: 1.0e-9,
            drain_counterbore_radius_m: 7.0e-10,
            drain_hydroxyl_every: 0,
        }
    }
}

/// A built hull-wall section.
#[derive(Clone, Debug)]
pub struct HullWall {
    /// The wall, end plates and seals, one body.
    pub part: Part,
    /// The chamber radius, in metres.
    pub chamber_radius_m: f64,
    /// The chamber half-height (the plate inner faces), in metres.
    pub chamber_half_m: f64,
    /// The half-angle of each seal arc about 90 and 270 degrees, in radians.
    pub seal_half_rad: f64,
    /// The y of the drain bore axis (it runs along x at z = 0), in metres.
    pub drain_y_m: f64,
    /// The number of hydroxyl groups lining the drain bore.
    pub drain_lining: usize,
    /// The inner radius of the face seal ring, in metres.
    pub face_seal_inner_radius_m: f64,
}

/// Builds the hull-wall section around `set` (ADR-0117).
pub fn hull_wall(set: &HemeRotorSet, spec: &HullWallSpec) -> Result<HullWall, PartError> {
    let chamber_r = set.disc_radius_m + spec.rim_gap_m;
    let chamber_half = set.half_thickness_m + spec.face_gap_m;
    let [mut top, mut bottom] = set.cam_tracks(chamber_half);
    if set.arms.is_empty() {
        // No arms, no cam: empty grooves and a plain plate.
        for track in [&mut top, &mut bottom] {
            track.z_min_m = 1.0;
            track.z_max_m = 0.0;
        }
    }
    let plate_top = if set.arms.is_empty() {
        chamber_half + 0.6e-9 + spec.plate_backing_m
    } else {
        top.z_max_m + spec.plate_backing_m
    };
    let plate_r = chamber_r + spec.plate_overhang_m;
    let half_t = 0.5 * spec.wall_thickness_m;
    let seal_half_rad = (half_t / chamber_r).asin();
    // The cam ramps must sit in the right places: the arm goes in inside the
    // tank-to-blood seal at 270 degrees, and comes out only once the pocket is
    // in the tank, past the seal at 90 degrees.
    let wrap = |a: f64| (a + PI).rem_euclid(2.0 * PI) - PI;
    // The O2 rotor: the arm goes in inside the tank-to-blood seal and comes
    // out in the tank. A plunger rotor (no sites): both ramps outside the
    // seals.
    let on_seal = |a: f64| {
        let d = wrap(a - 0.5 * PI).abs().min(wrap(a - 1.5 * PI).abs());
        d < seal_half_rad + set.cam_ramp_half_rad
    };
    if set.sites.is_empty() && !set.arms.is_empty() {
        if on_seal(set.cam_inlet_rad) || on_seal(set.cam_outlet_rad) {
            return Err(PartError::InvalidGeometry(
                "a plunger ramp overlaps a seal".to_owned(),
            ));
        }
    } else if set.sites.is_empty() {
    } else if wrap(set.cam_inlet_rad - 1.5 * PI).abs() + CAM_RAMP_HALF_RAD > seal_half_rad {
        return Err(PartError::InvalidGeometry(
            "the arm-in ramp is not inside the seal".to_owned(),
        ));
    }
    if !set.sites.is_empty() && set.cam_outlet_rad - CAM_RAMP_HALF_RAD < 0.5 * PI + seal_half_rad {
        return Err(PartError::InvalidGeometry(
            "the arm-out ramp starts before the pocket is in the tank".to_owned(),
        ));
    }
    // The chamber: full height inside the face seal ring, contact height
    // (the face seal) from the ring to the rim gap. A CSG face does not land
    // on a lattice plane, so the seal is set by measured layers: the first
    // plate layer sits a whole number of (001) layers beyond the outermost
    // rotor carbon layer in the seal band. The rotor and the wall are cut
    // from the same lattice, so the gap is exact.
    let seal_inner_r = set.disc_radius_m - spec.face_seal_width_m;
    let layer = 0.25 * crate::diamond_solid::DIAMOND_LATTICE_CONSTANT_M;
    let (mut face_top, mut face_bottom) = (f64::NEG_INFINITY, f64::INFINITY);
    for i in 0..set.rotor.atom_count() {
        if set.rotor.topology.element(i) != Some(Element::CARBON) {
            continue;
        }
        if let Some(p) = set.rotor.topology.position_m(i) {
            if p[0].hypot(p[1]) >= seal_inner_r {
                face_top = face_top.max(p[2]);
                face_bottom = face_bottom.min(p[2]);
            }
        }
    }
    let layers = (spec.face_seal_gap_m / layer).round();
    let seal_top = face_top + (layers - 0.3) * layer;
    let seal_bottom = face_bottom - (layers - 0.3) * layer;
    let wall_y = plate_r + spec.wall_margin_m;
    let drain_y = plate_r + 0.5 * spec.wall_margin_m;
    let body = Union {
        a: Box::new(Box3 {
            min_m: [-half_t, -wall_y, -plate_top],
            max_m: [half_t, wall_y, plate_top],
        }),
        b: Box::new(Union {
            // The plates reach down to the face seal; the chamber cut opens
            // them again inside the seal ring.
            a: Box::new(Cylinder {
                center_m: [0.0, 0.0, 0.5 * (seal_top + plate_top)],
                radius_m: plate_r,
                height_m: plate_top - seal_top,
            }),
            b: Box::new(Cylinder {
                center_m: [0.0, 0.0, 0.5 * (seal_bottom - plate_top)],
                radius_m: plate_r,
                height_m: plate_top + seal_bottom,
            }),
        }),
    };
    let narrow_half = 0.5 * spec.drain_narrow_length_m;
    let bore = |radius_m: f64, from_m: f64, to_m: f64| -> Box<dyn Solid> {
        Box::new(Placed {
            inner: Box::new(RadialCylinder {
                azimuth_rad: 0.0,
                center_radius_m: 0.5 * (from_m + to_m),
                radius_m,
                half_length_m: 0.5 * (to_m - from_m),
                center_z_m: 0.0,
            }),
            rotation_rad: 0.0,
            offset_m: [0.0, drain_y, 0.0],
        })
    };
    let reach = half_t + 1.0e-9;
    let drain = Union {
        a: bore(spec.drain_radius_m, -reach, reach),
        b: Box::new(Union {
            a: bore(spec.drain_counterbore_radius_m, narrow_half, reach),
            b: bore(spec.drain_counterbore_radius_m, -reach, -narrow_half),
        }),
    };
    let cuts = Union {
        a: Box::new(Union {
            a: Box::new(Cylinder {
                center_m: [0.0; 3],
                radius_m: seal_inner_r,
                height_m: 2.0 * chamber_half,
            }),
            b: Box::new(Cylinder {
                center_m: [0.0, 0.0, 0.5 * (seal_top + seal_bottom)],
                radius_m: chamber_r,
                height_m: seal_top - seal_bottom,
            }),
        }),
        b: Box::new(Union {
            a: Box::new(Union {
                a: Box::new(top),
                b: Box::new(bottom),
            }),
            b: if spec.drain {
                Box::new(drain) as Box<dyn Solid>
            } else {
                // An empty box: no drain.
                Box::new(Box3 {
                    min_m: [1.0; 3],
                    max_m: [0.0; 3],
                })
            },
        }),
    };
    let solid = Difference {
        outer: Box::new(body),
        inner: Box::new(cuts),
    };
    let mut topology = capped_lattice(fill_solid(&solid))?;

    // Line the drain: each cap inside the bore becomes a hydroxyl group.
    let in_bore = |p: [f64; 3]| {
        let r = (p[1] - drain_y).hypot(p[2]);
        let x = p[0].abs();
        x < half_t + 2.0e-10
            && (r < spec.drain_radius_m
                || (x > narrow_half - 1.0e-10 && r < spec.drain_counterbore_radius_m))
    };
    let caps: Vec<(usize, usize)> = topology
        .bonds()
        .filter_map(|b| {
            let (u, v) = (b.u as usize, b.v as usize);
            match (topology.element(u), topology.element(v)) {
                (Some(Element::CARBON), Some(Element::HYDROGEN)) => Some((u, v)),
                (Some(Element::HYDROGEN), Some(Element::CARBON)) => Some((v, u)),
                _ => None,
            }
        })
        .filter(|&(_, h)| topology.position_m(h).is_some_and(in_bore))
        .collect();
    // Thin the lining: keep every n-th cap along the bore (none when 0).
    let mut caps = caps;
    caps.sort_by(|a, b| {
        let key = |h: usize| {
            let p = topology.position_m(h).unwrap_or([0.0; 3]);
            (p[0], (p[2]).atan2(p[1] - drain_y))
        };
        key(a.1)
            .partial_cmp(&key(b.1))
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let caps: Vec<(usize, usize)> = match spec.drain_hydroxyl_every {
        0 => Vec::new(),
        n => caps.into_iter().step_by(n).collect(),
    };
    let lined = hydroxylate(&mut topology, &caps)?;
    crate::pocket::charge_wall(&mut topology, crate::group::FunctionalGroup::Hydroxyl);
    Ok(HullWall {
        part: Part::new("hull_wall", topology).with_material("diamond"),
        chamber_radius_m: chamber_r,
        chamber_half_m: chamber_half,
        seal_half_rad,
        drain_y_m: drain_y,
        drain_lining: lined,
        face_seal_inner_radius_m: seal_inner_r,
    })
}

/// The geometry of a passive carrier heme site (ADR-0129).
///
/// One heme-like Fe-porphine sits in a hydrophobic pocket cut into a diamond
/// scaffold. A snug mouth passes O2. A distal hydroxyl guard hangs over the
/// bound O2 and hydrogen-bonds it, as the distal histidine does in hemoglobin,
/// so the iron resists autoxidation. The carrier has no arm and no cam: it
/// loads and releases by the oxygen partial pressure alone.
#[derive(Clone, Copy, Debug)]
pub struct CarrierSiteSpec {
    /// Radius of the ring slot that holds the porphine.
    pub slot_radius_m: f64,
    /// Half-thickness of the ring slot along z.
    pub slot_half_m: f64,
    /// Radius of the O2 pocket above the iron.
    pub pocket_radius_m: f64,
    /// Height of the pocket above the ring slot.
    pub pocket_height_m: f64,
    /// Radius of the mouth in the top face.
    pub mouth_radius_m: f64,
    /// Length of the mouth channel.
    pub mouth_length_m: f64,
    /// Diamond below the ring slot.
    pub floor_m: f64,
    /// Half-width of the scaffold slab.
    pub half_width_m: f64,
    /// Azimuth of the distal guard group, in radians.
    pub guard_azimuth_rad: f64,
    /// Whether the site shows the bound O2.
    pub loaded: bool,
    /// Radius of the proximal pocket that houses the repair donor.
    pub proximal_radius_m: f64,
    /// Depth of the proximal pocket below the ring slot.
    pub proximal_depth_m: f64,
    /// Azimuth of the distal catalyst site, in radians.
    pub catalyst_azimuth_rad: f64,
    /// Whether to build the repair donor (ADR-0108).
    pub donor: bool,
    /// Whether to build the dismutase/catalase catalyst (ADR-0139).
    pub catalyst: bool,
}

impl Default for CarrierSiteSpec {
    fn default() -> Self {
        Self {
            slot_radius_m: 7.0e-10,
            slot_half_m: 2.0e-10,
            pocket_radius_m: 3.6e-10,
            pocket_height_m: 3.5e-10,
            mouth_radius_m: 2.6e-10,
            mouth_length_m: 2.0e-10,
            floor_m: 2.4e-10,
            half_width_m: 1.2e-9,
            guard_azimuth_rad: 0.0,
            loaded: true,
            proximal_radius_m: 3.0e-10,
            proximal_depth_m: 3.0e-10,
            catalyst_azimuth_rad: PI,
            donor: true,
            catalyst: true,
        }
    }
}

/// The repair donor (ADR-0108): a proximal group that can re-reduce the iron
/// to Fe(II) after an oxidation. It sits in a pocket on the far side of the
/// ring from the O2, so it never meets a bound O2. Modelled as a proximal
/// hydroxyl (a tyrosine-like one-electron donor); the identity and the redox
/// potential are a model, per ADR-0108.
#[derive(Clone, Debug)]
pub struct DonorSite {
    /// The donor oxygen atom index inside the scaffold.
    pub oxygen: usize,
    /// The donor hydrogen atom index inside the scaffold.
    pub hydrogen: usize,
    /// Distance from the donor oxygen to the iron, in metres.
    pub fe_distance_m: f64,
}

/// The dismutase/catalase catalyst (ADR-0139): a manganese centre in the
/// pocket wall with an open site that faces the pocket. It dismutates the
/// superoxide before it can attack the ring. Modelled as a substitutional Mn
/// that keeps three wall bonds and leaves one open position; the turnover is
/// a model, per ADR-0139.
#[derive(Clone, Debug)]
pub struct CatalystSite {
    /// The manganese atom index inside the scaffold.
    pub manganese: usize,
    /// Distance from the manganese to the iron, in metres.
    pub fe_distance_m: f64,
}

/// A carrier site: the diamond scaffold (with the distal guard) and the heme.
#[derive(Clone, Debug)]
pub struct CarrierSite {
    /// The diamond scaffold, with the guard hydroxyl.
    pub scaffold: Part,
    /// The heme-like site: Fe-porphine, its tethers, and the O2 when loaded.
    pub site: Part,
    /// The iron atom index inside the site part.
    pub site_iron: usize,
    /// The (site atom, scaffold atom) pairs of the covalent tethers.
    pub site_anchor: Vec<(usize, usize)>,
    /// The tethers (ADR-0137), with scaffold indexes after the replaced
    /// hydrogens were removed.
    pub tethers: Vec<Tether>,
    /// The guard oxygen atom index inside the scaffold, when a guard was made.
    pub guard_oxygen: Option<usize>,
    /// The guard hydrogen atom index inside the scaffold.
    pub guard_hydrogen: Option<usize>,
    /// The two O2 atom indices inside the site, when loaded.
    pub o2: Option<(usize, usize)>,
    /// Radius of the O2 pocket, in metres.
    pub pocket_radius_m: f64,
    /// Radius of the mouth, in metres.
    pub mouth_radius_m: f64,
    /// The repair donor (ADR-0108), when built.
    pub donor: Option<DonorSite>,
    /// The dismutase/catalase catalyst (ADR-0139), when built.
    pub catalyst: Option<CatalystSite>,
}

/// Builds a passive carrier heme site.
#[derive(Clone, Copy, Debug, Default)]
pub struct CarrierSiteGenerator;

/// Turns the cap carbon nearest `q` into a hydroxyl whose oxygen reaches
/// toward `target`. Returns the (oxygen, hydrogen) atom indices. The hydrogen
/// sits on the C-O-H cone at 108.5 degrees, turned toward the target so it can
/// hydrogen-bond it.
fn add_guard_oh(
    top: &mut Topology,
    azimuth_rad: f64,
    target: [f64; 3],
) -> Result<Option<(usize, usize)>, PartError> {
    // A C-O-H group is fixed by C-O 1.43 A, O-H 0.96 A, and the C-O-H angle
    // 108.5 degrees. Place the oxygen so that the O-H bond points at `target`:
    // then the angle at the oxygen between C and the target is 108.5 degrees.
    // That fixes the O-target distance for a given carbon-target distance.
    let k = 1.43e-10;
    let (s108, c108) = 108.5_f64.to_radians().sin_cos();
    let mut best: Option<(f64, usize, usize, f64, f64)> = None;
    for b in top.bonds() {
        let (u, v) = (b.u as usize, b.v as usize);
        let (c, h) = match (top.element(u), top.element(v)) {
            (Some(Element::CARBON), Some(Element::HYDROGEN)) => (u, v),
            (Some(Element::HYDROGEN), Some(Element::CARBON)) => (v, u),
            _ => continue,
        };
        let Some(cp) = top.position_m(c) else {
            continue;
        };
        if cp[2] <= 0.0 {
            continue;
        }
        let az = cp[1].atan2(cp[0]);
        let daz = ((az - azimuth_rad + PI).rem_euclid(2.0 * PI) - PI).abs();
        if daz > 1.6 {
            continue;
        }
        let d = dist(cp, target);
        let disc = d * d - k * k * s108 * s108;
        if disc <= 0.0 {
            continue;
        }
        let o_target = k * c108 + disc.sqrt();
        if !(2.2e-10..=3.3e-10).contains(&o_target) {
            continue;
        }
        let score = (o_target - 2.9e-10).abs();
        if best.is_none_or(|(bs, _, _, _, _)| score < bs) {
            best = Some((score, c, h, d, o_target));
        }
    }
    let Some((_, c, h, d, o_target)) = best else {
        return Ok(None);
    };
    let cp = top
        .position_m(c)
        .ok_or_else(|| PartError::InvalidGeometry("no guard carbon".to_owned()))?;
    let u = normalize([target[0] - cp[0], target[1] - cp[1], target[2] - cp[2]]);
    let rad = normalize([cp[0], cp[1], 0.0]);
    let mut n = [
        rad[0] - (rad[0] * u[0] + rad[1] * u[1] + rad[2] * u[2]) * u[0],
        rad[1] - (rad[0] * u[0] + rad[1] * u[1] + rad[2] * u[2]) * u[1],
        rad[2] - (rad[0] * u[0] + rad[1] * u[1] + rad[2] * u[2]) * u[2],
    ];
    if norm(n) < 1.0e-12 {
        n = normalize([0.0, 0.0, 1.0]);
    } else {
        n = normalize(n);
    }
    let cos_theta = ((k - c108 * o_target) / d).clamp(-1.0, 1.0);
    let sin_theta = (1.0 - cos_theta * cos_theta).max(0.0).sqrt();
    // Bend the oxygen toward the axis, into the pocket void.
    let o = [
        cp[0] + k * (cos_theta * u[0] - sin_theta * n[0]),
        cp[1] + k * (cos_theta * u[1] - sin_theta * n[1]),
        cp[2] + k * (cos_theta * u[2] - sin_theta * n[2]),
    ];
    top.set_element(h, Element::OXYGEN)?;
    top.set_atom_type(h, "O_hydroxyl")?;
    top.set_position_m(h, o)?;
    let hd = normalize([target[0] - o[0], target[1] - o[1], target[2] - o[2]]);
    let hp = [
        o[0] + 0.96e-10 * hd[0],
        o[1] + 0.96e-10 * hd[1],
        o[2] + 0.96e-10 * hd[2],
    ];
    let hi = top.add_atom(Atom::new(Element::HYDROGEN, hp, 0.0, "H"));
    top.add_bond(Bond::new(h as u32, hi, 1, BondType::Single))?;
    Ok(Some((h, hi as usize)))
}

fn norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Places an oxygen on the cap carbon `c` (its hydrogen `h`) along the
/// direction to `target`, with the hydrogen turned away. The oxygen can then
/// bond the target. Returns the (oxygen, hydrogen) indices.
fn place_ligand_oh(
    top: &mut Topology,
    c: usize,
    h: usize,
    target: [f64; 3],
) -> Result<(usize, usize), PartError> {
    let cp = top
        .position_m(c)
        .ok_or_else(|| PartError::InvalidGeometry("no ligand carbon".to_owned()))?;
    let d = normalize([target[0] - cp[0], target[1] - cp[1], target[2] - cp[2]]);
    let o = [
        cp[0] + 1.43e-10 * d[0],
        cp[1] + 1.43e-10 * d[1],
        cp[2] + 1.43e-10 * d[2],
    ];
    top.set_element(h, Element::OXYGEN)?;
    top.set_atom_type(h, "O_hydroxyl")?;
    top.set_position_m(h, o)?;
    let mut perp = cross(d, [0.0, 0.0, 1.0]);
    if norm(perp) < 1.0e-9 {
        perp = cross(d, [1.0, 0.0, 0.0]);
    }
    perp = normalize(perp);
    let (s, co) = 108.5_f64.to_radians().sin_cos();
    let hdir = normalize([
        -co * d[0] + s * perp[0],
        -co * d[1] + s * perp[1],
        -co * d[2] + s * perp[2],
    ]);
    let hp = [
        o[0] + 0.96e-10 * hdir[0],
        o[1] + 0.96e-10 * hdir[1],
        o[2] + 0.96e-10 * hdir[2],
    ];
    let hi = top.add_atom(Atom::new(Element::HYDROGEN, hp, 0.0, "H"));
    top.add_bond(Bond::new(h as u32, hi, 1, BondType::Single))?;
    Ok((h, hi as usize))
}

/// Adds the repair donor (ADR-0108): a proximal hydroxyl on the floor of the
/// proximal pocket, with the oxygen pointing up at the iron. It sits on the
/// far side of the ring from the O2, so it never meets a bound O2.
fn add_donor_oh(top: &mut Topology) -> Result<Option<DonorSite>, PartError> {
    let fe = [0.0, 0.0, 0.0];
    let mut best: Option<(f64, usize, usize)> = None;
    for b in top.bonds() {
        let (u, v) = (b.u as usize, b.v as usize);
        let (c, h) = match (top.element(u), top.element(v)) {
            (Some(Element::CARBON), Some(Element::HYDROGEN)) => (u, v),
            (Some(Element::HYDROGEN), Some(Element::CARBON)) => (v, u),
            _ => continue,
        };
        let Some(cp) = top.position_m(c) else {
            continue;
        };
        // Below the ring, and pointing up at the iron.
        if cp[2] >= -0.5e-10 {
            continue;
        }
        let d = dist(cp, fe);
        if !(2.6e-10..=6.0e-10).contains(&d) {
            continue;
        }
        if best.is_none_or(|(bs, _, _)| d < bs) {
            best = Some((d, c, h));
        }
    }
    let Some((_, c, h)) = best else {
        return Ok(None);
    };
    let (o, hi) = place_ligand_oh(top, c, h, fe)?;
    let op = top.position_m(o).unwrap_or([0.0; 3]);
    Ok(Some(DonorSite {
        oxygen: o,
        hydrogen: hi,
        fe_distance_m: dist(op, fe),
    }))
}

/// Finds the catalyst carbon (ADR-0139): the pocket-wall cap carbon near
/// `azimuth_rad` whose open side faces the pocket. Returns the (carbon,
/// hydrogen) pair; the caller replaces the carbon with manganese and removes
/// the hydrogen, leaving the open site into the pocket.
fn find_catalyst_carbon(
    top: &Topology,
    target: [f64; 3],
    max_radius_m: f64,
) -> Option<(usize, usize)> {
    // Hydrogens per carbon, so the chosen carbon has exactly one and the
    // removed hydrogen leaves a clean open site on the manganese.
    let mut hydrogens: Vec<usize> = vec![0; top.atom_count()];
    for b in top.bonds() {
        if top.element(b.u as usize) == Some(Element::HYDROGEN) {
            hydrogens[b.v as usize] += 1;
        }
        if top.element(b.v as usize) == Some(Element::HYDROGEN) {
            hydrogens[b.u as usize] += 1;
        }
    }
    let mut best: Option<(f64, usize, usize)> = None;
    for b in top.bonds() {
        let (u, v) = (b.u as usize, b.v as usize);
        let (c, h) = match (top.element(u), top.element(v)) {
            (Some(Element::CARBON), Some(Element::HYDROGEN)) => (u, v),
            (Some(Element::HYDROGEN), Some(Element::CARBON)) => (v, u),
            _ => continue,
        };
        if hydrogens[c] != 1 {
            continue;
        }
        let Some(cp) = top.position_m(c) else {
            continue;
        };
        // On the pocket wall, not the outer slot wall (which holds the tethers).
        if cp[2] <= 0.0 || cp[0].hypot(cp[1]) > max_radius_m {
            continue;
        }
        let Some(hp) = top.position_m(h) else {
            continue;
        };
        // The hydrogen points into the pocket (toward the axis).
        let hd = normalize([hp[0] - cp[0], hp[1] - cp[1], hp[2] - cp[2]]);
        let rad = normalize([cp[0], cp[1], 0.0]);
        let inward = -(hd[0] * rad[0] + hd[1] * rad[1]);
        if inward < 0.0 {
            continue;
        }
        // Nearest to the intended wall position.
        let score = dist(cp, target);
        if best.is_none_or(|(bs, _, _)| score < bs) {
            best = Some((score, c, h));
        }
    }
    best.map(|(_, c, h)| (c, h))
}

/// Adds the repair donor to a rotor pocket (ADR-0108, ADR-0141): the nearest
/// scaffold cap carbon on the proximal side of the ring (opposite the O2)
/// becomes a hydroxyl with its oxygen pointing at the iron. The pocket axis is
/// `u`, so proximal means a negative projection on `u`. Returns `None` when no
/// carbon is in reach.
fn add_rotor_donor(
    top: &mut Topology,
    fe: [f64; 3],
    u: [f64; 3],
) -> Result<Option<DonorSite>, PartError> {
    let mut best: Option<(f64, usize, usize)> = None;
    for b in top.bonds() {
        let (a, c) = (b.u as usize, b.v as usize);
        let (c_i, h_i) = match (top.element(a), top.element(c)) {
            (Some(Element::CARBON), Some(Element::HYDROGEN)) => (a, c),
            (Some(Element::HYDROGEN), Some(Element::CARBON)) => (c, a),
            _ => continue,
        };
        let Some(cp) = top.position_m(c_i) else {
            continue;
        };
        let along = (cp[0] - fe[0]) * u[0] + (cp[1] - fe[1]) * u[1] + (cp[2] - fe[2]) * u[2];
        // Proximal: on the far side of the ring from the O2.
        if along >= -0.5e-10 {
            continue;
        }
        let d = dist(cp, fe);
        if !(2.6e-10..=ROTOR_DONOR_REACH_M).contains(&d) {
            continue;
        }
        if best.is_none_or(|(bs, _, _)| d < bs) {
            best = Some((d, c_i, h_i));
        }
    }
    let Some((_, c_i, h_i)) = best else {
        return Ok(None);
    };
    let (o, hi) = place_ligand_oh(top, c_i, h_i, fe)?;
    let op = top.position_m(o).unwrap_or([0.0; 3]);
    Ok(Some(DonorSite {
        oxygen: o,
        hydrogen: hi,
        fe_distance_m: dist(op, fe),
    }))
}

/// Finds the catalyst carbon for a rotor pocket (ADR-0139, ADR-0141): the
/// distal pocket-wall cap carbon nearest the iron, whose single hydrogen
/// points into the pocket. Returns the (carbon, hydrogen) pair; the caller
/// replaces the carbon with manganese and removes the hydrogen.
fn find_rotor_catalyst_carbon(
    top: &Topology,
    fe: [f64; 3],
    u: [f64; 3],
    max_radius_m: f64,
) -> Option<(usize, usize)> {
    let mut hydrogens: Vec<usize> = vec![0; top.atom_count()];
    for b in top.bonds() {
        if top.element(b.u as usize) == Some(Element::HYDROGEN) {
            hydrogens[b.v as usize] += 1;
        }
        if top.element(b.v as usize) == Some(Element::HYDROGEN) {
            hydrogens[b.u as usize] += 1;
        }
    }
    let mut best: Option<(f64, usize, usize)> = None;
    for b in top.bonds() {
        let (a, c) = (b.u as usize, b.v as usize);
        let (c_i, h_i) = match (top.element(a), top.element(c)) {
            (Some(Element::CARBON), Some(Element::HYDROGEN)) => (a, c),
            (Some(Element::HYDROGEN), Some(Element::CARBON)) => (c, a),
            _ => continue,
        };
        if hydrogens[c_i] != 1 {
            continue;
        }
        let Some(cp) = top.position_m(c_i) else {
            continue;
        };
        let rel = [cp[0] - fe[0], cp[1] - fe[1], cp[2] - fe[2]];
        let along = rel[0] * u[0] + rel[1] * u[1] + rel[2] * u[2];
        // Distal: on the pocket side of the ring, where the O2 binds.
        if along <= 0.0 {
            continue;
        }
        let perp = [
            rel[0] - along * u[0],
            rel[1] - along * u[1],
            rel[2] - along * u[2],
        ];
        if norm(perp) > max_radius_m {
            continue;
        }
        let Some(hp) = top.position_m(h_i) else {
            continue;
        };
        // The hydrogen points into the pocket (toward the pocket axis).
        let hd = normalize([hp[0] - cp[0], hp[1] - cp[1], hp[2] - cp[2]]);
        let rad = normalize(perp);
        if -(hd[0] * rad[0] + hd[1] * rad[1] + hd[2] * rad[2]) < 0.0 {
            continue;
        }
        let d = dist(cp, fe);
        if best.is_none_or(|(bs, _, _)| d < bs) {
            best = Some((d, c_i, h_i));
        }
    }
    best.map(|(_, c_i, h_i)| (c_i, h_i))
}

impl CarrierSiteGenerator {
    /// Builds the carrier site from an explicit spec.
    pub fn build_spec(&self, spec: &CarrierSiteSpec) -> Result<CarrierSite, PartError> {
        let z_top = spec.slot_half_m + spec.pocket_height_m + spec.mouth_length_m;
        let z_floor = -(spec.floor_m + spec.proximal_depth_m);
        let slab = Box3 {
            min_m: [-spec.half_width_m, -spec.half_width_m, z_floor],
            max_m: [spec.half_width_m, spec.half_width_m, z_top],
        };
        let slot = Cylinder {
            center_m: [0.0, 0.0, 0.0],
            radius_m: spec.slot_radius_m,
            height_m: 2.0 * spec.slot_half_m,
        };
        let pocket = Cylinder {
            center_m: [0.0, 0.0, spec.slot_half_m + 0.5 * spec.pocket_height_m],
            radius_m: spec.pocket_radius_m,
            height_m: spec.pocket_height_m,
        };
        let mouth = Cylinder {
            center_m: [
                0.0,
                0.0,
                spec.slot_half_m + spec.pocket_height_m + 0.5 * spec.mouth_length_m,
            ],
            radius_m: spec.mouth_radius_m,
            height_m: spec.mouth_length_m + 1.0e-9,
        };
        // The proximal pocket, on the far side of the ring, houses the repair
        // donor (ADR-0108).
        let proximal = Cylinder {
            center_m: [0.0, 0.0, 0.5 * z_floor],
            radius_m: spec.proximal_radius_m,
            height_m: -z_floor,
        };
        let cut: Box<dyn Solid> = Box::new(Union {
            a: Box::new(slot),
            b: Box::new(Union {
                a: Box::new(pocket),
                b: Box::new(Union {
                    a: Box::new(mouth),
                    b: Box::new(proximal),
                }),
            }),
        });
        let solid = Difference {
            outer: Box::new(slab),
            inner: cut,
        };
        let mut scaffold = capped_lattice(fill_solid(&solid))?;
        let (sn60, cs60) = 60.0_f64.to_radians().sin_cos();
        let o2 = [O_O_M * sn60, 0.0, FE_O_M + O_O_M * cs60];
        let target = if spec.loaded { o2 } else { [0.0, 0.0, FE_O_M] };
        let guard = add_guard_oh(&mut scaffold, spec.guard_azimuth_rad, target)?;
        let donor = if spec.donor {
            add_donor_oh(&mut scaffold)?
        } else {
            None
        };
        // The catalyst: a pocket-wall carbon becomes manganese, with its
        // hydrogen removed so the open site faces the pocket (ADR-0139).
        let mut catalyst_c: Option<usize> = None;
        let mut catalyst_h: Option<usize> = None;
        if spec.catalyst {
            let (sc, cc) = spec.catalyst_azimuth_rad.sin_cos();
            let cat_target = [
                0.7 * spec.pocket_radius_m * cc,
                0.7 * spec.pocket_radius_m * sc,
                spec.slot_half_m + 0.4 * spec.pocket_height_m,
            ];
            if let Some((c, h)) =
                find_catalyst_carbon(&scaffold, cat_target, spec.pocket_radius_m + 1.0e-10)
            {
                scaffold.set_element(c, Element::MANGANESE)?;
                scaffold.set_atom_type(c, "Mn")?;
                catalyst_c = Some(c);
                catalyst_h = Some(h);
            }
        }

        let u = [0.0, 0.0, 1.0];
        let t = [1.0, 0.0, 0.0];
        let (site, fe_i, mut tethers) =
            build_heme_site([0.0, 0.0, 0.0], u, t, &scaffold, spec.loaded, 0.0)?;
        // Each tether replaced a scaffold hydrogen, and the catalyst lost one:
        // remove those, and renumber the anchors, the guard, and the catalyst.
        let mut gone: Vec<usize> = tethers.iter().map(|x| x.wall_hydrogen).collect();
        if let Some(h) = catalyst_h {
            gone.push(h);
        }
        let (pruned, map) = scaffold.without_atoms(&gone);
        scaffold = pruned;
        let remap = |i: usize| -> Result<usize, PartError> {
            map.get(i)
                .copied()
                .flatten()
                .map(|v| v as usize)
                .ok_or_else(|| PartError::InvalidGeometry("a tether anchor was removed".to_owned()))
        };
        for tether in tethers.iter_mut() {
            tether.wall_atom = remap(tether.wall_atom)?;
        }
        let guard = match guard {
            Some((o, h)) => Some((remap(o)?, remap(h)?)),
            None => None,
        };
        let donor = match donor {
            Some(d) => Some(DonorSite {
                oxygen: remap(d.oxygen)?,
                hydrogen: remap(d.hydrogen)?,
                fe_distance_m: d.fe_distance_m,
            }),
            None => None,
        };
        let catalyst = match catalyst_c {
            Some(c) => {
                let mn = remap(c)?;
                let p = scaffold.position_m(mn).unwrap_or([0.0; 3]);
                Some(CatalystSite {
                    manganese: mn,
                    fe_distance_m: dist(p, [0.0, 0.0, 0.0]),
                })
            }
            None => None,
        };
        let anchors: Vec<(usize, usize)> =
            tethers.iter().map(|x| (x.site_atom, x.wall_atom)).collect();
        let o2_idx = if spec.loaded {
            let oxy: Vec<usize> = (0..site.atom_count())
                .filter(|&i| site.element(i) == Some(Element::OXYGEN))
                .collect();
            if oxy.len() == 2 {
                Some((oxy[0], oxy[1]))
            } else {
                None
            }
        } else {
            None
        };
        Ok(CarrierSite {
            scaffold: Part::new("carrier_scaffold", scaffold).with_material("diamond"),
            site: Part::new("carrier_site", site).with_material("Fe-porphine (model)"),
            site_iron: fe_i,
            site_anchor: anchors,
            tethers,
            guard_oxygen: guard.map(|g| g.0),
            guard_hydrogen: guard.map(|g| g.1),
            o2: o2_idx,
            pocket_radius_m: spec.pocket_radius_m,
            mouth_radius_m: spec.mouth_radius_m,
            donor,
            catalyst,
        })
    }

    /// Merges the scaffold and the site into one part, adding the tether bonds.
    pub fn merged(&self, spec: &CarrierSiteSpec) -> Result<Part, PartError> {
        self.build_spec(spec)?.merged_part()
    }
}

impl CarrierSite {
    /// Merges the scaffold and the site into one part, adding the tether bonds.
    pub fn merged_part(&self) -> Result<Part, PartError> {
        let mut t = self.scaffold.topology.clone();
        let base = t.atom_count();
        for i in 0..self.site.atom_count() {
            let e = self.site.topology.element(i).unwrap_or(Element::CARBON);
            let p = self.site.topology.position_m(i).unwrap_or([0.0; 3]);
            let ty = self.site.topology.atom_type(i).unwrap_or("X").to_owned();
            t.add_atom(Atom::new(e, p, 0.0, ty));
        }
        for b in self.site.topology.bonds() {
            t.add_bond(Bond::new(
                b.u + base as u32,
                b.v + base as u32,
                b.order,
                b.bond_type,
            ))?;
        }
        for (s, w) in &self.site_anchor {
            t.add_bond(Bond::new((s + base) as u32, *w as u32, 1, BondType::Single))?;
        }
        Ok(Part::new("carrier_site", t).with_material("diamond + Fe-porphine (model)"))
    }

    /// The iron index inside the merged part.
    pub fn merged_iron(&self) -> usize {
        self.scaffold.atom_count() + self.site_iron
    }
}

impl PartGenerator for CarrierSiteGenerator {
    fn id(&self) -> &'static str {
        "carrier_site"
    }

    fn name(&self) -> &'static str {
        "Passive carrier heme site"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        CARRIER_SITE_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        let mut spec = CarrierSiteSpec::default();
        if let Some(v) = parameters.get("pocket_radius_m") {
            spec.pocket_radius_m = v;
        }
        if let Some(v) = parameters.get("mouth_radius_m") {
            spec.mouth_radius_m = v;
        }
        if let Some(v) = parameters.get("slot_radius_m") {
            spec.slot_radius_m = v;
        }
        if let Some(v) = parameters.get("guard_azimuth_rad") {
            spec.guard_azimuth_rad = v;
        }
        self.merged(&spec)
    }
}

impl PartGenerator for HemeRotorGenerator {
    fn id(&self) -> &'static str {
        "heme_rotor"
    }

    fn name(&self) -> &'static str {
        "Heme-site oxygen rotor"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        HEME_ROTOR_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        Ok(self.build(parameters)?.rotor)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn set() -> HemeRotorSet {
        HemeRotorGenerator
            .build(&ParameterSet::new())
            .expect("builds")
    }

    #[test]
    fn every_site_has_iron_with_four_nitrogens_and_the_arm_tip_where_intended() {
        let s = set();
        for (k, site) in s.sites.iter().enumerate() {
            let fe = site.topology.position_m(s.site_iron[k]).expect("fe");
            let n4 = (0..site.atom_count())
                .filter(|&i| site.topology.element(i) == Some(Element::NITROGEN))
                .filter(|&i| {
                    (dist(site.topology.position_m(i).expect("n"), fe) - 2.0e-10).abs() < 1e-13
                })
                .count();
            assert_eq!(n4, 4);
            let n = s.arms[k].topology.position_m(s.arm_tip[k]).expect("tip");
            let want = if s.loaded[k] {
                LIGAND_IN_M
            } else {
                LIGAND_OUT_M + s.release_margin_m
            };
            assert!((dist(n, fe) - want).abs() < 1e-12, "{}", dist(n, fe));
        }
    }

    #[test]
    fn an_empty_pocket_can_show_the_domed_five_coordinate_iron() {
        let s = HemeRotorGenerator
            .build(&ParameterSet::new().with("iron_dome_m", IRON_DOME_M))
            .expect("builds");
        // The dome lifts the iron out of the ring plane by IRON_DOME_M, so the
        // four iron-nitrogen bonds lengthen from 2.00 A to sqrt(2.00^2 + dome^2)
        // (M28-34). The loaded pockets stay in-plane.
        for (k, site) in s.sites.iter().enumerate() {
            let want = (2.00e-10_f64.powi(2) + s.site_dome_m[k].powi(2)).sqrt();
            let fe = site.topology.position_m(s.site_iron[k]).expect("fe");
            let ns: Vec<[f64; 3]> = (0..site.atom_count())
                .filter(|&i| site.topology.element(i) == Some(Element::NITROGEN))
                .filter_map(|i| site.topology.position_m(i))
                .collect();
            assert_eq!(ns.len(), 4, "the four porphine nitrogens");
            let center = [
                ns.iter().map(|p| p[0]).sum::<f64>() / 4.0,
                ns.iter().map(|p| p[1]).sum::<f64>() / 4.0,
                ns.iter().map(|p| p[2]).sum::<f64>() / 4.0,
            ];
            assert!((dist(fe, center) - s.site_dome_m[k]).abs() < 1e-12);
            for n in &ns {
                assert!((dist(*n, fe) - want).abs() < 1e-12, "{}", dist(*n, fe));
            }
            // The four tether meso carbons are unchanged, so the tethers take
            // the dome without a rebuild.
            assert_eq!(s.site_tethers[k].len(), 4);
        }
        assert!(s.site_dome_m.contains(&IRON_DOME_M));
        assert!(s
            .loaded
            .iter()
            .enumerate()
            .all(|(k, &l)| !l || s.site_dome_m[k] == 0.0));
    }

    #[test]
    fn the_rotor_lattice_clears_every_site_and_arm() {
        let s = set();
        let rotor: Vec<[f64; 3]> = (0..s.rotor.atom_count())
            .filter_map(|i| s.rotor.topology.position_m(i))
            .collect();
        // The designed tether atoms reach the wall (ADR-0113). Walk out from
        // each anchor carbon over the tether, stop at the porphyrin ring.
        // A site is the iron and 32 porphine atoms, then the O2 (two atoms)
        // when loaded, then the tethers.
        let mut tether = std::collections::HashSet::new();
        for (k, anchors) in s.site_anchor.iter().enumerate() {
            let ring_atoms = if s.loaded[k] { 35usize } else { 33 };
            for &(si, _) in anchors {
                let mut stack = vec![si];
                while let Some(a) = stack.pop() {
                    if !tether.insert((k, a)) {
                        continue;
                    }
                    for b in s.sites[k].topology.bonds() {
                        let nb = if b.u as usize == a {
                            b.v as usize
                        } else if b.v as usize == a {
                            b.u as usize
                        } else {
                            continue;
                        };
                        if nb >= ring_atoms {
                            stack.push(nb);
                        }
                    }
                }
            }
        }
        let mut worst = f64::INFINITY;
        for (k, part) in s.sites.iter().chain(s.arms.iter()).enumerate() {
            for i in 0..part.atom_count() {
                if k < s.sites.len() && tether.contains(&(k, i)) {
                    continue;
                }
                let p = part.topology.position_m(i).expect("p");
                for q in &rotor {
                    if (p[0] - q[0]).abs() < 3e-10 && (p[1] - q[1]).abs() < 3e-10 {
                        worst = worst.min(dist(p, *q));
                    }
                }
            }
        }
        assert!(worst > 1.6e-10, "closest rotor atom {worst:e} m");
    }

    /// The arm cam over a full turn (ADR-0116). The Fe-N bond pulls each arm
    /// toward its iron, so the arm rests at the furthest position the groove
    /// allows toward the iron. On the loading arc that position must let the
    /// nitrogen reach 2.1 A; on the return arc it must keep the nitrogen at
    /// 4.0 A or more. The pin must fit the groove everywhere.
    #[test]
    fn the_arm_cam_holds_every_arm_over_a_turn() {
        let s = set();
        let stroke = s.pin_loaded_radius_m - s.pin_released_radius_m;
        let want = LIGAND_OUT_M + s.release_margin_m - LIGAND_IN_M;
        assert!((stroke - want).abs() < 2e-12, "{stroke:e}");
        let wall = hull_wall(&s, &HullWallSpec::default()).expect("wall");
        let chamber_half = wall.chamber_half_m;
        let housing = wall.part;
        let contact = 2.2e-10;
        let step = 0.05e-10;
        // One loaded arm per side, in its pocket frame (radial, tangential,
        // z) less its tip radius. Keep only the part above the disc face.
        for k in [0usize, 1] {
            assert!(s.loaded[k]);
            let side = s.pin_side[k];
            let (sn, cs) = s.azimuth_rad[k].sin_cos();
            let tip_r = s.iron_radius_m - LIGAND_IN_M;
            let local: Vec<[f64; 3]> = (0..s.arms[k].atom_count())
                .filter_map(|i| s.arms[k].topology.position_m(i))
                .filter(|p| side * p[2] > s.half_thickness_m)
                .map(|p| [p[0] * cs + p[1] * sn - tip_r, -p[0] * sn + p[1] * cs, p[2]])
                .collect();
            let wall: Vec<[f64; 3]> = (0..housing.atom_count())
                .filter_map(|i| housing.topology.position_m(i))
                .filter(|q| side * q[2] > chamber_half - 5e-10)
                .collect();
            let (mut worst_load, mut worst_release) = (f64::INFINITY, f64::NEG_INFINITY);
            let mut ramp_worst = f64::INFINITY;
            for deg in (0..360).step_by(2) {
                let th = (deg as f64).to_radians();
                let f = ArmCamTrack::loaded_fraction(th);
                if f != 0.0 && f != 1.0 {
                    // On a ramp the arm only has to fit the groove.
                    let ligand =
                        LIGAND_IN_M + (1.0 - f) * (LIGAND_OUT_M + s.release_margin_m - LIGAND_IN_M);
                    let (sn, cs) = th.sin_cos();
                    let base = s.iron_radius_m - ligand;
                    let mut nearest = f64::INFINITY;
                    for l in &local {
                        let r = base + l[0];
                        let m = [r * cs - l[1] * sn, r * sn + l[1] * cs, l[2]];
                        for q in &wall {
                            if (m[0] - q[0]).abs() < 3e-10 && (m[1] - q[1]).abs() < 3e-10 {
                                nearest = nearest.min(dist(m, *q));
                            }
                        }
                    }
                    ramp_worst = ramp_worst.min(nearest);
                    continue;
                }
                let ligand = if f == 1.0 {
                    LIGAND_IN_M
                } else {
                    LIGAND_OUT_M + s.release_margin_m
                };
                let (sn, cs) = th.sin_cos();
                let base = s.iron_radius_m - ligand;
                let near: Vec<[f64; 3]> = wall
                    .iter()
                    .copied()
                    .filter(|q| {
                        let (x, y) = (base * cs, base * sn);
                        (q[0] - x).hypot(q[1] - y) < 2.5e-9
                    })
                    .collect();
                // Free at a shift `d` toward the iron.
                let free = |d: f64| {
                    local.iter().all(|l| {
                        let r = base + d + l[0];
                        let m = [r * cs - l[1] * sn, r * sn + l[1] * cs, l[2]];
                        near.iter().all(|q| dist(m, *q) >= contact)
                    })
                };
                assert!(free(0.0), "side {side}: the pin does not fit at {deg} deg");
                let mut up = 0.0;
                while up < 4e-10 && free(up + step) {
                    up += step;
                }
                // The nitrogen rests at `ligand - up` from the iron.
                if f == 1.0 {
                    worst_load = worst_load.min(up);
                } else {
                    worst_release = worst_release.max(up);
                }
            }
            eprintln!(
                "side {side}: loading arc free travel toward the iron >= {:.2} A; return arc nitrogen >= {:.2} A from the iron",
                worst_load * 1e10,
                (LIGAND_OUT_M + s.release_margin_m - worst_release) * 1e10
            );
            eprintln!(
                "side {side}: ramps, nearest wall atom {:.2} A",
                ramp_worst * 1e10
            );
            assert!(
                ramp_worst > 1.9e-10,
                "side {side}: the arm hits a ramp at {ramp_worst:e}"
            );
            assert!(worst_load >= 0.0);
            assert!(
                LIGAND_OUT_M + s.release_margin_m - worst_release >= LIGAND_OUT_M - 1e-12,
                "side {side}: the arm reaches {:e} m from the iron on the return arc",
                LIGAND_OUT_M + s.release_margin_m - worst_release
            );
        }
    }

    /// The hull wall (ADR-0117): the rim is closed only on the two seal arcs,
    /// open to the blood (+x) and the tank (-x) elsewhere; the drain bore is
    /// open and hydroxyl-lined.
    #[test]
    fn the_hull_wall_seals_the_rim_only_where_it_crosses() {
        let s = set();
        let wall = hull_wall(&s, &HullWallSpec::default()).expect("wall");
        let atoms: Vec<[f64; 3]> = (0..wall.part.atom_count())
            .filter_map(|i| wall.part.topology.position_m(i))
            .collect();
        // Wall atoms within 1 nm outside the chamber, at the rim mid-plane.
        let covered = |theta: f64| {
            let r = wall.chamber_radius_m + 0.3e-9;
            let p = [r * theta.cos(), r * theta.sin(), 0.0];
            atoms.iter().any(|q| dist(p, *q) < 2.0e-10)
        };
        for deg in [90.0_f64, 270.0] {
            assert!(covered(deg.to_radians()), "no seal at {deg}");
        }
        for deg in [0.0_f64, 45.0, 135.0, 180.0, 225.0, 315.0] {
            assert!(!covered(deg.to_radians()), "rim closed at {deg}");
        }
        assert!(
            wall.seal_half_rad.to_degrees() > 12.0,
            "{}",
            wall.seal_half_rad.to_degrees()
        );
        // The default wall has no drain bore (ADR-0122). With one, the
        // H-capped bore leaves 2.8 A to its cap hydrogens.
        assert_eq!(wall.drain_lining, 0);
        let drained = hull_wall(
            &s,
            &HullWallSpec {
                drain: true,
                ..HullWallSpec::default()
            },
        )
        .expect("wall with drain");
        let mut nearest = f64::INFINITY;
        for i in 0..drained.part.atom_count() {
            let q = drained.part.topology.position_m(i).expect("p");
            if q[0].abs() < 0.4e-9 {
                nearest = nearest.min((q[1] - drained.drain_y_m).hypot(q[2]));
            }
        }
        eprintln!(
            "seal half-angle {:.1} deg, drain bore nearest atom {:.2} A, atoms {}",
            wall.seal_half_rad.to_degrees(),
            nearest * 1e10,
            wall.part.atom_count()
        );
        assert!((2.7e-10..2.9e-10).contains(&nearest), "{nearest:e}");
        // Every pocket's state follows the cam profile.
        for (k, &az) in s.azimuth_rad.iter().enumerate() {
            let f = ArmCamTrack::loaded_fraction(az);
            assert!(f == 0.0 || f == 1.0, "pocket {k} is mid-ramp");
            assert_eq!(s.loaded[k], f == 1.0);
        }
        assert_eq!(s.loaded.iter().filter(|&&l| l).count(), 5);
    }

    /// The water-return plunger rotor (ADR-0122): flush plungers on the
    /// blood arc, withdrawn ones on the tank arc; every plunger clears the
    /// rotor and the hull wall; the pocket of a withdrawn plunger is open
    /// 0.8 nm deep.
    #[test]
    fn the_water_rotor_plungers_fit() {
        let w = water_rotor(&WaterRotorSpec::default()).expect("water rotor");
        assert_eq!(w.arms.len(), 6);
        assert_eq!(w.loaded.iter().filter(|&&f| f).count(), 3);
        // A pocket crosses the blood-to-tank seal (90 deg) flush, so it
        // carries no water into the tank, and the tank-to-blood seal (270
        // deg) withdrawn, full of sump water.
        let f = |deg: f64| {
            ArmCamTrack::fraction_with(
                deg.to_radians(),
                w.cam_inlet_rad,
                w.cam_outlet_rad,
                w.cam_ramp_half_rad,
            )
        };
        for deg in 70..=110 {
            assert_eq!(f(deg as f64), 1.0, "not flush at {deg}");
        }
        for deg in 250..=290 {
            assert_eq!(f(deg as f64), 0.0, "not withdrawn at {deg}");
        }
        let wall = hull_wall(&w, &HullWallSpec::default()).expect("wall");
        let fixed: Vec<([f64; 3], &str)> = (0..w.rotor.atom_count())
            .filter_map(|i| w.rotor.topology.position_m(i).map(|p| (p, "rotor")))
            .chain(
                (0..wall.part.atom_count())
                    .filter_map(|i| wall.part.topology.position_m(i).map(|p| (p, "wall"))),
            )
            .collect();
        for (k, arm) in w.arms.iter().enumerate() {
            let mut near = (f64::INFINITY, [0.0; 3], "");
            for i in 0..arm.atom_count() {
                let p = arm.topology.position_m(i).expect("p");
                for (q, what) in &fixed {
                    if (p[0] - q[0]).abs() < 3e-10 && (p[1] - q[1]).abs() < 3e-10 {
                        let d = dist(p, *q);
                        if d < near.0 {
                            near = (d, p, what);
                        }
                    }
                }
            }
            let (d, p, what) = near;
            eprintln!(
                "plunger {k}: nearest {what} {:.2} A at r {:.2} nm z {:.2} nm",
                d * 1e10,
                p[0].hypot(p[1]) * 1e9,
                p[2] * 1e9
            );
            assert!(d > 1.9e-10, "plunger {k} touches at {d:e}");
        }
    }

    /// Over a full turn in 2 degree steps, ramps included, each plunger
    /// (one per pin side) placed by the cam clears the hull wall, the end
    /// plates and the grooves by more than 1.9 A (ADR-0122).
    #[test]
    fn the_water_rotor_cam_holds_every_plunger_over_a_turn() {
        let w = water_rotor(&WaterRotorSpec::default()).expect("w");
        let wall = hull_wall(&w, &HullWallSpec::default()).expect("wall");
        let wall_atoms: Vec<[f64; 3]> = (0..wall.part.atom_count())
            .filter_map(|i| wall.part.topology.position_m(i))
            .filter(|q| q[2].abs() > w.half_thickness_m)
            .collect();
        let stroke = w.pin_loaded_radius_m - w.pin_released_radius_m;
        let mut worst = f64::INFINITY;
        for k in [0usize, 1] {
            let side = w.pin_side[k];
            let (sn, cs) = w.azimuth_rad[k].sin_cos();
            let tip0 = if w.loaded[k] {
                w.iron_radius_m
            } else {
                w.iron_radius_m - stroke
            };
            // The plunger above the disc face, in its pocket frame.
            let local: Vec<[f64; 3]> = (0..w.arms[k].atom_count())
                .filter_map(|i| w.arms[k].topology.position_m(i))
                .filter(|p| side * p[2] > w.half_thickness_m)
                .map(|p| [p[0] * cs + p[1] * sn - tip0, -p[0] * sn + p[1] * cs, p[2]])
                .collect();
            for deg in (0..360).step_by(2) {
                let th = (deg as f64).to_radians();
                let f = ArmCamTrack::fraction_with(
                    th,
                    w.cam_inlet_rad,
                    w.cam_outlet_rad,
                    w.cam_ramp_half_rad,
                );
                let tip = w.iron_radius_m - (1.0 - f) * stroke;
                let (s2, c2) = th.sin_cos();
                for l in &local {
                    let r = tip + l[0];
                    let m = [r * c2 - l[1] * s2, r * s2 + l[1] * c2, l[2]];
                    for q in &wall_atoms {
                        if (m[0] - q[0]).abs() < 3e-10 && (m[1] - q[1]).abs() < 3e-10 {
                            let d = dist(m, *q);
                            if d < 1.9e-10 && d < worst {
                                eprintln!(
                                    "k {k} deg {deg} f {f:.2}: {:.2} A at r {:.2} nm z {:.2} nm",
                                    d * 1e10,
                                    m[0].hypot(m[1]) * 1e9,
                                    m[2] * 1e9
                                );
                            }
                            worst = worst.min(d);
                        }
                    }
                }
            }
        }
        eprintln!("plunger cam sweep: nearest wall atom {:.2} A", worst * 1e10);
        assert!(worst > 1.9e-10, "{worst:e}");
    }
}

#[cfg(test)]
mod carrier_tests {
    use super::*;

    fn site() -> CarrierSite {
        CarrierSiteGenerator
            .build_spec(&CarrierSiteSpec::default())
            .expect("builds")
    }

    #[test]
    fn the_site_has_an_iron_with_four_nitrogens() {
        let c = site();
        let fe = c.site.topology.position_m(c.site_iron).expect("fe");
        let n4 = (0..c.site.atom_count())
            .filter(|&i| c.site.topology.element(i) == Some(Element::NITROGEN))
            .filter(|&i| {
                (dist(c.site.topology.position_m(i).expect("n"), fe) - 2.0e-10).abs() < 1e-13
            })
            .count();
        assert_eq!(n4, 4);
    }

    #[test]
    fn the_bound_o2_sits_on_the_sixth_site_at_the_scanned_distance() {
        let c = site();
        let fe = c.site.topology.position_m(c.site_iron).expect("fe");
        let (a, b) = c.o2.expect("loaded");
        let oa = c.site.topology.position_m(a).expect("a");
        let ob = c.site.topology.position_m(b).expect("b");
        assert!(
            (dist(oa, fe) - FE_O_M).abs() < 1e-12,
            "Fe-O {}",
            dist(oa, fe)
        );
        assert!((dist(oa, ob) - O_O_M).abs() < 1e-12, "O-O {}", dist(oa, ob));
    }

    #[test]
    fn the_distal_guard_reaches_the_bound_o2() {
        let c = site();
        let g = c.guard_oxygen.expect("guard");
        let go = c.scaffold.topology.position_m(g).expect("go");
        assert_eq!(c.scaffold.topology.atom_type(g), Some("O_hydroxyl"));
        let (a, b) = c.o2.expect("loaded");
        let fe = c.site.topology.position_m(c.site_iron).expect("fe");
        let pa = c.site.topology.position_m(a).expect("a");
        let pb = c.site.topology.position_m(b).expect("b");
        let oa = if dist(pa, fe) > dist(pb, fe) { pa } else { pb };
        let d = dist(go, oa);
        assert!((2.2e-10..=3.4e-10).contains(&d), "O(guard)...O(O2) {d}");
        let h = c.guard_hydrogen.expect("h");
        let gh = c.scaffold.topology.position_m(h).expect("gh");
        let dh = dist(gh, oa);
        assert!((1.5e-10..=2.6e-10).contains(&dh), "H...O {dh}");
    }

    #[test]
    fn the_mouth_passes_the_o2_molecule() {
        let c = site();
        assert!(
            c.mouth_radius_m > 2.0e-10,
            "mouth radius {} m is too tight for O2",
            c.mouth_radius_m
        );
    }

    fn degrees(t: &Topology) -> Vec<usize> {
        let mut d = vec![0usize; t.atom_count()];
        for b in t.bonds() {
            d[b.u as usize] += 1;
            d[b.v as usize] += 1;
        }
        d
    }

    #[test]
    fn no_carbon_in_the_carrier_site_has_more_than_four_bonds() {
        // ADR-0137: the old builder bonded each tether to a wall carbon
        // without removing its hydrogen, which made five-bonded carbons.
        let part = site().merged_part().expect("merge");
        let d = degrees(&part.topology);
        for (i, &n) in d.iter().enumerate() {
            if part.topology.element(i) == Some(Element::CARBON) {
                assert!(n <= 4, "carbon {i} has {n} bonds");
            }
        }
    }

    #[test]
    fn the_carrier_tethers_have_tetrahedral_like_angles() {
        let c = site();
        assert_eq!(c.tethers.len(), 4);
        let part = c.merged_part().expect("merge");
        let t = &part.topology;
        let mut nbrs = vec![Vec::new(); t.atom_count()];
        for b in t.bonds() {
            nbrs[b.u as usize].push(b.v as usize);
            nbrs[b.v as usize].push(b.u as usize);
        }
        let base = c.scaffold.atom_count();
        for tether in &c.tethers {
            assert!((1..=3).contains(&tether.bridge_atoms));
            // Walk from the anchor into the site over the bridge carbons.
            let mut prev = tether.wall_atom;
            let mut cur = base + tether.site_atom;
            for _ in 0..tether.bridge_atoms {
                let heavy: Vec<usize> = nbrs[cur]
                    .iter()
                    .copied()
                    .filter(|&k| t.element(k) != Some(Element::HYDROGEN))
                    .collect();
                let next = heavy.iter().copied().find(|&k| k != prev).expect("chain");
                let p = |k: usize| t.position_m(k).expect("p");
                let ang = bend(p(prev), p(cur), p(next)).to_degrees();
                assert!((95.0..=140.0).contains(&ang), "bridge angle {ang}");
                prev = cur;
                cur = next;
            }
            let len = dist(
                t.position_m(base + tether.site_atom).expect("s"),
                t.position_m(tether.wall_atom).expect("w"),
            );
            assert!((1.45e-10..=1.75e-10).contains(&len), "anchor bond {len}");
        }
    }

    #[test]
    fn no_rotor_carbon_gains_a_fifth_bond_from_a_tether() {
        let set = HemeRotorGenerator
            .build(&ParameterSet::new())
            .expect("builds");
        let mut d = degrees(&set.rotor.topology);
        for anchors in &set.site_anchor {
            for &(_, w) in anchors {
                d[w] += 1;
            }
        }
        for (i, &n) in d.iter().enumerate() {
            if set.rotor.topology.element(i) == Some(Element::CARBON) {
                assert!(n <= 4, "rotor carbon {i} has {n} bonds with its tethers");
            }
        }
    }

    #[test]
    fn the_repair_donor_sits_below_the_ring_out_of_the_o2_path() {
        let c = site();
        let d = c.donor.expect("donor");
        let o = c.scaffold.topology.position_m(d.oxygen).expect("donor o");
        assert!(
            (3.0e-10..=6.5e-10).contains(&d.fe_distance_m),
            "Fe-O(donor) {} m",
            d.fe_distance_m
        );
        // Proximal: below the ring, on the far side from the distal pocket.
        assert!(o[2] < 0.0, "donor oxygen is not proximal: z={}", o[2]);
        let dir = normalize(o);
        assert!(dir[2] < -0.5, "the donor is not off the distal axis");
        assert_eq!(c.scaffold.topology.atom_type(d.oxygen), Some("O_hydroxyl"));
    }

    #[test]
    fn the_catalyst_is_manganese_with_an_open_site() {
        let c = site();
        let cat = c.catalyst.expect("catalyst");
        assert_eq!(
            c.scaffold.topology.element(cat.manganese),
            Some(Element::MANGANESE)
        );
        assert_eq!(c.scaffold.topology.atom_type(cat.manganese), Some("Mn"));
        // Three wall bonds and one open position into the pocket.
        let bonds = c
            .scaffold
            .topology
            .bonds()
            .filter(|b| b.u as usize == cat.manganese || b.v as usize == cat.manganese)
            .count();
        assert_eq!(bonds, 3, "the manganese has {bonds} bonds");
        assert!(cat.fe_distance_m > 2.5e-10, "Mn too close to Fe");
    }

    #[test]
    fn the_self_sustaining_groups_keep_every_carbon_four_bonded() {
        let part = site().merged_part().expect("merge");
        let d = degrees(&part.topology);
        for (i, &n) in d.iter().enumerate() {
            if part.topology.element(i) == Some(Element::CARBON) {
                assert!(n <= 4, "carbon {i} has {n} bonds");
            }
        }
    }

    #[test]
    fn every_rotor_pocket_has_a_donor_and_a_catalyst() {
        let set = HemeRotorGenerator
            .build(&ParameterSet::new())
            .expect("builds");
        assert_eq!(set.site_donor.len(), set.azimuth_rad.len());
        assert_eq!(set.site_catalyst.len(), set.azimuth_rad.len());
        assert!(
            set.site_donor.iter().all(|d| d.is_some()),
            "a pocket has no repair donor"
        );
        assert!(
            set.site_catalyst.iter().all(|c| c.is_some()),
            "a pocket has no catalyst"
        );
    }

    #[test]
    fn the_rotor_repair_donor_is_proximal_and_in_reach() {
        // ADR-0141: the arm collar occupies the proximal axis, so the donor
        // sits on the collar wall, further from the iron than the free site.
        let set = HemeRotorGenerator
            .build(&ParameterSet::new())
            .expect("builds");
        for (k, donor) in set.site_donor.iter().enumerate() {
            let d = donor.as_ref().expect("donor");
            assert!(
                (3.0e-10..=ROTOR_DONOR_REACH_M).contains(&d.fe_distance_m),
                "Fe-O(donor) {} m is out of reach",
                d.fe_distance_m
            );
            let az = set.azimuth_rad[k];
            let (s, c) = az.sin_cos();
            let u = [c, s, 0.0];
            let fe = set.sites[k]
                .topology
                .position_m(set.site_iron[k])
                .expect("fe");
            let o = set.rotor.topology.position_m(d.oxygen).expect("donor o");
            let along = (o[0] - fe[0]) * u[0] + (o[1] - fe[1]) * u[1] + (o[2] - fe[2]) * u[2];
            assert!(along < 0.0, "the donor is not proximal: {along} m");
            assert_eq!(set.rotor.topology.atom_type(d.oxygen), Some("O_hydroxyl"));
        }
    }

    #[test]
    fn the_rotor_catalyst_is_manganese_with_an_open_site() {
        let set = HemeRotorGenerator
            .build(&ParameterSet::new())
            .expect("builds");
        for cat in set.site_catalyst.iter().flatten() {
            assert_eq!(
                set.rotor.topology.element(cat.manganese),
                Some(Element::MANGANESE)
            );
            assert_eq!(set.rotor.topology.atom_type(cat.manganese), Some("Mn"));
            // Three wall bonds and one open position into the pocket.
            let bonds = set
                .rotor
                .topology
                .bonds()
                .filter(|b| b.u as usize == cat.manganese || b.v as usize == cat.manganese)
                .count();
            assert_eq!(bonds, 3, "the manganese has {bonds} bonds");
            assert!(cat.fe_distance_m > 2.5e-10, "Mn too close to Fe");
        }
    }

    #[test]
    fn the_rotor_repair_groups_can_be_disabled() {
        let set = HemeRotorGenerator
            .build(
                &ParameterSet::new()
                    .with("repair_donor", 0.0)
                    .with("catalyst", 0.0),
            )
            .expect("builds");
        assert!(set.site_donor.iter().all(|d| d.is_none()));
        assert!(set.site_catalyst.iter().all(|c| c.is_none()));
    }

    #[test]
    fn the_registry_builds_the_carrier_site() {
        let part =
            crate::registry::generate("carrier_site", &ParameterSet::new()).expect("generate");
        assert_eq!(part.name, "carrier_site");
    }
}
