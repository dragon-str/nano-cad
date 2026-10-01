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

static HEME_ROTOR_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new("disc_radius_m", Some(Unit::Metre), 4.0e-9, 2.5e-9, 2.0e-8, false, "rim radius of the rotor disc"),
    ParameterSpec::new("thickness_m", Some(Unit::Metre), 2.4e-9, 1.6e-9, 1.0e-8, false, "axial thickness of the disc"),
    ParameterSpec::new("pocket_count", None, 12.0, 1.0, 64.0, true, "number of heme pockets on the rim"),
    ParameterSpec::new("bore_radius_m", Some(Unit::Metre), 5.6e-10, 4.0e-10, 2.0e-9, false, "radius of the bore from the ring face to the rim (carbon centres)"),
    ParameterSpec::new("bore_length_m", Some(Unit::Metre), 3.0e-10, 0.0, 2.0e-9, false, "length of the bore above the ring slot"),
    ParameterSpec::new("slot_radius_m", Some(Unit::Metre), 9.0e-10, 7.0e-10, 2.0e-9, false, "radius of the slot that holds the ring"),
    ParameterSpec::new("slot_half_m", Some(Unit::Metre), 3.6e-10, 3.0e-10, 1.0e-9, false, "half-thickness of the ring slot along the radius"),
    ParameterSpec::new("arm_radius_m", Some(Unit::Metre), 2.6e-10, 1.5e-10, 1.0e-9, false, "radius of the ligand arm lattice"),
    ParameterSpec::new("arm_clearance_m", Some(Unit::Metre), 5.5e-10, 3.0e-10, 1.0e-9, false, "carbon gap between the arm and its channel"),
    ParameterSpec::new("arm_length_m", Some(Unit::Metre), 1.4e-9, 5.0e-10, 3.0e-9, false, "length of the ligand arm"),
    ParameterSpec::new("loaded_pockets", None, 6.0, 0.0, 64.0, true, "pockets 0..n-1 show the loaded state (ligand in, O2 bound); the rest show release (ligand out, empty)"),
    ParameterSpec::new("cam_pin_offset_m", Some(Unit::Metre), 3.0e-10, 2.6e-10, 1.0e-9, false, "distance of the cam pin axis from the back end of the arm"),
    ParameterSpec::new("cam_pin_top_m", Some(Unit::Metre), 2.4e-9, 1.4e-9, 5.0e-9, false, "height of the cam pin end above the disc mid-plane"),
    ParameterSpec::new("cam_release_margin_m", Some(Unit::Metre), 2.0e-10, 0.0, 1.0e-9, false, "extra withdrawal of the release track beyond the 4.0 A design, taken up by the lattice steps of the groove wall"),
];

/// The carbon-centre gap between the cam pin and a groove wall: the two
/// hydrogen caps at van der Waals contact (1.09 + 2.4 + 1.09 A, rounded up to
/// the running clearance of the module).
pub const CAM_CONTACT_M: f64 = 4.6e-10;
/// The azimuth where the cam pushes the arm in (the inlet), in radians.
pub const CAM_INLET_RAD: f64 = 345.0 * PI / 180.0;
/// The azimuth where the cam pulls the arm out (the outlet), in radians.
pub const CAM_OUTLET_RAD: f64 = 165.0 * PI / 180.0;
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
}

impl ArmCamTrack {
    /// The fraction of the stroke pushed in at azimuth `theta_rad`: 1 on the
    /// loading arc from the inlet to the outlet (rotor turning +z), 0 on the
    /// return arc, with a half-cosine ramp across each end.
    pub fn loaded_fraction(theta_rad: f64) -> f64 {
        let wrap = |a: f64| (a + PI).rem_euclid(2.0 * PI) - PI;
        let d_in = wrap(theta_rad - CAM_INLET_RAD);
        if d_in.abs() < CAM_RAMP_HALF_RAD {
            return 0.5 - 0.5 * (PI * (d_in + CAM_RAMP_HALF_RAD) / (2.0 * CAM_RAMP_HALF_RAD)).cos();
        }
        let d_out = wrap(theta_rad - CAM_OUTLET_RAD);
        if d_out.abs() < CAM_RAMP_HALF_RAD {
            return 0.5
                + 0.5 * (PI * (d_out + CAM_RAMP_HALF_RAD) / (2.0 * CAM_RAMP_HALF_RAD)).cos();
        }
        let from_inlet = (theta_rad - CAM_INLET_RAD).rem_euclid(2.0 * PI);
        let span = (CAM_OUTLET_RAD - CAM_INLET_RAD).rem_euclid(2.0 * PI);
        if from_inlet < span {
            1.0
        } else {
            0.0
        }
    }

    /// The centreline radius at azimuth `theta_rad`, in metres.
    pub fn centre_radius_m(&self, theta_rad: f64) -> f64 {
        let f = Self::loaded_fraction(theta_rad);
        self.released_radius_m + f * (self.loaded_radius_m - self.released_radius_m)
    }
}

impl Solid for ArmCamTrack {
    fn contains_m(&self, p: [f64; 3]) -> bool {
        if p[2] < self.z_min_m || p[2] > self.z_max_m {
            return false;
        }
        let d = p[0].hypot(p[1]) - self.centre_radius_m(p[1].atan2(p[0]));
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
    /// The azimuth of each pocket, in radians.
    pub azimuth_rad: Vec<f64>,
    /// Whether each pocket shows the loaded state.
    pub loaded: Vec<bool>,
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
        let loaded_n = r.require("loaded_pockets")?.round() as usize;
        let release_margin_m = r.require("cam_release_margin_m")?;
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

        // The disc and its cuts.
        let mut cuts: Option<Box<dyn Solid>> = None;
        let azimuths: Vec<f64> = (0..count)
            .map(|k| 2.0 * PI * k as f64 / count as f64)
            .collect();
        for &az in &azimuths {
            let outer = disc_m + 1.0e-9;
            let pieces: [Box<dyn Solid>; 3] = [
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
                Box::new(RadialCylinder {
                    azimuth_rad: az,
                    center_radius_m: iron_r_m - 0.5 * chan_len_m,
                    radius_m: arm_chan_r_m,
                    half_length_m: 0.5 * chan_len_m,
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
        let rotor = capped_lattice(fill_solid(solid.as_ref()))?;

        let ring = porphine_plane();
        let mut sites = Vec::with_capacity(count);
        let mut arms = Vec::with_capacity(count);
        let mut site_iron = Vec::with_capacity(count);
        let mut site_anchor: Vec<Vec<(usize, usize)>> = Vec::with_capacity(count);
        let mut loaded = Vec::with_capacity(count);
        let mut pin_side = Vec::with_capacity(count);
        for (k, &az) in azimuths.iter().enumerate() {
            let (s, c) = az.sin_cos();
            let u = [c, s, 0.0];
            let t = [-s, c, 0.0];
            let fe = [iron_r_m * c, iron_r_m * s, 0.0];
            let is_loaded = k < loaded_n;
            loaded.push(is_loaded);
            let mut site = Topology::new();
            let fe_i = site.add_atom(Atom::new(Element::IRON, fe, 0.0, "Fe"));
            // The porphine, less its four meso hydrogens: each meso carbon
            // takes a two-carbon tether to the slot wall instead (below).
            let mut meso_c: Vec<usize> = Vec::new();
            for (idx, (e, [a, b])) in ring.iter().enumerate() {
                if idx % 9 == 8 {
                    continue;
                }
                let p = [fe[0] + a * t[0], fe[1] + a * t[1], *b];
                let ai = site.add_atom(Atom::new(*e, p, 0.0, e.symbol().unwrap_or("X")));
                if idx % 9 == 7 {
                    meso_c.push(ai as usize);
                }
            }
            if is_loaded {
                let o1 = [fe[0] + FE_O_M * u[0], fe[1] + FE_O_M * u[1], 0.0];
                let (sn, cs) = 60.0_f64.to_radians().sin_cos();
                let o2 = [
                    o1[0] + O_O_M * (cs * u[0] + sn * t[0]),
                    o1[1] + O_O_M * (cs * u[1] + sn * t[1]),
                    0.0,
                ];
                site.add_atom(Atom::new(Element::OXYGEN, o1, 0.0, "O_O2"));
                site.add_atom(Atom::new(Element::OXYGEN, o2, 0.0, "O_O2"));
            }
            bond_by_distance(&mut site)?;
            // A two-carbon tether at each meso carbon, to the nearest rotor
            // wall carbon. The chain is a symmetric zigzag in the plane of the
            // radial and the pocket axis, so its hydrogens stay clear of the
            // ring. The end carbon is the anchor, bonded to the rotor atom when
            // the unit is assembled (ADR-0113).
            let mut anchors: Vec<(usize, usize)> = Vec::new();
            for &mi in &meso_c {
                let m = site.position_m(mi).unwrap_or([0.0; 3]);
                let mut best = (f64::INFINITY, [0.0; 3], 0usize);
                for j in 0..rotor.atom_count() {
                    if rotor.element(j) != Some(Element::CARBON) {
                        continue;
                    }
                    let q = rotor.position_m(j).unwrap_or([0.0; 3]);
                    let d = dist(m, q);
                    if d < best.0 {
                        best = (d, q, j);
                    }
                }
                let (l, w, wj) = best;
                let b = 1.54e-10;
                let mut nb = (l / b).ceil() as usize;
                if nb < 2 {
                    nb = 2;
                }
                if nb % 2 == 1 {
                    nb += 1;
                }
                if nb > 6 {
                    continue;
                }
                let uu = normalize([w[0] - m[0], w[1] - m[1], w[2] - m[2]]);
                let dot = uu[2];
                let nn = normalize([-dot * uu[0], -dot * uu[1], 1.0 - dot * dot]);
                let cos_a = (l / (nb as f64 * b)).clamp(-1.0, 1.0);
                let sin_a = (1.0 - cos_a * cos_a).sqrt();
                let (mut axial, mut perp) = (0.0_f64, 0.0_f64);
                let mut chain: Vec<usize> = Vec::new();
                for i in 1..nb {
                    axial += b * cos_a;
                    perp += if i % 2 == 1 { b * sin_a } else { -b * sin_a };
                    let cp = [
                        m[0] + axial * uu[0] + perp * nn[0],
                        m[1] + axial * uu[1] + perp * nn[1],
                        m[2] + axial * uu[2] + perp * nn[2],
                    ];
                    chain.push(site.add_atom(Atom::new(Element::CARBON, cp, 0.0, "C")) as usize);
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
                        w
                    };
                    let a1 = normalize([n1[0] - p[0], n1[1] - p[1], n1[2] - p[2]]);
                    let a2 = normalize([n2[0] - p[0], n2[1] - p[1], n2[2] - p[2]]);
                    let s = normalize([-(a1[0] + a2[0]), -(a1[1] + a2[1]), -(a1[2] + a2[2])]);
                    let q = normalize(cross(a1, a2));
                    for sign in [1.0, -1.0] {
                        let hdir =
                            normalize([s[0] + sign * q[0], s[1] + sign * q[1], s[2] + sign * q[2]]);
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
                anchors.push((*chain.last().unwrap_or(&mi), wj));
            }
            site_anchor.push(anchors);
            site_iron.push(fe_i as usize);
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

        Ok(HemeRotorSet {
            rotor: Part::new("heme_rotor", rotor).with_material("diamond"),
            sites,
            arms,
            site_iron,
            site_anchor,
            arm_tip: vec![cap; count],
            azimuth_rad: azimuths,
            loaded,
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
        })
    }
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

/// Builds the module housing: a ring around the rotor with two end plates and
/// one straight channel through the ring at `channel_angle_rad`, which opens
/// the ring on both sides (inlet and outlet, 180 degrees apart), plus an
/// optional `extra_cut` (the arm-cam grooves in the end plates). Loose
/// lattice carbons are removed before capping, so no two caps share a site.
pub fn module_housing(
    chamber_radius_m: f64,
    wall_m: f64,
    inner_half_height_m: f64,
    plate_m: f64,
    channel_width_m: f64,
    channel_angle_rad: f64,
    extra_cut: Option<Box<dyn Solid>>,
) -> Result<Part, PartError> {
    let outer_r = chamber_radius_m + wall_m;
    let shell = Cylinder {
        center_m: [0.0; 3],
        radius_m: outer_r,
        height_m: 2.0 * (inner_half_height_m + plate_m),
    };
    let chamber = Cylinder {
        center_m: [0.0; 3],
        radius_m: chamber_radius_m,
        height_m: 2.0 * inner_half_height_m,
    };
    let reach = outer_r + 1.0e-9;
    let channel = Placed {
        inner: Box::new(Box3 {
            min_m: [-reach, -0.5 * channel_width_m, -inner_half_height_m],
            max_m: [reach, 0.5 * channel_width_m, inner_half_height_m],
        }),
        rotation_rad: channel_angle_rad,
        offset_m: [0.0; 3],
    };
    let mut inner: Box<dyn Solid> = Box::new(Union {
        a: Box::new(chamber),
        b: Box::new(channel),
    });
    if let Some(cut) = extra_cut {
        inner = Box::new(Union { a: inner, b: cut });
    }
    let solid = Difference {
        outer: Box::new(shell),
        inner,
    };
    let topology = capped_lattice(fill_solid(&solid))?;
    Ok(Part::new("module_housing", topology).with_material("diamond"))
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
    fn the_rotor_lattice_clears_every_site_and_arm() {
        let s = set();
        let rotor: Vec<[f64; 3]> = (0..s.rotor.atom_count())
            .filter_map(|i| s.rotor.topology.position_m(i))
            .collect();
        // The designed tether atoms reach the wall (ADR-0113). Walk out from
        // each anchor carbon over the tether, stop at the porphyrin ring.
        let ring_atoms = 36usize;
        let mut tether = std::collections::HashSet::new();
        for (k, anchors) in s.site_anchor.iter().enumerate() {
            for &(si, _) in anchors {
                let mut stack = vec![si];
                while let Some(a) = stack.pop() {
                    if !tether.insert(a) {
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
                if k < s.sites.len() && tether.contains(&i) {
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
        // The module geometry of `nanocad_jigs::module_scene`: rim gap 4.6 A,
        // face gap 6.0 A, 0.5 nm of plate behind each groove.
        let chamber_half = s.half_thickness_m + 6.0e-10;
        let [top, bottom] = s.cam_tracks(chamber_half);
        let housing = module_housing(
            s.disc_radius_m + 4.6e-10,
            0.8e-9,
            chamber_half,
            top.z_max_m - chamber_half + 0.5e-9,
            1.6e-9,
            -15.0_f64.to_radians(),
            Some(Box::new(Union {
                a: Box::new(top.clone()),
                b: Box::new(bottom),
            })),
        )
        .expect("housing");
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
            for deg in (0..360).step_by(2) {
                let th = (deg as f64).to_radians();
                let f = ArmCamTrack::loaded_fraction(th);
                if f != 0.0 && f != 1.0 {
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
            assert!(worst_load >= 0.0);
            assert!(
                LIGAND_OUT_M + s.release_margin_m - worst_release >= LIGAND_OUT_M - 1e-12,
                "side {side}: the arm reaches {:e} m from the iron on the return arc",
                LIGAND_OUT_M + s.release_margin_m - worst_release
            );
        }
    }
}
