//! The oxygen-loading rotor of the respirocyte design (ADR-0111).
//!
//! A diamond disc carries `pocket_count` radial pockets on its rim. Each
//! pocket holds one heme-like site: an idealized Fe-porphine whose plane is
//! normal to the radius, with the iron on the pocket axis. A snug bore runs
//! from the ring face out to the rim; O2 binds on the outer face of the iron.
//! Below the ring, a radial channel carries an arm: a diamond rod whose tip is
//! an NH2 group, the axial ligand. A cam (not built) slides the arm: tip in
//! at 2.1 A from the iron loads O2, tip out at 4.0 A releases it (ADR-0106).
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
    ParameterSpec::new("arm_length_m", Some(Unit::Metre), 1.2e-9, 5.0e-10, 3.0e-9, false, "length of the ligand arm"),
    ParameterSpec::new("loaded_pockets", None, 6.0, 0.0, 64.0, true, "pockets 0..n-1 show the loaded state (ligand in, O2 bound); the rest show release (ligand out, empty)"),
];

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
        let half_m = 0.5 * thick_m;
        if slot_r_m > half_m {
            return Err(PartError::InvalidGeometry(
                "the ring slot is wider than the disc is thick".to_owned(),
            ));
        }
        // The iron plane sits so that the slot and the bore end at the rim.
        let iron_r_m = disc_m - slot_h_m - bore_len_m;
        let arm_chan_r_m = arm_r_m + arm_clr_m;
        let chan_len_m = LIGAND_OUT_M + arm_len_m + 1.0e-9;

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
        let rotor = capped_lattice(fill_solid(solid.as_ref()))?;

        // The arm is a diamond rod along the crystal [111] axis: a carbon on
        // its (111) end face has one bond straight along the axis, so the
        // C-N bond and the nitrogen lone pair face the iron. The rod is then
        // turned so [111] becomes +x.
        let a111 = [1.0 / 3f64.sqrt(); 3];
        let arm_solid = AxisRod {
            axis: a111,
            radius_m: arm_r_m,
            length_m: arm_len_m,
        };
        let mut arm = Topology::new();
        let arm_c: Vec<u32> = keep_well_bonded(fill_solid(&arm_solid))
            .into_iter()
            .map(|p| arm.add_atom(Atom::new(Element::CARBON, p, 0.0, "C")))
            .collect();
        diamond_solid::bond_and_cap(&mut arm, &arm_c)?;
        let turn = rotation_onto_x(a111);
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

        let ring = porphine_plane();
        let mut sites = Vec::with_capacity(count);
        let mut arms = Vec::with_capacity(count);
        let mut site_iron = Vec::with_capacity(count);
        let mut loaded = Vec::with_capacity(count);
        for (k, &az) in azimuths.iter().enumerate() {
            let (s, c) = az.sin_cos();
            let u = [c, s, 0.0];
            let t = [-s, c, 0.0];
            let fe = [iron_r_m * c, iron_r_m * s, 0.0];
            let is_loaded = k < loaded_n;
            loaded.push(is_loaded);
            let mut site = Topology::new();
            let fe_i = site.add_atom(Atom::new(Element::IRON, fe, 0.0, "Fe"));
            for (e, [a, b]) in &ring {
                let p = [fe[0] + a * t[0], fe[1] + a * t[1], *b];
                site.add_atom(Atom::new(*e, p, 0.0, e.symbol().unwrap_or("X")));
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
            site_iron.push(fe_i as usize);
            sites.push(
                Part::new(format!("heme_site_{k}"), site).with_material("Fe-porphine (model)"),
            );

            // The arm: +x of the arm frame points along +u; its N at the
            // ligand distance below the iron.
            let tip_r_m = iron_r_m - if is_loaded { LIGAND_IN_M } else { LIGAND_OUT_M };
            let mut placed = arm.clone();
            for i in 0..placed.atom_count() {
                if let Some(p) = placed.position_m(i) {
                    let along = tip_r_m + p[0];
                    placed.set_position_m(i, [along * c - p[1] * s, along * s + p[1] * c, p[2]])?;
                }
            }
            arms.push(Part::new(format!("ligand_arm_{k}"), placed).with_material("diamond + NH2"));
        }

        Ok(HemeRotorSet {
            rotor: Part::new("heme_rotor", rotor).with_material("diamond"),
            sites,
            arms,
            site_iron,
            arm_tip: vec![cap; count],
            azimuth_rad: azimuths,
            loaded,
            iron_radius_m: iron_r_m,
            disc_radius_m: disc_m,
            half_thickness_m: half_m,
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
    axis: [f64; 3],
    radius_m: f64,
    length_m: f64,
}

impl Solid for AxisRod {
    fn contains_m(&self, p: [f64; 3]) -> bool {
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
        Bounds {
            min_m: [e[0].min(0.0) - r, e[1].min(0.0) - r, e[2].min(0.0) - r],
            max_m: [e[0].max(0.0) + r, e[1].max(0.0) + r, e[2].max(0.0) + r],
        }
    }

    fn clone_box(&self) -> Box<dyn Solid> {
        Box::new(self.clone())
    }
}

/// The rotation that turns the unit `a` onto +x (Rodrigues).
fn rotation_onto_x(a: [f64; 3]) -> [[f64; 3]; 3] {
    let k = [0.0, a[2], -a[1]]; // a x (1,0,0)
    let s = (k[1] * k[1] + k[2] * k[2]).sqrt();
    let c = a[0];
    if s < 1e-12 {
        return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    }
    let k = [0.0, k[1] / s, k[2] / s];
    let t = 1.0 - c;
    [
        [
            t * k[0] * k[0] + c,
            t * k[0] * k[1] - s * k[2],
            t * k[0] * k[2] + s * k[1],
        ],
        [
            t * k[0] * k[1] + s * k[2],
            t * k[1] * k[1] + c,
            t * k[1] * k[2] - s * k[0],
        ],
        [
            t * k[0] * k[2] - s * k[1],
            t * k[1] * k[2] + s * k[0],
            t * k[2] * k[2] + c,
        ],
    ]
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
/// the ring on both sides (inlet and outlet, 180 degrees apart). Loose
/// lattice carbons are removed before capping, so no two caps share a site.
pub fn module_housing(
    chamber_radius_m: f64,
    wall_m: f64,
    inner_half_height_m: f64,
    plate_m: f64,
    channel_width_m: f64,
    channel_angle_rad: f64,
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
    let solid = Difference {
        outer: Box::new(shell),
        inner: Box::new(Union {
            a: Box::new(chamber),
            b: Box::new(channel),
        }),
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
                LIGAND_OUT_M
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
        let mut worst = f64::INFINITY;
        for part in s.sites.iter().chain(s.arms.iter()) {
            for i in 0..part.atom_count() {
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
}
