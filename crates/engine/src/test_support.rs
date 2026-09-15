//! Helpers shared by force-term tests. Compiled only for tests.

use crate::angle_bend::{AngleBendParams, AngleBendTerm};
use crate::bond_stretch::BondStretchTerm;
use crate::electrostatic::ElectrostaticTerm;
use crate::nonbonded::{Cutoff, PeriodicBox};
use crate::out_of_plane::{OutOfPlaneParams, OutOfPlaneTerm};
use crate::system::System;
use crate::torsion::{TorsionParams, TorsionTerm};
use crate::van_der_waals::{VanDerWaalsTerm, VdwParams};

/// The mass of one carbon-12 atom in kilograms.
pub(crate) const CARBON_MASS_KG: f64 = 1.992_646_879_92e-26;

/// A small deterministic xorshift random generator.
pub(crate) struct Rng(u64);

impl Rng {
    /// Builds a generator from a seed.
    pub(crate) const fn new(seed: u64) -> Self {
        Self(seed)
    }

    /// Returns the next value in `[0, 1)`.
    pub(crate) fn next_f64(&mut self) -> f64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        (x >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Returns a value in `[-scale, scale)`.
    pub(crate) fn symmetric(&mut self, scale_m: f64) -> f64 {
        (self.next_f64() - 0.5) * 2.0 * scale_m
    }
}

/// Returns the central finite-difference gradient of a scalar function.
///
/// Entry `index` is `(f(x + h) - f(x - h)) / (2 h)`. The step is in the same
/// unit as the coordinates.
pub(crate) fn central_difference<F>(mut energy: F, positions: &[f64], step: f64) -> Vec<f64>
where
    F: FnMut(&[f64]) -> f64,
{
    let mut gradient = vec![0.0; positions.len()];
    for index in 0..positions.len() {
        let mut forward = positions.to_vec();
        forward[index] += step;
        let mut backward = positions.to_vec();
        backward[index] -= step;
        gradient[index] = (energy(&forward) - energy(&backward)) / (2.0 * step);
    }
    gradient
}

const BONDS: [[u32; 2]; 5] = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5]];
const ANGLES: [[u32; 3]; 4] = [[0, 1, 2], [1, 2, 3], [2, 3, 4], [3, 4, 5]];
const TORSIONS: [[u32; 4]; 3] = [[0, 1, 2, 3], [1, 2, 3, 4], [2, 3, 4, 5]];
const IMPROPERS: [[u32; 4]; 2] = [[0, 1, 2, 3], [2, 3, 4, 5]];
const EXCLUSIONS: [[u32; 2]; 12] = [
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 4],
    [4, 5],
    [0, 2],
    [1, 3],
    [2, 4],
    [3, 5],
    [0, 3],
    [1, 4],
    [2, 5],
];

/// Returns the flat positions of a small six-atom chain molecule, in metres.
///
/// The geometry is a slightly jittered chain with about 1.5e-10 m spacings.
/// The jitter is deterministic.
pub(crate) fn small_molecule_positions_m() -> Vec<f64> {
    let mut rng = Rng::new(0x51CE);
    let base_m = [
        [0.0, 0.0, 0.0],
        [1.5e-10, 0.0, 0.0],
        [2.2e-10, 1.3e-10, 0.0],
        [3.6e-10, 1.5e-10, 0.2e-10],
        [4.3e-10, 2.8e-10, -0.3e-10],
        [5.7e-10, 3.0e-10, 0.1e-10],
    ];
    let mut positions_m = Vec::with_capacity(18);
    for atom in base_m {
        positions_m.push(atom[0] + rng.symmetric(1.0e-12));
        positions_m.push(atom[1] + rng.symmetric(1.0e-12));
        positions_m.push(atom[2] + rng.symmetric(1.0e-12));
    }
    positions_m
}

/// Builds a system with bond-stretch and angle-bend terms only.
///
/// The system has no non-bonded terms, so it needs no neighbor list. Every atom
/// takes the mass given in `masses_kg`.
pub(crate) fn bonded_system(positions_m: &[f64], masses_kg: &[f64]) -> System {
    let atom_count = positions_m.len() / 3;
    let mut system = System::with_masses_kg(atom_count, masses_kg).expect("valid masses");

    let mut bonds = BondStretchTerm::new();
    for [u, v] in BONDS {
        bonds.add_bond(u, v, 300.0, 1.5e-10).expect("valid bond");
    }
    system.set_bond_stretch(bonds).expect("bond term fits");

    let angle_params: Vec<AngleBendParams> = ANGLES
        .iter()
        .map(|_| AngleBendParams::new(5.0e-19, 1.9))
        .collect();
    let angles = AngleBendTerm::from_triples(&ANGLES, &angle_params).expect("valid angles");
    system.set_angle_bend(angles).expect("angle term fits");

    system
}

/// Builds a system with every force term active.
///
/// Bond stretch, angle bend, torsion, out-of-plane, van der Waals, and
/// electrostatic terms act at once. Bonded and 1-3 and 1-4 pairs are excluded
/// from the non-bonded terms.
pub(crate) fn full_system(positions_m: &[f64], masses_kg: &[f64]) -> System {
    let atom_count = positions_m.len() / 3;
    let mut system = bonded_system(positions_m, masses_kg);

    let torsion_params: Vec<TorsionParams> = TORSIONS
        .iter()
        .map(|_| TorsionParams::new(2.0e-20, -3.0e-20, 5.0e-20))
        .collect();
    let torsions = TorsionTerm::from_quads(&TORSIONS, &torsion_params).expect("valid torsions");
    system.set_torsion(torsions).expect("torsion term fits");

    let improper_params: Vec<OutOfPlaneParams> = IMPROPERS
        .iter()
        .map(|_| OutOfPlaneParams::new(4.0e-19, 0.0))
        .collect();
    let impropers =
        OutOfPlaneTerm::from_quads(&IMPROPERS, &improper_params).expect("valid impropers");
    system
        .set_out_of_plane(impropers)
        .expect("improper term fits");

    let cutoff = Cutoff::new(1.0e-9, 0.8e-9).expect("valid cutoff");
    let vdw_params: Vec<VdwParams> = (0..atom_count)
        .map(|atom| VdwParams::new(2.2e-19 + 1.0e-21 * atom as f64, 4.2e10, 8.0e-79))
        .collect();
    let vdw = VanDerWaalsTerm::from_params(&vdw_params, cutoff, PeriodicBox::non_periodic())
        .expect("valid van der Waals term");
    system
        .set_van_der_waals(vdw)
        .expect("van der Waals term fits");

    let charges_c: Vec<f64> = (0..atom_count)
        .map(|atom| {
            let magnitude = 0.5e-19 + 0.1e-19 * atom as f64;
            if atom % 2 == 0 {
                magnitude
            } else {
                -magnitude
            }
        })
        .collect();
    let electrostatic =
        ElectrostaticTerm::from_charges(&charges_c, cutoff, PeriodicBox::non_periodic())
            .expect("valid electrostatic term");
    system
        .set_electrostatic(electrostatic)
        .expect("electrostatic term fits");

    for [i, j] in EXCLUSIONS {
        system.add_exclusion(i, j).expect("valid exclusion");
    }
    system
}
