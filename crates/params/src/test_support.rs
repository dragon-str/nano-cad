//! Helpers shared by extraction tests. Compiled only for tests.

use nanocad_engine::{BondStretchTerm, Cutoff, PeriodicBox, System, VanDerWaalsTerm, VdwParams};

/// The mass of one carbon-12 atom in kilograms.
pub(crate) const CARBON_MASS_KG: f64 = 1.992_646_879_92e-26;

/// The harmonic bond force constant of the test chain, in newtons per metre.
pub(crate) const K_N_PER_M: f64 = 300.0;

/// The equilibrium bond length of the test chain, in metres.
pub(crate) const R0_M: f64 = 1.5e-10;

/// The cross-sectional area used in the stiffness and failure tests, in square
/// metres. It is a stated model choice, not a measured value.
pub(crate) const AREA_M2: f64 = 1.0e-20;

/// The number of bonds in the test chain.
pub(crate) const BOND_COUNT: usize = 20;

/// The atom count of the test chain.
pub(crate) const CHAIN_ATOM_COUNT: usize = BOND_COUNT + 1;

/// The analytic elastic modulus of the test chain at the stated area.
///
/// A uniform harmonic chain under a uniaxial strain stores
/// `U = (1/2) (N-1) k r0^2 eps^2`. With volume `A L` and `L = (N-1) r0`, the
/// stress is `k r0 eps / A`, so the modulus is `k r0 / A`.
pub(crate) fn true_modulus_pa() -> f64 {
    K_N_PER_M * R0_M / AREA_M2
}

/// Builds a straight harmonic chain along x with one bond per atom pair.
pub(crate) fn uniform_chain() -> (System, Vec<f64>) {
    let mut system = System::with_masses_kg(CHAIN_ATOM_COUNT, &[CARBON_MASS_KG; CHAIN_ATOM_COUNT])
        .expect("valid masses");
    let mut bonds = BondStretchTerm::new();
    for atom in 0..BOND_COUNT {
        bonds
            .add_bond(atom as u32, atom as u32 + 1, K_N_PER_M, R0_M)
            .expect("valid bond");
    }
    system.set_bond_stretch(bonds).expect("bond term fits");

    let mut positions_m = Vec::with_capacity(3 * CHAIN_ATOM_COUNT);
    for atom in 0..CHAIN_ATOM_COUNT {
        positions_m.push(atom as f64 * R0_M);
        positions_m.push(0.0);
        positions_m.push(0.0);
    }
    (system, positions_m)
}

/// The number of fixed atoms in the sliding test, spaced like a graphene row.
pub(crate) const SUBSTRATE_COUNT: usize = 41;

/// The substrate row spacing in metres.
pub(crate) const SUBSTRATE_SPACING_M: f64 = 2.46e-10;

/// The atom index of the slider in the sliding test.
pub(crate) const SLIDER_ATOM: u32 = SUBSTRATE_COUNT as u32;

/// Builds a fixed substrate row with one mobile slider atom above it.
///
/// The interaction is Buckingham van der Waals. The parameters have a minimum
/// near 3.4e-10 m with a well depth near 2.0e-21 J, so the contact is stable.
/// The slider is atom [`SLIDER_ATOM`].
pub(crate) fn sliding_system() -> (System, Vec<f64>) {
    let atom_count = SUBSTRATE_COUNT + 1;
    let params: Vec<VdwParams> = (0..atom_count)
        .map(|_| VdwParams::new(1.2e-15, 4.0e10, 5.5e-78))
        .collect();
    let cutoff = Cutoff::new(1.2e-9, 1.0e-9).expect("valid cutoff");
    let vdw = VanDerWaalsTerm::from_params(&params, cutoff, PeriodicBox::non_periodic())
        .expect("valid term");
    let mut system = System::with_masses_kg(atom_count, &vec![CARBON_MASS_KG; atom_count])
        .expect("valid masses");
    system.set_van_der_waals(vdw).expect("term fits");

    let mut positions_m = Vec::with_capacity(3 * atom_count);
    let midpoint = (SUBSTRATE_COUNT as f64 - 1.0) / 2.0;
    for atom in 0..SUBSTRATE_COUNT {
        positions_m.push((atom as f64 - midpoint) * SUBSTRATE_SPACING_M);
        positions_m.push(0.0);
        positions_m.push(0.0);
    }
    positions_m.push(0.0);
    positions_m.push(0.0);
    positions_m.push(3.8e-10);
    (system, positions_m)
}
