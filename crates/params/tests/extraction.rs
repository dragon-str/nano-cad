//! The three extractions fill the PartRecord fields of `PARAMETERS.md`.
//!
//! The systems are synthetic. The test checks the dimension, the uncertainty,
//! and the provenance contract, not a physical value.

mod common;

use common::good_record;
use nanocad_engine::{BondStretchTerm, System};
use nanocad_params::{
    check, extract_failure_stress, extract_stiffness, summarize_friction, BondLimit, FailureConfig,
    StiffnessConfig,
};

const CARBON_MASS_KG: f64 = 1.992_646_879_92e-26;
const K_N_PER_M: f64 = 300.0;
const R0_M: f64 = 1.5e-10;
const AREA_M2: f64 = 1.0e-20;
const BOND_COUNT: usize = 20;
const ATOM_COUNT: usize = BOND_COUNT + 1;

fn chain() -> (System, Vec<f64>) {
    let mut system =
        System::with_masses_kg(ATOM_COUNT, &[CARBON_MASS_KG; ATOM_COUNT]).expect("valid masses");
    let mut bonds = BondStretchTerm::new();
    for atom in 0..BOND_COUNT {
        bonds
            .add_bond(atom as u32, atom as u32 + 1, K_N_PER_M, R0_M)
            .expect("valid bond");
    }
    system.set_bond_stretch(bonds).expect("bond term fits");
    let mut positions_m = Vec::with_capacity(3 * ATOM_COUNT);
    for atom in 0..ATOM_COUNT {
        positions_m.push(atom as f64 * R0_M);
        positions_m.push(0.0);
        positions_m.push(0.0);
    }
    (system, positions_m)
}

fn bond_limits(rupture_strain: f64) -> Vec<BondLimit> {
    (0..BOND_COUNT)
        .map(|atom| BondLimit {
            u: atom as u32,
            v: atom as u32 + 1,
            r0_m: R0_M,
            rupture_strain,
        })
        .collect()
}

#[test]
fn extracted_quantities_fill_a_consistent_part_record() {
    let stiffness_config = StiffnessConfig {
        reference_length_m: BOND_COUNT as f64 * R0_M,
        cross_section_area_m2: AREA_M2,
        ..StiffnessConfig::default()
    };
    let (mut system, positions_m) = chain();
    let stiffness =
        extract_stiffness(&mut system, &positions_m, &stiffness_config).expect("stiffness");

    let failure_config = FailureConfig {
        reference_length_m: BOND_COUNT as f64 * R0_M,
        cross_section_area_m2: AREA_M2,
        ..FailureConfig::default()
    };
    let (mut system, positions_m) = chain();
    let failure = extract_failure_stress(
        &mut system,
        &positions_m,
        &failure_config,
        &bond_limits(0.05),
    )
    .expect("failure");

    let lateral_forces_n: Vec<f64> = (0..500)
        .map(|index| 2.0e-11 * (2.0 * std::f64::consts::PI * index as f64 / 500.0).sin())
        .collect();
    let friction = summarize_friction(&lateral_forces_n, 1.0e-10).expect("friction");

    assert!(stiffness.elastic_modulus_pa.uncertainty_si().is_some());
    assert!(failure.failure_stress_pa.uncertainty_si().is_some());
    assert!(friction.friction_coefficient.uncertainty_si().is_some());

    let mut record = good_record();
    record.elastic_modulus_pa = stiffness.elastic_modulus_pa.clone();
    record.failure_stress_pa = failure.failure_stress_pa.clone();
    record.friction_coefficient = friction.friction_coefficient.clone();

    let problems = check(&record);
    assert!(problems.is_empty(), "unexpected problems: {problems:?}");
    assert_eq!(record.elastic_modulus_pa.unit(), "Pa");
    assert_eq!(record.failure_stress_pa.unit(), "Pa");
    assert_eq!(record.friction_coefficient.unit(), "1");
}
