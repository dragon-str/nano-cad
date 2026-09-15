//! Prints a fixed periodic SPC-like water box as JSON for an OpenMM cross-check.
//!
//! The output holds every parameter, the box, the positions, the potential
//! energy, and the per-atom forces. The comparison driver in
//! `benchmarks/openmm/` rebuilds the identical model in OpenMM from these
//! fields. All values are SI: metres, joules, and newtons.

use std::error::Error;

use nanocad_engine::{
    spc_water_box, COULOMB_CONSTANT_N_M2_PER_C2, ELEMENTARY_CHARGE_C, HYDROGEN_MASS_KG,
    OXYGEN_MASS_KG, SPC_ANGLE_K_J_PER_RAD2, SPC_ANGLE_RAD, SPC_BOND_K_N_PER_M, SPC_BOND_LENGTH_M,
    SPC_HYDROGEN_CHARGE_C, SPC_OXYGEN_CHARGE_C, SPC_VDW_A_J, SPC_VDW_B_PER_M, SPC_VDW_C_J_M6,
};

const MOLECULE_COUNT: usize = 200;
const TEMPERATURE_K: f64 = 300.0;
const SEED: u64 = 0x5A17E2;

fn box_length_m(molecule_count: usize) -> f64 {
    (molecule_count as f64 * 33.4e-30).cbrt()
}

fn join(values: &[f64]) -> String {
    values
        .iter()
        .map(|value| format!("{value:e}"))
        .collect::<Vec<String>>()
        .join(",")
}

fn main() -> Result<(), Box<dyn Error>> {
    let box_length_m = box_length_m(MOLECULE_COUNT);
    let cutoff_m = 0.35 * box_length_m;
    let switch_on_m = 0.9 * cutoff_m;
    let skin_m = 0.1 * box_length_m;

    let mut water = spc_water_box(MOLECULE_COUNT, box_length_m, TEMPERATURE_K, SEED)?;
    let positions_m = water.positions_m().to_vec();
    let energy_j = water.potential_energy_j()?;
    let forces_n = water.system_mut().forces_n(&positions_m)?;

    let mut out = String::new();
    out.push_str("{\n");
    out.push_str(&format!("\"molecule_count\":{MOLECULE_COUNT},\n"));
    out.push_str(&format!("\"atom_count\":{},\n", water.atom_count()));
    out.push_str(&format!("\"seed\":{SEED},\n"));
    out.push_str(&format!("\"temperature_k\":{TEMPERATURE_K:e},\n"));
    out.push_str(&format!("\"box_length_m\":{box_length_m:e},\n"));
    out.push_str(&format!("\"cutoff_m\":{cutoff_m:e},\n"));
    out.push_str(&format!("\"switch_on_m\":{switch_on_m:e},\n"));
    out.push_str(&format!("\"skin_m\":{skin_m:e},\n"));
    out.push_str(&format!("\"bond_k_n_per_m\":{SPC_BOND_K_N_PER_M:e},\n"));
    out.push_str(&format!("\"bond_r0_m\":{SPC_BOND_LENGTH_M:e},\n"));
    out.push_str(&format!(
        "\"angle_k_j_per_rad2\":{SPC_ANGLE_K_J_PER_RAD2:e},\n"
    ));
    out.push_str(&format!("\"angle_theta0_rad\":{SPC_ANGLE_RAD:e},\n"));
    out.push_str(&format!("\"vdw_a_j\":{SPC_VDW_A_J:e},\n"));
    out.push_str(&format!("\"vdw_b_per_m\":{SPC_VDW_B_PER_M:e},\n"));
    out.push_str(&format!("\"vdw_c_j_m6\":{SPC_VDW_C_J_M6:e},\n"));
    out.push_str(&format!("\"oxygen_charge_c\":{SPC_OXYGEN_CHARGE_C:e},\n"));
    out.push_str(&format!(
        "\"hydrogen_charge_c\":{SPC_HYDROGEN_CHARGE_C:e},\n"
    ));
    out.push_str(&format!(
        "\"elementary_charge_c\":{ELEMENTARY_CHARGE_C:e},\n"
    ));
    out.push_str(&format!(
        "\"coulomb_constant_n_m2_per_c2\":{COULOMB_CONSTANT_N_M2_PER_C2:e},\n"
    ));
    out.push_str(&format!("\"oxygen_mass_kg\":{OXYGEN_MASS_KG:e},\n"));
    out.push_str(&format!("\"hydrogen_mass_kg\":{HYDROGEN_MASS_KG:e},\n"));
    out.push_str(&format!("\"energy_j\":{energy_j:e},\n"));
    out.push_str(&format!("\"positions_m\":[{}],\n", join(&positions_m)));
    out.push_str(&format!("\"forces_n\":[{}]\n", join(&forces_n)));
    out.push_str("}\n");

    print!("{out}");
    Ok(())
}
