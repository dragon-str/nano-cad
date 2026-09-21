//! Characterise a model iron-oxygen coordination contact and print it as JSON.
//!
//! Oxygen binds a ferrous iron by charge transfer. Carbon dioxide does not form
//! that bond: it approaches a metal only by dispersion, which the classical
//! pocket term already covers. So this example gives the iron-oxygen contact a
//! well depth, and it gives carbon dioxide no coordination term.
//!
//! The well depth is assigned from the binding enthalpy of an oxygen carrier,
//! about `-40 kJ/mol`. This is an estimate, not a measurement.
//!
//! The free energy that sets `P50` is smaller than the well depth, because
//! binding removes the translational and rotational entropy of the gas. A
//! static potential well does not carry that entropy. The reconciliation of the
//! well depth with the delivery model is open. See ADR-0083.

use nanocad_engine::{MetalLigandParams, MetalLigandTerm};

/// The Boltzmann constant, in joules per kelvin.
const BOLTZMANN_J_PER_K: f64 = 1.380_649e-23;

/// The body temperature, in kelvin.
const BODY_TEMPERATURE_K: f64 = 310.15;

/// The equilibrium iron-oxygen distance, in metres.
const IRON_OXYGEN_M: f64 = 1.90e-10;

/// The iron-oxygen stretching force constant, in newtons per metre.
///
/// A heme iron-oxygen stretch is a soft mode of about one hundred and fifty
/// reciprocal centimetres. This is a model value.
const IRON_OXYGEN_K_N_PER_M: f64 = 2.0e2;

/// The iron-oxygen well depth in joules. About `-40 kJ/mol`, an estimate.
const IRON_OXYGEN_WELL_DEPTH_J: f64 = 6.642e-20;

/// The stroke of the ejector rod, in metres. It matches the rotor scene.
const EJECTION_STROKE_M: f64 = 1.5e-9;

/// Returns the maximum attractive force of one contact, in newtons.
///
/// The value is the largest `-dU/dr` over a scan outside the equilibrium
/// length. The ejector rod must exceed it to pull the guest off the metal.
fn maximum_attractive_force_n(params: MetalLigandParams) -> f64 {
    let mut term = MetalLigandTerm::new();
    term.add_pair(0, 1, params).expect("valid pair");
    let mut maximum_n = 0.0_f64;
    let steps = 400;
    for step in 0..=steps {
        let distance_m = params.equilibrium_m * (1.0 + 1.5 * step as f64 / steps as f64);
        let force_n = term.radial_force_n(0, distance_m);
        if -force_n > maximum_n {
            maximum_n = -force_n;
        }
    }
    maximum_n
}

fn main() {
    let thermal_j = BOLTZMANN_J_PER_K * BODY_TEMPERATURE_K;
    let params = MetalLigandParams::from_force_constant(
        IRON_OXYGEN_WELL_DEPTH_J,
        IRON_OXYGEN_M,
        IRON_OXYGEN_K_N_PER_M,
    );
    let maximum_force_n = maximum_attractive_force_n(params);
    let average_force_n = params.well_depth_j / EJECTION_STROKE_M;
    println!(
        "{{\"contact\":\"iron_oxygen\",\"equilibrium_m\":{:.4e},\"well_depth_j\":{:.4e},\
         \"well_depth_kj_per_mol\":{:.2},\"well_depth_over_kt\":{:.1},\
         \"force_constant_n_per_m\":{:.4e},\"maximum_attractive_force_n\":{:.4e},\
         \"average_ejection_force_n\":{:.4e},\"ejection_stroke_m\":{:.4e},\
         \"note\":\"assigned well depth, no carbon dioxide coordination term, entropy reconciliation open\"}}",
        params.equilibrium_m,
        params.well_depth_j,
        params.well_depth_j * 6.022_140_76e23 / 1000.0,
        params.well_depth_j / thermal_j,
        params.force_constant_n_per_m(),
        maximum_force_n,
        average_force_n,
        EJECTION_STROKE_M
    );
}
