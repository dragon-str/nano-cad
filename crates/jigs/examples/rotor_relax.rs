//! Relaxes a representative sub-assembly of the sorting rotor and prints the
//! stability facts.
//!
//! The sub-assembly holds the whole ejection rod body (the shaft, the follower
//! pin and the leaf spring) and the whole cam hub. These two bodies meet at the
//! pin and the hub surface, and they are the parts that move against each
//! other, so this pair is the representative test of the device.
//!
//! Each body is relaxed whole. A cut-out of a solid has broken bonds at the
//! cut, and those boundary atoms relax against nothing, so a cut-out would
//! report a false collapse.
//!
//! The metric relaxes the covalent lattice to its ideal bond lengths, then
//! reports the energy drop, the largest force on one atom, the largest bond
//! strain, the number of clashes and the largest displacement. A stable
//! assembly does not explode and does not fall apart.

use nanocad_jigs::{build_rotor_scene, CAM_BODY, ROD_BODY_FIRST, ROD_COUNT};
use nanocad_meter::{relax_subassembly, RelaxAtom, RelaxTarget};
use nanocad_model::Element;

fn element_of(atomic_number: u8) -> Element {
    Element::from_atomic_number(atomic_number).unwrap_or(Element::CARBON)
}

fn main() {
    let scene = match build_rotor_scene() {
        Ok(scene) => scene,
        Err(error) => {
            eprintln!("the rotor scene failed: {error}");
            std::process::exit(1);
        }
    };
    if ROD_COUNT == 0 {
        eprintln!("the rotor has no rod");
        std::process::exit(1);
    }
    let target = RelaxTarget::default();
    let rod_body = ROD_BODY_FIRST;

    let atoms: Vec<RelaxAtom> = scene
        .atomistic
        .atoms
        .iter()
        .filter(|atom| atom.body == rod_body || atom.body == CAM_BODY)
        .map(|atom| RelaxAtom {
            position_m: atom.position_m,
            element: element_of(atom.atomic_number),
            body: atom.body,
        })
        .collect();

    let report = match relax_subassembly(&atoms, &target) {
        Some(report) => report,
        None => {
            eprintln!("the sub-assembly did not relax ({} atoms)", atoms.len());
            std::process::exit(1);
        }
    };

    println!(
        "{{\"subassembly_atoms\":{},\"subassembly_bonds\":{},\
         \"initial_energy_j\":{:e},\"final_energy_j\":{:e},\"energy_drop_j\":{:e},\
         \"peak_force_n\":{:e},\"max_bond_strain\":{:e},\
         \"initial_clash_count\":{},\"clash_count\":{},\
         \"max_displacement_m\":{:e},\"converged\":{}}}",
        report.atom_count,
        report.bond_count,
        report.initial_energy_j,
        report.final_energy_j,
        report.energy_drop_j,
        report.peak_force_n,
        report.max_bond_strain,
        report.initial_clash_count,
        report.clash_count,
        report.max_displacement_m,
        report.converged,
    );
}
