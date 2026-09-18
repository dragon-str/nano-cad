//! Prints the binding pocket facts and the fit of every guest as one JSON line.
//!
//! The line reports the pocket size, the number of decorated wall sites, the
//! free cavity volume, and for each guest the clearance and the wall well
//! depth. Every value is a model value. See ADR-0063 and ADR-0064.

use std::f64::consts::PI;

use nanocad_model::Element;
use nanocad_parts::{
    guests, pocket_radius_for, wall_contact_distance_m, wall_well_depth_j, BindingPocketGenerator,
    ParameterSet, PartGenerator,
};

fn main() {
    let mut parameters = ParameterSet::new();
    if let Some(text) = std::env::args().nth(1) {
        if let Ok(radius_m) = text.parse::<f64>() {
            parameters.set("pocket_radius_m", radius_m);
        }
    }
    let pocket = match BindingPocketGenerator.generate(&parameters) {
        Ok(part) => part,
        Err(error) => {
            println!("{{\"error\":\"{error}\"}}");
            return;
        }
    };
    let radius_m = resolved_m(&parameters, "pocket_radius_m");
    let depth_m = resolved_m(&parameters, "pocket_depth_m");
    let wall_m = resolved_m(&parameters, "wall_m");

    let carbons = (0..pocket.atom_count())
        .filter(|index| pocket.topology.element(*index) == Some(Element::CARBON))
        .count();
    let decorated = (0..pocket.atom_count())
        .filter(|index| {
            matches!(
                pocket.topology.element(*index),
                Some(Element::OXYGEN)
                    | Some(Element::NITROGEN)
                    | Some(Element::FLUORINE)
                    | Some(Element::CHLORINE)
                    | Some(Element::SULFUR)
            )
        })
        .count();

    print!("{{\"pocket_radius_m\":{radius_m}");
    print!(",\"pocket_depth_m\":{depth_m}");
    print!(",\"wall_m\":{wall_m}");
    print!(
        ",\"wall_group\":\"{}\"",
        BindingPocketGenerator.wall_group().label()
    );
    print!(",\"pocket_atoms\":{}", pocket.atom_count());
    print!(",\"pocket_carbons\":{carbons}");
    print!(",\"decorated_sites\":{decorated}");
    print!(
        ",\"cavity_volume_m3\":{}",
        PI * radius_m * radius_m * depth_m
    );
    print!(
        ",\"wall_contact_distance_m\":{}",
        wall_contact_distance_m().expect("carbon is in the table")
    );
    print!(",\"guests\":[");
    for (slot, molecule) in guests().iter().enumerate() {
        if slot > 0 {
            print!(",");
        }
        let clearance_m = radius_m - molecule.extent_m();
        print!(
            "{{\"id\":\"{}\",\"extent_m\":{},\"clearance_m\":{clearance_m}",
            molecule.id,
            molecule.extent_m()
        );
        print!(
            ",\"required_radius_m\":{}",
            pocket_radius_for(molecule, 0.0)
        );
        print!(
            ",\"wall_well_depth_j\":{}",
            wall_well_depth_j(molecule).expect("the guest elements are in the table")
        );
        print!(",\"fits\":{}}}", clearance_m > 0.0);
    }
    println!("]}}");
}

fn resolved_m(parameters: &ParameterSet, name: &str) -> f64 {
    BindingPocketGenerator
        .resolve(parameters)
        .expect("the parameters resolve")
        .require(name)
        .expect("the parameter is declared")
}
