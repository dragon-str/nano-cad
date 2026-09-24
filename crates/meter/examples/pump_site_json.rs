//! The binding search and the free volume of one gear-pump tooth gap (M27-04).
//!
//! The example builds the default pump, takes the iron site of gear A that
//! faces away from the mesh, and searches oxygen, water and carbon dioxide at
//! it with the metal-ligand term. The best pose energies feed the competition
//! model at alveolar partial pressures and body temperature. It also measures
//! the free volume of that tooth gap, from the root to the housing bore, and
//! states it as a count of liquid-water molecules.
//!
//! Every number is a rigid, solvent-free estimate. The coordination wells
//! come from the Fe(NH3)4 proxy, not from the diamond root site.

use std::f64::consts::PI;

use nanocad_meter::{
    competition_report, iron_coordination, metal_site_binding, Competitor, MetalSiteTarget,
    AVOGADRO_PER_MOL,
};
use nanocad_model::Element;
use nanocad_parts::{guest, GearPumpGenerator, ParameterSet};

const BODY_TEMPERATURE_K: f64 = 310.15;
const OXYGEN_PA: f64 = 13_300.0;
const WATER_VAPOUR_PA: f64 = 6_280.0;
const CARBON_DIOXIDE_PA: f64 = 5_330.0;
// Standard gas entropies at 298 K, one bar, J/(mol K). The guest loses all of
// it on binding (ADR-0086).
const OXYGEN_S: f64 = 205.2;
const WATER_S: f64 = 188.8;
const CARBON_DIOXIDE_S: f64 = 213.8;
/// The volume of one molecule of liquid water at 310 K, in cubic metres
/// (18.015 g/mol over 0.9934 g/cm^3).
const WATER_MOLECULE_M3: f64 = 3.011e-29;
/// The grid step of the free-volume count, in metres.
const GRID_M: f64 = 0.2e-10;
/// The gap of gear A whose root iron is searched. Gap 3 of eight faces -x,
/// away from the mesh.
const GAP: usize = 3;

/// The Bondi van der Waals radius of an element, in metres. Iron has no Bondi
/// value; 2.0e-10 m is a common choice for a transition metal.
fn vdw_radius_m(element: Element) -> f64 {
    match element.atomic_number() {
        1 => 1.20e-10,
        6 => 1.70e-10,
        7 => 1.55e-10,
        8 => 1.52e-10,
        _ => 2.0e-10,
    }
}

fn main() {
    let set = match GearPumpGenerator.build(&ParameterSet::new()) {
        Ok(set) => set,
        Err(error) => {
            eprintln!("the pump failed: {error}");
            std::process::exit(1);
        }
    };
    let design = &set.design;
    let iron = set.metal_sites[GAP];

    let target = MetalSiteTarget::default();
    let mut rows = Vec::new();
    let mut competitors = Vec::new();
    for (id, pressure_pa, entropy) in [
        ("o2", OXYGEN_PA, OXYGEN_S),
        ("water", WATER_VAPOUR_PA, WATER_S),
        ("co2", CARBON_DIOXIDE_PA, CARBON_DIOXIDE_S),
    ] {
        let Some(molecule) = guest(id) else {
            eprintln!("no guest {id}");
            std::process::exit(1);
        };
        let coordination = iron_coordination(id);
        let Some(report) =
            metal_site_binding(&set.part, iron, molecule, coordination.as_ref(), &target)
        else {
            eprintln!("the site search failed for {id}");
            std::process::exit(1);
        };
        let per_mol = |j: f64| j * AVOGADRO_PER_MOL / 1.0e3;
        let enthalpy_kj = if report.minimum_energy_j.is_finite() {
            per_mol(report.minimum_energy_j)
        } else {
            0.0
        };
        competitors.push(Competitor::from_molar(
            id,
            pressure_pa,
            enthalpy_kj,
            entropy,
        ));
        rows.push(format!(
            "{{\"guest\":\"{id}\",\"energy_kj_per_mol\":{:.2},\"coordination_kj_per_mol\":{:.2},\
             \"wall_kj_per_mol\":{:.2},\"donor_distance_m\":{:.3e},\"poses\":{},\"blocked\":{}}}",
            enthalpy_kj,
            per_mol(report.coordination_energy_j),
            per_mol(report.wall_energy_j),
            report.donor_distance_m,
            report.poses,
            report.blocked
        ));
    }
    let competition = competition_report(&competitors, BODY_TEMPERATURE_K);
    let occupancy: Vec<String> = competition
        .guests
        .iter()
        .map(|g| format!("\"{}\":{:.4e}", g.id, g.occupancy))
        .collect();

    // Free volume of the gap: grid points in the gap sector of gear A, from
    // the root circle to the housing bore, within the cavity height, that lie
    // outside every van der Waals sphere.
    let center_m = design.gear_a_center_m();
    let pitch_rad = 2.0 * PI / design.teeth() as f64;
    let gap_rad = (GAP as f64 + 0.5) * pitch_rad;
    let inner_m = design.root_radius_m() - 1.0e-10;
    let outer_m = design.cavity_radius_m();
    let half_height_m = design.cavity_half_height_m();
    let atoms: Vec<([f64; 3], f64)> = (0..set.part.atom_count())
        .filter_map(|i| {
            let p = set.part.topology.position_m(i)?;
            let e = set.part.topology.element(i)?;
            let r = (p[0] - center_m[0]).hypot(p[1] - center_m[1]);
            (r > inner_m - 4.0e-10).then_some((p, vdw_radius_m(e)))
        })
        .collect();
    let mut free = 0usize;
    let reach = (outer_m / GRID_M).ceil() as i64;
    let height = (half_height_m / GRID_M).ceil() as i64;
    for ix in -reach..=reach {
        for iy in -reach..=reach {
            let x = ix as f64 * GRID_M;
            let y = iy as f64 * GRID_M;
            let r = x.hypot(y);
            if r < inner_m || r > outer_m {
                continue;
            }
            let delta = (y.atan2(x) - gap_rad + PI).rem_euclid(2.0 * PI) - PI;
            if delta.abs() > 0.5 * pitch_rad {
                continue;
            }
            for iz in -height..=height {
                let p = [center_m[0] + x, center_m[1] + y, iz as f64 * GRID_M];
                if p[2].abs() > half_height_m {
                    continue;
                }
                let clear = atoms.iter().all(|(q, radius)| {
                    let d2 = (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2);
                    d2 >= radius * radius
                });
                if clear {
                    free += 1;
                }
            }
        }
    }
    let free_m3 = free as f64 * GRID_M * GRID_M * GRID_M;

    println!(
        "{{\"iron_atom\":{iron},\"gap\":{GAP},\"temperature_k\":{BODY_TEMPERATURE_K},\
         \"binding\":[{}],\"occupancy\":{{{},\"empty\":{:.4e}}},\
         \"gap_free_volume_m3\":{free_m3:.3e},\"gap_free_volume_waters\":{:.1},\
         \"note\":\"rigid, solvent-free estimates; wells from the Fe(NH3)4 proxy\"}}",
        rows.join(","),
        occupancy.join(","),
        competition.empty,
        free_m3 / WATER_MOLECULE_M3
    );
}
