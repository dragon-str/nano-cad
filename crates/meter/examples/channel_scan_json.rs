//! Drain-channel kill test, part 1 (M28-09): a rigid scan of water and oxygen
//! through a narrow lined channel in diamond.
//!
//! The channel is a `binding_pocket` whose well runs through the whole block,
//! so it is open at both ends. Its wall carries one functional group kind.
//! Each guest moves along the channel axis; at each step the scan tries
//! every orientation (a sphere of directions times spins) and small offsets
//! from the axis, and keeps the lowest energy (UFF Lennard-Jones plus point
//! charges, no solvent). The barrier is the highest of those minima along the
//! axis, relative to the guest far outside. A blocked step (every pose
//! overlaps) is reported as a barrier of infinity.
//!
//! The rigid scan measures steric exclusion. It does not include the loss of
//! hydrogen bonds for water leaving bulk liquid; the MD test does that.
//!
//! Arguments: `<wall_group index> <radius_m>...`

use nanocad_model::{pair_params, Element};
use nanocad_parts::{guest, BindingPocketGenerator, Guest, ParameterSet, PartGenerator};

const COULOMB: f64 = 8.987_551_792_3e9;
const AVOGADRO: f64 = 6.022_140_76e23;
const THICKNESS_M: f64 = 1.0e-9;
const OVERLAP_RATIO: f64 = 0.85;

fn sphere(count: usize) -> Vec<[f64; 3]> {
    let golden = std::f64::consts::PI * (3.0 - 5.0_f64.sqrt());
    (0..count)
        .map(|i| {
            let z = 1.0 - 2.0 * (i as f64 + 0.5) / count as f64;
            let r = (1.0 - z * z).sqrt();
            [r * (golden * i as f64).cos(), r * (golden * i as f64).sin(), z]
        })
        .collect()
}

/// Rotation matrix taking the z axis to `d`, then spun by `spin` about `d`.
fn frame(d: [f64; 3], spin: f64) -> [[f64; 3]; 3] {
    let helper = if d[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
    let mut u = [
        helper[1] * d[2] - helper[2] * d[1],
        helper[2] * d[0] - helper[0] * d[2],
        helper[0] * d[1] - helper[1] * d[0],
    ];
    let n = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
    u = [u[0] / n, u[1] / n, u[2] / n];
    let v = [
        d[1] * u[2] - d[2] * u[1],
        d[2] * u[0] - d[0] * u[2],
        d[0] * u[1] - d[1] * u[0],
    ];
    let (s, c) = spin.sin_cos();
    let a = [c * u[0] + s * v[0], c * u[1] + s * v[1], c * u[2] + s * v[2]];
    let b = [-s * u[0] + c * v[0], -s * u[1] + c * v[1], -s * u[2] + c * v[2]];
    // Columns: a, b, d. Body x -> a, body y -> b, body z -> d.
    [[a[0], b[0], d[0]], [a[1], b[1], d[1]], [a[2], b[2], d[2]]]
}

fn energy_at(
    wall: &[([f64; 3], f64, Element)],
    guest: &Guest,
    centre: [f64; 3],
    rotation: &[[f64; 3]; 3],
) -> f64 {
    let c = guest.centroid_m();
    let mut e = 0.0;
    for atom in guest.atoms {
        let l = [
            atom.position_m[0] - c[0],
            atom.position_m[1] - c[1],
            atom.position_m[2] - c[2],
        ];
        let p = [
            centre[0] + rotation[0][0] * l[0] + rotation[0][1] * l[1] + rotation[0][2] * l[2],
            centre[1] + rotation[1][0] * l[0] + rotation[1][1] * l[1] + rotation[1][2] * l[2],
            centre[2] + rotation[2][0] * l[0] + rotation[2][1] * l[1] + rotation[2][2] * l[2],
        ];
        for &(w, q, element) in wall {
            let r2 = (w[0] - p[0]).powi(2) + (w[1] - p[1]).powi(2) + (w[2] - p[2]).powi(2);
            if r2 > 1.0e-18 {
                continue;
            }
            let Some((eps, sigma)) = pair_params(atom.element, element) else {
                continue;
            };
            let floor = OVERLAP_RATIO * sigma;
            if r2 < floor * floor {
                return f64::INFINITY;
            }
            let x = sigma * sigma / r2;
            let x6 = x * x * x;
            e += 4.0 * eps * (x6 * x6 - x6) + COULOMB * atom.charge_c * q / r2.sqrt();
        }
    }
    e
}

fn main() {
    let mut args = std::env::args().skip(1);
    let group: f64 = args.next().and_then(|v| v.parse().ok()).unwrap_or(0.0);
    let radii: Vec<f64> = args.filter_map(|v| v.parse().ok()).collect();
    let directions = sphere(40);
    let spins = 8;
    let mut rows = Vec::new();
    for radius_m in radii {
        let part = match BindingPocketGenerator.generate(
            &ParameterSet::new()
                .with("pocket_radius_m", radius_m)
                .with("pocket_depth_m", THICKNESS_M)
                .with("thickness_m", THICKNESS_M)
                .with("wall_m", 0.8e-9)
                .with("wall_group", group),
        ) {
            Ok(part) => part,
            Err(error) => {
                rows.push(format!("{{\"radius_m\":{radius_m:.3e},\"error\":\"{error}\"}}"));
                continue;
            }
        };
        let wall: Vec<([f64; 3], f64, Element)> = (0..part.atom_count())
            .filter_map(|i| {
                Some((
                    part.topology.position_m(i)?,
                    part.topology.charge_c(i).unwrap_or(0.0),
                    part.topology.element(i)?,
                ))
            })
            .collect();
        // The narrowest free radius: distance from the axis to the nearest
        // atom inside the block height, minus that atom's UFF radius.
        let free_m = wall
            .iter()
            .filter(|(p, _, _)| p[2].abs() < 0.4 * THICKNESS_M)
            .map(|(p, _, e)| {
                let r = if *e == Element::HYDROGEN { 1.2e-10 } else { 1.52e-10 };
                p[0].hypot(p[1]) - r
            })
            .fold(f64::INFINITY, f64::min);
        let mut cells = Vec::new();
        for id in ["water", "o2"] {
            let g = guest(id).expect("guest");
            let mut barrier = f64::NEG_INFINITY;
            let mut steps = 0;
            let mut z = -0.9 * THICKNESS_M;
            while z <= 0.9 * THICKNESS_M + 1.0e-15 {
                let mut best = f64::INFINITY;
                for ox in [-0.4e-10, 0.0, 0.4e-10] {
                    for oy in [-0.4e-10, 0.0, 0.4e-10] {
                        for d in &directions {
                            for s in 0..spins {
                                let rot = frame(*d, std::f64::consts::TAU * s as f64 / spins as f64);
                                best = best.min(energy_at(&wall, g, [ox, oy, z], &rot));
                            }
                        }
                    }
                }
                barrier = barrier.max(best);
                steps += 1;
                z += 0.25e-10;
            }
            let kj = barrier * AVOGADRO / 1.0e3;
            cells.push(format!(
                "\"{id}_barrier_kj_per_mol\":{}",
                if kj.is_finite() { format!("{kj:.1}") } else { "null".to_string() }
            ));
            let _ = steps;
        }
        rows.push(format!(
            "{{\"radius_m\":{radius_m:.3e},\"free_radius_m\":{free_m:.3e},{}}}",
            cells.join(",")
        ));
    }
    println!("[{}]", rows.join(","));
}
