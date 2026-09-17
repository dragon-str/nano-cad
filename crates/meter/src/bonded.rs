//! The bonded terms beyond bond stretch: torsion and out-of-plane.
//!
//! A diamond lattice is staggered, so the torsion potential is already at its
//! minimum at the generated geometry. The term adds stiffness without a net
//! force. An sp3 carbon has four neighbours, so its geometry is not planar and
//! an out-of-plane term would fight the lattice. The out-of-plane term is
//! therefore added only to a carbon with exactly three neighbours, where the
//! pucker is real. A perfect lattice has no such carbon.

use nanocad_engine::{EngineError, OutOfPlaneTerm, TorsionTerm};

/// The threefold torsion barrier of a carbon-carbon single bond, in joules.
///
/// The alkane barrier is about 12 kJ per mole.
pub const TORSION_V3_J: f64 = 2.0e-20;

/// The out-of-plane force constant of a trigonal carbon, in joules per radian
/// squared.
pub const IMPROPER_K_J_PER_RAD2: f64 = 2.0e-18;

/// Adds one torsion for every four-atom path inside a cluster.
///
/// The neighbour lists hold local indices. Each dihedral `i-j-k-l` is added
/// once, with `i < k`, so the count is the number of unique paths. The return
/// value is the number of torsions added.
pub fn add_torsions(term: &mut TorsionTerm, neighbours: &[Vec<u32>]) -> Result<usize, EngineError> {
    let mut added = 0;
    for j in 0..neighbours.len() {
        for &i in &neighbours[j] {
            for &k in &neighbours[j] {
                if k <= i {
                    continue;
                }
                for &l in &neighbours[k as usize] {
                    if l == j as u32 {
                        continue;
                    }
                    term.add_torsion(i, j as u32, k, l, 0.0, 0.0, TORSION_V3_J)?;
                    added += 1;
                }
            }
        }
    }
    Ok(added)
}

/// The out-of-plane angle below which a three-coordinate atom counts as
/// trigonal, in radians.
///
/// A tetrahedral centre sits at about 0.616 rad, so this threshold keeps it
/// out. An sp2 defect centre sits near zero and is kept.
pub const TRIGONAL_TOLERANCE_RAD: f64 = 0.5;

/// Adds one out-of-plane term for every TRIGONAL three-coordinate atom.
///
/// A cut-out cluster has boundary carbons with three neighbours inside the
/// cluster, but their geometry stays tetrahedral. The term is added only when
/// the measured out-of-plane angle is already close to planar, so a pyramidal
/// sp3 centre does not fight the lattice. The return value is the number of
/// impropers added.
pub fn add_impropers(
    term: &mut OutOfPlaneTerm,
    neighbours: &[Vec<u32>],
    positions_m: &[[f64; 3]],
) -> Result<usize, EngineError> {
    let mut added = 0;
    for (central, list) in neighbours.iter().enumerate() {
        if list.len() != 3 {
            continue;
        }
        let (a, b, c) = (list[0], list[1], list[2]);
        let chi = out_of_plane_rad(positions_m, a, central as u32, b, c);
        if chi.abs() >= TRIGONAL_TOLERANCE_RAD {
            continue;
        }
        term.add_improper(a, central as u32, c, b, IMPROPER_K_J_PER_RAD2, 0.0)?;
        added += 1;
    }
    Ok(added)
}

/// Returns the out-of-plane angle of atom `out` over the plane `i, centre, k`.
fn out_of_plane_rad(positions_m: &[[f64; 3]], i: u32, centre: u32, out: u32, k: u32) -> f64 {
    let origin = positions_m[centre as usize];
    let first = subtract(positions_m[i as usize], origin);
    let second = subtract(positions_m[k as usize], origin);
    let normal = cross(first, second);
    let length = norm(normal);
    if length == 0.0 {
        return 0.0;
    }
    let unit = scale(normal, 1.0 / length);
    let arm = normalize(subtract(positions_m[out as usize], origin));
    dot(arm, unit).clamp(-1.0, 1.0).asin()
}

fn subtract(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn scale(a: [f64; 3], factor: f64) -> [f64; 3] {
    [a[0] * factor, a[1] * factor, a[2] * factor]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn normalize(a: [f64; 3]) -> [f64; 3] {
    let length = norm(a);
    if length == 0.0 {
        a
    } else {
        scale(a, 1.0 / length)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A chain of four atoms gives exactly one torsion.
    #[test]
    fn a_four_atom_chain_has_one_torsion() {
        let neighbours = vec![vec![1], vec![0, 2], vec![1, 3], vec![2]];
        let mut term = TorsionTerm::new();
        assert_eq!(add_torsions(&mut term, &neighbours).unwrap(), 1);
    }

    /// A tetrahedral centre with four arms gives six torsions.
    #[test]
    fn a_tetrahedral_centre_has_six_torsions() {
        let neighbours = vec![
            vec![1, 2, 3, 4],
            vec![0, 5],
            vec![0, 6],
            vec![0, 7],
            vec![0, 8],
            vec![1],
            vec![2],
            vec![3],
            vec![4],
        ];
        let mut term = TorsionTerm::new();
        assert_eq!(add_torsions(&mut term, &neighbours).unwrap(), 6);
    }

    /// A planar three-coordinate centre gives one improper.
    #[test]
    fn a_trigonal_centre_has_one_improper() {
        let neighbours = vec![vec![1, 2, 3], vec![0], vec![0], vec![0]];
        let positions_m = [
            [0.0, 0.0, 0.0],
            [1.0e-10, 0.0, 0.0],
            [0.0, 1.0e-10, 0.0],
            [1.0e-10, 1.0e-10, 0.0],
        ];
        let mut term = OutOfPlaneTerm::new();
        assert_eq!(
            add_impropers(&mut term, &neighbours, &positions_m).unwrap(),
            1
        );
    }

    /// A four-coordinate centre gives no improper.
    #[test]
    fn a_tetrahedral_centre_has_no_improper() {
        let neighbours = vec![vec![1, 2, 3, 4], vec![0], vec![0], vec![0], vec![0]];
        let positions_m = [[0.0; 3]; 5];
        let mut term = OutOfPlaneTerm::new();
        assert_eq!(
            add_impropers(&mut term, &neighbours, &positions_m).unwrap(),
            0
        );
    }

    /// A pyramidal three-coordinate centre is kept out of the improper term.
    #[test]
    fn a_pyramidal_centre_has_no_improper() {
        let neighbours = vec![vec![1, 2, 3], vec![0], vec![0], vec![0]];
        let positions_m = [
            [0.0, 0.0, 0.0],
            [1.0e-10, 0.0, 0.0],
            [0.0, 1.0e-10, 0.0],
            [0.0, 0.0, 3.0e-10],
        ];
        let mut term = OutOfPlaneTerm::new();
        assert_eq!(
            add_impropers(&mut term, &neighbours, &positions_m).unwrap(),
            0
        );
    }
}

#[cfg(test)]
mod relaxation {
    use super::TORSION_V3_J;
    use nanocad_engine::{
        hessian_finite_difference, mass_weighted, minimize_with, spectrum, symmetric_eigenvalues,
        BondStretchTerm, MinimizeMethod, MinimizeOptions, System, TorsionTerm,
    };

    /// A relaxed torsion has no negative curvature.
    ///
    /// The hessian is a finite difference of the gradient, so it is only
    /// positive at a stationary point. The fixture is relaxed first, then the
    /// spectrum is checked. This proves the analytic torsion gradient agrees
    /// with the torsion energy.
    #[test]
    fn a_relaxed_torsion_has_no_negative_curvature() {
        const CARBON_MASS_KG: f64 = 1.992_646_879_92e-26;
        let a = 1.544e-10;
        let positions = [
            [0.0, 0.0, 0.0],
            [a, 0.0, 0.0],
            [a, a, 0.0],
            [a * 0.5, a, a * 0.866],
        ];
        let n = 4;
        let flat: Vec<f64> = positions.iter().flat_map(|p| p.iter().copied()).collect();
        let masses = vec![CARBON_MASS_KG; n];
        let mut system = System::with_masses_kg(n, &masses).unwrap();
        let mut bonds = BondStretchTerm::new();
        for (u, v) in [(0usize, 1usize), (1, 2), (1, 3)] {
            let d = ((positions[u][0] - positions[v][0]).powi(2)
                + (positions[u][1] - positions[v][1]).powi(2)
                + (positions[u][2] - positions[v][2]).powi(2))
            .sqrt();
            bonds.add_bond(u as u32, v as u32, 300.0, d).unwrap();
        }
        system.set_bond_stretch(bonds).unwrap();
        let mut torsions = TorsionTerm::new();
        torsions
            .add_torsion(2, 1, 0, 3, 0.0, 0.0, TORSION_V3_J)
            .unwrap();
        system.set_torsion(torsions).unwrap();

        let mut relaxed = flat.clone();
        let options = MinimizeOptions {
            max_iterations: 5000,
            gradient_tolerance_n: 1.0e-16,
            initial_step_m: 1.0e-14,
        };
        let outcome = minimize_with(
            &mut system,
            &mut relaxed,
            &options,
            MinimizeMethod::ConjugateGradientThenLbfgs,
        )
        .unwrap();
        assert!(outcome.converged, "the fixture did not relax");

        let hessian = hessian_finite_difference(&mut system, &relaxed, 1.0e-13).unwrap();
        let weighted = mass_weighted(&hessian, &masses).unwrap();
        let eigenvalues = symmetric_eigenvalues(&weighted, 3 * n).unwrap();
        let found = spectrum(eigenvalues);
        assert_eq!(
            found.unstable_count, 0,
            "the relaxed torsion has negative curvature: {:?}",
            found.eigenvalues
        );
    }
}
