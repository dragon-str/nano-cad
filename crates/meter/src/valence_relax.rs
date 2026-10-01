//! A valence force-field relaxation of one generated body.
//!
//! The generator builds a body on ideal lattice sites, then reconstructs the
//! surface (dimers, filled notches). The reconstruction moves a few atoms by
//! hand and leaves some bonds strained. This relaxation lets every atom move
//! under a harmonic valence force field, so the strain spreads into the
//! lattice as it does in a real crystal.
//!
//! The terms are harmonic bond stretch and angle bend on the topology bonds,
//! and a Lennard-Jones 12-6 term (UFF parameters from `nanocad_model`) for
//! every pair that is not 1-2, 1-3 or 1-4. The force constants are MM3-scale
//! values for alkanes (Allinger, Yuh and Lii, J. Am. Chem. Soc. 111, 8551,
//! 1989): C-C 4.49 mdyn/A, C-H 4.74 mdyn/A, C-C-C 0.67, C-C-H 0.59 and H-C-H
//! 0.55 mdyn A/rad^2. The rest angle is tetrahedral. The rest length of C-C is
//! the diamond bond, 1.544e-10 m. There is no torsion, no cross term and no
//! solvent, so the result is an estimate of the relaxed shape, not a
//! quantum-chemistry geometry.

use std::collections::HashSet;
use std::ops::Range;

use nanocad_engine::{
    minimize_with, AngleBendTerm, BondStretchTerm, Cutoff, EngineError, LennardJonesTerm,
    MinimizeMethod, MinimizeOptions, PeriodicBox, System,
};
use nanocad_model::{pair_params, Element, Topology};

/// C-C bond stiffness, in newtons per metre (MM3 4.49 mdyn/A).
const C_C_K_N_PER_M: f64 = 449.0;
/// C-H bond stiffness, in newtons per metre (MM3 4.74 mdyn/A).
const C_H_K_N_PER_M: f64 = 474.0;
/// The diamond C-C bond, in metres.
const C_C_R0_M: f64 = 1.544e-10;
/// The C-H bond, in metres.
const C_H_R0_M: f64 = 1.09e-10;
/// Fe-C bond stiffness, in newtons per metre. It is the iron-ligand value of
/// the metal-ligand example (ADR-0083), an estimate.
const FE_C_K_N_PER_M: f64 = 200.0;
/// The Fe-C sigma bond of a root site, in metres.
const FE_C_R0_M: f64 = 1.95e-10;
/// The stiffness of a bond with a rest length taken from the built
/// geometry, in newtons per metre (an aromatic-bond scale).
const FOREIGN_K_N_PER_M: f64 = 500.0;

fn is_ch(element: Element) -> bool {
    element == Element::CARBON || element == Element::HYDROGEN
}

fn is_fe_c(a: Element, b: Element) -> bool {
    (a == Element::IRON && b == Element::CARBON) || (a == Element::CARBON && b == Element::IRON)
}

fn distance(flat: &[f64], u: usize, v: usize) -> f64 {
    ((flat[3 * u] - flat[3 * v]).powi(2)
        + (flat[3 * u + 1] - flat[3 * v + 1]).powi(2)
        + (flat[3 * u + 2] - flat[3 * v + 2]).powi(2))
    .sqrt()
}

fn angle(flat: &[f64], a: usize, b: usize, c: usize) -> f64 {
    let d1 = [
        flat[3 * a] - flat[3 * b],
        flat[3 * a + 1] - flat[3 * b + 1],
        flat[3 * a + 2] - flat[3 * b + 2],
    ];
    let d2 = [
        flat[3 * c] - flat[3 * b],
        flat[3 * c + 1] - flat[3 * b + 1],
        flat[3 * c + 2] - flat[3 * b + 2],
    ];
    let dot = d1[0] * d2[0] + d1[1] * d2[1] + d1[2] * d2[2];
    let n = (d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2]).sqrt()
        * (d2[0] * d2[0] + d2[1] * d2[1] + d2[2] * d2[2]).sqrt();
    (dot / n).clamp(-1.0, 1.0).acos()
}

/// C-C-C bend constant, in joules per radian squared (0.67 mdyn A/rad^2).
const CCC_K_J_PER_RAD2: f64 = 0.67e-18;
/// C-C-H bend constant, in joules per radian squared.
const CCH_K_J_PER_RAD2: f64 = 0.59e-18;
/// H-C-H bend constant, in joules per radian squared.
const HCH_K_J_PER_RAD2: f64 = 0.55e-18;
/// The mass the minimizer needs. The minimum does not depend on it.
const ATOM_MASS_KG: f64 = 2.0e-26;

/// The settings of a valence relaxation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ValenceRelaxOptions {
    /// The Lennard-Jones cutoff, in metres.
    pub cutoff_m: f64,
    /// The iteration budget of the minimizer.
    pub max_iterations: usize,
    /// The gradient tolerance, in newtons.
    pub gradient_tolerance_n: f64,
    /// The first line-search step, in metres.
    pub initial_step_m: f64,
    /// The first `frozen_atoms` atoms of the range stay fixed: they act on
    /// the others (bonds, angles, contacts) but never move. A fixed wall.
    pub frozen_atoms: usize,
}

impl Default for ValenceRelaxOptions {
    fn default() -> Self {
        Self {
            cutoff_m: 5.0e-10,
            max_iterations: 300,
            gradient_tolerance_n: 1.0e-11,
            initial_step_m: 1.0e-12,
            frozen_atoms: 0,
        }
    }
}

/// The result of a valence relaxation.
#[derive(Clone, Debug, PartialEq)]
pub struct ValenceRelaxReport {
    /// The number of atoms relaxed.
    pub atom_count: usize,
    /// The number of bonds.
    pub bond_count: usize,
    /// The number of angles.
    pub angle_count: usize,
    /// The energy before the relaxation, in joules.
    pub initial_energy_j: f64,
    /// The energy after the relaxation, in joules.
    pub final_energy_j: f64,
    /// The largest remaining force component, in newtons.
    pub gradient_norm_n: f64,
    /// The minimizer iterations.
    pub iterations: usize,
    /// Whether the gradient reached the tolerance.
    pub converged: bool,
    /// The longest C-C bond before the relaxation, in metres.
    pub initial_max_c_c_m: f64,
    /// The longest C-C bond after the relaxation, in metres.
    pub max_c_c_m: f64,
    /// The largest atom displacement, in metres.
    pub max_displacement_m: f64,
    /// The largest carbon displacement, in metres. A hydrogen can swing
    /// freely about its bond, because the field has no torsion term, so the
    /// carbon value measures the framework.
    pub max_carbon_displacement_m: f64,
    /// The root-mean-square atom displacement, in metres.
    pub rms_displacement_m: f64,
    /// The shortest distance between two relaxed atoms that are not 1-2 or
    /// 1-3, in metres. Infinity when no pair is within the cutoff.
    pub min_nonbonded_m: f64,
    /// The relaxed positions, in the order of the atom range.
    pub positions_m: Vec<[f64; 3]>,
}

/// Relaxes the atoms `range` of a topology and returns the relaxed positions.
///
/// Only bonds with both ends in the range count. An element without UFF data
/// is an error.
pub fn relax_valence(
    topology: &Topology,
    range: Range<usize>,
    options: &ValenceRelaxOptions,
) -> Result<ValenceRelaxReport, EngineError> {
    let start = range.start;
    let count = range.len();
    let local = |index: u32| -> Option<usize> {
        let index = index as usize;
        range.contains(&index).then(|| index - start)
    };
    let mut elements = Vec::with_capacity(count);
    let mut flat = Vec::with_capacity(3 * count);
    for index in range.clone() {
        elements.push(topology.element(index).unwrap_or(Element::CARBON));
        let position_m = topology.position_m(index).unwrap_or([0.0; 3]);
        flat.extend_from_slice(&position_m);
    }

    let mut neighbours: Vec<Vec<usize>> = vec![Vec::new(); count];
    // A carbon with three bonds is sp2 (a porphine carbon): its bonds and
    // angles keep their built geometry. Every diamond carbon has four.
    let mut degree = vec![0usize; count];
    for bond in topology.bonds() {
        if let (Some(u), Some(v)) = (local(bond.u), local(bond.v)) {
            degree[u] += 1;
            degree[v] += 1;
        }
    }
    let sp2: Vec<bool> = (0..count)
        .map(|a| elements[a] == Element::CARBON && degree[a] == 3)
        .collect();
    let mut bonds = BondStretchTerm::new();
    let mut c_c_pairs = Vec::new();
    for bond in topology.bonds() {
        let (Some(u), Some(v)) = (local(bond.u), local(bond.v)) else {
            continue;
        };
        let hydrogen = elements[u] == Element::HYDROGEN || elements[v] == Element::HYDROGEN;
        let iron = elements[u] == Element::IRON || elements[v] == Element::IRON;
        // A bond to an atom other than C or H (the porphine N, a bound O2, an
        // amine tip) keeps its built length: those geometries come from the
        // quantum-chemistry models, and this field has no parameters for them.
        let foreign = [u, v].iter().any(|&a| !is_ch(elements[a]) || sp2[a]);
        let (k_n_per_m, r0_m) = if foreign && !(iron && is_fe_c(elements[u], elements[v])) {
            (FOREIGN_K_N_PER_M, distance(&flat, u, v))
        } else if iron {
            (FE_C_K_N_PER_M, FE_C_R0_M)
        } else if hydrogen {
            (C_H_K_N_PER_M, C_H_R0_M)
        } else {
            c_c_pairs.push((u, v));
            (C_C_K_N_PER_M, C_C_R0_M)
        };
        bonds.add_bond(u as u32, v as u32, k_n_per_m, r0_m)?;
        neighbours[u].push(v);
        neighbours[v].push(u);
    }
    let bond_count = bonds.bond_count();

    let tetrahedral_rad = (-1.0_f64 / 3.0).acos();
    let mut angles = AngleBendTerm::new();
    let mut excluded: HashSet<(usize, usize)> = HashSet::new();
    for (centre, list) in neighbours.iter().enumerate() {
        for &other in list {
            excluded.insert((centre.min(other), centre.max(other)));
        }
        for (i, &left) in list.iter().enumerate() {
            for &right in &list[i + 1..] {
                let hydrogens = [left, right]
                    .iter()
                    .filter(|&&atom| elements[atom] == Element::HYDROGEN)
                    .count();
                let k = match hydrogens {
                    0 => CCC_K_J_PER_RAD2,
                    1 => CCH_K_J_PER_RAD2,
                    _ => HCH_K_J_PER_RAD2,
                };
                let rest_rad = if [left, centre, right]
                    .iter()
                    .any(|&a| !is_ch(elements[a]) || sp2[a])
                {
                    angle(&flat, left, centre, right)
                } else {
                    tetrahedral_rad
                };
                // A straight angle (trans N-Fe-N in a porphine) has no bend
                // direction; the other angles hold that geometry.
                if rest_rad > 172.0_f64.to_radians() {
                    excluded.insert((left.min(right), left.max(right)));
                    continue;
                }
                angles.add_angle(left as u32, centre as u32, right as u32, k, rest_rad)?;
                excluded.insert((left.min(right), left.max(right)));
            }
        }
    }
    let angle_count = angles.angle_count();
    let close_pairs = excluded.clone();
    // The force field has no torsion term, so the 1-4 pairs carry no
    // nonbonded term either. UFF Lennard-Jones between 1-4 carbons, 2.96e-10 m
    // apart in diamond with sigma 3.43e-10 m, would push the lattice apart.
    for (b, list) in neighbours.iter().enumerate() {
        for &c in list {
            if c < b {
                continue;
            }
            for &a in &neighbours[b] {
                if a == c {
                    continue;
                }
                for &d in &neighbours[c] {
                    if d != b && d != a {
                        excluded.insert((a.min(d), a.max(d)));
                    }
                }
            }
        }
    }

    let cutoff = Cutoff::new(options.cutoff_m, 0.9 * options.cutoff_m)?;
    let mut lennard_jones = LennardJonesTerm::new(cutoff, PeriodicBox::non_periodic())?;
    for &element in &elements {
        let (epsilon_j, sigma_m) =
            pair_params(element, element).ok_or(EngineError::InvalidVdwParameters { atom: 0 })?;
        lennard_jones.add_atom(epsilon_j, sigma_m)?;
    }

    let mut system = System::with_masses_kg(count, &vec![ATOM_MASS_KG; count])?;
    system.set_bond_stretch(bonds)?;
    system.set_angle_bend(angles)?;
    system.set_lennard_jones(lennard_jones)?;
    if options.frozen_atoms > 0 {
        system.set_frozen_atoms((0..count).map(|a| a < options.frozen_atoms).collect())?;
    }
    // Inside one diamond lattice the bond rest length is already the crystal
    // value, so a Lennard-Jones term between its carbons double counts and
    // contracts the lattice by about 4 percent (ADR-0111). Exclude pairs of
    // four-bonded carbons of the same connected piece within the cutoff.
    {
        let mut parent: Vec<usize> = (0..count).collect();
        fn find(parent: &mut [usize], mut x: usize) -> usize {
            while parent[x] != x {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            x
        }
        for (a, list) in neighbours.iter().enumerate() {
            for &b in list {
                let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
                if ra != rb {
                    parent[ra] = rb;
                }
            }
        }
        let lattice: Vec<usize> = (0..count)
            .filter(|&a| elements[a] == Element::CARBON && degree[a] == 4)
            .collect();
        let cut = options.cutoff_m;
        let key = |a: usize| {
            (
                (flat[3 * a] / cut).floor() as i64,
                (flat[3 * a + 1] / cut).floor() as i64,
                (flat[3 * a + 2] / cut).floor() as i64,
            )
        };
        let mut grid: std::collections::HashMap<(i64, i64, i64), Vec<usize>> =
            std::collections::HashMap::new();
        for &a in &lattice {
            grid.entry(key(a)).or_default().push(a);
        }
        for &a in &lattice {
            let (x, y, z) = key(a);
            let root_a = find(&mut parent, a);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        for &b in grid.get(&(x + dx, y + dy, z + dz)).into_iter().flatten() {
                            if b <= a || distance(&flat, a, b) > cut {
                                continue;
                            }
                            if find(&mut parent, b) == root_a {
                                excluded.insert((a, b));
                            }
                        }
                    }
                }
            }
        }
    }
    let mut pairs: Vec<_> = excluded.into_iter().collect();
    pairs.sort_unstable();
    for (i, j) in pairs {
        system.add_exclusion(i as u32, j as u32)?;
    }
    system.rebuild_neighbors(&flat)?;
    let initial_energy_j = system.energy_j(&flat)?;

    let max_c_c = |positions: &[f64]| {
        c_c_pairs.iter().fold(0.0_f64, |max, &(u, v)| {
            let d = ((positions[3 * u] - positions[3 * v]).powi(2)
                + (positions[3 * u + 1] - positions[3 * v + 1]).powi(2)
                + (positions[3 * u + 2] - positions[3 * v + 2]).powi(2))
            .sqrt();
            max.max(d)
        })
    };
    let initial_max_c_c_m = max_c_c(&flat);

    let mut relaxed = flat.clone();
    let result = minimize_with(
        &mut system,
        &mut relaxed,
        &MinimizeOptions {
            max_iterations: options.max_iterations,
            gradient_tolerance_n: options.gradient_tolerance_n,
            initial_step_m: options.initial_step_m,
        },
        MinimizeMethod::ConjugateGradientThenLbfgs,
    )?;

    let mut max_displacement_m = 0.0_f64;
    let mut max_carbon_displacement_m = 0.0_f64;
    let mut sum_squared_m2 = 0.0_f64;
    let mut positions_m = Vec::with_capacity(count);
    for atom in 0..count {
        let p = [
            relaxed[3 * atom],
            relaxed[3 * atom + 1],
            relaxed[3 * atom + 2],
        ];
        let d = ((p[0] - flat[3 * atom]).powi(2)
            + (p[1] - flat[3 * atom + 1]).powi(2)
            + (p[2] - flat[3 * atom + 2]).powi(2))
        .sqrt();
        max_displacement_m = max_displacement_m.max(d);
        sum_squared_m2 += d * d;
        if elements[atom] == Element::CARBON {
            max_carbon_displacement_m = max_carbon_displacement_m.max(d);
        }
        positions_m.push(p);
    }

    let min_nonbonded_m = min_nonbonded_m(&positions_m, &close_pairs, options.cutoff_m);

    Ok(ValenceRelaxReport {
        atom_count: count,
        bond_count,
        angle_count,
        initial_energy_j,
        final_energy_j: result.energy_j,
        gradient_norm_n: result.gradient_norm_n,
        iterations: result.iterations,
        converged: result.converged,
        initial_max_c_c_m,
        max_c_c_m: max_c_c(&relaxed),
        max_displacement_m,
        max_carbon_displacement_m,
        rms_displacement_m: (sum_squared_m2 / count.max(1) as f64).sqrt(),
        min_nonbonded_m,
        positions_m,
    })
}

/// The shortest distance between two atoms that are not in `skip`, found on
/// a grid of cells one cutoff wide.
fn min_nonbonded_m(positions_m: &[[f64; 3]], skip: &HashSet<(usize, usize)>, cutoff_m: f64) -> f64 {
    let cell = |p: [f64; 3]| {
        (
            (p[0] / cutoff_m).floor() as i64,
            (p[1] / cutoff_m).floor() as i64,
            (p[2] / cutoff_m).floor() as i64,
        )
    };
    let mut grid: std::collections::HashMap<(i64, i64, i64), Vec<usize>> =
        std::collections::HashMap::new();
    for (index, p) in positions_m.iter().enumerate() {
        grid.entry(cell(*p)).or_default().push(index);
    }
    let mut best = f64::INFINITY;
    for (i, p) in positions_m.iter().enumerate() {
        let (cx, cy, cz) = cell(*p);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    for &j in grid.get(&(cx + dx, cy + dy, cz + dz)).into_iter().flatten() {
                        if j <= i || skip.contains(&(i, j)) {
                            continue;
                        }
                        let q = positions_m[j];
                        let d =
                            ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2))
                                .sqrt();
                        best = best.min(d);
                    }
                }
            }
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;
    use nanocad_model::{Atom, Bond, BondType};

    /// A stretched ethane-like C-C pair with capped ends relaxes back toward
    /// the diamond bond length.
    #[test]
    fn a_stretched_bond_relaxes() {
        let mut topology = Topology::new();
        let a = topology.add_atom(Atom::new(Element::CARBON, [0.0, 0.0, 0.0], 0.0, "C"));
        let b = topology.add_atom(Atom::new(Element::CARBON, [2.0e-10, 0.0, 0.0], 0.0, "C"));
        topology
            .add_bond(Bond::new(a, b, 1, BondType::Single))
            .expect("bond");
        let report =
            relax_valence(&topology, 0..2, &ValenceRelaxOptions::default()).expect("relax");
        assert!((report.initial_max_c_c_m - 2.0e-10).abs() < 1.0e-15);
        assert!(
            (report.max_c_c_m - C_C_R0_M).abs() < 1.0e-12,
            "{}",
            report.max_c_c_m
        );
        assert!(report.final_energy_j < report.initial_energy_j);
    }

    /// A frozen atom holds its place; the free end of the bond relaxes onto
    /// it alone.
    #[test]
    fn a_frozen_atom_does_not_move() {
        let mut topology = Topology::new();
        let a = topology.add_atom(Atom::new(Element::CARBON, [0.0, 0.0, 0.0], 0.0, "C"));
        let b = topology.add_atom(Atom::new(Element::CARBON, [2.0e-10, 0.0, 0.0], 0.0, "C"));
        topology
            .add_bond(Bond::new(a, b, 1, BondType::Single))
            .expect("bond");
        let options = ValenceRelaxOptions {
            frozen_atoms: 1,
            ..ValenceRelaxOptions::default()
        };
        let report = relax_valence(&topology, 0..2, &options).expect("relax");
        assert_eq!(report.positions_m[0], [0.0, 0.0, 0.0]);
        assert!((report.positions_m[1][0] - C_C_R0_M).abs() < 1.0e-12);
    }
}
