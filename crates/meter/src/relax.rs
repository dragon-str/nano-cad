//! The bonded-lattice relaxation of a representative sub-assembly.
//!
//! The metric takes a set of atoms and relaxes the covalent network to the
//! nearest potential minimum. The bonds come from the geometry: two atoms of
//! the same body closer than the bond cutoff are one covalent bond, and the
//! rest length is the current distance.
//!
//! The metric checks that the covalent lattice is a stable minimum. It does
//! not apply the van der Waals term between bodies: a mechanism holds its
//! bodies by joints, not by chemistry, and the stiff Buckingham term against a
//! soft lattice has no stable equilibrium. The separate clash count reports
//! any non-bonded overlap in the geometry, so the metric still detects a bad
//! placement.
//!
//! The report answers the question "does this assembly fall apart?". A stable
//! assembly drops little energy, keeps every bond strain small, moves every
//! atom by a small distance, and has no clash. The metric does not model the
//! solvent and it does not include the torsion or the out-of-plane terms.

use nanocad_engine::{minimize_with, BondStretchTerm, MinimizeMethod, MinimizeOptions, System};
use nanocad_model::Element;

/// The mass of one atom in kilograms. The minimizer needs masses, but the
/// gradient does not, so the value only has to be positive.
const ATOM_MASS_KG: f64 = 1.992_646_879_92e-26;

/// The bond stiffness of the lattice, in newtons per metre.
const BOND_STIFFNESS_N_PER_M: f64 = 300.0;

/// The ideal carbon-carbon bond length, in metres.
const C_C_BOND_M: f64 = 1.544e-10;

/// The ideal carbon-hydrogen bond length, in metres.
const C_H_BOND_M: f64 = 1.09e-10;

/// The ideal hydrogen-hydrogen separation when it is a bond, in metres.
const H_H_BOND_M: f64 = 7.4e-11;

/// One atom of the sub-assembly.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RelaxAtom {
    /// The position at the zero configuration, in metres.
    pub position_m: [f64; 3],
    /// The element, which sets the ideal bond rest length.
    pub element: Element,
    /// The index of the device body that owns this atom.
    ///
    /// The metric bonds two atoms only when they share a body, because each
    /// body is one covalent part.
    pub body: usize,
}

/// The settings of a sub-assembly relaxation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RelaxTarget {
    /// The largest atom count the metric accepts.
    pub max_atoms: usize,
    /// The separation below which two atoms of one body are one bond, in
    /// metres.
    pub bond_cutoff_m: f64,
    /// The separation below which a non-bonded pair is a clash, in metres.
    pub clash_m: f64,
    /// The largest number of minimizer iterations.
    pub max_iterations: usize,
    /// The gradient norm at which the minimizer stops, in newtons.
    pub gradient_tolerance_n: f64,
    /// The initial minimizer step, in metres.
    pub initial_step_m: f64,
}

impl Default for RelaxTarget {
    fn default() -> Self {
        Self {
            max_atoms: 4000,
            bond_cutoff_m: 1.62e-10,
            clash_m: 1.6e-10,
            max_iterations: 600,
            gradient_tolerance_n: 1.0e-13,
            initial_step_m: 1.0e-12,
        }
    }
}

/// The result of a sub-assembly relaxation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RelaxReport {
    /// The number of atoms.
    pub atom_count: usize,
    /// The number of perceived bonds.
    pub bond_count: usize,
    /// The potential energy at the zero configuration, in joules.
    pub initial_energy_j: f64,
    /// The potential energy at the relaxed configuration, in joules.
    pub final_energy_j: f64,
    /// The energy drop, `initial - final`, in joules.
    pub energy_drop_j: f64,
    /// The largest force on one atom at the relaxed configuration, in newtons.
    pub peak_force_n: f64,
    /// The largest absolute bond strain at the relaxed configuration.
    pub max_bond_strain: f64,
    /// The number of non-bonded pairs closer than the clash distance at the
    /// zero configuration.
    pub initial_clash_count: usize,
    /// The number of non-bonded pairs closer than the clash distance at the
    /// relaxed configuration.
    pub clash_count: usize,
    /// The largest displacement of one atom, in metres.
    pub max_displacement_m: f64,
    /// True when the minimizer reached the gradient tolerance.
    pub converged: bool,
}

/// Relaxes a sub-assembly and returns the stability report.
///
/// Returns `None` when the atom count is outside the accepted range, when a
/// term is inconsistent, or when the minimizer fails.
pub fn relax_subassembly(atoms: &[RelaxAtom], target: &RelaxTarget) -> Option<RelaxReport> {
    let atom_count = atoms.len();
    if atom_count < 2 || atom_count > target.max_atoms {
        return None;
    }
    let positions_m: Vec<[f64; 3]> = atoms.iter().map(|atom| atom.position_m).collect();
    let flat: Vec<f64> = positions_m.iter().flat_map(|p| p.iter().copied()).collect();

    let bond_cutoff_squared = target.bond_cutoff_m * target.bond_cutoff_m;
    let mut bonds = BondStretchTerm::new();
    let mut neighbours: Vec<Vec<u32>> = vec![Vec::new(); atom_count];
    for i in 0..atom_count {
        for j in (i + 1)..atom_count {
            if atoms[i].body != atoms[j].body {
                continue;
            }
            let squared = distance_squared(positions_m[i], positions_m[j]);
            if squared > bond_cutoff_squared {
                continue;
            }
            let r0_m = rest_length_m(atoms[i].element, atoms[j].element);
            bonds
                .add_bond(i as u32, j as u32, BOND_STIFFNESS_N_PER_M, r0_m)
                .ok()?;
            neighbours[i].push(j as u32);
            neighbours[j].push(i as u32);
        }
    }
    let bond_count = bonds.bond_count();

    let mut system = System::with_masses_kg(atom_count, &vec![ATOM_MASS_KG; atom_count]).ok()?;
    system.set_bond_stretch(bonds).ok()?;
    system.rebuild_neighbors(&flat).ok()?;
    let initial_energy_j = system.energy_j(&flat).ok()?;

    let mut relaxed_m = flat.clone();
    let options = MinimizeOptions {
        max_iterations: target.max_iterations,
        gradient_tolerance_n: target.gradient_tolerance_n,
        initial_step_m: target.initial_step_m,
    };
    let result = minimize_with(
        &mut system,
        &mut relaxed_m,
        &options,
        MinimizeMethod::ConjugateGradientThenLbfgs,
    )
    .ok()?;
    system.rebuild_neighbors(&relaxed_m).ok()?;
    let final_energy_j = system.energy_j(&relaxed_m).ok()?;

    let forces = system.forces_n(&relaxed_m).ok()?;
    let mut peak_force_n: f64 = 0.0;
    for atom in 0..atom_count {
        let (fx, fy, fz) = (forces[3 * atom], forces[3 * atom + 1], forces[3 * atom + 2]);
        peak_force_n = peak_force_n.max(fx.hypot(fy).hypot(fz));
    }

    let strains = system.bond_strains(&relaxed_m).ok()?;
    let max_bond_strain = strains
        .iter()
        .fold(0.0_f64, |max, strain| max.max(strain.abs()));

    let initial_clash_count = count_clashes(&flat, &neighbours, atom_count, target);
    let clash_count = count_clashes(&relaxed_m, &neighbours, atom_count, target);

    let mut max_displacement_m: f64 = 0.0;
    for atom in 0..atom_count {
        let (dx, dy, dz) = (
            relaxed_m[3 * atom] - flat[3 * atom],
            relaxed_m[3 * atom + 1] - flat[3 * atom + 1],
            relaxed_m[3 * atom + 2] - flat[3 * atom + 2],
        );
        max_displacement_m = max_displacement_m.max(dx.hypot(dy).hypot(dz));
    }

    Some(RelaxReport {
        atom_count,
        bond_count,
        initial_energy_j,
        final_energy_j,
        energy_drop_j: initial_energy_j - final_energy_j,
        peak_force_n,
        max_bond_strain,
        initial_clash_count,
        clash_count,
        max_displacement_m,
        converged: result.converged,
    })
}

/// Relaxes the atoms of one whole part, in the part frame.
pub fn relax_part(part: &nanocad_model::Part, target: &RelaxTarget) -> Option<RelaxReport> {
    let topology = &part.topology;
    let mut atoms = Vec::with_capacity(topology.atom_count());
    for index in 0..topology.atom_count() {
        atoms.push(RelaxAtom {
            position_m: topology.position_m(index)?,
            element: topology.element(index)?,
            body: 0,
        });
    }
    relax_subassembly(&atoms, target)
}

/// Counts the non-bonded pairs closer than the clash distance.
fn count_clashes(
    relaxed_m: &[f64],
    neighbours: &[Vec<u32>],
    atom_count: usize,
    target: &RelaxTarget,
) -> usize {
    let clash_squared = target.clash_m * target.clash_m;
    let mut bonded = vec![false; atom_count * atom_count];
    for (u, list) in neighbours.iter().enumerate() {
        for &v in list {
            bonded[u * atom_count + v as usize] = true;
        }
    }
    let mut clashes = 0;
    for i in 0..atom_count {
        for j in (i + 1)..atom_count {
            if bonded[i * atom_count + j] {
                continue;
            }
            let squared = squared_of_flat(relaxed_m, i, j);
            if squared < clash_squared {
                clashes += 1;
            }
        }
    }
    clashes
}

/// Returns the squared distance between two atoms of a flat position buffer.
fn squared_of_flat(positions_m: &[f64], i: usize, j: usize) -> f64 {
    let (dx, dy, dz) = (
        positions_m[3 * i] - positions_m[3 * j],
        positions_m[3 * i + 1] - positions_m[3 * j + 1],
        positions_m[3 * i + 2] - positions_m[3 * j + 2],
    );
    dx * dx + dy * dy + dz * dz
}

/// Returns the squared distance between two positions.
fn distance_squared(a: [f64; 3], b: [f64; 3]) -> f64 {
    let (dx, dy, dz) = (a[0] - b[0], a[1] - b[1], a[2] - b[2]);
    dx * dx + dy * dy + dz * dz
}

/// Returns the ideal bond rest length of an element pair, in metres.
fn rest_length_m(a: Element, b: Element) -> f64 {
    match (a, b) {
        (Element::HYDROGEN, Element::HYDROGEN) => H_H_BOND_M,
        (Element::HYDROGEN, _) | (_, Element::HYDROGEN) => C_H_BOND_M,
        _ => C_C_BOND_M,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nanocad_parts::{LeafSpringGenerator, PartGenerator};

    #[test]
    fn a_whole_leaf_spring_relaxes_without_a_clash() {
        let part = LeafSpringGenerator
            .generate_with_defaults()
            .expect("the spring builds");
        let report = relax_part(&part, &RelaxTarget::default()).expect("the spring relaxes");
        assert_eq!(report.atom_count, part.atom_count());
        assert!(report.bond_count > 0);
        assert!(report.energy_drop_j >= -1.0e-22, "unexpected energy rise");
        assert!(report.max_bond_strain < 0.1);
        assert_eq!(report.initial_clash_count, 0);
        assert!(report.converged);
    }

    #[test]
    fn a_small_cluster_stays_bonded() {
        let bond_m = 1.544e-10;
        let atoms = vec![
            RelaxAtom {
                position_m: [0.0, 0.0, 0.0],
                element: Element::CARBON,
                body: 0,
            },
            RelaxAtom {
                position_m: [1.04 * bond_m, 0.0, 0.0],
                element: Element::CARBON,
                body: 0,
            },
        ];
        let report = relax_subassembly(&atoms, &RelaxTarget::default()).expect("the pair relaxes");
        assert_eq!(report.bond_count, 1);
        assert_eq!(report.clash_count, 0);
        assert!(report.converged);
    }

    #[test]
    fn an_atom_count_outside_the_range_is_refused() {
        let target = RelaxTarget::default();
        assert!(relax_subassembly(&[], &target).is_none());
        let one = [RelaxAtom {
            position_m: [0.0, 0.0, 0.0],
            element: Element::CARBON,
            body: 0,
        }];
        assert!(relax_subassembly(&one, &target).is_none());
    }
}
