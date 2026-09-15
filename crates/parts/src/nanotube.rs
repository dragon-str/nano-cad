use std::collections::BTreeSet;
use std::f64::consts::PI;

use nanocad_model::{Atom, Element, Part, Topology};

use crate::error::PartError;
use crate::generator::PartGenerator;
use crate::geometry::first_shell_bonds;
use crate::parameter::{ParameterSet, ParameterSpec};

/// The graphene in-plane lattice constant at 300 K, in metres.
///
/// Source: P. Trucano and R. Chen, "Structure of graphite by neutron
/// diffraction", Nature 258, 136 (1975).
const GRAPHENE_LATTICE_CONSTANT_M: f64 = 2.461e-10;

/// The graphene carbon-carbon bond length `a / sqrt(3)`, in metres.
///
/// Source: same lattice constant. The ideal value is 1.4209e-10 m.
const GRAPHENE_BOND_M: f64 = 1.420_9e-10;

/// The number of decimal places used to key a chiral-equivalent atom.
const CHIRAL_KEY_SCALE: f64 = 1.0e9;

/// The modulus of the wrapped chiral key. It merges the two faces of the seam.
const CHIRAL_KEY_MOD: i64 = 1_000_000_000;

static NANOTUBE_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new("chiral_n", None, 5.0, 0.0, 40.0, true, "chiral index n"),
    ParameterSpec::new("chiral_m", None, 5.0, 0.0, 40.0, true, "chiral index m"),
    ParameterSpec::new(
        "length_cells",
        None,
        1.0,
        1.0,
        60.0,
        true,
        "translational unit cells along the tube axis",
    ),
];

/// Builds a single-walled carbon nanotube from chiral indices `(n, m)`.
///
/// The chiral vector is `C_h = n a1 + m a2` on the graphene sheet. The
/// translation vector `T` is the shortest lattice vector perpendicular to
/// `C_h`. Atoms in the fundamental parallelogram are rolled onto a cylinder of
/// radius `|C_h| / (2 pi)`. The roll is an isometry, so every graphene bond
/// length is preserved.
#[derive(Clone, Copy, Debug, Default)]
pub struct NanotubeGenerator;

impl NanotubeGenerator {
    /// Returns the analytic cylinder radius `|C_h| / (2 pi)`, in metres.
    pub fn radius_m(n: usize, m: usize) -> f64 {
        let a_m = GRAPHENE_LATTICE_CONSTANT_M;
        let squared = (n * n + n * m + m * m) as f64;
        a_m * squared.sqrt() / (2.0 * PI)
    }

    /// Returns the number of atoms in one translational unit cell.
    ///
    /// It is `4 (n^2 + n m + m^2) / d_R`, where `d_R = gcd(2m+n, 2n+m)`.
    pub fn atoms_per_cell(n: usize, m: usize) -> Result<usize, PartError> {
        let d_r = gcd(2 * m + n, 2 * n + m);
        if d_r == 0 {
            return Err(PartError::InvalidGeometry(
                "chiral indices (0, 0) name no tube".to_owned(),
            ));
        }
        let numerator = 4 * (n * n + n * m + m * m);
        if !numerator.is_multiple_of(d_r) {
            return Err(PartError::InvalidGeometry(
                "chiral atom count is not integral".to_owned(),
            ));
        }
        Ok(numerator / d_r)
    }
}

impl PartGenerator for NanotubeGenerator {
    fn id(&self) -> &'static str {
        "nanotube"
    }

    fn name(&self) -> &'static str {
        "Single-walled carbon nanotube"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        NANOTUBE_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        let resolved = self.resolve(parameters)?;
        let chiral_n = resolved.require("chiral_n")? as usize;
        let chiral_m = resolved.require("chiral_m")? as usize;
        let length_cells = resolved.require("length_cells")? as usize;

        if chiral_n == 0 && chiral_m == 0 {
            return Err(PartError::InvalidGeometry(
                "chiral indices (0, 0) name no tube".to_owned(),
            ));
        }

        let d_r = gcd(2 * chiral_m + chiral_n, 2 * chiral_n + chiral_m);
        let t1 = ((2 * chiral_m + chiral_n) / d_r) as isize;
        let t2_magnitude = (2 * chiral_n + chiral_m) / d_r;
        let t2 = -(t2_magnitude as isize);
        let determinant = chiral_n as isize * t2 - t1 * chiral_m as isize;

        let basis = [(0.0, 0.0), (1.0 / 3.0, 1.0 / 3.0)];
        let search_range = chiral_n + chiral_m + t1.unsigned_abs() + t2_magnitude + 4;
        let mut seen = BTreeSet::new();
        let mut fractional = Vec::new();
        for i in 0..=search_range {
            for j in 0..=search_range {
                for (basis_u, basis_v) in basis {
                    let lattice_u = i as f64 + basis_u;
                    let lattice_v = j as f64 + basis_v;
                    let alpha =
                        (lattice_u * t2 as f64 - t1 as f64 * lattice_v) / determinant as f64;
                    let beta = (chiral_n as f64 * lattice_v - chiral_m as f64 * lattice_u)
                        / determinant as f64;
                    let wrapped_u = alpha - alpha.floor();
                    let wrapped_v = beta - beta.floor();
                    let key = (
                        ((wrapped_u * CHIRAL_KEY_SCALE).round() as i64).rem_euclid(CHIRAL_KEY_MOD),
                        ((wrapped_v * CHIRAL_KEY_SCALE).round() as i64).rem_euclid(CHIRAL_KEY_MOD),
                    );
                    if seen.insert(key) {
                        fractional.push((wrapped_u, wrapped_v));
                    }
                }
            }
        }

        let expected = Self::atoms_per_cell(chiral_n, chiral_m)?;
        if fractional.len() != expected {
            return Err(PartError::InvalidGeometry(format!(
                "chiral search found {} atoms, expected {expected}",
                fractional.len()
            )));
        }

        let radius_m = Self::radius_m(chiral_n, chiral_m);
        let chiral_length_m = GRAPHENE_LATTICE_CONSTANT_M
            * ((chiral_n * chiral_n + chiral_n * chiral_m + chiral_m * chiral_m) as f64).sqrt();
        let translation_length_m = 3.0_f64.sqrt() * chiral_length_m / d_r as f64;

        let mut topology = Topology::new();
        for cell in 0..length_cells {
            for (wrapped_u, wrapped_v) in &fractional {
                let angle_rad = 2.0 * PI * wrapped_u;
                let axial_m = (cell as f64 + wrapped_v) * translation_length_m;
                let position_m = [
                    radius_m * angle_rad.cos(),
                    radius_m * angle_rad.sin(),
                    axial_m,
                ];
                topology.add_atom(Atom::new(Element::CARBON, position_m, 0.0, "C"));
            }
        }

        first_shell_bonds(&mut topology, GRAPHENE_BOND_M)?;

        let name = format!("nanotube-{chiral_n}-{chiral_m}x{length_cells}");
        Ok(Part::new(name, topology).with_material("carbon"))
    }
}

fn gcd(mut a: usize, mut b: usize) -> usize {
    while b != 0 {
        let remainder = a % b;
        a = b;
        b = remainder;
    }
    a
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use super::*;
    use crate::geometry::distance_m;

    fn generate(n: f64, m: f64, cells: f64) -> Part {
        NanotubeGenerator
            .generate(
                &ParameterSet::new()
                    .with("chiral_n", n)
                    .with("chiral_m", m)
                    .with("length_cells", cells),
            )
            .expect("valid nanotube")
    }

    fn max_bond_error(part: &Part, expected_m: f64) -> f64 {
        let mut worst = 0.0_f64;
        for bond in part.topology.bonds() {
            let u_m = part.topology.position_m(bond.u as usize).expect("atom u");
            let v_m = part.topology.position_m(bond.v as usize).expect("atom v");
            let delta_m = (distance_m(u_m, v_m) - expected_m).abs();
            if delta_m > worst {
                worst = delta_m;
            }
        }
        worst
    }

    #[test]
    fn armchair_five_five_has_twenty_atoms_per_cell() {
        assert_eq!(NanotubeGenerator::atoms_per_cell(5, 5), Ok(20));
        let part = generate(5.0, 5.0, 1.0);
        assert_eq!(part.atom_count(), 20);
        assert_eq!(part.material, "carbon");
        assert_eq!(part.name, "nanotube-5-5x1");
    }

    #[test]
    fn zigzag_six_zero_has_twenty_four_atoms_per_cell() {
        assert_eq!(NanotubeGenerator::atoms_per_cell(6, 0), Ok(24));
        let part = generate(6.0, 0.0, 1.0);
        assert_eq!(part.atom_count(), 24);
    }

    #[test]
    fn atom_count_grows_with_length() {
        let short = generate(5.0, 5.0, 1.0);
        let long = generate(5.0, 5.0, 3.0);
        assert!(long.atom_count() > short.atom_count());
        assert_eq!(long.atom_count(), 3 * short.atom_count());
        assert_eq!(long.atom_count() - short.atom_count(), 2 * 20);
    }

    #[test]
    fn every_bond_matches_the_graphene_bond_length() {
        // A large radius keeps the curvature compression of circumferential
        // bonds below the 1e-3 tolerance.
        let part = generate(20.0, 20.0, 2.0);
        assert!(part.bond_count() > 0, "expected bonds");
        let expected_m = GRAPHENE_LATTICE_CONSTANT_M / 3.0_f64.sqrt();
        let worst_m = max_bond_error(&part, expected_m);
        assert!(
            worst_m <= expected_m * 1.0e-3,
            "worst bond error {worst_m:e} m"
        );
    }

    #[test]
    fn a_small_radius_tube_keeps_bonds_within_one_percent() {
        let part = generate(8.0, 4.0, 2.0);
        let expected_m = GRAPHENE_LATTICE_CONSTANT_M / 3.0_f64.sqrt();
        let worst_m = max_bond_error(&part, expected_m);
        assert!(
            worst_m <= expected_m * 1.0e-2,
            "worst bond error {worst_m:e} m"
        );
    }

    #[test]
    fn every_atom_sits_on_the_analytic_radius() {
        for (n, m) in [(5.0, 5.0), (10.0, 0.0), (8.0, 4.0)] {
            let part = generate(n, m, 1.0);
            let radius_m = NanotubeGenerator::radius_m(n as usize, m as usize);
            for atom in part.topology.atoms() {
                let radial_m = (atom.position_m[0] * atom.position_m[0]
                    + atom.position_m[1] * atom.position_m[1])
                    .sqrt();
                assert!(
                    (radial_m - radius_m).abs() <= radius_m * 1.0e-12,
                    "radius {radial_m:e} m differs from {radius_m:e} m"
                );
            }
        }
    }

    #[test]
    fn the_radius_matches_the_chiral_formula() {
        let a_m = GRAPHENE_LATTICE_CONSTANT_M;
        let expected_m = a_m * (64.0_f64 + 32.0 + 16.0).sqrt() / (2.0 * PI);
        assert!((NanotubeGenerator::radius_m(8, 4) - expected_m).abs() <= expected_m * 1.0e-12);
    }

    #[test]
    fn invalid_chiral_indices_are_an_error_not_a_panic() {
        let zero = ParameterSet::new()
            .with("chiral_n", 0.0)
            .with("chiral_m", 0.0);
        assert!(NanotubeGenerator.generate(&zero).is_err());
        let fraction = ParameterSet::new().with("chiral_n", 2.5);
        assert!(NanotubeGenerator.generate(&fraction).is_err());
        let too_long = ParameterSet::new().with("length_cells", 0.0);
        assert!(NanotubeGenerator.generate(&too_long).is_err());
    }

    #[test]
    fn a_generator_reports_its_identity() {
        assert_eq!(NanotubeGenerator.id(), "nanotube");
        assert_eq!(NanotubeGenerator.name(), "Single-walled carbon nanotube");
        assert_eq!(NanotubeGenerator.parameters().len(), 3);
    }
}
