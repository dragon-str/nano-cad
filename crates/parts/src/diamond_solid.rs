//! Solid hydrogen-capped diamond filling for a gear profile.
//!
//! A gear is not a wire outline. It is a solid block of diamond cubic
//! carbon, cut to the involute profile. This module fills the profile with
//! the real diamond lattice, bonds every carbon to its first-shell
//! neighbours, and caps each surface carbon with hydrogen so that every
//! carbon keeps a valence of four.
//!
//! The geometry follows the diamond cubic structure. The conventional
//! lattice constant is `a = 3.567e-10 m`. The first-shell bond length is
//! `a * sqrt(3) / 4 = 1.544e-10 m`, and the tetrahedral bond angle is
//! `acos(-1/3) = 109.4712 deg`.

use std::collections::HashMap;

use nanocad_model::{Atom, Bond, BondType, Element, Topology};

use crate::error::PartError;

/// The diamond conventional cubic lattice constant at 300 K.
///
/// Source: CRC Handbook of Chemistry and Physics; N. W. Ashcroft and
/// N. D. Mermin, Solid State Physics (1976), chapter 4.
pub(crate) const DIAMOND_LATTICE_CONSTANT_M: f64 = 3.567e-10;

/// The first-shell carbon-carbon bond length in diamond.
pub(crate) fn diamond_bond_length_m() -> f64 {
    DIAMOND_LATTICE_CONSTANT_M * (3.0_f64).sqrt() / 4.0
}

/// The carbon-hydrogen bond length of a hydrogen-capped surface.
///
/// Source: the C-H bond length in methane, 1.09e-10 m.
const C_H_BOND_M: f64 = 1.09e-10;

/// The spacing between adjacent diamond (001) atomic planes, `a / 4`.
pub(crate) fn diamond_plane_spacing_m() -> f64 {
    DIAMOND_LATTICE_CONSTANT_M / 4.0
}

/// The eight basis atoms of the diamond cubic cell, in fractional
/// coordinates. Two interpenetrating FCC lattices, offset by (1/4,1/4,1/4).
pub(crate) const BASIS: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0],
    [0.0, 0.5, 0.5],
    [0.5, 0.0, 0.5],
    [0.5, 0.5, 0.0],
    [0.25, 0.25, 0.25],
    [0.25, 0.75, 0.75],
    [0.75, 0.25, 0.75],
    [0.75, 0.75, 0.25],
];

/// The four tetrahedral neighbour directions of a sublattice-A carbon.
/// A sublattice-B carbon uses the negative of each direction.
const TETRAHEDRAL_A: [[f64; 3]; 4] = [
    [1.0, 1.0, 1.0],
    [1.0, -1.0, -1.0],
    [-1.0, 1.0, -1.0],
    [-1.0, -1.0, 1.0],
];

/// Returns the diamond carbon positions for one gear.
///
/// The outline is a closed polygon in the local `xy` plane, in metres. For an
/// external gear the solid is the interior of the outline. For an internal
/// ring the solid is the annulus between the outline and `outer_radius_m`.
///
/// The gear has `layers` atomic (001) planes. Consecutive planes are spaced
/// `a / 4`, and the slab is centred on `z = 0`. Each returned atom sits on a
/// real diamond site.
pub(crate) fn fill_profile(
    outline_m: &[[f64; 2]],
    internal: bool,
    outer_radius_m: Option<f64>,
    layers: usize,
) -> Vec<[f64; 3]> {
    if layers == 0 || outline_m.len() < 3 {
        return Vec::new();
    }
    let a = DIAMOND_LATTICE_CONSTANT_M;
    let quarter_m = a / 4.0;

    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for point in outline_m {
        min_x = min_x.min(point[0]);
        min_y = min_y.min(point[1]);
        max_x = max_x.max(point[0]);
        max_y = max_y.max(point[1]);
    }
    if let Some(radius_m) = outer_radius_m {
        min_x = min_x.min(-radius_m);
        min_y = min_y.min(-radius_m);
        max_x = max_x.max(radius_m);
        max_y = max_y.max(radius_m);
    }

    let min_m = [min_x, min_y, 0.0];
    let max_m = [max_x, max_y, (layers - 1) as f64 * quarter_m];
    crate::lattice_fill::lattice_atoms(min_m, max_m, &mut |point_m| {
        let radius_m = (point_m[0] * point_m[0] + point_m[1] * point_m[1]).sqrt();
        if internal {
            let within_rim = outer_radius_m.is_none_or(|rim| radius_m <= rim);
            within_rim && !point_in_polygon(point_m[0], point_m[1], outline_m)
        } else {
            point_in_polygon(point_m[0], point_m[1], outline_m)
        }
    })
}

/// One free tetrahedral direction of a surface carbon.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct FreeDirection {
    /// The surface carbon.
    pub host: u32,
    /// The slot of the direction, from 0 to 3.
    pub slot: usize,
    /// The unit direction from the carbon outwards.
    pub direction_m: [f64; 3],
}

/// A group at one free direction, keyed by host and slot.
///
/// A free direction outside the plan takes one hydrogen.
pub(crate) type CapPlan = HashMap<(u32, usize), crate::group::FunctionalGroup>;

/// Bonds every first-shell carbon pair in `atoms` and caps the surface.
///
/// The bonds use the diamond first-shell length. A carbon with fewer than
/// four bonds receives one hydrogen for each missing tetrahedral direction,
/// so every carbon reaches a valence of four.
pub(crate) fn bond_and_cap(topology: &mut Topology, atoms: &[u32]) -> Result<Vec<u32>, PartError> {
    bond_and_cap_limited(topology, atoms, MAX_STRAINED_DIMER_M)
}

/// As [`bond_and_cap`], with the longest strained dimer it will pair.
pub(crate) fn bond_and_cap_limited(
    topology: &mut Topology,
    atoms: &[u32],
    max_strained_m: f64,
) -> Result<Vec<u32>, PartError> {
    bond_atoms(topology, atoms)?;
    let mut all: Vec<u32> = atoms.to_vec();
    let mut added: Vec<u32> = Vec::new();
    // Close the notches that no dimer can close. Each round fills the shared
    // sites that would still hold two caps, then plans again.
    for _ in 0..MAX_FILL_ROUNDS {
        let free = free_directions(topology, &all);
        let plan = plan_dimers_limited(topology, &free, max_strained_m);
        if plan.clash_sites_m.is_empty() {
            break;
        }
        let mut fresh = Vec::with_capacity(plan.clash_sites_m.len());
        for site_m in plan.clash_sites_m {
            fresh.push(topology.add_atom(Atom::new(Element::CARBON, site_m, 0.0, "C")));
        }
        all.extend(fresh.iter().copied());
        bond_new_atoms(topology, &fresh, &all)?;
        added.extend(fresh);
    }
    let free = free_directions(topology, &all);
    let plan = plan_dimers_limited(topology, &free, max_strained_m);
    let free = apply_dimers(topology, &free, &plan)?;
    added.extend(cap_free(topology, &free, &CapPlan::new())?);
    Ok(added)
}

/// The C-C length of a surface dimer after the 2x1 reconstruction.
///
/// Source: the hydrogenated diamond C(100)-(2x1):H dimer is about 1.6e-10 m
/// long in DFT and LEED studies (for example Furthmuller et al., Phys. Rev. B
/// 53, 7334, 1996). The value stays below the 1.62e-10 m bond cutoff of the
/// sense audit, so the audit perceives the dimer as a bond.
pub(crate) const SURFACE_DIMER_M: f64 = 1.60e-10;

/// The longest back bond a dimer shift may leave, in metres. A longer bond
/// would fall out of the audit's bond perception, so such a pair is closed by
/// a filled site instead.
const MAX_BACK_BOND_M: f64 = 1.615e-10;

/// The most notch-filling rounds before capping.
const MAX_FILL_ROUNDS: usize = 4;

/// The longest strained dimer the plan accepts before a relaxation, in
/// metres. A longer pair is not paired.
const MAX_STRAINED_DIMER_M: f64 = 2.6e-10;

/// A missing site is filled with carbon only when this many caps would share
/// it. A site with two caps on a face is a dimer site; filling it would grow
/// the face one row per round.
const MIN_FILL_CAPS: usize = 3;

/// One planned dimer: the two hosts, the free-direction index of each into
/// the shared site, and the two moved positions.
pub(crate) type PlannedDimer = (u32, u32, usize, usize, [f64; 3], [f64; 3]);

/// One candidate partner of a host: the partner, the two free-direction
/// indices, and the moved positions of the host and the partner.
type Candidate = (u32, usize, usize, [f64; 3], [f64; 3]);

/// A surface reconstruction plan.
#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct DimerPlan {
    /// The dimers: two hosts, the free-direction index of each into the
    /// shared site, and the two moved positions.
    pub pairs: Vec<PlannedDimer>,
    /// The missing sites that would still take two or more caps.
    pub clash_sites_m: Vec<[f64; 3]>,
}

/// Plans the 2x1 reconstruction of a capped diamond surface.
///
/// On a (001)-like face, one missing lattice site sits one bond from two
/// surface carbons. Capping both would put two hydrogen atoms 0.74e-10 m
/// apart in that site (ADR-0075). Real diamond reconstructs this face 2x1:
/// the two carbons bond to each other, and each keeps one hydrogen.
///
/// The plan builds a graph of hosts that share a site, takes a leaf-first
/// matching (maximum on the chains that a face gives), and keeps each pair
/// whose move to [`SURFACE_DIMER_M`] leaves every back bond below
/// [`MAX_BACK_BOND_M`]. A host with more than two free directions is not
/// paired. A site that would still take two caps is a clash site; the caller
/// fills it with carbon.
#[allow(dead_code)]
pub(crate) fn plan_dimers(topology: &Topology, free: &[FreeDirection]) -> DimerPlan {
    plan_dimers_limited(topology, free, MAX_STRAINED_DIMER_M)
}

/// As [`plan_dimers`], with the longest strained dimer it will pair.
pub(crate) fn plan_dimers_limited(
    topology: &Topology,
    free: &[FreeDirection],
    max_strained_m: f64,
) -> DimerPlan {
    let bond_m = diamond_bond_length_m();
    let site_of = |direction: &FreeDirection| {
        topology.position_m(direction.host as usize).map(|host_m| {
            [
                host_m[0] + bond_m * direction.direction_m[0],
                host_m[1] + bond_m * direction.direction_m[1],
                host_m[2] + bond_m * direction.direction_m[2],
            ]
        })
    };

    let mut free_count: HashMap<u32, usize> = HashMap::new();
    for direction in free {
        *free_count.entry(direction.host).or_insert(0) += 1;
    }
    let mut by_site: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
    for (index, direction) in free.iter().enumerate() {
        if let Some(site_m) = site_of(direction) {
            by_site.entry(site_key(site_m)).or_default().push(index);
        }
    }
    let mut sites: Vec<_> = by_site.into_iter().collect();
    sites.sort_unstable_by_key(|(site, _)| *site);

    // The host graph: an edge joins two eligible hosts that share a site.
    let eligible = |host: u32| matches!(free_count.get(&host), Some(1) | Some(2));
    let mut edges: HashMap<u32, Vec<(u32, usize, usize)>> = HashMap::new();
    for (_, members) in &sites {
        for (i, &left) in members.iter().enumerate() {
            for &right in &members[i + 1..] {
                let (u, v) = (free[left].host, free[right].host);
                if u == v || !eligible(u) || !eligible(v) {
                    continue;
                }
                edges.entry(u).or_default().push((v, left, right));
                edges.entry(v).or_default().push((u, right, left));
            }
        }
    }

    // Neighbour lists of the graph hosts.
    let mut neighbours: HashMap<u32, Vec<u32>> = HashMap::new();
    for bond in topology.bonds() {
        if edges.contains_key(&bond.u) {
            neighbours.entry(bond.u).or_default().push(bond.v);
        }
        if edges.contains_key(&bond.v) {
            neighbours.entry(bond.v).or_default().push(bond.u);
        }
    }
    // A host moves toward its partner along the pair axis. It moves at most
    // half the gap to a full dimer, and less when a back bond would pass
    // MAX_BACK_BOND_M. A step carbon whose back bond lies along the axis moves
    // little, so the dimer there stays strained until a relaxation settles it.
    let allowed_shift = |host: u32, host_m: [f64; 3], unit: [f64; 3], full_m: f64| {
        let back = |shift_m: f64| {
            let moved = [
                host_m[0] + shift_m * unit[0],
                host_m[1] + shift_m * unit[1],
                host_m[2] + shift_m * unit[2],
            ];
            neighbours.get(&host).is_none_or(|list| {
                list.iter().all(|&n| {
                    topology
                        .position_m(n as usize)
                        .is_none_or(|p| distance_m(p, moved) <= MAX_BACK_BOND_M)
                })
            })
        };
        if back(full_m) {
            return full_m;
        }
        let (mut low, mut high) = (0.0_f64, full_m);
        for _ in 0..30 {
            let mid = 0.5 * (low + high);
            if back(mid) {
                low = mid;
            } else {
                high = mid;
            }
        }
        low
    };
    let feasible = |u: u32, v: u32| -> Option<([f64; 3], [f64; 3])> {
        let u_m = topology.position_m(u as usize)?;
        let v_m = topology.position_m(v as usize)?;
        let length_m = distance_m(u_m, v_m);
        if length_m <= SURFACE_DIMER_M {
            return None;
        }
        let unit = [
            (v_m[0] - u_m[0]) / length_m,
            (v_m[1] - u_m[1]) / length_m,
            (v_m[2] - u_m[2]) / length_m,
        ];
        let full_m = 0.5 * (length_m - SURFACE_DIMER_M);
        let shift_u = allowed_shift(u, u_m, unit, full_m);
        let shift_v = allowed_shift(v, v_m, [-unit[0], -unit[1], -unit[2]], full_m);
        if length_m - shift_u - shift_v > max_strained_m {
            return None;
        }
        let mut new_u = u_m;
        let mut new_v = v_m;
        for axis in 0..3 {
            new_u[axis] += shift_u * unit[axis];
            new_v[axis] -= shift_v * unit[axis];
        }
        Some((new_u, new_v))
    };

    // Leaf-first greedy matching over the feasible edges only.
    let mut open: HashMap<u32, Vec<Candidate>> = HashMap::new();
    for (&host, list) in &edges {
        for &(other, own, theirs) in list {
            if host < other {
                if let Some((new_host, new_other)) = feasible(host, other) {
                    open.entry(host)
                        .or_default()
                        .push((other, own, theirs, new_host, new_other));
                    open.entry(other)
                        .or_default()
                        .push((host, theirs, own, new_other, new_host));
                }
            }
        }
    }
    let mut hosts: Vec<u32> = open.keys().copied().collect();
    hosts.sort_unstable();
    let mut matched: std::collections::HashSet<u32> = std::collections::HashSet::new();
    let degree = |host: u32, matched: &std::collections::HashSet<u32>| {
        open.get(&host)
            .map(|list| list.iter().filter(|e| !matched.contains(&e.0)).count())
            .unwrap_or(0)
    };
    let mut plan = DimerPlan::default();
    loop {
        let mut best: Option<(usize, u32)> = None;
        for &host in &hosts {
            if matched.contains(&host) {
                continue;
            }
            let d = degree(host, &matched);
            if d > 0 && best.is_none_or(|(b, _)| d < b) {
                best = Some((d, host));
                if d == 1 {
                    break;
                }
            }
        }
        let Some((_, host)) = best else {
            break;
        };
        let Some(&(other, own, theirs, new_host, new_other)) = open[&host]
            .iter()
            .filter(|e| !matched.contains(&e.0))
            .min_by_key(|e| (degree(e.0, &matched), e.0))
        else {
            break;
        };
        matched.insert(host);
        matched.insert(other);
        plan.pairs
            .push((host, other, own, theirs, new_host, new_other));
    }

    // The sites that still take two or more caps after the dimers.
    let mut dropped: std::collections::HashSet<usize> = std::collections::HashSet::new();
    for pair in &plan.pairs {
        dropped.insert(pair.2);
        dropped.insert(pair.3);
    }
    for (site, members) in &sites {
        let remaining = members
            .iter()
            .filter(|&&index| !dropped.contains(&index) && !matched.contains(&free[index].host))
            .count();
        if remaining >= MIN_FILL_CAPS {
            let quarter_m = DIAMOND_LATTICE_CONSTANT_M / 4.0;
            plan.clash_sites_m.push([
                site.0 as f64 * quarter_m,
                site.1 as f64 * quarter_m,
                site.2 as f64 * quarter_m,
            ]);
        }
    }
    plan
}

/// Applies a dimer plan: moves and bonds each pair, drops the directions into
/// the shared sites, and re-aims the last free direction of each dimer carbon
/// so it completes the tetrahedron of its three carbon bonds.
///
/// Returns the free directions that remain to be capped.
pub(crate) fn apply_dimers(
    topology: &mut Topology,
    free: &[FreeDirection],
    plan: &DimerPlan,
) -> Result<Vec<FreeDirection>, PartError> {
    let mut dropped: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let mut hosts: std::collections::HashSet<u32> = std::collections::HashSet::new();
    for &(u, v, own, theirs, new_u, new_v) in &plan.pairs {
        topology.set_position_m(u as usize, new_u)?;
        topology.set_position_m(v as usize, new_v)?;
        topology.add_bond(Bond::new(u, v, 1, BondType::Single))?;
        dropped.insert(own);
        dropped.insert(theirs);
        hosts.insert(u);
        hosts.insert(v);
    }
    let mut neighbours: HashMap<u32, Vec<u32>> = HashMap::new();
    for bond in topology.bonds() {
        if hosts.contains(&bond.u) {
            neighbours.entry(bond.u).or_default().push(bond.v);
        }
        if hosts.contains(&bond.v) {
            neighbours.entry(bond.v).or_default().push(bond.u);
        }
    }
    let mut out = Vec::with_capacity(free.len());
    for (index, direction) in free.iter().enumerate() {
        if dropped.contains(&index) {
            continue;
        }
        let mut direction = *direction;
        if hosts.contains(&direction.host) {
            if let Some(host_m) = topology.position_m(direction.host as usize) {
                let mut sum = [0.0_f64; 3];
                for &n in neighbours.get(&direction.host).into_iter().flatten() {
                    if let Some(p) = topology.position_m(n as usize) {
                        let d = distance_m(p, host_m);
                        for axis in 0..3 {
                            sum[axis] += (p[axis] - host_m[axis]) / d;
                        }
                    }
                }
                let norm = (sum[0] * sum[0] + sum[1] * sum[1] + sum[2] * sum[2]).sqrt();
                if norm > 1.0e-12 {
                    direction.direction_m = [-sum[0] / norm, -sum[1] / norm, -sum[2] / norm];
                }
            }
        }
        out.push(direction);
    }
    Ok(out)
}

/// The integer key of a diamond site on the quarter-cell grid.
pub(crate) fn site_key(point_m: [f64; 3]) -> (i64, i64, i64) {
    let quarter_m = DIAMOND_LATTICE_CONSTANT_M / 4.0;
    (
        (point_m[0] / quarter_m).round() as i64,
        (point_m[1] / quarter_m).round() as i64,
        (point_m[2] / quarter_m).round() as i64,
    )
}

/// Bonds each new carbon to every first-shell carbon of `all`.
fn bond_new_atoms(topology: &mut Topology, fresh: &[u32], all: &[u32]) -> Result<(), PartError> {
    let bond_m = diamond_bond_length_m();
    let tolerance_m = bond_m * 1.0e-3;
    let fresh_set: std::collections::HashSet<u32> = fresh.iter().copied().collect();
    let mut grid: HashMap<(i64, i64, i64), Vec<u32>> = HashMap::new();
    for &index in all {
        if let Some(position_m) = topology.position_m(index as usize) {
            grid.entry(cell_of(position_m, bond_m))
                .or_default()
                .push(index);
        }
    }
    for &u in fresh {
        let Some(u_m) = topology.position_m(u as usize) else {
            continue;
        };
        let (cx, cy, cz) = cell_of(u_m, bond_m);
        let mut partners = Vec::new();
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    for &v in grid.get(&(cx + dx, cy + dy, cz + dz)).into_iter().flatten() {
                        // A pair of two new atoms is bonded once, from the lower index.
                        if v == u || (fresh_set.contains(&v) && v < u) {
                            continue;
                        }
                        if topology
                            .position_m(v as usize)
                            .is_some_and(|v_m| (distance_m(u_m, v_m) - bond_m).abs() <= tolerance_m)
                        {
                            partners.push(v);
                        }
                    }
                }
            }
        }
        for v in partners {
            topology.add_bond(Bond::new(u, v, 1, BondType::Single))?;
        }
    }
    Ok(())
}

/// Adds one bond for every first-shell carbon pair in `atoms`.
///
/// The pass adds no cap. Call [`free_directions`] and then [`cap_free`] to
/// fill the surface, or use [`bond_and_cap`] to do both.
pub(crate) fn bond_atoms(topology: &mut Topology, atoms: &[u32]) -> Result<(), PartError> {
    let bond_m = diamond_bond_length_m();
    let tolerance_m = bond_m * 1.0e-3;
    let cell_m = bond_m;

    let mut grid: HashMap<(i64, i64, i64), Vec<u32>> = HashMap::new();
    for &index in atoms {
        if let Some(position_m) = topology.position_m(index as usize) {
            grid.entry(cell_of(position_m, cell_m))
                .or_default()
                .push(index);
        }
    }

    for &u in atoms {
        let position_u_m = match topology.position_m(u as usize) {
            Some(position) => position,
            None => continue,
        };
        let (cx, cy, cz) = cell_of(position_u_m, cell_m);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let Some(neighbours) = grid.get(&(cx + dx, cy + dy, cz + dz)) else {
                        continue;
                    };
                    for &v in neighbours {
                        if v <= u {
                            continue;
                        }
                        let Some(position_v_m) = topology.position_m(v as usize) else {
                            continue;
                        };
                        if (distance_m(position_u_m, position_v_m) - bond_m).abs() <= tolerance_m {
                            topology.add_bond(Bond::new(u, v, 1, BondType::Single))?;
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

/// Returns every free tetrahedral direction of the atoms in `atoms`.
///
/// A direction is free when the carbon has fewer than four bonds inside
/// `atoms`, and when no atom of `atoms` already sits one bond along that
/// direction.
pub(crate) fn free_directions(topology: &Topology, atoms: &[u32]) -> Vec<FreeDirection> {
    let bond_m = diamond_bond_length_m();
    let tolerance_m = bond_m * 1.0e-3;
    let cell_m = bond_m;
    let members: std::collections::HashSet<u32> = atoms.iter().copied().collect();

    let mut grid: HashMap<(i64, i64, i64), Vec<u32>> = HashMap::new();
    for &index in atoms {
        if let Some(position_m) = topology.position_m(index as usize) {
            grid.entry(cell_of(position_m, cell_m))
                .or_default()
                .push(index);
        }
    }

    // One pass over the bonds gives the degree of every member. A per-atom
    // scan would be quadratic, and a gear has more than a hundred thousand
    // atoms.
    let mut degree: HashMap<u32, u32> = HashMap::new();
    for bond in topology.bonds() {
        if members.contains(&bond.u) && members.contains(&bond.v) {
            *degree.entry(bond.u).or_insert(0) += 1;
            *degree.entry(bond.v).or_insert(0) += 1;
        }
    }

    let mut free = Vec::new();
    for &index in atoms {
        let Some(position_m) = topology.position_m(index as usize) else {
            continue;
        };
        if *degree.get(&index).unwrap_or(&0) >= 4 {
            continue;
        }
        for (slot, direction) in tetrahedral_directions(position_m, DIAMOND_LATTICE_CONSTANT_M)
            .iter()
            .enumerate()
        {
            // An existing carbon neighbour sits at the C-C bond length along
            // this direction, so probe there, not at the C-H length.
            let neighbour_m = [
                position_m[0] + bond_m * direction[0],
                position_m[1] + bond_m * direction[1],
                position_m[2] + bond_m * direction[2],
            ];
            // The probe point can sit exactly on a cell boundary, where
            // floating error moves it into a neighbouring cell. Scan the
            // surrounding cells, as the bonding pass does.
            let (px, py, pz) = cell_of(neighbour_m, cell_m);
            let mut occupied = false;
            'cells: for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let Some(neighbours) = grid.get(&(px + dx, py + dy, pz + dz)) else {
                            continue;
                        };
                        for &other in neighbours {
                            if other == index {
                                continue;
                            }
                            if topology.position_m(other as usize).is_some_and(|position| {
                                distance_m(position, neighbour_m) <= tolerance_m
                            }) {
                                occupied = true;
                                break 'cells;
                            }
                        }
                    }
                }
            }
            if occupied {
                continue;
            }
            free.push(FreeDirection {
                host: index,
                slot,
                direction_m: *direction,
            });
        }
    }
    free
}

/// Fills every free direction with a group from `plan`, or one hydrogen.
///
/// The returned indices are the atoms that fill the directions, in the order
/// of `free`.
pub(crate) fn cap_free(
    topology: &mut Topology,
    free: &[FreeDirection],
    plan: &CapPlan,
) -> Result<Vec<u32>, PartError> {
    let mut added = Vec::new();
    for direction in free {
        let Some(host_position_m) = topology.position_m(direction.host as usize) else {
            continue;
        };
        if let Some(group) = plan.get(&(direction.host, direction.slot)) {
            let group_atoms = crate::group::add_group(
                topology,
                direction.host,
                host_position_m,
                direction.direction_m,
                *group,
            )?;
            added.extend(group_atoms);
            continue;
        }
        let candidate_m = [
            host_position_m[0] + C_H_BOND_M * direction.direction_m[0],
            host_position_m[1] + C_H_BOND_M * direction.direction_m[1],
            host_position_m[2] + C_H_BOND_M * direction.direction_m[2],
        ];
        let hydrogen = topology.add_atom(Atom::new(Element::HYDROGEN, candidate_m, 0.0, "H"));
        topology.add_bond(Bond::new(direction.host, hydrogen, 1, BondType::Single))?;
        added.push(hydrogen);
    }
    Ok(added)
}

/// The four tetrahedral directions of the carbon at `position_m`.
fn tetrahedral_directions(position_m: [f64; 3], a: f64) -> [[f64; 3]; 4] {
    let inv = 1.0 / 3.0_f64.sqrt();
    let fx = position_m[0] / a;
    let fy = position_m[1] / a;
    let fz = position_m[2] / a;
    let sum = fx + fy + fz;
    let residual = sum - sum.floor();
    let is_sublattice_a = residual < 0.125 || residual > 0.875;
    let mut directions = [[0.0; 3]; 4];
    for (slot, direction) in directions.iter_mut().zip(TETRAHEDRAL_A) {
        let sign = if is_sublattice_a { 1.0 } else { -1.0 };
        *slot = [
            sign * direction[0] * inv,
            sign * direction[1] * inv,
            sign * direction[2] * inv,
        ];
    }
    directions
}

fn point_in_polygon(x: f64, y: f64, polygon: &[[f64; 2]]) -> bool {
    let mut inside = false;
    let count = polygon.len();
    let mut j = count - 1;
    for i in 0..count {
        let (xi, yi) = (polygon[i][0], polygon[i][1]);
        let (xj, yj) = (polygon[j][0], polygon[j][1]);
        let crosses = (yi > y) != (yj > y) && x < (xj - xi) * (y - yi) / (yj - yi) + xi;
        if crosses {
            inside = !inside;
        }
        j = i;
    }
    inside
}

fn cell_of(position_m: [f64; 3], cell_m: f64) -> (i64, i64, i64) {
    (
        (position_m[0] / cell_m).floor() as i64,
        (position_m[1] / cell_m).floor() as i64,
        (position_m[2] / cell_m).floor() as i64,
    )
}

fn distance_m(a_m: [f64; 3], b_m: [f64; 3]) -> f64 {
    let dx = a_m[0] - b_m[0];
    let dy = a_m[1] - b_m[1];
    let dz = a_m[2] - b_m[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice_fill::fill_solid;
    use crate::shape::Box3;

    /// Builds a capped diamond block and returns its topology.
    fn capped_block(size_cells: [f64; 3]) -> Topology {
        let a = DIAMOND_LATTICE_CONSTANT_M;
        // A carbon with one lattice bond is a loose corner group. The pass
        // does not pair it, so remove it first, as the gear pump does.
        let sites = crate::gear_pump::keep_well_bonded(fill_solid(&Box3 {
            min_m: [0.0, 0.0, 0.0],
            max_m: [size_cells[0] * a, size_cells[1] * a, size_cells[2] * a],
        }));
        let mut topology = Topology::new();
        let atoms: Vec<u32> = sites
            .into_iter()
            .map(|site| topology.add_atom(Atom::new(Element::CARBON, site, 0.0, "C")))
            .collect();
        bond_and_cap(&mut topology, &atoms).expect("the block caps");
        topology
    }

    fn close_pairs(topology: &Topology, limit_m: f64) -> usize {
        let count = topology.atom_count();
        let mut close = 0;
        for i in 0..count {
            for j in (i + 1)..count {
                let (Some(p), Some(q)) = (topology.position_m(i), topology.position_m(j)) else {
                    continue;
                };
                if distance_m(p, q) < limit_m {
                    close += 1;
                }
            }
        }
        close
    }

    /// A (001) block has flat faces where two caps share one site. The 2x1
    /// reconstruction leaves no pair closer than the 0.9e-10 m audit limit.
    #[test]
    fn a_reconstructed_block_has_no_cap_clash() {
        let topology = capped_block([3.0, 3.0, 2.0]);
        assert_eq!(close_pairs(&topology, 0.9e-10), 0);
    }

    /// Every carbon still has four bonds after the reconstruction.
    #[test]
    fn every_carbon_keeps_four_bonds() {
        let topology = capped_block([3.0, 3.0, 2.0]);
        let mut degree = vec![0usize; topology.atom_count()];
        for bond in topology.bonds() {
            degree[bond.u as usize] += 1;
            degree[bond.v as usize] += 1;
        }
        for (index, element) in topology.elements().iter().enumerate() {
            if *element == Element::CARBON {
                assert_eq!(degree[index], 4, "carbon {index}");
            }
        }
    }

    /// The flat faces carry surface dimers near the 1.6e-10 m target.
    #[test]
    fn the_faces_carry_dimers() {
        let topology = capped_block([3.0, 3.0, 2.0]);
        let dimers = topology
            .bonds()
            .filter(|bond| {
                let (Some(p), Some(q)) = (
                    topology.position_m(bond.u as usize),
                    topology.position_m(bond.v as usize),
                ) else {
                    return false;
                };
                (distance_m(p, q) - SURFACE_DIMER_M).abs() < 1.0e-13
            })
            .count();
        assert!(dimers > 0);
    }
}
