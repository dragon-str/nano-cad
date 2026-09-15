use crate::error::EngineError;
use crate::nonbonded::PeriodicBox;

/// A Verlet neighbor list on a grid cell list.
///
/// The grid cell size is the cutoff plus the skin. Two atoms closer than the
/// cutoff always land in the same cell or in adjacent cells, so a scan of the
/// 27 neighboring cells finds every pair.
///
/// The pair list is flat. Each consecutive pair of `u32` entries, at offsets
/// `2k` and `2k + 1`, is one pair. The first index is smaller than the second,
/// and every unordered pair appears exactly once.
///
/// The list is non-periodic by default. Give it a [`PeriodicBox`] to apply the
/// minimum-image convention in the pair search. The box cutoff plus skin must
/// not exceed half of any periodic length. The non-bonded terms apply the same
/// minimum image, so the list and the energy stay consistent.
#[derive(Clone, Debug, PartialEq)]
pub struct VerletList {
    cutoff_m: f64,
    skin_m: f64,
    cell_size_m: f64,
    periodic_box: PeriodicBox,
    nx: usize,
    ny: usize,
    nz: usize,
    origin_m: [f64; 3],
    axis_cell_size_m: [f64; 3],
    cell_starts: Vec<u32>,
    cell_atoms: Vec<u32>,
    pairs: Vec<u32>,
    reference_positions_m: Vec<f64>,
    valid: bool,
}

impl VerletList {
    /// Creates an empty list for a cutoff and a skin, both in SI metres.
    ///
    /// The cutoff must be positive and finite. The skin must be non-negative
    /// and finite. The list is not periodic. An invalid input returns an error
    /// and never panics.
    pub fn new(cutoff_m: f64, skin_m: f64) -> Result<Self, EngineError> {
        Self::with_periodic_box(cutoff_m, skin_m, PeriodicBox::non_periodic())
    }

    /// Creates an empty list for a cutoff, a skin, and a periodic box.
    ///
    /// The cutoff and skin are in SI metres. The cutoff plus the skin must not
    /// exceed half of any periodic box length, because the minimum image is
    /// ambiguous beyond that.
    pub fn with_periodic_box(
        cutoff_m: f64,
        skin_m: f64,
        periodic_box: PeriodicBox,
    ) -> Result<Self, EngineError> {
        if !cutoff_m.is_finite() || cutoff_m <= 0.0 {
            return Err(EngineError::NonPositiveCutoff { cutoff_m });
        }
        if !skin_m.is_finite() || skin_m < 0.0 {
            return Err(EngineError::InvalidSkin { skin_m });
        }
        let reach_m = cutoff_m + skin_m;
        for length_m in periodic_box.lengths_m() {
            if length_m > 0.0 && reach_m > 0.5 * length_m {
                return Err(EngineError::CutoffExceedsHalfBox {
                    cutoff_m: reach_m,
                    length_m,
                });
            }
        }
        let cell_size_m = reach_m;
        Ok(Self {
            cutoff_m,
            skin_m,
            cell_size_m,
            periodic_box,
            nx: 1,
            ny: 1,
            nz: 1,
            origin_m: [0.0, 0.0, 0.0],
            axis_cell_size_m: [cell_size_m; 3],
            cell_starts: vec![0],
            cell_atoms: Vec::new(),
            pairs: Vec::new(),
            reference_positions_m: Vec::new(),
            valid: false,
        })
    }

    /// Returns the periodic box. It is non-periodic by default.
    pub fn periodic_box(&self) -> &PeriodicBox {
        &self.periodic_box
    }

    /// Builds the pair list for a flat position buffer.
    ///
    /// The buffer holds three `f64` values per atom. Every value must be
    /// finite. An empty buffer is valid and yields an empty pair list.
    pub fn build(&mut self, positions_m: &[f64]) -> Result<(), EngineError> {
        let atom_count = self.validate_positions(positions_m)?;
        if atom_count == 0 {
            self.pairs.clear();
            self.cell_starts = vec![0];
            self.cell_atoms.clear();
            self.nx = 1;
            self.ny = 1;
            self.nz = 1;
            self.origin_m = [0.0, 0.0, 0.0];
            self.axis_cell_size_m = [self.cell_size_m; 3];
            self.reference_positions_m.clear();
            self.valid = true;
            return Ok(());
        }

        self.build_grid(atom_count, positions_m)?;
        self.build_pairs(positions_m);

        self.reference_positions_m.clear();
        self.reference_positions_m.extend_from_slice(positions_m);
        self.valid = true;
        Ok(())
    }

    /// Reports whether the atoms have moved far enough to need a rebuild.
    ///
    /// A rebuild is needed when any atom moves more than half the skin from
    /// its position at the last build. An unbuilt or mismatched list always
    /// needs a rebuild.
    pub fn needs_rebuild(&self, positions_m: &[f64]) -> bool {
        if !self.valid || positions_m.len() != self.reference_positions_m.len() {
            return true;
        }
        let half_skin_m = 0.5 * self.skin_m;
        let limit_sq_m2 = half_skin_m * half_skin_m;
        let atom_count = positions_m.len() / 3;
        for atom in 0..atom_count {
            let base = atom * 3;
            let delta_m = self.periodic_box.minimum_image([
                positions_m[base] - self.reference_positions_m[base],
                positions_m[base + 1] - self.reference_positions_m[base + 1],
                positions_m[base + 2] - self.reference_positions_m[base + 2],
            ]);
            let dx = delta_m[0];
            let dy = delta_m[1];
            let dz = delta_m[2];
            if dx * dx + dy * dy + dz * dz > limit_sq_m2 {
                return true;
            }
        }
        false
    }

    /// Returns the flat pair list.
    pub fn pairs(&self) -> &[u32] {
        &self.pairs
    }

    /// Returns the number of pairs.
    pub fn pair_count(&self) -> usize {
        self.pairs.len() / 2
    }

    /// Returns the cutoff in SI metres.
    pub fn cutoff_m(&self) -> f64 {
        self.cutoff_m
    }

    /// Returns the skin in SI metres.
    pub fn skin_m(&self) -> f64 {
        self.skin_m
    }

    /// Returns the grid cell size in SI metres.
    pub fn cell_size_m(&self) -> f64 {
        self.cell_size_m
    }

    /// Reports whether a pair list has been built.
    pub fn is_valid(&self) -> bool {
        self.valid
    }

    fn validate_positions(&self, positions_m: &[f64]) -> Result<usize, EngineError> {
        let atom_count = positions_m.len() / 3;
        if atom_count * 3 != positions_m.len() {
            return Err(EngineError::PositionBufferSizeMismatch {
                len: positions_m.len(),
                atom_count,
            });
        }
        for (index, value) in positions_m.iter().enumerate() {
            if !value.is_finite() {
                return Err(EngineError::NonFinitePosition { index });
            }
        }
        Ok(atom_count)
    }

    fn build_grid(&mut self, atom_count: usize, positions_m: &[f64]) -> Result<(), EngineError> {
        let (min_m, extent_m) = bounding_box(positions_m);
        let lengths_m = self.periodic_box.lengths_m();

        let (nx, origin_x, cell_x) =
            axis_grid(lengths_m[0], min_m[0], extent_m[0], self.cell_size_m)?;
        let (ny, origin_y, cell_y) =
            axis_grid(lengths_m[1], min_m[1], extent_m[1], self.cell_size_m)?;
        let (nz, origin_z, cell_z) =
            axis_grid(lengths_m[2], min_m[2], extent_m[2], self.cell_size_m)?;
        let cell_count = nx
            .checked_mul(ny)
            .and_then(|value| value.checked_mul(nz))
            .ok_or(EngineError::GridTooLarge)?;
        let starts_len = cell_count.checked_add(1).ok_or(EngineError::GridTooLarge)?;
        if atom_count > u32::MAX as usize {
            return Err(EngineError::TooManyAtoms { atom_count });
        }

        self.nx = nx;
        self.ny = ny;
        self.nz = nz;
        self.origin_m = [origin_x, origin_y, origin_z];
        self.axis_cell_size_m = [cell_x, cell_y, cell_z];

        self.cell_starts = vec![0u32; starts_len];
        let mut counts = vec![0u32; cell_count];
        for atom in 0..atom_count {
            let cell = self.cell_of(atom, positions_m);
            counts[cell] += 1;
        }
        let mut running = 0u32;
        for (start, count) in self.cell_starts.iter_mut().zip(counts.iter()) {
            *start = running;
            running += *count;
        }
        if let Some(last) = self.cell_starts.last_mut() {
            *last = running;
        }

        self.cell_atoms = vec![0u32; atom_count];
        let mut cursor = self.cell_starts[..cell_count].to_vec();
        for atom in 0..atom_count {
            let cell = self.cell_of(atom, positions_m);
            let slot = cursor[cell] as usize;
            self.cell_atoms[slot] = atom as u32;
            cursor[cell] += 1;
        }
        Ok(())
    }

    fn build_pairs(&mut self, positions_m: &[f64]) {
        self.pairs.clear();
        let lengths_m = self.periodic_box.lengths_m();
        let cutoff_sq_m2 = self.cutoff_m * self.cutoff_m;
        for cz in 0..self.nz {
            for cy in 0..self.ny {
                for cx in 0..self.nx {
                    let center = self.flat_cell(cx, cy, cz);
                    let center_atoms = &self.cell_atoms
                        [self.cell_starts[center] as usize..self.cell_starts[center + 1] as usize];
                    let mut neighbors = [0usize; 27];
                    let mut neighbor_count = 0;
                    for dx in -1i64..=1 {
                        for dy in -1i64..=1 {
                            for dz in -1i64..=1 {
                                let (ax, ay, az) = (
                                    neighbor_cell(cx, dx, self.nx, lengths_m[0] > 0.0),
                                    neighbor_cell(cy, dy, self.ny, lengths_m[1] > 0.0),
                                    neighbor_cell(cz, dz, self.nz, lengths_m[2] > 0.0),
                                );
                                let (Some(ax), Some(ay), Some(az)) = (ax, ay, az) else {
                                    continue;
                                };
                                let neighbor = self.flat_cell(ax, ay, az);
                                if neighbors[..neighbor_count].contains(&neighbor) {
                                    continue;
                                }
                                neighbors[neighbor_count] = neighbor;
                                neighbor_count += 1;
                            }
                        }
                    }
                    for &neighbor in &neighbors[..neighbor_count] {
                        let neighbor_atoms = &self.cell_atoms[self.cell_starts[neighbor] as usize
                            ..self.cell_starts[neighbor + 1] as usize];
                        for &i in center_atoms {
                            for &j in neighbor_atoms {
                                if i >= j {
                                    continue;
                                }
                                if self.pair_distance_sq_m2(positions_m, i, j) <= cutoff_sq_m2 {
                                    self.pairs.push(i);
                                    self.pairs.push(j);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    fn cell_of(&self, atom: usize, positions_m: &[f64]) -> usize {
        let base = atom * 3;
        let lengths_m = self.periodic_box.lengths_m();
        let ix = cell_coordinate(
            positions_m[base],
            self.origin_m[0],
            self.axis_cell_size_m[0],
            self.nx,
            lengths_m[0] > 0.0,
        );
        let iy = cell_coordinate(
            positions_m[base + 1],
            self.origin_m[1],
            self.axis_cell_size_m[1],
            self.ny,
            lengths_m[1] > 0.0,
        );
        let iz = cell_coordinate(
            positions_m[base + 2],
            self.origin_m[2],
            self.axis_cell_size_m[2],
            self.nz,
            lengths_m[2] > 0.0,
        );
        self.flat_cell(ix, iy, iz)
    }

    fn pair_distance_sq_m2(&self, positions_m: &[f64], i: u32, j: u32) -> f64 {
        let a = i as usize * 3;
        let b = j as usize * 3;
        let delta_m = self.periodic_box.minimum_image([
            positions_m[a] - positions_m[b],
            positions_m[a + 1] - positions_m[b + 1],
            positions_m[a + 2] - positions_m[b + 2],
        ]);
        delta_m[0] * delta_m[0] + delta_m[1] * delta_m[1] + delta_m[2] * delta_m[2]
    }

    fn flat_cell(&self, ix: usize, iy: usize, iz: usize) -> usize {
        (iz * self.ny + iy) * self.nx + ix
    }
}

fn bounding_box(positions_m: &[f64]) -> ([f64; 3], [f64; 3]) {
    let mut min_m = [f64::INFINITY; 3];
    let mut max_m = [f64::NEG_INFINITY; 3];
    let atom_count = positions_m.len() / 3;
    for atom in 0..atom_count {
        let base = atom * 3;
        for axis in 0..3 {
            let value = positions_m[base + axis];
            if value < min_m[axis] {
                min_m[axis] = value;
            }
            if value > max_m[axis] {
                max_m[axis] = value;
            }
        }
    }
    let extent_m = [
        max_m[0] - min_m[0],
        max_m[1] - min_m[1],
        max_m[2] - min_m[2],
    ];
    (min_m, extent_m)
}

fn grid_extent(extent_m: f64, cell_size_m: f64) -> Result<usize, EngineError> {
    let cells = (extent_m / cell_size_m).floor();
    if !cells.is_finite() || cells < 0.0 {
        return Err(EngineError::GridTooLarge);
    }
    let count = cells as usize;
    count.checked_add(1).ok_or(EngineError::GridTooLarge)
}

/// Returns the cell count, the grid origin, and the axis cell size.
///
/// A periodic axis spans the box and wraps. A non-periodic axis spans the atom
/// bounding box.
fn axis_grid(
    length_m: f64,
    min_m: f64,
    extent_m: f64,
    cell_size_m: f64,
) -> Result<(usize, f64, f64), EngineError> {
    if length_m > 0.0 {
        let count = ((length_m / cell_size_m).floor() as usize).max(1);
        return Ok((count, 0.0, length_m / count as f64));
    }
    let count = grid_extent(extent_m, cell_size_m)?;
    Ok((count, min_m, cell_size_m))
}

fn cell_coordinate(
    value_m: f64,
    origin_m: f64,
    cell_size_m: f64,
    count: usize,
    periodic: bool,
) -> usize {
    let raw = ((value_m - origin_m) / cell_size_m).floor();
    if periodic {
        return (raw as i64).rem_euclid(count as i64) as usize;
    }
    let upper = (count - 1) as f64;
    raw.clamp(0.0, upper) as usize
}

/// Returns the neighbor cell index, wrapped on a periodic axis and `None` off
/// the grid on a non-periodic axis.
fn neighbor_cell(index: usize, delta: i64, count: usize, periodic: bool) -> Option<usize> {
    let raw = index as i64 + delta;
    if periodic {
        return Some(raw.rem_euclid(count as i64) as usize);
    }
    if raw < 0 || raw as usize >= count {
        return None;
    }
    Some(raw as usize)
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Rng(u64);

    impl Rng {
        fn next_f64(&mut self) -> f64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            (x >> 11) as f64 / (1u64 << 53) as f64
        }

        fn in_range(&mut self, low_m: f64, high_m: f64) -> f64 {
            low_m + (high_m - low_m) * self.next_f64()
        }
    }

    fn random_positions(count: usize, side_m: f64, seed: u64) -> Vec<f64> {
        let mut rng = Rng(seed);
        let mut positions_m = Vec::with_capacity(count * 3);
        for _ in 0..count {
            positions_m.push(rng.in_range(0.0, side_m));
            positions_m.push(rng.in_range(0.0, side_m));
            positions_m.push(rng.in_range(0.0, side_m));
        }
        positions_m
    }

    fn plain_distance_sq(positions_m: &[f64], i: u32, j: u32) -> f64 {
        let a = i as usize * 3;
        let b = j as usize * 3;
        let dx = positions_m[a] - positions_m[b];
        let dy = positions_m[a + 1] - positions_m[b + 1];
        let dz = positions_m[a + 2] - positions_m[b + 2];
        dx * dx + dy * dy + dz * dz
    }

    fn brute_force_pairs(positions_m: &[f64], cutoff_m: f64) -> Vec<u32> {
        let atom_count = positions_m.len() / 3;
        let cutoff_sq_m2 = cutoff_m * cutoff_m;
        let mut pairs = Vec::new();
        for i in 0..atom_count {
            for j in (i + 1)..atom_count {
                if plain_distance_sq(positions_m, i as u32, j as u32) <= cutoff_sq_m2 {
                    pairs.push(i as u32);
                    pairs.push(j as u32);
                }
            }
        }
        pairs
    }

    fn normalize_pairs(pairs: &[u32]) -> Vec<(u32, u32)> {
        let mut normalized: Vec<(u32, u32)> = Vec::with_capacity(pairs.len() / 2);
        let mut offset = 0;
        while offset + 1 < pairs.len() {
            let a = pairs[offset];
            let b = pairs[offset + 1];
            normalized.push((a.min(b), a.max(b)));
            offset += 2;
        }
        normalized.sort_unstable();
        normalized
    }

    #[test]
    fn a_zero_cutoff_is_rejected() {
        assert!(matches!(
            VerletList::new(0.0, 1.0e-10),
            Err(EngineError::NonPositiveCutoff { .. })
        ));
    }

    #[test]
    fn a_negative_or_non_finite_skin_is_rejected() {
        assert!(matches!(
            VerletList::new(1.0e-9, -1.0),
            Err(EngineError::InvalidSkin { .. })
        ));
        assert!(matches!(
            VerletList::new(1.0e-9, f64::NAN),
            Err(EngineError::InvalidSkin { .. })
        ));
    }

    #[test]
    fn an_empty_topology_yields_an_empty_pair_list() {
        let mut list = VerletList::new(1.0e-9, 2.0e-10).expect("valid cutoff");
        list.build(&[]).expect("empty positions are valid");
        assert_eq!(list.pair_count(), 0);
        assert!(list.pairs().is_empty());
        assert!(list.is_valid());
    }

    #[test]
    fn a_position_buffer_that_is_not_a_multiple_of_three_is_rejected() {
        let mut list = VerletList::new(1.0e-9, 2.0e-10).expect("valid cutoff");
        assert!(matches!(
            list.build(&[0.0, 0.0]),
            Err(EngineError::PositionBufferSizeMismatch { .. })
        ));
    }

    #[test]
    fn a_non_finite_position_is_rejected() {
        let mut list = VerletList::new(1.0e-9, 2.0e-10).expect("valid cutoff");
        assert!(matches!(
            list.build(&[0.0, 0.0, f64::INFINITY]),
            Err(EngineError::NonFinitePosition { .. })
        ));
    }

    #[test]
    fn the_pair_list_matches_a_brute_force_reference() {
        let cutoff_m = 1.0e-9;
        let skin_m = 2.0e-10;
        let positions_m = random_positions(200, 5.0e-9, 0x9E3779B97F4A7C15);
        let mut list = VerletList::new(cutoff_m, skin_m).expect("valid cutoff");
        list.build(&positions_m).expect("valid positions");

        let expected = normalize_pairs(&brute_force_pairs(&positions_m, cutoff_m));
        let found = normalize_pairs(list.pairs());
        assert_eq!(found, expected);
        assert_eq!(list.pair_count(), expected.len());
    }

    #[test]
    fn the_pair_list_holds_each_pair_once_with_the_smaller_index_first() {
        let cutoff_m = 1.5e-9;
        let positions_m = random_positions(120, 4.0e-9, 0xDEADBEEFCAFEF00D);
        let mut list = VerletList::new(cutoff_m, 1.0e-10).expect("valid cutoff");
        list.build(&positions_m).expect("valid positions");
        let mut offset = 0;
        while offset + 1 < list.pairs().len() {
            assert!(list.pairs()[offset] < list.pairs()[offset + 1]);
            offset += 2;
        }
        let mut normalized = normalize_pairs(list.pairs());
        normalized.dedup();
        assert_eq!(normalized.len(), list.pair_count());
    }

    #[test]
    fn a_large_skin_can_only_add_pairs_not_remove_them() {
        let cutoff_m = 1.0e-9;
        let positions_m = random_positions(80, 3.0e-9, 0x123456789ABCDEF);
        let mut list = VerletList::new(cutoff_m, 0.0).expect("valid cutoff");
        list.build(&positions_m).expect("valid positions");
        let expected = normalize_pairs(list.pairs());

        let mut with_skin = VerletList::new(cutoff_m, 5.0e-10).expect("valid cutoff");
        with_skin.build(&positions_m).expect("valid positions");
        assert_eq!(normalize_pairs(with_skin.pairs()), expected);
    }

    #[test]
    fn a_rebuild_is_needed_when_an_atom_moves_more_than_half_the_skin() {
        let cutoff_m = 1.0e-9;
        let skin_m = 4.0e-10;
        let positions_m = random_positions(10, 2.0e-9, 0xABCDEF0123456789);
        let mut list = VerletList::new(cutoff_m, skin_m).expect("valid cutoff");
        list.build(&positions_m).expect("valid positions");

        assert!(!list.needs_rebuild(&positions_m));

        let mut small_move = positions_m.clone();
        small_move[0] += 1.0e-10;
        assert!(!list.needs_rebuild(&small_move));

        let mut large_move = positions_m.clone();
        large_move[0] += 3.0e-10;
        assert!(list.needs_rebuild(&large_move));
    }

    #[test]
    fn an_unbuilt_list_always_needs_a_rebuild() {
        let list = VerletList::new(1.0e-9, 1.0e-10).expect("valid cutoff");
        assert!(list.needs_rebuild(&[0.0, 0.0, 0.0]));
    }

    #[test]
    fn a_pair_across_a_cell_boundary_is_found() {
        let cutoff_m = 1.0e-9;
        let skin_m = 0.0;
        let cell_size_m = cutoff_m;
        let positions_m = vec![0.0, 0.0, 0.0, 0.5 * cell_size_m, 0.0, 0.0];
        let mut list = VerletList::new(cutoff_m, skin_m).expect("valid cutoff");
        list.build(&positions_m).expect("valid positions");
        assert_eq!(normalize_pairs(list.pairs()), vec![(0, 1)]);
    }

    #[test]
    fn a_pair_beyond_the_cutoff_is_excluded() {
        let cutoff_m = 1.0e-9;
        let positions_m = vec![0.0, 0.0, 0.0, 2.0e-9, 0.0, 0.0];
        let mut list = VerletList::new(cutoff_m, 1.0e-10).expect("valid cutoff");
        list.build(&positions_m).expect("valid positions");
        assert_eq!(list.pair_count(), 0);
    }

    #[test]
    fn a_reach_longer_than_half_the_periodic_box_is_rejected() {
        let periodic_box = PeriodicBox::new([1.0e-9, 1.0e-9, 1.0e-9]).expect("valid box");
        assert!(matches!(
            VerletList::with_periodic_box(0.6e-9, 0.0, periodic_box),
            Err(EngineError::CutoffExceedsHalfBox { .. })
        ));
    }

    #[test]
    fn a_periodic_pair_across_the_boundary_is_found() {
        let cutoff_m = 0.4e-9;
        let periodic_box = PeriodicBox::new([1.0e-9, 1.0e-9, 1.0e-9]).expect("valid box");
        let positions_m = vec![0.05e-9, 0.0, 0.0, 0.95e-9, 0.0, 0.0];
        let mut list =
            VerletList::with_periodic_box(cutoff_m, 0.0, periodic_box).expect("valid cutoff");
        list.build(&positions_m).expect("valid positions");
        assert_eq!(normalize_pairs(list.pairs()), vec![(0, 1)]);
    }

    #[test]
    fn the_periodic_pair_list_matches_a_minimum_image_brute_force_reference() {
        let cutoff_m = 1.0e-9;
        let skin_m = 2.0e-10;
        let periodic_box = PeriodicBox::new([5.0e-9, 5.0e-9, 5.0e-9]).expect("valid box");
        let positions_m = random_positions(200, 5.0e-9, 0xB5297A4D3C1E5F27);
        let mut list =
            VerletList::with_periodic_box(cutoff_m, skin_m, periodic_box).expect("valid cutoff");
        list.build(&positions_m).expect("valid positions");

        let cutoff_sq_m2 = cutoff_m * cutoff_m;
        let atom_count = positions_m.len() / 3;
        let mut expected = Vec::new();
        for i in 0..atom_count {
            for j in (i + 1)..atom_count {
                let a = i * 3;
                let b = j * 3;
                let delta_m = periodic_box.minimum_image([
                    positions_m[a] - positions_m[b],
                    positions_m[a + 1] - positions_m[b + 1],
                    positions_m[a + 2] - positions_m[b + 2],
                ]);
                let r_sq_m2 =
                    delta_m[0] * delta_m[0] + delta_m[1] * delta_m[1] + delta_m[2] * delta_m[2];
                if r_sq_m2 <= cutoff_sq_m2 {
                    expected.push(i as u32);
                    expected.push(j as u32);
                }
            }
        }
        assert_eq!(normalize_pairs(list.pairs()), normalize_pairs(&expected));
    }

    #[test]
    fn a_periodic_rebuild_is_not_needed_when_an_atom_wraps_the_boundary() {
        let cutoff_m = 0.4e-9;
        let skin_m = 1.0e-9;
        let periodic_box = PeriodicBox::new([3.0e-9, 3.0e-9, 3.0e-9]).expect("valid box");
        let positions_m = vec![2.9e-9, 0.0, 0.0];
        let mut list =
            VerletList::with_periodic_box(cutoff_m, skin_m, periodic_box).expect("valid cutoff");
        list.build(&positions_m).expect("valid positions");

        let wrapped_m = vec![0.05e-9, 0.0, 0.0];
        assert!(!list.needs_rebuild(&wrapped_m));
    }
}
