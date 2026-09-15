use crate::error::EngineError;

/// Validates a flat position buffer against an atom count.
pub(crate) fn validate_positions(
    positions_m: &[f64],
    atom_count: usize,
) -> Result<(), EngineError> {
    if positions_m.len() != 3 * atom_count {
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
    Ok(())
}

/// Reads one atom position from a flat buffer.
pub(crate) fn atom_position_m(positions_m: &[f64], atom: u32) -> [f64; 3] {
    let base = atom as usize * 3;
    [
        positions_m[base],
        positions_m[base + 1],
        positions_m[base + 2],
    ]
}

pub(crate) fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub(crate) fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

pub(crate) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub(crate) fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

pub(crate) fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

pub(crate) fn scale(a: [f64; 3], factor: f64) -> [f64; 3] {
    [a[0] * factor, a[1] * factor, a[2] * factor]
}

pub(crate) fn add_assign(target: &mut [f64], atom: u32, value: [f64; 3]) {
    let base = atom as usize * 3;
    target[base] += value[0];
    target[base + 1] += value[1];
    target[base + 2] += value[2];
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vector_helpers_are_consistent() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        assert_eq!(sub(b, a), [3.0, 3.0, 3.0]);
        assert_eq!(add(a, b), [5.0, 7.0, 9.0]);
        assert_eq!(dot(a, b), 32.0);
        assert_eq!(cross(a, b), [-3.0, 6.0, -3.0]);
        assert!((norm([3.0, 4.0, 0.0]) - 5.0).abs() < 1.0e-15);
        assert_eq!(scale(a, 2.0), [2.0, 4.0, 6.0]);
        let mut buffer = vec![0.0; 6];
        add_assign(&mut buffer, 1, [1.0, 2.0, 3.0]);
        assert_eq!(buffer, vec![0.0, 0.0, 0.0, 1.0, 2.0, 3.0]);
    }
}
