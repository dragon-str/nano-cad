//! Criterion benchmark harness for the units crate.
//!
//! The first benchmark is trivial. It proves the harness compiles and runs.
//! Real timing benchmarks arrive with their tasks in `TASKS.md`.

use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;

/// Convert a length in Angstrom to a length in metre.
fn angstrom_to_metre(length_angstrom: f64) -> f64 {
    length_angstrom * 1.0e-10
}

fn bench_angstrom_to_metre(c: &mut Criterion) {
    c.bench_function("angstrom_to_metre", |b| {
        b.iter(|| angstrom_to_metre(black_box(1.0)))
    });
}

criterion_group!(benches, bench_angstrom_to_metre);
criterion_main!(benches);
