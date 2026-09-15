//! Criterion benchmarks for the per-step force and energy cost.
//!
//! The reference case is a periodic SPC-like water box. `spc_water_box` draws
//! velocities from a seed, so the setup is reproducible. The benchmark runs on
//! the release profile, so the timings are not the debug-profile timings.

use criterion::{criterion_group, criterion_main, Criterion};
use nanocad_engine::spc_water_box;
use std::hint::black_box;

const MOLECULE_COUNT: usize = 200;
const WATER_VOLUME_M3_PER_MOLECULE: f64 = 33.4e-30;
const SEED: u64 = 0xB3C0;

fn box_length_m(molecule_count: usize) -> f64 {
    (molecule_count as f64 * WATER_VOLUME_M3_PER_MOLECULE).cbrt()
}

fn bench_water_energy_and_gradient(c: &mut Criterion) {
    let mut water = spc_water_box(MOLECULE_COUNT, box_length_m(MOLECULE_COUNT), 300.0, SEED)
        .expect("valid water box");
    let positions_m = water.positions_m().to_vec();
    c.bench_function("water_200_energy_and_gradient", |b| {
        b.iter(|| {
            water
                .system_mut()
                .energy_and_gradient_j(black_box(&positions_m))
                .expect("valid")
        })
    });
}

fn bench_water_serial_forces(c: &mut Criterion) {
    let mut water = spc_water_box(MOLECULE_COUNT, box_length_m(MOLECULE_COUNT), 300.0, SEED)
        .expect("valid water box");
    let positions_m = water.positions_m().to_vec();
    c.bench_function("water_200_forces_serial", |b| {
        b.iter(|| {
            water
                .system_mut()
                .forces_n(black_box(&positions_m))
                .expect("valid")
        })
    });
}

fn bench_water_parallel_forces(c: &mut Criterion) {
    let mut water = spc_water_box(MOLECULE_COUNT, box_length_m(MOLECULE_COUNT), 300.0, SEED)
        .expect("valid water box");
    let positions_m = water.positions_m().to_vec();
    c.bench_function("water_200_forces_parallel_4", |b| {
        b.iter(|| {
            water
                .system_mut()
                .forces_n_parallel(black_box(&positions_m), 4)
                .expect("valid")
        })
    });
}

criterion_group!(
    benches,
    bench_water_energy_and_gradient,
    bench_water_serial_forces,
    bench_water_parallel_forces
);
criterion_main!(benches);
