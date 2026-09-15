# Benchmarks

This page reports the engine timings for nano-cad. It states the hardware, the
toolchain, the build profile, the exact command, and what each metric measures.
You can reproduce every number on your own machine. Your numbers will differ.

## Scope and caveat

- All numbers below were measured on the machine that the record names. They
  are **not** universal. Do not compare them across machines.
- `cargo bench` is **timing only**. It is not a merge gate.
- The parallel force path is **opt-in**. The default engine path is serial.
- The build profile is `bench`, which is the release profile with thin LTO and
  one codegen unit.

## Machine and toolchain

The recorded run is `benchmarks/results/2026-09-14-engine-timings.txt`. It is
the M2-12 engine timing record.

| Item | Value | Source |
|---|---|---|
| Date (UTC) | 2026-09-15T00:40:09Z | `2026-09-14-engine-timings.txt` |
| Host | Apple M4, 10 logical cores | `2026-09-14-engine-timings.txt` |
| OS | macOS 26.6.2 (build 25G83) | `2026-09-14-engine-timings.txt` |
| Build profile | `bench` (release, LTO thin, codegen-units 1) | `2026-09-14-engine-timings.txt` |
| Bench tag | `nanocad-engine v0.1.0` | `2026-09-14-engine-timings.txt` |
| Toolchain | cargo 1.98.1; rustc 1.98.1 | `20260915T010005Z-engine-timings.txt` |

The M2-12 record did not capture the toolchain. The reproduction record
`benchmarks/results/20260915T010005Z-engine-timings.txt` captured it. That
reproduction ran on the same M4 host and reports cargo 1.98.1 (797e8a9bc
2026-08-05) and rustc 1.98.1 (48a229cea 2026-09-01).

## How to reproduce

Run the whole workspace benchmark set:

```
just bench
```

Run only the engine benchmark:

```
cargo bench -p nanocad-engine
```

Run the engine benchmark with the recorded settings and write a timestamped
record under `benchmarks/results/`:

```
./scripts/bench.sh
```

The script uses POSIX shell. It stops with an error when `cargo` is not on the
PATH. It writes one file per run. The file name starts with the UTC timestamp
`YYYYMMDDTHHMMSSZ`.

Criterion writes its full report to `target/criterion/`. Open
`target/criterion/report/index.html` in a browser. The record files hold the
summary only, because the raw report is large and machine-specific.

## Reference case

Every engine benchmark uses the same reference case.

- A periodic SPC-like water box.
- 200 molecules, 600 atoms.
- Box side 1.884e-9 m.
- Temperature 300 K.
- Seed `0xB3C0` for the velocity draw, so the setup is reproducible.
- The neighbor pair list holds 30827 pairs.

Source: `2026-09-14-engine-timings.txt` and
`crates/engine/benches/engine_bench.rs`.

## Metric definitions

Each benchmark times one full call. The timer covers the call only. It does not
cover the box construction.

| Benchmark | Function | What it times |
|---|---|---|
| `water_200_energy_and_gradient` | `System::energy_and_gradient_j` | One pass that returns the total potential energy in joules and the full gradient in newtons. It sums the bond, angle, torsion, out-of-plane, van der Waals, and electrostatic terms. |
| `water_200_forces_serial` | `System::forces_n` | The total force vector in newtons, the negative of the gradient, on the default serial path. |
| `water_200_forces_parallel_4` | `System::forces_n_parallel` | The same total force vector, with a deterministic parallel reduction over 4 workers. The result is bit-identical run to run and independent of the thread count. |

`forces_n` computes the gradient and then negates it. Therefore the serial
forces time and the energy-and-gradient time are almost the same number. The
parallel benchmark measures the neighbor-pair-list reduction; the bonded terms
stay serial.

## Measured numbers

Machine: Apple M4, 10 logical cores. Profile: `bench`. Command:

```
cargo bench -p nanocad-engine --bench engine_bench -- --warm-up-time 0.5 --measurement-time 1.0 --sample-size 20
```

Criterion reports a mean and a 95% confidence interval `[lower mean upper]`.
The 95% level is Criterion's default. All times are microseconds (us).

| Benchmark | Mean | Lower bound | Upper bound |
|---|---|---|---|
| `water_200_energy_and_gradient` | 680.97 us | 677.92 us | 685.88 us |
| `water_200_forces_serial` | 679.83 us | 678.45 us | 681.14 us |
| `water_200_forces_parallel_4` | 370.97 us | 369.07 us | 373.05 us |

Source: `benchmarks/results/2026-09-14-engine-timings.txt`.

The 4-worker deterministic reduction is about 1.83x faster than the serial path
on this host. Source: the note in `2026-09-14-engine-timings.txt`.

## Variance

- Criterion runs 20 samples in the recorded run. The 95% confidence intervals
  above show the spread.
- The intervals are narrow for the energy and serial benchmarks, about 1% of
  the mean (computed from the intervals above).
- The parallel benchmark varies more, because thread scheduling and the
  machine load change between runs.
- A second run on the same host gave means of 653.02 us (energy), 660.21 us
  (serial), and 303.08 us (parallel). Source:
  `benchmarks/results/20260915T010005Z-engine-timings.txt`. The parallel mean
  moved from 370.97 us to 303.08 us, a decrease of about 18% (computed from the
  two cited means). Treat the parallel number as the least stable.
- Close other programs before you compare two runs.

## Traceability

| Number | Value | File |
|---|---|---|
| Energy mean | 680.97 us | `benchmarks/results/2026-09-14-engine-timings.txt` |
| Energy interval | 677.92 to 685.88 us | `benchmarks/results/2026-09-14-engine-timings.txt` |
| Serial mean | 679.83 us | `benchmarks/results/2026-09-14-engine-timings.txt` |
| Serial interval | 678.45 to 681.14 us | `benchmarks/results/2026-09-14-engine-timings.txt` |
| 4-worker mean | 370.97 us | `benchmarks/results/2026-09-14-engine-timings.txt` |
| 4-worker interval | 369.07 to 373.05 us | `benchmarks/results/2026-09-14-engine-timings.txt` |
| Speedup | about 1.83x | `benchmarks/results/2026-09-14-engine-timings.txt` |
| Host and OS | Apple M4, 10 cores, macOS 26.6.2 | `benchmarks/results/2026-09-14-engine-timings.txt` |
| Atoms, box, pairs | 600 atoms, 1.884e-9 m, 30827 pairs | `benchmarks/results/2026-09-14-engine-timings.txt` |
| Toolchain | cargo 1.98.1, rustc 1.98.1 | `benchmarks/results/20260915T010005Z-engine-timings.txt` |
| Reproduction means | 653.02, 660.21, 303.08 us | `benchmarks/results/20260915T010005Z-engine-timings.txt` |

**Your numbers will differ.** This page reports one machine and one toolchain.
