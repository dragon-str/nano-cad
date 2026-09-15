# Benchmarks

This directory holds the benchmark harness notes and the recorded results.

The public benchmark page with hardware, timing, and caveats is
[`docs/benchmarks.md`](../docs/benchmarks.md). A reproducible runner is
[`scripts/bench.sh`](../scripts/bench.sh).

## Run the benchmarks

Run every benchmark in the workspace:

```
just bench
```

`just bench` calls `cargo bench`. To run one crate alone:

```
cargo bench -p nanocad-units
```

To compile the benchmarks without running the timing:

```
cargo bench -p nanocad-units --no-run
```

## Where the results go

Criterion writes its report and its raw data to `target/criterion/`. Open
`target/criterion/report/index.html` in a browser to read the report.
Criterion stores the previous run and reports the change.

`benchmarks/results/` holds durable timing records that we keep in the
repository, for example the results of task `M2-12`. Criterion output does not
go here, because it is large and machine-specific. Record only the summary
numbers, the hardware, and the date.
