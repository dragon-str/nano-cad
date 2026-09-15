# Reproduce the demo

This page tells a stranger how to reproduce the nano-cad demo from a clean
checkout. It lists the prerequisites, the exact commands, what each step
proves, and what the demo does **not** reproduce.

The demo is the planetary gearbox. An agent builds it through the MCP tool
surface only. The measured result is the sun-to-carrier gear ratio.

## Prerequisites

| Item | Requirement |
|---|---|
| Operating system | Linux or macOS with a POSIX `/bin/sh`. |
| Rust | `cargo` and `rustc` on `PATH`. Install through <https://rustup.rs>. The workspace needs Rust 1.85 or newer (MSRV in `Cargo.toml`). |
| Python | `python3`, version 3.10 or newer. |
| Network | Needed on the first run only. The script downloads `maturin` and `pytest` from PyPI. |
| Disk | About 1 GB for the Rust `target/` directory and the virtual environment. |

You need no prebuilt artifact. The script builds the extension and the Rust
crates from source.

## Commands

From the repository root:

```sh
just reproduce
```

The recipe calls `scripts/reproduce.sh`. You can run the script directly:

```sh
sh scripts/reproduce.sh
```

The script uses the virtual environment at `<repo>/.venv` by default. Set
`NCAD_VENV` to use another location:

```sh
NCAD_VENV=/tmp/ncvenv sh scripts/reproduce.sh
```

The script exits nonzero when any step fails. It prints one clear line for each
step.

## What each step proves

| Step | Action | It proves |
|---|---|---|
| 1 | Check `cargo`, `rustc`, `python3` | The host has the toolchains. A missing tool stops the run with exit 127. |
| 2 | Create or reuse the Python virtual environment | The demo runs in an isolated environment and does not touch the system Python. |
| 3 | Install `maturin` and `pytest` | The build and test tools are present. This step needs network access on the first run. |
| 4 | Build the extension with `maturin develop` from `python/` | The PyO3 bindings compile and install as `nanocad._core`. The build must run from `python/`, because the `pyproject.toml` there sets `module-name` and `python-source`. |
| 5 | `cargo test --workspace` | Every Rust gate passes: units, model, formats, engine, jigs, parts, params, CLI. This includes the finite-difference gradient tests and the gear-ratio test. |
| 6 | `pytest mcp/tests/test_agent_demo.py` | The agent builds the gearbox through `call_tool` only, with no direct core import. It builds the part, assembles the L2 device, measures the ratio, relaxes a lattice, saves and loads NCZ, and exports URDF. |
| 7 | Print the key result | The measured gear ratio is inside tolerance. |

## Expected result

The last line is the measured result. The numbers below came from a run on
Apple M4 with Python 3.14 and cargo 1.98.1:

```
==> Step 7/7: key measured result
REPRODUCE OK: gearbox demo: analytic_ratio 3.500000000 measured_ratio 3.500000000 relative_error 1.146e-10 tolerance 1.0e-03
The agent built the planetary gearbox through the MCP tool surface only.
```

The demo also exports a URDF with 7 links and 6 joints. Your timing will
differ. The ratio and the URDF structure must not.

## The guard test

The script must not pass silently. `mcp/tests/test_reproduce_script.py` checks
two things:

- `sh -n` accepts the script (valid POSIX syntax).
- The script exits nonzero and reports `cargo` when `cargo` is absent from
  `PATH`.

Run it with:

```sh
python -m pytest mcp/tests/test_reproduce_script.py
```

## What is NOT reproduced

- **The 90-second video (M8-02).** The video and its script are a separate
  artifact. This script reproduces the computation, not the recording.
- **The written report (M8-03).** The report is a separate document.
- **The three-scale export (M8-01).** It is not done yet.
- **The benchmark numbers.** See `docs/benchmarks.md`. Every timing number is
  machine-specific. This script does not run `cargo bench`.
- **The OpenMM cross-check (M2-11).** OpenMM is not installed. The water-box
  self-check runs, but the cross-check waits on an OpenMM install.
- **The public benchmark page.** It reports one machine and one toolchain.

## Known caveats

- **The first run needs network access.** The script downloads `maturin` and
  `pytest`. A later run reuses the virtual environment and skips the download.
- **The `maturin` build must run from `python/`.** A build from the repository
  root with `--manifest-path` misnames the module. The script changes to
  `python/` before it calls maturin.
- **The build is a release build.** `maturin develop --release` takes longer to
  compile than a debug build. The demo then runs fast.
- **`cargo test --workspace` does not build the Python feature.** The `pyo3`
  dependency is optional and stays off for the default Rust gates (ADR-0012).
  The extension is built only in step 4, by maturin.
- **The tree must compile when you run the script.** If another change is
  mid-edit, step 5 fails. Wait for the tree to become green and retry.
- **The script writes `<repo>/.venv` and `target/`.** Both are gitignored.
