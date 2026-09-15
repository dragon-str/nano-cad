# ARCHITECTURE.md

## Overview

The system has one headless core and five model levels. Every level reads and
writes the parameter schema in `PARAMETERS.md`. The core serves three clients:
the GUI, the CLI, and the agent.

```
L0 Quantum adapter      PySCF / xTB via ASE
L1 Atomistic engine     our Rust engine (parts, forces, jigs)
L2 Device layer         rigid bodies, joints, parameter handoff
L3 Continuum adapter    external FEA / flow solvers
L4 System layer         lumped ODE models (CellML / SBML interop)

Core: document model, part generators, constraints, units, formats,
      parameter store, Python API, CLI, MCP server
```

## Repo layout

```
nano-cad/
  AGENTS.md README.md PLAN.md ARCHITECTURE.md PARAMETERS.md TASKS.md DECISIONS.md
  Cargo.toml                 # workspace
  justfile                   # task runner: just verify, just bench
  crates/
    units/                   # SI units and conversions
    model/                   # document model: atoms, bonds, topology, parts
    format/                  # NCZ, MMP, XYZ, PDB, URDF
    engine/                  # L1: cell list, pair list, forces, integrators
    jigs/                    # anchors, motors, springs, rigid bodies, device, URDF, scene, nanomedicine
    parts/                   # parametric generators: gear, nanotube, lattice, respirocyte
    params/                  # parameter store, provenance, verification
    python/                  # PyO3 bindings -> the nanocad package
    cli/                     # the nanocad command line
  python/
    nanocad/                 # the Python package and the L0/L3/L4 adapters
    tests/
  mcp/                       # stdio JSON-RPC MCP server wrapping the API
  scripts/                   # reproducible shell entry points (bench, ...)
  tests/
    golden/                  # regression inputs and outputs
  benchmarks/
  docs/
```

## Data model

- `Document`: versioned, serializable design graph. Immutable updates give free
  undo and redo. It is separate from simulation state.
- `Part`: a named collection of atoms and bonds, plus a material and metadata.
- `Assembly`: parts placed by transforms, with constraints and mates.
- `Atom`: element, position, charge, type. Stored in struct-of-arrays form.
- `Bond`: two atom indices, order, type.
- `Constraint` / `Jig`: anchor, motor, spring, thermostat, rigid body, joint.
- `Quantity`: value plus unit. Internally SI.

Topology uses indices, never pointers. Positions and forces are contiguous
arrays. Pair lists are flat 32-bit index arrays.

## L1 engine

- Cell list on a grid, plus a Verlet skin. Rebuild the pair list when any atom
  moves more than half the skin.
- Flat pair list of `u32` indices for cache locality.
- Force terms: bond stretch, angle bend, torsion, out-of-plane, van der Waals
  (Buckingham or MM4), electrostatic with a cutoff.
- Integrators: velocity Verlet. Thermostats: Langevin and Berendsen.
- Minimizer: conjugate gradient, then L-BFGS near the minimum.
- Parallel: multi-threaded force accumulation with per-thread buffers and a
  reduction. Provide a deterministic reduction mode for reproducible runs.
- SIMD: start scalar. Add portable SIMD after correctness.

## L2 device layer

- Rigid-body dynamics with position and orientation integration.
- Joints: revolute, prismatic, fixed, gear coupling.
- Parameters come from L1 by the extraction pipeline in `PARAMETERS.md`.
- Export a multibody model to URDF for external tools.

## Stack and dependencies

| Component | Choice | License | Why |
|---|---|---|---|
| Core language | Rust | Apache-2.0 (ours) | safety, speed, threads |
| Python bindings | PyO3 / maturin | Apache-2.0 / MIT | scripting and science interop |
| Standard MD and GPU | OpenMM | MIT | solved problem, GPU on Apple |
| QM and file interop | ASE | LGPL-2.1+ | one interface to many codes |
| Quantum | PySCF or xTB | BSD / LGPL-3.0 | L0 work |
| Chemistry and formats | RDKit | BSD | molecules and file formats |
| Linear algebra | nalgebra / faer | Apache-2.0 / MIT | dense and sparse |
| Task runner | just | MIT | simple commands |
| Python lint | ruff, black | MIT | fast, standard |
| Python build | maturin | MIT/Apache-2.0 | PyO3 packaging |
| Rust errors | thiserror | MIT/Apache-2.0 | library error types |
| Benchmarks | criterion | MIT/Apache-2.0 | statistics for timings |
| Serialization | serde, serde_json | MIT/Apache-2.0 | parameter and NCZ metadata |
| NCZ container | zip 7.2.0 | MIT | NCZ is a stored zip archive |
| PyO3 | pyo3 0.29 (optional, `python` feature) | MIT/Apache-2.0 | M7-01 bindings, abi3-py3.10 |

`zip` is pinned to 7.x for MSRV 1.83; zip 8 needs Rust 1.88. `pyo3` is an
optional dependency behind the `python` feature, so the default Rust gates do
not need a Python interpreter (ADR-0016). No linear-algebra dependency is used
yet: the jigs device code hand-rolls `[f64; 3]` and a quaternion.

Check every new dependency for its license. Record it here. Keep LGPL use
dynamic-link only.

## Interfaces

One headless API. Three clients call it.

Agent tools exposed through MCP (all typed):

- `create_part`, `generate_gear`, `generate_nanotube`, `generate_lattice`
- `relax`, `run_md`, `minimize`
- `add_jig`, `add_constraint`, `assemble`
- `measure`, `extract_parameters`
- `convert_units`, `save`, `load`, `list_parts`, `export`

The agent calls these tools. It never edits raw coordinates.

## Formats

- `NCZ`: native, versioned. A zip of JSON metadata and binary arrays. Little
  endian. A header carries the schema version.
- `MMP`: import and export for NanoEngineer interop.
- `XYZ`, `PDB`: through ASE.
- `URDF`: export for the device level.
- `CellML`, `SBML`: later, for L4.

## Compute and performance strategy

- Correct scalar code first, with a benchmark and reference outputs.
- Then flat arrays and structure-of-arrays.
- Then multi-thread the force loop.
- Then add SIMD.
- Then consider the GPU through OpenMM for standard systems.
- Numerical type: `f64` for accumulation and correctness. `f32` only behind a
  measured flag.

## Agent interface

- The Python package exposes the tools.
- An MCP server wraps the package.
- The agent gives an intent in words. The tool layer converts it to typed
  operations. The core validates every result.
- The agent must show its work and its uncertainty. Treat output as a proposal.

## Verification architecture

- Finite-difference gradient test per force term, in `cargo test`.
- Energy-conservation test for each integrator, NVE drift within a set bound.
- Reference benchmarks against OpenMM and GROMACS CPU.
- Golden regression files with hashes for import and export.
- Sanitizers for any C or FFI boundary.
- A public benchmark page with timing and hardware.
