# TASKS.md — The Iterative Backlog

This is the backlog. A coding agent picks the first unchecked task in the
current milestone. When the user names a task, do that one.

Work loop: mark the line with ` <!-- WIP -->`, implement with tests, run
`just verify`, then change `[ ]` to `[x]` and append a result note.

## Ready now

Start here on a clean session. `M0-01` is the first task.

---

## M0 — Skeleton

- [x] **M0-01** Initialize the Rust workspace and the crate list in
  `ARCHITECTURE.md`. Each crate builds and has one placeholder test.
  - Done when: `cargo build --workspace` and `cargo test --workspace` succeed.
  - Result: 9 crates (`units`, `model`, `format`, `engine`, `jigs`, `parts`,
    `params`, `python`, `cli`). `cargo test --workspace` runs 68 tests green.
- [x] **M0-02** Add CI (GitHub Actions): fmt, clippy, test on Linux and macOS.
  - Result: `.github/workflows/ci.yml` runs fmt, clippy `-D warnings`, and test
    on `ubuntu-latest` and `macos-latest`.
- [x] **M0-03** Add `justfile` with `verify`, `fmt`, `lint`, `test`, `bench`.
  - Done when: `just verify` runs every gate and returns zero.
  - Result: `just verify` = `fmt-check` + `lint` + `test`; also `bench` and
    `py-test`. Green.
- [x] **M0-04** Add `LICENSE` (canonical Apache-2.0 text) and a license-header
  policy. Never copy GPL code.
  - Result: canonical Apache-2.0 vendored in `LICENSE` (202 lines), `NOTICE`
    added, policy in `docs/license-header.md`.
- [x] **M0-05** Add `docs/glossary.md`: atom, bond, part, assembly, jig, jig
  types, L0 to L4, pair list, skin, NCZ, MMP.
  - Result: `docs/glossary.md` defines all required terms.
- [x] **M0-06** Build the `units` crate: SI base units, a `Quantity` type, and
  conversions for Angstrom, nanometre, kilocalorie per mole, femtosecond,
  picosecond, atomic mass unit, and electron charge. Round-trip tests.
  - Result: `Dimension`, `Unit`, `Quantity` in `crates/units`. 7 SI base units
    plus 7 non-SI units. 17 tests pass, including round trips and a
    dimension-mismatch error.
- [x] **M0-07** Add the benchmark harness (criterion) with one trivial
  benchmark, and a results directory.
  - Result: Criterion bench `crates/units/benches/units_bench.rs`,
    `benchmarks/README.md`, `benchmarks/results/`. Compiles clean.
- [x] **M0-08** Add the Python package skeleton (pyproject, maturin) with an
  empty module and one pytest.
  - Result: `python/pyproject.toml` (maturin), `python/nanocad/__init__.py`
    with typed `version()`, and `python/tests/test_smoke.py` (passes under
    pytest). pyo3 binding is deferred to M7-01.

## M1 — Document model and formats

- [x] **M1-01** `model`: `Atom`, `Bond`, `Topology` in struct-of-arrays form,
  with index-based access. Unit tests.
  - Result: `crates/model` stores elements, a flat three-per-atom position
    buffer, charges, and types separately. Index access returns `Option` or
    `Result`; no path panics. 43 tests pass.
- [x] **M1-02** `model`: `Part` and `Document`, versioned. Serialization to a
  byte buffer.
  - Result: `Part` and `Document` (`SCHEMA_VERSION = 1`) with a little-endian
    codec. Header carries the `NCAD` magic and the schema version. Decode
    rejects bad magic, a foreign version, truncation, and trailing bytes.
    Round-trip tests pass.
- [x] **M1-03** `format`: NCZ writer and reader (zip of JSON metadata and
  binary arrays), with a header schema version. Round-trip test.
  - Result: `write_ncz`/`read_ncz` in `crates/format/src/ncz.rs`. ZIP (Stored),
    `header.json` carries `schema_version = 1`. f64 coordinates, bit-exact.
- [x] **M1-04** `format`: XYZ import and export.
  - Result: `import_xyz`/`export_xyz`. Angstrom at the boundary. Approximate,
    not bit-exact.
- [x] **M1-05** `format`: MMP import. Parse atom lines, bond lines, atom types,
  and coordinates. Handle the NanoEngineer conventions.
  - Result: `import_mmp` in `crates/format/src/mmp.rs`. MMP units are 1e-13 m
    (0.001 Angstrom). Known gap: `bondg`/`bondc` import as `BondType::Unknown`.
- [x] **M1-06** `format`: MMP export.
  - Result: `export_mmp`, round-trips topology and coordinates.
- [x] **M1-07** Golden round-trip test: import a real MMP, save to NCZ, reload,
  and reproduce the coordinates bit for bit. Use the universal joint file from
  the old tree at `/tmp/nanoengineer/cad/partlib/couplings/`.
  - Result: `tests/golden/golden_roundtrip.rs` plus an original synthetic
    fixture. Real file test ran: 6828 atoms and 10452 bonds reproduced bit for
    bit through NCZ. The GPL source file was read, never copied.
- [x] **M1-08** `cli`: `nanocad info <file>` and `nanocad convert <in> <out>`.
  - Result: `nanocad` in `crates/cli` supports `info`, `convert`, `--help`,
    `--version`. Format from the extension. Bad input exits 1 or 2 with a clear
    stderr line; no panic. 13 tests pass.

## M2 — Forces and minimization

- [x] **M2-01** Neighbor search: grid cell list plus a Verlet skin, producing a
  flat `u32` pair list. Test against a brute-force O(n^2) reference.
  - Result: `VerletList` in `crates/engine/src/neighbor.rs`. Flat `u32` pairs,
    one per unordered pair. Rebuild when an atom moves more than half the skin.
    Brute-force reference test passes. 21 engine tests total.
- [x] **M2-02** Bond-stretch force term plus a finite-difference gradient test.
  - Result: `BondStretchTerm` in `crates/engine/src/bond_stretch.rs`. Central
    finite difference, step 1e-14 m. Measured max error 1.28e-17 N.
- [x] **M2-03** Angle-bend term plus a finite-difference gradient test.
  - Result: `angle_bend.rs`, `U = 0.5 k (theta - theta0)^2`. Max FD error
    8.3e-18 N.
- [x] **M2-04** Torsion term plus a finite-difference gradient test.
  - Result: `torsion.rs`, three-term cosine series. Max FD error 7.5e-18 N.
- [x] **M2-05** Out-of-plane term plus a finite-difference gradient test.
  - Result: `out_of_plane.rs`, harmonic in the signed angle. Max FD error
    9.9e-17 N.
- [x] **M2-06** van der Waals term (Buckingham or MM4 table) plus a
  finite-difference gradient test.
  - Result: `van_der_waals.rs`, Buckingham with a switching cutoff. Max FD
    error 7.0e-20 N.
- [x] **M2-07** Electrostatic term with a cutoff plus a finite-difference test.
  - Result: `electrostatic.rs`, Coulomb with a switch. Max FD error 1.7e-18 N.
  - Note: non-bonded terms scan all pairs; wiring them to the Verlet list is
    part of M2-08.
- [x] **M2-08** Assemble total potential and force; finite-difference test on a
  small molecule.
  - Result: `System` in `crates/engine/src/system.rs` owns every term and now
    feeds the non-bonded terms from the Verlet list. Max FD error 5.7e-17 N on a
    6-atom chain with six terms.
- [x] **M2-09** Conjugate-gradient minimizer with a convergence test.
  - Result: `minimize` in `crates/engine/src/minimize.rs`. Polak-Ribiere+ with a
    strong-Wolfe line search. Converged to a gradient of 5.9e-15 N in 53
    iterations.
- [x] **M2-10** NVE energy-conservation test with a drift bound.
  - Result: `VelocityVerlet` in `crates/engine/src/integrator.rs`. Relative
    drift 1.278e-5 over 4000 steps, below the 1e-4 bound.
- [x] **M2-11** Reference water-box benchmark compared to OpenMM within
  tolerance for energy and forces.
  - Result: PASS against OpenMM 8.6.1 on the Reference platform. 200 molecules,
    600 atoms, seed 0x5A17E2, cubic box 1.8833259369189425e-9 m. Energy
    relative difference 2.341e-14, maximum per-atom force relative difference
    1.242e-14, tolerance 1.0e-6. `crates/engine/examples/water_box.rs`,
    `benchmarks/openmm/compare_water_box.py`, and
    `benchmarks/results/openmm-water-crosscheck.txt`. The CPU platform agrees
    to about 4e-6 (single precision).
  - Note: a coding check, not a physics validation. Both codes share the cutoff
    and the missing long-range correction. The model is SPC-like, not SPC.
- [x] **M2-12** Record timings in `benchmarks/results/`.
  - Result: `crates/engine/benches/engine_bench.rs` plus
    `benchmarks/results/2026-09-14-engine-timings.txt`. Apple M4 release: energy
    680.97 us, serial forces 679.83 us, 4-worker 370.97 us.

## M3 — Dynamics and jigs

- [x] **M3-01** Velocity Verlet integrator plus energy conservation.
  - Result: `VelocityVerlet` in `crates/engine/src/integrator.rs`. See M2-10 for
    the drift measurement. Masses live on the `System`.
- [x] **M3-02** Langevin thermostat plus a temperature-distribution test.
  - Result: `LangevinThermostat` in `crates/engine/src/langevin.rs`, with a
    seeded in-crate `Rng` (`rng.rs`). Mean 300.58 K at a 300 K target.
- [x] **M3-03** Berendsen thermostat.
  - Result: `BerendsenThermostat` in `crates/engine/src/berendsen.rs`. Relaxed
    10.23 K to 300.00 K. Documented as a relaxation thermostat, not a correct
    canonical sampler.
- [x] **M3-04** Anchor jig.
  - Result: `AnchorJig` in `crates/jigs/src/anchor.rs`. Max FD error 4.8e-22 N.
- [x] **M3-05** Motor jig: constant angular velocity or torque.
  - Result: `MotorJig` in `crates/jigs/src/motor.rs`. A torque source, not a
    potential. Tests check zero net force and exact axial torque.
- [x] **M3-06** Spring jig.
  - Result: `SpringJig` in `crates/jigs/src/spring.rs`. Atom-atom and
    atom-to-fixed-point. Max FD error 1.0e-16 N.
- [x] **M3-07** Test machine: a rotor driven by a motor at 300 K runs stable
  for many frames.
  - Result: `crates/jigs/src/rotor.rs`. 200000 frames: mean temperature
    309.57 K, max bond strain 1.965e-2, max force residual 8.564e-26 N.
- [x] **M3-08** Deterministic multi-threaded force reduction with a
  reproducibility test.
  - Result: `System::forces_n_parallel` and `gradient_j_per_m_parallel` using
    `std::thread::scope`, no new dependency. 1-thread and N-thread results are
    bit-identical.

## M4 — Part generation

- [x] **M4-01** Define the part-generator trait and the part schema.
  - Result: `PartGenerator` trait, `PartSchema`, `ParameterSpec`, `ParameterSet`
    in `crates/parts`. Inputs are resolved, defaulted, and range-checked.
- [x] **M4-02** Lattice builder for diamond and graphite.
  - Result: `DiamondGenerator` and `GraphiteGenerator` in
    `crates/parts/src/lattice.rs`. Diamond a = 3.567e-10 m, graphite a =
    2.461e-10 m, layers 3.354e-10 m. Bond-length and coordination tests pass.
- [x] **M4-03** Nanotube builder from chiral indices.
  - Result: `NanotubeGenerator` in `crates/parts/src/nanotube.rs`. Chiral
    vector roll-up, analytic radius, graphene bond 1.4209e-10 m. `(5,5)` gives
    20 atoms per cell, `(6,0)` gives 24.
- [x] **M4-04** Gear-tooth profile generator (a diamondoid analog).
  - Result: `GearProfile` in `crates/parts/src/gear_profile.rs`. Involute
    flanks, ISO full-depth form, undercut limit. Pitch diameter = module x teeth.
- [x] **M4-05** Spur-gear generator: teeth, thickness, bore, hub.
  - Result: `SpurGearGenerator` in `crates/parts/src/spur_gear.rs`. Rim,
    spokes, hub, empty bore. Skeletal, as the upstream gear moieties are.
- [x] **M4-06** Planetary-set generator: sun, planets, ring, carrier.
  - Result: `PlanetaryDesign` and `PlanetaryGenerator` in
    `crates/parts/src/planetary.rs`. Ring teeth = sun + 2*planet. Assembly and
    spacing constraints checked. Default sun 24, ring 60 gives ratio 3.5.
- [x] **M4-07** Valence and strain validation checker.
  - Result: `validate_part` in `crates/parts/src/validation.rs` returns every
    `Violation` in one pass: `OverCoordinated`, `BondStrain`, `UnknownValence`.
- [x] **M4-08** Agent-facing `generate_gear(name, specs)`, exposed in the
  Python API.
  - Result: `generate_gear(name, specs)` in `crates/parts/src/registry.rs`, plus
    `gear_generator` and `gear_generators`. Python exposure stays open under
    M7-01.

## M5 — Parameter extraction (L1 to L2)

- [x] **M5-01** `params`: the `PartRecord` type and its serialization.
  - Result: `PartRecord`, `Provenance`, `Method`, `Validation` in
    `crates/params`. Versioned JSON, header `schema`/`version`. Null means an
    unknown uncertainty. 38 tests pass.
- [x] **M5-02** Stiffness extraction: strain sweep, fit, uncertainty.
  - Result: `extract_stiffness` in `crates/params/src/stiffness.rs`. Synthetic
    chain k = 300 N/m gave 4.49999999999e12 Pa against an analytic 4.5e12 Pa
    (relative error 2.32e-12). Noisy run agrees within 3 sigma.
- [x] **M5-03** Friction extraction: sliding test, mean and spread.
  - Result: `extract_friction` and `summarize_friction` in
    `crates/params/src/friction.rs`. Sliding contact mu = 0.4433 with spread
    3.37e-12 N. A simulated estimate, not a material constant.
- [x] **M5-04** Failure-stress extraction.
  - Result: `extract_failure_stress` in `crates/params/src/failure.rs`. Rupture
    at strain 0.05 gave 2.25e11 Pa against an analytic 2.25e11 Pa.
- [x] **M5-05** Thermal-property extraction under NVT.
  - Result: `crates/params/src/thermal.rs`. Specific heat from kinetic-energy
    fluctuations (40-atom ideal gas: 1.0414e3 J/(kg*K) vs analytic 1.0393e3,
    relative error 2.0e-3). Conductivity from a two-bath direct estimate
    (harmonic chain: 8.34 +/- 0.93 W/(m*K)). Both `Method::Md`,
    `Validation::Unverified`. 6 tests.
- [x] **M5-06** Versioned parameter library with provenance.
  - Result: `ParameterLibrary` in `crates/params/src/library.rs`. Keyed by
    `part_id`, several revisions each, geometry-hash check, JSON persistence,
    revision counter.
- [x] **M5-07** Consistency-check tool for a parameter record.
  - Result: `check(&PartRecord) -> Vec<Problem>` in `crates/params/src/check.rs`.
    Checks the version, required fields, dimensions, nonnegative uncertainty,
    and physical plausibility. Returns all problems in one pass.

## M6 — Device level (L2)

- [x] **M6-01** Rigid-body representation and integrator.
  - Result: `RigidBody` and `RigidBodySystem` in `crates/jigs/src/device.rs`.
    Semi-implicit Euler, world angular momentum as the primary state. No
    nalgebra added.
- [x] **M6-02** Revolute joint.
  - Result: `RevoluteJoint` in `crates/jigs/src/device.rs`. Max anchor error
    0 m, axis error 9.9956e-13 rad over 2000 steps.
- [x] **M6-03** Prismatic joint.
  - Result: `PrismaticJoint` in `crates/jigs/src/device.rs`. Max error 0 m,
    orientation 5.42e-20 rad over 2000 steps.
- [x] **M6-04** Gear-coupling constraint.
  - Result: `GearCoupling` in `crates/jigs/src/device.rs`. Max error 5.551e-17 m.
    Ratio -0.500000000 for r_a = 0.05 m and r_b = 0.10 m.
  - Note: the constraint pass is position-level only. A velocity-level pass
    stays open.
- [x] **M6-05** Assemble the planetary gearbox from `PartRecord`s.
  - Result: `assemble_planetary`, `assemble_from_records`, `planetary_records`
    in `crates/jigs/src/assembly.rs`. Builds the device from the parts
    generators and the parameter library.
- [x] **M6-06** Gear-ratio test within tolerance.
  - Result: measured ratio 3.500000000, analytic 3.500000000, relative error
    1.146e-10, tolerance 1.0e-3.
- [x] **M6-07** URDF export.
  - Result: `export_urdf`, `parse_urdf`, `format_urdf` in
    `crates/jigs/src/urdf.rs`. Emitted XML round-trips bodies, joints, and axes.

## M7 — Agent interface

- [x] **M7-01** Python API surface for every tool in `ARCHITECTURE.md`.
  - Result: pyo3 bindings in `crates/python/src/bindings/**` behind an optional
    `python` feature (ADR-0012 keeps the default build a plain rlib). The
    extension `nanocad._core` exposes units, model, formats, parts, params,
    engine, and jigs. Built with maturin. `python/tests` passes 20 tests against
    a fresh abi3 build. `extract_parameters` stays out until M5-05 lands.
- [x] **M7-02** MCP server wrapping the API.
  - Result: `mcp/nanocad_mcp/` (stdlib JSON-RPC 2.0 over stdio, no third-party
    dependency). 18 tools. Subprocess test passes `initialize`, `tools/list`,
    and tool calls.
- [x] **M7-03** Tool schemas and usage docs.
  - Result: `mcp/tools.json` plus `mcp/README.md` with a worked example.
- [x] **M7-04** Scripted end-to-end agent demo test: build the gearbox through
  tools alone.
  - Result: `mcp/tests/test_agent_demo.py`. Builds the gearbox through
    `call_tool` only. Measured gear ratio 3.500000000, relative error
    1.146e-10, tolerance 1.0e-3. URDF 7 links, 6 joints.

## M8 — Demo and release

- [x] **M8-01** Three-scale export for the visualization.
  - Result: `crates/jigs/src/scene.rs`. `build_scene(assembly, set)` produces
    atomistic, device, and coarse layers. `scene_to_json`,
    `scene_to_json_pretty`, `scene_from_json`, `write_scene_json`. Schema
    `nanocad.scene` version 1. Default design: 7 bodies, 6 joints, 3 planets.
    See `docs/three-scale.md`.
- [ ] **M8-02** The 90-second gearbox video and its script.
  - Result: the script is written at `docs/video-script.md` (10 shots, 90 s
    budget, commands, asset list). The video is not recorded. This stays open
    until the recording exists.
- [x] **M8-03** Written report with the measured numbers and the caveats.
  - Result: `docs/report.md`. Every number cites its source. Caveats and
    limitations section. Marked simulated, no medical claim.
- [x] **M8-04** Public benchmark page with hardware and timing.
  - Result: `docs/benchmarks.md` plus `scripts/bench.sh`. Apple M4, cargo
    1.98.1: energy 680.97 us, serial forces 679.83 us, 4-worker 370.97 us.
    Every number is machine-specific.
- [x] **M8-05** Clean-clone reproduction test for the whole demo.
  - Result: `scripts/reproduce.sh` plus `just reproduce`. Creates or reuses a
    venv, builds the extension from `python/`, runs the Rust suite and the MCP
    demo test, prints the measured gear ratio. First run needs network access.
    See `docs/reproduce.md`. Guard test in
    `mcp/tests/test_reproduce_script.py`.

## Phase 6 — Later

- [ ] **L0-01** Quantum adapter through ASE and PySCF or xTB.
- [ ] **L3-01** Continuum adapter for structural and flow.
- [ ] **L4-01** Lumped system model with CellML or SBML interop.
- [ ] **FL-01** Respirocyte rotor and bearing subsystem.
- [ ] **FL-02** Respirocyte pump and gas tank subsystem.
- [ ] **NM-01** One nanomedicine slice: a machine coupled to a coarse host.
