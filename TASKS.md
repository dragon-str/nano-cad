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
- [x] **M8-02** The 90-second gearbox video and its script.
  - Result: the generated artifact is `docs/media/nano-cad-gearbox.mp4` with
    the sidecar `docs/media/nano-cad-gearbox.srt`. H.264, 1920x1080, 30 fps,
    90.000 s, AAC audio, burned-in captions. `scripts/make_video.sh` regenerates
    both from scratch with `/opt/homebrew/bin/python3` and ffmpeg;
    `scripts/render_video.py` draws every frame. The narration is synthetic
    text-to-speech (macOS voice Samantha), and the visuals are schematic and
    labelled simulated. Every on-screen number cites its source file.
  - Fix: the gear scenes now use the fixed-ring kinematics. The carrier rate is
    `w_c = w_s * 2/7` and the planet absolute spin is `w_p = -(2/3) * w_s`, so
    the sun and the planets turn in opposite directions. The planet tooth
    phase keeps the sun, planet, and ring teeth meshed. See `docs/video.md`.
    `docs/video-script.md` holds the storyboard and the real repository URL.
    See `docs/video.md`. The video is not part of `just verify`.
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

- [x] **L0-01** Quantum adapter through ASE and PySCF or xTB.
  - Result: `python/nanocad/quantum_adapter.py`. PySCF 2.14.0 from
    `/tmp/nc-qm-venv` ran RHF/STO-3G on H2 and gave -1.116999 Ha against the
    cited Hehre-Stewart-Pople 1969 value -1.1167 Ha. ASE EMT is a labelled
    classical fallback. Backends import lazily; the test skips when absent.
  - Note: single-point energy only. No forces, so no finite-difference test.
    xTB path is written but not exercised.
- [x] **L3-01** Continuum adapter for structural and flow.
  - Result: `python/nanocad/continuum_adapter.py`. Axial bar and
    Euler-Bernoulli beam, plus Hagen-Poiseuille tube and plane-channel flow.
    Relative error 1.4e-16 (bar), 2.4e-12 (beam), 5.8e-11 (tube), 1.5e-11
    (channel). `voxelize_part` maps a plain atom list to density.
- [x] **L4-01** Lumped system model with CellML or SBML interop.
  - Result: `python/nanocad/lumped_adapter.py`. A two-compartment respirocyte
    gas-transport model, fixed-step RK4, and an SBML Level 3 subset round trip.
    Analytic check reaches 6.55e-12 relative error.
- [x] **FL-01** Respirocyte rotor and bearing subsystem.
  - Result: `RespirocyteRotorGenerator` in `crates/parts/src/respirocyte.rs`.
    Rotor and bearing from the diamond lattice, a stated radial gap, and a
    rotational degree of freedom. No bonds cross the gap.
- [x] **FL-02** Respirocyte pump and gas tank subsystem.
  - Result: `RespirocytePumpGenerator` and `RespirocyteTankGenerator` in
    `crates/parts/src/respirocyte.rs`. Cylinder and piston, plus a hollow
    spherical shell. The tank interior is empty.
- [x] **NM-01** One nanomedicine slice: a machine coupled to a coarse host.
  - Result: `crates/jigs/src/nanomedicine.rs`. `HostEnvironment` (Stokes drag,
    optional flow, tether, drive) and `NanoMachine`. Terminal velocity matches
    the analytic value to 2.0e-9 relative; the energy balance closes to 3.4e-7.
    All geometry is a skeletal lattice model, not a validated device.

## Phase 7 — Depth, interop, device and media

The depth pass after M8. Each task below deepens a layer that was a stub in
Phase 6. All result notes cite the test that measures the number.

- [x] **ENG-01** Add an L-BFGS minimizer stage next to the conjugate gradient.
  - Result: `crates/engine/src/minimize.rs` adds `MinimizeMethod`
    (`ConjugateGradient`, `Lbfgs`, `ConjugateGradientThenLbfgs`) and
    `minimize_with`. `minimize` stays conjugate-gradient only. On a stiff chain
    of 40 iterations, the gradient reaches 2.5945e-7 N with CG and 1.00997e-7 N
    with L-BFGS. The hybrid method converges in 53 iterations to 5.9049e-15 N.
- [x] **ENG-02** Make the parallel force reduction deterministic and faster.
  - Result: `crates/engine/src/system.rs` partitions the pair list into
    canonical blocks and sums the same per-block partials in serial and in
    parallel. `forces_n_parallel` is bit-identical to `forces_n` for one, two,
    four and eight workers (400 atoms, 12414 pairs, 7 blocks). Serial path
    changed to sum the same partials. Serial 423.00 us, 8 workers 214.18 us.
- [x] **ENG-03** Add a portable nonbonded kernel with a scalar fallback.
  - Result: `crates/engine/src/pair_kernel.rs` holds a four-lane unrolled
    Buckingham plus electrostatic kernel and a scalar fallback. Both paths are
    bit-identical (maximum error 0.0). `std::simd` is deferred because it needs
    nightly; the MSRV stays 1.85 and the rustc is 1.98.1.
- [x] **IO-01** Import and export PDB.
  - Result: `crates/format/src/pdb.rs` adds `import_pdb` and `export_pdb`.
    ATOM, HETATM, element, CONECT bonds, CRYST1 box. The boundary is Angstrom;
    PDB carries no charge field, so charge is lost. Multi-part documents do not
    round-trip. New `FormatError::PdbParse`; the CLI accepts `FileFormat::Pdb`.
    8 tests in the module.
- [x] **IO-02** Add an optional RDKit molecule adapter.
  - Result: `python/nanocad/rdkit_adapter.py` converts to and from an RDKit
    molecule and reads and writes SMILES and SDF. RDKit 2026.03.6 runs from
    `/tmp/nc-qm-venv`. 10 tests; they skip when RDKit is absent.
- [x] **IO-03** Export the lumped model as CellML.
  - Result: `python/nanocad/lumped_adapter.py` adds `to_cellml` and
    `from_cellml`. The subset is honest and stated. It has no reaction element.
    The units attribute is dimensionless and a `nanocad:unit` attribute carries
    the true unit. 19 tests in `python/tests/test_lumped.py`.
- [x] **DEV-01** Add a fixed joint and a velocity-level constraint pass.
  - Result: `crates/jigs/src/device.rs` adds `FixedJoint` and
    `ConstraintOptions { velocity_pass }` (on by default). The velocity pass
    projects the constraint velocity to zero for revolute, prismatic, gear and
    fixed joints. Measured revolute maximum constraint velocity 1.118 to 0;
    gear 0.1 to 2.78e-17; weld anchor 0 m and orientation 2.23e-14 rad. The
    Python binding exposes the new methods.
- [x] **L3-02** Solve two L3 problems with an external library.
  - Result: `python/nanocad/continuum_adapter.py` adds
    `solve_plane_stress_cantilever` (Q4 finite elements through scipy) and
    `solve_poiseuille_2d` (finite differences). The cantilever tip is
    1.5893e-5 m against a refined 1.5986e-5 m, relative 5.78e-3. The flow
    maximum velocity matches the analytic value to 2.44e-4 relative.
- [x] **ENG-04** Compute the electrostatic energy with a particle-mesh method.
  - Result: `python/nanocad/pme_adapter.py` adds `electrostatic_energy_pme`
    through the OpenMM PME reference. A Na-Cl box of 3 nm gives -3.3587e-19 J
    against a cut-off -3.2627e-19 J, relative -2.94%. At 10 nm the relative
    difference is 2.42e-4. OpenMM is a cross-check, not a second engine.
- [x] **FL-03** Add the respirocyte seal clearance and leakage model.
  - Result: `crates/parts/src/respirocyte.rs` adds `AnnularGapFlow`. The
    pressure-driven flow is `Q_p = pi R h^3 / (6 mu) * dp/L`; the Couette flow
    is `Q_c = pi R h U` (Bird, Stewart and Lightfoot). The bearing leakage is
    9.114807e-20 m^3/s and matches the hand value. The seal leakage is
    1.675516e-21 m^3/s. New parameters `seal_clearance_m` and `seal_land_m`.
- [x] **NM-02** Add Brownian motion to the host environment.
  - Result: `crates/jigs/src/nanomedicine.rs` adds
    `HostEnvironment::with_temperature_kbt_j`, `temperature_k`,
    `NanoMachine::with_seed`, and a Brownian force. The velocity variance is
    3.9988e-3 against the kBT/m value 4.1419e-3, relative 3.46e-2. The zero
    temperature terminal velocity matches to 2.04e-9 relative.
- [x] **VAL-01** Cross-check the engine against GROMACS.
  - Result: `benchmarks/gromacs/` holds `compare_water_box.py`,
    `selftest_switch.py`, `README.md` and `results/`. GROMACS 2026.3 does not
    use the Buckingham potential with the Verlet cutoff scheme. The driver
    records SKIP, writes no number, and exits 0. `docs/validation.md` states
    this limitation. The Coulomb constant check matches to 6.75e-15 relative.
- [x] **M9-01** Build a static viewer and a documentation site.
  - Result: `crates/jigs/examples/scene_json.rs` writes `site/scene.json` and
    `site/scene.data.js`: 7 bodies, 6 joints, 3 planets, ratio 3.5, 2202 atoms.
    `site/index.html`, `site/viewer.js` and `site/style.css` draw three layers
    on a canvas and work from `file://`. `scripts/build_docs_site.py` writes 17
    pages to `site/docs`. `site/check.py` checks the site. See `docs/viewer.md`.
- [x] **M9-02** Make the atom layer of the video a real bonded lattice.
  - Result: `scene_json.rs` also writes `site/scene.bonds.json`.
    `scripts/check_atom_geometry.py` is the gate. `scripts/render_video.py`
    draws atoms and bonds in shots 1 and 8. `scripts/make_video.sh` runs the
    gate. Measured 216 atoms and 333 bonds: mean C-C distance 1.544556e-10 m
    against 1.544e-10 m (relative 3.603e-4) and mean bond angle 109.471221
    degrees against 109.4712 degrees (error 2.06e-5 degrees). 64 interior
    carbons all have 4 bonds.
  - Limitation: `PlanetaryGenerator` returns a skeletal gear outline, not a
    diamondoid solid. The gear shape is therefore not yet chemically dense. The
    `DiamondGenerator` block in `crates/parts/src/lattice.rs` is the verified
    material basis. See ADR-0038. The visuals were corrected in M9-03.
- [x] **M9-03** Make the schematic match the atomic layer and draw the gears as
  atoms.
  - Result: `scripts/gear_profile.py` is a port of
    `crates/parts/src/gear_profile.rs`. The renderer and the checker import it,
    so one formula drives both. The schematic now uses the exact full-depth
    involute profile with addendum `1.0 m` and dedendum `1.25 m`, and the ring
    is drawn as an internal gear with inward teeth. `scripts/render_video.py`
    adds a `GearLattice` class that draws the real `PlanetaryGenerator` atoms
    and bonds from `site/scene.json` and `site/scene.bonds.json`, animated with
    the fixed-ring rates.
  - Verification: `scripts/check_atom_geometry.py` now checks the gear layer.
    Measured tip and root radii equal the involute profile exactly: sun
    5.375e-9 to 6.500e-9 m, planet 3.875e-9 to 5.000e-9 m about its center,
    ring internal 1.450e-8 to 1.5625e-8 m. It also checks that the shared
    schematic formula gives the same extremes. The diamond-block checks stay.
    The video is regenerated: 1920x1080, 30 fps, 90.000 s. At that time the
    gear layer was a single 2D layer of 2202 atoms and 2205 bonds; M9-04
    extrudes it to four axial layers. See ADR-0039.
- [x] **M9-04** Give the gears an axial thickness and make the atom layer mesh
  with the schematic.
  - Result: `crates/parts/src/planetary.rs` adds `layers` and
    `layer_spacing_m` parameters and an `add_extruded_loop` helper. The sun,
    the planets, and the ring are extruded into centered axial layers with
    vertical bonds; the part metadata carries `layers`, `layer_spacing_m`, and
    `thickness_m`. `crates/jigs/src/scene.rs` exports the three values in
    `SceneDesign`. The planet rotation now carries the half-tooth offset
    `pi/N_p`, so a sun tooth enters a planet space and a planet tooth enters a
    ring space. Without it the generator placed tooth on tooth and the gears
    collided. The renderer now uses one shared rate function: the planet spins
    about its center at the relative rate `-(N_s/N_p)*(w_s - w_c)` and the
    carrier orbit supplies the rest.
  - Verification: `cargo test -p nanocad-parts planetary` -> 13 pass,
    including `the_planet_presents_a_tooth_space_at_each_mesh_line` and
    `the_gear_thickness_is_the_axial_layer_span`.
    `scripts/check_atom_geometry.py` checks the mesh phase at every planet, the
    layer count and the thickness, and that the schematic and the atoms share
    one rate. Default scene: 8340 atoms and 14481 bonds, 4 layers x 1.544e-10 m
    = 4.632e-10 m. The video is regenerated at 1920x1080, 30 fps, 90.000 s.
    See ADR-0040.
- [x] **M9-05** Add a live parameter panel and a chat interface.
  - Result: `app/chat.py` is a rule-based parser for short edit commands
    ("one atomic layer thicker", "6 atoms thick", "add more teeth", "4
    planets", "move the layers 20 pm apart", "reset"). It clamps every value to
    a stated limit and reports the new state. `app/server.py` is a standard
    library web server that runs the Rust generator for each change and serves
    the viewer. `site/index.html`, `site/viewer.js`, and `site/style.css` add
    the numeric panel, the chat box, and the live scene update. A change the
    engine rejects keeps the old parameters and reports the engine error.
  - Verification: `app/tests/test_chat_parser.py` -> 12 pass. The server was
    smoke-tested on `127.0.0.1:8123`: `/api/meta`, `/api/build?layers=6`, and a
    `/api/chat` message all returned 200 and a new scene. See ADR-0041.
- [x] **M9-06** Fill the gears with solid hydrogen-capped diamond and switch the
  viewer between atoms and the schematic by zoom.
  - Result: `crates/parts/src/diamond_solid.rs` fills each gear with a real
    diamond lattice cut to the involute profile, bonds it at the C-C length, and
    caps every under-coordinated carbon with a hydrogen at the C-H length. One
    atomic layer is one diamond (001) plane, so the plane spacing is fixed at
    `a/4` and `layer_spacing_m` is no longer a parameter. `site/viewer.js` shows
    depth-shaded atoms when the projected C-C spacing is six pixels or more and
    the exact involute schematic below that; both share one profile function and
    one world fit. `app/chat.py` no longer edits the spacing and explains that
    the crystal fixes it.
  - Verification: `cargo test -p nanocad-parts` -> 88 pass, including
    `every_solid_gear_atom_has_the_exact_diamond_valence`. `just verify` -> all
    gates passed. `/opt/homebrew/bin/python3 scripts/check_atom_geometry.py`
    -> PASS: 46618 C-C and 37284 C-H bonds, 0 outside tolerance, 0 carbon not
    4-bonded, 0 hydrogen not 1-bonded, mean C-C-C angle 109.471221 deg. The app
    smoke test on `127.0.0.1` returned 200 for `/api/meta`, `/api/build`, and
    `/api/chat`. See ADR-0043 and ADR-0044.

- [x] **M9-07** Give the atom layer a lit WebGL render, remove the scale slider,
  add panning, and add the viewer features that a solid-gear viewer needs.
  - Result: `site/render_atoms.js` is a WebGL2 renderer with instanced sphere
    impostors, key/fill/rim lighting, specular, a contact-shadow bias, a
    half-resolution screen-space ambient-occlusion pass, and a vignette. The
    viewer adds pan (shift-drag, middle-drag, shift-wheel), preset views, a
    turntable, motion playback of the planetary kinematics, fullscreen, PNG
    export, an element legend with per-element visibility, an atom-size control
    (space-filling to ball), a z-section clip, a hover pick, a two-click
    distance measure, a nanometre scale bar, and an orientation triad. The scale
    slider is gone; the layer checkboxes remain. When WebGL2 is absent the flat
    2D dots are the fallback. See ADR-0045.
  - Verification: `node --check site/render_atoms.js` and
    `node --check site/viewer.js` pass. `node site/render_atoms.test.js` -> 11
    passed. `/opt/homebrew/bin/python3 site/check.py` -> all checks passed.

- [x] **M9-08** Show only the atoms or the schematic, outline the atoms, and
  let the user hide the panel.
  - Result: The **Atoms** layer is on by default; the **Device markers** and
    **Coarse** layers are off, and the device body marker no longer draws a
    translucent sphere over the atoms. Each atom sphere gets a dark silhouette
    outline in the WebGL shader, and the flat fallback dots get a dark stroke.
    A **Hide** button collapses the panel. See ADR-0045.
  - Verification: `node --check` on both viewer files passes,
    `node site/render_atoms.test.js` -> 11 passed, and
    `/opt/homebrew/bin/python3 site/check.py` -> all checks passed.

- [x] **M9-09** Make every viewer control work in a real browser and soften the
  atom outlines.
  - Result: Two defects stopped the WebGL path. The AO, blur and composite
    programs were built with the vertex and fragment shaders swapped, so the
    renderer failed to build and the viewer silently used the flat fallback. The
    sphere impostor quads were then back-face culled because the NDC y-flip
    reverses the winding. With both fixed the lit spheres, the SSAO pass, the
    atom-size control, the clip and the legend behave. The gear spin now also
    animates the 2D schematic, which is what shows at the default zoom where the
    atom spacing is below six pixels. The outlines are softer: the WebGL edge
    mixes to 45 percent and the flat dots use a 0.45-alpha stroke.
  - Verification: `/tmp/nc-browser/verify.js` against `app/server.py` in headless
    Chrome (swiftshader) -> play advances, turntable, presets, atom size, clip,
    ambient occlusion, legend hide, panel hide, hover pick and measure all pass.
    `node --check` on both files, `node site/render_atoms.test.js` -> 11 passed,
    `site/check.py` -> all checks passed, `just verify` -> all gates passed.

- [x] **M9-10** Fix the dark diagonal band in the atom render.
  - Result: The band was the ambient-occlusion pass. The fragment shader
    projected each occlusion sample back to the screen with an inverted
    factor, so all sixteen samples read the centre depth and a large
    spurious region turned black. The projection now inverts the atom
    camera correctly. The result is a soft contact shadow of about
    25 percent, not a black band.
  - Verification: A raw-pixel probe (read the GL canvas into a 2D canvas)
    shows no dark rows with ambient occlusion on or off; the on/off
    luminance differs by about 25 percent. Screenshots at device pixel
    ratios 1 and 2 look the same. `node --check`, `node
    site/render_atoms.test.js` -> 11 passed, and the headless-browser
    feature check -> all features pass.

- [x] **M9-11** Fix the black mask over the top-right of the atom view.
  - Result: The final composite, the ambient-occlusion pass and the blur pass
    all draw one full-screen quad. The quad helper drew three vertices as
    triangles, so each pass covered one half of the screen. The final image
    showed atoms in only one half-plane and left the other half transparent:
    the reported black mask. The quad now draws four vertices as a triangle
    strip.
  - Verification: The composite, the direct atom pass and the composite again
    give the same pixel grid at zoom 4. A vision check reports no black mask.
    `node --check`, `node site/render_atoms.test.js` -> 11 passed, and the
    headless-browser feature check -> all features pass.

- [x] **M9-12** Remove the carrier atom plate.
  - Result: The carrier is a scene body only. The part carries no
    carrier atoms, the `carrier_offset_m` parameter is gone, and the
    generator has eight parameters. The assembly synthesises the
    carrier mass properties from a carbon ring at the carrier radius.
  - Verification: The scene has no atom at the carrier offset.
    `cargo test -p nanocad-parts` -> 88 passed. `site/check.py` ->
    all checks passed.

- [x] **M9-13** Give the gear teeth clearance so no atoms overlap.
  - Result: The addendum coefficient is 0.5, the dedendum coefficient is
    1.7, and the backlash is 9.0e-10 m. `site/viewer.js` mirrors the
    three constants.
  - Verification: An all-atom sweep over one carrier revolution gives
    sun-planet 3.0561 A, planet-ring 2.9265 A and planet-ring 2.8114 A.
    The minimum is above the 2.52 A target. 56160 atoms, 25928 carbon
    at degree 4 and 30232 hydrogen at degree 1. `check_atom_geometry.py`
    -> PASS. `just verify` -> all gates passed.

- [x] **M9-14** Add the metric crate and the clearance metric.
  - Result: New crate `nanocad-meter` with `Fidelity`, `MetricValue`,
    `Score` and a `Clearance` metric over the relative motion of the
    planetary bodies. `crates/meter/examples/planetary_score.rs` prints
    the score of the default set.
  - Verification: `cargo run -p nanocad-meter --example planetary_score`
    -> least distance 2.7753e-10 m between planet_1 and ring at 1.3475 s,
    which passes the 2.52e-10 m target and agrees with the off-line sweep.
    6 unit tests pass. `just verify` -> all gates passed. See ADR-0047.

- [x] **M9-15** Show a scorecard in the app.
  - Result: Geometric metrics live in `nanocad-meter`: atom count, contact
    ratio and clearance. The server route `GET /api/score` returns them as
    JSON. The app shows one row per metric with its value, its unit and its
    fidelity badge. `crates/meter/examples/score_json.rs` writes the JSON.
  - Verification: the browser shows atom count 56,160 atoms, contact ratio
    0.860 and clearance 0.278 nm with a pass mark. `just verify` -> all
    gates passed; 477 Rust tests. See ADR-0048.

- [x] **M9-16** Rebuild the default set at the recommended nanoscale size.
  - Result: The default is now module 1.5 nm, 12 sun teeth, 9 planet teeth,
    30 ring teeth, a 30 degree pressure angle, addendum coefficient 0.8,
    dedendum coefficient 1.25 and backlash 1.0e-9 m. The mesh phase follows
    the tooth-count parity: a planet presents a space toward the sun when
    its tooth count is odd, and the ring takes a half-pitch rotation when
    the planet tooth count is odd. The scene has 142091 atoms, 68197 carbon
    and 73894 hydrogen.
  - Verification: `cargo run -p nanocad-meter --example score_json` ->
    atom count 142091, contact ratio 1.0028, clearance 2.9968e-10 m
    (passes). `check_atom_geometry.py` -> PASS. The renderer draws 142091
    atoms in the browser. `just verify` -> all gates passed; 477 Rust
    tests. See ADR-0049.

- [x] **M9-17** Add the quasi-static slip-barrier metric.
  - Result: New metric `slip_barrier` in `crates/meter/src/slip.rs`. It
    sweeps one sun tooth pitch with the real gear kinematics and sums a
    shifted Lennard-Jones interaction between the sun and the first
    planet. The barrier is the largest energy minus the smallest. The
    verdict reads the ratio to the thermal energy `k T`.
  - Verification: `cargo run -p nanocad-meter --example score_json` ->
    slip barrier 5.7328e-20 J, which is 13.8 kT at 300 K ("the slip is
    marginal") with 27531 pairs at the peak. 9 unit tests pass. See
    ADR-0050.

- [x] **M9-18** Make the default drag rotate the scene with the pointer.
  - Result: The default drag rotates, so a point follows the pointer: a drag
    right turns the scene right and a drag up turns it up. Pan moves to
    Shift-drag and middle-drag. The hint text and `docs/viewer.md` follow.
  - Verification: A headless drag right gives yaw -0.650 -> -1.650 with panX
    unchanged, a drag up gives tilt 1.020 -> 1.820 with panY unchanged, and
    a Shift-drag right gives panX 0 -> 100 with yaw unchanged. `node --check
    site/viewer.js` passes, the JS suite is 11 passed, and `site/check.py`
    -> all checks passed.

- [x] **M9-19** Add a relaxed slip barrier that lowers the rigid bound.
  - Result: New metric `relaxed_slip_barrier` in
    `crates/meter/src/relaxed.rs`. It selects the atoms that face the mesh,
    frees the planet contact atoms, and relaxes them against the sun on a
    300 N/m tether. The energy is the interaction energy plus the tether
    energy, so the rigid path is a feasible point. The report also carries
    `rigid_barrier_j` on the same samples and atom sets.
  - Fix: the sample count was too low. One sun tooth pitch spans about
    fifty carbon bonds, so 24 samples aliased the barrier by 70 percent.
    The default is now 240 samples for the rigid metric and 120 for the
    relaxed metric. The rigid value at 60 and 240 samples agrees.
  - Verification: `cargo run -p nanocad-meter --example score_json` ->
    relaxed barrier 5.7371e-20 J against the rigid 5.7328e-20 J, so
    relaxation changes the bound by +0.07 percent and both read 13.8 and
    13.9 kT at 300 K. 12 unit tests pass. The scorecard shows 5 rows in
    the browser. `just verify` -> all gates passed; 483 Rust tests. See
    ADR-0051.

- [x] **M9-20** Add the harmonic analyze stage with the mesh normal modes.
  - Result: `crates/engine/src/hessian.rs` adds a central-difference
    Hessian, a Jacobi symmetric eigensolver, a mass weighting and a
    `Spectrum`. `crates/meter/src/harmonic.rs` adds the `mesh_mode`
    metric at `Fidelity::Harmonic`. It cuts a cluster of up to 180 atoms
    from the sun and the first planet around their closest approach,
    builds a `System` with a 300 N/m bond stretch on the real bonds and a
    Buckingham van der Waals term, relaxes the cluster, and reports the
    softest internal mode. The scorecard now shows six metrics.
  - Verification: `cargo run -p nanocad-meter --example score_json` ->
    mesh mode 1.184e13 Hz (395.1 per cm) over 148 atoms and 157 bonds,
    with 0 unstable modes. `cargo test -p nanocad-engine` -> 112 pass;
    `cargo test -p nanocad-meter` -> 15 pass. `just verify` -> all gates
    passed; 489 Rust tests. See ADR-0052.

- [x] **M9-21** Add the optimize stage with a parameter search.
  - Result: New crate `nanocad-opt` with a self-contained, deterministic
    CMA-ES and a planetary objective. `ParameterBounds` clamps each
    parameter and keeps an integer parameter whole. `Objective` maps a
    parameter vector to a cost, and returns `None` for a design the
    generator refuses. `PlanetarySearch` searches the sun teeth, the
    planet teeth and the planet count; the objective is the slip barrier
    in kT plus a clearance-shortfall penalty. The app route
    `/api/optimize` runs the search, and the panel has an Optimize
    section that shows the best design and applies it.
  - Verification: `cargo run --release -p nanocad-opt --example
    optimize_planetary -- population=8 generations=10 seed=5` -> best 12
    sun teeth, 10 planet teeth, 4 planets at 8.92 kT from a start of
    13.84 kT (36 percent lower), clearance 2.8865e-10 m passes, 178073
    atoms, 80 evaluations and 35 rejected. `cargo test --workspace` ->
    496 passed (7 new opt tests, 1 ignored). The browser clicks the
    Optimize button and shows the best design with an apply button. See
    ADR-0053.

## M10 — Part vocabulary and general molecular shapes

M9 designs one part: an involute gear. A molecular machine needs many
parts that connect. M10 adds a general solid, a port, and a lattice
defect. Every later part generator uses them.

- [x] **M10-01** Add a solid-shape trait with primitives and boolean
  operations. A `Solid` answers one question: is a point inside? The
  shapes are a box, a hex prism, a cylinder, and a placed (rotated)
  shape. The operations are difference, intersection, and union. Every
  shape returns its bounding box, so the lattice filler can enumerate
  sites without a scan of the whole cell.
  - Result: `crates/parts/src/shape.rs` adds the `Solid` trait with
    `contains_m` and `bounds_m`, the `Bounds` box, and the shapes `Box3`,
    `HexPrism`, `Cylinder`, `Placed`, `Difference`, `Intersection`, `Union`
    and `Profile`. The trait carries `clone_box`, because a boxed trait
    object is not `Clone`. A hex prism states its circumradius and its
    rotation, and a profile states a polygon plus an internal-rim mode.
  - Verification: 8 unit tests pass, including the hex six-side test, the
    placed rotation test and the difference hole test. `cargo clippy -p
    nanocad-parts --all-targets -- -D warnings` is clean.
- [x] **M10-02** Fill any solid with the diamond lattice. Generalise
  `diamond_solid::fill_profile` to a solid. The gear profile becomes one
  `Solid` implementation, so the gear set and the new parts share one
  filler and one bond-and-cap pass.
  - Result: `crates/parts/src/lattice_fill.rs` adds `fill_solid` and the
    shared `lattice_atoms` enumeration. `diamond_solid::fill_profile` now
    calls `lattice_atoms` with the polygon test, so the gear path and the
    general path place atoms on the same diamond sites. A solid with an
    unbounded axis gives an empty list instead of a hang.
  - Verification: 4 unit tests pass: the box fills, the nearest site is
    the diamond bond length (relative error below 1e-6), a difference
    removes the inner sites, and an unbounded profile returns nothing.
    The 13 existing diamond-solid tests still pass after the refactor.
- [x] **M10-03** Add a port model and a part connector. A port is a
  frame on a part: an origin, an axis, one allowed degree of freedom, and
  a stated interface gap. `connect` returns the transform that mates two
  ports. This is how an axle meets a bearing and a clutch meets a shaft.
  - Result: `crates/parts/src/port.rs` adds `Dof`, `Port`, `PortFrame`,
    `connect` and `mate_offset_m`. A plug meets a socket at the sum of the
    two gaps, and the plug axis turns to face the socket axis. The map is
    planar, so a vertical axis takes no rotation.
  - Verification: 6 unit tests pass, including the gap sum, the extra
    separation, the x-axis yaw flip and the vertical no-rotation case.
    Clippy is clean.
- [x] **M10-04** Add a wedge disclination for the 5-7 lattice defect.
  A wedge disclination removes a 60 degree wedge from the lattice. It
  turns two six-membered rings into a five-membered ring and a
  seven-membered ring. This defect curves a surface, so a tooth flank or
  a bearing race is not a staircase.
  - Result: `crates/parts/src/dislocation.rs` adds `WedgeDisclination`,
    `five_seven_wedge_rad`, `displace_point_m`, `displace` and
    `core_atom_mask`. The angular map is `theta / (1 - strength)`; the
    radius and z stay. A singular wedge and a core atom pass through
    unchanged, so no value becomes a NaN.
  - Verification: 7 unit tests pass, including the radius-and-z
    preservation, the angle stretch, the singular refusal and the core
    mask. Clippy is clean.
- [x] **M10-05** Generate a hex axle and a plain shaft. The cross section
  is a hexagonal prism, because a hex axle passes torque with a flat
  face. The generator states the across-flats size, the length, and the
  end chamfer.
  - Result: `crates/parts/src/axle.rs` adds `HexAxleGenerator` (id
    `hex_axle`) and `PlainShaftGenerator` (id `plain_shaft`). The hex
    prism circumradius is `across_flats_m / sqrt(3)`. The chamfer is an
    axial shortening, stated in the module doc as not a conical cut. The
    default hex axle has 19537 atoms and 33140 bonds; the shaft has 7309
    atoms and 11828 bonds.
  - Verification: 5 unit tests pass. The radius and length tests read
    the filled carbon atoms, because a hydrogen cap sits up to the C-H
    length outside the surface. `cargo test --workspace` -> 542 passed.
    Clippy is clean.
- [x] **M10-06** Generate structural blocks: a plate, a beam, and a
  bracket. A beam states its section and its length. A bracket is an
  L shape with a stated leg length and thickness.
  - Result: `crates/parts/src/block.rs` adds `PlateGenerator` (`plate`),
    `BeamGenerator` (`beam`) and `BracketGenerator` (`bracket`). One
    private `part_from_solid` helper fills the shape. A plate hole is a
    subtracted cylinder. Default counts: plate 6097 atoms and 9524
    bonds, beam 3677 and 5476, bracket 2825 and 4180.
  - Verification: 5 unit tests pass. The hole test reads carbons only,
    because the caps point into the void. Clippy is clean.
- [x] **M10-07** Generate a radial bearing and a bushing. A bearing has
  an inner race, an outer race, and a stated radial gap. The rolling
  elements are a stated count of cylinders. A bushing is a plain sleeve
  with a stated clearance.
  - Result: `crates/parts/src/bearing.rs` adds `RadialBearingGenerator`
    (`radial_bearing`) and `BushingGenerator` (`bushing`), plus the
    `inner_port`, `outer_port` and `shaft_port` helpers. A bearing with
    an outer radius below `inner + 2 * gap` is refused with
    `PartError::InvalidGeometry`. Default counts: bearing 19480 atoms
    and 26876 bonds, bushing 9668 and 13996.
  - Verification: 6 unit tests pass, including the port degree of
    freedom and the refusal. Clippy is clean.
- [x] **M10-08** Generate a clutch plate and a ratchet. The clutch is a
  stack of toothed plates with a stated gap. The ratchet is a pawl and a
  toothed wheel with a stated tooth count.
  - Result: `crates/parts/src/clutch.rs` adds `ClutchPlateGenerator`
    (`clutch_plate`) and `RatchetGenerator` (`ratchet`). The plate is a
    square-tooth annular outline filled through a `Profile` clipped by a
    `Box3`, with the centre removed. The ratchet is a one-way saw-tooth
    wheel plus a rigid block pawl. Default counts: clutch plate 3760
    atoms, ratchet 13538 atoms.
  - Verification: 5 unit tests pass. The radius test reads carbons only,
    because the caps sit past the surface. Clippy is clean.

## M11 — Loaded contact and driven motion

- [x] **M11-01** Build the full bonded terms in the meter potentials. The
  engine has torsion and out-of-plane terms. The meter builds bonds and
  van der Waals only. Add the missing terms with real parameters.
  - Result: New module `crates/meter/src/bonded.rs` adds `add_torsions`
    (every three-bond dihedral path, with `TORSION_V3_J = 2.0e-20` about
    12 kJ per mole) and `add_impropers` (an out-of-plane term with
    `IMPROPER_K_J_PER_RAD2 = 2.0e-18`). The improper is added only at a
    three-coordinate centre whose measured out-of-plane angle is within
    `TRIGONAL_TOLERANCE_RAD = 0.5`, so a tetrahedral centre stays out.
  - Verification: `cargo test -p nanocad-meter` -> 21 passed, including the
    relaxation test `a_relaxed_torsion_has_no_negative_curvature`, which
    relaxes a small torsion system and finds 0 unstable modes, so the
    engine torsion energy and its analytic gradient agree. `cargo test
    --workspace` -> 548 passed. `just verify` -> all gates passed. See
    ADR-0054.
  - Limit: the terms are built and tested, but they are not wired into the
    `mesh_mode` Hessian yet. The `mesh_mode` cluster is a cut-out of a
    larger gear, so its boundary carbons have fewer neighbours. The 536
    torsions over that cut-out add about -8e-6 N/m of curvature, which
    overlaps the genuine soft-mode band. A whole part is the correct host;
    see task M11-04.
- [x] **M11-02** Add a loaded-contact friction metric. The slip barrier
  has no normal load. Press two surfaces together, slide them, and report
  the friction force and the wear.
  - Result: New metric `contact_friction` in `crates/meter/src/loaded.rs`. A
    bisection finds the press offset that makes the normal force equal the
    stated load, then the sweep slides the pair through one sun tooth pitch.
    The lateral force is the energy fall over the arc length. The report holds
    the mean force, the peak force, the coefficient, the wear pair count, the
    contact pair count and the pressed separation. The coherence count comes
    from the same potential terms as the slip barrier. The bodies are rigid, so
    the value is a screening estimate.
  - Verification: `cargo test -p nanocad-meter` -> 26 passed, including
    `a_pressed_pair_has_friction` and
    `a_heavier_load_presses_the_bodies_together`. The default set reports a
    mean friction force of 1.6766e-10 N at a 1.0e-11 N load, 10 wear pairs and
    a pressed separation of about 3e-10 m. The app scorecard shows seven
    metrics. `just verify` -> all gates passed; 553 Rust tests. See ADR-0055.
- [x] **M11-03** Add a steered joint drive with a torque report. Hold one
  port and drive another at a stated rate. Report the torque, the energy
  loss, and the temperature rise.
  - Result: New metric `steered_drive` in `crates/meter/src/drive.rs`, at
    `Fidelity::Dynamics`. It builds a contact cluster, holds the planet atoms
    at their relaxed positions, and writes the sun velocity from the drive
    rate at every velocity-Verlet step. The report holds the mean and peak
    resisting torque, the drive work, the energy loss, the free-atom count,
    the final temperature and its rise. The example
    `crates/meter/examples/drive_json.rs` prints the report as JSON. The
    metric is NOT in the app scorecard, because it is a dynamics metric and
    the fit is a screening estimate.
  - Verification: `cargo test -p nanocad-meter` -> 35 passed, including
    `the_work_is_the_torque_times_the_angle` and
    `a_closed_body_has_no_torque`. On the default set the cluster holds 66
    atoms and the mean torque is 1.99e-45 N m, so the non-bonded contact at
    the relaxed mesh is nearly torque-free. The cluster holds every atom (the
    ring is about 1e-8 m away and out of the 6e-10 m collection radius), so
    the note states that the temperature rise is not measurable. `just verify`
    -> all gates passed; 562 Rust tests. See ADR-0056.
- [x] **M11-04** Report the normal modes of a whole part. The harmonic
  metric uses a 148-atom cluster. Extend it to a complete part.
  - Result: `HarmonicMesh::measure_part(&self, part: &nanocad_model::Part)`
    in `crates/meter/src/harmonic.rs` measures a whole part instead of a
    cut-out cluster. It returns the same `mesh_mode` metric at the same
    harmonic fidelity: the softest internal mode above 1e9 Hz. Because a whole
    part has no dangling crystal boundary, the potential now includes the
    torsion and out-of-plane terms from `crates/meter/src/bonded.rs`. The van
    der Waals term excludes every bonded pair and every pair inside the
    `COVALENT_EXCLUSION_M = 3.1e-10` m covalent shell, and the neighbour list
    is rebuilt after the relaxation and before the Hessian. A part above
    `target.max_atoms` is refused with an empty report.
  - Verification: `cargo test -p nanocad-meter` -> 35 passed, including
    `a_small_part_has_a_spectrum`, `a_part_above_the_cap_is_refused` and
    `the_whole_part_mode_is_at_the_harmonic_fidelity`. The smallest plate
    (477 atoms, 612 bonds) gives `lowest_hz = 5.733e10 Hz` with 1 unstable
    mode in a release run of about 215 s, against the 21 unstable modes of the
    cut-out cluster; this confirms ADR-0054. The tests use small uncapped
    diamond blocks because the diagonalisation is cubic in the atom count.
    `just verify` -> all gates passed; 562 Rust tests. See ADR-0057.

## M12 — Design workflow

- [x] **M12-01** Add a node selection language for atoms, bonds and parts.
  - Result: `crates/model/src/selection.rs` adds `parse_selection(text) ->
    Result<Selection, ModelError>` and `Selection::evaluate(&Document) ->
    SelectionResult { parts, atoms, bonds }`. The grammar is a boolean
    expression over predicates: `and`, `or`, `not`, parentheses, `all`,
    `none`, `atom`, `bond`, `part`, the element names, and the comparisons
    `charge`, `degree`, `atom_index`, `order`, `bond_length`, `bond_index`,
    `part_index`, `name == "..."` and `within <d> of atom [i]`. Keywords are
    case-insensitive. An expression evaluates to one set per node kind, so
    `and`, `or` and `not` stay per kind. A parse error is
    `ModelError::SelectionParse { position, message }`.
  - Verification: `cargo test -p nanocad-model` -> 54 passed, including
    `carbon_and_degree_selects_the_middle_carbon`,
    `bond_length_selects_by_length`,
    `parse_errors_carry_a_position_and_do_not_panic` and
    `within_and_of_atom_stays_in_one_part`. 11 new tests. `just verify` -> all
    gates passed; 581 Rust tests. See ADR-0058.
- [x] **M12-02** Add a part library and an interactive placement tool.
  - Result: `crates/parts/src/registry.rs` is now a general library. It
    registers all 18 generators in four categories (`gears`, `lattice`,
    `structure`, `device`) and exposes `library()`, `generator(name)`,
    `generate(name, specs)` and `library_categories()`. The gear-only
    functions stay. `crates/parts/src/placement.rs` adds `place(part, plug,
    socket, separation_m) -> PlacedPart` and
    `assembly_document(&[PlacedPart]) -> Document`. A placed part carries the
    port frame, and `transformed_part` moves every atom into the socket frame.
  - Verification: `cargo test -p nanocad-parts` -> 142 passed, including
    `the_library_lists_every_generator_with_a_unique_id`,
    `every_library_entry_names_a_generator`,
    `a_placed_part_moves_to_the_socket` and
    `an_assembly_holds_every_placed_part`. 8 new tests. `just verify` -> all
    gates passed; 581 Rust tests. See ADR-0059.
- [ ] **M12-03** Store a design document with the parameters, the parts
  and the measurements, and add multi-level undo and redo.
