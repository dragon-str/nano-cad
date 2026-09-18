# Decision Log (ADRs)

Record every architectural decision here. One entry per decision. Do not delete
old entries. Supersede them with a new entry.

Template:

```
## ADR-NNNN: Title
- Status: proposed | accepted | superseded by ADR-NNNN
- Date: YYYY-MM-DD
- Context: why this decision was needed
- Decision: what we decided
- Consequences: what this enables and what it costs
```

## ADR-0001: Fresh engine, Apache-2.0
- Status: accepted
- Date: 2026-09-13
- Context: The upstream NanoEngineer C engine is GPLv2-or-later. A fork would be
  GPL. We want the option to build a commercial tier and to partner freely.
- Decision: Write a fresh engine in clean code under Apache-2.0. Read the old C
  engine for algorithms only. Do not copy its code.
- Consequences: We avoid the copyleft constraint. We lose the existing working
  code and must rebuild it. We gain a clean, modern codebase.

## ADR-0002: Rust core with Python bindings
- Status: accepted
- Date: 2026-09-13
- Context: We need speed, memory safety, threads, and an easy scripting path.
- Decision: Write the core in Rust. Expose a Python package through PyO3. Put
  the CLI in Rust.
- Consequences: Speed and safety. The scientific ecosystem is reachable from
  Python. Rust packages must be built for each platform.

## ADR-0003: SI units internally
- Status: accepted
- Date: 2026-09-13
- Context: The old code mixed 0.01 Angstrom, picometre, and nanometre units.
  That caused bugs.
- Decision: Store all quantities in SI inside the core. Convert at the
  boundary. Use a units type. Never mix units in one expression.
- Consequences: Fewer unit bugs. Some conversion overhead at the boundary.

## ADR-0004: Reuse OpenMM and ASE; do not rebuild them
- Status: accepted
- Date: 2026-09-13
- Context: Standard MD forces, GPU support, and file interop are solved.
- Decision: Use OpenMM (MIT) for GPU and standard force fields. Use ASE (LGPL)
  for QM and file interop. Use PySCF (BSD) or xTB (LGPL) for quantum work. Keep
  our engine for atomically precise parts and jigs.
- Consequences: We move fast and reach many formats. We depend on their
  licenses and their release cycles. LGPL use must stay dynamic-link only.

## ADR-0005: Open document format, with NanoEngineer interop
- Status: accepted
- Date: 2026-09-13
- Context: The tool needs a native format and must trade parts with old tools.
- Decision: Define `NCZ`, an open, versioned document format (metadata plus
  binary arrays). Import and export MMP (NanoEngineer). Export XYZ, PDB, and
  URDF.
- Consequences: No lock-in. Users can migrate. The format spec is a published
  asset.

## ADR-0006: The multiscale parameter schema is the central contract
- Status: accepted
- Date: 2026-09-13
- Context: The unique value is the bridge between scales.
- Decision: Define one versioned parameter schema. Each level reads and writes
  it. Every parameter carries provenance and uncertainty. See `PARAMETERS.md`.
- Consequences: Levels are decoupled. The schema must be versioned carefully.

## ADR-0007: Agent interface through a headless API and an MCP server
- Status: accepted
- Date: 2026-09-13
- Context: Users want to design by talking to the tool.
- Decision: Expose one headless API. The GUI, the CLI, and the agent all call
  it. Publish an MCP server as the agent bridge. The agent calls generators and
  solvers. It never edits coordinates directly.
- Consequences: One code path for all clients. The agent stays reliable because
  it calls typed tools.

## ADR-0008: Verification is a gate, not a phase
- Status: accepted
- Date: 2026-09-13
- Context: Scientific trust is the product.
- Decision: Finite-difference gradient checks and energy-conservation runs are
  merge gates. Keep reference benchmarks against OpenMM and GROMACS.
- Consequences: Slower early development. Higher trust later.

## ADR-0009: Build order is a vertical slice; the gearbox is the first demo
- Status: accepted
- Date: 2026-09-13
- Context: The original project failed from breadth.
- Decision: Build one planetary gearbox through every level before going broad.
- Consequences: Early end-to-end proof. Some layers stay thin at first.

## ADR-0010: The verified-data library is governed, not proprietary
- Status: accepted
- Date: 2026-09-13
- Context: Facts are not copyrightable. Curated, provenance-tracked data is the
  durable asset.
- Decision: Publish parameters under a governed open data license (CC-BY) with
  provenance. Sell certification and hosting, not the raw facts.
- Consequences: Trust and reuse grow. Revenue comes from services.

## ADR-0011: The model byte encoding is private; NCZ is the public format
- Status: accepted
- Date: 2026-09-14
- Context: `Part` and `Document` need a fast in-memory round trip. The public
  file format (NCZ) is a separate, versioned concern (M1-03). If the model
  codec and NCZ were the same bytes, every internal layout change would break
  the public format.
- Decision: Give the model a private little-endian codec with an `NCAD` magic
  and a schema version. Keep NCZ as the published format in the `format` crate.
  A migration function converts between them.
- Consequences: Internal layout can change without a public format break. We
  must write and test the conversion. The `NCAD` magic never appears in a file
  a user shares.

## ADR-0012: Python bindings stay a plain crate until M7
- Status: accepted
- Date: 2026-09-14
- Context: The `python` crate is listed with PyO3 bindings, but the bindings
  have no API to expose until the tool surface exists (M7-01). Adding PyO3 now
  would force a Python development environment into every `cargo test` run.
- Decision: Keep `crates/python` a plain library crate for now. Ship the
  Python package as a pure-Python module with a typed `version()`. Add PyO3 at
  M7-01, when there are tools to bind.
- Consequences: `cargo test --workspace` stays free of the Python interpreter.
  M7-01 must add the `pyo3` dependency and the extension module.

## ADR-0013: Periodic minimum image lives in the neighbor list and the non-bonded terms
- Status: accepted
- Date: 2026-09-14
- Context: `VerletList` and the non-bonded terms originally ignored a periodic
  box. `System` was correct only for a non-periodic cell. The water-box
  benchmark (M2-11) needs a periodic box.
- Decision: Give `VerletList` an optional `PeriodicBox`. Apply the minimum
  image in the neighbor search and in the per-pair non-bonded distance. Keep
  the non-periodic behavior the default. Propagate the box to the non-bonded
  terms when `System::set_periodic_box` is called.
- Consequences: The force terms and the neighbor list use one box convention.
  A term replaced after the box is set is not re-checked; the caller must set
  the box again.

## ADR-0014: Device constraints use position-level projection; gear coupling is a constraint
- Status: accepted
- Date: 2026-09-14
- Context: M6 needs revolute, prismatic, and gear constraints on rigid bodies.
  A velocity-level solver is more accurate but larger. The gear coupling is not
  a potential, so a force term cannot express it.
- Decision: Use semi-implicit Euler with world angular momentum as the primary
  angular state. Enforce constraints with a position-level Gauss-Seidel,
  mass-weighted projection. Model the gear coupling as an exact linearized
  angle constraint `r_a dtheta_a + r_b dtheta_b = -error`. There is no
  velocity-level correction pass yet.
- Consequences: Constraint error stays near machine precision for the tested
  cases. A residual constraint velocity can persist. A velocity-level pass
  stays open for later.

## ADR-0015: URDF export lives in the jigs crate
- Status: accepted
- Date: 2026-09-14
- Context: M6-07 exports the device as URDF. The device model (`RigidBodySystem`
  and its joints) lives in `crates/jigs`. The `format` crate holds file formats
  for molecules, not devices.
- Decision: Put the URDF writer and a small reader in `crates/jigs/src/urdf.rs`.
  Do not add an XML dependency.
- Consequences: The exporter travels with the device model it serializes. If a
  later device format appears, this placement is revisited.

## ADR-0016: Python bindings are feature-gated
- Status: accepted
- Date: 2026-09-14
- Context: ADR-0012 deferred PyO3 to M7-01. This is M7-01. The workspace
  `cargo test` and `cargo clippy` gates must still run without a Python
  interpreter.
- Decision: Add `pyo3` as an optional dependency behind a `python` feature in
  `crates/python`. The default build stays a plain rlib. Maturin enables the
  feature through `python/pyproject.toml`.
- Consequences: The Rust gates stay interpreter-free. The extension builds only
  through maturin with the feature on. CI must keep both paths green.

## ADR-0017: Deterministic parallel reduction uses std threads
- Status: accepted
- Date: 2026-09-14
- Context: M3-08 needs a multi-threaded force reduction that is reproducible
  bit for bit. A work-stealing pool would not be deterministic.
- Decision: Use `std::thread::scope` with fixed chunks and a deterministic
  reduction order. Do not add a dependency. Keep the serial path the default.
- Consequences: The result is independent of the thread count. The parallel
  path is a correctness path, not the fastest path.

## ADR-0018: Thermal extraction uses a kinetic estimate and a two-bath direct method
- Status: accepted
- Date: 2026-09-14
- Context: M5-05 extracts `specific_heat_J_kg_K` and
  `thermal_conductivity_W_m_K`. A full Green-Kubo heat flux needs the per-term
  virial, which the engine does not expose.
- Decision: Estimate the specific heat from kinetic-energy fluctuations in the
  NVT ensemble, `C_v = Var(E_k) / (k_B T^2)`. Estimate the conductivity with a
  two-bath direct method, `q * distance / (area * delta_T)`. Mark both
  `Method::Md` and `Validation::Unverified`.
- Consequences: The specific heat is a kinetic estimate; it omits `Var(U)`.
  The conductivity is an effective value for the stated geometry; a harmonic
  chain shows anomalous transport. Both are estimates, not material constants.

## ADR-0019: The MCP server uses the Python standard library only
- Status: accepted
- Date: 2026-09-14
- Context: M7-02 wraps the Python API as an MCP server. The official `mcp`
  package would add a dependency and a supply-chain surface for a small
  JSON-RPC protocol.
- Decision: Hand-roll a JSON-RPC 2.0 stdio server in `mcp/nanocad_mcp/`. Keep
  the tool schemas in `mcp/tools.json`. Add no third-party package.
- Consequences: No license change. State is in memory, so a restart loses
  handles. The server emits a fixed tool list and no list-changed
  notification.

## ADR-0020: The three-scale scene export lives in `nanocad-jigs`
- Status: accepted
- Date: 2026-09-14
- Context: M8-01 needs one export that a visualization can animate across the
  atomistic (L1), device (L2), and coarse scales. The device, the assembly, and
  the URDF export already live in `crates/jigs`.
- Decision: Put `scene.rs` in `crates/jigs`. Serialize a `nanocad.scene` JSON
  document with three layers. Zero positions are metres. Angles are radians.
  Quaternions are `[w, x, y, z]`. The frame is right-handed with z as the gear
  axis. Add `serde` and `serde_json`, both already used in the workspace.
- Consequences: `crates/format` stays free of device dependencies. The atomistic
  layer is a static snapshot with no bonds. A prismatic joint reports a null
  world anchor. Every coarse cylinder aligns with the z axis.

## ADR-0021: The clean-clone reproduction is a shell script plus a `just` recipe
- Status: accepted
- Date: 2026-09-14
- Context: PLAN.md M8 requires that a stranger reproduce the demo from the
  repository. The demo spans Rust, the maturin extension, and the MCP test.
- Decision: `scripts/reproduce.sh` runs the seven steps and exits nonzero on
  any failure. `just reproduce` calls it. `docs/reproduce.md` documents the
  prerequisites. A guard test asserts a missing `cargo` fails the script.
- Consequences: The first run needs network access to install maturin and
  pytest. A later run reuses the venv. The video and the written report stay
  outside the script.

## ADR-0022: Velocity Verlet takes an optional external-force buffer; bonds are readable
- Status: accepted
- Date: 2026-09-14
- Context: The friction extractor held its own velocity-Verlet loop to add a
  drive spring, a load, and drag. The failure extractor required the caller to
  build bond limits by hand, because `System` exposed no bond list.
- Decision: `VelocityVerlet::step_with_external_forces` adds a constant
  per-atom external force over one step. `System::bond_count`, `System::bond`,
  `System::bond_forces_n`, and `System::bond_strains` expose the bonds and
  their scalar quantities. `BondStretchTerm::bond_info` is the source.
  `bond_limits_from_system` derives the failure limits from the system.
- Consequences: The friction and failure extractors use the public engine path.
  The external force is constant over the step, so the method is exact for it,
  and no finite-difference test applies to that step. Bond forces are
  magnitudes, not signed vectors.

## ADR-0023: Python exposes the parameter extractors as properties, not one call
- Status: accepted
- Date: 2026-09-14
- Context: M5-02 to M5-05 produce five property results. M7-01 left them out
  of the Python surface. One `extract_parameters` call would need a config
  object on the Python side.
- Decision: Expose per-property functions: `extract_stiffness`,
  `extract_friction`, `summarize_friction`, `extract_failure_stress`,
  `extract_failure_stress_from_system`, `extract_specific_heat`,
  `extract_thermal_conductivity`, and `extract_thermal`. Each returns a result
  class with `Quantity` fields and provenance.
- Consequences: The functions take many keyword arguments and no Python config
  class exists yet. The `python` feature still gates every binding (ADR-0016).

## ADR-0024: L0 quantum adapter is a lazy Python slice over PySCF or ASE
- Status: accepted
- Date: 2026-09-14
- Context: L0 needs a quantum energy, but a hard PySCF dependency would break
  the default install and the test suite.
- Decision: `python/nanocad/quantum_adapter.py` imports every backend inside
  the call and takes plain atom data. PySCF is preferred, ASE is a fallback,
  and the result states `is_quantum`, `method` and `validation`.
- Consequences: The default test run skips the backend tests. The adapter is a
  single-point energy only; it returns no forces and no gradient.

## ADR-0025: L3 continuum adapter uses closed-form checks, numpy only
- Status: accepted
- Date: 2026-09-14
- Context: L3 needs structural and flow results that a reviewer can verify.
- Decision: `python/nanocad/continuum_adapter.py` solves an axial bar, an
  Euler-Bernoulli beam, and Poiseuille tube and channel flow, and checks each
  against the analytic closed form. It uses numpy and the standard library.
- Consequences: The models are 1D, steady, and laminar. numpy is a test extra
  in `python/pyproject.toml` and the reproduce script installs it.

## ADR-0026: L4 lumped model emits an honest SBML subset
- Status: accepted
- Date: 2026-09-14
- Context: L4 needs a lumped model and a standard interchange format. Full SBML
  support is large and would need a dependency.
- Decision: `python/nanocad/lumped_adapter.py` emits and parses an SBML Level 3
  subset with `xml.etree.ElementTree`. It supports compartments, species,
  parameters, and reactions with a MathML kinetic law.
- Consequences: Foreign SBML that uses assignment rules, events, or function
  definitions is rejected with a typed error. Units ride in a
  `nanocad:unit` extension attribute.

## ADR-0027: Respirocyte subsystems live in nanocad-parts as lattice generators
- Status: accepted
- Date: 2026-09-14
- Context: FL-01 and FL-02 need rotor, bearing, pump, and tank geometry that
  the part machinery can generate and validate.
- Decision: Add `crates/parts/src/respirocyte.rs` with three
  `PartGenerator`s over the diamond lattice, each returning Parts plus a
  metadata struct with SI values.
- Consequences: The geometry is skeletal, not a validated device. No force term
  was added, so no finite-difference gradient test applies.

## ADR-0028: NM-01 couples a rigid body to a coarse Stokes host
- Status: accepted
- Date: 2026-09-14
- Context: NM-01 needs a machine in a host environment without a full fluid
  solver.
- Decision: Add `crates/jigs/src/nanomedicine.rs` with `HostEnvironment`
  (Stokes drag, optional flow and tether) and `NanoMachine`, integrated by the
  existing `RigidBodySystem`. The slice is checked against terminal velocity.
- Consequences: The host is a single-sphere Stokes model with no wall
  correction, no Brownian noise, and no hydrodynamic interaction. The energy
  balance closes only to the integrator order.

## ADR-0029: L-BFGS is a second minimization stage with explicit selection
- Status: accepted
- Date: 2026-09-14
- Context: Conjugate gradient alone stalls on stiff systems. ENG-01 needs a
  tighter gradient without changing the default behavior of `minimize`.
- Decision: Add `MinimizeMethod` and `minimize_with` in
  `crates/engine/src/minimize.rs`. `ConjugateGradientThenLbfgs` is the new
  default for the new entry point. The old `minimize` stays conjugate-gradient.
- Consequences: The callers choose the method. The hybrid reaches 5.9049e-15 N
  on the stiff test. The library keeps one behavior change, at a named call.

## ADR-0030: The parallel force reduction uses a canonical pair-block partition
- Status: accepted
- Date: 2026-09-14
- Context: A naive parallel reduction was non-deterministic, so tests and
  simulations could not repeat. ENG-02 needs both determinism and speed.
- Decision: Partition the pair list into canonical blocks and sum the same
  per-block partials in serial and in parallel. The serial path also sums the
  partials, so the two paths agree bit for bit.
- Consequences: The result is deterministic for one, two, four and eight
  workers. The serial path changed, so its old accumulation order is gone. The
  speedup is about 1.98x on eight workers on an Apple M4.

## ADR-0031: The nonbonded kernel is scalar plus a portable multi-lane path
- Status: accepted
- Date: 2026-09-14
- Context: ENG-03 needs a faster inner loop. `std::simd` is the natural tool
  but it is nightly-only, and the MSRV is 1.85 on a stable rustc.
- Decision: Write a four-lane unrolled Buckingham plus electrostatic kernel and
  keep a scalar fallback. Both paths are bit-identical.
- Consequences: The crate builds on stable. The `std::simd` version is deferred
  until it reaches stable or the MSRV moves. The hand-written lanes are more
  code to maintain.

## ADR-0032: nanocad-format imports and exports PDB
- Status: accepted
- Date: 2026-09-14
- Context: IO-01 needs to exchange structures with the wider molecular tooling.
- Decision: Add `crates/format/src/pdb.rs` with `import_pdb` and `export_pdb`
  over ATOM, HETATM, CONECT and CRYST1. Map lengths at the Angstrom boundary.
- Consequences: Charge is lost because PDB has no charge field. Multi-part
  documents do not round-trip. A malformed line raises `FormatError::PdbParse`
  with the line number.

## ADR-0033: RDKit is an optional Python adapter
- Status: accepted
- Date: 2026-09-14
- Context: IO-02 needs cheminformatics interchange without making RDKit a core
  dependency.
- Decision: Put the RDKit conversion in `python/nanocad/rdkit_adapter.py`. The
  import is lazy and the tests skip when RDKit is absent.
- Consequences: The core stays free of the dependency. The adapter is tested in
  `/tmp/nc-qm-venv`, which has RDKit 2026.03.6, and skipped in the plain venv.

## ADR-0034: The lumped model exports an honest CellML subset
- Status: accepted
- Date: 2026-09-14
- Context: IO-03 needs a second standard interchange format for L4.
- Decision: Add `to_cellml` and `from_cellml` in
  `python/nanocad/lumped_adapter.py`. The subset has no reaction element. The
  units attribute is dimensionless and a `nanocad:unit` attribute carries the
  true unit.
- Consequences: Foreign CellML that uses unsupported elements is rejected with
  a typed error. The unit carries through the round trip only under our
  extension attribute.

## ADR-0035: A fixed joint and a velocity-level constraint pass
- Status: accepted
- Date: 2026-09-14
- Context: DEV-01 needs a weld and a drift-free constraint solution.
- Decision: Add `FixedJoint` and `ConstraintOptions { velocity_pass }` in
  `crates/jigs/src/device.rs`, on by default. The pass projects the constraint
  velocity to zero for revolute, prismatic, gear and fixed joints.
- Consequences: The revolute constraint velocity falls from 1.118 to 0. The
  gear residual is 2.78e-17. The pass adds one linear solve for each step.

## ADR-0036: OpenMM and GROMACS are cross-checks, not engine backends
- Status: accepted
- Date: 2026-09-14
- Context: ENG-04 and VAL-01 need an independent number for the electrostatics
  and a second opinion on the potential. The engine stays the one core.
- Decision: Add `python/nanocad/pme_adapter.py` for OpenMM PME and
  `benchmarks/gromacs/` for the GROMACS water-box check. Neither is a backend.
  GROMACS 2026.3 cannot use Buckingham with the Verlet cutoff scheme, so that
  driver records SKIP and writes no number.
- Consequences: OpenMM PME agrees at 2.42e-4 relative at 10 nm. The GROMACS
  check is honestly open. No false number enters the docs.

## ADR-0037: The respirocyte seal is an annular-gap leakage model
- Status: accepted
- Date: 2026-09-14
- Context: FL-03 and NM-02 need a leak estimate and thermal motion without a
  full fluid solver.
- Decision: Add `AnnularGapFlow` to `crates/parts/src/respirocyte.rs` and
  Brownian motion to `HostEnvironment`. Use the pressure-driven and Couette
  closed forms from Bird, Stewart and Lightfoot.
- Consequences: The bearing leakage is 9.114807e-20 m^3/s and matches the hand
  value. The model is one-dimensional and laminar. The Brownian variance is
  within 3.46e-2 of kBT/m at the test temperature.

## ADR-0038: The video atom layer uses DiamondGenerator; the gear stays skeletal
- Status: superseded by ADR-0039
- Date: 2026-09-14
- Context: M9-02 needs a physically accurate atom layer. `PlanetaryGenerator`
  returns a skeletal gear outline with nearest-neighbour distances from
  5.28e-11 m to 3.72e-10 m and polygon angles, so it is not a bonded solid.
- Decision: Build the atom layer from `DiamondGenerator` in
  `crates/parts/src/lattice.rs`, and gate it with
  `scripts/check_atom_geometry.py`. State the limitation in `docs/video.md`.
- Consequences: The atom layer is real diamond: 216 atoms, 333 bonds, mean C-C
  distance 1.544556e-10 m and mean bond angle 109.471221 degrees. The gear shape
  is still not chemically accurate. A diamondoid gear generator with surface
  passivation is a separate research task.

## ADR-0039: The schematic draws the exact involute profile and the gears are atoms
- Status: accepted
- Date: 2026-09-14
- Context: The video's schematic gears and the underlying atomic layer did not
  agree. The drawing used a crude four-point trapezoid tooth with the tip at the
  pitch radius and the root too deep, and it drew the ring as an external gear
  with outward teeth. The atomic layer was a separate diamond block, not the
  gears. M9-02 and the user asked that the gears be made of atoms and that the
  schematic match the dimensions.
- Decision: Add `scripts/gear_profile.py`, a port of
  `crates/parts/src/gear_profile.rs`. The renderer and the checker import it, so
  one formula drives both. Draw the ring as an internal gear with the reflection
  `2 p - r`. Render the gears as their real `PlanetaryGenerator` atoms and bonds
  from `site/scene.json` and `site/scene.bonds.json`, animated with the same
  fixed-ring rates as the schematic.
- Consequences: The schematic tip and root radii now equal the atomic layer
  exactly: sun `5.375e-9` to `6.500e-9 m`, planet `3.875e-9` to `5.000e-9 m`,
  ring internal `1.450e-8` to `1.5625e-8 m`. `scripts/check_atom_geometry.py`
  asserts this and `scripts/make_video.sh` fails when it fails. The
  `PlanetaryGenerator` stays skeletal, so its bond lengths and angles are not
  diamond values. The `DiamondGenerator` block is kept as the verified material
  basis. Supersedes ADR-0038.

## ADR-0040: The gear generator extrudes the profile into atomic layers
- Status: accepted
- Date: 2026-09-14
- Context: The gears were a single plane of atoms, so they had no thickness and
  the video could not show a solid. The generator also placed the planet teeth
  on the sun and ring teeth, so the atom gear set collided.
- Decision: Add `layers` and `layer_spacing_m` to `PLANETARY_PARAMETERS` and an
  `add_extruded_loop` helper in `crates/parts/src/planetary.rs`. The sun, the
  planets, and the ring are copied into centered axial layers with vertical
  bonds. The planet rotation carries a half-tooth offset `pi/N_p`.
- Consequences: The default gear set has four layers and is `4.632e-10 m` thick.
  The atom layer meshes: a sun tooth enters a planet space and a planet tooth
  enters a ring space. The thickness is a stated number of atoms, not a render
  effect. The carrier stays on one plane.

## ADR-0041: One rate function drives the schematic and the atoms
- Status: accepted
- Date: 2026-09-14
- Context: The schematic and the atom layer each carried their own hand-coded
  kinematics. The atom layer double-counted the carrier, so the schematic and
  the atoms turned the planets at different rates.
- Decision: Derive both from one phase and rate function. With sun rate `w_s`,
  `w_c = w_s*N_s/(N_s+N_r)`, the planet relative rate is
  `w_rel = -(N_s/N_p)*(w_s - w_c)`, and the planet absolute rate is
  `w_p = w_c + w_rel`. The planet placement phase is `2*pi*k/P + pi/N_p`.
  `scripts/check_atom_geometry.py` measures the atom planet orientation rate
  and compares it to the schematic rate.
- Consequences: The schematic and the atom layer share one rate, checked to
  `1e-6 rad/s`. The renderer no longer keeps a second set of rates. A change to
  the kinematics must change one function.

## ADR-0042: The interactive app is a standard-library server with a rule-based parser
- Status: accepted
- Date: 2026-09-14
- Context: The user asked to change the design in words and with side-panel
  parameters. A language model in the loop would add a dependency and could
  invent a value.
- Decision: `app/server.py` is a `http.server` app that runs the Rust generator
  for each change. `app/chat.py` is a rule-based parser. Both use only the
  Python standard library. The page talks to the server under `/api/`.
- Consequences: The app needs no external service and no network. It accepts a
  stated set of commands only; other text returns help. Every value is clamped,
  and the engine error is shown verbatim. The panel and the chat show engine
  output, never an invented number. The parser has its own tests, and
  `just py-test` runs them.

## ADR-0043: The gears are solid hydrogen-capped diamond
- Status: accepted
- Date: 2026-09-15
- Context: The gears were a hollow extruded profile. A reader could not see the
  diamond material, and the axial spacing was a free parameter, so the layers
  were not a crystal.
- Decision: Fill each gear with a real diamond lattice cut to the involute
  profile in `crates/parts/src/diamond_solid.rs`. One atomic layer is one diamond
  (001) plane, so the plane spacing is fixed at `a/4 = 8.9175e-11 m` and only the
  layer count is a parameter. Bond the lattice at the C-C length, then cap every
  carbon that has fewer than four neighbours with a hydrogen at the C-H length
  (`1.09e-10 m`). The capping probe tests occupancy at the C-C length and scans
  the 27 neighbouring grid cells, because a neighbour can sit on a cell
  boundary.
- Consequences: Every diamond carbon is exactly four-bonded and every hydrogen
  exactly one-bonded; the default gear set has 32786 carbons and 37284 hydrogens
  and obeys the diamond bond angle `109.4712 deg`. The layer spacing is no longer
  a parameter. The mechanical carrier stays a degree-2 ring and is excluded from
  the chemistry checks.

## ADR-0044: The viewer level of detail switches on the atom spacing in pixels
- Status: accepted
- Date: 2026-09-15
- Context: The user asked to see real atoms when zoomed in and the schematic when
  zoomed out, with the two agreeing in size.
- Decision: `site/viewer.js` measures the projected carbon-carbon spacing in
  pixels. Above six pixels it draws depth-shaded atoms with the true radii; below
  six it draws the exact involute schematic. Both use the same profile function
  and the same world fit, so the switch does not move the geometry.
- Consequences: The user sees atoms at close zoom and a clean profile at far
  zoom. The schematic and the atom outline share one code path, so a change to
  the profile shows in both.

## ADR-0045: The atom layer renders with WebGL sphere impostors and screen-space ambient occlusion
- Status: accepted
- Date: 2026-09-15
- Context: The user asked for a render close to Philip Turner's published
  atom-level images and asked to remove the scale slider, add panning, and add
  the features that a viewer of this kind needs.
- Decision: `site/render_atoms.js` adds a WebGL2 renderer with no external
  library. It draws one instanced sphere impostor for each atom, with key, fill
  and rim light, a specular term, a contact-shadow bias, a screen-space ambient
  occlusion pass at half resolution, and a vignette. `site/viewer.js` drives it:
  pan, preset views, a turntable, a PNG export, fullscreen, an element legend
  with per-element visibility, an atom-size control, a z-section clip, a hover
  pick, a two-click distance measure, a nanometre scale bar and an orientation
  triad. The scale slider is removed because the atom layer already switches to
  the schematic by zoom.
- Consequences: The viewer looks closer to a real atomistic render and stays
  offline and dependency-free. The occlusion, the lighting and the sphere radii
  are display effects, so they carry no physical result. When WebGL2 is absent
  the viewer falls back to the flat 2D dots, and the page still works. The
  banner stays on the page.

## ADR-0046: The carrier carries no atoms, and the gears get tooth clearance

Status: accepted.

Context. The carrier was a ring of carbon atoms at the axial offset. The
scene showed it as a plane of atoms that floated above the gear. The
carrier ring also overlapped the planets in the plane. The solid gear
teeth also met too closely. A discrete diamond lattice cannot follow the
smooth involute curve, so atoms stand proud of the true profile. Two
meshing gears then overlapped.

Decision. The carrier is a scene body, and the part carries no carrier
atoms. The generator has eight parameters. The mass properties of the
carrier come from the assembly, not from atoms.

The gear profile gets tooth clearance. The addendum coefficient is 0.5,
the dedendum coefficient is 1.7, and the backlash is 9.0e-10 m. The
backlash thins each flank and the small addendum shortens each tooth
tip. `site/viewer.js` mirrors all three constants, because the
schematic must match the atom layer.

Consequences. The minimum atom-to-atom distance between two bodies is
2.8114 A over one carrier revolution. This is more than the 2.52 A
target, the nearest non-bonded diamond spacing. The teeth are shorter
than a standard gear. The contact ratio is lower, and the set is a
display model, not a working gearbox.

The clearance test is an off-line all-atom sweep, not a unit test. The
generator places the planet atoms with the assembly, so a unit test in
the parts crate cannot reproduce the body transforms. See
`/tmp/overlap.py` and `TASKS.md` M9-13.

## ADR-0047: Metrics live in a separate crate and carry a fidelity label

Status: accepted.

Context. The user wants the software to guide the design of nanoscale
parts and to optimise them. A design needs scores, and the scores come
from many methods at different cost and accuracy. A single number hides
that difference.

Decision. Add the crate `nanocad-meter`. A metric returns a
`MetricValue` with a name, a value, a unit, a `Fidelity` and a note.
`Fidelity` has four levels in cost order: `Geometric`, `QuasiStatic`,
`Harmonic` and `Dynamics`. A `Score` collects metric values.

The first metric is `Clearance`. It sweeps the relative motion of the
bodies and reports the least atom-to-atom distance between two bodies.
The default target is 2.52e-10 m, the nearest non-bonded diamond
spacing.

Consequences. A score is comparable only at one fidelity level. The
optimiser must not mix levels. Every metric needs a reference check, as
every force term needs a finite-difference gradient test.

The clearance metric found a real error on its first use. The planet
motion applies the spin in the world frame, so the absolute planet rate
is the carrier rate plus the rate relative to the carrier. The earlier
off-line sweep had this right; the first version of the metric did not.

## ADR-0048: The app shows a scorecard with a fidelity badge on each metric

Status: accepted.

Context. The metric layer must be visible to the user. A raw number hides
how it was measured, and cheap and expensive metrics must not be confused.

Decision. The app shows a Score panel with one row for each metric. Each
row shows the name, the value in a readable unit, the fidelity badge
(`geometric`, `quasi_static`, `harmonic`, `dynamics`) and the note. The
server route `GET /api/score` returns the score as JSON, and the metric
crate writes that JSON through `crates/meter/examples/score_json.rs`.

The first metrics are all `geometric`, because they read the generated
geometry. Later metrics carry a higher fidelity.

Consequences. A low-fidelity value is always labelled, so the user cannot
compare it with a high-fidelity value by accident. The first real finding
is that the default contact ratio is 0.86, below 1. The short addendum of
0.5, which the clearance work chose, reduces the contact ratio. The drive
can skip a tooth. This is a design trade-off to present to the user.

## ADR-0049: The default gear set uses the recommended nanoscale size

Status: accepted.

Context. The physics gives a smallest reasonable gear set. A tooth needs
at least four diamond lattice rows, so the module needs about 1.5 nm. A
small pressure angle undercuts a small planet, so 30 degrees is used.

Decision. The default parameters are module `1.5e-9 m`, 12 sun teeth, 9
planet teeth, 30 ring teeth, a 30 degree pressure angle, an addendum
coefficient of 0.8, a dedendum coefficient of 1.25 and a backlash of
`1.0e-9 m`. The set has 142091 atoms.

The mesh phase depends on the tooth-count parity. A planet presents a
space toward the sun when the planet tooth count is odd, and a tooth when
it is even. The ring takes a half-pitch rotation `pi/N_r` when the planet
tooth count is odd. Exactly one gear of the planet and ring pair takes the
half-pitch offset, so the three gears interleave.

The backlash of `1.0e-9 m` is larger than the van der Waals estimate of
0.6 nm. A backlash of 0.6 nm gave a dynamic minimum of `1.3e-10 m`, which
is an overlap for the 12 and 9 tooth pair. The clearance metric sets the
real value.

Consequences. The default set is 2.5 times the previous atom count, so a
build is slower. The contact ratio is 1.00, which is a working value. The
teeth are thin, because the backlash is large.

## ADR-0050: A quasi-static slip barrier measures the friction at the mesh

Status: accepted.

Context. The clearance metric proves the atoms do not overlap. It says
nothing about the friction. Two meshing surfaces can be clear and still
jam, if the potential energy swings by much more than the thermal energy
as one tooth passes. Stick-slip friction is the sign of that barrier.

Decision. The `nanocad-meter` crate gets a `slip_barrier` metric. The
metric sweeps one sun tooth pitch with the real gear kinematics. At every
step it sums a shifted Lennard-Jones interaction between the sun and the
first planet. The barrier is the largest energy minus the smallest.

The barrier is compared with the thermal energy `k T`. Below 5 `k T` the
thermal motion smooths the slip. Between 5 and 20 `k T` the slip is
marginal. Above 20 `k T` the surfaces jam. The metric carries the
`quasi_static` fidelity label.

The model is rigid. The atoms follow the gear motion, and the lattice
does not relax. A relaxed lattice has a lower barrier, so this value is an
upper bound. The interaction uses one carbon-like well for every atom, so
the value compares designs. It is not an absolute friction coefficient.

Consequences. The metric adds about 4 seconds to a score in release mode.
It is the first quasi-static metric, so the scorecard now shows two
fidelity labels. A relaxed variant can use the engine `System` and
`minimize`, and it will lower the bound.

## ADR-0051: A relaxed sweep bounds the slip barrier from below

Status: accepted.

Context. ADR-0050 gives a rigid slip barrier. The atoms follow the gear
path and the lattice never relaxes. That value is an upper bound. The
real barrier depends on how much the lattice gives way at the mesh.

Decision. The `nanocad-meter` crate also gets a `relaxed_slip_barrier`
metric. At every time sample it selects the atoms of the sun and the
first planet that face each other inside a contact radius. The planet
contact atoms are free. Each free atom feels the shifted Lennard-Jones
interaction of the sun and a harmonic tether to its own rigid position.
The tether stiffness is 300 N/m, the value the parameter crate extracts
for a carbon bond. The free atoms relax by damped gradient descent, with
the displacement of one step capped at 4e-12 m.

The reported energy is the interaction energy plus the tether energy, so
the rigid configuration is a feasible point of the relaxation. The
metric also reports the rigid barrier on the same samples and atom sets,
so the two numbers compare directly.

The samples must resolve the atom spacing. One sun tooth pitch of 9 nm
spans about fifty carbon bonds, so the default is 240 samples for the
rigid metric and 120 for the relaxed metric. A coarse sweep of 24 samples
raised the barrier by 70 percent through aliasing. The value at 60 and
240 samples agrees to all printed digits, which shows that 240 samples
resolve the curve.

Consequences. On the default set the relaxed barrier is 5.7371e-20 J
against a rigid 5.7328e-20 J. Relaxation changes the bound by +0.07
percent. The lattice is stiff, so the rigid value is already a tight
bound and the two metrics agree on the verdict. The relaxed metric is a
check on the rigid one, not a replacement. A true relaxed value needs the
real bond topology and a full minimization. It can use the engine
`System` and `minimize` later.

## ADR-0052: A finite-difference Hessian gives the mesh normal modes

Status: accepted.

Context. The scorecard has geometric and quasi-static metrics. The next
fidelity rung is `Harmonic`. No Hessian code existed in the workspace.

Decision. The `nanocad-engine` crate gets a `hessian` module. The Hessian
is a central finite difference of the analytic gradient, and a Jacobi
rotation sweep gives the symmetric eigenvalues. The finite-difference
step must be much smaller than a bond: 1e-8 m crosses the bond minimum
and gives a 1.5 percent error, so the tests use 1e-12 m and the meter
uses 1e-13 m.

The `nanocad-meter` crate gets a `mesh_mode` metric with `Fidelity::Harmonic`.
It cuts a cluster of up to 180 atoms from the sun and the first planet
around their closest approach, builds a `System` with a 300 N/m bond
stretch on the real bonds and a Buckingham van der Waals term, relaxes
the cluster, and reports the softest internal mode.

The exclusion rules carry the physics. A bonded carbon pair sits at
1.54e-10 m, far below the van der Waals minimum, so the term repels it;
without exclusions the relaxation fails. The long-range van der Waals
tail has negative curvature and makes the cluster unstable, so the
contact shell is limited to 3.5e-10 m. Above that radius the tail is
excluded. With both rules the cluster has zero unstable modes.

Consequences. On the default set the softest mesh mode is 1.184e13 Hz
(395.1 per cm) over 148 atoms and 157 bonds, with zero unstable modes.
The mode is far above the thermal frequency, so the mesh is stiff. This
agrees with the slip-barrier verdict.

The cost is one relaxation and one 444x444 Hessian per score, about nine
seconds in release, in line with the other metrics. The cluster is a
local model, not the whole gear, so the value is a local mode. Six rigid
modes are always zero, and a cluster atom with no bond and no contact
adds three more, so the report skips every frequency below 1e9 Hz.

## ADR-0053

### A CMA-ES search over the tooth counts lowers the slip barrier

Status: accepted.

The design loop needs a search stage. The score already measures a design,
so the search only needs to propose designs and read the score. We add a
new crate, `nanocad-opt`, with a self-contained covariance matrix
adaptation evolution strategy (CMA-ES).

The strategy is deterministic. The workspace has no random crate, and a
reproducible run matters more than a fast one here. A SplitMix64 generator
seeds the Box-Muller normal draws from a stated seed, so one seed gives one
result. The population is `max(4, population)`, and the search stops early
when it finds no better design for eight generations.

The objective is the slip barrier in units of the thermal energy. An
infeasible design has no gradient, so the objective adds a penalty for the
clearance shortfall, scaled by 100 kT for every 1.0e-10 m. Without the
penalty 18 of 24 candidates were rejected and the search had no signal.

The search covers three integer parameters: the sun teeth (9 to 16), the
planet teeth (7 to 12) and the planet count (2 to 4). The module and the
layer count stay at their defaults, so the designs compare at one size.

Consequences. On the default set the search starts at 13.84 kT (12 sun
teeth, 9 planet teeth, 3 planets) and reaches 8.92 kT (12 sun teeth, 10
planet teeth, 4 planets) in 80 evaluations. The barrier falls by 36
percent. The clearance of the best design is 2.8865e-10 m, which passes
the 2.52e-10 m target. The best design has 178073 atoms, so it is larger
than the default.

The cost is about two seconds for each evaluation in release, because each
candidate builds a full gear set and runs the slip sweep. The app route
`/api/optimize` caps the population at 16 and the generations at 12. The
search is a local improvement, not a global optimum. The result is a
simulated estimate, and it is not a validated design.

## ADR-0054

### The bonded terms are built as reusable helpers, but the mesh Hessian keeps its cut-out limit

Status: accepted.

The meter potential needs every bonded term that the engine can provide. The
engine already has a torsion term and an out-of-plane term. The meter built
bonds and van der Waals only. The new module `crates/meter/src/bonded.rs` adds
both missing terms as reusable builders.

The torsion builder takes the neighbour list of every atom. For each central
bond `j->k` it pairs every `i` on `j` with every `l` on `k`, and adds one
torsion for each path. The amplitude is `V3` only, with
`TORSION_V3_J = 2.0e-20` J, which is about 12 kJ per mole. The diamond lattice
is staggered, so this term sits at its minimum in the ideal crystal. A test
confirms that the generated geometry has a torsion energy of exactly zero.

The out-of-plane builder adds an improper term at a three-coordinate centre.
The amplitude is `IMPROPER_K_J_PER_RAD2 = 2.0e-18` J per radian squared. The
builder reads the current positions and adds the term only when the measured
out-of-plane angle is within `TRIGONAL_TOLERANCE_RAD = 0.5` rad. A tetrahedral
centre sits at about 0.616 rad, so the gate keeps it out. The gate exists
because a cut-out cluster leaves boundary carbons with three neighbours. Their
geometry is still tetrahedral. Forcing them planar made the minimizer fail with
`LineSearchFailed`.

**The two terms are not wired into the `mesh_mode` Hessian.** The investigation
found the reason. The engine torsion gradient is correct: a test relaxes a
small torsion system to a stationary point and then finds no unstable mode. But
the `mesh_mode` cluster is a cut-out of a larger gear. Its 536 torsions add
about `-8e-6` N/m of negative curvature. That value overlaps the genuine soft
mode band, because a 1e9 Hz mode is `4e19` per second squared and a 1e10 Hz
mode is `3.9e21`. No threshold separates the artifact from a real soft mode.
With the torsion on, the softest mode fell from `1.184e13` Hz with 0 unstable
modes to `1.094e9` Hz with 21 unstable modes. The Hessian of a whole part has no
cut-out boundary, so it is the correct host for these terms. Task M11-04 does
that work.

The result is a simulated estimate of a model potential. The amplitudes are
representative values, not fitted parameters, and they are not validated
against an experiment.

## ADR-0055

### A loaded contact reports the friction force as a screening estimate

Status: accepted.

The slip barrier gives no normal load. A real tooth pair carries a load, so the
tangential force under that load is the useful number.

The metric is `contact_friction` in `crates/meter/src/loaded.rs`. It presses the
sun and the first planet together along the contact normal. A bisection finds
the offset that makes the normal force equal the stated load. It then slides the
pair through one sun tooth pitch. The lateral force is the fall of the
interaction energy over the arc length. The mean of the absolute value is the
friction force, and the peak is the static value. The wear count is the number
of unique sun-planet pairs that come within the wear distance at any sample.

The bodies are rigid, so the model has no elastic contact area and no plastic
flow. The normal load is a small model value, so the friction coefficient can
exceed one. That is a property of the model, not of a real material. The metric
is a screening estimate. It compares designs and it does not certify a friction
coefficient.

The default set gives a mean friction force of 1.677e-10 N at a 1.0e-11 N load
over a 60-step sweep. Ten pairs come within 3.0e-10 m during the slide.

## ADR-0056

### A steered drive reports the torque on a held joint

Status: accepted.

A design is not proven by a static score. The question is what torque the drive
must supply, and how much energy the contact loses.

The metric is `steered_drive` in `crates/meter/src/drive.rs`, at the dynamics
fidelity. It builds a contact cluster around the sun-planet mesh. It holds the
planet atoms at their relaxed positions. At every velocity-Verlet step it writes
the sun velocity from the stated drive rate. The resisting torque is the
negative z component of the torque on the sun atoms. The report holds the mean
and peak torque, the drive work, the energy loss, the free-atom count, the final
temperature, and the temperature rise.

The bodies are rigid and the cluster is a cut-out, so the value is a screening
estimate. The cluster rejects every atom that holds a bond across its own
boundary; a free atom is one that only moves through the non-bonded term.

On the default set the cluster holds 66 atoms and the mean torque is
1.99e-45 N m. That is effectively zero: the non-bonded contact at the relaxed
mesh is nearly torque-free. The ring sits about 1.0e-8 m away, so it never
enters the 6.0e-10 m collection radius, and the cluster has no free atom. The
temperature rise is therefore not measurable, and the note says so.

The metric is not in the app scorecard. A dynamics metric with 4000 steps is
too slow for the score route, and a near-zero torque is not a useful scorecard
row. It stays available as a metric and as
`crates/meter/examples/drive_json.rs`.

## ADR-0057

### A whole part carries the bonded terms in its Hessian

Status: accepted.

ADR-0054 kept the torsion and out-of-plane terms out of the mesh Hessian. The
reason was the cut-out cluster: its dangling crystal boundary turned the bonded
terms into spurious negative curvature. A complete part has no such boundary, so
the terms are physical there.

`HarmonicMesh::measure_part` in `crates/meter/src/harmonic.rs` measures a whole
part. It returns the same `mesh_mode` metric at the same harmonic fidelity. It
adds the torsion and out-of-plane terms from `crates/meter/src/bonded.rs`,
because a whole part has no dangling boundary. It excludes every pair inside the
`COVALENT_EXCLUSION_M = 3.1e-10` m covalent shell from the van der Waals term.
Without that exclusion the minimizer fails with `CoincidentNonbondedAtoms`: a
part is one covalent network, so its second and third shells must not act
through the non-bonded term. It rebuilds the neighbour list after the relaxation
and before the finite-difference Hessian, because the relaxation moves the
atoms.

A part above `target.max_atoms` is refused with an empty report, because the
Jacobi diagonalisation is cubic in the atom count.

The smallest plate, 477 atoms and 612 bonds, gives a softest internal mode of
5.733e10 Hz with 1 unstable mode in a release run of about 215 s. The cut-out
cluster of the same mesh gave 21 unstable modes under the same terms. This
confirms the ADR-0054 diagnosis: the boundary, not the terms, caused the
artifact.

## ADR-0058

### A small selection language chooses nodes

Status: accepted.

A design grows past the size where a caller can name nodes by index. A short
query text is the practical interface, and SAMSON shows the value of the idea
with its Node Specification Language.

The language lives in `crates/model/src/selection.rs`. It is a boolean
expression over predicates: `and`, `or`, `not`, parentheses, `all`, `none`,
`atom`, `bond`, `part`, the four element names, and the comparisons `charge`,
`degree`, `atom_index`, `order`, `bond_length`, `bond_index`, `part_index`,
`name == "..."` and `within <d> of atom [i]`. The keywords are
case-insensitive.

The language is a stated subset of the SAMSON language. It has no arithmetic,
no variables, no units, and no nested node sets. `within` works inside one part
only, because an atom index is part-local along with the rest of the model.

Every expression evaluates to one set per node kind. `and`, `or` and `not` act
per kind, so `carbon or bond` is well defined. A predicate that produces only
atoms leaves the bond and part sets empty.

The parser is a hand-written recursive descent. It returns
`ModelError::SelectionParse { position, message }` and it does not panic.

## ADR-0059

### The part library is one table and placement is a port frame

Status: accepted.

M10 added nine generators, and M10-05 to M10-08 left them out of the registry.
The registry was gear-only, so nothing could list or dispatch the full set.

`crates/parts/src/registry.rs` now holds one static table of 18 generators,
each with a category: `gears`, `lattice`, `structure` or `device`. The
`LibraryEntry` list is built from that table by calling `id()` and `name()` on
each generator, so a name cannot drift from its generator. The gear-only
functions stay, because the gear tests and the gear callers depend on them.

Placement is a small piece of the connector idea from M10-03.
`crates/parts/src/placement.rs` resolves a plug port and a socket port into a
`PortFrame`, then moves a whole part by that frame: a rotation about z through
the origin, then the frame offset. `assembly_document` collects the placed
parts into one document.

The frame is planar, because the port model is planar. A placement with a
three-dimensional rotation needs a full transform, and that is future work.

## ADR-0060

### The design document is a plain-data snapshot history

Status: accepted.

A design has parameters, parts and measurements. The user needs to walk the
history back and forward. The parts crate holds the generators, and a saved
document must load without them.

`crates/jigs/src/design.rs` holds the document. A `DesignSnapshot` carries
named parameters, part records and measurements. Every field is plain data, so
serde reads and writes a document with no generator present. `PartRecord`
reads a part schema, and `Measurement` reads a metric value, but neither value
holds a live part.

`DesignDocument` holds one array of snapshots and one cursor. Undo and redo
move the cursor, so the two operations share one code path. `MAX_HISTORY` is
64. A commit equal to the current snapshot changes nothing and returns false.
A new commit drops every snapshot after the cursor, because the redo tail is
no longer reachable.

The app keeps the same rule in `site/viewer.js`, with the same limit. Each
build commits a parameter snapshot, and the Undo and Redo buttons move the
cursor and rebuild the scene.

The document holds parameters, parts and measurements. It does not hold a
motion, a force field, or a camera. Those belong to a later document version.

## ADR-0061

### The sorting rotor is built as two diamond parts that mate on the axis

Status: accepted.

The sorting rotor is the device in Freitas, *Nanomedicine* Volume I, Section
3.4.2, after Drexler. A disk carries binding pockets along its rim, and each
pocket carries a bound molecule from the outside solution to an inner chamber
as the disk turns.

Two generators build the mechanics of the device. `SortingRotorGenerator` cuts
a disk with a `pocket_count` of cylindrical pockets on the rim and a central
bore. A pocket opens to the rim when its centre circle plus its radius reaches
past the outside radius, because that is the face the solution sees.
`RotorHousingGenerator` cuts an annular slab with a chamber for the disk and
two radial channels, one inlet and one outlet. Each channel runs from inside
the chamber wall to past the outside radius, so it opens the chamber instead of
forming a trapped pocket.

Both parts are cut from the diamond cubic lattice, so every bond is the diamond
first-shell length and every carbon keeps a valence of four. Both are registered
in the part library under the `device` category.

The rotor mates to the housing through the port model. The rotor `axis_port` is
revolute on `z` and the housing `chamber_port` is fixed on `z`, so the two parts
share one axis and the disk is free to turn.

Scope limit. The generators build the mechanics and the example computes the
geometry, the mass and the kinematics. The model does not include a binding
site, a target molecule, a solvent or an ion. It therefore does not simulate
selectivity: it does not show that a pocket binds one species and rejects
another. That needs a chemistry-specific force field and it is a separate
milestone. The drag power of 1e-16 W in the reference is likewise not
simulated; the example reports it only as a reference figure.

The measured agreement is close. At the default sizes the rotor is 53964 atoms
and 39848 atoms of housing, so the assembly is 93812 atoms and 1.46e-21 kg,
against the reference of about 1e5 atoms and 2e-21 kg. At a stated rate of
86000 rev/s the rim turns at 3.78e-3 m/s, against the reference 2.7e-3 m/s, and
twelve pockets pass 1.03e6 molecules per second, against the reference rate of
about 1e6.

## ADR-0062

### One scene schema carries both mechanisms

Status: accepted.

The sorting rotor needs a viewer, and the planetary gearbox already has one.
Two schemas would double the viewer and the exporter, so the rotor scene uses
the gearbox scene schema, and the schema states the mechanism in its data.

The housing is body 0 and it is fixed. The rotor is body 1 and it is free. The
device layer carries one revolute joint on the world `z` axis through the
origin, because the rotor axis port and the housing chamber port share that
axis. The gear couplings and the gear constraints are empty.

The design block describes a gear set, and a rotor has none. The rotor scene
therefore writes a zero module. A viewer reads a zero module as "not a
gearbox", and it then drives the bodies from the joint instead of from the gear
kinematics. This is a convention on an existing field, so the schema version
stays at 1. The alternative, a nullable design block or a new mechanism field,
would change the schema and every reader for one extra mechanism.

Two consequences follow from the reuse. The viewer must not score a rotor,
because the scorecard measures the gearbox, so it clears the scorecard and
states the reason. The viewer must also not show the real turn rate. The real
rate is 86000 revolutions per second, and no display shows that, so the viewer
turns the rotor at 2.0 rad/s and the readout names it a display rate.

The scene holds the geometry, the mass and the joint. It does not hold a
binding site, a target molecule, a solvent or an ion, so it does not show
molecular selectivity. That limit is unchanged from ADR-0061.

## ADR-0063

### The chemistry layer states its data and never invents it

Status: accepted.

The chemistry layer needs numbers that a real force field would supply: a
valence, an atomic mass, a contact distance and a partial charge. This
repository has no such table, and a fitted parameter set is only worth
something with a citation.

The rule is therefore: state the data, cite the source, and never invent a
number.

The element table in `crates/model/src/chemistry.rs` holds a valence, an
atomic mass and a covalent radius for nine elements. The valence is the usual
valence of the neutral element. The mass is the IUPAC standard atomic weight.
The radius is the Cordero 2008 single-bond radius. Each value is a published
fact and not a fitted constant. An element outside the table returns `None`
instead of a silent default, so a caller cannot use the table by accident.

The functional groups in `crates/parts/src/group.rs` take their bond lengths
from that table, so a group and the diamond lattice share one source of radii.
The remaining bonds take the tetrahedral angle, which is geometry and not a
parameter.

The guest molecules in `crates/parts/src/guest_data.rs` are frozen data. The
geometry comes from RDKit ETKDGv3 embedding and MMFF94 optimization. The
partial charges are Gasteiger-Marsili PEOE values. RDKit is BSD-3 licensed,
so the values are free to use, and the emitter script states the seed and the
iteration count, so anyone can reproduce the table.

This repository still has no nonbonded parameter set. `VanDerWaalsTerm` takes
raw per-atom `a`, `b` and `c` values and nothing maps an element to them. A
later task must add that table, and it must cite its source.

Consequence: every result that uses a charge from this layer is a model
result. The charges are not computed for the pocket, and they are not
validated against experiment. A report must say so.
