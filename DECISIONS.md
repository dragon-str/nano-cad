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
