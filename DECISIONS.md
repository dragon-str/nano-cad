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

## ADR-0064

### The nonbonded parameters are UFF, and the term is Lennard-Jones

Status: accepted.

ADR-0063 refused to invent a nonbonded parameter set. This decision supplies
one, so that a pocket can be sized from contact distances and a binding metric
can have a repulsion term.

Source. The values are the nonbonded parameters of the Universal Force Field,
from Rappe, Casewit, Colwell, Goddard and Skiff, *Journal of the American
Chemical Society* 1992, 114, 10024. We transcribed them through the RDKit
implementation, which is BSD-3 licensed. The transcription is reproducible:
the table is at `Code/ForceField/UFF/Params.cpp` in the RDKit repository.

The reduction to one row per element. UFF has one row for each atom type, and
it has several types for one element, for example `C_3`, `C_2` and `C_R`. For
every element in our table, all of its types share the same `x1` and `D1`, so
the table reduces to one row for each element. We checked this claim for all
eight elements rather than assuming it.

The form. UFF states its nonbonded term as
`U(r) = D_ij ((x_ij/r)^12 - 2 (x_ij/r)^6)` with `x_ij = sqrt(x_i x_j)` and
`D_ij = sqrt(D_i D_j)`. That is a Lennard-Jones 12-6 curve with
`epsilon = D_ij / 4` and `sigma = x_ij / 2^(1/6)`. We add a Lennard-Jones term
to the engine, in `crates/engine/src/lennard_jones.rs`, rather than fitting a
Buckingham curve to the same values. A fit would be a model of a model, and
the engine would then disagree with the source it cites.

The term carries its own finite-difference gradient test, as the ground rules
of this repository require for every analytic force term.

The limits. This table gives a size and a well depth. It does not give a
partial charge, a torsional barrier, an angle force constant or a hydrogen
bond. The charges remain the stated model charges of ADR-0063. The agreement
of a computed binding energy with experiment therefore follows from the whole
model, not from this table, and no result in this repository is validated
against experiment.

UFF is also known to be weak for metals and for inorganic compounds, and most
of its entries were never validated. Every result that uses this table is a
model result for that reason as well.

## ADR-0065

### The binding pocket is a decorated well in the diamond lattice

Status: accepted.

A pocket that sorts one molecule from another needs two properties. It needs a
cavity that fits one guest, and it needs a wall that holds one guest better
than another. `BindingPocketGenerator` in `crates/parts/src/pocket.rs` builds
both.

The cavity is a cylinder subtracted from a diamond block. The well opens at
the top face, so a guest can enter from the solution and leave to the chamber.
The block is large enough to hold a wall around the well, so the well is a hole
and not a through slot.

The wall chemistry is decoration, not geometry. `diamond_solid` already splits
a fill into bonding, free directions and capping. This generator reads the free
directions, and for each one it steps a short distance along the bond. When
that probe point lands inside the cavity, the generator attaches the chosen
`FunctionalGroup` there and lets the group's own hydrogens finish the valence.
Everywhere else it caps with hydrogen, so the surface is unchanged.

The probe rule is geometric on purpose. A direction that points into the well
is a direction that faces the guest, and the wall atoms that face the guest are
exactly the atoms that can hold it. No separate surface-area calculation is
needed, and the bond lengths and angles come from the chemistry table.

The group is a parameter index, not a name. A `ParameterSet` holds numbers, so
`wall_group` is an index into `FUNCTIONAL_GROUPS`. This keeps the generator on
the same interface as every other generator.

Scope limit. The pocket selects by shape, through the size of the well, and by
stated charge and dispersion, through the wall group. It does not model a
solvent, an ion, a hydrogen bond with an explicit geometry, or a
conformational change. It does not prove that a pocket binds one species in a
real solution. It measures the interaction of a frozen guest with a frozen
wall. Every number it gives is a model number.

Three limits follow from the lattice and must be read with every result.

1. The lattice is discrete. The well radius moves in lattice steps, so two
   radii a tenth of an Angstrom apart can give the same part. Any claim about
   a pocket that turns on a finer distinction is not supported.
2. The decorated groups reach into the nominal cavity. The effective opening
   is smaller than `pocket_radius_m`, so a guest that fits by the stated
   radius can still fail in a built part.
3. Size alone cannot separate two guests of the same size class. Ethanol and
   dimethyl ether share the formula C2H6O and their extents differ by one
   tenth of an Angstrom. The well that admits one admits the other, so only
   the wall chemistry can separate them. That is the measurement of M14-06,
   and it is the measurement that gives this layer its meaning.

## ADR-0066

### The binding metric is a rigid search in a frozen pocket

Status: accepted.

The pocket generator of ADR-0065 builds a cavity with a stated wall
chemistry. This decision adds the measurement that gives the wall chemistry a
meaning.

`crates/meter/src/binding.rs` states one model. The pocket is rigid. The guest
is rigid. The interaction is the Lennard-Jones sum from ADR-0064 plus the
Coulomb sum, over every atom pair inside a cutoff. The charges come from the
Gasteiger-Marsili model, stored as data by ADR-0063, and they are written onto
the wall by `charge_wall`, because a built part starts neutral.

The search moves the guest along the pocket axis and turns it about that axis.
The lowest energy of the search is the binding energy, and the ratio of two
binding energies is the selectivity.

A position where any pair comes closer than `overlap_ratio * sigma` is not a
bound state, and the search rejects it. This guard is necessary. Without it the
search finds interpenetrating positions and reports an energy that no real
contact can give. `overlap_ratio` is 0.85, and it is a stated parameter of the
model.

Three limits must be read with every result.

1. The search is coarse. It slides and turns the guest, and it does not relax
   the guest or the wall. The energy is an upper bound on a relaxed binding
   energy, so the value is comparable only within one target and one guest.
2. The model has no solvent, no ion and no explicit hydrogen bond. A hydrogen
   bond appears only as the sum of its Coulomb and dispersion parts, so the
   model cannot state a hydrogen-bond energy.
3. The pocket must be larger than a first estimate suggests. A bound guest
   needs the well radius to exceed its extent plus one contact distance, and
   the decorated groups themselves intrude into the cavity. A well of
   radius 3.6e-10 m rejects every position, and a well of radius 7.0e-10 m
   accepts every one.

The result of the measurement at a well radius of 7.0e-10 m, in kT at 300 K,
is the energy of each isomer, and then the ratio of ethanol to dimethyl ether:

| wall group | ethanol | dimethyl ether | ratio |
|---|---|---|---|
| hydroxyl | 26.37 | 22.21 | 1.188 |
| amino | 19.23 | 19.39 | 0.992 |
| methyl | 9.50 | 10.93 | 0.869 |
| fluoro | 6.44 | 10.77 | 0.598 |
| chloro | 6.05 | 6.26 | 0.967 |
| thiol | 18.47 | 9.34 | 1.978 |

The two isomers share the formula C2H6O, and they differ in shape by one tenth
of an Angstrom and in charge. The shape difference alone cannot separate them,
as ADR-0065 states. The table shows that the wall chemistry can. The choice
moves by a factor of 3.3 across the six groups, and it changes sign: the
hydroxyl wall takes ethanol, the hydrogen-bond donor, and the fluoro wall takes
dimethyl ether.

This is a model statement. It is not a claim that a real pocket separates
ethanol from dimethyl ether in water.

A unit error was found while this metric was built, and the fix matters for
every earlier result. Gasteiger returns a charge in units of the elementary
charge. `Atom::charge_c`, by contract, holds coulombs. The first data file
stored the raw Gasteiger number, which made every charge 1.6e19 times too
large and every Coulomb energy about 1e17 J. The emitters now multiply by
`ELEMENTARY_CHARGE_C`. The table above holds only after that fix.

## ADR-0067

### The ejection rod is a stepped column, and its work is a barrier

Status: accepted.

The reference states that the bound molecules "are forcibly ejected by rods
thrust outward by the cam surface". `EjectionRodGenerator` in
`crates/parts/src/rod.rs` builds the rod. The part is a stepped diamond
column: a shaft that runs in a bore, and a narrower tip that enters the pocket
and touches the guest. The step is the point of the shape, because the tip must
pass the pocket mouth while the shaft stays in its bore.

The rod carries two ports. `cam_port` is prismatic and it sits at the shaft
root, pointing back along the shaft. A prismatic joint is the honest model of a
rod in a bore: the rod has one freedom and no more. `tip_port` is fixed at the
tip, so a pocket or a guest can be attached to it.

The scope limit is the cam. The generator does not build a cam surface, and it
does not model the sliding fit, the shaft friction, or the force that the cam
applies. It builds the rod and states the direction in which the rod moves.

The work metric is in `crates/meter/src/ejection.rs`. The guest starts at the
bound position and the bound orientation that the binding metric found. The
rod then pushes the guest along the pocket axis in fixed steps. The interaction
energy rises while the guest passes the wall and falls when the guest leaves
the well. The work that the rod must supply is the height of the barrier above
the bound state. The peak sample and the last sample are both reported, so a
reader can see whether the guest really left.

Both metrics now share one kernel. `pocket_wall` returns the wall atoms that
face a guest, and `interaction_energy_j` is crate-visible. A second copy of the
pair loop would drift from the first.

The honest comparison. The reference gives 10 to 40 zJ per molecule at 310 K.
That figure is the reversible work at a stated concentration ratio. A
reversible process returns the binding energy, so the reference figure is a
lower bound and it is not the raw binding energy. The model does not recover
any energy, so a model work above the bound is expected and it is not a
contradiction. At a well radius of 7.0e-10 m the model work runs from 1.9e-20 J
(19 zJ, chloro wall, dimethyl ether) to 2.8e-19 J (283 zJ, hydroxyl wall,
methanol), and the best case is 1.9 times the reference bound.

The selectivity survives the change of metric. The wall group changes the work
by a factor of 14. The thiol wall takes 1.46e-19 J to eject ethanol and only
3.70e-20 J to eject dimethyl ether, while the chloro wall takes 1.94e-20 J to
eject dimethyl ether. The rod therefore sees the same preference that the
binding metric sees, which is the point of measuring both.

Limits on this metric. The pocket and the guest stay rigid. The search is
coarse, so the work is an upper bound for a given target. There is no solvent,
no ion and no explicit hydrogen bond. The rod is not simulated, so the metric
does not include the work that the rod itself returns, the work against the
bore, or any elastic energy in the rod. Every value is a model value.

## ADR-0068

### The rotor scene holds the whole device, and the panel measures the wall

Status: accepted.

M14 built a binding pocket, a rod and two metrics, but no scene joined them
to the rotor. A reader could see the rotor turn and could read a selectivity
table, and nothing tied the two together. M15 joins them.

The rod is body 2 of the rotor scene. `build_rotor_scene` in
`crates/jigs/src/rotor_scene.rs` builds the rod with `EjectionRodGenerator`,
places its tip at the centre of pocket 0 with `place` on the `tip_port`, and
declares a prismatic joint whose parent is the housing and whose axis is
radial. The choice of parent matters. The rod presses on the pocket wall, so
the housing must carry the reaction and not the rotor. The axis is radial,
because that is the direction that pushes a bound guest out of a pocket.

A prismatic joint is a translation without rotation. The rod therefore slides
along one line and it does not turn with the rotor. This is the behaviour that
the reference describes: rods thrust outward by a cam surface. The scene does
not build the cam. It builds the joint that a cam would drive, and the viewer
shows the travel that the joint permits.

The example and the panel report the whole device. `rotor_json` prints the
rod atom count, the rod mass and the totals of all three parts, so the facts
that the panel shows cover the housing, the rotor and the rod. Measured: the
rod holds 328 atoms and 3.6185271751219556e-24 kg, and the device holds 95444
atoms and 1.4668573965626264e-21 kg. The rotor holds 55268 atoms after the
pocket decoration, and the housing holds 39848.

The viewer stroke is a display convention. `ROTOR_ROD_STROKE_M` in
`site/viewer.js` is 1.5e-9 m, and it matches `EjectionTarget::travel_m` in
`crates/meter/src/ejection.rs`, so the picture and the metric agree on the
same distance. The rotor rate in the viewer is 2.0 rad/s and the real rate is
86000 rev/s, so the viewer shows one stroke for each turn and hides the rate.
The panel says so.

The panel reads the selectivity table from the metric and not from a copy.
`GET /api/selectivity` runs `crates/meter/examples/binding_json.rs` through
`_run_example` and returns its JSON, so one Rust program is the only source of
the numbers. A Python table would drift from the model. The route takes an
optional `radius_m`, so a caller can ask for another well radius.

The table is a model result. It states binding energy as a multiple of kT for
the ethanol and dimethyl-ether isomer pair, and it gives the ratio. The search
holds the pocket and the guest rigid, the model has no solvent and no ion, and
there is no explicit hydrogen bond. The panel prints that limit under the
table. The measured ratios at a well radius of 7.0e-10 m are hydroxyl 1.188,
amino 0.992, methyl 0.869, fluoro 0.598, chloro 0.967 and thiol 1.978. The
wall chemistry changes which isomer binds more strongly and it changes the
sign, so the ratio is a real output of the force model. It is still not a
prediction of separation in water.

Limits on this work. The scene has no cam surface, no bore, no shaft friction
and no elastic energy in the rod, so it shows the motion that the joint allows
and not the force that drives it. The rod is one of twelve, because the scene
places one rod at pocket 0 and the other eleven pockets have no rod. The
panel shows one more device than M14 did and it still does not show a guest, a
solvent or a molecule in a pocket.

## ADR-0069

### The cam ring stays in the rotor plane, and the rotor bore clears it

Status: accepted.

M15 placed one ejection rod in the scene, and the picture was wrong. The rod
stood vertically, in free space below the disk, because `place` rotates about
z only and the rod axis has no xy part. M16 houses the rods where they belong
and drives them, which is what the reference describes: "the bound molecules
are forcibly ejected by rods thrust outward by the cam surface."

A rod is a long thin column, so it must lie along the radius of its pocket and
slide along that radius. That poses a problem, because the rod axis is local
`+z`, and a `Placed` solid and the `place` mate can rotate about z only. The
scene therefore does not use `place` for the rod. It maps the rod atoms
directly. The map sends local `+z` to the radial direction at the pocket
azimuth, and it sends local `x` to world `-z`. The rod is symmetric about
local `x`, so the picture is the same either way. The tip lands on the inner
wall of the pocket and the inboard end lands on the base of the cam, so the
geometry states the stroke and the mount together.

The cam ring is a fixed flat ring, and its lobe is a radial key that rises
once around the ring. It must thrust a rod outward, so the lobe is taller than
the ring base and it must move into the region that the rods occupy. A rotating
slot cannot always hold a fixed key. The rotor therefore has a central bore of
3.5e-9 m, which is exactly the peak radius of the lobe. The whole region inside
that radius is open, so the fixed lobe never touches rotor material, and the
whole mechanism stays in the z plane. No follower pin leaves the plane, and no
face groove is needed. This is the reason for the change of the `bore_radius_m`
default, and it is why the rotor is now an annulus.

Each rod sits in a bore of its own, and each bore runs from the central bore
out to the inner wall of its pocket. The rod is captive in the rotor, so its
prismatic joint parent is the rotor and not the housing. This supersedes the
one-rod decision of ADR-0068, which made the housing the parent. The rotor
carries the rods around, and the fixed lobe pushes one rod at a time as that
pocket crosses the lobe.

The rod is 4.0e-9 m long and the cam stroke is 2.0e-9 m, which is the full
diameter of a pocket, so an extended rod sweeps the pocket and reaches the rim.
The rod count is 12, one for each pocket. The measured device holds 91553
atoms: the rotor 39368, the housing 39848, the cam 3793, and 712 atoms for
each of the 12 rods. The mass is 1.290735679733618e-21 kg, against the
reference 2e-21 kg for a rotor with its housing and a share of the drive.

The viewer follows the same model as one Rust program. It reads the prismatic
joints from the scene and gives each rod a smooth phase ratio against the
fixed lobe, so one rod extends as its pocket crosses the lobe and the others
stay home. The readout states the body count and the joint count. The
headless harness measured the radius gain of the extended rod as 1.747299e-9 m
against a model gain of 1.751806e-9 m, and it found exactly one extended rod
of twelve.

Limits on this work. The cam pushes outward and nothing pulls a rod back, so
the model omits the return spring; the panel says that the profile is stated
and the spring is not built. The lobe is a radial key and not a tuned motion
law, so the rod acceleration is not designed. The rod generator still builds
no follower pin and no sliding fit, so the model shows the motion that the
joint allows and not the contact between a pin and a profile. The work of
ejection remains the measured barrier in `crates/meter/src/ejection.rs` and
not a force on a simulated rod. The scene still holds no guest, no solvent and
no molecule in a pocket.

## ADR-0070

### The cam leaves the rotor plane, and a keyed drive shaft turns the rotor

Status: accepted. This supersedes ADR-0069.

ADR-0069 kept the cam in the rotor plane and opened the rotor bore to the peak
of the lobe. The user reported that the result does not make mechanical sense:
the cam ring cannot be fixed to the housing where it overlaps the turning rods,
and the ring should be held from above and below. The user also asked how the
rotor is powered and recalled an off-centre cam that pushes and pulls the rods
from the rotation alone. That memory is correct. This decision moves the cam out
of the rotor plane.

The cam is now a plate below the rotor. The rotor carries the 12 rods, one for
each pocket, and each rod carries a follower pin that hangs into a groove in the
plate face. The groove is a circle whose centre is offset from the rotor axis by
1.0e-9 m. The groove centreline radius therefore runs from 2.0e-9 m to 4.0e-9 m,
a stroke of 2.0e-9 m, and it returns with no spring: one wall of the groove
pushes a rod out as the offset grows and the opposite wall pulls it back as the
offset falls. The cam plate is fixed, so it cannot be struck by a rod that
turns; the rods turn above it, and only the pins enter the groove.

The cam plate rests under the housing ring. Its top face is flush with the
housing lower face, and its outer radius equals the housing outside radius, so
the plate laps under the ring and the housing holds it from below and around its
rim. The bore of the plate (1.5e-9 m) clears the drive shaft. Because the cam no
longer enters the rotor plane, the rotor central bore returns to 1.5e-9 m from
the 3.5e-9 m of ADR-0069, and the rotor can carry a hub again.

The rotor is powered through a drive shaft on the rotor axis. A key on the shaft
sits in a keyway that `SortingRotorGenerator` cuts into the rotor bore. The key
locks the shaft to the rotor rigidly, so the shaft and the rotor turn as one
body. The scene therefore puts the shaft atoms on the rotor body and gives the
shaft no joint of its own. The scene schema has no rigid coupling, and a
revolute joint between the shaft and the rotor would model a bearing, not a key.
The machine that turns the shaft is out of scope, because the user set it aside
for a later design; the panel says that the shaft is keyed to the rotor and that
the external drive is not shown.

The rod's prismatic joint still names the rotor as its first body, because the
rod remains captive in a radial ejection bore in the rotor. The rod atoms are
still placed by a direct atom map, because `place` rotates about z only and
cannot tilt the rod into the radial direction. The pin atoms ride the rod body,
so they translate with the rod and stay in the groove.

The measured scene now holds 15 bodies and 13 joints and 157637 atoms. The rotor
body holds 54844 atoms (49708 rotor and 5136 shaft), the housing 39848, the cam
plate 54953, and each rod 666 (602 rod and 64 pin). The mass is
2.2431991352736912e-21 kg. The rod is 3.5e-9 m long, and its tip runs from
5.5e-9 m to 7.5e-9 m, so an extended rod sweeps its pocket to the rim.

The viewer follows the same model. It drives every rod with the eccentric
formula `0.5 * (1 + cos(delta))`, where `delta` is the rod azimuth less the
groove centre angle. All 12 rods move smoothly; the model no longer needs a
one-lobe window. The headless harness measured the counts above, found 12 of 12
rods moving, and measured the worst deviation from the eccentric model as
3.9e-12 m, with the cam plate at rest.

Limits on this work. The pin slides in the groove; a roller follower would roll.
The eccentric circle is a radial key and not a tuned motion law, so the rod
acceleration is not designed and a profiled groove could add a dwell. The scene
models no contact force between the pin and the groove wall, no torque on the
drive shaft, and no external machine. The work of one ejection stays the
measured barrier in `crates/meter/src/ejection.rs`, not a force on a simulated
rod. The scene still holds no guest, no solvent and no molecule in a pocket.

## ADR-0071

### The cam groove dwells, and the cutaway shows the mechanism

Status: accepted. This supersedes ADR-0070.

ADR-0070 moved the cam below the rotor, but the groove stayed a circle. A
circular groove gives a sinusoid, so a rod moved for the whole turn. The user
reported three defects. First, an extended rod overlapped the housing atoms.
Second, a rod pushed the guest out too early, before its pocket reached the
outlet. Third, the mechanism below the housing was not visible.

The overlap came from the rod length. The groove centreline ran from 2.0e-9 m to
4.0e-9 m, and the rod was 3.5e-9 m long, so a rod tip reached 7.5e-9 m. The
housing chamber radius is 7.2e-9 m, so the tip crossed the housing ring at the
outlet. The groove rise is now 1.5e-9 m, so the centreline runs from 2.0e-9 m to
3.5e-9 m and a rod tip stops at 7.0e-9 m, which is the rotor rim.
`build_rotor_scene` refuses a design whose fully extended tip passes the rim.

The groove is now profiled. Its centreline dwells at `groove_base_radius_m`, and
it rises by `groove_rise_m` over a raised cosine of half-angle
`groove_ramp_half_angle_rad` about the azimuth `groove_angle_rad`. The default
ramp is 0.21 rad, about 12 degrees, which is close to the housing outlet window
of about 8 degrees at the chamber radius. A rod therefore pushes only as its
pocket meets the outlet, and it stays retracted for the rest of the turn.
`CamPlateGenerator` sweeps 240 disks of the groove width along the centreline
into one slot, so the groove keeps its width through the ramp. The eccentric
circle of ADR-0070 is gone.

The viewer uses the same profile. It gives a rod an offset of
`0.5 * (1 + cos(pi * delta / ramp))` inside the ramp and zero outside it, where
`delta` is the rod azimuth less the groove azimuth. The stroke is 1.5e-9 m.

The viewer also shows the mechanism. `site/index.html` adds a **Hide rotor top**
control, and `viewer.js` hides the rotor body atoms above the rotor mid-plane
when the control is on. The user asked to keep the housing in place, because the
outer housing gives context and the rotor disc is the part that covers the pins
and the groove. `uploadAtoms` and the 2D fallback both apply the cutaway, so the
radial bores, the follower pins and the cam groove below the rotor become
visible.

The measured scene now holds 157857 atoms and 2.2505103227300436e-21 kg. The cam
plate holds 55173 atoms, the housing 39848, the rotor body 54844 (49708 rotor and
5136 shaft), and each rod 666 (602 rod and 64 pin). The headless harness
measured the rod peak as 1.475e-9 m, which is 98.4 percent of the stroke and is
limited by the sample step. The worst deviation from the profiled model was
3.2e-12 m. The cam plate did not move. The largest rod radius over a turn was
6.91e-9 m, below the housing inner radius of 7.2e-9 m, so no rod crosses the
housing. The harness reported zero console errors.

Limits on this work. The ramp angle and the rise are chosen values, not a tuned
motion law, so a profiled groove could still change the acceleration. The pin
slides in the groove; a roller follower would roll. The scene models no contact
force between the pin and the groove wall, no torque on the drive shaft, and no
external machine. The scene still holds no guest, no solvent and no molecule in
a pocket.

## ADR-0072

### A leaf spring returns each rod, and the cam is a one-sided hub

Status: accepted. This supersedes ADR-0071.

ADR-0071 gave the cam a closed groove that both pushes and retracts a rod. A
closed groove needs a groove wall on each side, so the cam is a wide plate with a
long slot. The user asked for the part that actually pushes the rod and pulls it
back. A one-sided cam with a return spring is smaller and shows the load path.

The cam is now `CamHubGenerator`. It is a profiled disc, one-sided. Its surface
radius is `base_radius_m` outside a ramp of half-angle `ramp_half_angle_rad`
about `angle_rad`, and it rises by `rise_m` over the ramp with a raised cosine.
The disc carries a central bore of `bore_radius_m`. `CamPlateGenerator` and its
closed groove stay in the library, but the rotor scene no longer uses it.

Each rod carries a `LeafSpringGenerator` beside its follower pin. The spring is a
diamond cantilever of `SPRING_LENGTH_M` (0.6e-9 m), `SPRING_WIDTH_M` (0.4e-9 m)
and `SPRING_THICKNESS_M` (0.2e-9 m). Its tip stiffness is `3EI/L^3` with
`E = 1.05e12 Pa`, which is 3.89 N/m. The spring holds the pin against the hub, so
the hub pushes the rod out and the spring pulls it back. The spring is a
schematic of the load path. It is drawn rigid, it applies no force in the scene,
and a 1.5e-9 m stroke is large for a 0.6e-9 m cantilever, so the stiffness is
indicative.

A first build of this scene placed the parts wrong. The hub was generated at its
default thickness of 0.9e-9 m plus a 1.0e-9 m flange, but the code placed it as
if it were only as thick as the gap under the rotor. The hub top rose about
0.18e-9 m into the rotor and the hub bottom hung about 0.57e-9 m below the
housing. The rotor-rod, rotor-hub and rod-hub interfaces also overlapped, because
a carbon surface carries hydrogen caps that reach about 0.11e-9 m past it.

The scene now uses `INTERFACE_CLEARANCE_M` (0.4e-9 m) at every moving interface.
The housing is 4.6e-9 m thick, so the chamber is tall enough. The hub is
generated at the true gap thickness and it has no flange. The rotor bore is
widened by the clearance, and the follower pin top sits one clearance below the
rotor. The gap is measured between nominal carbon surfaces, so 0.4e-9 m keeps
every inter-part atom pair above the 0.16e-9 m clash threshold.

A new metric `nanocad_meter::relax_subassembly` confirms the assembly is stable.
It relaxes the covalent lattice only, with ideal rest lengths (C-C 1.544e-10 m,
C-H 1.09e-10 m, H-H 7.4e-11 m), and it counts non-bonded pairs closer than
1.6e-10 m as clashes. It does not use van der Waals: the reaction forces in a
mechanism come from the joints, not from chemistry, and a stiff Buckingham term
against a soft lattice has no stable equilibrium. The metric reports the energy
drop, the peak force, the largest bond strain, the clash count and the largest
displacement.

The relaxation of the rod and the cam hub (1692 atoms, 2259 bonds) converges. The
energy falls by 2.63e-21 J, the largest bond strain is 6.6e-6, the largest
displacement is 2.75e-13 m, and the clash count is zero before and after. The
scene holds 127330 atoms and 1.8243665191584498e-21 kg.

Limits on this work. The relaxation covers the rod and the hub, not the whole
scene; the rotor and the housing are large and their cut-out boundaries carry
broken bonds that give a false collapse. The spring makes no contact and applies
no force. The scene still models no contact force between the pin and the hub, no
torque on the drive shaft, no guest, no solvent and no molecule in a pocket.

## ADR-0073: The rod shaft is the wide piston face

Status: accepted. It refines ADR-0072.

Context. A guest can bind anywhere in the pocket, not only on the pocket axis.
The old ejection rod was a wide shaft that stepped down to a narrow tip. The
narrow tip was the part that entered the pocket. A new geometric metric,
`nanocad_meter::capture`, asks how much of the pocket this face can cover.

Finding. The narrow tip (0.2 nm radius) covers only 0.52-0.54 of the pocket. The
shaft face (0.4 nm radius) covers 0.77-0.79. A guest against the sidewall beyond
the reach of the face is not pushed. A face of radius 0.583e-9 (dimethyl ether)
or 0.563e-9 (ethanol) covers every allowed offset.

Decision. The rod shaft becomes the pushing face. The default
`shaft_radius_m` is 0.65e-9 and the default `tip_radius_m` is 0.4e-9. The step
guard stays: the tip is narrower than the shaft. The rotor ejection bore is
`rod_radius_m + EJECTION_BORE_CLEARANCE_M` (0.94e-9). The bore cannot use
`INTERFACE_CLEARANCE_M`, because it must stay below the 1.0e-9 pocket radius.

Consequence. The shaft face covers every allowed offset for ethanol and dimethyl
ether (coverage 1.0, dead zone 0). The cost is a wider bore and a heavier rod.
The scene now holds 127036 atoms and 1.787369134208518e-21 kg. The rod and hub
relaxation still converges with 0 clashes (2142 atoms, 2993 bonds, bond strain
1.9e-6).

Limits. The metric is a rigid geometric coverage model. It uses one sphere
radius for the guest. It does not model the guest shape, the contact force, or
the guest motion. A face that covers the pocket can still miss a guest that
leaves the pocket before the rod arrives.

## ADR-0074: Capped surfaces need a wider clearance than the atom radius

Status: accepted. Refines ADR-0073.

Context. A user reported that the rod parts below the rotor bottom plate
overlapped the rotor during the in-and-out motion. A swept check measured the
least distance between each rod and the rotor over the full 1.5e-9 stroke. At
every stroke fraction the least distance was 0.107-0.130e-9, and 68 to 100 atom
pairs were closer than the 0.16e-9 clash threshold.

Finding. The closest pairs were rod-shaft hydrogen atoms against rotor-bore-wall
hydrogen atoms. The bore clearance was `EJECTION_BORE_CLEARANCE_M` = 0.29e-9.
The rod surface and the bore wall are both capped with hydrogen. Each cap
reaches about 0.11e-9 past the carbon surface, so the cap-to-cap gap was only
`0.29 - 0.22` = 0.07e-9. `INTERFACE_CLEARANCE_M` is 0.40e-9 for the same
reason on planar faces. The bore used less because it must stay below the
1.0e-9 pocket.

Decision. The bore clearance `EJECTION_BORE_CLEARANCE_M` is 0.39e-9. The
default `shaft_radius_m` is 0.60e-9. The bore is then 0.99e-9, strictly below
the pocket. The tip radius stays 0.4e-9.

Consequence. The least rod-to-rotor distance is 2.5e-10 to 2.6e-10 at every
stroke fraction, and no pair is closer than the clash threshold. The shaft face
still covers every allowed offset for ethanol and dimethyl ether (coverage 1.0
for both, though the dimethyl ether margin is only 0.017e-9). The scene now
holds 120884 atoms and 1.715827677935783e-21 kg. The rod and hub relaxation
converges with 0 clashes (1962 atoms, 2709 bonds, bond strain 4.4e-6). A
permanent test, `the_rods_clear_the_rotor_through_the_stroke`, holds the limit.

Limits. The check covers the rod and the rotor only. It does not cover the
guard, the leaf spring bond to the rotor, or a guest in the pocket. The capture
margin for dimethyl ether is small.

## ADR-0075: Surface clashes are fixed in the lattice, not in the cap

Context. The physical-sense audit (`sense_report`) found hydrogen-hydrogen
pairs 0.74e-10 apart on curved diamond surfaces, in every scene. Each pair is
two hydrogen caps that point into the same missing lattice site. The site
separates two surface carbons.

Finding. A minimum-separation rule in `cap_free` removes every such pair. It
also removes 11 750 hydrogen caps from the rotor scene. The carbon valence then
falls below four. The parts tests `every_carbon_keeps_a_valid_valence` and
`every_solid_gear_atom_has_a_valid_diamond_valence` then fail with carbon
degrees of three, two and one. A degree-one carbon cannot come from a refused
cap alone. The source lattice therefore contains isolated and under-coordinated
surface sites.

Decision. `cap_free` does not refuse a cap. A cap is always added, and a
hydrogen-hydrogen clash stays visible. The audit reports the clash. The root
cause is in lattice generation and surface reconstruction, not in capping.
`spur_gear` has a known duplicate: the rim, the hub and the spokes are three
independent atom sets, so coincident sites duplicate exactly.

Consequence. The tree stays green, and both valence tests still pass. The audit
tool `sense_report`, its unit tests and the `sense_json` example are merged. The
open work is the lattice dedup and surface-reconstruction task.

Limits. The audit threshold `MIN_SEPARATION_M` is 0.9e-10. A real hydrogen
bond is near 0.74e-10, so the audit flags a clash that a relaxed structure
could still tolerate. The finding is not yet traced to a single lattice
function.

## ADR-0076: The pin and the leaf spring are separate bodies on fixed joints

Context. The user reported that the follower pin and the leaf spring are not
attached to the rod. They are not bonds, and their atoms sit in the rod body.
The pin hangs 0.74e-9 m below the rod, and the spring hangs 0.60e-9 m below it.
No perceived bond joins them, so the pin cannot push the rod and the spring
cannot return it.

Finding. A bonded pin is not possible. The rotor bore is a radial cylinder with
`y^2 + z^2 < ejection_bore_radius^2`. Its radius is 0.99e-9 m. The rotor
half-thickness is 1.0e-9 m. The generator also refuses a bore radius at or above
the pocket radius of 1.0e-9 m. So the bore cannot cut the bottom face, and a
0.01e-9 m skin stays. A vertical pin at radius 2.75e-9 m must cross that skin.

Decision. The pin and the leaf spring are their own bodies. A fixed joint
attaches the pin to its rod, and a fixed joint attaches the spring to the rotor.
The rotor scene holds 39 bodies and 37 joints.

Consequence. The linkage is attached, and `sense_report` finds no floating body
and no split linkage body. The joint is the attachment, not a lattice bond, so a
small gap stays in the model. The viewer draws the joints and the body markers in
the 3D pass at high zoom, so the atoms occlude them.

Limits. The audit still reports 84 cross-body perceived bonds and two split
bodies (the housing into 2 parts, the rotor into 9). These are older defects and
they are not from this change. A real lattice bond needs a through-bore, a
thinner rotor, or a counterbore.

## ADR-0077: The pin is machined into the rod through a rotor clearance slot

Context. ADR-0076 attached the follower pin to its rod with a fixed joint,
because the rotor skin blocked a bond. The user asked for a real bond, not a
joint. The ejection bore is a radial cylinder with a radius below the pocket
radius of 1.0e-9 m. The rotor half-thickness is 1.0e-9 m. So the bore cannot cut
the bottom face, and the pin cannot cross the skin (ADR-0076).

Finding. A second void can cut the skin. The rotor generator now takes three
parameters: `pin_clearance_orbit_m`, `pin_clearance_half_length_m` and
`pin_clearance_half_width_m`. It cuts one radial `Box3` slot for each pocket,
through the full thickness, at the pin orbit. The slot is not the bore.

Decision. The sorting rotor cuts a pin clearance slot. The rod generator
`EjectionRodGenerator::generate_with_follower_pin` machines the pin into the rod
inner end from a `FollowerPinFeature`. The scene sets the orbit to
`cam_retract + rise/2`, the half-length to `rise/2 + head_radius + margin`, and
the half-width to `head_radius + margin`, with `margin = 3.5e-10 m`. The pin is
now one lattice with the rod. The rotor scene holds 27 bodies and 25 joints.
ADR-0076 is superseded.

Consequence. The stroke test passes with no rod-to-rotor clash. `sense_report`
shows no rod or spring body floating or split. The spring keeps its fixed joint
to the rotor, because the rotor still blocks a bonded spring.

Limits. The margin of 3.5e-10 m is large, because the pin cap extends about
0.11e-9 m past the carbon surface. A margin of 1.5e-10 m left a 1.19e-10 m clash
at rod 6. The scene still reports 84 cross-body perceived bonds and two split
bodies (housing into 2 parts, rotor into 9). These are older defects.

## ADR-0078: The cam pressure angle governs the drive, and it is too steep

Context. The user asked if the follower pin is angled enough to push the rod
smoothly. The pin is vertical. The cam is a radial cam: a disc whose outer
radius rises with azimuth (`cam_hub.rs`). The pin rides the outer wall as a
roller follower, and the rod slides along its own radius. The push is radial, so
the pin must not be angled. The governing angle is the pressure angle.

Finding. For a translating radial follower, `tan(angle) = (dr/dtheta) / r`,
where `r` is the follower-centre radius. The new `CamDrive` metric in
`crates/meter` samples one turn and reports the pressure angle, its position,
the mechanical advantage, the side-load ratio, and the lift derivatives. At the
current defaults the maximum pressure angle is `74.99 deg`, and the mechanical
advantage is `0.26`. A smooth cam wants at most `30 deg`. The ramp is capped at
15 deg, because the 12 pockets are 30 deg apart and only one pin may ride the
ramp. So the 12-pocket rotor cannot reach 30 deg with a 1.5 nm stroke on a
2-3 nm hub.

Decision. `crates/meter/src/cam_drive.rs` holds the `CamDrive` metric. It uses
`CamHubGenerator::surface_radius_with`, which is now public. The example
`cam_drive_json` prints the report. `crates/opt/src/cam_drive.rs` holds
`CamDriveSearch`, a CMA-ES search over `base_radius_m`, `ramp_half_angle_rad`
and `pin_radius_m`. The rise holds the stroke target, so every candidate is
feasible. The hub base is limited to `3.0e-9 m`, because a larger lobe would
breach the pocket inner radius. The objective is the maximum pressure angle plus
a jerk penalty.

Consequence. The search lowers the maximum pressure angle from `74.99 deg` to
`63.50 deg` at the hub limit, the ramp limit and a pin radius of `0.80e-9 m`.
The design stays feasible and not smooth.

Limits. This is a model search. It reads the cam profile and the pin radius
only. It models no friction, no spring preload, no contact force, and no
materials. A smooth angle needs a much larger hub (near `2.0e-8 m`), a shorter
stroke, a closed-groove cam, or fewer pockets. This ADR records the estimate; it
does not claim a validated design.

## ADR-0079: Face skins keep every rotor-scene body one piece

Context. `sense_report` on the rotor scene found the housing split into two
halves and the rotor split into 33 pieces. A body with more than one piece is
not a physical part. The cause was the through cuts. The two housing channels
ran through the full thickness, so they cut the ring into two arcs. The rotor
ejection bores and pin clearance slots also ran through the full thickness, so
they cut the disc into a hub and many sectors.

Decision. Each functional cut stops at both faces, and a solid skin stays.

- `housing.rs` keeps `CHANNEL_SKIN_M = 3.5e-10 m`. Each channel is a window in
  the mid-plane. The ring keeps a top and a bottom skin, and the guest still
  passes in the mid-plane.
- `rotor.rs` keeps `TOP_SKIN_M = 1.5e-10 m`. The 12 ejection bores are clipped
  under both faces. The 12 pin clearance slots are clipped under the top face
  only, because the pin exits the bottom face. The central bore, the keyway and
  the pockets stay through-thickness.

The rotor skin is thin because the rod slides under it. The rod top sits at
`z = 0.6e-9`, and the rotor half-thickness is `1.0e-9`. A `1.5e-10 m` skin puts
the plate at `z = 0.85e-9..1.0e-9`, and the rod keeps a clearance above the
`1.7e-10 m` clash threshold.

Consequence. `sense_report` reports no disjoint body and no floating body for the
rotor scene. Every one of the 27 bodies is one connected lattice. The rotor
holds 38 472 atoms and the housing 73 208.

Limits. The change does not remove the pre-existing `overlap_count` from
sub-lattice H-H clashes (ADR-0075), and it does not remove the 84 cross-body
perceived bonds. Those stay open under M21-03.

## ADR-0080: A wide cam ramp gives a gradual push

Context. The user said the cam lobe is too abrupt. It looks like a straight
bump, and the pin seems to crash instead of being pushed. The pressure angle
at the default ramp of `0.21 rad` is `75 deg`. The pin then has a mechanical
advantage of `0.26`, and the side load is four times the useful force.

Decision. The rotor scene uses `CAM_RAMP_HALF_ANGLE_RAD = 1.4`. The cam is
generated with that ramp, and both pin-centre computations call
`CamHubGenerator::surface_radius_with` with the same value, so the scene and
the pin geometry agree. `site/viewer.js` sets `ROTOR_CAM_RAMP_RAD = 1.4`.

Consequence. The pressure angle is near `29 deg`, and the push is gradual. The
cam surface grows from 1 184 to 1 822 atoms. The stroke stays at `1.5 nm`.

Trade-off. A ramp of `1.4 rad` covers `2.8 rad` of the turn, so about five pins
sit on the ramp at once, and several pockets open together. If the machine
must open one pocket at a time, then the ramp must be short, and the stroke
must fall to about `0.45 nm`, or the hub must grow to about `20 nm`. This ADR
keeps the full stroke and accepts the wide window. M23-03 records the choice.

## ADR-0081: The gas guests use point charges, and the loading model is kinetic

Context. The user asked for a binding site that distinguishes oxygen from
carbon dioxide, and for the binding dynamics at physiological concentrations.
The guest table held five organic molecules. It held no gas.

Decision. The guest table gains `o2` and `co2`. The RDKit emitter cannot run
here, so both entries are hand-authored from published force-field parameters,
not from a conformer search:

- `o2` is a rigid linear molecule with a `1.208e-10` bond, a molar mass of
  `3.199800e-2`, and zero atom charges. The neutral charge is correct for the
  ground-state triplet at the point-charge level.
- `co2` is a rigid linear molecule with a `1.16e-10` bond, a molar mass of
  `4.400900e-2`, a carbon charge of `+0.70e`, and an oxygen charge of `-0.35e`.
  These are the TraPPE charges of Potoff and Siepmann (2001).

The header of `guest_data.rs` states the source and the limits.

A new `nanocad-meter` module, `loading.rs`, estimates the loading time. The
arrival rate is the gas kinetic flux onto the pocket mouth, `n * v / 4`, times
the mouth area and a sticking probability. The loading angle is the rotor
angle that the saturation time spans, `omega * t`.

Consequence. The binding metric holds the carbon dioxide far more strongly
than the oxygen in every wall group. The gas ratio runs from `0.022` (amino)
to `0.071` (chloro). The carbon dioxide is worth `14` to `36 kT`, and the
oxygen is worth `0.5` to `1.2 kT`.

At alveolar partial pressures (`O2 1.33e4 Pa`, `CO2 5.33e3 Pa`), `310.15 K`, a
`0.6 nm` mouth and a `1e6 rad/s` rotor, one pocket fills from oxygen in
`1.58e-9 s` and from carbon dioxide in `4.63e-9 s`. The carbon dioxide sets the
loading angle at `0.266 deg`. One pocket spacing `30 deg` takes `524e-9 s`, so
the gas arrival is not the limit.

Trade-off. This ADR does `not` deliver an oxygen-selective site. The model
holds a neutral oxygen and a charged carbon dioxide, so the Coulomb term
always favors the carbon dioxide. An oxygen-selective site needs a model that
this metric does not have: a magnetic or charge-transfer center, or an
explicit quadrupole. The result is stated as a limit, not as a success. The
loading model states no liquid diffusion, no crowding, and a rigid pocket.

The gas guests carry no independent validation. Every value is a model value.

## ADR-0082: One rotor per target, and the delivered fraction sets the design

Context. The user asked for the optimal rotor for oxygen and for carbon dioxide
in blood. The user then clarified the architecture: the device uses one rotor
for each target molecule, and each rotor has a fixed binding site. A single
site does not switch between molecules. A pH-sensitive site is useful but
optional.

Decision. `crates/meter/src/delivery.rs` models the delivered fraction of one
site over one concentration swing. A site follows the Hill curve
`theta(P) = P^n / (P50^n + P^n)`. The delivered fraction is the occupancy at
the loading pressure minus the occupancy at the unloading pressure. The `P50`
that maximises the fraction is the geometric mean of the two pressures, for
every Hill coefficient.

The binding energy and the `P50` are the same statement:
`P50 = P0 * exp(dG / kT)`, where `P0` is one bar. A design names one of the two.

Consequence. The model sets the design target for each rotor:

- Oxygen, 13.3 kPa to 5.3 kPa: a flat site delivers `0.226` of its capacity
  per cycle, at `P50 = 8.40 kPa` and `dG = -6.39 kJ/mol`.
- Oxygen, 13.3 kPa to 3.0 kPa: a flat site delivers `0.356`. A cooperative
  site with a Hill coefficient of `2.8` delivers `0.779`, at `P50 = 6.32 kPa`.
- Carbon dioxide, 6.1 kPa to 5.3 kPa: a flat site delivers `0.035`. At the
  tissue swing of 6.7 kPa to 5.3 kPa it delivers `0.059`.

So the oxygen rotor is feasible with a fixed site and a large pressure swing.
The carbon dioxide rotor is not efficient as a free-gas pump, because the
physiological pressure swing is small. The carbon dioxide rotor must move a
larger pool (bicarbonate, about `24 mM`, against dissolved CO2 of about
`1.2 mM`) or use a chemical reaction (carbamate on an amine).

Trade-off. The model is a single equilibrium site. It states no kinetics, no
hysteresis, no competing ligand and no crowding. The delivered fraction is an
upper bound for a real site. The `P50` values are model targets, not measured
affinities. ADR-0081 still holds: the current point-charge binding metric
cannot prefer oxygen, so the oxygen site needs the metal-ligand term of M25-02.

## ADR-0083: The oxygen site gets a metal-ligand coordination term

Context. ADR-0081 measured that a wall of Lennard-Jones plus point charges
binds carbon dioxide more strongly than oxygen in every functional group. The
cause is physical: carbon dioxide is larger and more polarizable, and it carries
the TraPPE charges. Oxygen cannot win on dispersion. Its advantage is
electronic, because it has two unpaired electrons and empty pi-star orbitals
that back-donate to a metal. A point-charge model cannot see that term. Task
M24-04 asks for an oxygen-selective site.

Decision. `crates/engine/src/metal_ligand.rs` adds a Morse coordination term
between a metal atom and a ligand atom:
`U(r) = D * (1 - exp(-a * (r - r0)))^2`.
The term is zero at the equilibrium distance `r0` and equals `D` at infinite
separation. `MetalLigandParams::from_force_constant(D, r0, k)` sets
`a = sqrt(k / (2 D))`. A gradient test with a central finite difference passes
to one part in `1e7`.

The term is ligand-specific. Oxygen gets one coordination well at a metal
center. Carbon dioxide gets none, because its dispersion and electrostatics are
already in the classical pocket term. The site is a metal macrocycle with one
open axial position and a distal hydrogen-bond donor, in the spirit of a heme.

An example `crates/engine/examples/metal_ligand_json.rs` uses one iron-oxygen
contact at `1.90e-10 m`, a well depth of `6.642e-20 J` (`40 kJ/mol`, `15.5 kT`)
and a force constant of `200 N/m`. The largest attractive force is `1.29 nN`.
The average force over a `1.5 nm` ejection stroke is `44 pN`.

Consequence. A `40 kJ/mol` oxygen bond is larger than a typical physisorption
well, so a metal site can prefer oxygen over a bare cavity. The same number is
comparable to the carbon-dioxide binding that the charged pocket already
reported (14 to 36 kT), so the two models must be reconciled. Published heme
oxygen enthalpies are about `-40` to `-60 kJ/mol`, which supports the chosen
depth as a starting value.

Trade-off. The well depth is a zero-Kelvin potential energy. The delivery model
of ADR-0082 states a free energy near `-7 kJ/mol` against one bar, because gas
binding loses translational and rotational entropy. The reconciliation between
the two is open (M25-06). The depth here is an assigned starting value, not a
measured affinity. The term is not yet integrated into the binding search
(M25-02). The engine has no carbon-dioxide coordination term, by design.

## ADR-0084: Quantum chemistry calibrates the oxygen well depth

Context. ADR-0083 added a metal-ligand Morse term and assigned a `40 kJ/mol`
well depth from the literature. Task M25-06 asks for a calibration and a
reconciliation with the delivery model of ADR-0082.

Decision. Run a small quantum-chemistry calibration with pyscf 2.14.0 (PBE,
`def2-svp`, grid level 1) on a tractable proxy for the heme metal site: a
ferrous iron with four square-planar ammonia ligands at
`Fe-N = 2.15e-10 m`. The full iron-porphine model did not converge in a
reasonable time, so the ammonia proxy stands in. The ammonia ligands are pure
sigma donors, so the model likely understates oxygen back-donation.

Results. Oxygen approaching the iron axially, in the triplet state, gives
binding energies of `-115.7 kJ/mol` at `1.9e-10 m`, `-95.6 kJ/mol` at
`2.3e-10 m` (converged) and `-78.3 kJ/mol` at `2.5e-10 m`. Carbon dioxide, with
its oxygen nearest the iron, gives `+0.6 kJ/mol` at `3.4e-10 m` (essentially no
binding) and `+9.4 kJ/mol` at `3.0e-10 m` (repulsive). The singlet
iron-ammonia ground state is `-1488.878981` hartree.

Consequence. The metal center binds oxygen and rejects carbon dioxide. This is
the electronic selectivity that the classical point-charge metric of ADR-0081
cannot produce. The example now uses a well depth of `1.3284e-19 J`
(`80 kJ/mol`), the midpoint of the corrected estimate, in place of `40 kJ/mol`.

Reconciliation. Binding an oxygen molecule from the gas removes its
translational and rotational entropy, about `78 kJ/mol` at `310 K`
(`205.1` plus `46.0 J/mol/K`). The free energy is therefore
`dG = dH + 64..78 kJ/mol`. The delivery model wants `dG` near `-7.25 kJ/mol`
for a target `P50` of `6 kPa`, which needs `dH` near `-71` to `-85 kJ/mol`. The
calibrated well depth sits in that window. The reconciliation is approximate.

Trade-off and limits. The quantum chemistry has no dispersion module, so carbon
dioxide physisorption is absent from it. The classical metric gives `14` to
`36 kT` for that wall term, so the true comparison is oxygen coordination
against carbon-dioxide physisorption, and oxygen still wins by a factor. There
is no counterpoise correction, several open-shell states did not converge, and
only one functional, one basis and one model were used. Every number here is an
estimate, not a measurement. The rotor remains a simulated design, and no
medical claim is made. The term is still not integrated into the binding search
(M25-02).

## ADR-0085: The viewer builds a solid model from the atoms

Status: accepted.

Context. At low zoom, the atoms became too small to draw, and the viewer showed
nothing for the rotor. The user wants each component to look like a 3D CAD part
at that scale, and wants to show or hide each component.

Decision. The browser builds one surface mesh for each rigid body from the
atoms of that body. `site/solid_mesh.js` sums a smooth field around each atom
(radius `3.2e-10 m`, hydrogen weight `0.6`) and extracts the surface with
surface nets on a grid of `0.6e-10` to `1.6e-10 m`. A lone carbon gives a
closed sphere of `1.73e-10 m`. `site/render_atoms.js` draws the meshes with the
same projection as the atoms, writes a normal and a body id to a second target,
and draws outline, crease and depth edges in the composite pass. Each body moves
as a rotation about z and a shift, the same as the atom motion. The scene schema
does not change.

Consequence. The Auto model mode shows atoms near and the solid model far. The
Components panel hides bodies in both modes. The mesh is an envelope of the
atoms for display. It is not a measured surface and carries no physical result.

Trade-off and limits. The mesh is rebuilt once for each scene, from the rest
positions. A body that deforms would need a rebuild. The upload took about 4.6 s
for the gearbox in the headless software renderer; it is not measured on a
hardware GPU (M26-04).

## ADR-0086: A selective gear pump replaces the rod rotor

Status: accepted. Supersedes the rod rotor as the target design (ADR-0078 to
ADR-0080 stay as the Freitas reference). Corrects the entropy bound in ADR-0084.

Context. The user asked for the minimum machine that binds one molecule on one
side of a membrane and releases it on the other side. The rod rotor has faults.
Its pocket (radius `1.0e-9 m`, depth `2.0e-9 m`) holds about `5 nm^3`, which is
150 to 200 water molecules. Its rim is one long sliding seal between the inlet
and the outlet. Its twelve rod bores leak. An empty pocket returns open. It
needs 27 bodies, 25 joints, a cam with a `2.8 rad` ramp, and a rod return that
does not exist.

Decision. The target design is a selective gear pump: two meshing spur gears in
a figure-eight housing, and the housing is the membrane. Each tooth gap is a
one-molecule pocket with a metal site at the tooth root. A gap loads at the
suction port, carries its guest along the sealed casing arc, and the mating
tooth fills the gap at the mesh. That is forced release by positive
displacement. The meshed gap has no volume, so nothing returns. The machine has
three bodies and two revolute joints with one gear coupling.

Kill test. Water is smaller than oxygen (`2.65e-10` against `3.46e-10 m` kinetic
diameter), so size cannot exclude it. The test compares binding free energies
with each guest at its partial pressure; water enters at its vapour pressure
(`6.28 kPa` at `310 K`), which is exact for a site in the liquid.
`crates/meter/src/competition.rs` computes the multi-component Langmuir
occupancy. The same quantum-chemistry model as ADR-0084 (pyscf PBE/def2-svp,
Fe(NH3)4 proxy) binds water by `-10.3 kJ/mol` at `2.7e-10 m`, against `-95.6`
for oxygen. At alveolar pressures the metal site holds oxygen `0.99997` of the
time, water `1.5e-14` and carbon dioxide `9e-18`. A bare methyl pocket of radius
`4.5e-10 m` without a metal stays empty (`1 - 2e-8`) and prefers water to oxygen
20 to 1. So the metal site is mandatory, and it passes the kill test by about
`14` orders of magnitude. A +-30 kJ/mol error leaves more than `1e8`.

Entropy correction. The standard entropy of oxygen gas, `205.2 J/(mol K)`, is
the total of translation, rotation and electronic terms. ADR-0084 added rotation
twice for its upper bound of `78 kJ/mol`. The correct cost at `310 K` is about
`64 kJ/mol`, so a passive site with `P50` near `6 kPa` needs `dH` near
`-71 kJ/mol`. The gear pump does not need that compromise: its site can bind
strongly (`dG` near `-32 kJ/mol`, `P50` well under `1 Pa`), because the mesh
forces the release.

Trade-off and limits. The kill test assumes the guest loses all gas entropy, no
dispersion, no counterpoise correction, one functional, one basis and one small
model. The seal gaps (backlash, tip, face) must be smaller than a water
molecule, and tooth stiffness at a module near `3e-10 m` is unproven. Every
number here is an estimate from a simulated design. No medical claim is made.

## ADR-0087: Gear-pump geometry: eight teeth, housing posts, tuned gaps

Context. M27-02 builds the gear pump of ADR-0086 as parts and a scene. The
existing `spur_gear` part is a skeletal rim with spokes, so it cannot seal. The
planetary set already fills a solid hydrogen-capped diamond gear from an
involute outline.

Decision. `crates/parts/src/gear_pump.rs` (`GearPumpGenerator`, id `gear_pump`)
builds two equal solid gears and one housing in a single part with three atom
ranges. The gears use module `0.6e-9 m`, eight teeth and a `30 deg` pressure
angle, which is the smallest undercut-free count at that angle. Each gear has
eight (001) layers and a bore of radius `1.0e-9 m`. The housing is one lattice:
two side plates, a figure-eight wall, a suction port at `-y` and a discharge
port at `+y` at the mesh waist, and one post through each gear bore. The posts
join the plates, so no body floats. The gear is held by its post and the plates,
not by its tips.

Every gap is a carbon-to-carbon distance near `4.5e-10 m`. Two hydrogen caps
then leave about `2.3e-10 m`, less than the `2.8e-10 m` diameter of water. The
housing clearance is `4.6e-10 m`. The gear-to-gear gap was tuned by a scan over
one tooth pitch: backlash `4e-10 m` and a centre offset `4e-10 m` give a worst
gap of `4.50e-10 m`. Backlash `2e-10 m` gave `2.8e-10 m`, where the caps touch.
A pruning pass removes lattice carbons with fewer than two carbon neighbours and
keeps the largest connected piece.

`crates/jigs/src/gear_pump_scene.rs` gives three bodies, two revolute joints at
the post centres and one gear coupling. The scene catalog lists it as a machine
with the id `gear_pump`. The viewer reads a zero design with a gear coupling as
joint-driven motion: each gear turns about its own joint anchor, and the
coupling sets the opposite rate.

Fix. `GearProfile::minimum_tooth_count` returned nine at `30 deg`, because
`2 / sin^2(30 deg)` rounds above eight in floating point. A `1e-9` tolerance
now gives eight.

Limits. The sense audit shows no cross-body overlap, no floating body, no
disjoint body and no cross-body bond. It still flags the same-body hydrogen cap
clash at a shared missing site (ADR-0075), which is in every scene. The scene
has 53 366 atoms, 48 206 of them in the housing. The part holds no binding site
(M27-04). No drive torque, contact force or solvent is modelled. The gaps are
from the rigid lattice, not from a relaxed structure.

## ADR-0088: 2x1 surface dimers remove the cap clash; a valence field checks it

Context. ADR-0075 left a defect in every scene: two hydrogen caps 0.74e-10 m
apart in one missing lattice site. On a (001)-like face each missing site sits
one bond from two surface carbons. The gear pump had 7 026 such pairs.

Decision. `diamond_solid::bond_and_cap` now reconstructs the surface before it
caps it, as real diamond C(100)-(2x1):H does. `plan_dimers` builds a graph of
carbons that share a site and takes a leaf-first matching. `apply_dimers` bonds
each pair, drops the two directions into the shared site, and re-aims the last
cap of each dimer carbon to complete its tetrahedron. A pair moves together to
the target dimer length `1.60e-10 m` (Furthmuller et al. 1996), as far as its
back bonds stay under `1.615e-10 m`. At a step, one back bond lies along the
pair axis (it would reach 1.94e-10 m), so that pair stays strained, up to
`2.6e-10 m`, until a relaxation. A site that three or more caps share is filled
with carbon (at most four rounds). Filling a two-cap site grew each face one
row per round and put carbon 5e-10 m into a gear, so it is not done.

Results. Gear pump: 7 026 overlaps to 0, clean audit before relaxation. Rotor:
13 649 overlaps to 6 053. The rest of the rotor overlaps are older defects in the
rotor body (5 466 pairs, worst 0.2e-10 m) and its cross-body contacts.

Relaxation. `nanocad_meter::relax_valence` is a harmonic valence field on the
topology bonds: MM3-scale bond and angle constants, tetrahedral rest angle,
diamond C-C rest length, and UFF Lennard-Jones on every pair that is not 1-2,
1-3 or 1-4. It has no torsion term. The pump gear converges in 605 iterations
(1 s). Its RMS displacement is 1.0e-10 m and its largest carbon moves
3.0e-10 m. The thin gear warps as a whole. The housing does not converge in 300
iterations (55 s). In the relaxed scene the housing comes within 0.71e-10 m of
a gear at 8 pairs. So the rigid 4.6e-10 m clearance does not survive this
relaxation. Without a torsion term the bending stiffness is uncertain, so the
size of the warp is an estimate.

Other changes. The sense audit perceives a bond below `1.62e-10 m` or below
`1.15` times the sum of the covalent radii (C-C 1.75e-10 m), so a strained
relaxed bond still counts. The rod bore clearance is `4.5e-10 m` and the rod
radius `0.54e-9 m` (bore still 0.99e-9 m): re-aimed caps point straight
across the gap, and the old `3.9e-10 m` left 1.44e-10 m between hydrogen atoms.
The rotor rod-and-hub check uses `relax_valence`. The old bond-only
`relax_subassembly` has no angle stiffness, so it folds a strained surface.

Limits. A carbon with one lattice bond (a corner group) is not paired. The
gear pump prunes such carbons; other parts keep them, and a plain box keeps
38 clash pairs at its corners. The topology cannot remove atoms, so the pass
cannot either.

## ADR-0089: An iron site in each tooth root; the gap is too large

Context. M27-04 puts the metal binding site of ADR-0086 into the gear pump
and integrates the metal-ligand term into the binding search (was M25-02).

Decision. `GearPumpGenerator` turns the root cap nearest each gap centre into
an iron atom, 1.95e-10 m from its carbon along the cap direction, open toward
the gap. Gear B takes the same local atoms. The model tables now hold iron:
Cordero low-spin radius 1.32e-10 m and UFF Fe (2.912e-10 m, 0.013 kcal/mol).
`nanocad_meter::metal_site_binding` places the donor oxygen on the open side
of the iron over 12 distances, 48 directions and 8 spins. A pose energy is the
UFF wall term plus a Morse coordination term, which replaces the Lennard-Jones
pair between the donor and the iron. The wells are from the Fe(NH3)4 proxy:
oxygen 80 kJ/mol at 1.90e-10 m, water 10.3 kJ/mol at 2.70e-10 m, carbon
dioxide none. One width, `a = 2.744e10 1/m`, matches k = 200 N/m for oxygen.

Results (`pump_site_json`, gap 3 of gear A, rigid, no solvent). Oxygen binds
at -82.6 kJ/mol (-79.5 coordination, -3.1 wall), water at -13.3 and carbon
dioxide at -5.5 (wall only). At alveolar pressures and 310 K the site holds
oxygen 0.995 of the time, water 7e-12 and carbon dioxide 1.4e-14. During the
mesh the nearest carbon of the other gear comes no closer than 5.46e-10 m to
the iron. A bound oxygen reaches only about 3.1e-10 m from the iron, and a
carbon excludes about 1.7e-10 m around itself, so the space in front of the
iron stays free. The mesh does not release the oxygen. (An earlier text of
this ADR said the reverse; that was wrong, see ADR-0090.)

Finding. The free volume of one gap, from the root to the housing bore and
outside every Bondi sphere, is 3.05e-27 m^3: about 100 molecules of liquid
water. ADR-0086 called the gap a one-molecule pocket; at module 0.6e-9 m it is
not. The site selects oxygen, but each gap would carry about 100 water
molecules with each oxygen. A volumetric gear pump in a liquid pumps the
liquid. The gap volume must fall to about one molecule, or the carried water
must be returned, before the pump can deliver oxygen without water.

Limits. The wells come from a proxy, not from the diamond root site. The
search is rigid, has no solvent and ignores the entropy of the bound pose.
The free volume uses zero probe radius, so the water count is an upper-range
estimate.

## ADR-0090: The site-blocking rod pocket: the block works, the pocket is wet

Context. ADR-0089 showed that the gear mesh never touches a bound oxygen and
that a gear gap carries about 100 water molecules. The replacement idea is a
rotor pocket with an iron in the bore wall and a rod that slides past the iron.
The rod fills the space in front of the iron, so no bound oxygen fits there.
Two tests could reject it: the rod must block the site (test 1), and the loaded
pocket must stay dry in liquid water (test 2).

Test 1 (`rod_pocket` part, `rod_block_json`). The iron replaces a bore-wall
carbon, so it keeps three lattice bonds and sits level with the wall, not
standing out. `metal_site_binding` now takes the open side from all the iron
bonds and searches tilted donor positions too. With a 4.5e-10 m rod-to-wall gap
the rod does not block the site: a bound oxygen fits between the rod caps at
-74 to -80 kJ/mol. Pressing on it costs only a few kJ/mol. A 3.6e-10 m gap with
a 2.6e-10 m rod blocks the site, but the narrow bore crowds the loaded iron to
-40 kJ/mol, too weak to load. A 3.5e-10 m rod in a 7.1e-10 m bore (gap
3.6e-10 m) does both. Rod pulled back: oxygen -81 kJ/mol, water -15. Face
1.5e-10 m below the iron: oxygen -6. Face at or past the iron: no bound pose.
Quasi-static peak push 2.5e-9 N (estimate). The rod stays at least 2.5e-10 m
from the block.

Test 2 (`scripts/pocket_water_md.py`, OpenMM 8.6, SPC/E water, frozen solid,
UFF with Lorentz-Berthelot mixing, iron-water Morse from the proxy, 310 K,
1 ns after 50 ps). The pocket (the bore above the pulled-back rod, 1.37 nm^3
cylinder) holds 20.3 water molecules on average from a wet start and 19.7
from a dry start, 25 at most, and is never empty. A water touches the iron
49 to 66 percent of the time. Test 2 fails.

Finding. A pocket wide enough for a rod to slide into is too wide to stay
dry. Water fills a 14e-10 m open bore; a hydrophobic wall does not hold it
out at that size. The carried water per oxygen is the water that fits in the
pocket, so the pocket volume must be about one bound oxygen (about 5e-11 m^3
scale). The blocker must then be thin enough to enter a pocket of that size,
for example a carbon spike on a pin wheel, not a 3.5e-10 m rod.

Limits. Rigid, frozen solid; no polarization; mixing differs from UFF; the
iron-water and iron-oxygen wells come from the Fe(NH3)4 proxy. The results
are simulated estimates.

## ADR-0091: A snug rod pocket still carries about four water molecules

Context. ADR-0090 found that the 7.1e-10 m rod pocket holds about 20 water
molecules. M28-03 makes the pocket as small as the loading and the block allow.

Decision. The `rod_pocket` defaults are now a 2.6e-10 m rod, a 3.0e-10 m gap
(bore 5.6e-10 m) and a mouth 3.5e-10 m above the block origin (3.2e-10 m
above the iron). The rod is pulled back 1.5e-10 m below the iron to load.
The earlier -40 kJ/mol at the 6.2e-10 m bore was a lattice effect at that
radius, not a limit of narrow bores: the 5.6e-10 m bore loads well.

Test 1 (`rod_block_json`). Loaded oxygen -74.5 kJ/mol, water -16. Face
0.5e-10 m below the iron: oxygen -18. Face level with the iron or past it: no
bound oxygen pose. Peak push 2.4e-9 N (estimate). Rod to block at least
2.27e-10 m. A meter test now guards the load and the block.

Test 2 (`pocket_water_md.py`, 1 ns each). Empty site: 5.2 water molecules
from a wet start and 5.1 from a dry start, never empty. Oxygen bound: 3.8
water molecules beside it (at most 7). The large pocket held 20.

Finding. Shrinking the pocket cuts the carried water by a factor of five, but
not to zero. Water (2.8e-10 m) is smaller than oxygen, and a pocket that an
oxygen can enter and bind in keeps crescents beside the bound oxygen that
water fills. With rigid hydrogen-capped diamond walls the pump carries about
four water molecules for each oxygen. Into an oxygen tank, that is a failure:
the tank fills with water.

Option for the next step (not tested). Accept the water and separate it: the
rotor delivers oxygen and about four water molecules into a small collection
chamber, where oxygen leaves the liquid as gas; a second, dry stage feeds the
tank from that gas; a plain liquid pump returns the water to the blood. The
returned water carries dissolved oxygen back. With Henry's law at 310 K
(about 1.07e-3 mol/(L atm), an extrapolation at high pressure) the loss is
about 7e-5 of the delivered oxygen per atmosphere of chamber pressure: 0.07
percent at 10 atm. This is an estimate and needs its own kill tests.

Limits. As ADR-0090: frozen solid, UFF with Lorentz-Berthelot mixing, proxy
wells, no polarization. The count uses the bore radius at the carbon centres
up to the carbon plane of the mouth.

## ADR-0092: Option A, separate the water: the gas-side stage stays nearly dry

Context. ADR-0091 left the snug rod pocket carrying about four water
molecules with each oxygen. Option A accepts that water at the first stage.
Stage 1 pumps oxygen with its water from the blood into a small collection
chamber. In the chamber, oxygen leaves the liquid into a permanent gas space.
Stage 2 is a second rod-pocket rotor that loads from that gas and fills the
tank. A liquid pump returns the chamber water to the blood.

Kill test 1, thermodynamics. The chamber gas is saturated with water vapour
(6.28 kPa at 310 K), so its water chemical potential equals that of the
liquid. Stage 2 is dry only if a pocket that faces gas holds far less water
than a pocket in liquid. The difference comes from the missing bulk water at
the mouth, so it must be measured.

Kill test 2, MD (`pocket_water_md.py --environment gas`). The pocket is filled
from a 50 ps liquid run. Every water outside the bore is then made
non-interacting, and a water that rises 3e-10 m past the mouth is removed.
This makes a perfect sink. The pocket emptied in every run: after 76, 168 and
460 ps with the site empty, and after 54, 82 and 90 ps beside a bound oxygen.
The mean residence of one water is 38 to 494 ps. Kinetic theory gives
condensation onto the whole bore mouth at 0.22 per ns at most (sticking
coefficient 1). So a stage-2 pocket in saturated vapour holds at most about
0.1 water with the site empty and 0.03 beside a bound oxygen. In liquid it
held 3.8. The test passes: stage 2 carries about 1/100 of the water of
stage 1.

Mass balance (estimate). The returned water carries dissolved oxygen back
to the blood. Henry's law at 310 K gives 1.9e-5 oxygen per water per atm.
With four waters for each oxygen, the loss is 7.7e-5 per atm of chamber
pressure: 0.008 percent at 1 atm. Stage 2 then compresses to the tank. The
pressure-volume work against 1000 atm is about 2 kJ/mol, small beside the
80 kJ/mol release.

Open risks. (1) The gas space must be permanent. A meniscus pinned at a
2e-9 m hydrophobic slot holds about 2 gamma / w = 70 MPa (estimate). Start-up
needs a first bubble. (2) Water enters the tank at up to 0.03 per oxygen. It
builds up unless the release path removes it; with dense oxygen at 1000 atm
this is not modelled. (3) Every energy uses the Fe(NH3)4 proxy wells, UFF
walls, a frozen solid and no polarization. (4) The residence times come from
six short runs and are uncertain by a factor of a few. The conclusion holds
unless the uncertainty exceeds about a factor of 30.

## ADR-0093: Option B, a strain switch on the iron: large effect, weak evidence

Context. Option B releases oxygen by a change of chemistry, not by a push. A
cam strains the atoms that hold the iron so that the Fe-O2 bond changes
between the inlet and the outlet. Kill test: some strain must change the
binding by at least 45-50 kJ/mol.

Method. `scripts/fe_o2_switch_scan.py`, Fe(NH3)4 proxy, UKS PBE/def2-SVP,
rigid fragments, three frame spins, O2 at three distances (triplet) and one
quintet. Modes: the iron domed out of the N4 plane (0.3, 0.6 A), the Fe-N
bonds stretched (2.35, 2.55 A), and an NH3 pushed onto the back of the iron
(trans, 2.10, 1.95 A). 50 calculations, 3 workers, peak 4.6 GB. 23 of the 50
did not converge.

Results (binding enthalpy against the strained frame, kJ/mol; converged
calculations only where a pair exists):

| State | All points | Converged only |
|---|---|---|
| unstrained | -129 | -93 |
| dome 0.3 A | -135 | -93 |
| dome 0.6 A | -94 | no converged pair |
| stretch 2.35 A | -176 | -153 |
| stretch 2.55 A | -156 | -167 |
| trans NH3 2.10 A | -203 | -203 |
| trans NH3 1.95 A | -225 | -202 |

Finding. Stretch and trans strain make the bond much stronger: a trans
ligand adds about 110 kJ/mol against the unstrained site (converged points).
Doming weakens it by about 35 kJ/mol, and only at 0.6 A, where no pair
converged. So the switch works in the reverse sense: a weak resting site,
and a cam that pushes a ligand onto the back of the iron at the inlet to load,
then withdraws at the outlet to release. The amplitude passes the 45-50
kJ/mol test by about a factor of two.

Limits. The evidence is weak. Half the calculations did not converge. PBE is
known to overbind O2 to iron by tens of kJ/mol, so the absolute energies are
too negative. The unstrained site at -93 kJ/mol is too strong to release into
a tank, so a real design needs a weaker resting site, which this proxy does
not have. Rigid fragments, no relaxation of the O2 pose, no solvent. The next
test must use a better method (a hybrid functional or a spin-state-checked
calculation) on a site that can exist on diamond.

## ADR-0094: Hybrid functionals confirm a large strain switch on the iron

Context. ADR-0093 found a strain switch with PBE, but 23 of 50 calculations
did not converge and PBE overbinds O2 on iron. This rerun uses two hybrid
functionals and a stronger SCF procedure.

Method. `fe_o2_switch_scan.py --xc tpssh b3lyp`, Fe(NH3)4 proxy, def2-SVP,
grid level 2, two initial guesses per point, damping and a level shift,
then a clean pass from the result, lowest converged energy kept, <S^2>
recorded. States: unstrained, dome 0.6 A, trans NH3 at 2.10 and 1.95 A. O2
at 1.85 and 2.05 A (triplet) and 2.05 A (quintet). 50 calculations, all
converged; 3 workers, peak 3.9 GB, about 3 hours.

Results (binding enthalpy against the frame at the same strain, kJ/mol):

| State | TPSSh | B3LYP |
|---|---|---|
| unstrained | -158 | -136 |
| dome 0.6 A | -51 | -37 |
| trans 2.10 A | -199 | -163 |
| trans 1.95 A | -234 | -202 |

Doming by 0.6 A weakens the bond by 108 (TPSSh) and 99 (B3LYP) kJ/mol. A
trans ligand strengthens it by 28 to 75 kJ/mol. The two functionals agree
within 20 kJ/mol on every switch amplitude. The kill test (45-50 kJ/mol)
passes by a factor of two.

Strain cost. The empty site is lower in energy domed, by 39 (TPSSh) and 37
(B3LYP) kJ/mol, and it turns high-spin (S=2). The loaded site costs 69 and
61 kJ/mol to dome. This is the hemoglobin pattern: the empty iron rests
domed, the bound iron sits flat. So the cam flattens the empty site at the
inlet (costs about 38 kJ/mol) so oxygen loads, and lets it dome at the
outlet, pushing it with about 65 kJ/mol so that oxygen leaves.

Release against a tank. With the gas entropy (about 64 kJ/mol at 310 K,
ADR-0086) the domed site has a binding free energy of about +13 (TPSSh) and
+27 (B3LYP) kJ/mol at 1 bar. Release into a tank at 1000 atm needs more than
RT ln(1000) = 17.8 kJ/mol. That is marginal for TPSSh, so the release side
wants a deeper dome or a weaker resting site.

Limits. Rigid fragments, no geometry relaxation, O2 pose fixed at 60 degrees
tilt, no dispersion, def2-SVP (basis-set superposition error overbinds), the
triplet complexes have <S^2> near 3 (open-shell Fe(III)-superoxide
character, normal for Fe-O2 but a sign of multireference character). It is
the Fe(NH3)4 proxy, not a site built on diamond. The switch amplitudes are
estimates good to perhaps 30 kJ/mol; the absolute energies are less
certain.

## ADR-0095: Option A fails its gas-space test: water floods a diamond slot

Context. Option A (ADR-0092) needs a permanent gas space in the collection
chamber, held in a slot between hydrogen-capped diamond walls. Kill test: a
slot that starts dry must stay dry in contact with liquid water.

Method. `slot_json` writes two parallel 3.2 nm plates with reconstructed
(001) faces; `scripts/slot_water_md.py` puts them in SPC/E water with a free
liquid surface (near zero pressure), OpenMM, 310 K, 1 ns, dry and wet starts,
gaps 1.0, 1.5 and 2.0 nm between the facing carbon planes. One process,
340 MB.

Results (water in the slot, mean of the last 500 ps):

| Gap | Dry start | Wet start | Fill of a liquid-density slot |
|---|---|---|---|
| 1.0 nm | 171.7 | 172.5 | 0.75 |
| 1.5 nm | 299.1 | 298.5 | 0.80 |
| 2.0 nm | 440.0 | 441.6 | 0.85 |

Every dry slot filled within a few hundred ps to the same count as the wet
start. The fill is below 1 only because the layer next to each wall is
excluded. Hydrogen-capped diamond with UFF walls is wetted, so capillary
action floods the chamber. Option A fails as designed: its gas space cannot
form in a diamond slot.

What would be needed. A strongly water-repelling wall (for example a
fluorinated surface, contact angle above about 110 degrees) or a positive gas
pressure that holds the meniscus. Either must be tested the same way. The
measured wetting also depends on the carbon well depth of UFF, which is not
fitted to the diamond-water contact angle, so the model may overstate the
wetting. The experimental contact angle of hydrogenated diamond is near 90
degrees, where small slots are known to fill.

## ADR-0096: Squeezing the pocket does not remove the carried water

Context. A cascade for a pure oxygen tank first needs a step that removes the
water that stage 1 carries (3.8 per oxygen, ADR-0091). The cheapest candidate:
raise the rod to just below the bound oxygen before the pocket closes, so
the smaller pocket pushes water out of the mouth.

Test. `rod_block_json` with the rod face 1.0e-10 m below the iron (oxygen
still bound at -73 kJ/mol; at 0.5e-10 m it falls to -18). `pocket_water_md.py
--o2 bound`, 1 ns, two seeds: 3.60 and 3.71 water molecules in the pocket,
at most 7, empty 0.5-1 percent of the time. With the rod at 1.5e-10 m it was
3.8. The water sits beside and above the bound oxygen, toward the mouth, not
under it. Squeezing from below fails.

Consequence. A chain of liquid-phase stages does not purify either: every
stage that loads from liquid carries about four water molecules per oxygen,
whatever the oxygen concentration of its feed. The size relation points the
other way: water (kinetic diameter 2.65e-10 m) is smaller than oxygen
(3.46e-10 m). So water cannot be kept out by size, but it can be let out by
size. Zeolite 3A dries gases on this principle.

## ADR-0097: Drain channel: no size window, but the O2 loss is tolerable

Context. ADR-0096 proposed a tank that drains its carried water out through
channels too narrow for oxygen. Kill test M28-09.

Rigid scan (`channel_scan_json`, hydroxyl-lined through-channels in a 1 nm
diamond slab, every orientation). The lattice gives only discrete sizes. At a
carbon radius of 3.6e-10 m both water and oxygen are blocked; at 3.8e-10 m
(free radius 1.36e-10 m) both pass with no barrier. There is no size window:
end-on O2 is no wider than the oxygen atom of water. The kinetic diameters
(2.65 vs 3.46 A) come from rotating gas molecules and do not apply in a
channel that aligns them.

The real requirement is milder. Drained water carries the oxygen that is
dissolved in it. Plain water at 1000 atm holds about 1.9e-2 O2 per water
(Henry, extrapolated), so with four waters per pumped oxygen the drain would
return 7.7 percent of the oxygen. A channel that excludes dissolved oxygen by
dG divides that by exp(dG / RT).

MD with Widom insertion (`channel_water_md.py`, SPC/E, charged hydroxyl
lining, 310 K, 1 bar). Bulk excess chemical potential of O2: +9.0 kJ/mol in
three 1 ns runs, against about +8.4 measured, which checks the method. In the
channel the Widom average is not converged: at 4.5e-10 m it gave +16, +38
and +63 kJ/mol in three runs, so exclusion 7 to 54 kJ/mol; at 3.8e-10 m, +8.
Widom overestimates in a dense region, so the lowest value is the safest:
exclusion about 7 kJ/mol, a factor of about 15.

Consequence. The drain returns about 0.5 percent of the pumped oxygen at 1000
atm (estimate), and less at lower tank pressure. That is tolerable. The
remaining risk is not the loss while draining but the leak while idle: if
oxygen can pass the water in the channel, the full tank bleeds by diffusion
(estimate: about 5e6 O2/s per channel at 4.5e-10 m, which empties a
respirocyte-scale tank in minutes). Only a strictly single-file channel,
where oxygen cannot overtake water, or a valve that opens the drain only
while it drains, stops that. The 3.8e-10 m channel is single-file.

Next test: a pressure-driven MD run through the single-file channel with
oxygen at high concentration on the tank side, counting water and oxygen
transits directly, with and without flow.

## ADR-0098: The single-file drain passes water and holds oxygen (estimate)

Context. ADR-0097 left one risk for the drain: oxygen leaking through the
channel while the tank is idle. Test M28-10 on the hydroxyl-lined 3.8e-10 m
channel (free radius 1.36e-10 m), `scripts/drain_transit_md.py`.

Transit counts (2 ns each, SPC/E, about 1 M dissolved O2 on the tank side,
pressure by a constant force on each water, Zhu et al. 2002). The channel
holds 4.0 to 4.7 waters in single file.

| dP | Water up | Water down | Net water per ns | O2 transits |
|---|---|---|---|---|
| 0 | 3 | 9 | -3 | 0 |
| 200 MPa | 38 | 1 | 18.5 | 0 |
| 500 MPa | 76 | 0 | 38 | 0 |

The flow is roughly linear in pressure: about 9e9 water/s per channel at
100 MPa (1000 atm). A tank that takes 3e8 O2/s brings 1.2e9 water/s, so one
channel drains it. No oxygen crossed in 6 ns, but that bound is weak.

Umbrella sampling of one O2 along the axis (41 windows 0.06 nm apart, 150 ps
each, WHAM). The PMF inside the channel averages +8.9 kJ/mol, in line with
the Widom value of +8 (ADR-0097). The peak is 21 kJ/mol above the mean bulk
level. The two bulk ends differ by about 9 kJ/mol, which measures how far
the profile is from convergence. A solubility-diffusion estimate with an
assumed D of 1e-9 m^2/s gives an idle leak of 7.6e5 O2/s per channel at 1 M.
A tank of 3e9 O2 would lose half in about an hour. At one fill and discharge
per circulation (about a minute) that is about 1.5 percent per cycle. With
the 9 kJ/mol uncertainty the range is about 0.05 to 45 percent.

Self-sealing. The channel wets, so gas oxygen at 1000 atm cannot push the
water out of it: the capillary entry pressure 2 gamma / r is about 7e8 Pa
(continuum extrapolation). The blood side keeps it full.

Conclusion. The drain passes the transit test as an estimate: water flows,
no oxygen was seen to cross, and the estimated idle leak is tolerable at the
central value. Tightening the PMF (longer windows, both pull directions) is
needed before the number is used in a design budget.

## ADR-0099: The drain PMF is reproducible; idle leak about 1-4 % per cycle

Context. ADR-0098 had one 150 ps-per-window PMF whose two bulk ends differed
by 9 kJ/mol. M28-11 reruns it: 41 windows of 300 ps pulled from the tank side
(`--save fwd.npz`), 41 windows of 300 ps pulled back from the far side with
a new seed (`--reverse`), and WHAM over both (`--combine`).

Results. The forward and reverse profiles agree within about 3 kJ/mol at
every point, so the profile is reproducible. Peak 22.6 (forward), 26.2
(reverse), 24.4 kJ/mol (both) against the mean of the two ends. Channel
interior +10.6, +13.4, +12.0. The ends still differ by 6.7 to 9.2 kJ/mol in
both directions, so this is not sampling noise: the windows reach only
0.7 nm past the slab, and the far face holds an O2 adsorption well of about
-6 kJ/mol. The tank-side end sits about 3.5 kJ/mol above the mean.

Leak. Solubility-diffusion with an assumed D of 1e-9 m^2/s, 1 M O2 on the
tank side: 1.2e5 to 4.0e5 O2/s per channel (2.3e5 combined), and about 4x
more when the barrier is taken from the tank-side level. For the Freitas tank
(1.51e9 O2) that is about 1 to 4 percent of a full tank per 60 s cycle. The
leaked oxygen enters the blood, so it is lost to control, not to the body.

Limits. D in a single-file channel is assumed, not measured; the diffusion
coefficient of O2 in a single-file water column is likely lower, which lowers
the leak. The windows should reach 1.5 nm past each face for a true bulk
reference.

## ADR-0100: Budget of the drained-tank respirocyte against Freitas

Reference figures (Freitas, "Respirocytes", 1996-1999 revision on
foresight.org, via web.archive.org). O2 tank 1.19e-19 m^3 holding 1.51e9 O2
at 1000 atm. Fill or discharge in about 10 s. Circulation about 60 s. Sorting
rotors: 1e6 molecules/s and 1e-22 J/molecule sorted, up to 1e-19 J/molecule
against 30,000 atm. Glucose engines 3e-13 W each, 7.3e-12 W with all 12
running. Water ballast tank 6.07e9 molecules. Rotor turn rate and the 236x
derivation are not given there.

Budget of this design (every entry an estimate):

| Item | Value | Against Freitas |
|---|---|---|
| Fill rate | 1.5e8 O2/s | same |
| Water carried in | 4 per O2, 6.0e9 per fill | equals his whole ballast tank |
| Drain capacity, one channel | 9e9 water/s at 1000 atm; blood osmotic pull 7.6 atm at start | enough |
| O2 lost while draining | 0.5 % | none in his design |
| Idle leak, one channel | 1-4 % of the tank per 60 s cycle | none in his design |
| N2 carried in | 6e4 per fill (vent needed) | none in his design |
| Energy per O2, rod push | 1.45e-19 J | about his 1e-19 J head-pressure cost |
| Energy per O2, ideal switch | 5.0e-20 J (reversible compression 3.8e-20 plus water) | half his head-pressure cost |
| Peak power, 10 s fill, rod push | 2.2e-11 W | 3x his 7.3e-12 W maximum |
| Peak power, 10 s fill, ideal switch | 7.6e-12 W | equal to his maximum |
| Whole-body power, 5.36e12 devices, one fill per minute | 20 W (push), 7 W (switch) | inside his 17-340 W |

Observations. (1) The only budget line that fails is peak power with the
rod push: the push dissipates the full 80 kJ/mol binding energy on every
molecule. The strain switch, if its cam stores and returns the strain energy,
fits Freitas's power. (2) Freitas's quoted engine output (3e-13 W at 1e8
molecules/s, 3e-21 J/molecule) is below the reversible work to compress O2
from 0.13 to 1000 atm (3.8e-20 J/molecule). That is consistent only if the
device recovers the compression work at discharge; his pages do not say so.
(3) Start-up: a tank full of liquid holds only 7e7 O2 dissolved at 1000 atm,
5 percent of the design load, so an O2 gas phase must nucleate in the tank.
That has not been tested.

## ADR-0101: A heme-like Fe-N4 site switches O2 when the iron domes 0.3 A

Context. The strain switch (ADR-0094) had only the free Fe(NH3)4 proxy.
Diamond cannot form a flat Fe-N4 ring itself, but Fe-N4 sites exist in
carbon materials and in heme. The site is an idealized Fe-porphine held
rigid by the rotor, with an axial NH3 below the iron. The cam pulls the
axial ligand, so the iron domes h out of the ring plane away from the O2
(the hemoglobin T-state motion).

Method. `scripts/heme_switch_scan.py`, UKS PBE, density fitting, def2-SVP on
Fe, N and O, 6-31G on ring carbons, STO-3G on hydrogens, second-order SCF
finish, rigid geometry, O2 end-on bent 60 degrees at 1.8 and 2.0 A, site
spins S = 0, 1, 2, complex spins S = 0, 1. 22 calculations, all converged,
3 workers, peak 4.0 GB, about 2 hours. A single calculation peaks at 1.7 GB.

Results (kJ/mol).

| Dome h | E_bind | dG at 1 bar (+64 entropy) | Empty-site strain | Loaded-site strain |
|---|---|---|---|---|
| 0.0 A | -82.6 | -18.6 | 0 | 0 |
| 0.3 A | -32.5 | +31.5 | 7.9 | 58.0 |
| 0.5 A | +19.7 (unbound) | +83.7 | 56.0 | 158.3 |

The flat site loads at the lung: occupancy 0.994 at 0.13 bar. The site domed
0.3 A releases into a 1000 atm tank (dG +31.5 against the +17.8 needed).
Doming 0.3 A switches the bond by 50 kJ/mol, as the proxy predicted.

Cycle. Load flat; dome the loaded site (58 kJ/mol); O2 leaves into the tank;
flatten the empty site (returns 7.9). Cam work 50 kJ/mol against a
reversible minimum of 23 (0.13 to 1000 atm). Per O2 9.5e-20 J with the water
work; peak power for a 10 s fill 1.4e-11 W, 2x the Freitas maximum of
7.3e-12 W (the rod push was 3x). Mean over a 60 s cycle 2.4e-12 W, inside it.

Difference from the proxy. PBE puts the empty site in S = 1 at every dome,
so the empty site does not rest domed; the proxy with hybrids rested domed
and high-spin. PBE is known to over-stabilize low spin in iron, so a hybrid
rerun is the next check: a high-spin domed empty site would cut the cam work.

Limits. Rigid idealized porphine, two O2 distances only, mixed basis, PBE,
no dispersion. The iron-porphine is attached to diamond only by assumption;
the attachment chemistry is not modelled.

## ADR-0102: Tank start-up: keep a permanent gas cushion

Context. M28-13. A tank full of liquid holds only 7e7 O2 dissolved at 1000
atm (5 percent of capacity), so the gas phase must form.

Classical nucleation theory (water surface tension 0.070 N/m, 310 K, kinetic
prefactor 1e39 per m^3 per s, one event per second in the 1.19e-19 m^3 tank,
which allows a barrier of 46 kT) gives the supersaturation needed: 1,680 atm
homogeneous, 1,190 atm on a 90-degree wall, 850 atm at 110 degrees. The tank
works at 1000 atm, and ADR-0095 found capped diamond wetted, so a bubble does
not form reliably at operating pressure.

Decision. The tank is built gas-filled and never discharges below a cushion
of about 10 atm. The cushion costs 1 percent of capacity. A bubble that fills
the 0.61 um tank has a Laplace pressure of only 4.5 atm, so it is stable.
The drain keeps the liquid out while the tank pressure exceeds the blood
pressure. Loss of the cushion (a full discharge or a leak) is a failure mode
that the controller must prevent. This is an estimate; nucleation rates
carry orders-of-magnitude uncertainty, which does not change the decision.

## ADR-0103: B3LYP confirms the switch size; experiment fixes the absolute level

Context. M28-15 reruns the heme-like site of ADR-0101 with B3LYP at dome 0
and 0.3 A (15 calculations, all converged, 3 workers, peak 5.8 GB).

Results (kJ/mol). Flat: E_bind -19.9 (dG at 1 bar +44). Domed 0.3 A: +37.5
(unbound). Switch 57.4 (PBE: 50.1). Loaded-site strain 60.0 (PBE 58.0),
empty-site strain 2.6 (PBE 7.9). The empty domed site stays S = 1 in both
methods (S = 2 lies 13 above), so the hoped-for high-spin saving does not
appear, and the cam work stays about 50-57 kJ/mol.

Two methods, one conclusion. The switch amplitude agrees within 7 kJ/mol.
The absolute binding differs by 63 kJ/mol: B3LYP is known to underbind O2 to
heme, and generalized-gradient functionals are closer. Experiment sets the
level: myoglobin at 37 C has a P50 near 0.3-0.4 kPa, so its flat (R-like)
site has dG near -14 kJ/mol at 1 bar. PBE gives -18.6 for the flat site here.
Taking the experimental level and the computed switch: the flat site loads at
the lung, and the domed site sits near +36 to +43 kJ/mol, above the +17.8
needed to release into 1000 atm. The switch needs at least 32 kJ/mol; the
two methods give 50 and 57.

Hemoglobin's own T/R switch is about 12-15 kJ/mol, so this design forces
about four times more strain than hemoglobin uses. That forced strain is what
the cam supplies.

Memory. Hybrid runs now plan 2.8 GB per worker (2 workers on this machine).

## ADR-0104: The axial lever fails: the Fe-N bond is weaker than the dome

Context. M28-17. The switch of ADR-0101 needs a cam that domes the loaded
iron 0.3 A. The plan was to pull the axial ligand below the iron.

Force needed. The loaded site costs 58 (PBE) to 60 (B3LYP) kJ/mol to dome
0.3 A: a mean force of 3.2-3.3 nN and, for a harmonic mode, 6.4-6.6 nN at
0.3 A.

Force available (`heme_switch_scan.py --pull`, loaded site, O2 at 1.85 A,
flat ring, axial NH3 pulled from 2.1 to 4.0 A, PBE, 15 calculations, all
converged, peak 4.0 GB). Energy along the pull: 0, 9.3, 24.7, 39.9, 45.7
(spin change), 75.0, 88.1 kJ/mol at 2.1, 2.3, 2.5, 2.7, 3.0, 3.4, 4.0 A. The
steepest slope is about 1.3 nN. The bond holds 88 kJ/mol, more than the dome
needs, but it cannot hold the force: with the dome stiffness of about
1300 kJ/mol/A^2, 1.3 nN domes the iron only about 0.06 A before the ligand
lets go. The lever fails by a factor of 3 to 5.

Consistency with biology. Hemoglobin pulls on its iron through the same
single proximal Fe-N bond, and its T/R switch is only 12-15 kJ/mol. A single
dative bond limits how hard any mechanism can pull on the iron.

Next candidate. Put the load on covalent bonds. The cam squeezes the
porphyrin core through its meso anchors (C-C bonds carry 4-6 nN). A core
too small for Fe(II) pushes the iron out of the plane, and the Fe-N bonds
work in compression, which is stiff. Kill test: scan core size against dome
and O2 binding.

## ADR-0105: Squeezing the porphyrin core does not dome the loaded iron

Context. M28-18 tested the lever proposed in ADR-0104: squeeze the ring so a
core too small for the iron pushes it out of plane.

Method. `heme_switch_scan.py --cores 2.06 2.00 1.92 1.85 --domes 0 0.25 0.5`,
PBE, ring scaled uniformly, empty site S = 1, 2 and loaded site S = 0, 1,
O2 at 1.85 A. 49 calculations, all converged, peak 4.1 GB.

Results. At every core size the loaded iron is lowest flat (h = 0). The
ring absorbs the squeeze: a 4 percent smaller core costs 380 kJ/mol, 7.5
percent costs 1300, almost all in its own C-C and C-N bonds. Binding at
h = 0 hardly moves (-76, -81, -85, -90 kJ/mol from large to small core).
Doming relieves an Fe-N bond only by sqrt(core^2 + h^2), so a 1.85 A core
would need 0.76 A of dome to reach 2.00 A. The lever fails. A few points are
not monotonic (for example the empty 2.06 A site at h = 0.25) and probably
landed on excited SCF states; they do not change the conclusion.

Spin check (the question from the hemoglobin comparison). With the core
expanded to 2.06 A and the iron domed 0.5 A, the empty site is high-spin
(S = 2), as in deoxy-heme. At the 2.00 A core it stays S = 1, so the fixed
2.00 A core of the earlier scans biased the resting state.

## ADR-0106: Withdraw the axial ligand instead of forcing the dome

Context. The dome is too stiff to force (ADR-0104, ADR-0105), but the axial
ligand itself changes O2 binding (the trans effect, ADR-0094). A cam that
moves the ligand in and out does not fight the ring. Pulling the ligand off
is meant to break its weak bond; pushing it back is compression.

Method. `heme_switch_scan.py --pull 2.1 3.0 4.0 --empty` (empty flat site,
S = 0, 1, 2, 10 calculations, all converged) with the loaded pull of
ADR-0104.

Results (PBE, kJ/mol, relative to the ligand at 2.1 A):

| Axial Fe-N | Empty site | Loaded site | Switch in O2 binding |
|---|---|---|---|
| 2.1 A | 0 | 0 | 0 |
| 3.0 A | 19.6 | 45.7 | 26.2 weaker |
| 4.0 A | 45.9 | 88.1 | 42.1 weaker |

The four-coordinate empty site is S = 1, as measured for four-coordinate
Fe(II) porphyrins. Withdrawing the ligand 1.9 A weakens O2 binding by 42
kJ/mol. With the flat, ligated site at the myoglobin level (-14 kJ/mol at 1
bar), the withdrawn site sits near +28, above the +17.8 needed to release
into 1000 atm. The margin is about 10 kJ/mol, so the test passes narrowly.

Force and work. The cam pulls a ligand bond that holds about 1.3 nN, within
what a diamond arm carries; the bond breaking is intended. Cam work per
cycle = 88.1 - 45.9 = 42 kJ/mol, less than the forced dome (50-57). Per O2
8.2e-20 J with the water work; peak power at a 10 s fill 1.2e-11 W, 1.7x the
Freitas maximum (the dome was 2x, the rod 3x).

Limits. PBE only; rigid flat ring; NH3 as the ligand (histidine in biology);
the ligand must travel about 2 A. A B3LYP check follows.

## ADR-0107: B3LYP confirms the axial-withdrawal switch

M28-20. `heme_switch_scan.py --xc b3lyp --pull 2.1 4.0` for the empty and
the loaded site (12 calculations, all converged, 2 workers). Withdrawing the
axial NH3 from 2.1 to 4.0 A costs 28.5 kJ/mol on the empty site and 91.5 on
the loaded site: the O2 bond weakens by 63.0 kJ/mol (PBE 42.1). Both methods
exceed the 32 kJ/mol needed. As for the dome (ADR-0103), B3LYP gives the
larger switch; the two bracket it at 42-63 kJ/mol. Cam work equals the
switch, so B3LYP raises the work to 63 kJ/mol and the peak power to about
2.4x Freitas; PBE gives 1.7x. The empty site is S = 1 at both distances in
both methods.

Monitoring note. During this run a progress check fell back to PID 0 when
the scan had ended and summed every process on the machine (12 GB). The scan
itself peaked below 6 GB. Progress checks must treat "no scan process" as
finished, not measure a default PID.

## ADR-0108: Keep iron; oxidation is a lifetime problem, solved by a repair station

Context. M28-21. Fe(II) sites autoxidize to Fe(III), which cannot bind O2.
Hemoglobin loses about 3 percent per day and survives by enzymatic repair
(cytochrome b5 reductase). The question: does cobalt resist this better, and
what does repair cost?

Scan (`heme_switch_scan.py --redox Fe Co`, PBE, heme-like site, gas phase
plus C-PCM water, 26 calculations, 2 workers at 3 GB planned each). All 26
gas-phase calculations converged; 11 of the PCM single points did not, so the
PCM numbers only rank.

| | Fe | Co |
|---|---|---|
| O2 binding, NH3 at 2.1 A (gas) | -80.9 | -70.6 |
| O2 binding, NH3 withdrawn to 4.0 A | -38.8 | +24.8 (unbound) |
| Withdrawal switch | 42.1 (reproduces ADR-0106) | 95.4 |
| M(II)-O2 -> M(III) + O2-, relative, gas | 0 | -22 |
| same, relative, C-PCM water | 0 | -24 |

Cobalt. The large Co switch matches the known chemistry: cobalt porphyrins
need an axial base to bind O2. But cobalt is not easier to protect. In this
model its autoxidation is about 22-24 kJ/mol more favourable than iron's,
consistently in gas and in water, and its bound O2 is about 10 kJ/mol weaker.
Measured Co-heme proteins bind O2 more weakly than their Fe forms, so loading
at the lung would be marginal (occupancy about 0.66 at the PBE level). The
absolute autoxidation energies (about 270-740 kJ/mol) are not physical: the
charge-separated products are poorly described, and real autoxidation is
proton-coupled. Only the Fe-Co difference is used. Decision: keep iron.

Budget. At a 20 s fill a device needs about 900 sites (76 rotors of 12 at
1e6 O2/s; Freitas lists 5,184 O2 rotors). Oxidation per site does not depend
on the throughput. At a myoglobin-like rate (about 0.06 per hour) with each
site loaded half the time, a device suffers about 0.008 oxidations per
second: 3e-10 per delivered O2. It releases about 0.008 superoxide per
second, against about 94 per second from one red cell. Without repair, half
the sites die in about 23 hours; with the four waters beside the bound O2 the
rate may be 10-100x higher. So repair is mandatory, but the electrons are
trivial: 0.01-1 per second, about 1e-20 W.

Repair station. Each site passes a fixed point of the housing about 8e4
times per second. Put a one-electron donor there, on the arc where the site
is empty (after release, before loading), so the donor can only reduce
Fe(III) and never meets a bound O2. A donor near the bound O2 would itself
drive autoxidation. Moser-Dutton tunnelling (lambda 0.8 eV, driving force
0.3 eV) with a contact of 1 percent of a turn: donor-iron distance 22 A
repairs 0.1 per second per site, 18 A repairs 25 per second, against an
oxidation rate of 2e-5 to 2e-3 per second. The station is recharged from the
device power supply. The donor potential must lie below the Fe(III)/Fe(II)
couple of the site, not yet computed.

What is left. Reversible oxidation is handled on paper. Irreversible damage
is not: superoxide and hydrogen peroxide degrade the porphyrin ring (heme
degradation), and the device cannot repair a broken ring. That rate, and
whether the superoxide stays in the pocket water long enough to attack it,
set the device lifetime. Not modelled.

## ADR-0109: O2 loading rate, not the rotor, limits throughput

Context. Every earlier budget assumed, as Freitas does (1e6 molecules/s per
12-site rotor, about 8e4 rev/s), that a site loads on every pass. A heme
site loads by diffusion and binding: rate = k_on c.

Estimate. With a myoglobin-like k_on of about 1.5e7 per M per s
(literature order; our narrow pocket may be slower) and 130 uM dissolved O2
at arterial pO2 100 mmHg, an exposed empty site loads at about 2e3 per
second (time constant 0.5 ms). A Freitas-rate rotor gives each site about
4 us on a one-third inlet arc: 0.8 percent loading per pass.

Consequence. Throughput per site is at most about 650 O2/s with a one-third
inlet arc. A 20 s fill needs about 1.2e5 sites: about 9,700 rotors of 12
(Freitas lists 5,184 O2 rotors), about 2.5 percent of the device volume. The
rotors turn at about 300 rev/s for 90 percent loading per pass, 300x slower
than the Freitas rate, which also cuts friction losses. Cam work per O2 is
unchanged. Oxidation grows with the site count, to about 1 per second per
device (repair about 1.6e-19 W; superoxide 1 per second against 94 from a
red cell), still negligible (ADR-0108).

Limits. k_on in the snug diamond pocket is not computed; a slower on-rate
raises the site count in proportion.

## ADR-0110: O2 arrival at the snug pocket is fast; bond formation sets the rate

M28-24. `scripts/pocket_entry_md.py`: the snug pocket (5.6 A bore, rod pulled
back), frozen solid, SPC/E water, five mobile O2 at about 0.48 M, 2 ns per
seed, 310 K, one process of 0.22 GB. O2 is kept 0.5 nm below the free water
surface; without that wall it collected at the surface in some seeds.

Seeds (entries, arrivals within 4 A of the iron): (66, 41), (17, 0), (0, 0),
(42, 52). The spread is large: the hydrophobic O2 sits on the diamond faces,
and whether it sits near the mouth decides the seed. Pooled over 8 ns: 93
arrivals, k_arrival about 2.4e10 per M per s, about 3e6 per second at 130 uM.
Even the pooled lower bound from the slow seeds stays above 1e5 per second.

Conclusion. Arrival is about 1000x faster than the myoglobin-based load rate
of ADR-0109 (2e3 per second), so diffusion into the pocket does not limit
throughput. What limits it is bond formation: displacing the water on the
iron and forming the Fe-O2 bond, which the MD does not model. Myoglobin's
k_on (about 1.5e7 per M per s) includes that step, so ADR-0109 stays the
conservative design basis: about 1.2e5 sites. If bond formation in this
pocket is faster than in myoglobin, the site count falls in proportion.

## ADR-0111: The respirocyte O2 module as a relaxed CAD assembly

Context. The oxygen chain passes its simulated tests (ADR-0086 to ADR-0110).
This ADR builds the loading module as one audited, relaxed atomistic model.

Parts. `HemeRotorGenerator` (`crates/parts/src/heme_rotor.rs`): a diamond
disc (radius 4.0 nm, thickness 2.4 nm) with 12 radial pockets. Each pocket
has a ring slot, a heme-like site (Fe-porphine, Fe-N 2.00 A) normal to the
radius, a 5.6 A bore from the ring face to the rim, and a channel below the
ring for a ligand arm. The arm is a diamond rod along [111], turned onto the
pocket axis, so its end carbon bonds straight along the axis to an NH2 tip;
the nitrogen lone pair faces the iron. A [100] rod could not do this: its
end bond lies 55 degrees off the axis, and the NH2 hydrogens came within
1.85 A of the iron. `module_housing` is a ring with end plates and one
channel (inlet 345, outlet 165 degrees). Pockets 0-5 show the loaded state
(tip 2.1 A from the iron, O2 end-on at 1.85 A); 6-11 show release (tip at
4.0 A, empty).

Surface. The module parts pair surface dimers only up to 2.1 A
(`MODULE_MAX_DIMER_M`) and remove one carbon of any pair that still shares
a cap site (a CH2 step), repeating until no two caps are closer than 0.9 A.
With the global 2.6 A limit, closing 1.8-2.45 A step dimers bent the relaxed
rim by up to 5.6 A; with 1.9 A the removal eroded the disc and the rotor
collapsed in relaxation; 2.1 A settles in a few rounds.

Relaxation (`module_scene::build_relaxed_module_scene`). The housing relaxes
alone; the rotor, the 12 sites and the 12 arms relax together with the six
Fe-N ligand bonds included, so the pockets stay fitted. `relax_valence`
now keeps the built bond lengths and angles for any bond or angle with an
atom other than C and H, or with a three-bonded (sp2) carbon, so the
porphine keeps its quantum-chemistry geometry; a straight angle is skipped;
Lennard-Jones between four-bonded carbons of one lattice is excluded (it
double counts the crystal bond length). Converged: housing 867 iterations,
RMS 0.77 A; rotor unit 939 iterations, RMS 0.98 A, max 2.8 A. Fe-N 2.00 A
and Fe-O 1.85 A are kept; the rings buckle 0.11-0.18 A RMS.

Audit of the relaxed model (55,550 atoms, 26 bodies, 25 joints): no overlap,
no floating or disconnected body. Contacts under 2.2 A between any two
bodies, joints included: the six Fe-N ligand bonds (2.10 A) and one H-H pair
at 2.19 A. The unrelaxed model still has 12 cap pairs at 0.74 A on filled
notch carbons; the relaxation removes them.

Not in the model. The porphine-to-diamond attachment (a fixed joint stands
in for it), the cam that slides the arms, the repair donor, the tank and its
drain channels, and any solvent. The pocket above the ring (slot plus bore)
is larger than the 5.6 A test pocket of ADR-0091, so its water count must be
measured again. The viewer shows the module under Machines as
"Respirocyte O2 module (relaxed)"; it is built on request in about 50 s.

## ADR-0112: The as-built module pocket carries about 1-2 waters

Context. ADR-0091 measured the snug test pocket of ADR-0086 at 5.2 waters,
and ADR-0095 showed that diamond slots flood. The module pocket is neither:
it is the ring slot (0.90 nm radius, 0.72 nm along the axis) plus the bore
to the rim (0.56 nm radius) plus the arm channel behind it (0.81 nm radius,
2.24 nm long). This ADR measures the water carried by that geometry.

Method. `scripts/module_pocket_water.py` takes the relaxed module scene,
cuts a cylinder of the rotor around pocket k (radius 1.6 nm, from 2.9 nm
behind the iron to 1.0 nm past the rim), keeps that pocket's heme site and
ligand arm, and freezes every solid atom. SPC/E water fills a periodic box
around the piece. The iron keeps its Morse well for water (ADR-0086). Two
starts bracket a slow filling or drying: `wet` fills the pocket and the arm
at the start, `dry` leaves them empty. One job, about 9,000 particles,
below 0.5 GB.

Result (2 ns, 310 K). Loaded pocket (k=0): wet start 2.03 waters, dry start
0.99. Released pocket (k=6): wet 1.99, dry 1.56 (1 ns). The iron is never
wetted in any run (zero contacts under 3.2 A). So a pocket carries about one
or two waters, not the four assumed in ADR-0086. The two starts still differ
at 2 ns: the pocket is a tight hydrophobic cavity and exchanges water slowly.
The design should use two waters per pocket.

Carry is not a sizing problem. Twelve pockets at two waters per revolution,
at 300 rev/s for 20 s, give about 144,000 waters per fill. The tank holds
1.5e9 oxygen molecules. The water is 0.01 percent by count and 0.005 percent
by mass.

Open. The cut leaves the back of the arm channel open to the bath, so the
arm count (0 to 60 waters, depending on the start) is not physical. In the
device the cam seals that end. A sealed-back model is needed before the arm
water is quoted. The result is a simulated estimate: frozen solid, UFF
walls, and no water-carbon well fitted.

## ADR-0113: The porphine mount is a flexible two-carbon tether

Context. ADR-0111 left the porphine-to-diamond attachment as a placeholder
fixed joint. This ADR measures the mount geometry and designs the link.

Measured gap. In the as-built slot, the porphine meso carbons sit at radius
0.344 nm and the nearest rotor carbon (the slot wall) at 0.52 to 0.61 nm. The
gap is about 0.55 nm. A direct C-C bond (0.154 nm) and a single methylene
bridge (0.31 nm) are both too short.

A rigid aryl strap does not fit. A meso-phenyl reaches the wall, but its
ortho and meta hydrogens clash with the ring's outer beta hydrogens: the
closest H-H contact is 0.084 nm in a model build. Real meso-aryl porphyrins
twist the aryl ring out of the porphyrin plane to avoid this clash, and that
twist points the strap away from the wall, along the pocket axis. So an
in-plane meso-aryl mount is not usable here.

The mount. Each site carries four two-carbon tethers, one at each meso
carbon. Each tether is a chain of four C-C bonds (meso-C, C-C, C-C, C-wall)
built as a symmetric zigzag of 1.54 A bonds in the plane of the pocket
radius and the pocket axis, so its hydrogens point away from the ring. The
end carbon is bonded to the nearest rotor wall carbon. The chains bow
(about 32 degrees off the direct line) and relax to bonds of 1.42 to
1.55 A.

Result. The relaxed module (55,934 atoms) has no overlap and no loose body.
Contacts below 2.2 A between bodies are the twelve tether bonds (1.42 to
1.55 A) and the six Fe-N ligand bonds (2.10 A). The rotor unit converges in
2,481 iterations, RMS 1.2 A, max 5.1 A.

Limits. A four-bond chain is flexible: this holds the ring near its place
but does not fix its orientation. The production design should size the
pocket (a counterbore matched to the ring plus one bond) so a short rigid
link fits. The model draws the bond; it does not model the reaction that
forms it. This is a simulated estimate.

## ADR-0114: The tank wall carries hydroxyl-lined drain channels

Context. ADR-0096 to ADR-0100 proposed that the tank drains its carried
water through a single-file channel. The channel was tested as an isolated
lined slab (M28-09, M28-10). The channel was not part of the tank CAD. This
ADR adds it.

Design. The tank generator takes two new parameters. `drain_channel_count`
sets the number of channels (default 1). `drain_channel_radius_m` sets the
cut radius to the wall carbon centres (default 3.4e-10 m). Each channel is a
radial hole in the equatorial plane, spaced evenly in azimuth. The cut
removes every lattice carbon within the cut radius of the channel axis. A
free tetrahedral direction that probes into a channel takes a hydroxyl
group; every other free direction stays bare, as the shell was before.

Result. At the defaults the channel carries 72 lining oxygens. The nearest
atom to the axis is a hydrogen at 2.29e-10 m, so the free radius to atom
centres is about 2.3e-10 m and the accessible radius is about 1.2e-10 m.
That is close to the 1.36e-10 m channel of the MD test (ADR-0098). A
straight probe along the axis crosses the wall with more than 1.0e-10 m
clearance. A tank with zero channels carries no lining.

Limits. The channel is a geometric model. The flow and leak numbers stay
from ADR-0098, which is a simulated estimate. The lining groups carry model
charges (ADR-0063). The channel is equatorial in the model; the real device
routes the drain to the blood side. The wall is left bare away from the
channels, as before.

## ADR-0115: The assembled device scene

Context. The O2 module (ADR-0111), the gas tank with drain channels (ADR-0114)
and the pump were separate parts. This ADR puts them in one scene.

Design. Body 0 is the module housing, fixed. Body 1 is the rotor on a revolute
joint about z. Bodies 2 to 25 are the 12 heme sites and the 12 ligand arms.
Body 26 is the tank shell, which encloses the module: its inner radius is the
module outer radius plus 0.6e-9 m, and its wall is 0.6e-9 m. Bodies 27 and
28 are the pump cylinder and piston, parked above the module at z = 3.6e-9 m.
The tank and the pump are fixed to the housing. Each subsystem relaxes on its
own, because no two subsystems touch.

Result. The scene holds 122,893 atoms, 29 bodies and 28 joints. It has no
overlap, no loose body and no disjoint body. The only cross-body contacts are
the 18 inside the module: the 12 tethers and the 6 Fe-N bonds.

Limits. The layout is schematic. The tank and the pump are drawn at the scale
of the module, not at the scale of a real respirocyte hull. The pump seal ring
is left out: its radial thickness, 0.1e-9 m, is below one bond, so the
lattice cut is not one connected body. The scene builds in about 100 s in a
release build. This is a simulated estimate.

## ADR-0116: The arm cam is a closed groove in the housing end plates

Context. M28-29. Each ligand arm sat in a closed channel in the disc with
nothing to drive it. The Fe-N bond already joins each loaded arm to its iron
in the relaxation (`relax_module_parts`).

Design. Each arm carries a follower pin along z, a [1-10] diamond rod of the
arm radius: up on even pockets, down on odd ones, so the slots in the disc
skins stay apart. The pin passes through a radial slot in the skin and runs in
a groove cut in the housing end plate (`ArmCamTrack`). Both groove walls stand
at H-H contact from the measured pin carbons (`CAM_CONTACT_M`, 4.6 A between
carbon centres), so the cam drives the arm both ways. The arm grows from 1.2
to 1.4 nm so the pin slot clears the ring slot.

Finding: the groove, not the pin, sets the release distance. The Fe-N bond
still pulls the arm at 3-4 A (ADR-0106: 19.6 to 45.9 kJ/mol over 1 A, about
0.4 nN), so a released arm rests on the groove wall nearest its iron. Diamond
cut surfaces step by about 0.9 A and the step pattern changes with azimuth,
so the free travel toward the iron was 0.75-1.4 A; the nitrogen would have
rested at 2.6-3.25 A, where the switch is about 26 kJ/mol (32 needed). The
release track is therefore withdrawn by a 2.0 A margin.

Results. A 2 degree sweep over a full turn (`the_arm_cam_holds_every_arm_over
_a_turn`): on the loading arc the arm can always reach 2.1 A (free travel
toward the iron >= 0.2 A); on the return arc the nitrogen stays >= 4.9 A. The
relaxed module: 0 overlaps, no disjoint body, only the 18 intended contacts,
loaded Fe-N 2.10 A, released 5.4-7.1 A (the valence model has no Fe-N pull, so
released arms rest anywhere in the play). The relaxation runs in two stages:
the housing alone, then the rotor unit against the housing atoms within 1 nm
held fixed (`ValenceRelaxOptions::frozen_atoms`; pairs of two frozen atoms are
dropped). Relaxed apart, the pins bent 4 A into the groove walls. The face gap
grows from 4.6 to 6.0 A because the relaxed disc bulges about 1 A. Module
relax 149 s, device 144 s; the server build limit is now 600 s.

Limits. The cam phases (arm in at 345 deg, out at 165 deg) follow the old
housing, which is wrong (ADR-0117); M28-32 re-phases them to the hull-wall
seals. Cam stroke 3.9 A; cam work is not recomputed past 4.0 A. Simulated
estimate.

## ADR-0117: The device layout was wrong; the piston pump is removed

Context. The device scene of ADR-0115 put a gas-tank shell around the whole
O2 module. Every pocket then saw only tank gas: no pocket met blood, and the
housing inlet and outlet both opened into the tank. The drain channels sat in
that shell, not between the tank and the blood. Calling the layout schematic
hid a topological fault.

Decision. The hull wall is the boundary between blood and tank (Freitas). The
rotor sits in the wall with its axis parallel to the wall surface: half the
rim in blood, half in the tank, with the housing sealing the rim at the two
crossings. A pocket loads with the arm in on the blood side, crosses the seal
with the arm in, releases with the arm out in the tank, and returns through
the second seal with the arm out; the arm goes in inside that seal. The water
carried in leaves through OH-lined drain channels in the hull wall. A flat
wall slab is used: over a few nanometres the curvature of a 1 um hull is about
0.1 percent. Built in M28-32. Discharge at the tissues needs a second rotor
with reversed cam phases (M28-33).

The piston pump is removed from the device: no step of the drained-tank design
(ADR-0100) uses it.

## ADR-0118: The O2 rotor sits in a hull-wall section with contact seals

Context. ADR-0117 fixed the layout: the hull wall separates blood from tank,
and the rotor sits in it with its axis in the wall plane. M28-32 builds it and
asks whether it seals.

Design. `hull_wall` builds one fixed body: a 2.4 nm diamond slab on x = 0
(blood on +x, tank on -x), two end plates over the rotor faces carrying the
cam grooves, and one OH-lined drain bore through the slab beside the plates.
The slab closes the rim only where the rim crosses it: two seal arcs of +-15.6
degrees about 90 and 270 degrees, rim gap 4.2 A. The cam is re-phased to the
seals: the arm goes in at 270 degrees (ramp 258-282, inside the seal) and comes
out at 120 degrees (ramp 108-132, in the tank). Pockets sit at 15 + 30k
degrees, so none is mid-ramp, and the cam profile sets each state (7 loaded,
5 released). The separate device scene and the tank shell are removed; the
module scene is the section.

Seal audit. `scripts/seal_audit.py` voxelizes the relaxed scene (0.5 A),
blocks every voxel within a Bondi radius plus a probe radius of an atom,
masks the drain, and flood-fills from the blood face. The first section
leaked O2: the axis lies in the wall, so every closed path around it crosses
the wall, and the 7 A face gap and the cam grooves joined the blood window to
the tank window. Masking showed the disc interior (pockets, arm channels) was
already sealed. Fix: a face-seal ring on each plate at the outer edge of the
rotor face. A CSG face does not land on a lattice plane, so the first ring
left up to 6.4 A between carbons and passed a 2.0 A probe; the ring is now set
by measured layers, five (001) layers (4.46 A) beyond the outermost rotor
layer. The rotor and the wall share one lattice, so the gap is exact.

Results (relaxed). 77,353 atoms. No overlap, no loose or disjoint body. The
wall touches the rotor only on the face-seal rings, at 2.01-2.14 A. O2 probe
(1.73 A): no path from blood to tank, also with every radius scaled by 0.9.
Water probe (1.40 A): passes the rim seal; closed at 1.50 A. The as-built
drain bore does not pass a rigid water probe.

Limits. Static, rigid-sphere, one rotor pose. A path found is real; a pass is
necessary, not sufficient: leak rates need MD (M28-36). The drain's water
transport rests on ADR-0098's MD of a wider channel (M28-35). The face seals
rub at contact all round; their friction is not estimated. Simulated
estimate.

## ADR-0119: The as-built drain does not carry water; two MD faults fixed

Context. M28-35. ADR-0118 left the hull-wall drain unproven: a 3.4 A bore
whose hydroxyl hydrogens narrowed it to 1.6 A. The water transport rests on
ADR-0098's MD of a 3.8 A, 1 nm channel.

Geometry changes. The bore cut radius is 4.2 A, and the narrow section is
1.0 nm long, the length of the tested channel; 7 A counterbores lead to it
from both wall faces. Each hydroxyl hydrogen starts at its clearest turn
about C-O. Pointing it away from the bore had pressed it to 1.9 A from the
wall carbons. Lining oxygens sit 3.08 A from the axis (tested channel 2.91).
`drain_patch_json` cuts the as-built bore out of the wall for the MD.

Two faults in `drain_transit_md.py`, both fixed:
1. Under-filled box. The water grid floored its row count on each axis, so
   a 1.8 nm box lost up to 14 percent of its water. A vapour bubble formed
   at every pressure, and the push force moved it to the low-pressure face.
   The rows are now rounded and stretched to fill the box.
2. Low side below zero. The constant force splits dP about the mean
   pressure. The water is now packed denser by 1 + kappa dP / 2.
After the fixes the ADR-0098 channel carries 12-13.5 waters/ns at 200 MPa
(was 18.5), with no cavity. The 200 and 500 MPa numbers of ADR-0098 came
from a box with a bubble; they are now a rough estimate.

Results, as-built drain, 2 ns each:

| Lining | dP | Net water | Narrow section |
|---|---|---|---|
| frozen | 100 MPa | 0 (101 up, 101 down) | 3 fixed water sites, 3 A apart |
| frozen | 200 MPa | 0 (104 / 104) | same |
| mobile OH | 100 MPa | 0 (0 / 0) | dry |
| mobile OH | 200 MPa | 0 (14 / 12) | dry |
| mobile OH | 0 | 0 (2 / 2) | dry; no O2 crossed |

Mobile OH (`--mobile-oh`): each lining hydrogen keeps a stiff O-H bond and a
C-O-H angle and turns about C-O. It has mass 4 amu, and the step is 1 fs,
because a 2 fs step blew up after 50-100 ps. 113 of 134 hydroxyls turn; 21
start within 2 A of a carbon or a cap hydrogen and stay frozen.

Finding. Neither lining model carries water at the tank pressure. A frozen
lining traps single waters in fixed pockets. A mobile lining turns its
hydrogens into the bore and fills it. The tested channel conducted only with
a frozen lining: all 16 of its hydroxyls start near cap hydrogens, so none
turned in the control. No channel has yet been shown to conduct with a
lining that can move. The drain, which ADR-0100 relies on, is unproven.

Also: the rotor relaxation needed 3934 iterations after this change. The
budget is now 6000. The rim seal stops a probe of 1.55 A or more, so water
passes it (M28-36). The O2 seal holds at 1.0x and 0.9x radii.

Next (M28-38). Scan the drain with a mobile lining: lining density (every
cap, every other cap, H-capped), bore radius 3.8-5 A, and length 0.5-1 nm.
Measure water flux at 100 MPa and O2 entry at idle. Simulated estimates.

## ADR-0120: No drain channel; every bore tested fails one test

Context. M28-38. ADR-0119 found the OH-lined drain carries no water when its
hydroxyls can turn. This ADR screens the alternatives with
`drain_transit_md.py` on bores cut from the hull wall (`drain_patch_json`,
raw, 1.0 nm narrow section, 7 A counterbores), 0.5 ns at 200 MPa.

| Bore (cut radius) | Lining | Net water | O2 in the channel |
|---|---|---|---|
| 4.6 A | H-capped | 38 /ns | 0 % of samples |
| 4.2 A | H-capped | 36 /ns | 0 % |
| 3.8 A | H-capped | 18 /ns | 16 % |
| 3.6 A | H-capped | 10 /ns | 41 % |
| 3.4 A | H-capped | 0 | 0 % |
| 4.2 A | OH, every 2nd cap, mobile | 0 | 0 % |
| 4.6, 5.0 A | OH, mobile | 0 | 1-3 % |

A mobile hydroxyl lining always turns into the bore and blocks it. An
H-capped bore that carries water also takes up O2: the umbrella PMF of O2 in
the 3.6 A bore (40 windows, 150 ps, forward only, ends mismatched by 6 kJ/mol)
puts the channel 10 kJ/mol below bulk. That is an idle leak of about 2.4e8
O2/s per channel from a 1 M tank, the whole tank in seconds. Wider bores carry
water without O2 in 0.5 ns, but any continuous wet path also passes
dissolved O2 by diffusion: about D c A / L = 4e8 /s for a 1 nm pore 2.4 nm
long at 1 M (an estimate).

Decision. The hull wall has no drain (`HullWallSpec::drain` is false). The
tank side must not hold a continuous wet path to the blood. Water leaves in
sealed packets instead (ADR-0122). Simulated estimates.

## ADR-0121: The O2 rotor: arm collar, eight pockets, a 4 A bore

Context. The rigid-probe audit of ADR-0118 is blind to a path that a real
molecule finds through soft contacts. `seal_audit.py --energy` computes the
Lennard-Jones energy of a one-site O2 or water probe on a 0.5 A grid (UFF
solid, Lorentz-Berthelot) and bisects the lowest barrier on any blood-to-tank
path, or between two points (`--pair`).

Findings and fixes.
1. Every pocket opened into the rotor interior through its arm channel, and
   the interior joins all pockets: no barrier from blood to tank through the
   rotor. The arm now enters the ring slot through a collar, a contact-gap bore
   (arm radius + 4.6 A, 0.4 nm long). At 12 pockets the ring slots of
   neighbours overlapped, so the rotor has 8.
2. The bore from the iron to the rim shrinks from 5.6 to 4.0 A. A loaded pocket
   coming from blood carries 0.41 waters into the tank (1 ns MD; 1.38 with the
   5.6 A bore). A free O2 still reaches the iron with no barrier (the pair
   path from the mouth ends at -12 kJ/mol).
3. The bore above the iron does not connect to the slot below the ring (no
   path below 200 kJ/mol for O2 or water): the ring and its tethers close it.

Results (relaxed, 80,620 atoms). No overlap, no loose body. Blood to tank:
O2 barrier 146 kJ/mol, water 57 kJ/mol (22 kT). Rigid O2 and water probes
find no path. The ramps of the cam groove are now checked as well (nearest
wall atom 1.99 A), and the groove width is measured normal to its
centreline (ADR-0122). Simulated estimates; one rotor pose.

## ADR-0122: Water returns by a plunger rotor

Context. The rotor carries water into the tank and no drain works (ADR-0120).
A closed hydroxyl-lined cup carries about 18 waters and no O2 from O2-rich
water (pocket MD), but it carries the same 18 back from the blood: a cup moves
no net water. Net transport needs a pocket whose volume changes with phase.

Design (`water_rotor`, `WaterRotorSpec`). The disc of the O2 rotor with 6
radial bores, each holding a plunger on the arm-cam mechanism (ADR-0116). A
plunger is a [111] diamond rod, radius 2.6 A, 1.9 nm long, with a cam pin at
its back. In the tank sump the cam withdraws it (at 165 degrees) and the
0.6 nm pocket fills with sump water; it crosses the 270 degree seal full;
in the blood it is pushed flush (at 345 degrees) and the water leaves; it
crosses the 90 degree seal flush, with no volume. Ramps are +-20 degrees,
midway between pockets: a pin spans about +-11 degrees, and a pocket 3
degrees past a ramp touched it.

Faults found and fixed on the way. (1) The skin slot of a pin reached the
withdrawn pocket and opened it to the face cavity: a guard now requires the
slot to end 2 A behind the withdrawn tip. (2) A 4.6 A plunger gap leaked:
the [111] plunger and a bore along an arbitrary lattice direction do not
match, so the real gap varies round the bore. At 4.0 A the section seals.
(3) On a steep ramp a radial-width groove pinches a round pin; the groove
width is now measured normal to its centreline, and the cam sweep checks the
ramps (nearest wall atom 2.16 A over a full turn in 2 degree steps).

Results. Relaxed: no overlap, no loose body. Blood to tank: no O2 path below
200 kJ/mol, water 149 kJ/mol. Pocket MD as built (1 ns each): withdrawn in
tank-sump water with dissolved O2, 7.35 waters and no O2 in 1000 samples;
flush in blood, 1.42 waters counted in front of the tip. Net 5.9 waters per
pocket pass, tank to blood.

Limits and open points. The dissolved O2 in the sump run settled at about
0.03 M away from the walls (O2 adsorbs on diamond), not the intended 1 M; the
zero is at that concentration. The pocket must meet liquid: in 1000 atm gas a
lined cup took 4.7 O2. The sump therefore needs a hydrophilic wick between
tank gas and sump water, whose 1 nm pores hold liquid against gas up to
about 2 gamma / r = 2400 atm (Laplace, an estimate). The wick and sump are not
built. Simulated estimates.

## ADR-0123: The pump balance; the top stage is the problem

Context. M28-41. Does the O2 rotor (ADR-0121) raise O2 from blood to 1000
atm? `pocket_carry_md.py` puts one pocket of the relaxed scene in O2 gas and
counts the free O2 in the pocket; `pump_balance.py` combines the runs. The
O2 equation of state of this force field (barostat MD, 310 K): 3.1, 6.8 and
14.1 molecules/nm^3 at 100, 300 and 1000 atm.

Balance. One pass from a low side to a high side moves
net = 1 + n_loaded(low) - n_released(high). From blood (liquid, no free gas)
a released pocket carries back n = 0.106 rho (11 runs): a single stage stalls
near 550 atm. In a cascade the loaded pocket also carries gas forward.

Same-pocket comparison (0.6 ns each):

| Pocket | loaded, rho 9.5 | released, rho 14 | net 600->1000 atm |
|---|---|---|---|
| 0 | 0.43 | 2.51 | about -1.0 |
| 1 | 0.39 | 1.46 | about -0.1 |
| 5 | 1.18 | 1.53 | about +0.6 |

Finding. Below about 500 atm every stage pumps (0.6-0.8 net per pass). At
the top stage the cavity that the released O2 leaves fills with free gas,
and the net per pass is near zero or negative and depends on each pocket's
atomic surface. Single short runs, several restarted after NaN starts; the
trend is clear, the numbers are not converged. Simulated estimates.

## ADR-0124: A real-size system model fed by the atomistic parameters

Context. M29-01. The atomistic work answers questions about one rotor. The
device questions (how long to fill, how much power, how many rotors, does a
cascade reach 1000 atm) need a model at the next level of abstraction, one
that behaves like a CAD package and can be changed in seconds.

Decision. `python/nanocad/respirocyte_system.py`: components with
parameters, each parameter with a value, a range, a unit, a tag (measured,
estimated, assumed, reference) and a source. A 1 um hull (Freitas); rotor
banks as counts of measured sections (13.7 x 6.9 nm on the hull, 10.5 nm
deep); a cascade of O2 rotor stages with small intermediate chambers; a
water-return plunger bank; a discharge bank; a power supply whose peak caps
the rotor speeds (duty cycle). Pocket balance: net = occupancy + s_fwd
rho_low - s_back rho_high, with the force-field O2 equation of state. The
pocket fits are refit with all runs: s_back 0.109 (range 0.076-0.142, 16
runs), s_fwd 0.081 (0.033-0.141, 11 runs); the ranges are the per-run
scatter, which is mostly pocket-to-pocket. Binding from gas depends on the
gas density (assumed collision-limited), so an empty chamber does not feed
the next stage. Backward Euler with Newton steps, no external solver, so the
server runs it on any Python. The page `/system` edits and draws it.

Results (default device: 10,000 loading, 6,000 per cascade stage, 2,000
water and 4,000 discharge sections; 48 % of the hull; Freitas power).

| Stages | central | pessimistic chemistry, Freitas power |
|---|---|---|
| 1 | stalls at 470 atm | stalls at 100 atm |
| 2 | 1000 atm in 139 s | stalls at 326 atm |
| 3 | 1000 atm in 111 s | stalls at 419 atm |
| 4 | 1000 atm in 139 s | 435 atm at 600 s |
| 6 | 1000 atm in 212 s | 437 atm at 600 s |

Findings. (1) A cascade answers the top-stage problem only if the free-gas
carry falls in the favourable half of its scatter: the limit of many close
stages is net = 1 - (s_back - s_fwd) rho, +0.61 per pass at 1000 atm in the
central case and -0.53 in the pessimistic one. (2) The fill is
power-limited: the design speeds ask 1.7e-11 W, against Freitas's 7.3e-12
W, so the rotors run at about half speed; a 10 s fill would need about 10
times his power. (3) Three stages are best: more stages add pocket passes and
power. Energy per delivered O2 5.5e-19 J, about 4 times the single-pass
figure of ADR-0100. (4) The water bank keeps the chamber-1 sump below 1e4
waters; the discharge runs in 28-34 s.

Limits. Lumped, well-mixed chambers; no transport inside the device; no
lungs or tissue model (blood is an infinite reservoir at fixed state); many
assumed values, tagged. Simulated estimates.


## ADR-0125: The respirocyte is sized by the circulation, not by the tank

Context. M29-04. Freitas fixed the hull at 1 um and the tank at 1000 atm and
asked whether such a device could be built. That leaves the actual design
variables free: the diameter, the tank pressure, the rotor count and the rotor
speed. The question is what values replace the red-cell gas capacity.

Decision. `python/nanocad/respirocyte_circulation.py`. The blood is a loop of
`blood_volume / cardiac_output` (60 s at rest, 12 s at exercise). A device
loads O2 during the pulmonary capillary transit `t_lung` and releases it during
the tissue transit `t_tissue` (0.75 s at rest, 0.30 s at exercise). The target
is the arterio-venous O2 difference, 50 mL per litre at rest and 150 at peak
exercise; the red-cell capacity is 201 mL per litre (15 g/dL Hb, 1.34 mL/g).
A hull of diameter D has `pi D^2` of area, a measured rotor section takes
`13.7 x 6.9 nm`, and the tank must hold one cycle of O2 inside the rotor shell.
The engines supply `2.32 W/m^2` of peak mechanical power (Freitas, 12 engines
at 1 um). Because each bank uses the whole surface power at its own time, the
rotor speed is set by the power and the throughput is
`P_peak * occupancy / work_per_o2`, independent of the rotor count. The work
per O2 is the cam work plus the compression work from the blood partial
pressure through the O2 equation of state (ADR-0123).

Results (lowest feasible pressure for each diameter; one stage, which stalls
at 224 atm; the tank never needs a cascade):

| D (um) | P (atm) | sections | load rotors | rev/s | nanocrit rest | nanocrit exercise | chem. power exercise |
|---|---|---|---|---|---|---|---|
| 0.10 | 40 | 199 | 90 | 973 | 0.24 % | 0.71 % | 630 W |
| 0.20 | 15 | 798 | 359 | 1004 | 0.47 % | 1.4 % | 619 W |
| 0.30 | 10 | 1,795 | 808 | 1017 | 0.70 % | 2.1 % | 615 W |
| 0.50 | 10 | 4,985 | 2,243 | 1017 | 1.2 % | 3.5 % | 615 W |
| 1.00 | 10 | 19,940 | 8,973 | 1017 | 2.3 % | 7.0 % | 615 W |
| 2.00 | 10 | 79,761 | 35,893 | 1017 | 4.6 % | 14 % | 615 W |

The lowest nanocrit is 0.39 % at a 0.05 um device (50 sections, tank-limited),
but that needs about 1e17 devices in the body. The pressure optimum is the
cushion, 10-40 atm, not 1000 atm.

Findings. (1) The device is rotor-area and power limited, never tank-volume
limited at a practical size, so the pressure should be low. 1000 atm adds
compression work and, above about 200 atm, a cascade, for no gain. (2) The
rotor speed is about 1000 rev/s (60,000 rpm) at every size, set by the power
density and the section footprint. (3) The nanocrit is 0.4-14 %, far below the
45 % hematocrit; Freitas's nanocrit limit is not binding. (4) The binding
constraint is energy: moving each O2 with the cam costs 650 W of chemical
power at peak exercise, about half the whole-body exercise metabolism. A
reversible cam (23 kJ/mol) costs 342 W; compression alone, with a passive
release, costs about 75 W. The cam irreversibility, not the hull strength, is
the real limit.

Limits. Lumped, well-mixed; a single transit window per bed; the power density
and the cam work are Freitas's and the atomistic estimates; the device count at
small diameters is not bounded by clearance or fabrication. Simulated
estimates.

## ADR-0126: A passive carrier beats the pump for gas transport

Context. M29-04 sized the active respirocyte and found the cam, not the hull or
the nanocrit, is the binding cost: 650 W of chemical power at peak exercise,
342 W even with a reversible switch. The red cell does not pump. It carries O2
on a site that loads in the lung and releases in the tissue by the pO2
difference alone. So the question is whether the pump is needed at all.

Decision. `python/nanocad/respirocyte_design.py`. Two architecture families,
compared on the two currencies that matter: the volume fraction (nanocrit) and
the metabolic power at rest and at peak exercise. A carrier has a Hill O2
binding curve `S(p) = p^n / (p50^n + p^n)`; the delivered O2 per pass is
`capacity * (S(lung) - S(tissue))`; the carrier mass per litre is the red-cell
flux divided by that; the nanocrit is the mass over the density; the power is
zero. A pump reuses `design_point` with the cam work as the only variable
(50, 23, 0 kJ/mol).

Results (target: 50 mL O2/L rest, 150 mL O2/L exercise):

| architecture | nanocrit | chem. W rest | chem. W exercise | provenance |
|---|---|---|---|---|
| heme, bare | 0.42 % | 0 | 0 | estimated |
| dendrimer-heme | 1.0 % | 0 | 0 | assumed |
| heme-cyclodextrin | 2.8 % | 0 | 0 | estimated |
| hemoglobin | 12.3 % | 0 | 0 | reference |
| perfluorocarbon | 316 % | 0 | 0 | reference |
| pump, ideal (no cam) | 2.6 % | 4.8 | 72 | measured |
| pump, reversible cam | 3.0 % | 22 | 329 | measured |
| pump, cam 50 kJ/mol | 3.5 % | 42 | 630 | measured |

Findings. (1) The carrier model reproduces the real red cell: free hemoglobin
needs a 12.3 % volume fraction to carry the flux, close to the real 15 g/dL.
(2) A passive carrier with a red-cell-like site (p50 about 27 mmHg,
cooperative) meets the flux at 1-3 % nanocrit and no engine, so it beats the
pump on both nanocrit and power. (3) A perfluorocarbon cannot do it at room
air, because it is linear and needs a high pO2. (4) The pump is only needed
when O2 must be stored or delivered against a gradient. Its whole cost is the
cam: without the cam it fits the 100 W budget, with the cam it does not.

Limits. The carrier p50 is assumed tunable to a red-cell-like value; a real
heme needs a protective pocket and a large enough scaffold to stay in blood.
The carrier is a transport molecule, not a reservoir, so it holds no O2 beyond
the equilibrium. The pump numbers inherit the M29-04 limits. Simulated
estimates.

## ADR-0127: Device specifications, counts, and the site lifetime

Context. M29-05 chose the passive carrier over the pump. The remaining
questions are buildable numbers: how many devices, how much material, how much
O2 per volume, and how long the binding sites last. The site lifetime is the
least certain number in the program. A heme loses capacity when its Fe(II)
oxidizes to Fe(III). A red cell lasts 120 days because it carries a
methemoglobin reductase; a synthetic carrier has no repair unless one is
engineered in.

Decision. `python/nanocad/respirocyte_specs.py`. It turns each design into a
specification sheet and models the site lifetime from the autoxidation rate
``k``: half-life ``ln2 / k``, daily replacement ``mass * (1 - 2^(-1 day /
t_half))``. The rate is a parameter with four measured reference points: free
heme, hemoglobin without a reductase, myoglobin, and an engineered pocket. The
pump lifetime is the shorter of the rotor wear (an assumed 1e11 cycles) and
the heme autoxidation.

Specifications (target: the red-cell flux, 50 mL O2/L rest, 150 exercise):

| design | volume | devices/L | mass body | O2 per L of device | power body | site life |
|---|---|---|---|---|---|---|
| red cell | 45 % | 5.0e12 | 2.5e13 cells | 447 L | 0 | 120 d |
| heme, bare | 0.42 % | 6.1e21 | 31 g | 48,000 L | 0 | rate-set |
| dendrimer-heme | 1.0 % | 8.8e20 | 73 g | 19,400 L | 0 | rate-set |
| heme-cyclodextrin | 2.8 % | 8.3e21 | 200 g | 7,300 L | 0 | rate-set |
| hemoglobin | 12 % | 1.5e21 | 820 g | 1,600 L | 0 | rate-set |
| pump, cam 50 | 3.5 % | 5.4e14 | n/a | 5,700 L | 630 W | rate-set |
| pump, ideal | 2.6 % | 4.0e14 | n/a | 7,700 L | 72 W | rate-set |

Site lifetime against the autoxidation rate (dendrimer-heme, 73 g body):

| case | k (1/s) | life | replacement |
|---|---|---|---|
| free heme | 1e-3 | 12 min | 73 g/day |
| hemoglobin, no reductase | 1.4e-5 | 14 h | 51 g/day |
| myoglobin | 1.2e-7 | 67 days | 0.75 g/day |
| engineered pocket | 1e-8 | 800 days | 0.06 g/day |

Findings. (1) The carrier is 4 to 100 times denser in O2 than a red cell, so
it carries the same gas flux at 1-3 % by volume instead of 45 %. (2) The pump
needs 1e14 to 1e15 devices, each a machine; the carrier needs 1e20 to 1e21
molecules, which bulk chemistry can make. (3) The site lifetime is not known.
It spans 12 minutes to 800 days across the measured references, and the
synthetic pocket is unmeasured. The engineering target is k below about
1e-6 /s (life over 8 days) for a practical daily dose; a myoglobin-like pocket
gives 67 days.

Limits. The autoxidation rates are for other environments; the pocket rate is
the open question (the same risk as the ring damage in M28-22). The pump mass
and the rotor wear are not measured. Simulated estimates.

## ADR-0128: A self-contained 3D model page for the respirocyte

Context. M29-06 gives the device numbers. A colleague or a friend cannot read
a table of nanocrit and rotor counts. The program needs one artifact that
shows the design and its scale in a form that travels on its own.

Decision. `site/models.html`, linked as "3D models". It is one HTML file with a
small software 3D renderer: it projects triangle meshes to Canvas 2D. Opaque
parts rasterize through a depth buffer, so their surfaces have no depth-sort
artifacts; translucent and wire parts use the painter's algorithm. Curved
meshes carry per-vertex normals and shade smoothly (diffuse plus specular). It
has no dependency and no network call, so the file opens from a disk or an
email. The
shapes are schematic (a hull sphere, rotor rings, a tank, a biconcave red
cell, a carrier molecule), not atomistic, because a 1 um hull and an 8 um cell
are far too large for atoms. Three scenes: (1) the respirocyte next to a red
blood cell at true relative scale, with a close up inset; (2) the passive
oxygen rotor in detail; (3) the passive carrier molecule. Drag rotates, scroll
zooms, a button toggles auto-spin, and a scale bar and a key explain each
scene.

Findings. The true-scale scene makes the point: the 1 um device is one eighth
of the red cell across and about one five-hundredth of its volume. The close
up and the rotor scene carry the parts the scale hides. The page is a
schematic, not a physical render, and it says so.

Limits. The shapes are schematic and not to the CAD geometry; the renderer has
no shadows or reflections. It is a communication artifact, not a design tool.
Simulated.

## ADR-0129: The passive carrier site as CAD

Context. ADR-0126 and ADR-0127 choose a passive heme carrier: no engine and no
cam. The iron carries the oxygen and the oxygen partial pressure releases it.
M29-07 shows a schematic of it. The program needs the site as real geometry, so
a later stage can test the pocket, the mouth, and the guard.

Decision. `CarrierSiteGenerator` (id `carrier_site`) in `heme_rotor.rs` builds
one carrier site from a `CarrierSiteSpec`. It cuts a diamond slab with a ring
slot (radius 7.0 A), an O2 pocket above the slot (radius 3.6 A, height 3.5 A),
and a mouth in the top face (radius 2.6 A). It drops the existing heme site
(the Fe-porphine with its four two-carbon tethers) into the slot through the
shared `build_heme_site` helper. It then turns the pocket-wall cap carbon
nearest the bound O2 into a distal hydroxyl guard.

The guard is a real C-O-H group: C-O 1.43 A, O-H 0.96 A, C-O-H 108.5 degrees.
The oxygen is placed so that the O-H bond points at the O2. This makes the
angle at the oxygen between the carbon and the O2 equal to 108.5 degrees, which
fixes the O...O2 distance from the carbon...O2 distance. The generator picks
the wall carbon that gives an O...O2 distance nearest 2.9 A.

Findings. The default site gives O(guard)...O(O2) 2.99 A and H...O2 2.03 A, a
normal hydrogen bond. The mouth (2.6 A radius) is wider than the O2 kinetic
radius (1.73 A), so O2 can enter. The merged scaffold-plus-site part is 1252
atoms. The site is registered, so an agent builds it by name. The pocket is
snug: it is a pocket below the ring, not a wide cavity.

Limits. The guard is one static conformer; the real distal group moves. The
pocket-water MD kill test is not run yet, so the water exclusion and the
autoxidation rate are not measured here. The ring-to-diamond joint is still a
placeholder. Simulated geometry.

## ADR-0130: The carrier pocket excludes water

Context. The rod pocket of ADR-0090 held about 20 water molecules and failed
the dry-pocket test. A carrier must keep water off the iron, because water and
its ions drive the autoxidation of the heme.

Decision. Test the carrier pocket with the explicit-water model of ADR-0090
(`scripts/pocket_water_md.py`). Two small changes make it fit the carrier: the
geometry JSON now carries an explicit pocket `axis` (the carrier iron has four
nitrogen neighbours, not three carbons, so the test cannot infer the open
side), and the UFF table gains nitrogen. A new exporter
(`crates/parts/examples/carrier_site_json.rs`) writes the merged scaffold and
site.

Findings. The pocket stays dry. Over 500 ps in liquid water, with both a dry
start and a wet start, the pocket holds 0 waters and the iron makes 0 water
contacts (573 waters in the box). The pocket cylinder is 0.224 nm3, but the
bound O2, the distal guard, and the pocket wall leave under 0.01 nm3 of volume
with a 0.26 nm water clearance, so less than one water molecule fits. The
mouth (2.6 A radius) still passes O2, whose kinetic radius is 1.73 A. The
carrier is the opposite of the rod pocket: the pocket is snug, not a cavity.

Limits. The wet start could not place a water molecule in the pocket, so it
equals the dry start; the result rests on the free-volume measurement, not on
an observed drying. The autoxidation rate is not measured, because a classical
MD model cannot measure it. The guard is one static conformer. Simulated.

## ADR-0131: The distal guard strengthens the O2 bond (QM)

Context. ADR-0129 builds a distal hydroxyl guard on the carrier pocket wall,
2.99 A from the bound O2. ADR-0130 shows the pocket excludes water but cannot
measure the guard. The guard must hold the oxygen: a stronger O2 bond has a
higher barrier to the superoxide release that autoxidizes the iron. M29-10
tests the guard with quantum chemistry.

Method. `scripts/heme_guard_scan.py`, UKS PBE/def2-SVP with density fitting,
grid level 1, no dispersion, a damped first pass and a Newton finish. The model
is the shared heme site of `heme_switch_scan.py` (Fe-porphine, one axial NH3)
with O2 at 1.85 A. The guard is a water molecule whose first O-H points at the
distal O2 oxygen, at the ADR-0129 geometry (O...O2 2.99 A, H...O2 2.03 A).
Eleven points: O2, water, the bare site (S=0,2,4), the site with water (S=0,2),
the O2 complex (S=0,2), and the complex with the guard (S=0,2). Three workers,
planned cap 6.4 GB; peak resident 3.1 GB. All points converged.

Findings. The guard makes a hydrogen bond of -21.8 kJ/mol to the complex. The
same water binds the bare site by only -6.9 kJ/mol, so most of the bond is to
the O2. The O2 binding energy is -80.9 kJ/mol without the guard and -95.8 kJ/mol
with it: the guard strengthens the O2 bond by 14.9 kJ/mol. This is the size of
the distal-histidine effect in myoglobin, where the same H-bond slows
autoxidation by one to two orders of magnitude. The guard also polarizes the
bound O2 (Mulliken charge -0.18 to -0.22; triplet spin on O2 1.24 to 1.15). The
singlet and the triplet stay nearly degenerate (0.4 kJ/mol apart without the
guard, 2.4 with it), so the spin state is not resolved.

Decision. Keep the guard. It forms a normal H-bond to the distal O2 and
strengthens the O2 bond by about 15 kJ/mol, the mechanism that a distal H-bond
uses to slow autoxidation. Together with the dry pocket of ADR-0130, this is
the second protection of the carrier site.

Limits. A proxy, not the diamond site: the guard is a water molecule, not the
real C-O-H group; the geometry is rigid; PBE overbinds O2 on iron by tens of
kJ/mol, so the absolute energies are too negative. The calculation gives the
O2 bond strength, not the autoxidation rate. A stronger O2 bond also makes the
tissue release harder, so the resting-site affinity must absorb the change.
Simulated.

## ADR-0132: The circulation cycle limits the device to a flow-through pump

Context. The M29-01 time course (ADR-0124) filled the tank from an infinite
blood reservoir and then discharged it. The real device rides the blood: it
loads for the pulmonary capillary transit and discharges for the systemic
transit, and it stores nothing in between. M29-02 replaces the reservoir with
that cycle and asks what the tank and the delivery really are.

Method. `cycle_steady_state` in `python/nanocad/respirocyte_circulation.py`
reuses the compartment rates and the backward-Euler stepper of
`respirocyte_system.py`. Each loop the device loads for `t_lung` (0.75 s at
rest) and discharges for `t_tissue` (0.75 s); the rest of the 60 s loop is
arterial and venous storage, with no transfer. The loop repeats until the tank
reaches a periodic steady state. Tests in
`python/tests/test_respirocyte_circulation.py`.

Findings. The loop converges in about six loops (drift < 0.001 percent). The
tank swings 0.02 to 3.7 atm, far below the 1000 atm target. The 0.75 s lung
window loads 1.28e7 O2, which fills the tank to only 3.7 atm. To reach 1000 atm
would need 1.67e9 O2, about 130 loops (2.2 h), and the discharge empties the
tank every loop, so it never charges. The device delivers 1.28e7 O2 per loop,
set by the loading rate; the load and the delivery balance. The tank drains
below the 10 atm cushion every loop, so the no-nucleation rule of ADR-0102
cannot be kept in a real cycle. The loading is power limited: the two cascade
banks together demand more than the 7.3 pW peak of ADR-0100, so every rotor
runs at about 73 percent. The default 1 um device needs 5.5 percent nanocrit at
rest and 16 percent at peak exercise to meet the red-cell target.

Decision. Model the circulation as a real cycle and report the per-loop
delivery, not the tank capacity. The Freitas 1000 atm tank is a capacity
ceiling that the pulmonary transit cannot fill; the device is a flow-through
pump near ambient pressure. The M29-04 top-down design already agrees, because
it also uses `t_lung` as the loading window, and its optimum is a small device
at low pressure.

Limits. The discharge has no tissue-demand feedback, so the tank drains to the
blood pressure. The blood pO2 is held at its arterial and venous values. The
rotor speeds are the design speeds. The cascade banks are over-provisioned, so
they waste the fill power; a matched design would load faster. Simulated.

## ADR-0133: One device, every scale

Context. The site had a page for each scale: the atomistic viewer, the system
model and the 3D models page. The user asked for one view that zooms from the
whole respirocyte to the atomistic detail without a page change.

Method. `site/device.html` and `site/device.js`, one canvas, no dependency.
Two data sources: `/api/system` (the L4 layout: hull, bands, tank, chambers)
and `/api/part?id=carrier_site` (the L1 CAD, 1252 atoms). One log-scale camera,
`logWidth` from -4.7 m (20 um) to -9.7 m (0.2 nm). The hull is a sphere with
the three bank bands plus a shell for the unmodelled south; the carrier site
sits at the north pole. The camera pans from the hull centre to the site
between -5.2 and -7.0. The hull cross-fades out (-8.4 to -6.4) as the site
cross-fades in (-6.2 to -7.4). Atoms are ball-and-stick; bonds are computed by
distance (1.75 A cutoff). A scale bar, a level label and an element legend
track the view. Drag rotates, scroll or the buttons zoom, and the two buttons
jump to either end. Linked from every page header.

Decision. Keep one device page and link it from the atomistic viewer, the
system model and the 3D models page. It is the single place to see the
respirocyte and its carrier site together.

Limits. The hull is schematic: three coloured bands and a shell, not the
33,233 measured sections. The atomistic target is the 1252-atom carrier site,
not the 80,945-atom module, which is too large for a Canvas 2D redraw. The
atoms are drawn without bonds from the part scene, so the bonds are inferred
by distance. Simulated schematic.

## ADR-0134: The top-stage carry is not converged, and the top stage is off the critical path

Context. M29-03. A cascade stage pumps only while the free O2 a released pocket
carries back stays below one molecule per pass. ADR-0124 fit the carry from 16
single-seed runs as s_back 0.109, which stalls a single stage near 530 atm, and
called the result unconverged. M29-03 asks for longer, repeated same-pocket MD
at 600-1000 atm for every pocket of the O2 rotor.

Method. `scripts/top_stage_sweep.py` drives `scripts/pocket_carry_md.py` over
all eight O2-rotor pockets. Released pockets (the clearance case) at 600 and
1000 atm, and loaded pockets (a bound O2 on the iron) at 600 atm, 600 ps of
production each, 4 OpenCL workers, one seed. Sixteen released runs and eight
loaded runs.

Results. The released pocket carries n_c = 0.30 O2 on average at the top stage
(per pocket 0.04 to 0.87), nearly independent of pressure. A single stage then
pumps at 1000 atm, net 1 - 0.30 = 0.70 per pass. Fit against the far gas
density: s_back 0.060 +- 0.013 (16 runs), s_fwd 0.065 +- 0.022 (8 runs). Both
fits are through the origin and dominated by the pocket-to-pocket scatter.

But the number is not converged. The pocket load rises with production time:
one pocket reads 0.04 at 100 ps and 0.38 at 600 ps at the same density. The
pocket region (r < 0.9 nm, z 0 to 0.66 nm) is 1.7 nm3 and would hold several O2
at the far density if it were open and at equilibrium, so the 600 ps value is a
lower bound. The earlier ADR-0124 runs found 1.5 to 2.5 at the same density.

Finding. The sweep does not settle the top stage. It brackets the released
carry between about 0.3 (600 ps) and about 2 (the earlier runs); the top-stage
stall lies between about 500 and about 1000 atm. This agrees with ADR-0124: the
top stage depends on the favourable half of a wide scatter.

Decision. Do not rely on the 1000 atm top stage. ADR-0132 shows the real
circulation never charges the tank past a few atm, so the device is a
flow-through pump near ambient and the top stage is off the critical path. Keep
the cascade as design margin, not a requirement.

Limits. 600 ps production, one seed per pocket, no repeated runs. The
far-density slab (z 0.9-1.3 nm) is depleted by the rotor, so the fit x-axis is
a local density, not the bulk; the pocket region also contains rotor atoms, so
its free volume is below its geometric volume. Simulated estimates, not
converged.

Revision. ADR-0138 found a protocol error in the runs above (the bound O2 was
not removed from the released pockets) and gives the corrected carry. Read
ADR-0138 for the current number.

## ADR-0135: One header and an overview page for every view

Context. The site grew one page at a time: the atomistic viewer, the system
model, the 3D models page, the device explorer and the docs. Each page had its
own header and its own list of links. Three pages lacked a link that another
page had, the 3D models page used a different header style, and no page said
how the views relate. The user asked for a better interface, kept to plain
HTML, JS and Canvas.

Decision. One header for every page, from one source. `scripts/site_nav.py`
holds the nav list and writes the header between two comment markers in each
app page; `scripts/build_docs_site.py` imports the same function for the docs.
`site/nav.css` styles it with fixed colours and a fixed 46 px height, which the
canvas pages already subtract. The header is static HTML, not injected by a
script, so it works without JavaScript and on `file://`. A new entry page,
`site/overview.html`, explains each view, the level it covers and the data
flow between them, and lists where the study stands with links to the ADRs.
The docs builder gives each heading an anchor (`#adr-NNNN` for decisions), so
a page can link to one decision. `site/check.py` fails on a stale header and
on a broken local link or anchor.

Limits. `/` still serves the atomistic viewer, to keep existing links; the
brand mark leads to the overview. The status list on the overview page is
written by hand and must be updated with new findings. The system model page
still needs about 735 px of width because of its three fixed columns.


## ADR-0136: The site lifetime is about a day without repair; the reference rates were wrong

Context. ADR-0127 modelled the site lifetime from four autoxidation rates
with no sources. A check against the literature found two of them wrong. The
"myoglobin" rate, 1.2e-7 /s (67 days), is about 130 times too slow. The
"engineered pocket" rate, 1e-8 /s (800 days), matches no measured system; it
is a hope, not a reference. The carrier lifetime in ADR-0127 used the myoglobin
value.

Decision. `respirocyte_specs.py` now holds a sourced table
(`AUTOXIDATION_REFERENCES`, each with its PubMed id) and a mechanistic estimate
of our pocket (`pocket_autoxidation`). Measured references (37 C, pH 7, air):

| reference | half-life | k (1/s) | PMID |
|---|---|---|---|
| myoglobin H64A, no distal H-bond | 43 s | 1.6e-2 | 8463233 |
| albumin-heme, covalent, red-cell p50 | 5 h | 3.9e-5 | 15656571 |
| myoglobin, native | 12.6 h | 1.5e-5 | 8463233 |
| hemoglobin A | 24 h | 8.1e-6 | 3422420 |
| albumin-heme FecycP(His), best synthetic | 25 h | 7.7e-6 | 12009926 |
| heme-cyclodextrin (25 C) | 30 h | 6.4e-6 | 19178464 |
| myoglobin L29F, crowded pocket | 5.8 d | 1.4e-6 | 8463233 |

The model is the two-channel mechanism of Brantley et al. (1993, PMID 8463233):
superoxide (HO2) dissociation from the oxy site, k_ox Y(p), and the
bimolecular reaction of O2 with a deoxy site that holds a water,
k_bi p (1 - Y(p)). It is fitted to hemoglobin in plasma at three pO2 (Furlan
and Bachofen 1976, PMID 1248570): worst residual 6 %. The fitted oxy channel
is 0.70 times the independent HbA-in-air rate, a fair cross-check. Two pocket
factors scale the channels: f_ox (the distal H-bond and pocket polarity;
assumed 0.3 / 1 / 30) and f_bi (water access; estimated 0.01 / 0.1 / 1 from
the dry pocket of ADR-0130 and the crowded L29F pocket). The rate is
time-weighted over arterial (95 mmHg, 20 %) and venous (40 mmHg, 80 %) blood.

Results (dendrimer-heme, p50 26.8 mmHg, 73 g in the body):

| case | f_ox | f_bi | half-life | replacement |
|---|---|---|---|---|
| hemoglobin-like pocket | 1 | 1 | 9.5 h | 60 g/day |
| optimistic | 0.3 | 0.01 | 5.5 d | 8.6 g/day |
| central | 1 | 0.1 | 31 h | 30 g/day |
| pessimistic | 30 | 1 | 1.3 h | 73 g/day |

Findings. (1) No measured heme carrier reaches the 8-day target (k below
1e-6 /s); the best is myoglobin L29F at 5.8 days, the best synthetic about
1 day. (2) Our pocket's central estimate is about 31 h, so a carrier without
repair needs about 30 g of new carrier each day. (3) Once the pocket is dry,
the superoxide channel dominates (67 % central), so the lever is f_ox, the
distal H-bond, not the affinity: moving the p50 from 10 to 40 mmHg changes the
half-life only from 33 to 29 h. (4) The guard is a hydroxyl, a proton donor;
Brantley found that a more polar pocket speeds protonation and autoxidation.
The guard may help or hurt, and only a QM or experimental test will decide it.
This is the M28-22 risk, now with a number.

Limits. The two pocket factors f_ox and f_bi are assumptions, not measurements.
The central case is a judgement, not a fit. The two-channel fit uses three
points, so it constrains two parameters with one degree of freedom. The rate
table mixes temperatures (37 C and 25 C) and solvent conditions; the myoglobin
H64A value is a lower bound on the distal effect, not our pocket. The model
gives the site lifetime, not the device lifetime; the ring damage from
superoxide (M28-22) and the repair path are still open.


## ADR-0137: The tether is a rigid ethynyl link; the old mount left five-bonded carbons

Context. The carrier site joined each porphine meso carbon to a diamond wall
carbon with a bridge of one to three sp3 carbons. A first strain check
(`scripts/tether_strain.py`, a harmonic bond-angle field) found that the
as-built three-carbon bridge carries about 7000 kJ/mol of strain. Worse, a
contact check found 90 non-bonded site-block contacts below 2.0 A, with a
hydrogen to hydrogen distance of 0.59 A. The cause was a defect, not tight
geometry: the builder bonded each tether to a wall carbon and did not remove
the wall hydrogen, so that carbon had five bonds.

Decision. Two changes. (1) `Topology::without_atoms` drops the replaced
hydrogens and remaps the site indices. The builder removes the anchor
hydrogen, so every carbon keeps four bonds, and it picks the wall carbon and
bridge length with the lowest strain. (2) The bridge is an ethynyl link,
meso-C#C-C(wall), in the CAD. A rigid, linear two-carbon link fixes the ring
to the wall and holds the porphine plane, and its sp carbons carry no
hydrogens, so it cannot form the close H...H contacts that the sp3 bridge
formed. The builder places the two sp carbons and removes the wall hydrogen.
`scripts/mount_uff.py` checks the result with an independent force field
(RDKit UFF): it relaxes the site in a frozen scaffold shell, reports the mount
strain and the close contacts, and runs a quasi-static hold test that moves
the ring core rigidly and fits the spring constant.

Results. The corrected sp3 build (slot radius 7.0 A, half-width 2.0 A) has
1177 scaffold atoms, four C-C anchor bonds at 1.55-1.60 A, and 36 site-block
contacts below 2.0 A (from 90). The UFF mount strain is 866 kJ/mol, and the
tethered ring is stiff: axial rms 0.03 A, in-plane 0.03 A, tilt 0.56 deg at
310 K. The control with the four tether bonds cut is soft (axial rms 0.09 A,
tilt 1.5 deg), so the tethers, not the slot, hold the ring. The ethynyl
screen (harmonic field) gives the lowest strain at slot radius 7.0-7.5 A,
half-width 3.0 A: the two independent links carry about 180 and 550 kJ/mol.

Limits. The ethynyl link is a design, not a synthesis. UFF is a crude force
field for a metal-porphyrin; it gives the mount stiffness and the contact
count, not the Fe-N or Fe-O energetics. The hold test freezes the scaffold
shell at 7.0 A, so a wall mode below that radius is not seen. The reaction
that forms the C#C bond, and the repair of a broken tether, are not modelled.

## ADR-0138: The corrected top-stage carry; a protocol error in ADR-0134

Context. ADR-0134 measured the top-stage carry from released pockets (0.30 O2)
and concluded that a single stage pumps at 1000 atm with net 0.70. A review
found a protocol error. The module scene keeps a bound O2 on the iron in the
five loaded pockets. The ADR-0134 released runs did not remove it, so most of
them measured a pocket that still held its bound O2 and had less room for free
gas; the loaded runs stacked a second O2 on top. ADR-0134 also used a short gas
box (z_front 1.4 nm) whose far slab the rotor surface depleted.

Method. `top_stage_sweep.py` now always passes `--remove-bound-o2`, and the
loaded runs use `--add-bound-o2` to place the bound O2 at the scanned geometry.
The gas box is 4.4 nm tall, and the start density is calibrated (13.3 at
600 atm, 17.6 at 1000 atm) so the near slab reaches the equation-of-state
density. `pocket_carry_md.py` gains `--prefill` and `--start-empty`, so each
released run starts empty or prefilled and brackets equilibrium from both
sides. 40 runs, 400 ps each, one seed: 16 released at 600 atm, 16 at 1000 atm,
8 loaded at 600 atm. The fit uses the last 200 ps.

Results (free O2 in the pocket, mean of the last 200 ps):

| condition | empty start | prefilled start |
|---|---|---|
| released, 600 atm | 0.40 | 0.64 / 1.43 |
| released, 1000 atm | 1.18 | 1.10 / 1.33 |
| loaded, 600 atm | 0.63 (fit) | - |

Calibration: the near slab holds 8.7 O2/nm3 at 600 atm and 12.2 at 1000 atm
(the equation of state gives 10.4 and 14.1). Fits of n = s * rho through the
origin: released s_back 0.076 +- 0.009 (32 runs), loaded s_fwd 0.064 +- 0.014
(8 runs).

Findings. (1) The corrected released carry is about 1.06 O2 at 1000 atm, not
0.30: ADR-0134 measured a pocket that still held its bound O2. (2) A single
stage from blood stalls at about 920 atm, not 530 atm. (3) The top stage
(600 to 1000 atm) still pumps, net 1 + 0.63 - 1.06 = 0.57 per pass, less than
the 0.70 of ADR-0134. (4) The margin is thin and pocket-dependent: seven
pockets pump (+0.25 to +0.98) and one, pocket 6, has a negative net (-0.57).
(5) At 1000 atm the two starts agree (1.18 empty, 1.10 prefilled), so that
number is converged; at 600 atm the empty start is a lower bound, because the
pocket fills slowly at the lower density.

Decision. Withdraw the ADR-0134 carry. The top stage pumps on average, but it
is marginal. Keep the cascade as design margin, not a requirement, as ADR-0134
decided, and do not count on every pocket. The device remains a flow-through
pump near ambient (ADR-0132).

Limits. One seed per pocket; the per-pocket scatter is wide. The pocket region
contains rotor atoms, so its free volume is below its geometric volume. The
near-slab density is not the bulk density. Simulated estimates.


## ADR-0139: The site can be self-sustaining; protect the ring, do not repair it

Context. M29-14. Can the heme site be repaired without new material each day,
and can a fully self-sustaining site last a very long time? ADR-0108 handled
the reversible iron oxidation with a donor station. ADR-0136 set the site
lifetime (central 31 h) and named the irreversible ring damage as the open
risk (M28-22).

Decision. The site degrades two ways. The iron oxidation is reversible and
needs only electrons. The ring damage is irreversible: a broken ring needs new
carbon, so no self-sustaining device can repair it. The strategy is to stop
the ring from breaking, with three permanent elements:

  1. the iron-repair donor on the empty arc (ADR-0108), so the site stays in
     the Fe(II) state;
  2. a distal pocket that suppresses autoxidation about 100-fold (a distal
     H-bond holds the bound O2 and slows the superoxide channel; a dry pocket
     keeps water off the iron and slows the water channel);
  3. a co-located catalyst that dismutates the superoxide before it reaches
     the ring, and then decomposes the peroxide.

None of the three is consumed. The inputs are electrons (a charged onboard
store) and protons (buffered blood). A consumable antioxidant would need new
material each day and is rejected.

Model. `respirocyte_repair.py` estimates the ring lifetime from the race
between the catalyst, escape, and the ring attack. The superoxide is born on
the iron; its irreversible fraction is phi = k_ring / (k_dismut + k_escape +
k_ring), and the ring lifetime is ln2 / (k_autox phi). The rates (k_escape
1e4 /s, k_ring 1e5 /s, catalyst 0 to 1e8 /s) are estimates. The rotor limit is
1e11 cycles at 1000 rev/s, 3.2 years.

Results (ring lifetime, years):

| pocket | no catalyst | Mn mimic | enzyme-like |
|---|---|---|---|
| central (1x) | 0.004 | 0.004 | 0.04 |
| H-bonded, dry (10x) | 0.04 | 0.04 | 0.4 |
| optimized (100x) | 0.4 | 0.4 | 3.9 |
| native-mutant span (1000x) | 3.9 | 4.2 | 39 |

Findings. (1) The iron repair has a 16,000-fold margin, so the iron is always
reduced (oxidized fraction 6e-5). (2) The optimized pocket alone gives 0.4
years and the catalyst alone 0.04 years; together they give 3.9 years, so the
site then outlives the 3.2-year rotor. (3) The electron budget is negligible:
5e-7 electrons per second for a dendrimer carrier, about 8e-26 A. (4) To reach
decades, use a stronger pocket (the 1000-fold native-to-mutant span of
Brantley) or a faster catalyst; the 1000-fold pocket with the enzyme-like
catalyst gives 39 years.

Decision. Build the site with the optimized distal pocket, the co-located
dismutase and catalase catalyst, and the iron-repair donor. Then the rotor, or
the scaffold, sets the device lifetime, not the site chemistry. For a passive
carrier, with no rotor, the site lifetime is the limit.

Sourced envelope (M29-15). The assumed rates are replaced by ranges from the
literature, in `FATE_REFERENCES`. The spontaneous dismutation is 0.1 to 1 /s
(Bielski 1985, at about 1 uM superoxide). A Mn-porphyrin SOD mimic turns over
1e6 to 1e8 /s. An enzyme-like catalyst reaches 1e5 to 1e9 /s. Escape by
diffusion is 1e6 to 1e10 /s. The ring attack is 1e2 to 1e6 /s; it is the least
certain value. `fate_scenarios` combines the worst and best ends. For the
optimized pocket (100x), the ring lifetime is 1.1 years pessimistic, 3.9 years
central, and 4e7 years optimistic. So the design outlives the 3.2-year rotor
unless the catalyst is slow (1e6 /s) and the ring attack is fast (1e6 /s) at
the same time. The catalyst turnover is the lever.

Limits. The superoxide fate rates are ranges, not measurements. The ring-attack
rate is the least certain, so the pessimistic case is a lower bound, not a
prediction. The 100-fold pocket suppression is a design goal, not a
demonstration. The catalyst is a design, not a synthesis. The M28-22 ring-damage
risk stays open until the ring-attack rate is measured.

## ADR-0140: The ring damage is mediated by H2O2; the rate is measured

Context. M28-22. ADR-0139 named the ring-attack rate as the largest uncertainty
and used an assumed range (1e2 to 1e6 /s). Can that rate come from measured
data instead of a guess?

Decision. The superoxide does not attack the ring directly. It converts almost
completely into hydrogen peroxide by dismutation (Wazawa 1992, PMID 1420896).
The H2O2 is the oxidant: it reacts with the Fe(II) heme to form a ferryl, which
opens the ring (Yusa & Shikama 1987, PMID 3427037). So the ring damage is
bounded by the heme-peroxide reaction, and that reaction is measured.

Model. `respirocyte_repair.py` now estimates the ring lifetime from the race
for the H2O2 between the heme and the co-located catalyst. The ring-damage
fraction is phi = k_heme / (k_heme + k_cat), where k_heme is the measured
deoxy-heme + H2O2 rate (3.6e3 /M/s) and k_cat is the catalyst's catalase-like
rate. The ring lifetime is ln2 / (k_autox phi). The no-catalyst case (phi = 1)
is an upper bound on the damage, because the ferryl can return to the resting
state and the superoxide can escape, both of which lower it.

Results (ring lifetime, years):

| pocket | no catalyst | Mn mimic (1e5) | catalase (1e7) |
|---|---|---|---|
| central (1x) | 0.004 | 0.10 | 9.8 |
| H-bonded, dry (10x) | 0.035 | 1.0 | 98 |
| optimized (100x) | 0.35 | 10.2 | 981 |
| native-mutant span (1000x) | 3.5 | 102 | 9800 |

Findings. (1) The measured heme-peroxide rate (3.6e3 /M/s) replaces the assumed
ring-attack range. (2) The optimized pocket with a Mn-porphyrin mimic gives a
10-year ring lifetime, so the site outlives the 3.2-year rotor with a catalyst
no faster than a known mimic. With a catalase it gives about 1000 years. (3)
The catalyst must decompose the H2O2, not only dismutate the superoxide: a pure
dismutase produces the peroxide that damages the ring. (4) The design needs a
catalase-like rate of about 1e5 /M/s or more.

Limits. The heme + H2O2 rate is measured for sperm-whale myoglobin, not for
this pocket, and it is an upper bound: not every ferryl opens the ring. The
superoxide-to-H2O2 conversion is taken as complete (conservative; escape only
lowers the damage). The catalase-like rate of a Mn-porphyrin is an estimate.
The 100-fold pocket suppression is a design goal. The catalyst is a design, not
a synthesis.
