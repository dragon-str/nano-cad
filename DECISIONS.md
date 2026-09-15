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
