# nano-cad planetary gearbox — written report

This report collects the measured numbers of the nano-cad project. It is the
M8-03 deliverable. Every number in this report cites the file or the test that
produced it. If no source gives a number, the report says that the number is
unknown.

The report uses SI units. It states the unit for every number.

> **Status note.** All physics in this report is **simulated**. No physical
> device exists. No result is validated for medicine. Several numbers are
> estimates and are marked as estimates.

## 1. What was implemented and simulated

nano-cad is a Rust workspace with nine crates (`units`, `model`, `format`,
`engine`, `jigs`, `parts`, `params`, `python`, `cli`). Source: `TASKS.md`
M0-01.

The workspace holds about 425 Rust tests and 40 Python and MCP tests. Source:
`README.md`. The M0-01 result note records 68 tests at the skeleton stage.
Source: `TASKS.md` M0-01.

The implemented parts are:

- The units crate with 7 SI base units and 7 non-SI units, with 17 passing
  tests. Source: `TASKS.md` M0-06.
- The document model with atoms, bonds, and topology. The model test count is
  43. Source: `TASKS.md` M1-01. The model byte codec carries the `NCAD` magic
  and `SCHEMA_VERSION = 1`. Source: `TASKS.md` M1-02.
- The NCZ, XYZ, and MMP formats. NCZ is a stored ZIP with `schema_version = 1`
  and f64 coordinates. Source: `TASKS.md` M1-03. MMP units are `1e-13 m`, which
  is 0.001 Angstrom. Source: `TASKS.md` M1-05.
- The CLI `nanocad` with `info` and `convert`, and 13 passing tests. Source:
  `TASKS.md` M1-08.
- The L1 force engine and the L2 device layer, described in sections 3 and 4.
- Six part generators: diamond, graphite, nanotube, gear profile, spur gear,
  and planetary set. Source: `TASKS.md` M4-02 to M4-06.
- The parameter record, the extraction pipeline, and the parameter library.
  Source: `TASKS.md` M5-01 to M5-07.
- The Python bindings. The `nanocad._core` extension passes 27 Python tests.
  Source: `TASKS.md` M7-01.
- The MCP server with 18 tools. Source: `TASKS.md` M7-02 and `mcp/tools.json`.
- The three-scale scene export. Source: `TASKS.md` M8-01 and
  `docs/three-scale.md`.

A golden round-trip test imported a real MMP file, saved it to NCZ, reloaded
it, and reproduced the coordinates bit for bit. The file had 6828 atoms and
10452 bonds. Source: `TASKS.md` M1-07.

The Apache-2.0 license text is 202 lines. Source: `TASKS.md` M0-04.

## 2. The multiscale design (L0 to L4)

The architecture defines five levels. Source: `ARCHITECTURE.md` and
`docs/glossary.md`.

| Level | Name | Role | State |
|---|---|---|---|
| L0 | Quantum adapter | PySCF or xTB through ASE | **Not implemented** |
| L1 | Atomistic engine | Rust engine for parts, forces, and jigs | Implemented |
| L2 | Device layer | Rigid bodies, joints, and the parameter handoff | Implemented |
| L3 | Continuum adapter | External FEA and flow solvers | **Not implemented** |
| L4 | System layer | Lumped ODE models with CellML or SBML interop | **Not implemented** |

Only L1 and L2 are implemented. L0, L3, and L4 are later work. The open tasks
are `L0-01`, `L3-01`, and `L4-01`. Source: `TASKS.md` (Phase 6).

The bridge between L1 and L2 is the parameter schema in `PARAMETERS.md`. Every
value carries a unit, a provenance record, and an uncertainty. An unknown
uncertainty is `null`. Source: `PARAMETERS.md`.

The L1 to L2 handoff is the extraction pipeline. It relaxes a part, measures
stiffness, friction, failure stress, thermal properties, and writes a
`PartRecord`. Source: `PARAMETERS.md` and `TASKS.md` M5-01 to M5-05.

The device layer reads the `PartRecord` and assembles rigid bodies and joints.
Source: `TASKS.md` M6-01 to M6-05.

The documentation number for this handoff is `nanocad.param.part`, version 1.
Source: `PARAMETERS.md`.

## 3. The force field and the finite-difference gradient errors

The L1 engine has six force terms. Every term has a central finite-difference
gradient test. Source: `TASKS.md` M2-01 to M2-08 and `ARCHITECTURE.md`.

The test for a term is named
`the_analytic_gradient_matches_a_central_finite_difference`. Source: for
example `crates/engine/src/bond_stretch.rs` line 424 and
`crates/engine/src/system.rs` line 712.

| Force term | Maximum finite-difference error | Source |
|---|---|---|
| Bond stretch | 1.28e-17 N | `TASKS.md` M2-02 |
| Angle bend | 8.3e-18 N | `TASKS.md` M2-03 |
| Torsion (three-term cosine series) | 7.5e-18 N | `TASKS.md` M2-04 |
| Out-of-plane (harmonic in the signed angle) | 9.9e-17 N | `TASKS.md` M2-05 |
| van der Waals (Buckingham with a switching cutoff) | 7.0e-20 N | `TASKS.md` M2-06 |
| Electrostatic (Coulomb with a switch) | 1.7e-18 N | `TASKS.md` M2-07 |
| Total potential on a 6-atom chain with six terms | 5.7e-17 N | `TASKS.md` M2-08 |

The bond-stretch test uses a central finite difference with a step of
`1e-14 m`. Source: `TASKS.md` M2-02. The unit of every error above is the
newton (N).

Two jig force terms also have finite-difference tests:

| Jig term | Maximum finite-difference error | Source |
|---|---|---|
| Anchor jig | 4.8e-22 N | `TASKS.md` M3-04 |
| Spring jig | 1.0e-16 N | `TASKS.md` M3-06 |

The non-bonded terms read the Verlet pair list after M2-08. Source: `TASKS.md`
M2-08. The collective test is
`the_nonbonded_pair_list_path_matches_the_all_pairs_reference` in
`crates/engine/src/system.rs` line 750.

The conjugate-gradient minimizer converged to a gradient of `5.9e-15 N` in 53
iterations. Source: `TASKS.md` M2-09. The minimizer test is
`a_diatomic_reaches_its_equilibrium_length` in
`crates/engine/src/minimize.rs` line 371.

### Caveat on the force field

- The non-bonded terms use a cutoff. Long-range electrostatics beyond a cutoff
  is a non-goal for the first release. Source: `PLAN.md`.
- The reference water-box cross-check against OpenMM is **done**. Source:
  `benchmarks/results/openmm-water-crosscheck.txt`. See the subsection below.
- The MCP `relax` and `run_md` tools assign bond-stretch terms only. They do
  not assign angles, torsions, van der Waals, or electrostatic terms. Source:
  `mcp/README.md`.

### The OpenMM cross-check

Task M2-11 compares the periodic SPC-like water box against OpenMM 8.6.1 on
the `Reference` (double precision) platform. Both codes use the same box, the
same parameters, and the same functional forms. The comparison uses 200
molecules (600 atoms), a cubic box of 1.8833259369189425e-9 m, and the fixed
seed `0x5A17E2`. Source:
`benchmarks/results/openmm-water-crosscheck.txt`.

| Quantity | nanocad | OpenMM | Relative difference |
|---|---|---|---|
| Potential energy | -2.737699062416127e-17 J | -2.7376990624161912e-17 J | 2.341e-14 |
| Maximum per-atom force difference | — | max magnitude 1.4621685148513376e-08 N | 1.242e-14 |

The declared tolerance is `1.0e-6` relative. The result is **PASS**. The
`CPU` platform agrees to 4.32e-6 on the energy and 3.99e-6 on the force, which
is single-precision evaluation, not a model mismatch. Source:
`benchmarks/results/openmm-water-crosscheck.txt`.

**This is a coding check, not a physics validation.** Both codes use the same
cutoff and no long-range correction. The model is flexible and SPC-like, so it
is not comparable to a standard SPC or TIP3P result. One fixed configuration
is checked. A geometry-symmetric model error could pass.

## 4. The integrator energy drift

The integrator is Velocity Verlet. Source: `TASKS.md` M2-10 and M3-01. The
integrator test is `the_nve_total_energy_drift_stays_below_the_bound` in
`crates/engine/src/integrator.rs` line 216.

| Quantity | Value | Bound | Source |
|---|---|---|---|
| NVE relative energy drift, 4000 steps | 1.278e-5 | 1.0e-4 | `TASKS.md` M2-10 |
| Periodic water-box NVE relative drift | 1.274e-6 | 1.0e-3 | `TASKS.md` M2-11 |

The drift is dimensionless. The second bound is `1.0e-3` in the test source
`crates/engine/src/reference.rs` line 339.

Thermostat results:

| Thermostat | Result | Source |
|---|---|---|
| Langevin | mean 300.58 K at a 300 K target | `TASKS.md` M3-02 |
| Berendsen | relaxed from 10.23 K to 300.00 K | `TASKS.md` M3-03 |

The Berendsen thermostat is a relaxation thermostat. It is not a correct
canonical sampler. Source: `TASKS.md` M3-03 and `DECISIONS.md` ADR-0018.

The rotor test ran 200000 frames at a 300 K target. The mean temperature was
309.57 K. The maximum bond strain was 1.965e-2 (dimensionless). The maximum
force residual was 8.564e-26 N. Source: `TASKS.md` M3-07. The rotor test is
`a_motor_drives_the_rotor_and_it_stays_stable_at_300_k` in
`crates/jigs/src/rotor.rs`.

The parallel force reduction is deterministic. One-thread and N-thread results
are bit-identical. Source: `TASKS.md` M3-08. The test is
`the_deterministic_parallel_forces_match_for_one_and_many_threads`. Source:
`crates/engine/src/system.rs`.

## 5. The gearbox gear ratio result

The planetary set obeys `ring_teeth = sun_teeth + 2 * planet_teeth`. Source:
`TASKS.md` M4-06. The default design has 24 sun teeth, 18 planet teeth, 60
ring teeth, and 3 planets. Sources: `TASKS.md` M4-06 and `docs/three-scale.md`.

The default gear ratio is 3.5. Sources: `TASKS.md` M4-06 and
`docs/three-scale.md`.

The device test measured the sun-to-carrier ratio through the M6 gear-ratio
test and through the MCP tool demo:

| Quantity | Value | Source |
|---|---|---|
| Analytic ratio | 3.500000000 | `TASKS.md` M6-06, `docs/reproduce.md` |
| Measured ratio | 3.500000000 | `TASKS.md` M6-06, `docs/reproduce.md` |
| Relative error | 1.146e-10 | `TASKS.md` M6-06, `docs/reproduce.md` |
| Tolerance | 1.0e-3 | `TASKS.md` M6-06, `mcp/tests/test_agent_demo.py` |

The ratio is dimensionless. The M6 test is
`the_gear_ratio_matches_the_analytic_ratio_within_tolerance` in
`crates/jigs/src/assembly.rs` line 662. The MCP demo test is
`test_agent_builds_the_planetary_gearbox_through_tools` in
`mcp/tests/test_agent_demo.py`.

The MCP demo also reports a maximum gear-constraint error below `1.0e-9 m`.
Source: `mcp/tests/test_agent_demo.py`.

The URDF export has 7 links and 6 joints. Source: `TASKS.md` M7-04 and
`docs/reproduce.md`. The test asserts `link_count == 7` and `joint_count == 6`.
Source: `mcp/tests/test_agent_demo.py`.

The scene export has 7 bodies, 6 joints, and 3 planets. Source: `TASKS.md`
M8-01. The scene design states the module as `5e-10 m`, the sun pitch radius
as `6e-9 m`, the planet pitch radius as `4.5e-9 m`, the ring pitch radius as
`1.5e-8 m`, and the carrier radius as `1.05e-8 m`. Source:
`docs/three-scale.md`.

## 6. The device constraint errors

The L2 device layer uses semi-implicit Euler with a position-level
Gauss-Seidel projection. There is no velocity-level correction pass.
Source: `DECISIONS.md` ADR-0014.

| Joint or constraint | Measured error | Steps | Source |
|---|---|---|---|
| Revolute anchor | 0 m | 2000 | `TASKS.md` M6-02 |
| Revolute axis | 9.9956e-13 rad | 2000 | `TASKS.md` M6-02 |
| Prismatic position | 0 m | 2000 | `TASKS.md` M6-03 |
| Prismatic orientation | 5.42e-20 rad | 2000 | `TASKS.md` M6-03 |
| Gear coupling | 5.551e-17 m | many | `TASKS.md` M6-04 |

The gear coupling gave a ratio of `-0.500000000` for `r_a = 0.05 m` and
`r_b = 0.10 m`. Source: `TASKS.md` M6-04. The ratio is dimensionless.

The relevant tests are:

- `a_revolute_joint_starts_satisfied_and_rejects_bad_setup` in
  `crates/jigs/src/device.rs` line 1913.
- `a_prismatic_joint_constrains_the_perpendicular_motion` in
  `crates/jigs/src/device.rs` line 1996.
- `a_gear_coupling_holds_the_ratio_over_many_steps` in
  `crates/jigs/src/device.rs` line 2042.

### Caveat on the constraints

- The constraint pass is position-level only. A residual constraint velocity
  can persist. A velocity-level pass stays open. Source: `TASKS.md` M6-04 note
  and `DECISIONS.md` ADR-0014.
- The readings above are from the tested cases only. They are not a general
  accuracy bound.

## 7. The parameter extraction results

The L1 to L2 extraction pipeline produced these results. All are **simulated**
and all are **estimates**, not material constants.

| Property | Extracted value | Analytic value | Relative error | Source |
|---|---|---|---|---|
| Elastic modulus | 4.49999999999e12 Pa | 4.5e12 Pa | 2.32e-12 | `TASKS.md` M5-02 |
| Failure stress | 2.25e11 Pa | 2.25e11 Pa | not stated | `TASKS.md` M5-04 |
| Specific heat | 1.0414e3 J/(kg·K) | 1.0393e3 J/(kg·K) | 2.0e-3 | `TASKS.md` M5-05 |
| Thermal conductivity | 8.34 +/- 0.93 W/(m·K) | unknown | unknown | `TASKS.md` M5-05 |

The stiffness test used a synthetic chain with a stiffness of `300 N/m`. The
stiffness test is `a_harmonic_chain_recovers_its_known_modulus` in
`crates/params/src/stiffness.rs` line 203. Source of the `300 N/m` value:
`TASKS.md` M5-02.

The failure test used a rupture strain of `0.05`. The failure test is
`a_known_rupture_strain_gives_the_known_failure_stress` in
`crates/params/src/failure.rs` line 305. The rupture strain is also a constant
in `crates/params/src/failure.rs` line 281.

The friction extraction gave a sliding-contact coefficient of `0.4433` with a
spread of `3.37e-12 N`. The coefficient is dimensionless. The friction test is
`the_sliding_contact_repeats_and_reports_a_spread` in
`crates/params/src/friction.rs` line 359. Source: `TASKS.md` M5-03. The
friction value is a simulated estimate.

The specific-heat test is
`a_monatomic_ideal_gas_recovers_the_equipartition_specific_heat` in
`crates/params/src/thermal.rs` line 804. The conductivity test is
`a_chain_gives_a_two_bath_conductivity_that_repeats_and_is_positive` in
`crates/params/src/thermal.rs` line 842.

### Caveats on the extraction

- Both thermal properties are marked `Method::Md` and `Validation::Unverified`.
  Source: `TASKS.md` M5-05.
- The specific heat is a kinetic estimate. It omits `Var(U)`. The conductivity
  is an effective value for the stated geometry. A harmonic chain shows
  anomalous transport. Sources: `DECISIONS.md` ADR-0018 and `TASKS.md` M5-05.
- The conductivity has **no analytic reference value**. Its accuracy is
  unknown.
- The MCP server does not expose `extract_parameters`. Source: `mcp/README.md`
  and `TASKS.md` M7-03.
- The `PARAMETERS.md` example shows an elastic modulus of `1.05e12 Pa` and a
  friction coefficient of `0.05`. Source: `PARAMETERS.md`. Those are schema
  examples, not measured results. They must not be quoted as results.

## 8. The benchmark timings

The benchmark hardware is an Apple M4 with 10 logical cores, macOS 26.6.2
(build 25G83). Source: `benchmarks/results/2026-09-14-engine-timings.txt`. The
toolchain is cargo 1.98.1 and rustc 1.98.1. Source:
`benchmarks/results/20260915T010005Z-engine-timings.txt`.

The reference case is a periodic SPC-like water box with 200 molecules and 600
atoms. The box side is `1.884e-9 m`. The temperature is `300 K`. The seed is
`0xB3C0`. The neighbor pair list holds 30827 pairs. Source:
`benchmarks/results/2026-09-14-engine-timings.txt` and
`crates/engine/benches/engine_bench.rs`.

The first recorded run (`2026-09-14-engine-timings.txt`):

| Benchmark | Mean | Lower bound | Upper bound |
|---|---|---|---|
| `water_200_energy_and_gradient` | 680.97 µs | 677.92 µs | 685.88 µs |
| `water_200_forces_serial` | 679.83 µs | 678.45 µs | 681.14 µs |
| `water_200_forces_parallel_4` | 370.97 µs | 369.07 µs | 373.05 µs |

The 4-worker deterministic reduction is about 1.83 times faster than the serial
path on this host. Source:
`benchmarks/results/2026-09-14-engine-timings.txt`. The unit is the
microsecond (µs).

The second recorded run on the same host
(`20260915T010005Z-engine-timings.txt`):

| Benchmark | Criterion time |
|---|---|
| `water_200_energy_and_gradient` | [648.79 µs 653.02 µs 659.23 µs] |
| `water_200_forces_serial` | [656.13 µs 660.21 µs 665.04 µs] |
| `water_200_forces_parallel_4` | [284.23 µs 303.08 µs 317.76 µs] |

Source: `benchmarks/results/20260915T010005Z-engine-timings.txt`.

The parallel mean moved from 370.97 µs to 303.08 µs. That is a decrease of
about 18 percent, computed from the two cited means. Source:
`docs/benchmarks.md`.

The exact command is:

```
cargo bench -p nanocad-engine --bench engine_bench -- --warm-up-time 0.5 --measurement-time 1.0 --sample-size 20
```

Source: `benchmarks/results/2026-09-14-engine-timings.txt`.

### Caveats on the timings

- Every timing number is machine-specific. The numbers are not universal.
  Source: `docs/benchmarks.md`.
- `cargo bench` is timing only. It is not a merge gate. Source:
  `benchmarks/results/2026-09-14-engine-timings.txt`.
- The parallel path is opt-in. The default engine path is serial. Source:
  `benchmarks/results/2026-09-14-engine-timings.txt` and `DECISIONS.md`
  ADR-0017.
- The parallel number is the least stable. Thread scheduling changes it.
  Source: `docs/benchmarks.md`.
- Criterion reports a 95 percent confidence interval. Source:
  `docs/benchmarks.md`.
- No GPU benchmark exists. The GPU path is a later plan. Source
  `ARCHITECTURE.md`.

## 9. Reproduction

A clean checkout reproduces the demo with one command:

```
just reproduce
```

Source: `README.md` and `docs/reproduce.md`. The recipe calls
`scripts/reproduce.sh`. Source: `docs/reproduce.md`.

The script runs 7 steps: toolchain check, virtual environment, install of
maturin and pytest, extension build, `cargo test --workspace`, the MCP demo
test, and the key-result print. Source: `scripts/reproduce.sh` and
`docs/reproduce.md`.

The observed key result line is:

```
REPRODUCE OK: gearbox demo: analytic_ratio 3.500000000 measured_ratio 3.500000000 relative_error 1.146e-10 tolerance 1.0e-03
```

Source: `docs/reproduce.md`. The reference host for that line is an Apple M4
with Python 3.14 and cargo 1.98.1. Source: `docs/reproduce.md`.

The demo exports a URDF with 7 links and 6 joints. Source: `docs/reproduce.md`.

Prerequisites: Rust 1.85 or newer (MSRV), Python 3.10 or newer, a POSIX
`/bin/sh`, and network access on the first run. Source: `docs/reproduce.md`.
The disk estimate is about 1 GB for `target/` and the virtual environment.
Source: `docs/reproduce.md`. That disk number is an estimate.

The first run needs network access to download `maturin` and `pytest`. Source:
`docs/reproduce.md`.

### Caveats on reproduction

- The script does not run `cargo bench`. The benchmark numbers are not
  reproduced. Source: `docs/reproduce.md`.
- The script does not run the OpenMM cross-check. Source: `docs/reproduce.md`.
- The `maturin` build must run from `python/`. Source: `docs/reproduce.md` and
  `DECISIONS.md` ADR-0021.
- The tree must compile when the script runs. Source: `docs/reproduce.md`.

## 10. Caveats and limitations

### Scope limits

- Only L1 and L2 are implemented. L0, L3, and L4 are not implemented. Sources:
  `ARCHITECTURE.md`, `TASKS.md` (Phase 6).
- The gearbox is a design hypothesis. No physical device exists. Source:
  `PLAN.md` ("Any claim of a built or medically validated device" is a
  non-goal).
- No medical claim is made. Source: `README.md` doctrine.

### Physics limits

- The OpenMM water-box cross-check is done and passes at 1e-14 relative. It is
  a coding check, not a physics validation. Source:
  `benchmarks/results/openmm-water-crosscheck.txt`.
- No long-range electrostatics beyond a cutoff. PME is a later plan. Source:
  `PLAN.md` and `ARCHITECTURE.md`.
- The default engine path is serial. The parallel path is opt-in. Source:
  `TASKS.md` M3-08 and `DECISIONS.md` ADR-0017.
- No SIMD path is implemented. Source: `ARCHITECTURE.md`.
- The `run_md` tool uses one carbon mass for every atom. That is an estimate,
  not a per-element mass. Source: `mcp/README.md`.

### Device limits

- The constraint pass is position-level only. A velocity-level pass stays open.
  Source: `TASKS.md` M6-04 and `DECISIONS.md` ADR-0014.
- The MCP server keeps state in memory. A restart loses every handle. Source:
  `mcp/README.md` and `DECISIONS.md` ADR-0019.

### Parameter limits

- The extraction results are simulated estimates. The stiffness and failure
  tests use synthetic systems. Source: `crates/params/tests/extraction.rs`
  header and `TASKS.md` M5-02 to M5-05.
- The thermal properties are `Validation::Unverified`. Source: `TASKS.md`
  M5-05.
- The thermal conductivity has no analytic reference. Its accuracy is unknown.
- The `PARAMETERS.md` example values are schema examples, not measured
  results. Source: `PARAMETERS.md`.

### Visualization limits

- The atomistic scene layer is a static snapshot. It has no bonds and no atom
  type. Source: `docs/three-scale.md` and `DECISIONS.md` ADR-0020.
- The coarse layer uses the device axis only. A tilted body does not tilt its
  bounding cylinder. Source: `docs/three-scale.md`.
- A prismatic joint reports a null world anchor. Source: `docs/three-scale.md`.

### Benchmark limits

- The timings are from one host and one toolchain. They are not universal.
  Source: `docs/benchmarks.md`.
- `cargo bench` is timing only. It is not a merge gate. Source:
  `docs/benchmarks.md`.

### Open tasks

| Task | State | Source |
|---|---|---|
| M2-11 OpenMM water-box cross-check | Done, passes | `benchmarks/results/openmm-water-crosscheck.txt` |
| M8-02 90-second video | Script written, video not recorded | `docs/video-script.md` |
| M8-03 This report | Draft, needs review | `TASKS.md` M8-03 |
| L0-01 Quantum adapter | Not implemented | `TASKS.md` |
| L3-01 Continuum adapter | Not implemented | `TASKS.md` |
| L4-01 Lumped system model | Not implemented | `TASKS.md` |

### Numbers that are unknown

The following numbers have no source in the repository. The report does not
invent them.

- The absolute accuracy of the thermal conductivity. No analytic reference
  exists.
- The absolute accuracy of the friction coefficient. The source calls it a
  simulated estimate.
- Any GPU timing. No GPU benchmark exists.
- Any medical or biological outcome. No such measurement exists.
- The number of tests at the M8 stage beyond the README count. The README
  states 425 Rust tests and 40 Python and MCP tests. Source: `README.md`. A
  later exact count is not recorded in the sources.

## Appendix A — Number traceability index

| Number | Value | Source |
|---|---|---|
| Crates | 9 | `TASKS.md` M0-01 |
| Rust tests | 425 | `README.md` |
| Python and MCP tests | 40 | `README.md` |
| License lines | 202 | `TASKS.md` M0-04 |
| Golden atoms | 6828 | `TASKS.md` M1-07 |
| Golden bonds | 10452 | `TASKS.md` M1-07 |
| MMP unit | 1e-13 m | `TASKS.md` M1-05 |
| FD step | 1e-14 m | `TASKS.md` M2-02 |
| Bond-stretch FD error | 1.28e-17 N | `TASKS.md` M2-02 |
| Angle-bend FD error | 8.3e-18 N | `TASKS.md` M2-03 |
| Torsion FD error | 7.5e-18 N | `TASKS.md` M2-04 |
| Out-of-plane FD error | 9.9e-17 N | `TASKS.md` M2-05 |
| van der Waals FD error | 7.0e-20 N | `TASKS.md` M2-06 |
| Electrostatic FD error | 1.7e-18 N | `TASKS.md` M2-07 |
| Total FD error | 5.7e-17 N | `TASKS.md` M2-08 |
| Anchor FD error | 4.8e-22 N | `TASKS.md` M3-04 |
| Spring FD error | 1.0e-16 N | `TASKS.md` M3-06 |
| Minimizer gradient | 5.9e-15 N | `TASKS.md` M2-09 |
| Minimizer iterations | 53 | `TASKS.md` M2-09 |
| NVE drift | 1.278e-5 | `TASKS.md` M2-10 |
| Water-box NVE drift | 1.274e-6 | `TASKS.md` M2-11 |
| Langevin mean | 300.58 K | `TASKS.md` M3-02 |
| Berendsen start and end | 10.23 K to 300.00 K | `TASKS.md` M3-03 |
| Rotor frames | 200000 | `TASKS.md` M3-07 |
| Rotor mean temperature | 309.57 K | `TASKS.md` M3-07 |
| Rotor max bond strain | 1.965e-2 | `TASKS.md` M3-07 |
| Rotor max force residual | 8.564e-26 N | `TASKS.md` M3-07 |
| Diamond lattice | 3.567e-10 m | `TASKS.md` M4-02 |
| Graphite lattice | 2.461e-10 m | `TASKS.md` M4-02 |
| Graphite layer | 3.354e-10 m | `TASKS.md` M4-02 |
| Graphene bond | 1.4209e-10 m | `TASKS.md` M4-03 |
| (5,5) nanotube cell | 20 atoms | `TASKS.md` M4-03 |
| (6,0) nanotube cell | 24 atoms | `TASKS.md` M4-03 |
| Sun teeth | 24 | `TASKS.md` M4-06 |
| Planet teeth | 18 | `docs/three-scale.md` |
| Ring teeth | 60 | `TASKS.md` M4-06 |
| Planet count | 3 | `TASKS.md` M4-06 |
| Design gear ratio | 3.5 | `TASKS.md` M4-06 |
| Measured gear ratio | 3.500000000 | `TASKS.md` M6-06 |
| Gear-ratio relative error | 1.146e-10 | `TASKS.md` M6-06 |
| Gear-ratio tolerance | 1.0e-3 | `TASKS.md` M6-06 |
| Elastic modulus | 4.49999999999e12 Pa | `TASKS.md` M5-02 |
| Analytic modulus | 4.5e12 Pa | `TASKS.md` M5-02 |
| Modulus relative error | 2.32e-12 | `TASKS.md` M5-02 |
| Friction coefficient | 0.4433 | `TASKS.md` M5-03 |
| Friction spread | 3.37e-12 N | `TASKS.md` M5-03 |
| Failure stress | 2.25e11 Pa | `TASKS.md` M5-04 |
| Rupture strain | 0.05 | `TASKS.md` M5-04 |
| Specific heat | 1.0414e3 J/(kg·K) | `TASKS.md` M5-05 |
| Analytic specific heat | 1.0393e3 J/(kg·K) | `TASKS.md` M5-05 |
| Specific-heat relative error | 2.0e-3 | `TASKS.md` M5-05 |
| Thermal conductivity | 8.34 +/- 0.93 W/(m·K) | `TASKS.md` M5-05 |
| Revolute axis error | 9.9956e-13 rad | `TASKS.md` M6-02 |
| Prismatic orientation error | 5.42e-20 rad | `TASKS.md` M6-03 |
| Gear-coupling error | 5.551e-17 m | `TASKS.md` M6-04 |
| Gear-coupling ratio | -0.500000000 | `TASKS.md` M6-04 |
| Constraint steps | 2000 | `TASKS.md` M6-02, M6-03 |
| Python tests | 20 | `TASKS.md` M7-01 |
| MCP tools | 18 | `TASKS.md` M7-02 |
| URDF links | 7 | `TASKS.md` M7-04 |
| URDF joints | 6 | `TASKS.md` M7-04 |
| Scene bodies | 7 | `TASKS.md` M8-01 |
| Scene joints | 6 | `TASKS.md` M8-01 |
| Scene planets | 3 | `TASKS.md` M8-01 |
| Scene module | 5e-10 m | `docs/three-scale.md` |
| Sun pitch radius | 6e-9 m | `docs/three-scale.md` |
| Planet pitch radius | 4.5e-9 m | `docs/three-scale.md` |
| Ring pitch radius | 1.5e-8 m | `docs/three-scale.md` |
| Carrier radius | 1.05e-8 m | `docs/three-scale.md` |
| Energy mean, run 1 | 680.97 µs | `benchmarks/results/2026-09-14-engine-timings.txt` |
| Serial mean, run 1 | 679.83 µs | `benchmarks/results/2026-09-14-engine-timings.txt` |
| Parallel mean, run 1 | 370.97 µs | `benchmarks/results/2026-09-14-engine-timings.txt` |
| Parallel speedup | about 1.83x | `benchmarks/results/2026-09-14-engine-timings.txt` |
| Energy mean, run 2 | 653.02 µs | `benchmarks/results/20260915T010005Z-engine-timings.txt` |
| Serial mean, run 2 | 660.21 µs | `benchmarks/results/20260915T010005Z-engine-timings.txt` |
| Parallel mean, run 2 | 303.08 µs | `benchmarks/results/20260915T010005Z-engine-timings.txt` |
| Water box atoms | 600 | `benchmarks/results/2026-09-14-engine-timings.txt` |
| Water box molecules | 200 | `benchmarks/results/2026-09-14-engine-timings.txt` |
| Water box side | 1.884e-9 m | `benchmarks/results/2026-09-14-engine-timings.txt` |
| Water box temperature | 300 K | `benchmarks/results/2026-09-14-engine-timings.txt` |
| Pair count | 30827 | `benchmarks/results/2026-09-14-engine-timings.txt` |
| Benchmark cores | 10 | `benchmarks/results/2026-09-14-engine-timings.txt` |
| MSRV | Rust 1.85 | `docs/reproduce.md` |
| Python minimum | 3.10 | `docs/reproduce.md` |
| Disk estimate | about 1 GB | `docs/reproduce.md` (estimate) |
