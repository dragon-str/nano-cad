# PLAN.md — Program Plan

## Vision

Give anyone the power to design atomically precise machines and to carry those
designs across scales, from a single bond up to a working device and its
environment.

## Mission

Build an open, multiscale CAD and simulation platform for molecular
nanotechnology and nanomedicine. Make it trustworthy, scriptable, and agent
driven. Keep the designs and the formats open.

## Users

- Nanotechnology researchers who design molecular machines.
- Computational chemists who need a design front end.
- Biomedical researchers who study nanomedicine.
- Educators who teach molecular design.
- Coding agents that design and check parts through an API.

## Principles

1. One headless core serves the GUI, the CLI, and the agent.
2. The multiscale parameter schema is the central contract.
3. Vertical slices before breadth.
4. Verification is a merge gate.
5. SI units with explicit types.
6. Reuse OpenMM and ASE. Do not rebuild solved software.
7. Open formats and open designs.
8. No medical claims.

## Scope

In scope for the first release:

- Atomistic engine: bonds, angles, torsions, van der Waals, electrostatics.
- Minimization and molecular dynamics with thermostats.
- Jigs: anchors, motors, springs, thermostats, applied forces.
- Rigid bodies and joints at the device level.
- A parametric part library, starting with a gear generator.
- Parameter extraction from atomistic parts to the device level.
- A versioned open document format and MMP interop.
- A Python API, a CLI, and an MCP server for agents.
- One demo: an atomically precise planetary gearbox.

## Non-goals for the first release

- A whole-body human twin.
- Many force fields and many QM codes at once.
- Long-range electrostatics beyond a cutoff (add PME later through OpenMM).
- Production-scale MD of millions of atoms.
- Closed formats or exclusive designs.
- Any claim of a built or medically validated device.

## Phases and milestones

Each milestone has an exit criterion. A milestone is done only when the
criterion holds and its tasks in `TASKS.md` are checked.

### Phase 0 — Foundation (M0, M1)

- **M0 Skeleton.** Workspace, CI, docs, units module, license headers, benchmark
  harness stub.
  - Exit: `just verify` runs green on an empty workspace.
- **M1 Document model and formats.** Atoms, bonds, topology, a versioned NCZ
  format, XYZ and MMP import and export, round-trip tests.
  - Exit: A real MMP file imports, saves to NCZ, reloads, and reproduces the
    original coordinates bit for bit.

### Phase 1 — Atomistic core (M2, M3)

- **M2 Forces and minimization.** Cell list, flat pair list, bond, angle,
  torsion, van der Waals, and electrostatic terms. A minimizer. A
  finite-difference gradient test for every term. An energy-conservation run.
  - Exit: Gradients match finite differences to tolerance. NVE energy drifts
    below the set bound. A water-box benchmark runs and is compared to OpenMM.
- **M3 Dynamics and jigs.** Velocity Verlet, Langevin or Berendsen thermostat,
  anchors, motors, springs.
  - Exit: A test machine runs stable dynamics at 300 K. A motor drives a rotor.

### Phase 2 — Part generation (M4)

- **M4 Gear generator and part library.** A parametric gear generator, a
  nanotube builder, a lattice builder, and a part schema.
  - Exit: The agent builds a specified gear by name and specs. The part passes
    a strain and valence check.

### Phase 3 — Multiscale bridge (M5, M6)

- **M5 Parameter extraction.** Relax a part, measure stiffness, friction,
  failure stress, and thermal properties. Write a parameter record with
  provenance and uncertainty. See `PARAMETERS.md`.
  - Exit: A gear yields a parameter record that passes a consistency check.
- **M6 Device level.** Rigid-body dynamics with joints. A gearbox assembled
  from atomistic parts. Handoff from L1 parameters to L2 bodies.
  - Exit: The gearbox reproduces its designed gear ratio within tolerance.

### Phase 4 — Agent interface (M7)

- **M7 Headless API and MCP server.** Typed tools: create part, generate gear,
  relax, run MD, add jig, measure, convert units, save, load, list.
  - Exit: An agent completes the gearbox demo through tools alone.

### Phase 5 — Demo and release (M8)

- **M8 Demo.** The 90-second gearbox video and a written report. The three-scale
  handoff animation. A public benchmark.
  - Exit: A stranger reproduces the demo from the repository.

### Phase 6 — Flagship and upper levels (later)

- Respirocyte subsystems: rotor, bearing, pump, gas tank.
- L3 continuum adapters and L4 lumped system models.
- One nanomedicine slice: a machine coupled to a coarse host model.

## The vertical slice

The first end-to-end path is the planetary gearbox. It touches M1, M2, M3, M4,
M5, and M6 thinly. Build the slice first. Deepen each level after the slice
runs. This prevents the breadth failure of the original project.

## Success metrics

- The gearbox demo runs end to end from a clean clone.
- Every force term has a passing finite-difference test.
- The engine matches OpenMM on a reference water box within tolerance.
- An agent builds a specified part through tools alone.
- A parameter record carries full provenance and uncertainty.

## Risks

| Risk | Impact | Mitigation |
|---|---|---|
| Scope creep | project death | vertical slice, hard non-goals |
| GPL contamination | legal exposure | clean room, no copying, license checks |
| Wrong physics | lost trust | finite-difference gates, reference benchmarks |
| Overclaiming | reputation damage | label all output as simulated |
| Tiny market | no revenue | serve adjacent fields, seek grants |
| Long timeline | burnout | phase exits, visible demos |
| Agent errors | bad designs | typed tools, validation, human review |

## Funding path

- Near term: grants, consulting, and sponsored development. See the earlier
  analysis. National science agencies, private foundations, and long-termist
  philanthropies are the likely sources.
- Mid term: support, training, certification, and hosted runs.
- Long term: a verified parameter library and a certified platform.

## Honest note

This is a multi-year program. The first usable tool is 1 to 3 years out. The
respirocyte study is a research artifact, not a product. Plan for steady support
and visible milestones.
