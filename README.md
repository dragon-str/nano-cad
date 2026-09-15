# Nano-CAD

An open, multiscale, atomically precise CAD and simulation platform for
molecular nanotechnology and nanomedicine.

Status: **pre-alpha.** The core is built end to end and **440 Rust tests plus
75 Python and MCP tests pass**. The 9-crate workspace has: units, the model and
the NCZ/MMP/XYZ formats, a force engine (bond, angle, torsion, out-of-plane,
van der Waals, electrostatic, periodic minimum image, velocity Verlet,
conjugate gradient, Langevin and Berendsen thermostats), part generators
(diamond, graphite, nanotube, involute gear, spur gear, planetary set,
respirocyte rotor, pump and tank), parameter records and extraction (stiffness,
friction, failure, thermal), jigs, rigid-body device joints, URDF and
three-scale scene export, the CLI, PyO3 bindings, and an MCP server. An agent
builds the planetary gearbox through the tool surface alone and measures the
analytic 3.5 gear ratio. A clean checkout reproduces the demo with
`just reproduce`. The periodic water box matches OpenMM 8.6.1 to 1e-14
relative. The later-phase adapters are started: L0 quantum (PySCF), L3
continuum (structural and flow), L4 lumped system (SBML subset), FL-01, FL-02
and NM-01. M8-02 needs the video recording. See `docs/benchmarks.md` and
`docs/report.md`.

## Doctrine

- Open core under Apache-2.0. Open formats. A governed verified-data library.
- One headless core. The agent, the GUI, and scripts all call the same API.
- Vertical slices before breadth. One part, end to end, before a whole layer.
- Verification first. Every force term needs a finite-difference gradient test.
- SI units internally, with explicit unit types. Never mix Angstrom and nanometre.
- No medical claims. Output is a design hypothesis, not a validated device.

## Documents

Read in this order.

| File | Purpose |
|---|---|
| `AGENTS.md` | Operating manual for a coding agent. Start here. |
| `PLAN.md` | Vision, scope, phases, milestones, risks, funding. |
| `ARCHITECTURE.md` | Layers, repo layout, data model, stack, interfaces. |
| `PARAMETERS.md` | The multiscale parameter schema. The central contract. |
| `TASKS.md` | The iterative backlog. Pick the next unchecked task. |
| `DECISIONS.md` | The decision log (ADRs). Record decisions here. |
| `docs/benchmarks.md` | Public benchmark page: hardware, timing, caveats. |

## Reproduce the demo

A clean checkout reproduces the planetary gearbox demo with one command. The
first run needs network access to download `maturin` and `pytest`.

```sh
just reproduce
```

The script checks the toolchains, builds the extension with maturin from
`python/`, runs the Rust tests, runs the MCP demo test, and prints the measured
gear ratio. See `docs/reproduce.md` for the prerequisites, the exact commands,
what each step proves, and what is not reproduced.

## Start a clean session

Point a fresh agent at this repository and give it this prompt:

> Read AGENTS.md, then PLAN.md, ARCHITECTURE.md, PARAMETERS.md, TASKS.md, and
> DECISIONS.md. Work on the first unchecked task in TASKS.md. Follow the work
> loop in AGENTS.md. Do not copy GPL code. Run the verification commands before
> you report done.

## License

Apache-2.0. See `LICENSE`. The upstream NanoEngineer code is GPLv2. Do not copy
it. Read it for algorithms only, then reimplement.
