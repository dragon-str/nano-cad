# AGENTS.md — Operating Manual

This file tells a coding agent how to work in this repository. Read it first.

## Read order

1. `AGENTS.md` (this file)
2. `PLAN.md` — what we are building and why
3. `ARCHITECTURE.md` — how it is built
4. `PARAMETERS.md` — the central data contract
5. `TASKS.md` — the backlog
6. `DECISIONS.md` — the decision log

## Ground rules

1. **License.** All code is Apache-2.0. Never copy GPL code. The upstream
   NanoEngineer C engine is GPLv2. You may read it for algorithms. You must
   reimplement. During development it lives at `/tmp/nanoengineer`.
2. **Units.** Use SI internally. Use the units module. Never mix Angstrom and
   nanometre. Put the unit in the variable name when it is ambiguous, for
   example `length_m` or `time_s`.
3. **Verification.** Every analytic force term needs a finite-difference
   gradient test. A term does not merge without one.
4. **No overclaiming.** Write "simulated", never "built", "proven", or
   "validated for medicine". Mark estimates as estimates.
5. **Vertical slices.** Build a thin end-to-end path before a broad layer.
   The first slice is the planetary gearbox. See `PLAN.md`.
6. **Small commits.** One task per commit. Keep the tree green.

## The work loop

1. Open `TASKS.md`. Pick the first unchecked task in the current milestone.
   If the user named a task, do that one.
2. Mark it in progress. Append ` <!-- WIP -->` to the task line.
3. Implement the task. Add tests at the same time.
4. Run the verification commands below. Fix every failure.
5. Update `TASKS.md`: change `[ ]` to `[x]`, remove `<!-- WIP -->`, and append
   a one-line result note.
6. If you made an architectural decision, add an entry to `DECISIONS.md`.
7. Commit. The message says what changed and which task it closes.
8. Report to the user in Simplified Technical English. State the route.

## Verification commands

Create these as the tree grows. They are the gate.

```
cargo fmt --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo bench            # timing only, not a gate
pytest python/tests
just verify            # runs all gates
```

## Definition of done

A task is done only when all of these are true.

- The code compiles with no warnings.
- Tests are added and pass.
- A finite-difference gradient test exists for any new force term.
- The energy-conservation check passes for any new integrator or thermostat.
- The relevant docs are updated.
- `TASKS.md` shows the task checked with a result note.
- No GPL code was copied.

## Code conventions

- **Rust.** `rustfmt`, `clippy`. Library crates use `#![forbid(unsafe_code)]`
  except a single audited FFI crate. Use `thiserror` in libraries. Do not panic
  in library code.
- **Python.** `ruff` and `black`. Type hints on public functions.
- **Naming.** `snake_case`. Suffix quantities with their unit.
- **Comments.** Explain why, not what. Keep them short.
- **Tests.** Unit tests beside the code. Golden files under `tests/golden/`.
  Benchmarks under `benchmarks/`.

## What not to touch

- Do not add credentials, tokens, or personal data to this repository.
- Do not add a dependency without checking its license. Record it in
  `ARCHITECTURE.md`.
- Do not change the parameter schema without a version bump and an ADR.

## Reporting to the user

- Open every reply with one route line, for example `-> Claude Opus`.
- Write in Simplified Technical English. Short sentences. Active voice.
- Report what is done, what is verified, and what is still open.
- Never report a delegated change as done before you read the diff.
