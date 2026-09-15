# OpenMM water-box cross-check

This directory holds the cross-check for task **M2-11**: the periodic SPC-like
water box in `crates/engine/src/reference.rs` compared to OpenMM for the
potential energy and the per-atom forces.

The check runs one fixed, deterministic configuration. It is not a timing
benchmark and it is **not** part of `just verify`. When OpenMM is absent the
driver prints `SKIP` and exits zero.

## Install OpenMM

The system Python 3.14 has no OpenMM wheel. Use `uv` to make a 3.12
environment:

```
uv venv --python 3.12 /tmp/openmm-venv
uv pip install --python /tmp/openmm-venv/bin/python openmm numpy
```

Recorded versions for the result in `benchmarks/results/`: Python 3.12.9 and
OpenMM 8.6.1 (`openmm.version.version` reports `8.6.1.dev-b399af4`).

## Run the check

```
/tmp/openmm-venv/bin/python benchmarks/openmm/compare_water_box.py
```

The `Reference` platform runs in double precision and is the default. To run a
different platform, set `NANOCAD_OPENMM_PLATFORM`, for example
`NANOCAD_OPENMM_PLATFORM=CPU`. The tolerance is `1e-6` on `Reference` and `1e-5`
on the single-precision platforms.

## What it compares

1. It runs `cargo run -q --example water_box` in `crates/engine/examples/`.
   The example builds the box with seed `0x5A17E2` and prints JSON: the box, the
   cutoff, every physical parameter, the positions in metres, the potential
   energy in joules, and the per-atom forces in newtons.
2. It rebuilds the identical model in OpenMM from those fields:
   `CustomBondForce`, `CustomAngleForce`, and two `CustomNonbondedForce` objects
   for the Buckingham oxygen-oxygen term and the Coulomb term. Both non-bonded
   expressions carry the same polynomial switching function as
   `crates/engine/src/nonbonded.rs`. The three intra-molecular pairs are
   excluded. The method is `CutoffPeriodic` with no built-in switching and no
   long-range correction.
3. It converts OpenMM output from nanometre and kilojoule per mole back to SI,
   then prints the relative energy difference and the maximum per-atom force
   difference.

## Units

nanocad is SI. OpenMM is nanometre and kilojoule per mole. The driver converts
at the edges:

```
1 nm        = 1e-9 m
1 kJ/mol    = 1000 / N_A J
1 kJ/mol/nm = (1000 / N_A) / 1e-9 N
```

The Coulomb constant is computed from the nanocad value, not from the OpenMM
default. The two agree to the OpenMM print precision.

## Caveats

- The `step(ron - r)` guard in the switching expression makes the value exactly
  one at and below `switch_on`. nanocad uses the same branch. The configuration
  holds no pair exactly at `switch_on`, so the branch is not exercised.
- The OpenMM `CPU` and `OpenCL` platforms use single precision in parts of the
  force evaluation. Their agreement is about `4e-6`, not `1e-14`. Certify the
  result with `Reference`.
- The model has no long-range correction, no Ewald sum, and no explicit
  dispersion tail. Both codes use the same plain cutoff, so the comparison is
  self-consistent. It does not test either code against experiment.
