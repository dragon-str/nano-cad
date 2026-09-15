# Validation

This page describes the reference cross-checks that verify the nanocad atomistic
engine against independent molecular-dynamics codes. It states what each check
proves, what each check does not prove, and how to run each check.

`ARCHITECTURE.md` names two external cross-checks: OpenMM and GROMACS CPU.

## The common reference model

Both cross-checks use the same fixed case: a periodic SPC-like water box from
`crates/engine/src/reference.rs`, printed as JSON by
`crates/engine/examples/water_box.rs`.

- 200 molecules, 600 atoms.
- Cubic box, side `1.8833259369189425e-9 m`, periodic.
- Cutoff `0.35 * box`, switch-on `0.9 * cutoff`, no long-range correction.
- Harmonic O-H bond: `k = 500 N/m`, `r0 = 1e-10 m`.
- Harmonic H-O-H angle: `k = 6e-19 J/rad^2`, `theta0 = 109.47 deg`.
- Buckingham O-O term: `U = A exp(-B r) - C/r^6`, with `A = 1.151e-15 J`,
  `B = 4.2e10 1/m`, `C = 2.6645e-78 J m^6`. Only oxygen carries vdW parameters.
- Charges: `q_O = -0.82 e`, `q_H = +0.41 e`.
- The 3 intra-molecular pairs per molecule are excluded.
- Seed `0x5A17E2` fixes the orientations and the velocities.

The model is flexible and SPC-like. It is not standard SPC, because the bonded
terms are flexible and the dispersion term is Buckingham.

## OpenMM cross-check

Task M2-11. Driver: `benchmarks/openmm/compare_water_box.py`. Result:
`benchmarks/results/openmm-water-crosscheck.txt`.

The driver runs the Rust example, rebuilds the identical model in OpenMM from
the JSON fields, and compares the potential energy and the per-atom forces.

| Quantity | nanocad | OpenMM | Relative difference |
|---|---|---|---|
| Potential energy | -2.737699062416127e-17 J | -2.7376990624161912e-17 J | 2.341e-14 |
| Maximum per-atom force | — | max magnitude 1.4621685148513376e-08 N | 1.242e-14 |

The declared tolerance is `1.0e-6` relative. The result is **PASS**.

### What the OpenMM cross-check proves

- The bonded terms, the Buckingham term, the Coulomb term, the switching
  function, the exclusions, and the minimum-image convention in nanocad agree
  with an independent implementation to double precision.
- A coding error in any of those terms shows up as a difference far above the
  tolerance. A missing switching function on the dispersion term is about 5e-4.

### What the OpenMM cross-check does not prove

- It is a coding check, not a physics validation. Both codes use the same plain
  cutoff and no long-range correction. Agreement does not mean the physics is
  correct.
- It does not compare nanocad to a standard SPC or TIP3P result, because the
  model is flexible and Buckingham-based.
- It checks one fixed configuration. A model error that is symmetric under this
  geometry could pass.
- The switching branch at exactly `switch_on` is not exercised, because no pair
  sits exactly at that distance.
- The `CPU` and `OpenCL` OpenMM platforms are single precision in places and
  agree to about `4e-6`, not `1e-14`. Certify with the `Reference` platform.

### How to run the OpenMM cross-check

```
uv venv --python 3.12 /tmp/openmm-venv
uv pip install --python /tmp/openmm-venv/bin/python openmm numpy
/tmp/openmm-venv/bin/python benchmarks/openmm/compare_water_box.py
```

The `Reference` platform is the default. Set `NANOCAD_OPENMM_PLATFORM=CPU` for
the single-precision CPU platform. The driver prints `SKIP` and exits zero when
OpenMM is absent.

## GROMACS CPU cross-check

Task from `ARCHITECTURE.md`. Driver:
`benchmarks/gromacs/compare_water_box.py`. Result:
`benchmarks/gromacs/results/gromacs-water-crosscheck.txt`.

### Current result: SKIP

GROMACS 2026.3 is installed (Homebrew build, `2026.3-Homebrew`). It cannot run
the model. `gmx mdrun` reports:

```
Fatal error:
The Verlet cutoff-scheme does not (yet) support Buckingham
```

Therefore the driver prints `SKIP`, records no agreement, and exits zero. It
does not fabricate a comparison.

The reasons are:

- GROMACS removed the group cutoff scheme in 2020. Buckingham was only ever
  supported with that scheme.
- GROMACS 2026 implements only the potential-shift modifier, not a switching
  function. The nanocad custom switching polynomial has no GROMACS equivalent.
- A tabulated nonbonded potential (`vdwtype = User`) is also unavailable with
  the Verlet scheme in 2026. The GROMACS developers state that it returns in the
  2027 release.

### What the GROMACS driver does verify

- `gmx grompp` accepts the exact topology: `nbfunc = 2` (Buckingham), the same
  charges, bond, angle, exclusions, box, and cutoff. The SI-to-GROMACS input
  conversion is therefore correct.
- The Coulomb constant converts exactly. GROMACS prints `138.935458` and the
  nanocad value is `138.935458 kJ/mol/nm/e^2` at the printed precision.
- The switch-correction function matches an independent Python evaluation of
  the reference model. Bonded plus switched non-bonded energy equals the
  recorded nanocad energy to `6.75e-15` relative. Reproduce this with
  `python3 benchmarks/gromacs/selftest_switch.py`.

### What the GROMACS cross-check would prove

When GROMACS can run the model, the driver compares the potential energy and the
per-atom forces, exactly as the OpenMM driver does. It would give a second,
independent CPU implementation to confirm the nanocad terms. It is subject to
the same limits as the OpenMM check: a coding check, not a physics validation,
and one fixed configuration only.

The driver adds an exact analytic switch correction, because GROMACS has no
nanocad switch. It compares GROMACS to the corrected prediction. The correction
is validated by the independent evaluation above.

### How to run the GROMACS cross-check

```
brew install gromacs
python3 benchmarks/gromacs/compare_water_box.py
```

The driver prints `SKIP` and exits zero when GROMACS is absent, and also when
GROMACS is present but cannot express the model. Use `--workdir DIR --keep` to
inspect the generated `.g96`, `.top`, and `.mdp` input.

## Other verification mechanisms

The reference cross-checks are one part of the verification architecture in
`ARCHITECTURE.md`. The other parts run in the normal test suite.

- **Finite-difference gradient test per force term.** Every analytic force term
  is compared to a central finite difference in `cargo test`. A term does not
  merge without one.
- **Energy-conservation test per integrator.** An NVE run must keep the total
  energy within a set bound.
- **Golden regression files** for import and export, under `tests/golden/`.

These checks run with `just verify`. They are merge gates. The OpenMM and
GROMACS cross-checks need an external code, so they are not gates.

## What no check proves

- No check validates the model against experiment. All results are **simulated**.
- No check shows medical or clinical validity.
- Agreement between two codes proves that the two implementations match. It does
  not prove that the model is physically correct.
