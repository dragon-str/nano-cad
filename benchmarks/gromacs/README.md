# GROMACS water-box cross-check

This directory holds the GROMACS CPU cross-check named in `ARCHITECTURE.md`
("Reference benchmarks against OpenMM and GROMACS CPU"). It targets the same
periodic SPC-like water box as the OpenMM cross-check in `benchmarks/openmm/`
and the same fixed configuration (200 molecules, seed `0x5A17E2`).

The check runs one fixed, deterministic configuration. It is not a timing
benchmark and it is **not** part of `just verify`.

## Result on this machine

GROMACS 2026.3 is installed. It **cannot run the model**, so the driver prints
`SKIP` and exits zero. The recorded run is
[`results/gromacs-water-crosscheck.txt`](results/gromacs-water-crosscheck.txt).

## Run the check

```
brew install gromacs
python3 benchmarks/gromacs/compare_water_box.py
```

The driver writes `.g96`, `.top` and `.mdp` input to a temporary directory, runs
`gmx grompp` and `gmx mdrun`, and compares the energy and forces. Use
`--workdir DIR --keep` to inspect the generated input.

When GROMACS is absent the driver prints `SKIP` and exits zero. It is never a
hard gate.

## Why GROMACS cannot run the model

The nanocad vdW term is Buckingham: `U = A exp(-B r) - C / r^6`. The installed
GROMACS reports:

```
Fatal error:
The Verlet cutoff-scheme does not (yet) support Buckingham
```

The facts behind that message:

- GROMACS removed the group cutoff scheme in 2020. Buckingham was only ever
  supported with that scheme.
- GROMACS 2026 implements only the potential-shift modifier, not a switching
  function. The reference manual keeps the switch formula for a future release,
  but the code does not apply it. The nanocad switching polynomial therefore has
  no GROMACS equivalent.
- A tabulated nonbonded potential (`vdwtype = User`) is also unavailable with
  the Verlet scheme in 2026. The GROMACS developers state that it returns in the
  2027 release.

The driver cannot record a measured agreement, because GROMACS produced no
energy and no forces. It does **not** fabricate a number.

## What the driver still verifies

1. `gmx grompp` accepts the exact topology: `nbfunc = 2` (Buckingham), the same
   charges, the same harmonic bond and angle, 3 excluded intra-molecular pairs
   per molecule, and the same cubic box and cutoff. The input conversion is
   therefore correct.
2. The Coulomb constant converts exactly. GROMACS prints `138.935458` and the
   nanocad value is `138.935458 kJ/mol/nm/e^2` at the printed precision.
3. The switch correction is consistent with an independent Python evaluation of
   the reference model. Bonded plus switched non-bonded energy equals the
   recorded nanocad energy to `6.75e-15` relative. The correction equals the
   plain-cutoff minus switched value to machine precision. Reproduce this with
   `python3 benchmarks/gromacs/selftest_switch.py`. That script does not need
   GROMACS.

## What the driver compares when GROMACS can run the model

The mdp uses `vdw-modifier = None` and `coulomb-modifier = None`, because
GROMACS has no nanocad switch. The driver then adds the exact analytic switch
difference to the nanocad model:

```
E_gromacs = E_nanocad + sum over pairs in [switch_on, cutoff] of V_full(r) (1 - S(r))
F_gromacs = F_nanocad - grad(E_gromacs - E_nanocad)
```

The driver computes this correction in double precision from the same positions,
box, cutoff, charges and vdW parameters. It compares GROMACS to the corrected
prediction, so the comparison isolates coding errors from the one documented
model difference. If GROMACS gains Buckingham and switch support, the same
driver runs without the correction.

## Units

nanocad is SI. GROMACS is nanometre and kilojoule per mole. The driver converts
at the edges:

```
1 nm         = 1e-9 m
1 kJ/mol     = 1000 / N_A J
1 kJ/mol/nm  = (1000 / N_A) / 1e-9 N
```

The Buckingham parameters convert as `A [kJ/mol] = A [J] N_A/1000`,
`B [1/nm] = B [1/m] 1e-9`, and `C [kJ/mol/nm^6] = C [J m^6] 1e9^6 N_A/1000`.
The harmonic constants convert to `kJ/mol/nm^2` and `kJ/mol/rad^2` with the
same `N_A/1000` factor.

## Caveats

- The recorded result is a skip, not an agreement.
- GROMACS 2026.3 is single precision in places. If it could run the model, the
  expected agreement is about `1e-5`, not the `1e-14` that OpenMM `Reference`
  reaches. The declared tolerance is `1e-4`.
- The check compares model implementations. It does not test physics against
  experiment, and it uses one fixed configuration.
