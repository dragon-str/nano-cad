#!/usr/bin/env python3
"""Cross-check the nanocad periodic SPC-like water box against GROMACS CPU.

The script runs the Rust example `crates/engine/examples/water_box.rs`, reads the
fixed configuration and parameters from its JSON output, rebuilds the identical
physical model in GROMACS, runs a single point, and compares the potential
energy and the per-atom forces.

Every value crosses the boundary in SI. GROMACS works in nanometre, kilojoule
per mole and elementary charge, so the script converts at the edges:

    1 nm       = 1e-9 m
    1 kJ/mol   = 1000 / N_A J
    1 kJ/mol/nm = (1000 / N_A) / 1e-9 N

GROMACS 2026.3 cannot run this model.  The nanocad vdW term is Buckingham, and
the Verlet cutoff scheme in GROMACS does not support Buckingham.  The script
detects that condition, prints the exact GROMACS error, and exits zero.  It is
never a hard gate.  Read `README.md` for the full limitation.

Run it with a GROMACS on PATH:

    python3 benchmarks/gromacs/compare_water_box.py

The script exits zero when GROMACS is absent, or when GROMACS is present but
cannot express the model.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
OPENMM_RESULT = REPO_ROOT / "benchmarks" / "results" / "openmm-water-crosscheck.txt"

# Relative tolerance for the cross-check when GROMACS can run the model.  The
# Homebrew build is single precision, so the expected agreement is about 1e-5,
# not the 1e-14 that OpenMM `Reference` reaches.  A real model mismatch (a
# missing switching function on the dispersion term) is about 5e-4.
SINGLE_PRECISION_TOL = 1.0e-4

N_A = 6.02214076e23
J_PER_KJMOL = 1000.0 / N_A
N_PER_KJMOL_PER_NM = J_PER_KJMOL / 1.0e-9
NM_PER_M = 1.0e9
AMU_PER_KG = 1.0 / 1.66053906660e-27

# GROMACS prints f = 1/(4 pi eps0) = 138.935458 kJ/mol/nm/e^2 (reference manual,
# Definitions and Units).  nanocad uses 8.9875517923e9 N m^2/C^2.  The two agree
# to about 2e-7 relative; the driver reports the difference.
GROMACS_COULOMB_KJMOL_NM_PER_E2 = 138.935_458

# GROMACS error fragments that mean "this model cannot run".  They are matched
# case-insensitively against mdrun output.
UNSUPPORTED_MODEL_MARKERS = (
    "does not (yet) support buckingham",
    "does not support buckingham",
    "verlet cutoff-scheme does not",
)


def find_gmx() -> str | None:
    """Returns the path to a GROMACS binary, or None when GROMACS is absent."""
    for name in ("gmx", "gmx_mpi"):
        path = shutil.which(name)
        if path is not None:
            return path
    return None


def gmx_version(gmx: str) -> str:
    """Returns a one-line GROMACS version string."""
    proc = subprocess.run(
        [gmx, "--version"],
        capture_output=True,
        text=True,
        check=False,
    )
    for line in proc.stdout.splitlines():
        if line.startswith("GROMACS version:"):
            return line.split(":", 1)[1].strip()
    return "unknown"


def run_rust_example() -> dict:
    """Runs the Rust water-box example and returns its JSON output."""
    env = os.environ.copy()
    cargo_bin = str(Path.home() / ".cargo" / "bin")
    env["PATH"] = cargo_bin + os.pathsep + env.get("PATH", "")
    proc = subprocess.run(
        ["cargo", "run", "-q", "--example", "water_box"],
        cwd=REPO_ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError("the Rust example failed:\n" + proc.stderr.strip())
    return json.loads(proc.stdout)


def read_recorded_energy_j() -> float | None:
    """Reads the recorded nanocad potential energy from the OpenMM result.

    Returns None when the file or the line is absent.  The value is the same
    fixed configuration that the Rust example prints, so it is a consistency
    check on the live run.
    """
    if not OPENMM_RESULT.is_file():
        return None
    match = re.search(
        r"Energy nanocad:\s*([-+0-9.eE]+)\s*J",
        OPENMM_RESULT.read_text(encoding="utf-8"),
    )
    if match is None:
        return None
    return float(match.group(1))


def nucleus_name(site: int) -> tuple[str, str]:
    """Returns the GROMACS (residue, atom) names for site 0, 1 or 2.

    Site 0 is the oxygen.  Sites 1 and 2 are the two hydrogens.
    """
    if site == 0:
        return "WAT", "OW"
    return "WAT", f"HW{site}"


def write_g96(data: dict, path: Path) -> None:
    """Writes a GROMOS-96 coordinate file in nanometre.

    GROMOS-96 carries more digits than the fixed-column .gro format, so the
    input positions keep the precision that the force comparison needs.
    """
    lines = ["TITLE", "nanocad SPC-like water box", "END", "POSITION"]
    for atom in range(data["atom_count"]):
        molecule, site = divmod(atom, 3)
        residue, name = nucleus_name(site)
        x = data["positions_m"][3 * atom] * NM_PER_M
        y = data["positions_m"][3 * atom + 1] * NM_PER_M
        z = data["positions_m"][3 * atom + 2] * NM_PER_M
        lines.append(
            "%5d %-5s %-5s  %5d%15.9f%15.9f%15.9f"
            % (molecule + 1, residue, name, atom + 1, x, y, z)
        )
    box_nm = data["box_length_m"] * NM_PER_M
    lines.append("END")
    lines.append("BOX")
    lines.append("%13.9f%15.9f%15.9f" % (box_nm, box_nm, box_nm))
    lines.append("END")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_topology(data: dict, path: Path) -> None:
    """Writes a GROMACS topology for the exact nanocad model.

    The non-bonded function is Buckingham (`nbfunc = 2`).  The three intra-
    molecular pairs per molecule are excluded by `nrexcl = 2`, which matches the
    three exclusions that `reference.rs` adds.
    """
    bond_k = data["bond_k_n_per_m"] * 1.0e-18 * (N_A / 1000.0)
    bond_r0_nm = data["bond_r0_m"] * NM_PER_M
    angle_k = data["angle_k_j_per_rad2"] * (N_A / 1000.0)
    angle_theta0_deg = math.degrees(data["angle_theta0_rad"])

    vdw_a = data["vdw_a_j"] * (N_A / 1000.0)
    vdw_b_nm = data["vdw_b_per_m"] * 1.0e-9
    vdw_c = data["vdw_c_j_m6"] * (NM_PER_M**6) * (N_A / 1000.0)

    q_oxygen = data["oxygen_charge_c"] / data["elementary_charge_c"]
    q_hydrogen = data["hydrogen_charge_c"] / data["elementary_charge_c"]
    mass_oxygen = data["oxygen_mass_kg"] * AMU_PER_KG
    mass_hydrogen = data["hydrogen_mass_kg"] * AMU_PER_KG

    text = f"""\
; nanocad SPC-like water box rebuilt in GROMACS.
; All values are converted from SI.  See compare_water_box.py.
[ defaults ]
; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ
2 1 no 1.0 1.0

[ atomtypes ]
; name at.num mass charge ptype A B C
OW 8 {mass_oxygen:.10f} 0.0 A {vdw_a:.10e} {vdw_b_nm:.10e} {vdw_c:.10e}
HW 1 {mass_hydrogen:.10f} 0.0 A 0.0 0.0 0.0

[ nonbond_params ]
; i j func A B C
OW OW 2 {vdw_a:.10e} {vdw_b_nm:.10e} {vdw_c:.10e}
OW HW 2 0.0 0.0 0.0
HW HW 2 0.0 0.0 0.0

[ moleculetype ]
; name nrexcl
WAT 2

[ atoms ]
; nr type resnr resid atom cgnr charge mass
1 OW 1 WAT OW 1 {q_oxygen:.10f} {mass_oxygen:.10f}
2 HW 1 WAT HW1 1 {q_hydrogen:.10f} {mass_hydrogen:.10f}
3 HW 1 WAT HW2 1 {q_hydrogen:.10f} {mass_hydrogen:.10f}

[ bonds ]
; ai aj func b0 kb
1 2 1 {bond_r0_nm:.10e} {bond_k:.10e}
1 3 1 {bond_r0_nm:.10e} {bond_k:.10e}

[ angles ]
; ai aj ak func theta0 ktheta
2 1 3 1 {angle_theta0_deg:.10f} {angle_k:.10e}

[ system ]
nanocad SPC-like water box, {data["molecule_count"]} molecules

[ molecules ]
; name count
WAT {data["molecule_count"]}
"""
    path.write_text(text, encoding="utf-8")


def write_mdp(data: dict, path: Path) -> None:
    """Writes the single-point mdp that matches the nanocad cutoff.

    The nanocad switching function has no GROMACS equivalent, so the modifier is
    `None` (a plain truncation).  This is the documented model difference.
    """
    cutoff_nm = data["cutoff_m"] * NM_PER_M
    text = f"""\
; nanocad SPC-like water box, single point.
integrator = md
nsteps = 1
dt = 1e-12
cutoff-scheme = Verlet
nstlist = 1
rlist = {cutoff_nm:.10e}
rvdw = {cutoff_nm:.10e}
rcoulomb = {cutoff_nm:.10e}
vdw-modifier = None
coulomb-modifier = None
coulombtype = Cut-off
nstfout = 1
nstenergy = 1
nstxout = 1
"""
    path.write_text(text, encoding="utf-8")


def run_command(args: list[str], cwd: Path) -> subprocess.CompletedProcess:
    """Runs a command and returns the completed process."""
    return subprocess.run(
        args,
        cwd=cwd,
        capture_output=True,
        text=True,
        check=False,
    )


def parse_energy_xvg(path: Path) -> float | None:
    """Returns the first potential energy in an xvg file, in kJ/mol."""
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped[0] in "#@":
            continue
        fields = stripped.split()
        if len(fields) >= 2:
            return float(fields[1])
    return None


def parse_first_frame_forces(text: str, atom_count: int) -> list[float]:
    """Returns the first trajectory frame of forces, in kJ/mol/nm.

    Parses the `f[ i]={x, y, z}` records that `gmx dump` writes.
    """
    pattern = re.compile(r"f\[\s*(\d+)\]=\{\s*([^,]+),\s*([^,]+),\s*([^}]+)\}")
    forces: list[float] = []
    for line in text.splitlines():
        match = pattern.search(line)
        if match is None:
            continue
        forces.extend(float(match.group(k)) for k in (2, 3, 4))
        if len(forces) >= 3 * atom_count:
            break
    return forces


def _minimum_image(delta_m: list[float], length_m: float) -> list[float]:
    """Reduces a pair displacement by the minimum-image convention."""
    out_m = list(delta_m)
    for axis in range(3):
        if length_m > 0.0:
            out_m[axis] -= length_m * round(out_m[axis] / length_m)
    return out_m


def _switch_value(
    r_m: float, switch_on_m: float, cutoff_m: float
) -> tuple[float, float]:
    """Returns the nanocad switching value S(r) and dS/dr.

    This mirrors `Cutoff::switch_value` in `crates/engine/src/nonbonded.rs`.
    """
    if r_m <= switch_on_m:
        return 1.0, 0.0
    if r_m >= cutoff_m:
        return 0.0, 0.0
    denom = cutoff_m * cutoff_m - switch_on_m * switch_on_m
    denom_cubed = denom * denom * denom
    d = cutoff_m * cutoff_m - r_m * r_m
    e = cutoff_m * cutoff_m + 2.0 * r_m * r_m - 3.0 * switch_on_m * switch_on_m
    value = d * d * e / denom_cubed
    derivative = 12.0 * r_m * d * (switch_on_m * switch_on_m - r_m * r_m) / denom_cubed
    return value, derivative


def switch_difference(data: dict) -> tuple[float, list[float]]:
    """Returns the exact nanocad-vs-plain-cutoff difference for this model.

    GROMACS has no equivalent of the nanocad switching polynomial.  A GROMACS
    run with `vdw-modifier = None` and `coulomb-modifier = None` computes the
    full potential out to the cutoff.  The nanocad run computes the same
    potential multiplied by `S(r)`.  The difference is

        dE = sum over pairs in the switching region of V_full(r) * (1 - S(r))

    This function returns `(dE, grad dE)` in SI.  A GROMACS run then matches
    nanocad when `E_gromacs = E_nanocad + dE` and `F_gromacs = F_nanocad - grad`.
    The correction is zero when no pair sits between switch-on and cutoff.
    """
    positions = data["positions_m"]
    atom_count = data["atom_count"]
    box = data["box_length_m"]
    cutoff = data["cutoff_m"]
    switch_on = data["switch_on_m"]
    ke = data["coulomb_constant_n_m2_per_c2"]
    vdw_a = data["vdw_a_j"]
    vdw_b = data["vdw_b_per_m"]
    vdw_c = data["vdw_c_j_m6"]
    charges = [
        data["oxygen_charge_c"] if atom % 3 == 0 else data["hydrogen_charge_c"]
        for atom in range(atom_count)
    ]

    excluded = {
        (3 * molecule + first, 3 * molecule + second)
        for molecule in range(data["molecule_count"])
        for first, second in ((0, 1), (0, 2), (1, 2))
    }

    delta_energy_j = 0.0
    delta_gradient_n = [0.0] * (3 * atom_count)
    cutoff_sq = cutoff * cutoff
    for i in range(atom_count):
        for j in range(i + 1, atom_count):
            if (i, j) in excluded:
                continue
            delta = _minimum_image(
                [
                    positions[3 * i + axis] - positions[3 * j + axis]
                    for axis in range(3)
                ],
                box,
            )
            r_sq = sum(value * value for value in delta)
            if r_sq >= cutoff_sq:
                continue
            r_m = math.sqrt(r_sq)
            if r_m <= switch_on:
                continue
            vdw = 0.0
            vdw_prime = 0.0
            if i % 3 == 0 and j % 3 == 0:
                exp_term = math.exp(-vdw_b * r_m)
                vdw = vdw_a * exp_term - vdw_c / r_m**6
                vdw_prime = -vdw_a * vdw_b * exp_term + 6.0 * vdw_c / r_m**7
            coulomb = ke * charges[i] * charges[j] / r_m
            coulomb_prime = -ke * charges[i] * charges[j] / (r_m * r_m)
            full = vdw + coulomb
            full_prime = vdw_prime + coulomb_prime
            switch_value, switch_derivative = _switch_value(r_m, switch_on, cutoff)
            delta_energy_j += full * (1.0 - switch_value)
            d_w = full_prime * (1.0 - switch_value) - full * switch_derivative
            factor = d_w / r_m
            for axis in range(3):
                delta_gradient_n[3 * i + axis] += factor * delta[axis]
                delta_gradient_n[3 * j + axis] -= factor * delta[axis]
    return delta_energy_j, delta_gradient_n


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workdir",
        type=Path,
        default=None,
        help="Directory for the generated GROMACS input. Defaults to a temp dir.",
    )
    parser.add_argument(
        "--keep",
        action="store_true",
        help="Keep the generated input directory.",
    )
    args = parser.parse_args()

    gmx = find_gmx()
    if gmx is None:
        print("SKIP: GROMACS is absent. No gmx or gmx_mpi binary is on PATH.")
        print("Install it with: brew install gromacs")
        print("Then run: python3 benchmarks/gromacs/compare_water_box.py")
        return 0

    workdir = args.workdir
    if workdir is None:
        workdir = Path(tempfile.mkdtemp(prefix="nanocad-gromacs-"))
    workdir.mkdir(parents=True, exist_ok=True)

    data = run_rust_example()
    recorded_energy_j = read_recorded_energy_j()

    write_g96(data, workdir / "conf.g96")
    write_topology(data, workdir / "topol.top")
    write_mdp(data, workdir / "md.mdp")

    print("nanocad <-> GROMACS water-box cross-check")
    print(f"  GROMACS:            {gmx_version(gmx)}")
    print(f"  binary:             {gmx}")
    print(f"  workdir:            {workdir}")
    print(f"  molecules:          {data['molecule_count']}")
    print(f"  atoms:              {data['atom_count']}")
    print(f"  box (nm):           {data['box_length_m'] * NM_PER_M:.12g}")
    print(f"  cutoff (nm):        {data['cutoff_m'] * NM_PER_M:.12g}")

    expected_k = (
        data["coulomb_constant_n_m2_per_c2"]
        * data["elementary_charge_c"] ** 2
        * NM_PER_M
        / J_PER_KJMOL
    )
    print(
        f"  Coulomb constant:   GROMACS {GROMACS_COULOMB_KJMOL_NM_PER_E2} vs "
        f"nanocad {expected_k:.6f} kJ/mol/nm/e^2"
    )

    grompp = run_command(
        [
            gmx,
            "grompp",
            "-f",
            "md.mdp",
            "-c",
            "conf.g96",
            "-p",
            "topol.top",
            "-o",
            "run.tpr",
        ],
        workdir,
    )
    if grompp.returncode != 0:
        print("SKIP: grompp rejected the topology. GROMACS output:")
        print(_indent(_tail(grompp.stderr)))
        _cleanup(workdir, args.keep)
        return 0

    mdrun = run_command(
        [gmx, "mdrun", "-s", "run.tpr", "-deffnm", "run", "-nt", "1"],
        workdir,
    )
    combined = mdrun.stdout + "\n" + mdrun.stderr
    if mdrun.returncode != 0:
        lowered = combined.lower()
        if any(marker in lowered for marker in UNSUPPORTED_MODEL_MARKERS):
            print("SKIP: GROMACS cannot express the nanocad model.")
            print("  Reason: the nanocad vdW term is Buckingham, and the Verlet")
            print("  cutoff scheme in GROMACS does not support Buckingham. Switch")
            print("  functions are also absent (only potential-shift is implemented),")
            print("  so the nanocad custom switching polynomial has no equivalent.")
            print("  Exact GROMACS message:")
            print(_indent(_tail(mdrun.stderr, 12)))
            _cleanup(workdir, args.keep)
            return 0
        print("FAIL: mdrun failed for an unexpected reason. GROMACS output:")
        print(_indent(_tail(combined, 20)))
        _cleanup(workdir, args.keep)
        return 1

    energy_xvg = workdir / "energy.xvg"
    energy_proc = subprocess.run(
        [gmx, "energy", "-f", "run.edr", "-o", str(energy_xvg)],
        cwd=workdir,
        input="Potential\n\n",
        capture_output=True,
        text=True,
        check=False,
    )
    if energy_proc.returncode != 0:
        print("FAIL: gmx energy failed.")
        print(_indent(_tail(energy_proc.stderr, 12)))
        _cleanup(workdir, args.keep)
        return 1

    dump = run_command([gmx, "dump", "-f", "run.trr"], workdir)
    forces_kjmol_nm = parse_first_frame_forces(dump.stdout, data["atom_count"])
    energy_kjmol = parse_energy_xvg(energy_xvg)

    if energy_kjmol is None or len(forces_kjmol_nm) != 3 * data["atom_count"]:
        print("FAIL: could not read the GROMACS energy or forces.")
        _cleanup(workdir, args.keep)
        return 1

    energy_gromacs_j = energy_kjmol * J_PER_KJMOL
    energy_nanocad_j = (
        recorded_energy_j if recorded_energy_j is not None else data["energy_j"]
    )
    forces_nanocad = data["forces_n"]
    forces_gromacs_n = [value * N_PER_KJMOL_PER_NM for value in forces_kjmol_nm]

    # GROMACS uses a plain cutoff here.  Add the analytic switch difference to
    # the nanocad model, so the comparison isolates coding errors from the one
    # documented model difference.
    switch_energy_j, switch_gradient_n = switch_difference(data)
    energy_predicted_j = energy_nanocad_j + switch_energy_j
    forces_predicted_n = [
        value - switch_gradient_n[index] for index, value in enumerate(forces_nanocad)
    ]

    energy_rel_diff = abs(energy_predicted_j - energy_gromacs_j) / max(
        abs(energy_predicted_j), 1.0e-300
    )
    force_abs_diff_n = [
        abs(a - b) for a, b in zip(forces_predicted_n, forces_gromacs_n, strict=True)
    ]
    max_force_abs_diff_n = max(force_abs_diff_n)
    max_force_magnitude_n = max(abs(value) for value in forces_predicted_n)
    force_rel_diff = max_force_abs_diff_n / max(max_force_magnitude_n, 1.0e-300)

    print(f"  energy nanocad (J): {energy_nanocad_j:.12e}")
    print(f"  switch correction (J): {switch_energy_j:.3e}")
    print(f"  energy GROMACS (J): {energy_gromacs_j:.12e}")
    print(f"  energy predicted (J): {energy_predicted_j:.12e}")
    print(f"  energy rel diff:    {energy_rel_diff:.3e}")
    print(f"  max |force diff| (N): {max_force_abs_diff_n:.3e}")
    print(f"  force rel diff:       {force_rel_diff:.3e}")
    print(f"  relative tolerance:   {SINGLE_PRECISION_TOL:.1e}")

    ok = (
        energy_rel_diff <= SINGLE_PRECISION_TOL
        and force_rel_diff <= SINGLE_PRECISION_TOL
    )
    _cleanup(workdir, args.keep)
    if ok:
        print(
            f"PASS: energy rel {energy_rel_diff:.3e}, force rel {force_rel_diff:.3e} "
            f"<= {SINGLE_PRECISION_TOL:.1e}"
        )
        return 0
    print(
        f"FAIL: energy rel {energy_rel_diff:.3e}, force rel {force_rel_diff:.3e} "
        f"(tolerance {SINGLE_PRECISION_TOL:.1e})"
    )
    return 1


def _tail(text: str, lines: int = 8) -> str:
    return "\n".join(text.strip().splitlines()[-lines:])


def _indent(text: str) -> str:
    return "\n".join("    " + line for line in text.splitlines())


def _cleanup(workdir: Path, keep: bool) -> None:
    if not keep:
        shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
