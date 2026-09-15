#!/usr/bin/env python3
"""Cross-check the nanocad periodic SPC-like water box against OpenMM.

The script runs the Rust example `crates/engine/examples/water_box.rs`, reads the
fixed configuration and parameters from its JSON output, rebuilds the identical
physical model in OpenMM with custom forces, and compares the potential energy
and the per-atom forces.

Every value crosses the boundary in SI. OpenMM works in nanometre, kilojoule per
mole and elementary charge, so the script converts at the edges:

    1 nm      = 1e-9 m
    1 kJ/mol  = 1000 / N_A J
    1 kJ/mol/nm = (1000 / N_A) / 1e-9 N

Run it with an interpreter that has OpenMM:

    /tmp/openmm-venv/bin/python benchmarks/openmm/compare_water_box.py

The script exits zero when OpenMM is absent, so it is never a hard gate.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

# Relative tolerance for the cross-check. Reference runs in double precision and
# reaches about 1e-14. The CPU and OpenCL platforms use single precision in
# places, so they get a looser bound. Either bound catches a real model mismatch
# (a missing switching function on the dispersion term is about 5e-4).
REFERENCE_TOL = 1.0e-6
OTHER_PLATFORM_TOL = 1.0e-5

N_A = 6.02214076e23
J_PER_KJMOL = 1000.0 / N_A
N_PER_KJMOL_PER_NM = J_PER_KJMOL / 1.0e-9
NM_PER_M = 1.0e9


def run_rust_example() -> dict:
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
        raise RuntimeError(
            "the Rust example failed:\n" + proc.stderr.strip()
        )
    return json.loads(proc.stdout)


def build_openmm_model(data: dict):
    import openmm

    box_nm = data["box_length_m"] * NM_PER_M
    cutoff_nm = data["cutoff_m"] * NM_PER_M
    switch_on_nm = data["switch_on_m"] * NM_PER_M
    molecule_count = data["molecule_count"]
    atom_count = data["atom_count"]

    bond_k = data["bond_k_n_per_m"] * 1.0e-18 * (N_A / 1000.0)
    bond_r0_nm = data["bond_r0_m"] * NM_PER_M
    angle_k = data["angle_k_j_per_rad2"] * (N_A / 1000.0)
    angle_theta0 = data["angle_theta0_rad"]

    vdw_a = data["vdw_a_j"] * (N_A / 1000.0)
    vdw_b = data["vdw_b_per_m"] / NM_PER_M
    vdw_c = data["vdw_c_j_m6"] * (NM_PER_M**6) * (N_A / 1000.0)

    q_oxygen = data["oxygen_charge_c"] / data["elementary_charge_c"]
    q_hydrogen = data["hydrogen_charge_c"] / data["elementary_charge_c"]
    coulomb = (
        data["coulomb_constant_n_m2_per_c2"]
        * data["elementary_charge_c"] ** 2
        * NM_PER_M
        / J_PER_KJMOL
    )

    system = openmm.System()
    for molecule in range(molecule_count):
        oxygen = 3 * molecule
        system.addParticle(data["oxygen_mass_kg"])
        system.addParticle(data["hydrogen_mass_kg"])
        system.addParticle(data["hydrogen_mass_kg"])

    bond_force = openmm.CustomBondForce("0.5*k*(r-r0)^2")
    bond_force.addPerBondParameter("k")
    bond_force.addPerBondParameter("r0")
    angle_force = openmm.CustomAngleForce("0.5*k*(theta-theta0)^2")
    angle_force.addPerAngleParameter("k")
    angle_force.addPerAngleParameter("theta0")
    for molecule in range(molecule_count):
        oxygen = 3 * molecule
        hydrogen_one = oxygen + 1
        hydrogen_two = oxygen + 2
        bond_force.addBond(oxygen, hydrogen_one, [bond_k, bond_r0_nm])
        bond_force.addBond(oxygen, hydrogen_two, [bond_k, bond_r0_nm])
        angle_force.addAngle(hydrogen_one, oxygen, hydrogen_two, [angle_k, angle_theta0])

    switching = (
        "select(step(ron-r), 1, "
        "((rc*rc-r*r)*(rc*rc-r*r)*(rc*rc+2*r*r-3*ron*ron))"
        "/((rc*rc-ron*ron)*(rc*rc-ron*ron)*(rc*rc-ron*ron)))"
    )
    vdw_force = openmm.CustomNonbondedForce(
        f"(sqrt(A1*A2)*exp(-0.5*(B1+B2)*r) - sqrt(C1*C2)/r^6)*({switching})"
    )
    vdw_force.addPerParticleParameter("A")
    vdw_force.addPerParticleParameter("B")
    vdw_force.addPerParticleParameter("C")
    vdw_force.addGlobalParameter("rc", cutoff_nm)
    vdw_force.addGlobalParameter("ron", switch_on_nm)

    electrostatic_force = openmm.CustomNonbondedForce(
        f"K*q1*q2*({switching})/r"
    )
    electrostatic_force.addPerParticleParameter("q")
    electrostatic_force.addGlobalParameter("K", coulomb)
    electrostatic_force.addGlobalParameter("rc", cutoff_nm)
    electrostatic_force.addGlobalParameter("ron", switch_on_nm)

    for atom in range(atom_count):
        if atom % 3 == 0:
            vdw_force.addParticle([vdw_a, vdw_b, vdw_c])
            electrostatic_force.addParticle([q_oxygen])
        else:
            vdw_force.addParticle([0.0, 0.0, 0.0])
            electrostatic_force.addParticle([q_hydrogen])

    for molecule in range(molecule_count):
        oxygen = 3 * molecule
        pairs = [(oxygen, oxygen + 1), (oxygen, oxygen + 2), (oxygen + 1, oxygen + 2)]
        for first, second in pairs:
            vdw_force.addExclusion(first, second)
            electrostatic_force.addExclusion(first, second)

    for force in (vdw_force, electrostatic_force):
        force.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
        force.setCutoffDistance(cutoff_nm)
        force.setUseSwitchingFunction(False)
        force.setUseLongRangeCorrection(False)

    system.addForce(bond_force)
    system.addForce(angle_force)
    system.addForce(vdw_force)
    system.addForce(electrostatic_force)
    system.setDefaultPeriodicBoxVectors(
        openmm.Vec3(box_nm, 0.0, 0.0),
        openmm.Vec3(0.0, box_nm, 0.0),
        openmm.Vec3(0.0, 0.0, box_nm),
    )
    return system


def compute_openmm(system, data: dict, platform_name: str):
    import numpy as np
    import openmm
    from openmm import unit

    positions_nm = [
        openmm.Vec3(
            data["positions_m"][3 * atom] * NM_PER_M,
            data["positions_m"][3 * atom + 1] * NM_PER_M,
            data["positions_m"][3 * atom + 2] * NM_PER_M,
        )
        * unit.nanometer
        for atom in range(data["atom_count"])
    ]
    integrator = openmm.VerletIntegrator(0.001)
    platform = openmm.Platform.getPlatformByName(platform_name)
    context = openmm.Context(system, integrator, platform)
    context.setPositions(positions_nm)
    state = context.getState(getEnergy=True, getForces=True)
    energy_kjmol = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    forces = state.getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer
    )
    forces_n = np.asarray(forces, dtype=float) * N_PER_KJMOL_PER_NM
    return energy_kjmol, forces_n


def main() -> int:
    try:
        import openmm  # noqa: F401
    except ImportError:
        print("SKIP: OpenMM is not importable in this interpreter.")
        print("Install it with: uv venv --python 3.12 /tmp/openmm-venv")
        print("                 uv pip install --python /tmp/openmm-venv/bin/python openmm numpy")
        print("Then run: /tmp/openmm-venv/bin/python benchmarks/openmm/compare_water_box.py")
        return 0

    import numpy as np
    import openmm

    platform_name = os.environ.get("NANOCAD_OPENMM_PLATFORM", "Reference")
    tolerance = REFERENCE_TOL if platform_name == "Reference" else OTHER_PLATFORM_TOL
    data = run_rust_example()
    system = build_openmm_model(data)
    energy_kjmol, forces_n = compute_openmm(system, data, platform_name)

    energy_openmm_j = energy_kjmol * J_PER_KJMOL
    energy_nanocad_j = data["energy_j"]
    energy_abs_diff_j = abs(energy_nanocad_j - energy_openmm_j)
    energy_rel_diff = energy_abs_diff_j / max(abs(energy_openmm_j), 1.0e-300)

    forces_nanocad = np.asarray(data["forces_n"], dtype=float).reshape(-1, 3)
    force_component_abs_diff = np.abs(forces_nanocad - forces_n)
    per_atom_diff = np.linalg.norm(forces_nanocad - forces_n, axis=1)
    max_force_ref = float(np.max(np.linalg.norm(forces_n, axis=1)))
    max_force_abs_diff = float(np.max(per_atom_diff))
    force_rel_diff = max_force_abs_diff / max(max_force_ref, 1.0e-300)
    component_abs_diff = float(np.max(force_component_abs_diff))

    print("nanocad <-> OpenMM water-box cross-check")
    print(f"  host platform:      {platform_name}")
    print(f"  OpenMM version:     {openmm.version.version}")
    print(f"  Python:             {sys.version.split()[0]}")
    print(f"  molecules:          {data['molecule_count']}")
    print(f"  atoms:              {data['atom_count']}")
    print(f"  box (nm):           {data['box_length_m'] * NM_PER_M:.12g}")
    print(f"  cutoff (nm):        {data['cutoff_m'] * NM_PER_M:.12g}")
    print(f"  energy nanocad (J): {energy_nanocad_j:.12e}")
    print(f"  energy OpenMM (J):  {energy_openmm_j:.12e}")
    print(f"  energy abs diff (J):{energy_abs_diff_j:.3e}  rel {energy_rel_diff:.3e}")
    print(f"  max |force diff| (N): {max_force_abs_diff:.3e}")
    print(f"  max component diff (N): {component_abs_diff:.3e}")
    print(f"  max force magnitude (N): {max_force_ref:.3e}")
    print(f"  force rel diff:       {force_rel_diff:.3e}")
    print(f"  relative tolerance:   {tolerance:.1e}")

    ok = energy_rel_diff <= tolerance and force_rel_diff <= tolerance
    if ok:
        print(
            f"PASS: energy rel {energy_rel_diff:.3e} <= {tolerance:.1e}, "
            f"force rel {force_rel_diff:.3e} <= {tolerance:.1e}"
        )
        return 0
    print(
        f"FAIL: energy rel {energy_rel_diff:.3e} (tol {tolerance:.1e}), "
        f"force rel {force_rel_diff:.3e} (tol {tolerance:.1e})"
    )
    return 1


if __name__ == "__main__":
    sys.exit(main())
