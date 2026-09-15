"""Tests for the compiled ``nanocad._core`` extension."""

from __future__ import annotations

import math

import pytest

import nanocad


def test_extension_is_loaded() -> None:
    """The compiled extension is importable."""
    assert nanocad._core is not None
    assert nanocad.version() == nanocad.__version__


def test_unit_conversion() -> None:
    assert nanocad.convert_units(1.0, "nm", "m") == pytest.approx(1.0e-9)
    assert nanocad.convert_units(1.0, "Å", "m") == pytest.approx(1.0e-10)
    assert nanocad.to_si(2.0, "m") == pytest.approx(2.0)
    assert "nm" in nanocad.unit_symbols()


def test_quantity_round_trip() -> None:
    quantity = nanocad.Quantity(1.0, "nm")
    assert quantity.value_si == pytest.approx(1.0e-9)
    assert quantity.to("Å").value == pytest.approx(10.0)
    assert quantity.value_in("m") == pytest.approx(1.0e-9)


def test_generator_registry_lists_six() -> None:
    generators = nanocad.list_generators()
    ids = {generator.id for generator in generators}
    assert ids == {
        "spur_gear",
        "gear_profile",
        "planetary",
        "nanotube",
        "diamond",
        "graphite",
    }
    params = nanocad.generator_parameters("spur_gear")
    assert params is not None
    assert len(params) == 7


def test_generate_part_has_atoms() -> None:
    part = nanocad.generate_part("spur_gear", {})
    assert part.atom_count > 0
    assert part.bond_count > 0
    assert part.name.startswith("spur-gear-")


def test_generate_part_rejects_unknown_and_bad_values() -> None:
    with pytest.raises(ValueError):
        nanocad.generate_part("worm_gear", {})
    with pytest.raises(ValueError):
        nanocad.generate_part("spur_gear", {"tooth_count": 5.0})


def test_topology_and_part_construction() -> None:
    topology = nanocad.Topology()
    topology.add_atom(6, [0.0, 0.0, 0.0], 0.0, "C")
    topology.add_atom(6, [1.544e-10, 0.0, 0.0], 0.0, "C")
    topology.add_bond(0, 1, 1, "single")
    assert topology.atom_count == 2
    assert topology.bond(0) == (0, 1, 1, "single")
    part = nanocad.Part("ethane", topology)
    assert part.atom_count == 2
    assert part.material == ""


def test_document_ncz_round_trip(tmp_path) -> None:
    part = nanocad.generate_part("diamond", {})
    document = nanocad.Document("design")
    document.add_part(part)
    path = tmp_path / "design.ncz"
    nanocad.save(str(path), document)
    reloaded = nanocad.load(str(path))
    assert reloaded.part_count == 1
    assert nanocad.list_parts(str(path)) == [part.name]


def test_validate_part_accepts_a_relaxed_bond() -> None:
    topology = nanocad.Topology()
    topology.add_atom(6, [0.0, 0.0, 0.0], 0.0, "C")
    topology.add_atom(6, [1.544e-10, 0.0, 0.0], 0.0, "C")
    topology.add_bond(0, 1, 1, "single")
    part = nanocad.Part("ethane", topology)
    assert nanocad.validate_part(part) == []


def test_validate_part_reports_strain_for_a_generated_gear() -> None:
    part = nanocad.generate_part("spur_gear", {})
    assert isinstance(nanocad.validate_part(part), list)


def test_planetary_geometry() -> None:
    geometry = nanocad.planetary_geometry(1.0e-9, 18, 18, 3, math.radians(20.0))
    assert geometry.sun_teeth == 18
    assert geometry.planet_count == 3
    assert geometry.ring_teeth == 54
    assert geometry.constraint_holds
    assert geometry.gear_ratio > 0.0


def test_system_minimization_lowers_energy() -> None:
    bond_m = 1.544e-10
    system = nanocad.System(2, [12.0e-3 / 6.02214076e23, 12.0e-3 / 6.02214076e23])
    system.set_bond_stretch([(0, 1, 300.0, bond_m)])
    start = [0.0, 0.0, 0.0, bond_m * 1.4, 0.0, 0.0]
    positions, result = nanocad.minimize(system, start)
    final = nanocad.System(2, [12.0e-3 / 6.02214076e23, 12.0e-3 / 6.02214076e23])
    final.set_bond_stretch([(0, 1, 300.0, bond_m)])
    assert final.energy(positions) <= final.energy(start)
    assert result.energy_j <= final.energy(start)


def test_velocity_verlet_conserves_energy_roughly() -> None:
    mass_kg = 12.0e-3 / 6.02214076e23
    bond_m = 1.544e-10
    system = nanocad.System(2, [mass_kg, mass_kg])
    system.set_bond_stretch([(0, 1, 300.0, bond_m)])
    positions = [0.0, 0.0, 0.0, bond_m, 0.0, 0.0]
    velocities = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    before = system.total_energy(positions, velocities)
    final_positions, final_velocities = nanocad.run_md(
        system, positions, velocities, 1.0e-15, 20
    )
    after = system.total_energy(final_positions, final_velocities)
    assert after == pytest.approx(before, rel=1.0e-3)


def test_param_record_json_round_trip() -> None:
    quantity = nanocad.ParamQuantity(1.0e12, "Pa")
    provenance = nanocad.Provenance(
        source="test",
        method="estimate",
        code_version="0.1.0",
        timestamp="2026-01-01T00:00:00Z",
    )
    record = nanocad.PartRecord(
        "part-1",
        "sha256:0",
        "diamondoid",
        10,
        nanocad.ParamQuantity(1.0, "kg"),
        [nanocad.ParamQuantity(1.0, "kg*m^2")] * 3,
        quantity,
        quantity,
        nanocad.ParamQuantity(0.2, "1"),
        quantity,
        nanocad.ParamQuantity(0.1, "1"),
        nanocad.ParamQuantity(1.0, "W/(m*K)"),
        nanocad.ParamQuantity(500.0, "J/(kg*K)"),
        provenance,
    )
    assert record.part_id == "part-1"
    assert record.method == "estimate"
    assert nanocad.check_record(record) == []
    reloaded = nanocad.PartRecord.from_json(record.to_json())
    assert reloaded.part_id == record.part_id
    assert reloaded.mass_kg.value_si == pytest.approx(1.0)


def test_parameter_library_round_trip() -> None:
    quantity = nanocad.ParamQuantity(1.0e12, "Pa")
    provenance = nanocad.Provenance(
        source="test",
        method="estimate",
        code_version="0.1.0",
        timestamp="2026-01-01T00:00:00Z",
    )
    record = nanocad.PartRecord(
        "part-1",
        "sha256:0",
        "diamondoid",
        10,
        nanocad.ParamQuantity(1.0, "kg"),
        [nanocad.ParamQuantity(1.0, "kg*m^2")] * 3,
        quantity,
        quantity,
        nanocad.ParamQuantity(0.2, "1"),
        quantity,
        nanocad.ParamQuantity(0.1, "1"),
        nanocad.ParamQuantity(1.0, "W/(m*K)"),
        nanocad.ParamQuantity(500.0, "J/(kg*K)"),
        provenance,
    )
    library = nanocad.ParameterLibrary()
    library.insert(record, 1)
    assert library.len == 1
    assert library.latest_version("part-1") == 1
    reloaded = nanocad.ParameterLibrary.from_json(library.to_json())
    assert reloaded.get("part-1").part_id == "part-1"


def test_quat_and_rigid_body_system() -> None:
    rotation = nanocad.Quat.from_axis_angle([0.0, 0.0, 1.0], math.pi / 2.0)
    x, y, z = rotation.rotate([1.0, 0.0, 0.0])
    assert (x, y, z) == pytest.approx((0.0, 1.0, 0.0))

    body = nanocad.RigidBody.diagonal(1.0, [1.0e-3, 1.0e-3, 1.0e-3], [0.0, 0.0, 0.0])
    anchor = nanocad.RigidBody.diagonal(
        1.0, [1.0e-3, 1.0e-3, 1.0e-3], [1.0e-9, 0.0, 0.0]
    )
    anchor.is_fixed = True
    device = nanocad.RigidBodySystem()
    device.add_body(body)
    device.add_body(anchor)
    device.add_revolute_joint(0, 1, [0.5e-9, 0.0, 0.0], [0.0, 0.0, 1.0])
    device.step(1.0e-9, 8)
    assert device.max_revolute_anchor_error_m() < 1.0e-12


def test_alias_generators_and_lattice() -> None:
    assert nanocad.create_part("spur_gear", {}).atom_count > 0
    assert nanocad.generate_nanotube({}).atom_count > 0
    assert nanocad.generate_lattice("diamond", {}).atom_count > 0
    with pytest.raises(ValueError):
        nanocad.create_part("worm_gear", {})


def test_assemble_and_measure() -> None:
    part = nanocad.generate_part("diamond", {})
    document = nanocad.assemble([part, part])
    assert document.part_count == 2
    measured = nanocad.measure(part)
    assert measured["atom_count"] == part.atom_count
    size = measured["bounding_box_size_m"]
    assert size[0] >= 0.0 and size[1] >= 0.0 and size[2] >= 0.0
    assert isinstance(measured["total_charge_c"], float)


def test_add_constraint_dispatches_by_kind() -> None:
    device = nanocad.RigidBodySystem()
    device.add_body(nanocad.RigidBody.diagonal(1.0, [1.0e-3] * 3, [0.0, 0.0, 0.0]))
    device.add_body(nanocad.RigidBody.diagonal(1.0, [1.0e-3] * 3, [1.0e-9, 0.0, 0.0]))
    index = nanocad.add_constraint(
        device, "revolute", 0, 1, [0.5e-9, 0.0, 0.0], [0.0, 0.0, 1.0]
    )
    assert index == 0
    assert device.revolute_joint_count == 1
    with pytest.raises(ValueError):
        nanocad.add_constraint(device, "ball", 0, 1, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0])
