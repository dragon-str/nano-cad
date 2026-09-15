"""M7-04: the planetary gearbox built through the tool surface alone.

This test calls :func:`nanocad_mcp.tools.call_tool`, the same path the MCP
server uses. It does not import the low-level core. It is the M7 exit
criterion: an agent completes the gearbox demo through tools alone.
"""

from __future__ import annotations

import json

from nanocad_mcp.session import Session
from nanocad_mcp.tools import call_tool

_GEAR_RATIO_TOLERANCE = 1.0e-3


def _call(session: Session, name: str, arguments: dict) -> dict:
    result = call_tool(name, arguments, session)
    assert result["isError"] is False, result["structuredContent"]
    return result["structuredContent"]


def test_agent_builds_the_planetary_gearbox_through_tools(tmp_path) -> None:
    session = Session()

    # 1. Generate the gear set by name and specs.
    planetary = _call(
        session, "generate_gear", {"name": "planetary", "specs": {"planet_count": 3}}
    )
    assert planetary["atom_count"] > 100
    assert planetary["bond_count"] > 0

    # 2. Assemble the L2 device from the parameter map.
    assembly = _call(session, "assemble_planetary", {"specs": {"planet_count": 3}})
    assert assembly["gear_ratio"] == 3.5
    assert assembly["sun_teeth"] == 24
    assert assembly["ring_teeth"] == 60
    assert assembly["planet_count"] == 3

    # 3. Drive the device and measure the sun-to-carrier ratio.
    ratio = _call(
        session, "measure_gear_ratio", {"assembly_id": assembly["assembly_id"]}
    )
    assert ratio["within_tolerance"] is True
    assert (
        abs(ratio["measured_gear_ratio"] - ratio["analytic_gear_ratio"])
        / ratio["analytic_gear_ratio"]
        < _GEAR_RATIO_TOLERANCE
    )
    assert ratio["max_gear_constraint_error_m"] < 1.0e-9

    # 4. Relax a part.
    lattice = _call(session, "generate_lattice", {"kind": "diamond", "specs": {}})
    relaxed = _call(session, "relax", {"part_id": lattice["part_id"]})
    assert relaxed["energy_j"] >= 0.0
    assert relaxed["converged"] is True

    # 5. Assemble the relaxed part into a document.
    document = _call(
        session,
        "assemble",
        {"part_ids": [lattice["part_id"]], "name": "gearbox-demo"},
    )
    assert document["part_count"] == 1

    # 6. Write and read the NCZ.
    path = tmp_path / "gearbox.ncz"
    saved = _call(
        session, "save", {"path": str(path), "document_id": document["document_id"]}
    )
    assert saved["part_count"] == 1
    loaded = _call(session, "load", {"path": str(path)})
    assert loaded["part_count"] == document["part_count"]
    assert loaded["parts"][0]["name"] == document["part_names"][0]

    # 7. Export URDF and check the structure.
    urdf = _call(session, "export_urdf", {"assembly_id": assembly["assembly_id"]})
    assert urdf["urdf"].startswith("<?xml")
    assert '<robot name="planetary_gearbox"' in urdf["urdf"]
    assert 'link="sun"' in urdf["urdf"]
    assert urdf["link_count"] == 7
    assert urdf["joint_count"] == 6

    # 8. The tool results are structured JSON.
    assert json.loads(json.dumps(ratio)) == ratio

    print(
        "gearbox demo: analytic_ratio {analytic:.9f} measured_ratio {measured:.9f} "
        "relative_error {error:.3e} tolerance {tolerance:.1e}".format(
            analytic=ratio["analytic_gear_ratio"],
            measured=ratio["measured_gear_ratio"],
            error=ratio["relative_error"],
            tolerance=ratio["tolerance"],
        )
    )
