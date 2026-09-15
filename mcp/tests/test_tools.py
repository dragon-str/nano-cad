"""Tests for the nanocad MCP tool registry.

These call the tool layer in process. The JSON-RPC transport is tested
separately in ``test_server.py``.
"""

from __future__ import annotations

import json
from pathlib import Path

from nanocad_mcp.session import Session
from nanocad_mcp.tools import call_tool, list_tools, tool_names

_TOOLS_JSON = Path(__file__).resolve().parents[1] / "tools.json"

_EXPECTED_TOOLS = {
    "create_part",
    "generate_gear",
    "generate_nanotube",
    "generate_lattice",
    "list_generators",
    "relax",
    "run_md",
    "add_jig",
    "measure",
    "validate_part",
    "convert_units",
    "save",
    "load",
    "list",
    "assemble",
    "assemble_planetary",
    "measure_gear_ratio",
    "export_urdf",
}


def test_schema_file_matches_the_registry() -> None:
    """The committed schema equals the live ``tools/list`` output."""
    with _TOOLS_JSON.open(encoding="utf-8") as handle:
        schema = json.load(handle)
    assert schema["schema"] == "nanocad.mcp.tools"
    assert schema["tools"] == list_tools()["tools"]
    assert schema["tools"]
    for entry in schema["tools"]:
        assert entry["inputSchema"]["type"] == "object"


def test_every_tool_has_a_name_description_and_schema() -> None:
    assert set(tool_names()) == _EXPECTED_TOOLS
    for entry in list_tools()["tools"]:
        assert entry["name"] == entry["name"].strip()
        assert entry["description"]
        assert isinstance(entry["inputSchema"]["properties"], dict)


def test_unknown_tool_returns_an_error_object() -> None:
    result = call_tool("does_not_exist", {}, Session())
    assert result["isError"] is True
    assert "does_not_exist" in result["structuredContent"]["error"]


def test_missing_argument_returns_an_error_object() -> None:
    result = call_tool("measure", {}, Session())
    assert result["isError"] is True
    assert result["structuredContent"]["error_type"] == "ValueError"


def test_generate_measure_and_validate_a_part() -> None:
    session = Session()
    created = call_tool("create_part", {"name": "diamond", "specs": {}}, session)
    assert created["isError"] is False
    part_id = created["structuredContent"]["part_id"]
    assert created["structuredContent"]["atom_count"] > 0

    measured = call_tool("measure", {"part_id": part_id}, session)
    assert measured["isError"] is False
    assert (
        measured["structuredContent"]["atom_count"]
        == created["structuredContent"]["atom_count"]
    )
    assert "bounding_box_size_m" in measured["structuredContent"]

    validated = call_tool("validate_part", {"part_id": part_id}, session)
    assert validated["isError"] is False
    assert isinstance(validated["structuredContent"]["violations"], list)


def test_add_anchor_jig_reports_units() -> None:
    session = Session()
    part_id = call_tool("generate_lattice", {"kind": "diamond"}, session)[
        "structuredContent"
    ]["part_id"]
    result = call_tool(
        "add_jig",
        {"kind": "anchor", "part_id": part_id, "k_n_per_m": 1000.0, "atoms": [0, 1]},
        session,
    )
    assert result["isError"] is False
    assert result["structuredContent"]["count"] == 2
    assert "energy_j" in result["structuredContent"]


def test_relax_lowers_the_energy() -> None:
    session = Session()
    # A stretched two-atom bond. The minimizer relaxes toward the rest length.
    created = call_tool("create_part", {"name": "spur_gear", "specs": {}}, session)
    part_id = created["structuredContent"]["part_id"]
    before = call_tool("measure", {"part_id": part_id}, session)["structuredContent"]
    relaxed = call_tool(
        "relax",
        {"part_id": part_id},
        session,
    )
    assert relaxed["isError"] is False
    assert relaxed["structuredContent"]["energy_j"] >= 0.0
    after = call_tool("measure", {"part_id": part_id}, session)["structuredContent"]
    assert after["atom_count"] == before["atom_count"]


def test_convert_units_is_labeled() -> None:
    result = call_tool(
        "convert_units", {"value": 1.0, "from_unit": "nm", "to_unit": "m"}
    )
    assert result["isError"] is False
    payload = result["structuredContent"]
    assert payload["from_unit"] == "nm"
    assert payload["to_unit"] == "m"
    assert payload["result"] == 1.0e-9


def test_save_and_load_round_trip(tmp_path) -> None:
    session = Session()
    part_id = call_tool("create_part", {"name": "diamond"}, session)[
        "structuredContent"
    ]["part_id"]
    document_id = call_tool("assemble", {"part_ids": [part_id]}, session)[
        "structuredContent"
    ]["document_id"]
    path = tmp_path / "design.ncz"
    saved = call_tool("save", {"path": str(path), "document_id": document_id}, session)
    assert saved["isError"] is False
    assert path.exists()

    loaded = call_tool("load", {"path": str(path)}, session)
    assert loaded["isError"] is False
    assert loaded["structuredContent"]["part_count"] == 1

    listed = call_tool(
        "list", {"document_id": loaded["structuredContent"]["document_id"]}, session
    )
    assert listed["isError"] is False
    assert listed["structuredContent"]["count"] == 1
