"""The nanocad MCP tool registry.

Every tool wraps the public ``nanocad`` Python API. The MCP server and the
end-to-end demo test both call :func:`call_tool`, so there is one call path.

Two rules hold for every tool:

1. The result is structured JSON. Numeric fields carry their unit in the field
   name, for example ``energy_j`` or ``length_m``.
2. No exception crosses the tool boundary. A bad request returns an error
   object with ``isError`` set.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

import nanocad

from .session import Session, carbon_atom_mass_kg

Handler = Callable[[Session, dict[str, Any]], dict[str, Any]]


def _number(description: str, default: float | None = None) -> dict[str, Any]:
    schema: dict[str, Any] = {"type": "number", "description": description}
    if default is not None:
        schema["default"] = default
    return schema


def _integer(description: str, default: int | None = None) -> dict[str, Any]:
    schema: dict[str, Any] = {"type": "integer", "description": description}
    if default is not None:
        schema["default"] = default
    return schema


def _string(description: str, enum: list[str] | None = None) -> dict[str, Any]:
    schema: dict[str, Any] = {"type": "string", "description": description}
    if enum is not None:
        schema["enum"] = enum
    return schema


def _object(description: str) -> dict[str, Any]:
    return {"type": "object", "description": description}


def _schema(
    properties: dict[str, Any],
    required: list[str] | None = None,
) -> dict[str, Any]:
    return {
        "type": "object",
        "properties": properties,
        "required": required or [],
        "additionalProperties": False,
    }


@dataclass(frozen=True)
class Tool:
    """One MCP tool: its name, description, input schema, and handler."""

    name: str
    description: str
    input_schema: dict[str, Any]
    handler: Handler


_TOOLS: list[Tool] = []
_HANDLERS: dict[str, Handler] = {}


def tool(
    name: str, description: str, input_schema: dict[str, Any]
) -> Callable[[Handler], Handler]:
    """Register a tool handler under a name and an input schema."""

    def register(handler: Handler) -> Handler:
        if name in _HANDLERS:
            raise RuntimeError(f"duplicate tool {name!r}")
        _TOOLS.append(Tool(name, description, input_schema, handler))
        _HANDLERS[name] = handler
        return handler

    return register


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _jsonable(value: Any) -> Any:
    """Return a JSON-safe copy. A non-finite float becomes ``None``."""
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, bool) or value is None:
        return value
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    return value


def _part_summary(session: Session, handle: str, part: Any) -> dict[str, Any]:
    return {
        "part_id": handle,
        "name": part.name,
        "material": part.material,
        "atom_count": part.atom_count,
        "bond_count": part.bond_count,
    }


def _part_positions_m(part: Any) -> list[float]:
    return list(part.topology.positions_m())


def _bond_term(part: Any, stiffness_n_per_m: float):
    """Build bond-stretch entries from the current bond lengths."""
    topology = part.topology
    positions_m = topology.positions_m()
    term = []
    for u, v, _order, _bond_type in topology.bonds():
        du = [
            positions_m[3 * u + axis] - positions_m[3 * v + axis] for axis in range(3)
        ]
        rest_m = math.sqrt(sum(value * value for value in du))
        term.append((u, v, stiffness_n_per_m, rest_m))
    return term


def _system_for_part(part: Any, stiffness_n_per_m: float, with_masses: bool):
    atom_count = part.atom_count
    if with_masses:
        mass_kg = carbon_atom_mass_kg()
        system = nanocad.System(atom_count, [mass_kg] * atom_count)
    else:
        system = nanocad.System(atom_count)
    system.set_bond_stretch(_bond_term(part, stiffness_n_per_m))
    return system


def _require(args: dict[str, Any], key: str) -> Any:
    if key not in args:
        raise ValueError(f"missing required argument {key!r}")
    return args[key]


# ---------------------------------------------------------------------------
# Part creation
# ---------------------------------------------------------------------------


@tool(
    "create_part",
    "Create a part from a registered generator id and a parameter map.",
    _schema(
        {
            "name": _string("Generator id, for example 'spur_gear' or 'nanotube'."),
            "specs": _object("Generator parameters. Missing values use defaults."),
            "part_id": _string("Optional handle to reuse for this part."),
        },
        ["name"],
    ),
)
def _create_part(session: Session, args: dict[str, Any]) -> dict[str, Any]:
    specs = dict(args.get("specs") or {})
    part = nanocad.create_part(_require(args, "name"), specs)
    handle = session.store_part(part, args.get("part_id"))
    return _part_summary(session, handle, part)


@tool(
    "generate_gear",
    "Generate a gear from the gear registry: spur_gear, gear_profile, or planetary.",
    _schema(
        {
            "name": _string(
                "Gear generator id.", ["spur_gear", "gear_profile", "planetary"]
            ),
            "specs": _object("Generator parameters. Missing values use defaults."),
            "part_id": _string("Optional handle to reuse for this part."),
        },
        ["name"],
    ),
)
def _generate_gear(session: Session, args: dict[str, Any]) -> dict[str, Any]:
    specs = dict(args.get("specs") or {})
    part = nanocad.generate_gear(_require(args, "name"), specs)
    handle = session.store_part(part, args.get("part_id"))
    return _part_summary(session, handle, part)


@tool(
    "generate_nanotube",
    "Generate a nanotube from the default parameters, or from a parameter map.",
    _schema(
        {
            "specs": _object("Nanotube parameters, for example chiral indices."),
            "part_id": _string("Optional handle to reuse for this part."),
        }
    ),
)
def _generate_nanotube(session: Session, args: dict[str, Any]) -> dict[str, Any]:
    specs = dict(args.get("specs") or {})
    part = nanocad.generate_nanotube(specs)
    handle = session.store_part(part, args.get("part_id"))
    return _part_summary(session, handle, part)


@tool(
    "generate_lattice",
    "Generate a diamond or graphite lattice part.",
    _schema(
        {
            "kind": _string("Lattice kind.", ["diamond", "graphite"]),
            "specs": _object("Lattice parameters. Missing values use defaults."),
            "part_id": _string("Optional handle to reuse for this part."),
        },
        ["kind"],
    ),
)
def _generate_lattice(session: Session, args: dict[str, Any]) -> dict[str, Any]:
    specs = dict(args.get("specs") or {})
    part = nanocad.generate_lattice(_require(args, "kind"), specs)
    handle = session.store_part(part, args.get("part_id"))
    return _part_summary(session, handle, part)


@tool(
    "list_generators",
    "List every registered generator with its parameter schema.",
    _schema({}),
)
def _list_generators(session: Session, args: dict[str, Any]) -> dict[str, Any]:
    generators = []
    for generator in nanocad.list_generators():
        generators.append(
            {
                "id": generator.id,
                "name": generator.name,
                "parameters": [
                    {
                        "name": parameter.name,
                        "unit": parameter.unit,
                        "default_value": parameter.default_value,
                        "min": parameter.min,
                        "max": parameter.max,
                        "integer": parameter.integer,
                        "description": parameter.description,
                    }
                    for parameter in generator.parameters
                ],
            }
        )
    return {"generators": generators, "count": len(generators)}


# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------


@tool(
    "relax",
    "Minimize a part's energy with harmonic bond-stretch terms at current lengths.",
    _schema(
        {
            "part_id": _string("Handle of the part to relax."),
            "stiffness_n_per_m": _number(
                "Bond-stretch stiffness in newtons per metre.", 300.0
            ),
            "max_iterations": _integer("Maximum minimizer iterations.", 200),
            "gradient_tolerance_n": _number(
                "Convergence threshold on the gradient norm, in newtons.", 1.0e-14
            ),
        },
        ["part_id"],
    ),
)
def _relax(session: Session, args: dict[str, Any]) -> dict[str, Any]:
    part = session.require_part(_require(args, "part_id"))
    stiffness = float(args.get("stiffness_n_per_m", 300.0))
    system = _system_for_part(part, stiffness, with_masses=False)
    positions_m = _part_positions_m(part)
    new_positions_m, result = nanocad.relax(
        system,
        positions_m,
        int(args.get("max_iterations", 200)),
        float(args.get("gradient_tolerance_n", 1.0e-14)),
    )
    part.set_positions_m(list(new_positions_m))
    return {
        "part_id": args["part_id"],
        "energy_j": result.energy_j,
        "gradient_norm_n": result.gradient_norm_n,
        "iterations": result.iterations,
        "converged": result.converged,
    }


@tool(
    "run_md",
    "Run Velocity Verlet molecular dynamics on a part and update its positions.",
    _schema(
        {
            "part_id": _string("Handle of the part to integrate."),
            "dt_s": _number("Time step in seconds.", 1.0e-15),
            "steps": _integer("Number of steps.", 100),
            "stiffness_n_per_m": _number(
                "Bond-stretch stiffness in newtons per metre.", 300.0
            ),
        },
        ["part_id"],
    ),
)
def _run_md(session: Session, args: dict[str, Any]) -> dict[str, Any]:
    part = session.require_part(_require(args, "part_id"))
    stiffness = float(args.get("stiffness_n_per_m", 300.0))
    dt_s = float(args.get("dt_s", 1.0e-15))
    steps = int(args.get("steps", 100))
    system = _system_for_part(part, stiffness, with_masses=True)
    positions_m = _part_positions_m(part)
    velocities_m_per_s = [0.0] * len(positions_m)
    new_positions_m, new_velocities = nanocad.run_md(
        system, positions_m, velocities_m_per_s, dt_s, steps
    )
    part.set_positions_m(list(new_positions_m))
    return {
        "part_id": args["part_id"],
        "steps": steps,
        "dt_s": dt_s,
        "temperature_k": system.temperature(list(new_velocities)),
        "kinetic_energy_j": system.kinetic_energy(list(new_velocities)),
    }


@tool(
    "add_jig",
    "Add an anchor or spring jig to a part. The jig is stored under ``jig_id``.",
    _schema(
        {
            "kind": _string("Jig kind.", ["anchor", "spring"]),
            "part_id": _string("Handle of the part the jig acts on."),
            "k_n_per_m": _number(
                "Spring constant in newtons per metre, for an anchor jig.", 1000.0
            ),
            "atoms": {
                "type": "array",
                "items": {"type": "integer"},
                "description": "Atom indices to hold, for an anchor jig.",
            },
            "springs": {
                "type": "array",
                "description": "Spring definitions, for a spring jig.",
                "items": {
                    "type": "object",
                    "properties": {
                        "i": {"type": "integer"},
                        "j": {"type": "integer"},
                        "k_n_per_m": {"type": "number"},
                        "rest_length_m": {"type": "number"},
                    },
                    "required": ["i", "j", "k_n_per_m", "rest_length_m"],
                },
            },
        },
        ["kind", "part_id"],
    ),
)
def _add_jig(session: Session, args: dict[str, Any]) -> dict[str, Any]:
    kind = _require(args, "kind")
    part = session.require_part(_require(args, "part_id"))
    topology = part.topology
    positions_m = _part_positions_m(part)

    if kind == "anchor":
        jig = nanocad.AnchorJig(topology, float(args.get("k_n_per_m", 1000.0)))
        atoms = args.get("atoms") or []
        if not atoms:
            raise ValueError("an anchor jig needs at least one atom index")
        for atom in atoms:
            jig.hold_atom_current(int(atom))
        energy_j = jig.energy(positions_m)
        count = jig.held_count
    elif kind == "spring":
        jig = nanocad.SpringJig(topology)
        springs = args.get("springs") or []
        if not springs:
            raise ValueError("a spring jig needs at least one spring")
        for spring in springs:
            jig.add_spring(
                int(spring["i"]),
                int(spring["j"]),
                float(spring["k_n_per_m"]),
                float(spring["rest_length_m"]),
            )
        energy_j = jig.energy(positions_m)
        count = jig.spring_count
    else:
        raise ValueError(f"unknown jig kind {kind!r}; expected anchor or spring")

    handle = session.store_jig(jig)
    return {
        "jig_id": handle,
        "kind": kind,
        "part_id": args["part_id"],
        "count": count,
        "energy_j": energy_j,
    }


# ---------------------------------------------------------------------------
# Inspection
# ---------------------------------------------------------------------------


@tool(
    "measure",
    "Measure a part: atom and bond counts, bounding box, centroid, and charge.",
    _schema({"part_id": _string("Handle of the part.")}, ["part_id"]),
)
def _measure(session: Session, args: dict[str, Any]) -> dict[str, Any]:
    part = session.require_part(_require(args, "part_id"))
    result = dict(nanocad.measure(part))
    result["part_id"] = args["part_id"]
    return result


@tool(
    "validate_part",
    "Check a part for over-coordination, bond strain, and unknown valence.",
    _schema(
        {
            "part_id": _string("Handle of the part."),
            "reference_m": _number(
                "Reference bond length in metres. Defaults to the diamond C-C bond."
            ),
            "strain_tolerance_relative": _number("Relative strain tolerance.", 0.35),
        },
        ["part_id"],
    ),
)
def _validate_part(session: Session, args: dict[str, Any]) -> dict[str, Any]:
    part = session.require_part(_require(args, "part_id"))
    reference_m = args.get("reference_m")
    tolerance = args.get("strain_tolerance_relative")
    if tolerance is None:
        violations = nanocad.validate_part(part, reference_m)
    else:
        violations = nanocad.validate_part(part, reference_m, float(tolerance))
    return {
        "part_id": args["part_id"],
        "valid": not violations,
        "violation_count": len(violations),
        "violations": list(violations),
    }


@tool(
    "convert_units",
    "Convert a value between two units of the same dimension.",
    _schema(
        {
            "value": _number("The value to convert."),
            "from_unit": _string("Source unit symbol, for example 'nm'."),
            "to_unit": _string("Target unit symbol, for example 'm'."),
        },
        ["value", "from_unit", "to_unit"],
    ),
)
def _convert_units(session: Session, args: dict[str, Any]) -> dict[str, Any]:
    value = float(_require(args, "value"))
    from_unit = _require(args, "from_unit")
    to_unit = _require(args, "to_unit")
    result = nanocad.convert_units(value, from_unit, to_unit)
    return {
        "value": value,
        "from_unit": from_unit,
        "to_unit": to_unit,
        "result": result,
    }


# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------


@tool(
    "save",
    "Save a document. The format follows the path extension: ncz, xyz, or mmp.",
    _schema(
        {
            "path": _string("Output file path."),
            "document_id": _string("Handle of the document to save."),
        },
        ["path", "document_id"],
    ),
)
def _save(session: Session, args: dict[str, Any]) -> dict[str, Any]:
    path = _require(args, "path")
    document = session.require_document(_require(args, "document_id"))
    nanocad.save(path, document)
    return {
        "path": path,
        "document_id": args["document_id"],
        "part_count": document.part_count,
        "part_names": document.part_names(),
    }


@tool(
    "load",
    "Load a document from a file. The format follows the path extension.",
    _schema({"path": _string("Input file path.")}, ["path"]),
)
def _load(session: Session, args: dict[str, Any]) -> dict[str, Any]:
    path = _require(args, "path")
    document = nanocad.load(path)
    handle = session.store_document(document)
    parts = []
    for index in range(document.part_count):
        part = document.part(index)
        part_handle = session.store_part(part)
        parts.append(_part_summary(session, part_handle, part))
    return {
        "document_id": handle,
        "path": path,
        "part_count": document.part_count,
        "parts": parts,
    }


@tool(
    "list",
    "List the part handles in a document, or every part handle in the session.",
    _schema(
        {"document_id": _string("Optional document handle. Omit to list the session.")}
    ),
)
def _list(session: Session, args: dict[str, Any]) -> dict[str, Any]:
    document_id = args.get("document_id")
    if document_id is not None:
        document = session.require_document(document_id)
        parts = []
        for index in range(document.part_count):
            part = document.part(index)
            parts.append(
                {
                    "part_id": None,
                    "name": part.name,
                    "material": part.material,
                    "atom_count": part.atom_count,
                    "bond_count": part.bond_count,
                }
            )
        return {"source": document_id, "count": len(parts), "parts": parts}
    parts = [
        _part_summary(session, handle, part) for handle, part in session.parts.items()
    ]
    return {"source": "session", "count": len(parts), "parts": parts}


# ---------------------------------------------------------------------------
# Assembly and the L2 device
# ---------------------------------------------------------------------------


@tool(
    "assemble",
    "Assemble parts into a document and store it under a document handle.",
    _schema(
        {
            "part_ids": {
                "type": "array",
                "items": {"type": "string"},
                "description": "Part handles to assemble, in order.",
            },
            "name": _string(
                "Document name.",
            ),
        },
        ["part_ids"],
    ),
)
def _assemble(session: Session, args: dict[str, Any]) -> dict[str, Any]:
    part_ids = _require(args, "part_ids")
    if not part_ids:
        raise ValueError("assemble needs at least one part id")
    parts = [session.require_part(part_id) for part_id in part_ids]
    document = nanocad.assemble(parts)
    document.name = args.get("name", "assembly")
    handle = session.store_document(document)
    return {
        "document_id": handle,
        "name": document.name,
        "part_count": document.part_count,
        "part_names": document.part_names(),
    }


@tool(
    "assemble_planetary",
    "Assemble a planetary gearbox at the L2 device level from a parameter map.",
    _schema(
        {
            "specs": _object("Planetary parameters, for example {'planet_count': 3}."),
            "assembly_id": _string("Optional handle to reuse for this assembly."),
        }
    ),
)
def _assemble_planetary(session: Session, args: dict[str, Any]) -> dict[str, Any]:
    specs = dict(args.get("specs") or {})
    assembly = nanocad.assemble_planetary(specs)
    handle = session.store_assembly(assembly, args.get("assembly_id"))
    return {
        "assembly_id": handle,
        "gear_ratio": assembly.gear_ratio,
        "gear_ratio_unit": "1",
        "sun_teeth": assembly.sun_teeth,
        "planet_teeth": assembly.planet_teeth,
        "ring_teeth": assembly.ring_teeth,
        "planet_count": assembly.planet_count,
        "body_count": assembly.body_count,
    }


_GEAR_RATIO_TOLERANCE = 1.0e-3


@tool(
    "measure_gear_ratio",
    "Drive the sun and measure the sun-to-carrier gear ratio of an assembled gearbox.",
    _schema(
        {
            "assembly_id": _string("Handle of an assembled planetary gearbox."),
            "sun_angular_velocity_rad_per_s": _number(
                "Sun angular velocity in radians per second.", 0.5
            ),
            "dt_s": _number("Device time step in seconds.", 1.0e-3),
            "steps": _integer("Number of device steps.", 2000),
            "iterations": _integer("Constraint iterations per step.", 64),
            "tolerance": _number("Relative tolerance for the ratio check.", 1.0e-3),
        },
        ["assembly_id"],
    ),
)
def _measure_gear_ratio(session: Session, args: dict[str, Any]) -> dict[str, Any]:
    assembly = session.require_assembly(_require(args, "assembly_id"))
    dt_s = float(args.get("dt_s", 1.0e-3))
    steps = int(args.get("steps", 2000))
    iterations = int(args.get("iterations", 64))
    tolerance = float(args.get("tolerance", _GEAR_RATIO_TOLERANCE))

    assembly.set_consistent_rotation(
        float(args.get("sun_angular_velocity_rad_per_s", 0.5))
    )
    sun_start_rad = assembly.sun_angle_rad()
    carrier_start_rad = assembly.carrier_angle_rad()
    max_gear_error_m = 0.0
    for _ in range(steps):
        assembly.step(dt_s, iterations)
        max_gear_error_m = max(max_gear_error_m, assembly.max_gear_constraint_error_m())

    sun_delta_rad = assembly.sun_angle_rad() - sun_start_rad
    carrier_delta_rad = assembly.carrier_angle_rad() - carrier_start_rad
    if carrier_delta_rad == 0.0:
        raise ValueError("the carrier did not move; increase the rotation or the steps")
    measured = sun_delta_rad / carrier_delta_rad
    analytic = assembly.gear_ratio
    relative_error = abs(measured - analytic) / analytic
    return {
        "assembly_id": args["assembly_id"],
        "analytic_gear_ratio": analytic,
        "measured_gear_ratio": measured,
        "relative_error": relative_error,
        "tolerance": tolerance,
        "within_tolerance": relative_error < tolerance,
        "sun_delta_rad": sun_delta_rad,
        "carrier_delta_rad": carrier_delta_rad,
        "max_gear_constraint_error_m": max_gear_error_m,
        "steps": steps,
        "dt_s": dt_s,
    }


@tool(
    "export_urdf",
    "Export an assembled gearbox to URDF text, and optionally write it to a file.",
    _schema(
        {
            "assembly_id": _string("Handle of an assembled planetary gearbox."),
            "name": _string(
                "Robot name in the URDF.",
            ),
            "path": _string("Optional output file path."),
        },
        ["assembly_id"],
    ),
)
def _export_urdf(session: Session, args: dict[str, Any]) -> dict[str, Any]:
    assembly = session.require_assembly(_require(args, "assembly_id"))
    name = args.get("name")
    xml = assembly.export_urdf(name)
    path = args.get("path")
    if path:
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(xml)
    return {
        "assembly_id": args["assembly_id"],
        "path": path,
        "link_count": xml.count("<link"),
        "joint_count": xml.count("<joint "),
        "urdf": xml,
    }


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------


def list_tools() -> dict[str, Any]:
    """Return the MCP ``tools/list`` result."""
    return {
        "tools": [
            {
                "name": entry.name,
                "description": entry.description,
                "inputSchema": entry.input_schema,
            }
            for entry in _TOOLS
        ]
    }


def tool_names() -> list[str]:
    """Return every registered tool name, in registration order."""
    return [entry.name for entry in _TOOLS]


def call_tool(
    name: str,
    arguments: dict[str, Any] | None = None,
    session: Session | None = None,
) -> dict[str, Any]:
    """Call a tool and return a JSON-RPC tool result.

    The result has ``structuredContent`` and ``isError``. No exception escapes.
    """
    session = session if session is not None else Session()
    handler = _HANDLERS.get(name)
    if handler is None:
        return {
            "isError": True,
            "structuredContent": {
                "error": f"unknown tool {name!r}",
                "error_type": "UnknownTool",
                "available_tools": tool_names(),
            },
        }
    try:
        result = handler(session, dict(arguments or {}))
    except Exception as exc:  # noqa: BLE001 - the boundary must not raise
        return {
            "isError": True,
            "structuredContent": {
                "error": str(exc),
                "error_type": type(exc).__name__,
                "tool": name,
            },
        }
    return {"isError": False, "structuredContent": _jsonable(result)}
