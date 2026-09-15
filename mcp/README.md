# nanocad MCP server

This directory holds the MCP server that wraps the `nanocad` Python tools. An
agent uses the server to build and check parts. The server speaks JSON-RPC 2.0
over stdio.

The server uses the Python standard library only. It adds no dependency.

## Requirements

- Python 3.10 or newer.
- The built `nanocad` extension. Build it from `python/` with maturin.
- `pytest` to run the tests.

## Run

```sh
cd mcp
python -m nanocad_mcp
```

The server reads newline-delimited JSON-RPC messages on stdin and writes one
JSON response on stdout for each request. It writes no other text to stdout.

To use the server from an MCP client, point the client at the module. This is
the shape for a `stdio` server:

```json
{
  "mcpServers": {
    "nanocad": {
      "command": "python",
      "args": ["-m", "nanocad_mcp"],
      "cwd": "/path/to/nano-cad/mcp"
    }
  }
}
```

## Tools

The tool schemas are in `tools.json`. The same schemas are returned by
`tools/list`. Every tool returns structured JSON. A numeric field carries its
unit in the field name, for example `energy_j`, `dt_s`, or `bounding_box_size_m`.

| Tool | Purpose |
|---|---|
| `create_part` | Create a part from a generator id and parameters. |
| `generate_gear` | Generate a gear: `spur_gear`, `gear_profile`, or `planetary`. |
| `generate_nanotube` | Generate a nanotube. |
| `generate_lattice` | Generate a `diamond` or `graphite` lattice. |
| `list_generators` | List the generators and their parameter schemas. |
| `relax` | Minimize a part with harmonic bond-stretch terms. |
| `run_md` | Run Velocity Verlet dynamics on a part. |
| `add_jig` | Add an `anchor` or `spring` jig. |
| `measure` | Measure counts, the bounding box, the centroid, and charge. |
| `validate_part` | Check coordination, bond strain, and valence. |
| `convert_units` | Convert a value between two units. |
| `save` | Save a document as `ncz`, `xyz`, or `mmp`. |
| `load` | Load a document. |
| `list` | List part handles in a document or the session. |
| `assemble` | Group parts into a document. |
| `assemble_planetary` | Build the L2 planetary gearbox device. |
| `measure_gear_ratio` | Drive the sun and measure the gear ratio. |
| `export_urdf` | Export the gearbox to URDF text. |

Tool calls name handles from earlier calls. A part handle is `part-1`, a
document handle is `document-1`, and an assembly handle is `assembly-1`.

## Worked example

The example below builds one diamond part, measures it, and converts a length.
Each request is one line. Each response is one line. The blank lines are only
for reading.

Initialize:

```json
{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2025-06-18","capabilities":{},"clientInfo":{"name":"example","version":"1.0"}}}
```

Response:

```json
{"jsonrpc":"2.0","id":1,"result":{"protocolVersion":"2025-06-18","capabilities":{"tools":{"listChanged":false}},"serverInfo":{"name":"nanocad","version":"0.1.0"}}}
```

Create a part:

```json
{"jsonrpc":"2.0","id":2,"method":"tools/call","params":{"name":"create_part","arguments":{"name":"diamond","specs":{}}}}
```

Response, shortened:

```json
{"jsonrpc":"2.0","id":2,"result":{"content":[{"type":"text","text":"{\"atom_count\": 8, \"bond_count\": 16, \"material\": \"diamond\", \"name\": \"diamond-...\", \"part_id\": \"part-1\"}"}],"structuredContent":{"part_id":"part-1","name":"diamond-...","material":"diamond","atom_count":8,"bond_count":16},"isError":false}}
```

Measure the part:

```json
{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"measure","arguments":{"part_id":"part-1"}}}
```

Convert a length:

```json
{"jsonrpc":"2.0","id":4,"method":"tools/call","params":{"name":"convert_units","arguments":{"value":1.0,"from_unit":"nm","to_unit":"m"}}}
```

Response, shortened:

```json
{"jsonrpc":"2.0","id":4,"result":{"structuredContent":{"value":1.0,"from_unit":"nm","to_unit":"m","result":1e-9},"isError":false}}
```

A bad call returns `isError` true and an `error` field. The server does not
raise across the boundary:

```json
{"jsonrpc":"2.0","id":5,"method":"tools/call","params":{"name":"measure","arguments":{}}}
```

```json
{"jsonrpc":"2.0","id":5,"result":{"structuredContent":{"error":"missing required argument 'part_id'","error_type":"ValueError","tool":"measure"},"isError":true}}
```

## Test

Run from the repository root:

```sh
python -m pytest mcp/tests
```

The server test spawns `python -m nanocad_mcp` as a subprocess. The demo test
builds the planetary gearbox through the tool layer alone.

## Limitations

- The server keeps state in memory. A restart loses every handle.
- `extract_parameters` is not exposed. The thermal-property extraction (M5-05)
  is not done.
- `run_md` uses one carbon mass for every atom. It is an estimate, not a
  per-element mass.
- `relax` and `run_md` build bond-stretch terms only. They do not assign
  angles, torsions, van der Waals, or electrostatic terms.
- The tool list is fixed at startup. The server sends no
  `notifications/tools/list_changed`.
