"""Spawn the MCP server as a subprocess and validate the JSON-RPC exchange."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

_MCP_ROOT = Path(__file__).resolve().parents[1]


class ServerProcess:
    """A running ``python -m nanocad_mcp`` subprocess."""

    def __init__(self) -> None:
        environment = dict(os.environ)
        environment["PYTHONPATH"] = (
            str(_MCP_ROOT) + os.pathsep + environment.get("PYTHONPATH", "")
        )
        self.process = subprocess.Popen(
            [sys.executable, "-m", "nanocad_mcp"],
            cwd=str(_MCP_ROOT),
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
            env=environment,
        )

    def request(self, message: dict) -> dict:
        """Send one message and read one response line."""
        assert self.process.stdin is not None
        assert self.process.stdout is not None
        self.process.stdin.write(json.dumps(message) + "\n")
        self.process.stdin.flush()
        line = self.process.stdout.readline()
        assert line, "the server closed the connection"
        return json.loads(line)

    def notify(self, message: dict) -> None:
        """Send a notification. The server sends no response."""
        assert self.process.stdin is not None
        self.process.stdin.write(json.dumps(message) + "\n")
        self.process.stdin.flush()

    def close(self) -> None:
        if self.process.stdin is not None:
            self.process.stdin.close()
        self.process.wait(timeout=10)


def test_initialize_tools_list_and_tool_calls() -> None:
    server = ServerProcess()
    try:
        initialize = server.request(
            {
                "jsonrpc": "2.0",
                "id": 1,
                "method": "initialize",
                "params": {
                    "protocolVersion": "2025-06-18",
                    "capabilities": {},
                    "clientInfo": {"name": "pytest", "version": "1.0"},
                },
            }
        )
        assert initialize["jsonrpc"] == "2.0"
        assert initialize["id"] == 1
        assert initialize["result"]["serverInfo"]["name"] == "nanocad"
        assert "tools" in initialize["result"]["capabilities"]

        server.notify(
            {"jsonrpc": "2.0", "method": "notifications/initialized", "params": {}}
        )

        listing = server.request(
            {"jsonrpc": "2.0", "id": 2, "method": "tools/list", "params": {}}
        )
        assert listing["id"] == 2
        tools = listing["result"]["tools"]
        names = {entry["name"] for entry in tools}
        assert {"create_part", "measure", "convert_units"} <= names

        converted = server.request(
            {
                "jsonrpc": "2.0",
                "id": 3,
                "method": "tools/call",
                "params": {
                    "name": "convert_units",
                    "arguments": {"value": 1.0, "from_unit": "nm", "to_unit": "m"},
                },
            }
        )
        assert converted["id"] == 3
        assert converted["result"]["isError"] is False
        text = converted["result"]["content"][0]["text"]
        payload = json.loads(text)
        assert payload["result"] == 1.0e-9
        assert converted["result"]["structuredContent"]["to_unit"] == "m"

        created = server.request(
            {
                "jsonrpc": "2.0",
                "id": 4,
                "method": "tools/call",
                "params": {
                    "name": "create_part",
                    "arguments": {"name": "diamond", "specs": {}},
                },
            }
        )
        assert created["id"] == 4
        assert created["result"]["isError"] is False
        assert created["result"]["structuredContent"]["atom_count"] > 0

        unknown = server.request(
            {
                "jsonrpc": "2.0",
                "id": 5,
                "method": "tools/call",
                "params": {"name": "nope", "arguments": {}},
            }
        )
        assert unknown["result"]["isError"] is True
        assert "nope" in unknown["result"]["structuredContent"]["error"]

        missing = server.request(
            {"jsonrpc": "2.0", "id": 6, "method": "does/not/exist", "params": {}}
        )
        assert missing["error"]["code"] == -32601
    finally:
        server.close()
