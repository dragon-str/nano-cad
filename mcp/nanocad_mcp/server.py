"""A hand-rolled MCP server over stdio.

This module speaks JSON-RPC 2.0 as newline-delimited JSON, as the MCP stdio
transport requires. It uses only the Python standard library. It implements
``initialize``, ``tools/list``, ``tools/call``, and ``ping``, and it ignores
notifications. No other dependency is added.
"""

from __future__ import annotations

import json
import sys
from typing import Any, TextIO

from . import __version__
from .session import Session
from .tools import call_tool, list_tools

# The MCP protocol revision this server implements. The tool surface is stable
# across the 2024-11-05 and 2025-06-18 revisions.
PROTOCOL_VERSION = "2025-06-18"

SERVER_NAME = "nanocad"

# JSON-RPC 2.0 error codes.
PARSE_ERROR = -32700
INVALID_REQUEST = -32600
METHOD_NOT_FOUND = -32601
INVALID_PARAMS = -32602
INTERNAL_ERROR = -32603


def _tool_result_payload(result: dict[str, Any]) -> dict[str, Any]:
    """Wrap a tool result as an MCP tool result with JSON text content."""
    structured = result["structuredContent"]
    return {
        "content": [
            {
                "type": "text",
                "text": json.dumps(structured, sort_keys=True, allow_nan=False),
            }
        ],
        "structuredContent": structured,
        "isError": bool(result["isError"]),
    }


class MCPServer:
    """Dispatch JSON-RPC messages to the nanocad tool registry."""

    def __init__(self, session: Session | None = None) -> None:
        self.session = session if session is not None else Session()

    def _result(self, message_id: Any, result: dict[str, Any]) -> dict[str, Any]:
        return {"jsonrpc": "2.0", "id": message_id, "result": result}

    def _error(self, message_id: Any, code: int, message: str) -> dict[str, Any]:
        return {
            "jsonrpc": "2.0",
            "id": message_id,
            "error": {"code": code, "message": message},
        }

    def handle(self, message: dict[str, Any]) -> dict[str, Any] | None:
        """Handle one JSON-RPC message. Return a response, or ``None``."""
        if not isinstance(message, dict) or "method" not in message:
            if "id" not in message:
                return None
            return self._error(message.get("id"), INVALID_REQUEST, "invalid request")

        method = message.get("method")
        message_id = message.get("id")
        params = message.get("params") or {}
        is_request = "id" in message

        if method == "notifications/initialized" or not is_request:
            return None

        if method == "initialize":
            return self._result(
                message_id,
                {
                    "protocolVersion": PROTOCOL_VERSION,
                    "capabilities": {"tools": {"listChanged": False}},
                    "serverInfo": {"name": SERVER_NAME, "version": __version__},
                    "instructions": (
                        "nanocad tools build, simulate, and inspect molecular "
                        "parts. Numeric results carry units in the field name."
                    ),
                },
            )

        if method == "ping":
            return self._result(message_id, {})

        if method == "tools/list":
            return self._result(message_id, list_tools())

        if method == "tools/call":
            if not isinstance(params, dict):
                return self._error(
                    message_id, INVALID_PARAMS, "params must be an object"
                )
            name = params.get("name")
            if not isinstance(name, str):
                return self._error(
                    message_id, INVALID_PARAMS, "params.name must be a string"
                )
            arguments = params.get("arguments")
            if arguments is not None and not isinstance(arguments, dict):
                return self._error(
                    message_id, INVALID_PARAMS, "params.arguments must be an object"
                )
            result = call_tool(name, arguments, self.session)
            return self._result(message_id, _tool_result_payload(result))

        return self._error(
            message_id, METHOD_NOT_FOUND, f"method {method!r} is not supported"
        )

    def serve(self, stdin: TextIO, stdout: TextIO) -> None:
        """Read newline-delimited JSON-RPC messages until ``stdin`` closes."""
        for line in stdin:
            line = line.strip()
            if not line:
                continue
            try:
                message = json.loads(line)
            except json.JSONDecodeError:
                response = self._error(None, PARSE_ERROR, "invalid JSON")
                stdout.write(json.dumps(response) + "\n")
                stdout.flush()
                continue
            response = self.handle(message)
            if response is not None:
                stdout.write(json.dumps(response) + "\n")
                stdout.flush()


def main() -> None:
    """Run the server on the process stdin and stdout."""
    server = MCPServer()
    server.serve(sys.stdin, sys.stdout)


if __name__ == "__main__":  # pragma: no cover
    main()
