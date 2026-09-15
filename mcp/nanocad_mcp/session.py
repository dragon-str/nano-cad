"""In-memory state for one MCP session.

The MCP server is a long-running process. Each tool call names a handle that an
earlier call produced. This module maps those handles to live ``nanocad``
objects. Handles are stable strings: ``part-1``, ``document-2``, ``assembly-1``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

_CARBON_ATOMIC_MASS_U = 12.011


def carbon_atom_mass_kg() -> float:
    """Return the carbon-12 atomic mass in kilograms.

    This is an estimate. It uses the standard atomic weight of carbon, not the
    exact carbon-12 mass.
    """
    import nanocad

    return float(nanocad.to_si(_CARBON_ATOMIC_MASS_U, "u"))


@dataclass
class Session:
    """Live nanocad objects for one MCP connection."""

    parts: dict[str, Any] = field(default_factory=dict)
    documents: dict[str, Any] = field(default_factory=dict)
    assemblies: dict[str, Any] = field(default_factory=dict)
    jigs: dict[str, Any] = field(default_factory=dict)
    _counters: dict[str, int] = field(default_factory=dict)

    def next_id(self, kind: str) -> str:
        """Return the next handle for a kind of object."""
        counter = self._counters.get(kind, 0) + 1
        self._counters[kind] = counter
        return f"{kind}-{counter}"

    def store_part(self, part: Any, part_id: str | None = None) -> str:
        """Store a part and return its handle."""
        handle = part_id or self.next_id("part")
        self.parts[handle] = part
        return handle

    def store_document(self, document: Any, document_id: str | None = None) -> str:
        """Store a document and return its handle."""
        handle = document_id or self.next_id("document")
        self.documents[handle] = document
        return handle

    def store_assembly(self, assembly: Any, assembly_id: str | None = None) -> str:
        """Store a planetary assembly and return its handle."""
        handle = assembly_id or self.next_id("assembly")
        self.assemblies[handle] = assembly
        return handle

    def store_jig(self, jig: Any, jig_id: str | None = None) -> str:
        """Store a jig and return its handle."""
        handle = jig_id or self.next_id("jig")
        self.jigs[handle] = jig
        return handle

    def require_part(self, part_id: str) -> Any:
        """Return a stored part or raise ``ValueError``."""
        if part_id not in self.parts:
            raise ValueError(f"unknown part id {part_id!r}")
        return self.parts[part_id]

    def require_document(self, document_id: str) -> Any:
        """Return a stored document or raise ``ValueError``."""
        if document_id not in self.documents:
            raise ValueError(f"unknown document id {document_id!r}")
        return self.documents[document_id]

    def require_assembly(self, assembly_id: str) -> Any:
        """Return a stored assembly or raise ``ValueError``."""
        if assembly_id not in self.assemblies:
            raise ValueError(f"unknown assembly id {assembly_id!r}")
        return self.assemblies[assembly_id]
