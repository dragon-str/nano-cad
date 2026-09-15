"""Pytest path setup for the mcp tests.

The ``nanocad_mcp`` package lives in this directory and is not installed. Add
the directory to ``sys.path`` so the tests import it.
"""

from __future__ import annotations

import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
