"""Guard tests for ``scripts/reproduce.sh``.

The reproduction script must not pass silently. These tests check the POSIX
shell syntax and check that the script exits nonzero when a required tool is
absent from PATH.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
_SCRIPT = _REPO_ROOT / "scripts" / "reproduce.sh"
_SH = shutil.which("sh") or "/bin/sh"


def test_script_has_valid_posix_syntax() -> None:
    result = subprocess.run(
        [_SH, "-n", str(_SCRIPT)], capture_output=True, text=True, check=False
    )
    assert result.returncode == 0, result.stderr


def test_script_fails_when_cargo_is_missing(tmp_path: Path) -> None:
    env = dict(os.environ)
    env["PATH"] = str(tmp_path)
    result = subprocess.run(
        [_SH, str(_SCRIPT)],
        cwd=str(_REPO_ROOT),
        env=env,
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )
    assert result.returncode != 0
    assert "cargo" in (result.stdout + result.stderr).lower()
