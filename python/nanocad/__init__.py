"""nanocad: an open multiscale CAD and simulation platform.

The compiled extension module ``nanocad._core`` holds the Rust core. This
package imports it and re-exports its public names. When the extension is
absent the package still imports and reports a version string.
"""

from __future__ import annotations

__version__ = "0.1.0"


def version() -> str:
    """Return the nanocad version string."""
    return __version__


try:
    from . import _core
except ImportError:  # pragma: no cover - the extension is built by maturin
    _core = None

if _core is not None:
    __version__ = _core.version()
    for _name in dir(_core):
        if not _name.startswith("_"):
            globals()[_name] = getattr(_core, _name)
    del _name

__all__ = ["__version__", "version"] + (
    sorted(name for name in dir(_core) if not name.startswith("_"))
    if _core is not None
    else []
)
