"""Smoke test for the nanocad package skeleton."""

import nanocad


def test_version_is_non_empty_string() -> None:
    """The package reports a non-empty version string."""
    assert isinstance(nanocad.__version__, str)
    assert nanocad.__version__ != ""
