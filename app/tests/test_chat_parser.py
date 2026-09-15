"""Tests for the rule-based chat parser in `app/chat.py`.

The parser is deterministic and needs no third-party package.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import chat  # noqa: E402


def test_defaults_are_within_limits():
    for key, value in chat.DEFAULT_PARAMS.items():
        low, high = chat.LIMITS[key]
        assert low <= value <= high


def test_one_layer_thicker():
    command = chat.parse("make it one atomic layer thicker")
    assert command.intent == "layers_thicker"
    assert command.changes == {"layers": 5.0}


def test_bare_thicker_defaults_to_one_layer():
    command = chat.parse("thicker")
    assert command.changes == {"layers": 5.0}


def test_atoms_thick_sets_the_layer_count():
    command = chat.parse("6 atoms thick")
    assert command.intent == "layers_set"
    assert command.changes == {"layers": 6.0}


def test_layer_count_is_clamped():
    command = chat.parse("500 atoms thick")
    assert command.changes == {"layers": float(chat.LIMITS["layers"][1])}


def test_more_teeth_changes_sun_and_planet():
    command = chat.parse("add more teeth to the gears")
    assert command.intent == "teeth_more"
    assert command.changes == {"sun_teeth": 30.0, "planet_teeth": 21.0}


def test_planet_count():
    command = chat.parse("4 planets")
    assert command.intent == "planets_set"
    assert command.changes == {"planet_count": 4.0}


def test_layer_spacing_is_fixed_by_the_crystal():
    command = chat.parse("move the layers 20 pm apart")
    assert command.intent is None
    assert command.changes == {}
    assert "crystal" in command.reply


def test_layer_spacing_command_explains_without_a_change():
    command = chat.parse("layer spacing 200 pm")
    assert command.changes == {}
    assert "crystal" in command.reply


def test_reset_returns_the_defaults():
    command = chat.parse("reset")
    assert command.intent == "reset"
    assert command.changes == chat.DEFAULT_PARAMS


def test_unknown_text_is_reported_without_a_change():
    command = chat.parse("hello there")
    assert command.intent is None
    assert command.changes == {}
    assert "did not understand" in command.reply


def test_reply_reports_the_new_state():
    command = chat.parse("6 atoms thick")
    assert "6 atomic layers" in command.reply
    assert "ring 60t" in command.reply
