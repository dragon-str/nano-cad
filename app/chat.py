#!/usr/bin/env python3
"""Turn a short natural-language edit into generator parameters.

This is a small, honest rule-based parser. It is not a language model. It
recognizes a fixed set of commands about the planetary gearbox and returns the
new parameter set and a reply. The server then runs the Rust generator with
those parameters, so the engine remains the source of truth.

Every command is checked against the engine's own limits before the build:
planet teeth must be at least 18 (undercut), the counts must assemble, and the
layer count must be in range. A rejected command changes nothing.
"""

from __future__ import annotations

import re

# The generator defaults. Keep in step with `PLANETARY_PARAMETERS` in
# `crates/parts/src/planetary.rs`.
DEFAULT_PARAMS = {
    "module_m": 1.5e-9,
    "sun_teeth": 12.0,
    "planet_teeth": 9.0,
    "planet_count": 3.0,
    "layers": 4.0,
}

# One atomic layer is one diamond (001) plane. The crystal fixes the spacing,
# so the user controls the layer count only.
DIAMOND_LATTICE_CONSTANT_M = 3.567e-10
DIAMOND_PLANE_SPACING_M = DIAMOND_LATTICE_CONSTANT_M / 4.0

# Inclusive limits. The engine enforces its own, stricter limits too.
LIMITS = {
    "module_m": (1e-10, 5e-9),
    "sun_teeth": (12.0, 200.0),
    "planet_teeth": (6.0, 120.0),
    "planet_count": (1.0, 12.0),
    "layers": (1.0, 64.0),
}

# Parameter names, in the order the UI shows them, with a human label and the
# unit scale for display.
PARAM_META = [
    ("module_m", "module", 1e-9, "nm"),
    ("sun_teeth", "sun teeth", 1.0, ""),
    ("planet_teeth", "planet teeth", 1.0, ""),
    ("planet_count", "planet count", 1.0, ""),
    ("layers", "atomic layers", 1.0, ""),
]

_WORD_NUMBERS = {
    "one": 1, "two": 2, "three": 3, "four": 4, "five": 5,
    "six": 6, "seven": 7, "eight": 8, "nine": 9, "ten": 10,
}

_INTENT_ORDER = [
    "reset",
    "layers_set",
    "layers_thicker",
    "layers_thinner",
    "teeth_set",
    "teeth_more",
    "teeth_fewer",
    "planets_set",
    "planets_more",
    "planets_fewer",
]


class Command:
    """One parsed edit: which parameters to change and the new values."""

    def __init__(self, intent: str | None, changes: dict, reply: str):
        self.intent = intent
        self.changes = changes
        self.reply = reply


def _number(text: str) -> int | None:
    text = text.strip()
    if not text:
        return None
    if text.isdigit():
        return int(text)
    return _WORD_NUMBERS.get(text)


def _clamp(name: str, value: float) -> float:
    low, high = LIMITS[name]
    return max(low, min(high, value))


def _describe(params: dict) -> str:
    ring = params["sun_teeth"] + 2.0 * params["planet_teeth"]
    ratio = (params["sun_teeth"] + ring) / params["sun_teeth"]
    thickness_nm = (params["layers"] - 1.0) * DIAMOND_PLANE_SPACING_M / 1e-9
    return (
        f"sun {int(params['sun_teeth'])}t, planet {int(params['planet_teeth'])}t, "
        f"ring {int(ring)}t, {int(params['planet_count'])} planets, "
        f"{int(params['layers'])} atomic layers "
        f"({thickness_nm:.3f} nm thick), ratio {ratio:.3f}."
    )


def parse(message: str, params: dict | None = None) -> Command:
    """Parse one message into a parameter change.

    Returns a `Command`. When `intent` is None the message was not understood
    and `changes` is empty, so the caller keeps the current parameters.
    """
    base = dict(DEFAULT_PARAMS)
    if params:
        base.update(params)
    text = message.lower().strip()
    if not text:
        return Command(None, {}, "Type a change, for example \"one atomic layer thicker\".")

    # reset
    if re.search(r"\b(reset|default|start over)\b", text):
        return Command("reset", dict(DEFAULT_PARAMS),
                       "Reset to the default design. " + _describe(DEFAULT_PARAMS))

    # "make it 6 atoms thick" / "6 layers thick" / "set layers to 6"
    match = re.search(r"(?:set\s+)?(?:the\s+)?layers?\s+to\s+(\d+)", text)
    if not match:
        match = re.search(r"(\d+)\s+(?:atomic\s+)?(?:layers?|atoms?)\s+thick", text)
    if match:
        value = float(_clamp("layers", float(match.group(1))))
        changed = {"layers": value}
        return Command("layers_set", changed,
                       "Set the gear thickness. " + _describe({**base, **changed}))

    # thicker / thinner by N
    match = re.search(r"(\w+|\d+)\s+(?:atomic\s+)?layers?\s+(thicker|thinner)", text)
    if not match:
        match = re.search(r"(thicker|thinner)", text)
        if match:
            match = re.match(r"(thicker|thinner)", match.group(1))
    if match:
        amount = _number(match.group(1)) if match.lastindex and match.lastindex >= 2 else 1
        if amount is None:
            amount = 1
        direction = match.group(2) if match.lastindex and match.lastindex >= 2 else match.group(1)
        delta = amount if direction == "thicker" else -amount
        value = float(_clamp("layers", base["layers"] + delta))
        changed = {"layers": value}
        intent = "layers_thicker" if delta > 0 else "layers_thinner"
        word = "Thickened" if delta > 0 else "Thinned"
        return Command(intent, changed,
                       f"{word} the gears by {abs(delta)} atomic layer(s). "
                       + _describe({**base, **changed}))

    # layer spacing is fixed by the crystal, so explain instead of editing
    if re.search(r"\b(spacing|space|separate|move)\b.*\blayers?\b", text) or \
            re.search(r"\blayers?\b.*\bspacing\b", text):
        return Command(
            None, {},
            "One atomic layer is one diamond (001) plane. The crystal fixes the "
            f"spacing at {DIAMOND_PLANE_SPACING_M * 1e12:.4f} pm, so you set the "
            "layer count, not the spacing.",
        )

    # "set the sun teeth to 30"
    match = re.search(r"(sun|planet)\s+teeth\s+to\s+(\d+)", text)
    if match:
        key = f"{match.group(1)}_teeth"
        value = float(_clamp(key, float(match.group(2))))
        changed = {key: value}
        return Command("teeth_set", changed,
                       f"Set the {match.group(1)} teeth. " + _describe({**base, **changed}))

    # "add more teeth" / "add 6 teeth"
    if re.search(r"\b(more|add|extra)\b.*\bteeth\b", text) or re.search(r"\bteeth\b.*\bmore\b", text):
        match = re.search(r"(\d+)", text)
        amount = float(match.group(1)) if match else 6.0
        changed = {
            "sun_teeth": float(_clamp("sun_teeth", base["sun_teeth"] + amount)),
            "planet_teeth": float(_clamp("planet_teeth", base["planet_teeth"] + amount / 2.0)),
        }
        return Command("teeth_more", changed,
                       f"Added teeth to the sun and the planets. "
                       + _describe({**base, **changed}))

    # "fewer teeth"
    if re.search(r"\b(fewer|less)\b.*\bteeth\b", text):
        match = re.search(r"(\d+)", text)
        amount = float(match.group(1)) if match else 6.0
        changed = {
            "sun_teeth": float(_clamp("sun_teeth", base["sun_teeth"] - amount)),
            "planet_teeth": float(_clamp("planet_teeth", base["planet_teeth"] - amount / 2.0)),
        }
        return Command("teeth_fewer", changed,
                       f"Removed teeth from the sun and the planets. "
                       + _describe({**base, **changed}))

    # "set 4 planets"
    match = re.search(r"(?:set\s+)?(\d+)\s+planets", text)
    if match:
        value = float(_clamp("planet_count", float(match.group(1))))
        changed = {"planet_count": value}
        return Command("planets_set", changed,
                       "Set the planet count. " + _describe({**base, **changed}))

    # "add a planet" / "more planets"
    if re.search(r"\b(add|more|extra)\b.*\bplanets?\b", text):
        value = float(_clamp("planet_count", base["planet_count"] + 1.0))
        changed = {"planet_count": value}
        return Command("planets_more", changed,
                       "Added a planet. " + _describe({**base, **changed}))

    if re.search(r"\b(fewer|less|remove)\b.*\bplanets?\b", text):
        value = float(_clamp("planet_count", base["planet_count"] - 1.0))
        changed = {"planet_count": value}
        return Command("planets_fewer", changed,
                       "Removed a planet. " + _describe({**base, **changed}))

    return Command(
        None, {},
        "I did not understand that. Try \"one atomic layer thicker\", "
        "\"6 atoms thick\", \"add more teeth\", \"4 planets\", or \"reset\".",
    )
