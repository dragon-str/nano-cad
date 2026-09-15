"""Lumped-parameter ODE models with SBML interchange (task L4-01).

This module implements one thin vertical slice of the L4 system layer: a
physically meaningful lumped model, a fixed-step integrator, an exact analytic
cross-check, and a round trip through a documented SBML-like XML subset.

Model: two-compartment respirocyte gas transport
------------------------------------------------
A high-pressure source compartment drains through a resistive channel into a
low-pressure sink compartment.  Each compartment stores volume V [m^3] that is
proportional to its pressure p [Pa] through a capacitance C [m^3/Pa].  The
channel carries volumetric flow Q [m^3/s] proportional to the pressure drop.

Governing ODEs, with state vector (p_source, p_sink) in Pa::

    C_source * dp_source/dt = -(p_source - p_sink) / R
    C_sink   * dp_sink/dt   = +(p_source - p_sink) / R

Equivalently::

    dp_source/dt = -(p_source - p_sink) / (R * C_source)
    dp_sink/dt   = +(p_source - p_sink) / (R * C_sink)

The system is linear.  The total stored volume ``C_source * p_source +
C_sink * p_sink`` is conserved, and both pressures relax exponentially to the
volume-weighted mean pressure with rate constant
``k = 1/(R*C_source) + 1/(R*C_sink)``.  That closed form is the exact analytic
solution used for the cross-check.

SI units.  The default parameter names carry their units:
``resistance_pa_s_per_m3`` [Pa*s/m^3], ``capacitance_source_m3_per_pa`` and
``capacitance_sink_m3_per_pa`` [m^3/Pa].

SBML subset
-----------
``to_sbml`` emits SBML Level 3 Version 1 elements: ``sbml``, ``model``,
``listOfCompartments``/``compartment``, ``listOfSpecies``/``species``,
``listOfParameters``/``parameter``, ``listOfReactions``/``reaction`` with
``listOfReactants``, ``listOfProducts`` and ``kineticLaw``, plus a MathML
``math`` expression over ``cn``, ``ci``, ``apply``, ``plus``, ``minus``,
``times``, ``divide`` and ``power``.  Two extensions keep the document honest:

* each ``species`` and ``parameter`` carries a ``nanocad:unit`` attribute in
  the ``https://nanocad.org/sbml-ext`` namespace;
* a compartment is a normalised holder (``size="1"``), because this model is
  not a concentration model.  Pressure is stored in ``initialConcentration``.

The subset is deliberately small.  The document is well-formed XML and uses
only the elements listed above.  It is not schema-validated against the full
SBML schema, and it does not claim SBML compliance beyond that subset.

CellML subset
-------------
``to_cellml`` emits CellML 2.0 elements: ``model``, one ``component``, a
``variable`` for time, each state and each parameter, and one ``math`` with
``apply``/``eq``/``diff``/``bvar`` and MathML ``cn``/``ci``/``apply``.  CellML
has no reaction element, so the writer folds the reactions into one derivative
equation per state.  ``from_cellml`` rebuilds one reaction per state from that
derivative, so the ODE system round-trips but the original reaction split does
not.  The CellML ``units`` attribute is ``dimensionless``; the real unit
string rides in a ``nanocad:unit`` attribute, as in the SBML path.

The subset is deliberately small.  It is not schema-validated, and it does not
claim CellML compliance beyond the listed elements.
"""

from __future__ import annotations

import math
import re
import xml.etree.ElementTree as ET
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass

__all__ = [
    "BinOp",
    "ExpressionError",
    "LumpedError",
    "Model",
    "ModelFormatError",
    "Neg",
    "Num",
    "Parameter",
    "Reaction",
    "SimulationResult",
    "Species",
    "Var",
    "XmlParseError",
    "derivatives",
    "evaluate",
    "from_cellml",
    "from_sbml",
    "initial_state",
    "parse_expression",
    "respirocyte_gas_transport",
    "rk4_step",
    "simulate",
    "species_ids",
    "to_cellml",
    "to_sbml",
    "two_compartment_pressures",
    "variables",
]

_SBML_NS = "http://www.sbml.org/sbml/level3/version1/core"
_CELLML_NS = "http://www.cellml.org/cellml/2.0#"
_MATHML_NS = "http://www.w3.org/1998/Math/MathML"
_NANOCAD_NS = "https://nanocad.org/sbml-ext"
_CELLML_NANOCAD_NS = "https://nanocad.org/cellml-ext"
_OP_TO_TAG = {"+": "plus", "-": "minus", "*": "times", "/": "divide", "^": "power"}
_TAG_TO_OP = {tag: op for op, tag in _OP_TO_TAG.items()}


# ---------------------------------------------------------------------------
# Typed errors.  Library code never lets a parse failure escape as a bare
# exception such as ET.ParseError, KeyError or ZeroDivisionError.
# ---------------------------------------------------------------------------
class LumpedError(Exception):
    """Base class for every error raised by this module."""


class ExpressionError(LumpedError):
    """A rate expression is empty, malformed, or uses an unknown name."""


class XmlParseError(LumpedError):
    """The supplied text is not well-formed XML."""


class ModelFormatError(LumpedError):
    """The XML is well-formed but is not a supported model document."""


# ---------------------------------------------------------------------------
# Rate-law expression tree and a small recursive-descent parser.
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Num:
    """A numeric literal."""

    value: float

    def to_infix(self) -> str:
        return repr(self.value)


@dataclass(frozen=True)
class Var:
    """A reference to a species id or a parameter id."""

    name: str

    def to_infix(self) -> str:
        return self.name


@dataclass(frozen=True)
class Neg:
    """Arithmetic negation."""

    operand: Expr

    def to_infix(self) -> str:
        return f"(-{self.operand.to_infix()})"


@dataclass(frozen=True)
class BinOp:
    """A binary operation.  ``op`` is one of ``+ - * / ^``."""

    op: str
    left: Expr
    right: Expr

    def to_infix(self) -> str:
        return f"({self.left.to_infix()} {self.op} {self.right.to_infix()})"


Expr = Num | Var | Neg | BinOp

_TOKEN_RE = re.compile(
    r"\s*(?:(\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?|([A-Za-z_][A-Za-z0-9_]*)|([()+\-*/^]))"
)


def _tokenize(text: str) -> list[tuple[str, object]]:
    tokens: list[tuple[str, object]] = []
    position = 0
    length = len(text)
    while position < length:
        match = _TOKEN_RE.match(text, position)
        if match is None:
            if text[position].isspace():
                position += 1
                continue
            raise ExpressionError(f"unexpected character {text[position]!r}")
        position = match.end()
        number, name, symbol = match.groups()
        if number is not None:
            tokens.append(("num", float(number)))
        elif name is not None:
            tokens.append(("id", name))
        else:
            tokens.append((symbol, symbol))
    tokens.append(("eof", ""))
    return tokens


class _ExpressionParser:
    """Recursive-descent parser for the tiny rate-law grammar."""

    def __init__(self, tokens: Sequence[tuple[str, object]]) -> None:
        self._tokens = tokens
        self._position = 0

    def parse(self) -> Expr:
        expression = self._parse_sum()
        if self._peek()[0] != "eof":
            raise ExpressionError(f"unexpected trailing token {self._peek()[1]!r}")
        return expression

    def _peek(self) -> tuple[str, object]:
        return self._tokens[self._position]

    def _next(self) -> tuple[str, object]:
        token = self._tokens[self._position]
        self._position += 1
        return token

    def _parse_sum(self) -> Expr:
        expression = self._parse_product()
        while self._peek()[0] in ("+", "-"):
            op = str(self._next()[0])
            expression = BinOp(op, expression, self._parse_product())
        return expression

    def _parse_product(self) -> Expr:
        expression = self._parse_power()
        while self._peek()[0] in ("*", "/"):
            op = str(self._next()[0])
            expression = BinOp(op, expression, self._parse_power())
        return expression

    def _parse_power(self) -> Expr:
        expression = self._parse_unary()
        if self._peek()[0] == "^":
            self._next()
            expression = BinOp("^", expression, self._parse_power())
        return expression

    def _parse_unary(self) -> Expr:
        kind = self._peek()[0]
        if kind in ("+", "-"):
            self._next()
            operand = self._parse_unary()
            return operand if kind == "+" else Neg(operand)
        return self._parse_primary()

    def _parse_primary(self) -> Expr:
        kind, value = self._next()
        if kind == "num":
            return Num(float(value))
        if kind == "id":
            return Var(str(value))
        if kind == "(":
            expression = self._parse_sum()
            if self._peek()[0] != ")":
                raise ExpressionError("missing closing parenthesis")
            self._next()
            return expression
        raise ExpressionError(f"unexpected token {value!r}")


def parse_expression(text: str) -> Expr:
    """Parse a rate law over species and parameter ids.

    The grammar supports numbers, identifiers, ``+ - * / ^``, parentheses and
    unary signs.  Raise :class:`ExpressionError` on any other input.
    """
    if not isinstance(text, str) or not text.strip():
        raise ExpressionError("empty expression")
    return _ExpressionParser(_tokenize(text)).parse()


def variables(expr: Expr) -> frozenset[str]:
    """Return the set of every variable name that appears in ``expr``."""
    if isinstance(expr, Var):
        return frozenset({expr.name})
    if isinstance(expr, Neg):
        return variables(expr.operand)
    if isinstance(expr, BinOp):
        return variables(expr.left) | variables(expr.right)
    return frozenset()


def evaluate(expr: Expr, environment: Mapping[str, float]) -> float:
    """Evaluate ``expr`` against a name-to-value mapping.

    Raise :class:`ExpressionError` on an unknown name or a division by zero.
    """
    if isinstance(expr, Num):
        return expr.value
    if isinstance(expr, Var):
        try:
            return float(environment[expr.name])
        except KeyError as exc:
            raise ExpressionError(f"unknown variable {expr.name!r}") from exc
    if isinstance(expr, Neg):
        return -evaluate(expr.operand, environment)
    left = evaluate(expr.left, environment)
    right = evaluate(expr.right, environment)
    if expr.op == "+":
        return left + right
    if expr.op == "-":
        return left - right
    if expr.op == "*":
        return left * right
    if expr.op == "^":
        return left**right
    if expr.op == "/":
        if right == 0.0:
            raise ExpressionError("division by zero in rate law")
        return left / right
    raise ExpressionError(f"unknown operator {expr.op!r}")


# ---------------------------------------------------------------------------
# Model data contract.
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Species:
    """A lumped state variable.  ``initial_value_si`` carries SI units."""

    id: str
    name: str
    initial_value_si: float
    unit: str = "Pa"
    compartment: str = "environment"


@dataclass(frozen=True)
class Parameter:
    """A constant model parameter with an explicit SI unit."""

    id: str
    name: str
    value_si: float
    unit: str


@dataclass(frozen=True)
class Reaction:
    """A rate law plus the net stoichiometry it drives.

    ``stoichiometry`` holds ``(species_id, net_coefficient)`` pairs.  The ODE
    term is ``net_coefficient * rate_law`` for each species.  The tuple is
    normalised to sorted order so that a parse round trip compares equal.
    """

    id: str
    name: str
    stoichiometry: tuple[tuple[str, float], ...]
    rate_law: Expr

    def __post_init__(self) -> None:
        object.__setattr__(self, "stoichiometry", tuple(sorted(self.stoichiometry)))


@dataclass(frozen=True)
class Model:
    """A lumped model: species, parameters and reactions."""

    id: str
    name: str
    species: tuple[Species, ...]
    parameters: tuple[Parameter, ...]
    reactions: tuple[Reaction, ...]
    method: str = "rk4-fixed-step"
    validation: str = "unverified"


@dataclass(frozen=True)
class SimulationResult:
    """Time series with provenance.  ``method`` and ``validation`` are always set."""

    times_s: tuple[float, ...]
    states: tuple[tuple[float, ...], ...]
    species_ids: tuple[str, ...]
    method: str
    validation: str
    relative_error: float | None = None


def species_ids(model: Model) -> tuple[str, ...]:
    """Return the species ids in state-vector order."""
    return tuple(species.id for species in model.species)


def initial_state(model: Model) -> tuple[float, ...]:
    """Return the initial state vector in ``species_ids`` order."""
    return tuple(species.initial_value_si for species in model.species)


def derivatives(model: Model, state: Sequence[float]) -> list[float]:
    """Return ``d(state)/dt`` for the model at ``state``.

    Raise :class:`ModelFormatError` if a reaction names an unknown species.
    """
    ids = species_ids(model)
    if len(state) != len(ids):
        raise ValueError(f"state has {len(state)} values, expected {len(ids)}")
    environment: dict[str, float] = {
        species_id: float(value) for species_id, value in zip(ids, state, strict=True)
    }
    environment.update(
        {parameter.id: parameter.value_si for parameter in model.parameters}
    )
    index = {species_id: position for position, species_id in enumerate(ids)}
    rate_of_change = [0.0] * len(ids)
    for reaction in model.reactions:
        rate = evaluate(reaction.rate_law, environment)
        for species_id, coefficient in reaction.stoichiometry:
            if species_id not in index:
                raise ModelFormatError(
                    f"reaction {reaction.id!r} names unknown species {species_id!r}"
                )
            rate_of_change[index[species_id]] += coefficient * rate
    return rate_of_change


# ---------------------------------------------------------------------------
# Fixed-step RK4 and the exact two-compartment solution.
# ---------------------------------------------------------------------------
def rk4_step(
    function: Callable[[float, Sequence[float]], Sequence[float]],
    state: Sequence[float],
    time_s: float,
    step_s: float,
) -> list[float]:
    """Advance one classical Runge-Kutta 4 step of size ``step_s``."""
    k1 = function(time_s, state)
    k2 = function(time_s + 0.5 * step_s, _axpy(state, 0.5 * step_s, k1))
    k3 = function(time_s + 0.5 * step_s, _axpy(state, 0.5 * step_s, k2))
    k4 = function(time_s + step_s, _axpy(state, step_s, k3))
    return [
        value + step_s / 6.0 * (a + 2.0 * b + 2.0 * c + d)
        for value, a, b, c, d in zip(state, k1, k2, k3, k4, strict=True)
    ]


def _axpy(
    state: Sequence[float], scale: float, direction: Sequence[float]
) -> list[float]:
    return [
        value + scale * delta for value, delta in zip(state, direction, strict=True)
    ]


def simulate(
    model: Model,
    step_s: float,
    step_count: int,
    *,
    analytic: Callable[[float], Sequence[float]] | None = None,
) -> SimulationResult:
    """Integrate ``model`` with fixed-step RK4.

    When ``analytic`` is given it returns the exact state at a time in seconds.
    The result then carries the maximum relative error across all species and
    steps.  ``validation`` becomes ``"cross-checked"`` only when that error is
    at or below ``1e-6``; a larger error gives ``"analytic-mismatch"``.
    """
    if step_s <= 0.0:
        raise ValueError("step_s must be positive")
    if step_count < 0:
        raise ValueError("step_count must not be negative")

    def function(time_s: float, state: Sequence[float]) -> list[float]:
        return derivatives(model, state)

    times = [0.0]
    states = [initial_state(model)]
    time_s = 0.0
    state = list(states[0])
    for _ in range(step_count):
        state = rk4_step(function, state, time_s, step_s)
        time_s += step_s
        times.append(time_s)
        states.append(tuple(state))

    relative_error: float | None = None
    validation = "unverified"
    if analytic is not None:
        relative_error = _maximum_relative_error(times, states, analytic)
        validation = "cross-checked" if relative_error <= 1e-6 else "analytic-mismatch"
    return SimulationResult(
        times_s=tuple(times),
        states=tuple(states),
        species_ids=species_ids(model),
        method=model.method,
        validation=validation,
        relative_error=relative_error,
    )


def _maximum_relative_error(
    times: Sequence[float],
    states: Sequence[Sequence[float]],
    analytic: Callable[[float], Sequence[float]],
) -> float:
    worst = 0.0
    for time_s, state in zip(times, states, strict=True):
        exact = analytic(time_s)
        if len(exact) != len(state):
            raise ValueError("analytic solution has the wrong number of components")
        for numeric, reference in zip(state, exact, strict=True):
            scale = abs(reference) if abs(reference) > 1e-30 else 1.0
            worst = max(worst, abs(numeric - reference) / scale)
    return worst


def two_compartment_pressures(
    pressure_source_initial_pa: float,
    pressure_sink_initial_pa: float,
    capacitance_source_m3_per_pa: float,
    capacitance_sink_m3_per_pa: float,
    resistance_pa_s_per_m3: float,
    time_s: float,
) -> tuple[float, float]:
    """Exact pressures ``(p_source, p_sink)`` at ``time_s`` for the linear model.

    The linear system has eigenvalues ``0`` and
    ``-(1/(R*C_source) + 1/(R*C_sink))``.  The zero mode is the conserved
    volume-weighted mean pressure ``p_eq``.
    """
    rate_source = 1.0 / (resistance_pa_s_per_m3 * capacitance_source_m3_per_pa)
    rate_sink = 1.0 / (resistance_pa_s_per_m3 * capacitance_sink_m3_per_pa)
    decay_rate = rate_source + rate_sink
    equilibrium_pa = (
        capacitance_source_m3_per_pa * pressure_source_initial_pa
        + capacitance_sink_m3_per_pa * pressure_sink_initial_pa
    ) / (capacitance_source_m3_per_pa + capacitance_sink_m3_per_pa)
    deviation_source = pressure_source_initial_pa - equilibrium_pa
    deviation_sink = pressure_sink_initial_pa - equilibrium_pa
    coefficient = (deviation_source - deviation_sink) / decay_rate
    decay = math.exp(-decay_rate * time_s)
    pressure_source_pa = (
        equilibrium_pa + deviation_source - (coefficient * rate_source * (1.0 - decay))
    )
    pressure_sink_pa = (
        equilibrium_pa + deviation_sink + (coefficient * rate_sink * (1.0 - decay))
    )
    return pressure_source_pa, pressure_sink_pa


# ---------------------------------------------------------------------------
# The respirocyte factory.
# ---------------------------------------------------------------------------
def respirocyte_gas_transport(
    pressure_source_initial_pa: float = 1.0e6,
    pressure_sink_initial_pa: float = 1.0e5,
    capacitance_source_m3_per_pa: float = 1.0e-12,
    capacitance_sink_m3_per_pa: float = 5.0e-12,
    resistance_pa_s_per_m3: float = 2.0e14,
    model_id: str = "respirocyte_gas_transport",
) -> Model:
    """Build the two-compartment respirocyte gas-transport model.

    The source drains into the sink through one resistive channel.  Two
    reactions encode the one physical channel: one reaction removes pressure
    from the source, and the other adds pressure to the sink.  Both reactions
    carry the same flow, so the capacitance of each compartment appears as the
    explicit divisor in its own rate law.
    """
    source_id = "pressure_source_pa"
    sink_id = "pressure_sink_pa"
    resistance_id = "resistance_pa_s_per_m3"
    source_capacitance_id = "capacitance_source_m3_per_pa"
    sink_capacitance_id = "capacitance_sink_m3_per_pa"

    species = (
        Species(
            id=source_id,
            name="Source compartment pressure",
            initial_value_si=pressure_source_initial_pa,
            unit="Pa",
            compartment="source",
        ),
        Species(
            id=sink_id,
            name="Sink compartment pressure",
            initial_value_si=pressure_sink_initial_pa,
            unit="Pa",
            compartment="sink",
        ),
    )
    parameters = (
        Parameter(
            id=resistance_id,
            name="Channel hydraulic resistance",
            value_si=resistance_pa_s_per_m3,
            unit="Pa*s/m^3",
        ),
        Parameter(
            id=source_capacitance_id,
            name="Source compartment capacitance",
            value_si=capacitance_source_m3_per_pa,
            unit="m^3/Pa",
        ),
        Parameter(
            id=sink_capacitance_id,
            name="Sink compartment capacitance",
            value_si=capacitance_sink_m3_per_pa,
            unit="m^3/Pa",
        ),
    )
    flow = f"({source_id} - {sink_id}) / {resistance_id}"
    reactions = (
        Reaction(
            id="source_drains_to_sink_source",
            name="Pressure leaves the source compartment",
            stoichiometry=((source_id, -1.0),),
            rate_law=parse_expression(f"{flow} / {source_capacitance_id}"),
        ),
        Reaction(
            id="source_drains_to_sink_sink",
            name="Pressure enters the sink compartment",
            stoichiometry=((sink_id, 1.0),),
            rate_law=parse_expression(f"{flow} / {sink_capacitance_id}"),
        ),
    )
    return Model(
        id=model_id,
        name="Respirocyte two-compartment gas transport",
        species=species,
        parameters=parameters,
        reactions=reactions,
    )


# ---------------------------------------------------------------------------
# SBML serialisation.  See the module docstring for the supported subset.
# ---------------------------------------------------------------------------
def to_sbml(model: Model) -> str:
    """Serialize ``model`` to a documented SBML-like XML document."""
    root = ET.Element(
        "sbml",
        {
            "xmlns": _SBML_NS,
            "xmlns:nanocad": _NANOCAD_NS,
            "level": "3",
            "version": "1",
        },
    )
    model_element = ET.SubElement(root, "model", {"id": model.id, "name": model.name})

    compartments = sorted({species.compartment for species in model.species})
    compartment_list = ET.SubElement(model_element, "listOfCompartments")
    for compartment_id in compartments:
        ET.SubElement(
            compartment_list,
            "compartment",
            {
                "id": compartment_id,
                "size": "1",
                "constant": "true",
                "nanocad:unit": "dimensionless",
            },
        )

    species_list = ET.SubElement(model_element, "listOfSpecies")
    for species in model.species:
        ET.SubElement(
            species_list,
            "species",
            {
                "id": species.id,
                "name": species.name,
                "compartment": species.compartment,
                "initialConcentration": repr(species.initial_value_si),
                "constant": "false",
                "boundaryCondition": "false",
                "hasOnlySubstanceUnits": "false",
                "nanocad:unit": species.unit,
            },
        )

    parameter_list = ET.SubElement(model_element, "listOfParameters")
    for parameter in model.parameters:
        ET.SubElement(
            parameter_list,
            "parameter",
            {
                "id": parameter.id,
                "name": parameter.name,
                "value": repr(parameter.value_si),
                "constant": "true",
                "nanocad:unit": parameter.unit,
            },
        )

    reaction_list = ET.SubElement(model_element, "listOfReactions")
    for reaction in model.reactions:
        reaction_element = ET.SubElement(
            reaction_list,
            "reaction",
            {"id": reaction.id, "name": reaction.name, "reversible": "true"},
        )
        reactants = [
            (species_id, coefficient)
            for species_id, coefficient in reaction.stoichiometry
            if coefficient < 0.0
        ]
        products = [
            (species_id, coefficient)
            for species_id, coefficient in reaction.stoichiometry
            if coefficient > 0.0
        ]
        if reactants:
            reactant_list = ET.SubElement(reaction_element, "listOfReactants")
            for species_id, coefficient in reactants:
                ET.SubElement(
                    reactant_list,
                    "speciesReference",
                    {
                        "species": species_id,
                        "stoichiometry": repr(-coefficient),
                        "constant": "true",
                    },
                )
        if products:
            product_list = ET.SubElement(reaction_element, "listOfProducts")
            for species_id, coefficient in products:
                ET.SubElement(
                    product_list,
                    "speciesReference",
                    {
                        "species": species_id,
                        "stoichiometry": repr(coefficient),
                        "constant": "true",
                    },
                )
        kinetic_law = ET.SubElement(reaction_element, "kineticLaw")
        math_element = ET.SubElement(kinetic_law, "math", {"xmlns": _MATHML_NS})
        math_element.append(_to_mathml(reaction.rate_law))

    ET.indent(root, space="  ")
    body = ET.tostring(root, encoding="unicode")
    return '<?xml version="1.0" encoding="UTF-8"?>\n' + body


def from_sbml(xml: str) -> Model:
    """Parse an SBML-like document into a :class:`Model`.

    Raise :class:`XmlParseError` for malformed XML.  Raise
    :class:`ModelFormatError` or :class:`ExpressionError` when the document is
    well-formed but outside the supported subset.
    """
    try:
        root = ET.fromstring(xml)
    except ET.ParseError as exc:
        raise XmlParseError(f"malformed XML: {exc}") from exc

    if _local(root.tag) != "sbml":
        raise ModelFormatError(
            f"expected <sbml> root element, found <{_local(root.tag)}>"
        )
    model_element = _find_child(root, "model")
    if model_element is None:
        raise ModelFormatError("document has no <model> element")

    model_id = model_element.get("id", "model")
    name = model_element.get("name", model_id)

    species_element = _find_child(model_element, "listOfSpecies")
    if species_element is None:
        raise ModelFormatError("model has no <listOfSpecies>")
    species = tuple(
        _parse_species(element) for element in _children(species_element, "species")
    )

    parameters: tuple[Parameter, ...] = ()
    parameter_element = _find_child(model_element, "listOfParameters")
    if parameter_element is not None:
        parameters = tuple(
            _parse_parameter(element)
            for element in _children(parameter_element, "parameter")
        )

    reactions: tuple[Reaction, ...] = ()
    reaction_element = _find_child(model_element, "listOfReactions")
    if reaction_element is not None:
        reactions = tuple(
            _parse_reaction(element)
            for element in _children(reaction_element, "reaction")
        )

    model = Model(
        id=model_id,
        name=name,
        species=species,
        parameters=parameters,
        reactions=reactions,
    )
    _validate_model(model)
    return model


def _parse_species(element: ET.Element) -> Species:
    species_id = element.get("id")
    if not species_id:
        raise ModelFormatError("a <species> has no id")
    text = element.get("initialConcentration", element.get("initialAmount", "0"))
    try:
        initial_value = float(text)
    except ValueError as exc:
        raise ModelFormatError(
            f"species {species_id!r} has a non-numeric initial value {text!r}"
        ) from exc
    return Species(
        id=species_id,
        name=element.get("name", species_id),
        initial_value_si=initial_value,
        unit=_unit_of(element),
        compartment=element.get("compartment", "environment"),
    )


def _parse_parameter(element: ET.Element) -> Parameter:
    parameter_id = element.get("id")
    if not parameter_id:
        raise ModelFormatError("a <parameter> has no id")
    text = element.get("value", "0")
    try:
        value = float(text)
    except ValueError as exc:
        raise ModelFormatError(
            f"parameter {parameter_id!r} has a non-numeric value {text!r}"
        ) from exc
    return Parameter(
        id=parameter_id,
        name=element.get("name", parameter_id),
        value_si=value,
        unit=_unit_of(element),
    )


def _parse_reaction(element: ET.Element) -> Reaction:
    reaction_id = element.get("id")
    if not reaction_id:
        raise ModelFormatError("a <reaction> has no id")
    stoichiometry: dict[str, float] = {}
    for list_name, sign in (("listOfReactants", -1.0), ("listOfProducts", 1.0)):
        list_element = _find_child(element, list_name)
        if list_element is None:
            continue
        for reference in _children(list_element, "speciesReference"):
            species_id = reference.get("species")
            if not species_id:
                raise ModelFormatError(
                    f"reaction {reaction_id!r} has a speciesReference with no species"
                )
            text = reference.get("stoichiometry", "1")
            try:
                coefficient = float(text)
            except ValueError as exc:
                raise ModelFormatError(
                    f"reaction {reaction_id!r} has a non-numeric stoichiometry"
                ) from exc
            stoichiometry[species_id] = (
                stoichiometry.get(species_id, 0.0) + sign * coefficient
            )
    kinetic_law = _find_child(element, "kineticLaw")
    if kinetic_law is None:
        raise ModelFormatError(f"reaction {reaction_id!r} has no <kineticLaw>")
    math_element = _find_child(kinetic_law, "math")
    if math_element is None:
        raise ModelFormatError(f"reaction {reaction_id!r} has no <math>")
    return Reaction(
        id=reaction_id,
        name=element.get("name", reaction_id),
        stoichiometry=tuple(stoichiometry.items()),
        rate_law=_from_mathml(math_element),
    )


def _validate_model(model: Model) -> None:
    known = set(species_ids(model))
    for reaction in model.reactions:
        for species_id, _ in reaction.stoichiometry:
            if species_id not in known:
                raise ModelFormatError(
                    f"reaction {reaction.id!r} names unknown species {species_id!r}"
                )
        for name in variables(reaction.rate_law):
            if name not in known and name not in {
                parameter.id for parameter in model.parameters
            }:
                raise ModelFormatError(
                    f"reaction {reaction.id!r} uses unknown name {name!r}"
                )


# ---------------------------------------------------------------------------
# CellML serialisation.  See the module docstring for the supported subset.
# ---------------------------------------------------------------------------
def to_cellml(model: Model) -> str:
    """Serialize ``model`` to a documented CellML 2.0 subset.

    The output holds one ``<model>``, one ``<component>`` with one
    ``<variable>`` per state, parameter and time, and one ``<math>`` with one
    derivative equation per state.  CellML has no reaction element, so the
    writer folds the reactions into the derivative equations.  The real unit
    of each variable rides in a ``nanocad:unit`` attribute, and the CellML
    ``units`` attribute stays ``dimensionless``.  Name and compartment ride in
    ``nanocad:name`` and ``nanocad:compartment``.
    """
    root = ET.Element(
        "model",
        {
            "xmlns": _CELLML_NS,
            "xmlns:nanocad": _CELLML_NANOCAD_NS,
            "name": model.name or model.id,
        },
    )
    component = ET.SubElement(root, "component", {"name": model.id or "model"})
    ET.SubElement(
        component,
        "variable",
        {
            "name": "time",
            "units": "second",
            "nanocad:unit": "s",
        },
    )
    for species in model.species:
        ET.SubElement(
            component,
            "variable",
            {
                "name": species.id,
                "units": "dimensionless",
                "initial_value": repr(species.initial_value_si),
                "nanocad:unit": species.unit,
                "nanocad:name": species.name,
                "nanocad:compartment": species.compartment,
            },
        )
    for parameter in model.parameters:
        ET.SubElement(
            component,
            "variable",
            {
                "name": parameter.id,
                "units": "dimensionless",
                "initial_value": repr(parameter.value_si),
                "nanocad:unit": parameter.unit,
                "nanocad:name": parameter.name,
            },
        )

    math_element = ET.SubElement(component, "math", {"xmlns": _MATHML_NS})
    for species in model.species:
        equation = ET.SubElement(math_element, "apply")
        ET.SubElement(equation, "eq")
        derivative = ET.SubElement(equation, "apply")
        ET.SubElement(derivative, "diff")
        bvar = ET.SubElement(derivative, "bvar")
        ET.SubElement(bvar, "ci").text = "time"
        ET.SubElement(derivative, "ci").text = species.id
        equation.append(_to_mathml(_derivative_expression(model, species.id)))

    ET.indent(root, space="  ")
    body = ET.tostring(root, encoding="unicode")
    return '<?xml version="1.0" encoding="UTF-8"?>\n' + body


def _derivative_expression(model: Model, species_id: str) -> Expr:
    """Fold every reaction term of one species into a single expression."""
    terms: list[Expr] = []
    for reaction in model.reactions:
        for name, coefficient in reaction.stoichiometry:
            if name == species_id:
                terms.append(BinOp("*", Num(coefficient), reaction.rate_law))
    if not terms:
        return Num(0.0)
    result = terms[0]
    for term in terms[1:]:
        result = BinOp("+", result, term)
    return result


def from_cellml(xml: str) -> Model:
    """Parse a CellML 2.0 subset document into a :class:`Model`.

    Raise :class:`XmlParseError` for malformed XML.  Raise
    :class:`ModelFormatError` or :class:`ExpressionError` when the document is
    well-formed but outside the supported subset.

    CellML carries a derivative equation for each state and no reaction.  This
    reader rebuilds one reaction per state whose rate law is that derivative.
    The ODE system is therefore preserved, but the original reaction split is
    not.
    """
    try:
        root = ET.fromstring(xml)
    except ET.ParseError as exc:
        raise XmlParseError(f"malformed XML: {exc}") from exc

    if _local(root.tag) != "model":
        raise ModelFormatError(
            f"expected <model> root element, found <{_local(root.tag)}>"
        )
    components = _children(root, "component")
    if not components:
        raise ModelFormatError("document has no <component> element")

    model_id = root.get("name") or "model"
    name = root.get("name") or model_id

    variables: dict[str, ET.Element] = {}
    for component in components:
        for variable in _children(component, "variable"):
            variable_name = variable.get("name")
            if not variable_name:
                raise ModelFormatError("a <variable> has no name")
            if variable_name in variables:
                raise ModelFormatError(f"duplicate variable {variable_name!r}")
            variables[variable_name] = variable

    derivatives: dict[str, Expr] = {}
    for component in components:
        for math_element in _children(component, "math"):
            for equation in _children(math_element, "apply"):
                parts = list(equation)
                if not parts or _local(parts[0].tag) != "eq":
                    continue
                if len(parts) < 3:
                    raise ExpressionError("an <eq> needs two sides")
                target = _cellml_diff_target(parts[1])
                if target is None:
                    continue
                derivatives[target] = _from_mathml(parts[2])

    species: list[Species] = []
    parameters: list[Parameter] = []
    for variable_name, variable in variables.items():
        if variable_name == "time":
            continue
        if variable_name in derivatives:
            text = variable.get("initial_value")
            if text is None:
                raise ModelFormatError(
                    f"state variable {variable_name!r} has no initial_value"
                )
            try:
                initial_value = float(text)
            except ValueError as exc:
                raise ModelFormatError(
                    f"state variable {variable_name!r} has a non-numeric "
                    f"initial_value {text!r}"
                ) from exc
            species.append(
                Species(
                    id=variable_name,
                    name=_extension_of(variable, "name", variable_name),
                    initial_value_si=initial_value,
                    unit=_unit_of(variable),
                    compartment=_extension_of(variable, "compartment", "environment"),
                )
            )
        elif variable.get("initial_value") is not None:
            text = variable.get("initial_value", "0")
            try:
                value = float(text)
            except ValueError as exc:
                raise ModelFormatError(
                    f"parameter {variable_name!r} has a non-numeric value {text!r}"
                ) from exc
            parameters.append(
                Parameter(
                    id=variable_name,
                    name=_extension_of(variable, "name", variable_name),
                    value_si=value,
                    unit=_unit_of(variable),
                )
            )

    reactions = tuple(
        Reaction(
            id=f"rate_of_{species_state.id}",
            name=f"Rate of {species_state.name}",
            stoichiometry=((species_state.id, 1.0),),
            rate_law=derivatives[species_state.id],
        )
        for species_state in species
    )
    model = Model(
        id=model_id,
        name=name,
        species=tuple(species),
        parameters=tuple(parameters),
        reactions=reactions,
    )
    _validate_model(model)
    return model


def _cellml_diff_target(lhs: ET.Element) -> str | None:
    """Return the state name of a ``diff`` expression, or None."""
    if _local(lhs.tag) != "apply":
        return None
    parts = list(lhs)
    if not parts or _local(parts[0].tag) != "diff":
        return None
    for child in parts[1:]:
        if _local(child.tag) == "ci":
            text = (child.text or "").strip()
            if text:
                return text
    return None


def _extension_of(element: ET.Element, local_name: str, default: str) -> str:
    """Read a namespaced extension attribute such as ``nanocad:name``."""
    for key, value in element.attrib.items():
        if "}" in key and _local(key) == local_name:
            return value
    return default


# ---------------------------------------------------------------------------
# MathML helpers.
# ---------------------------------------------------------------------------
def _to_mathml(expr: Expr) -> ET.Element:
    if isinstance(expr, Num):
        element = ET.Element("cn")
        element.text = repr(expr.value)
        return element
    if isinstance(expr, Var):
        element = ET.Element("ci")
        element.text = expr.name
        return element
    if isinstance(expr, Neg):
        element = ET.Element("apply")
        element.append(ET.Element("minus"))
        element.append(_to_mathml(expr.operand))
        return element
    element = ET.Element("apply")
    element.append(ET.Element(_OP_TO_TAG[expr.op]))
    element.append(_to_mathml(expr.left))
    element.append(_to_mathml(expr.right))
    return element


def _from_mathml(node: ET.Element) -> Expr:
    tag = _local(node.tag)
    if tag == "math":
        children = [child for child in node if _local(child.tag) != "annotation"]
        if not children:
            raise ExpressionError("empty <math> element")
        return _from_mathml(children[0])
    if tag == "cn":
        text = (node.text or "").strip()
        try:
            return Num(float(text))
        except ValueError as exc:
            raise ExpressionError(f"invalid <cn> value {text!r}") from exc
    if tag == "ci":
        name = (node.text or "").strip()
        if not name:
            raise ExpressionError("empty <ci> element")
        return Var(name)
    if tag == "apply":
        children = list(node)
        if not children:
            raise ExpressionError("empty <apply> element")
        operator_tag = _local(children[0].tag)
        arguments = [_from_mathml(child) for child in children[1:]]
        if operator_tag == "minus" and len(arguments) == 1:
            return Neg(arguments[0])
        operator = _TAG_TO_OP.get(operator_tag)
        if operator is None:
            raise ExpressionError(f"unsupported MathML operator <{operator_tag}>")
        if not arguments:
            raise ExpressionError(f"operator <{operator_tag}> has no arguments")
        result = arguments[0]
        for argument in arguments[1:]:
            result = BinOp(operator, result, argument)
        return result
    raise ExpressionError(f"unsupported MathML element <{tag}>")


# ---------------------------------------------------------------------------
# Small XML helpers that ignore namespaces.
# ---------------------------------------------------------------------------
def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _find_child(element: ET.Element, local_name: str) -> ET.Element | None:
    for child in element:
        if _local(child.tag) == local_name:
            return child
    return None


def _children(element: ET.Element, local_name: str) -> list[ET.Element]:
    return [child for child in element if _local(child.tag) == local_name]


def _unit_of(element: ET.Element, default: str = "dimensionless") -> str:
    for key, value in element.attrib.items():
        if _local(key) == "unit":
            return value
    return default
