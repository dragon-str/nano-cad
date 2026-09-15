"""Tests for the L4-01 lumped adapter and its SBML round trip."""

from __future__ import annotations

import xml.etree.ElementTree as ET

import pytest

from nanocad import lumped_adapter as lumped

PRESSURE_SOURCE_INITIAL_PA = 1.0e6
PRESSURE_SINK_INITIAL_PA = 1.0e5
CAPACITANCE_SOURCE_M3_PER_PA = 1.0e-12
CAPACITANCE_SINK_M3_PER_PA = 5.0e-12
RESISTANCE_PA_S_PER_M3 = 2.0e14


def _model() -> lumped.Model:
    return lumped.respirocyte_gas_transport(
        pressure_source_initial_pa=PRESSURE_SOURCE_INITIAL_PA,
        pressure_sink_initial_pa=PRESSURE_SINK_INITIAL_PA,
        capacitance_source_m3_per_pa=CAPACITANCE_SOURCE_M3_PER_PA,
        capacitance_sink_m3_per_pa=CAPACITANCE_SINK_M3_PER_PA,
        resistance_pa_s_per_m3=RESISTANCE_PA_S_PER_M3,
    )


def _analytic(time_s: float) -> tuple[float, float]:
    return lumped.two_compartment_pressures(
        PRESSURE_SOURCE_INITIAL_PA,
        PRESSURE_SINK_INITIAL_PA,
        CAPACITANCE_SOURCE_M3_PER_PA,
        CAPACITANCE_SINK_M3_PER_PA,
        RESISTANCE_PA_S_PER_M3,
        time_s,
    )


def test_rk4_matches_two_compartment_analytic() -> None:
    result = lumped.simulate(_model(), step_s=1.0, step_count=500, analytic=_analytic)
    assert result.method == "rk4-fixed-step"
    assert result.validation == "cross-checked"
    assert result.relative_error is not None
    assert result.relative_error < 1.0e-9


def test_simulate_without_analytic_is_unverified() -> None:
    result = lumped.simulate(_model(), step_s=1.0, step_count=10)
    assert result.validation == "unverified"
    assert result.relative_error is None


def test_total_volume_is_conserved() -> None:
    result = lumped.simulate(_model(), step_s=1.0, step_count=500)
    initial = (
        CAPACITANCE_SOURCE_M3_PER_PA * PRESSURE_SOURCE_INITIAL_PA
        + CAPACITANCE_SINK_M3_PER_PA * PRESSURE_SINK_INITIAL_PA
    )
    final_source_pa, final_sink_pa = result.states[-1]
    final = (
        CAPACITANCE_SOURCE_M3_PER_PA * final_source_pa
        + CAPACITANCE_SINK_M3_PER_PA * final_sink_pa
    )
    assert final == pytest.approx(initial, rel=1.0e-9)


def test_pressures_relax_to_weighted_mean() -> None:
    equilibrium_pa = (
        CAPACITANCE_SOURCE_M3_PER_PA * PRESSURE_SOURCE_INITIAL_PA
        + CAPACITANCE_SINK_M3_PER_PA * PRESSURE_SINK_INITIAL_PA
    ) / (CAPACITANCE_SOURCE_M3_PER_PA + CAPACITANCE_SINK_M3_PER_PA)
    result = lumped.simulate(_model(), step_s=1.0, step_count=4000)
    final_source_pa, final_sink_pa = result.states[-1]
    assert final_source_pa == pytest.approx(equilibrium_pa, rel=1.0e-6)
    assert final_sink_pa == pytest.approx(equilibrium_pa, rel=1.0e-6)


def test_sbml_round_trip_preserves_model() -> None:
    model = _model()
    parsed = lumped.from_sbml(lumped.to_sbml(model))
    assert parsed == model
    assert parsed.species[0].initial_value_si == PRESSURE_SOURCE_INITIAL_PA
    assert parsed.parameters[0].value_si == RESISTANCE_PA_S_PER_M3
    assert parsed.reactions[0].rate_law.to_infix() == (
        model.reactions[0].rate_law.to_infix()
    )


def test_sbml_document_is_well_formed() -> None:
    document = lumped.to_sbml(_model())
    root = ET.fromstring(document)
    assert root.tag.rsplit("}", 1)[-1] == "sbml"
    assert document.startswith("<?xml")


def test_from_sbml_malformed_xml_raises_typed_error() -> None:
    with pytest.raises(lumped.XmlParseError):
        lumped.from_sbml("<sbml><model></sbml>")
    with pytest.raises(lumped.XmlParseError):
        lumped.from_sbml("not xml at all <<<")


def test_from_sbml_missing_model_raises_format_error() -> None:
    with pytest.raises(lumped.ModelFormatError):
        lumped.from_sbml("<sbml><foo/></sbml>")
    with pytest.raises(lumped.ModelFormatError):
        lumped.from_sbml("<html/>")


def test_from_sbml_unknown_species_raises_format_error() -> None:
    document = lumped.to_sbml(_model()).replace(
        'species="pressure_sink_pa"', 'species="ghost"'
    )
    with pytest.raises(lumped.ModelFormatError):
        lumped.from_sbml(document)


def test_from_sbml_unsupported_mathml_raises_expression_error() -> None:
    document = lumped.to_sbml(_model()).replace("<divide />", "<log />")
    with pytest.raises(lumped.ExpressionError):
        lumped.from_sbml(document)


def test_cellml_round_trip_preserves_the_dynamics() -> None:
    model = _model()
    parsed = lumped.from_cellml(lumped.to_cellml(model))
    assert [item.id for item in parsed.species] == [item.id for item in model.species]
    for original, restored in zip(model.species, parsed.species, strict=True):
        assert restored.initial_value_si == original.initial_value_si
        assert restored.unit == original.unit
        assert restored.name == original.name
        assert restored.compartment == original.compartment
    for original, restored in zip(model.parameters, parsed.parameters, strict=True):
        assert restored.id == original.id
        assert restored.value_si == original.value_si
        assert restored.unit == original.unit
    state = lumped.initial_state(model)
    assert lumped.derivatives(parsed, state) == pytest.approx(
        lumped.derivatives(model, state)
    )


def test_cellml_document_is_well_formed() -> None:
    document = lumped.to_cellml(_model())
    root = ET.fromstring(document)
    assert root.tag.rsplit("}", 1)[-1] == "model"
    assert document.startswith("<?xml")


def test_from_cellml_malformed_xml_raises_typed_error() -> None:
    with pytest.raises(lumped.XmlParseError):
        lumped.from_cellml("<model><component></model>")
    with pytest.raises(lumped.XmlParseError):
        lumped.from_cellml("not xml at all <<<")


def test_from_cellml_missing_component_raises_format_error() -> None:
    with pytest.raises(lumped.ModelFormatError):
        lumped.from_cellml(
            '<model xmlns="http://www.cellml.org/cellml/2.0#" name="m"/>'
        )
    with pytest.raises(lumped.ModelFormatError):
        lumped.from_cellml("<html/>")


def test_from_cellml_unknown_math_operator_raises_expression_error() -> None:
    document = lumped.to_cellml(_model()).replace("<divide />", "<log />")
    with pytest.raises(lumped.ExpressionError):
        lumped.from_cellml(document)


def test_from_cellml_unknown_name_raises_format_error() -> None:
    document = lumped.to_cellml(_model()).replace(
        "<ci>resistance_pa_s_per_m3</ci>", "<ci>ghost</ci>"
    )
    with pytest.raises(lumped.ModelFormatError):
        lumped.from_cellml(document)


def test_from_cellml_state_without_initial_value_raises() -> None:
    xml = (
        '<model xmlns="http://www.cellml.org/cellml/2.0#" name="m">'
        '<component name="c">'
        '<variable name="time" units="second"/>'
        '<variable name="x" units="dimensionless"/>'
        '<variable name="k" units="dimensionless" initial_value="1"/>'
        '<math xmlns="http://www.w3.org/1998/Math/MathML">'
        "<apply><eq/><apply><diff/><bvar><ci>time</ci></bvar><ci>x</ci></apply>"
        "<apply><times/><ci>k</ci><ci>x</ci></apply></apply>"
        "</math></component></model>"
    )
    with pytest.raises(lumped.ModelFormatError):
        lumped.from_cellml(xml)


def test_expression_parser_rejects_bad_text() -> None:
    with pytest.raises(lumped.ExpressionError):
        lumped.parse_expression("1 +")
    with pytest.raises(lumped.ExpressionError):
        lumped.parse_expression("(1 + 2")
    with pytest.raises(lumped.ExpressionError):
        lumped.parse_expression("")


def test_evaluate_unknown_name_raises_typed_error() -> None:
    expression = lumped.parse_expression("2 * missing")
    with pytest.raises(lumped.ExpressionError):
        lumped.evaluate(expression, {"present": 1.0})
