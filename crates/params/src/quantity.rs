use nanocad_units::{Dimension, Unit};
use serde::de::{self, Deserializer};
use serde::ser::{SerializeStruct, Serializer};
use serde::{Deserialize, Serialize};

use crate::error::ParamError;

/// A physical value in SI base units, plus a presentation unit and an optional
/// uncertainty.
///
/// The value is always stored in SI. The `unit` field is a label. On
/// serialization the value and the uncertainty are written in that unit. This
/// follows the `Quantity` form in `PARAMETERS.md`:
///
/// ```json
/// { "value": 1.05e12, "unit": "Pa", "uncertainty": 0.08e12 }
/// ```
///
/// The uncertainty is [`None`] when it is unknown. It is never zero for an
/// unknown uncertainty, per Rule 3 of `PARAMETERS.md`.
#[derive(Clone, Debug, PartialEq)]
pub struct Quantity {
    value_si: f64,
    unit: String,
    uncertainty_si: Option<f64>,
}

impl Quantity {
    /// Build a quantity from a value in a unit from `nanocad-units`.
    pub fn new(value: f64, unit: Unit) -> Self {
        Self {
            value_si: value * unit.factor_to_si(),
            unit: unit.symbol().to_string(),
            uncertainty_si: None,
        }
    }

    /// Build a quantity from a value that is already in SI base units.
    pub fn from_si(value_si: f64, unit: Unit) -> Self {
        Self {
            value_si,
            unit: unit.symbol().to_string(),
            uncertainty_si: None,
        }
    }

    /// Build a quantity with a derived unit such as `Pa` or `kg*m^2`.
    ///
    /// The value is already in SI base units. The unit must be in the registry.
    pub fn derived(value_si: f64, unit: &str) -> Result<Self, ParamError> {
        if unit_info(unit).is_none() {
            return Err(ParamError::UnknownUnit {
                symbol: unit.to_string(),
            });
        }
        Ok(Self {
            value_si,
            unit: unit.to_string(),
            uncertainty_si: None,
        })
    }

    /// Attach an uncertainty that is already in SI base units.
    pub fn with_uncertainty_si(mut self, uncertainty_si: f64) -> Self {
        self.uncertainty_si = Some(uncertainty_si);
        self
    }

    /// Attach an uncertainty expressed in the presentation unit.
    pub fn with_uncertainty(mut self, uncertainty: f64) -> Self {
        let factor = unit_info(&self.unit)
            .map(|(_, factor)| factor)
            .unwrap_or(1.0);
        self.uncertainty_si = Some(uncertainty * factor);
        self
    }

    /// The value in SI base units.
    pub fn value_si(&self) -> f64 {
        self.value_si
    }

    /// The uncertainty in SI base units, or [`None`] when unknown.
    pub fn uncertainty_si(&self) -> Option<f64> {
        self.uncertainty_si
    }

    /// The presentation unit symbol.
    pub fn unit(&self) -> &str {
        &self.unit
    }

    /// The dimension of the presentation unit, or [`None`] when unknown.
    pub fn dimension(&self) -> Option<Dimension> {
        unit_info(&self.unit).map(|(dimension, _)| dimension)
    }

    /// The value expressed in the presentation unit.
    pub fn value(&self) -> f64 {
        self.value_si / self.factor()
    }

    /// The uncertainty expressed in the presentation unit, or [`None`].
    pub fn uncertainty(&self) -> Option<f64> {
        self.uncertainty_si.map(|u| u / self.factor())
    }

    fn factor(&self) -> f64 {
        unit_info(&self.unit)
            .map(|(_, factor)| factor)
            .unwrap_or(1.0)
    }
}

/// The dimension and the factor-to-SI of a unit symbol.
///
/// The seven SI base units and the accepted non-SI units come from
/// `nanocad-units`. The derived units in `PARAMETERS.md` are listed here,
/// because `nanocad-units` holds base and accepted non-SI units only.
pub(crate) fn unit_info(symbol: &str) -> Option<(Dimension, f64)> {
    if let Some(unit) = known_unit(symbol) {
        return Some((unit.dimension(), unit.factor_to_si()));
    }
    Some(match symbol {
        "1" => (Dimension::new(0, 0, 0, 0, 0, 0, 0), 1.0),
        "Pa" => (Dimension::new(-1, 1, -2, 0, 0, 0, 0), 1.0),
        "kg*m^2" => (Dimension::new(2, 1, 0, 0, 0, 0, 0), 1.0),
        "W/(m*K)" => (Dimension::new(1, 1, -3, 0, -1, 0, 0), 1.0),
        "J/(kg*K)" => (Dimension::new(2, 0, -2, 0, -1, 0, 0), 1.0),
        _ => return None,
    })
}

fn known_unit(symbol: &str) -> Option<Unit> {
    Some(match symbol {
        "m" => Unit::Metre,
        "kg" => Unit::Kilogram,
        "s" => Unit::Second,
        "A" => Unit::Ampere,
        "K" => Unit::Kelvin,
        "mol" => Unit::Mole,
        "cd" => Unit::Candela,
        "Å" => Unit::Angstrom,
        "nm" => Unit::Nanometre,
        "kcal/mol" => Unit::KilocaloriePerMole,
        "fs" => Unit::Femtosecond,
        "ps" => Unit::Picosecond,
        "u" => Unit::AtomicMassUnit,
        "e" => Unit::ElectronCharge,
        _ => return None,
    })
}

impl Serialize for Quantity {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let factor = self.factor();
        let mut state = serializer.serialize_struct("Quantity", 3)?;
        state.serialize_field("value", &(self.value_si / factor))?;
        state.serialize_field("unit", &self.unit)?;
        state.serialize_field("uncertainty", &self.uncertainty_si.map(|u| u / factor))?;
        state.end()
    }
}

#[derive(Deserialize)]
struct QuantityRepr {
    value: f64,
    unit: String,
    #[serde(default)]
    uncertainty: Option<f64>,
}

impl<'de> Deserialize<'de> for Quantity {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let repr = QuantityRepr::deserialize(deserializer)?;
        let (_, factor) = unit_info(&repr.unit)
            .ok_or_else(|| de::Error::custom(format!("unknown unit symbol {:?}", repr.unit)))?;
        Ok(Self {
            value_si: repr.value * factor,
            unit: repr.unit,
            uncertainty_si: repr.uncertainty.map(|u| u * factor),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn base_unit_converts_to_si() {
        let q = Quantity::new(2.5, Unit::Nanometre);
        assert_eq!(q.unit(), "nm");
        assert!((q.value_si() - 2.5e-9).abs() < 1.0e-24);
        assert!((q.value() - 2.5).abs() < 1.0e-12);
    }

    #[test]
    fn derived_unit_keeps_si_value() {
        let q = Quantity::derived(1.05e12, "Pa").unwrap();
        assert_eq!(q.unit(), "Pa");
        assert_eq!(q.value_si(), 1.05e12);
        assert_eq!(q.value(), 1.05e12);
        assert_eq!(q.uncertainty_si(), None);
    }

    #[test]
    fn unknown_derived_unit_is_an_error() {
        assert!(matches!(
            Quantity::derived(1.0, "furlong"),
            Err(ParamError::UnknownUnit { .. })
        ));
    }

    #[test]
    fn uncertainty_is_null_when_unknown_and_kept_when_known() {
        let unknown = Quantity::derived(0.05, "1").unwrap();
        assert_eq!(unknown.uncertainty_si(), None);
        let known = Quantity::derived(0.05, "1").unwrap().with_uncertainty(0.02);
        assert_eq!(known.uncertainty_si(), Some(0.02));
        assert_eq!(known.uncertainty(), Some(0.02));
    }

    #[test]
    fn dimensions_match_the_schema_fields() {
        assert_eq!(
            Quantity::derived(1.0, "Pa").unwrap().dimension(),
            Some(Dimension::new(-1, 1, -2, 0, 0, 0, 0))
        );
        assert_eq!(
            Quantity::derived(1.0, "kg*m^2").unwrap().dimension(),
            Some(Dimension::new(2, 1, 0, 0, 0, 0, 0))
        );
        assert_eq!(
            Quantity::derived(1.0, "W/(m*K)").unwrap().dimension(),
            Some(Dimension::new(1, 1, -3, 0, -1, 0, 0))
        );
        assert_eq!(
            Quantity::derived(1.0, "J/(kg*K)").unwrap().dimension(),
            Some(Dimension::new(2, 0, -2, 0, -1, 0, 0))
        );
        assert_eq!(
            Quantity::derived(1.0, "1").unwrap().dimension(),
            Some(Dimension::new(0, 0, 0, 0, 0, 0, 0))
        );
    }

    #[test]
    fn quantity_json_round_trip() {
        let q = Quantity::derived(1.05e12, "Pa")
            .unwrap()
            .with_uncertainty(0.08e12);
        let json = serde_json::to_string(&q).unwrap();
        let back: Quantity = serde_json::from_str(&json).unwrap();
        assert_eq!(back, q);
    }

    #[test]
    fn unknown_unit_in_json_is_rejected() {
        let json = r#"{"value":1.0,"unit":"furlong","uncertainty":null}"#;
        let parsed: Result<Quantity, _> = serde_json::from_str(json);
        assert!(parsed.is_err());
    }
}
