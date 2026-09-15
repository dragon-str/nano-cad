use std::collections::BTreeMap;

use nanocad_units::Unit;

use crate::error::PartError;

/// The specification of one generator parameter.
///
/// A spec carries the name, the presentation unit, the default, and the valid
/// range. An integer spec reports that the value must be a whole number.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ParameterSpec {
    /// The parameter name. It is unique within one generator.
    pub name: &'static str,
    /// The presentation unit. `None` means dimensionless.
    pub unit: Option<Unit>,
    /// The default value, in SI base units or dimensionless.
    pub default: f64,
    /// The inclusive lower bound.
    pub min: f64,
    /// The inclusive upper bound.
    pub max: f64,
    /// Whether the value must be a whole number.
    pub integer: bool,
    /// A short description for the agent-facing schema.
    pub description: &'static str,
}

impl ParameterSpec {
    /// Builds a parameter spec.
    pub const fn new(
        name: &'static str,
        unit: Option<Unit>,
        default: f64,
        min: f64,
        max: f64,
        integer: bool,
        description: &'static str,
    ) -> Self {
        Self {
            name,
            unit,
            default,
            min,
            max,
            integer,
            description,
        }
    }

    /// Checks one value against the range and the integer rule.
    pub fn accepts(&self, value: f64) -> Result<(), PartError> {
        if !(self.min..=self.max).contains(&value) {
            return Err(PartError::ParameterOutOfRange {
                name: self.name.to_owned(),
                value,
                min: self.min,
                max: self.max,
            });
        }
        if self.integer && value.fract() != 0.0 {
            return Err(PartError::NonIntegerParameter {
                name: self.name.to_owned(),
                value,
            });
        }
        Ok(())
    }
}

/// A named set of parameter values.
///
/// The set holds only the values that a caller supplies. [`PartGenerator`]
/// fills the missing entries from the spec defaults during `resolve`.
///
/// [`PartGenerator`]: crate::PartGenerator
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParameterSet {
    values: BTreeMap<String, f64>,
}

impl ParameterSet {
    /// Creates an empty set.
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds one value and returns the set, for chaining.
    pub fn with(mut self, name: impl Into<String>, value: f64) -> Self {
        self.set(name, value);
        self
    }

    /// Inserts or replaces one value.
    pub fn set(&mut self, name: impl Into<String>, value: f64) {
        self.values.insert(name.into(), value);
    }

    /// Returns a value when it is present.
    pub fn get(&self, name: &str) -> Option<f64> {
        self.values.get(name).copied()
    }

    /// Returns a value or a [`PartError::MissingParameter`].
    pub fn require(&self, name: &str) -> Result<f64, PartError> {
        self.get(name)
            .ok_or_else(|| PartError::MissingParameter(name.to_owned()))
    }

    /// Reports whether a name is present.
    pub fn contains(&self, name: &str) -> bool {
        self.values.contains_key(name)
    }

    /// The number of stored values.
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// Reports whether the set is empty.
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// Iterates over the values in name order.
    pub fn iter(&self) -> impl Iterator<Item = (&str, f64)> + '_ {
        self.values
            .iter()
            .map(|(name, value)| (name.as_str(), *value))
    }

    /// Returns the underlying map.
    pub fn as_map(&self) -> &BTreeMap<String, f64> {
        &self.values
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_parameter_set_stores_and_requires_values() {
        let set = ParameterSet::new()
            .with("cells_x", 2.0)
            .with("cells_y", 3.0);
        assert_eq!(set.get("cells_x"), Some(2.0));
        assert_eq!(set.require("cells_y"), Ok(3.0));
        assert_eq!(set.len(), 2);
        assert!(!set.is_empty());
    }

    #[test]
    fn a_missing_parameter_is_an_error_not_a_panic() {
        let set = ParameterSet::new();
        assert_eq!(
            set.require("cells_x"),
            Err(PartError::MissingParameter("cells_x".to_owned()))
        );
        assert_eq!(set.get("cells_x"), None);
    }

    #[test]
    fn a_spec_accepts_a_value_inside_its_range() {
        let spec = ParameterSpec::new("cells_x", None, 1.0, 1.0, 8.0, true, "cells along x");
        assert_eq!(spec.accepts(1.0), Ok(()));
        assert_eq!(spec.accepts(8.0), Ok(()));
        assert_eq!(spec.accepts(4.0), Ok(()));
    }

    #[test]
    fn a_spec_rejects_a_value_outside_its_range() {
        let spec = ParameterSpec::new("cells_x", None, 1.0, 1.0, 8.0, true, "cells along x");
        assert_eq!(
            spec.accepts(0.0),
            Err(PartError::ParameterOutOfRange {
                name: "cells_x".to_owned(),
                value: 0.0,
                min: 1.0,
                max: 8.0,
            })
        );
        assert!(spec.accepts(9.0).is_err());
    }

    #[test]
    fn an_integer_spec_rejects_a_fraction() {
        let spec = ParameterSpec::new("cells_x", None, 1.0, 1.0, 8.0, true, "cells along x");
        assert_eq!(
            spec.accepts(2.5),
            Err(PartError::NonIntegerParameter {
                name: "cells_x".to_owned(),
                value: 2.5,
            })
        );
    }

    #[test]
    fn iteration_visits_values_in_name_order() {
        let set = ParameterSet::new().with("b", 2.0).with("a", 1.0);
        let names: Vec<&str> = set.iter().map(|(name, _)| name).collect();
        assert_eq!(names, vec!["a", "b"]);
    }
}
