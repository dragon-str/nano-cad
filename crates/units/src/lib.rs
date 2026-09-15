//! SI units, `Quantity`, and conversions.
//!
//! Every value inside a [`Quantity`] is stored in SI base units. The stored
//! [`Unit`] is a presentation label. This follows ADR-0003: SI internally, and
//! conversion at the boundary.
#![forbid(unsafe_code)]

use core::fmt;

/// Exponents of the seven SI base dimensions.
///
/// The order is metre, kilogram, second, ampere, kelvin, mole, candela.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Dimension {
    pub length: i8,
    pub mass: i8,
    pub time: i8,
    pub current: i8,
    pub temperature: i8,
    pub amount: i8,
    pub luminous_intensity: i8,
}

impl Dimension {
    /// Build a dimension from its seven exponents.
    pub const fn new(
        length: i8,
        mass: i8,
        time: i8,
        current: i8,
        temperature: i8,
        amount: i8,
        luminous_intensity: i8,
    ) -> Self {
        Self {
            length,
            mass,
            time,
            current,
            temperature,
            amount,
            luminous_intensity,
        }
    }
}

/// A unit of measure. The factor converts a value in this unit to SI.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Unit {
    // SI base units.
    Metre,
    Kilogram,
    Second,
    Ampere,
    Kelvin,
    Mole,
    Candela,
    // Accepted non-SI units.
    Angstrom,
    Nanometre,
    KilocaloriePerMole,
    Femtosecond,
    Picosecond,
    AtomicMassUnit,
    ElectronCharge,
}

impl Unit {
    /// The dimension of this unit.
    pub const fn dimension(self) -> Dimension {
        match self {
            Unit::Metre | Unit::Angstrom | Unit::Nanometre => Dimension::new(1, 0, 0, 0, 0, 0, 0),
            Unit::Kilogram | Unit::AtomicMassUnit => Dimension::new(0, 1, 0, 0, 0, 0, 0),
            Unit::Second | Unit::Femtosecond | Unit::Picosecond => {
                Dimension::new(0, 0, 1, 0, 0, 0, 0)
            }
            Unit::Ampere => Dimension::new(0, 0, 0, 1, 0, 0, 0),
            Unit::Kelvin => Dimension::new(0, 0, 0, 0, 1, 0, 0),
            Unit::Mole => Dimension::new(0, 0, 0, 0, 0, 1, 0),
            Unit::Candela => Dimension::new(0, 0, 0, 0, 0, 0, 1),
            Unit::ElectronCharge => Dimension::new(0, 0, 1, 1, 0, 0, 0),
            Unit::KilocaloriePerMole => Dimension::new(2, 1, -2, 0, 0, -1, 0),
        }
    }

    /// The multiplier that converts a value in this unit to SI.
    pub const fn factor_to_si(self) -> f64 {
        match self {
            Unit::Metre => 1.0,
            Unit::Kilogram => 1.0,
            Unit::Second => 1.0,
            Unit::Ampere => 1.0,
            Unit::Kelvin => 1.0,
            Unit::Mole => 1.0,
            Unit::Candela => 1.0,
            Unit::Angstrom => 1.0e-10,
            Unit::Nanometre => 1.0e-9,
            // Thermochemical kilocalorie, exact by definition: 4184 J/mol.
            Unit::KilocaloriePerMole => 4184.0,
            Unit::Femtosecond => 1.0e-15,
            Unit::Picosecond => 1.0e-12,
            // CODATA 2018 unified atomic mass unit.
            Unit::AtomicMassUnit => 1.660_539_066_60e-27,
            // SI 2019 elementary charge, exact by definition.
            Unit::ElectronCharge => 1.602_176_634e-19,
        }
    }

    /// The short symbol of this unit.
    pub const fn symbol(self) -> &'static str {
        match self {
            Unit::Metre => "m",
            Unit::Kilogram => "kg",
            Unit::Second => "s",
            Unit::Ampere => "A",
            Unit::Kelvin => "K",
            Unit::Mole => "mol",
            Unit::Candela => "cd",
            Unit::Angstrom => "Å",
            Unit::Nanometre => "nm",
            Unit::KilocaloriePerMole => "kcal/mol",
            Unit::Femtosecond => "fs",
            Unit::Picosecond => "ps",
            Unit::AtomicMassUnit => "u",
            Unit::ElectronCharge => "e",
        }
    }

    /// True when both units measure the same dimension.
    pub const fn is_compatible_with(self, other: Unit) -> bool {
        dimensions_equal(self.dimension(), other.dimension())
    }
}

const fn dimensions_equal(a: Dimension, b: Dimension) -> bool {
    a.length == b.length
        && a.mass == b.mass
        && a.time == b.time
        && a.current == b.current
        && a.temperature == b.temperature
        && a.amount == b.amount
        && a.luminous_intensity == b.luminous_intensity
}

impl fmt::Display for Unit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.symbol())
    }
}

/// An error from a unit conversion.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum UnitError {
    /// The two units measure different dimensions.
    DimensionMismatch { from: Unit, to: Unit },
}

impl fmt::Display for UnitError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            UnitError::DimensionMismatch { from, to } => {
                write!(f, "cannot convert from {from} to {to}: dimension mismatch")
            }
        }
    }
}

impl std::error::Error for UnitError {}

/// A physical value. The internal value is always in SI base units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Quantity {
    value_si: f64,
    unit: Unit,
}

impl Quantity {
    /// Build a quantity from a value in `unit`. The value is converted to SI.
    pub fn new(value: f64, unit: Unit) -> Self {
        Self {
            value_si: value * unit.factor_to_si(),
            unit,
        }
    }

    /// Build a quantity from a value already in SI base units.
    pub fn from_si(value_si: f64, unit: Unit) -> Self {
        Self { value_si, unit }
    }

    /// The value in SI base units.
    pub fn value_si(&self) -> f64 {
        self.value_si
    }

    /// The presentation unit.
    pub fn unit(&self) -> Unit {
        self.unit
    }

    /// Read the value expressed in `target`.
    ///
    /// Return an error when `target` measures a different dimension.
    pub fn value_in(&self, target: Unit) -> Result<f64, UnitError> {
        if !self.unit.is_compatible_with(target) {
            return Err(UnitError::DimensionMismatch {
                from: self.unit,
                to: target,
            });
        }
        Ok(self.value_si / target.factor_to_si())
    }

    /// Return the same physical value with `target` as the presentation unit.
    pub fn convert_to(&self, target: Unit) -> Result<Quantity, UnitError> {
        if !self.unit.is_compatible_with(target) {
            return Err(UnitError::DimensionMismatch {
                from: self.unit,
                to: target,
            });
        }
        Ok(Quantity {
            value_si: self.value_si,
            unit: target,
        })
    }
}

impl fmt::Display for Quantity {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let shown = self.value_si / self.unit.factor_to_si();
        write!(f, "{shown} {}", self.unit)
    }
}

/// Returns the crate version string.
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(got: f64, want: f64) {
        let tol = 1.0e-12 * want.abs().max(1.0);
        assert!(
            (got - want).abs() <= tol,
            "got {got:?}, want {want:?}, tolerance {tol:?}"
        );
    }

    fn round_trip(value: f64, unit: Unit) {
        let q = Quantity::new(value, unit);
        assert_close(q.value_in(unit).unwrap(), value);
        assert_close(q.value_si(), value * unit.factor_to_si());
    }

    #[test]
    fn version_is_not_empty() {
        assert!(!version().is_empty());
    }

    #[test]
    fn si_base_units_round_trip() {
        round_trip(3.5, Unit::Metre);
        round_trip(2.0, Unit::Kilogram);
        round_trip(7.25, Unit::Second);
        round_trip(0.5, Unit::Ampere);
        round_trip(300.0, Unit::Kelvin);
        round_trip(4.0, Unit::Mole);
        round_trip(1.0, Unit::Candela);
    }

    #[test]
    fn angstrom_round_trip() {
        round_trip(5.0, Unit::Angstrom);
        round_trip(1.234, Unit::Angstrom);
    }

    #[test]
    fn nanometre_round_trip() {
        round_trip(5.0, Unit::Nanometre);
        round_trip(0.75, Unit::Nanometre);
    }

    #[test]
    fn femtosecond_round_trip() {
        round_trip(200.0, Unit::Femtosecond);
        round_trip(1.0e6, Unit::Femtosecond);
    }

    #[test]
    fn picosecond_round_trip() {
        round_trip(200.0, Unit::Picosecond);
        round_trip(0.125, Unit::Picosecond);
    }

    #[test]
    fn atomic_mass_unit_round_trip() {
        round_trip(12.0, Unit::AtomicMassUnit);
        round_trip(1.007_276, Unit::AtomicMassUnit);
    }

    #[test]
    fn electron_charge_round_trip() {
        round_trip(1.0, Unit::ElectronCharge);
        round_trip(-2.0, Unit::ElectronCharge);
    }

    #[test]
    fn kilocalorie_per_mole_round_trip() {
        round_trip(1.0, Unit::KilocaloriePerMole);
        round_trip(0.239, Unit::KilocaloriePerMole);
    }

    #[test]
    fn angstrom_to_nanometre() {
        let q = Quantity::new(10.0, Unit::Angstrom);
        assert_close(q.value_in(Unit::Nanometre).unwrap(), 1.0);
    }

    #[test]
    fn si_values_are_exact_for_scaled_units() {
        assert_close(Quantity::new(1.0, Unit::Angstrom).value_si(), 1.0e-10);
        assert_close(Quantity::new(1.0, Unit::Nanometre).value_si(), 1.0e-9);
        assert_close(Quantity::new(1.0, Unit::Femtosecond).value_si(), 1.0e-15);
        assert_close(Quantity::new(1.0, Unit::Picosecond).value_si(), 1.0e-12);
        assert_close(
            Quantity::new(1.0, Unit::ElectronCharge).value_si(),
            1.602_176_634e-19,
        );
        assert_close(
            Quantity::new(1.0, Unit::AtomicMassUnit).value_si(),
            1.660_539_066_60e-27,
        );
        assert_close(
            Quantity::new(1.0, Unit::KilocaloriePerMole).value_si(),
            4184.0,
        );
    }

    #[test]
    fn convert_to_preserves_si_value() {
        let q = Quantity::new(2.5, Unit::Angstrom);
        let converted = q.convert_to(Unit::Nanometre).unwrap();
        assert_close(converted.value_si(), q.value_si());
        assert_eq!(converted.unit(), Unit::Nanometre);
        assert_close(converted.value_in(Unit::Angstrom).unwrap(), 2.5);
    }

    #[test]
    fn from_si_is_complete() {
        let q = Quantity::from_si(1.0e-10, Unit::Angstrom);
        assert_close(q.value_in(Unit::Angstrom).unwrap(), 1.0);
    }

    #[test]
    fn dimension_mismatch_is_an_error() {
        let q = Quantity::new(1.0, Unit::Angstrom);
        assert_eq!(
            q.value_in(Unit::Second),
            Err(UnitError::DimensionMismatch {
                from: Unit::Angstrom,
                to: Unit::Second,
            })
        );
        assert_eq!(
            q.convert_to(Unit::Second),
            Err(UnitError::DimensionMismatch {
                from: Unit::Angstrom,
                to: Unit::Second,
            })
        );
    }

    #[test]
    fn electron_charge_has_charge_dimension() {
        assert_eq!(
            Unit::ElectronCharge.dimension(),
            Dimension::new(0, 0, 1, 1, 0, 0, 0)
        );
        assert!(!Unit::ElectronCharge.is_compatible_with(Unit::Ampere));
    }

    #[test]
    fn kilocalorie_per_mole_dimension() {
        assert_eq!(
            Unit::KilocaloriePerMole.dimension(),
            Dimension::new(2, 1, -2, 0, 0, -1, 0)
        );
    }

    #[test]
    fn unit_symbols_and_display() {
        assert_eq!(Unit::Metre.symbol(), "m");
        assert_eq!(format!("{}", Unit::KilocaloriePerMole), "kcal/mol");
        assert_eq!(format!("{}", Quantity::new(2.0, Unit::Nanometre)), "2 nm");
    }
}
