use nanocad_model::Part;

use crate::error::PartError;
use crate::parameter::{ParameterSet, ParameterSpec};

/// A description of a generated part.
///
/// The schema names the generator, the parameters that produced the part, and
/// the resulting counts. It is the compact record that a caller or an agent
/// stores beside the geometry.
#[derive(Clone, Debug, PartialEq)]
pub struct PartSchema {
    /// The part name.
    pub name: String,
    /// The generator id.
    pub generator_id: String,
    /// The resolved parameters that produced the part.
    pub parameters: ParameterSet,
    /// The number of atoms.
    pub atom_count: usize,
    /// The number of bonds.
    pub bond_count: usize,
    /// The material label.
    pub material: String,
}

impl PartSchema {
    /// The schema version for the part description.
    pub const SCHEMA_VERSION: u32 = 1;

    /// Builds a schema from a generator id, the resolved parameters, and the
    /// generated part.
    pub fn from_part(
        generator_id: impl Into<String>,
        parameters: &ParameterSet,
        part: &Part,
    ) -> Self {
        Self {
            name: part.name.clone(),
            generator_id: generator_id.into(),
            parameters: parameters.clone(),
            atom_count: part.atom_count(),
            bond_count: part.bond_count(),
            material: part.material.clone(),
        }
    }
}

/// A parametric part generator.
///
/// Every generator declares its name, its parameter specs, and a `generate`
/// method. `generate` validates its inputs and returns a [`Part`] or a
/// [`PartError`]. No path panics.
pub trait PartGenerator {
    /// The stable generator id, for example `"diamond"`.
    fn id(&self) -> &'static str;

    /// The human-readable generator name.
    fn name(&self) -> &'static str;

    /// The parameter specs, in a stable order.
    fn parameters(&self) -> &'static [ParameterSpec];

    /// Validates the inputs and builds the part.
    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError>;

    /// Fills the defaults, rejects unknown names, and checks every range.
    ///
    /// A generator calls this first, then reads the resolved values.
    fn resolve(&self, supplied: &ParameterSet) -> Result<ParameterSet, PartError> {
        for (name, _) in supplied.iter() {
            if !self.parameters().iter().any(|spec| spec.name == name) {
                return Err(PartError::UnknownParameter(name.to_owned()));
            }
        }
        let mut resolved = ParameterSet::new();
        for spec in self.parameters() {
            let value = supplied.get(spec.name).unwrap_or(spec.default);
            spec.accepts(value)?;
            resolved.set(spec.name, value);
        }
        Ok(resolved)
    }

    /// Builds the part with every default parameter.
    fn generate_with_defaults(&self) -> Result<Part, PartError> {
        self.generate(&ParameterSet::new())
    }

    /// Builds the part schema for a result of [`PartGenerator::generate`].
    fn schema(&self, parameters: &ParameterSet, part: &Part) -> PartSchema {
        PartSchema::from_part(self.id(), parameters, part)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nanocad_model::Topology;
    use nanocad_units::Unit;

    const SPECS: &[ParameterSpec] = &[ParameterSpec::new(
        "count",
        None,
        2.0,
        1.0,
        4.0,
        true,
        "a test count",
    )];

    struct Stub;

    impl PartGenerator for Stub {
        fn id(&self) -> &'static str {
            "stub"
        }

        fn name(&self) -> &'static str {
            "Stub generator"
        }

        fn parameters(&self) -> &'static [ParameterSpec] {
            SPECS
        }

        fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
            let resolved = self.resolve(parameters)?;
            let count = resolved.require("count")? as usize;
            let mut topology = Topology::new();
            for index in 0..count {
                topology.add_atom(nanocad_model::Atom::new(
                    nanocad_model::Element::CARBON,
                    [index as f64 * 1.0e-10, 0.0, 0.0],
                    0.0,
                    "C",
                ));
            }
            Ok(Part::new("stub", topology).with_material("carbon"))
        }
    }

    #[test]
    fn resolve_fills_the_default_when_none_is_supplied() {
        let resolved = Stub.resolve(&ParameterSet::new()).expect("valid");
        assert_eq!(resolved.get("count"), Some(2.0));
    }

    #[test]
    fn resolve_applies_a_supplied_value() {
        let resolved = Stub
            .resolve(&ParameterSet::new().with("count", 3.0))
            .expect("valid");
        assert_eq!(resolved.get("count"), Some(3.0));
    }

    #[test]
    fn resolve_rejects_an_unknown_name() {
        assert_eq!(
            Stub.resolve(&ParameterSet::new().with("nope", 1.0)),
            Err(PartError::UnknownParameter("nope".to_owned()))
        );
    }

    #[test]
    fn resolve_rejects_an_out_of_range_value() {
        assert!(Stub
            .resolve(&ParameterSet::new().with("count", 0.0))
            .is_err());
    }

    #[test]
    fn generate_with_defaults_and_schema_agree() {
        let part = Stub.generate_with_defaults().expect("valid");
        assert_eq!(part.atom_count(), 2);
        let schema = Stub.schema(&ParameterSet::new().with("count", 2.0), &part);
        assert_eq!(schema.generator_id, "stub");
        assert_eq!(schema.name, "stub");
        assert_eq!(schema.atom_count, 2);
        assert_eq!(schema.bond_count, 0);
        assert_eq!(schema.material, "carbon");
        assert_eq!(PartSchema::SCHEMA_VERSION, 1);
        assert_eq!(SPECS[0].unit, None);
        assert_eq!(
            ParameterSpec::new("x", Some(Unit::Metre), 1.0, 0.0, 2.0, false, "x").unit,
            Some(Unit::Metre)
        );
    }
}
