use nanocad_params::{Method, PartRecord, Provenance, Quantity, Validation, SCHEMA, VERSION};

/// Build a quantity for a registry unit. An absent uncertainty stays `None`.
pub fn quantity(value_si: f64, unit: &str, uncertainty_si: Option<f64>) -> Quantity {
    let quantity = Quantity::derived(value_si, unit).expect("the unit must be in the registry");
    match uncertainty_si {
        Some(uncertainty) => quantity.with_uncertainty_si(uncertainty),
        None => quantity,
    }
}

/// A complete, physically plausible record with known values.
pub fn good_record() -> PartRecord {
    PartRecord {
        schema: SCHEMA.to_string(),
        version: VERSION,
        part_id: "gear.sun.diamond.2nm.12t".to_string(),
        geometry_ref: "sha256:deadbeef".to_string(),
        material: "diamond".to_string(),
        atoms: 1234,
        mass_kg: quantity(1.5e-21, "kg", Some(0.05e-21)),
        inertia_kg_m2: [
            quantity(1.0e-40, "kg*m^2", Some(1.0e-42)),
            quantity(1.1e-40, "kg*m^2", Some(1.0e-42)),
            quantity(1.2e-40, "kg*m^2", None),
        ],
        elastic_modulus_pa: quantity(1.05e12, "Pa", Some(0.08e12)),
        shear_modulus_pa: quantity(4.2e11, "Pa", Some(0.03e11)),
        poisson_ratio: quantity(0.2, "1", Some(0.01)),
        failure_stress_pa: quantity(2.0e10, "Pa", Some(0.5e10)),
        friction_coefficient: quantity(0.05, "1", Some(0.02)),
        thermal_conductivity_w_m_k: quantity(1000.0, "W/(m*K)", Some(50.0)),
        specific_heat_j_kg_k: quantity(500.0, "J/(kg*K)", Some(20.0)),
        method: Method::Md,
        validation: Validation::Unverified,
        provenance: Provenance {
            source: "nanocad-engine".to_string(),
            method: Method::Md,
            code_version: "0.1.0".to_string(),
            force_field: Some("MM4".to_string()),
            timestamp: "2026-09-13T00:00:00Z".to_string(),
            uncertainty: None,
            validation: Validation::Unverified,
            notes: "synthetic test record; not a measurement".to_string(),
        },
    }
}
