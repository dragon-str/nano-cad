//! The `nanocad._core` extension module.
//!
//! Every tool named in `ARCHITECTURE.md` is registered here. Each function
//! converts its errors to a Python exception; no call panics.

mod assembly;
mod device;
mod engine;
mod formats;
mod jigs;
mod model;
mod params;
mod parts;
mod units;
mod util;

use pyo3::prelude::*;

/// The crate version string.
#[pyfunction]
fn version() -> &'static str {
    crate::version()
}

#[pymodule]
fn _core(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(version, module)?)?;

    // Units.
    module.add_class::<units::Quantity>()?;
    module.add_function(wrap_pyfunction!(units::convert_units, module)?)?;
    module.add_function(wrap_pyfunction!(units::to_si, module)?)?;
    module.add_function(wrap_pyfunction!(units::from_si, module)?)?;
    module.add_function(wrap_pyfunction!(units::unit_symbols, module)?)?;
    module.add_function(wrap_pyfunction!(units::unit_dimension, module)?)?;

    // Model.
    module.add_class::<model::Atom>()?;
    module.add_class::<model::Topology>()?;
    module.add_class::<model::Part>()?;
    module.add_class::<model::Document>()?;
    module.add_function(wrap_pyfunction!(model::document_schema_version, module)?)?;

    // Formats.
    module.add_function(wrap_pyfunction!(formats::read_ncz, module)?)?;
    module.add_function(wrap_pyfunction!(formats::write_ncz, module)?)?;
    module.add_function(wrap_pyfunction!(formats::ncz_to_bytes, module)?)?;
    module.add_function(wrap_pyfunction!(formats::ncz_from_bytes, module)?)?;
    module.add_function(wrap_pyfunction!(formats::import_xyz, module)?)?;
    module.add_function(wrap_pyfunction!(formats::export_xyz, module)?)?;
    module.add_function(wrap_pyfunction!(formats::import_mmp, module)?)?;
    module.add_function(wrap_pyfunction!(formats::export_mmp, module)?)?;
    module.add_function(wrap_pyfunction!(formats::save, module)?)?;
    module.add_function(wrap_pyfunction!(formats::load, module)?)?;
    module.add_function(wrap_pyfunction!(formats::export, module)?)?;
    module.add_function(wrap_pyfunction!(formats::list_parts, module)?)?;

    // Parts.
    module.add_class::<parts::ParameterInfo>()?;
    module.add_class::<parts::GeneratorInfo>()?;
    module.add_class::<parts::PlanetaryGeometry>()?;
    module.add_function(wrap_pyfunction!(parts::list_generators, module)?)?;
    module.add_function(wrap_pyfunction!(parts::generator_parameters, module)?)?;
    module.add_function(wrap_pyfunction!(parts::generator_schema, module)?)?;
    module.add_function(wrap_pyfunction!(parts::generate_part, module)?)?;
    module.add_function(wrap_pyfunction!(parts::create_part, module)?)?;
    module.add_function(wrap_pyfunction!(parts::generate_gear, module)?)?;
    module.add_function(wrap_pyfunction!(parts::generate_nanotube, module)?)?;
    module.add_function(wrap_pyfunction!(parts::generate_lattice, module)?)?;
    module.add_function(wrap_pyfunction!(parts::assemble, module)?)?;
    module.add_function(wrap_pyfunction!(parts::measure, module)?)?;
    module.add_function(wrap_pyfunction!(parts::planetary_geometry, module)?)?;
    module.add_function(wrap_pyfunction!(parts::validate_part, module)?)?;
    module.add_function(wrap_pyfunction!(parts::element_valence, module)?)?;

    // Parameters.
    module.add_class::<params::ParamQuantity>()?;
    module.add_class::<params::Provenance>()?;
    module.add_class::<params::PartRecord>()?;
    module.add_class::<params::ParameterLibrary>()?;
    module.add_function(wrap_pyfunction!(params::check_record, module)?)?;
    module.add_function(wrap_pyfunction!(params::is_consistent, module)?)?;

    // Engine.
    module.add_class::<engine::System>()?;
    module.add_class::<engine::MinimizeResult>()?;
    module.add_class::<engine::VelocityVerlet>()?;
    module.add_function(wrap_pyfunction!(engine::minimize, module)?)?;
    module.add_function(wrap_pyfunction!(engine::relax, module)?)?;
    module.add_function(wrap_pyfunction!(engine::integrate, module)?)?;
    module.add_function(wrap_pyfunction!(engine::run_md, module)?)?;
    module.add_function(wrap_pyfunction!(engine::coulomb_constant, module)?)?;

    // Jigs.
    module.add_class::<jigs::AnchorJig>()?;
    module.add_class::<jigs::SpringJig>()?;

    // Device layer.
    module.add_class::<device::Quat>()?;
    module.add_class::<device::RigidBody>()?;
    module.add_class::<device::RigidBodySystem>()?;
    module.add_function(wrap_pyfunction!(device::add_constraint, module)?)?;

    // Planetary device assembly and URDF export.
    module.add_class::<assembly::PlanetaryAssembly>()?;
    module.add_function(wrap_pyfunction!(assembly::assemble_planetary, module)?)?;

    Ok(())
}
