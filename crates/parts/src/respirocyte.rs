//! Respirocyte machine parts: rotor and bearing, pump, and gas tank.
//!
//! A respirocyte is a hypothetical diamondoid blood-borne oxygen carrier. This
//! module builds three skeletal subsystems from a cut diamond cubic lattice:
//! a rotating disk with a concentric bearing sleeve, a cylinder and piston
//! pair, and a hollow spherical pressure shell. The atoms sit on the diamond
//! lattice, so every bond is the diamond first-shell distance. The generators
//! follow [`PartGenerator`] and never panic.
//!
//! Source of the respirocyte concept: R. A. Freitas Jr., "Exploratory Design in
//! Medical Nanorobotics: The Respirocyte", Artificial Cells, Blood Substitutes,
//! and Immobilization Biotechnology 26(4), 411 (1998). This code simulates a
//! skeletal geometry only. It is not a validated medical device.

use std::collections::HashMap;

use nanocad_model::{Atom, Bond, BondType, Element, Part, Topology};
use nanocad_units::Unit;

use crate::error::PartError;
use crate::generator::PartGenerator;
use crate::geometry::{distance_m, BOND_TOLERANCE_RELATIVE};
use crate::parameter::{ParameterSet, ParameterSpec};

/// The diamond conventional cubic lattice constant at 300 K, in metres.
///
/// Source: CRC Handbook of Chemistry and Physics; N. W. Ashcroft and N. D.
/// Mermin, Solid State Physics (1976), chapter 4. The same value is used by the
/// diamond and spur-gear generators in this crate.
const DIAMOND_LATTICE_CONSTANT_M: f64 = 3.567e-10;

/// The diamond carbon-carbon first-shell bond length `a sqrt(3) / 4`, in
/// metres. Source: the lattice constant above. The same value is used by the
/// spur-gear generator in this crate.
const DIAMONDOID_BOND_M: f64 = 1.544e-10;

/// The eight basis atoms of the diamond cubic cell, in fractional coordinates.
/// Two interpenetrating FCC lattices, offset by (1/4, 1/4, 1/4).
const DIAMOND_BASIS: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0],
    [0.0, 0.5, 0.5],
    [0.5, 0.0, 0.5],
    [0.5, 0.5, 0.0],
    [0.25, 0.25, 0.25],
    [0.25, 0.75, 0.75],
    [0.75, 0.25, 0.75],
    [0.75, 0.75, 0.25],
];

/// The smallest atom count that still describes a body. Below it the cut
/// region missed the lattice and the generator returns an error.
const MIN_BODY_ATOMS: usize = 4;

static ROTOR_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new(
        "rotor_outer_radius_m",
        Some(Unit::Metre),
        3.0e-9,
        0.5e-9,
        50.0e-9,
        false,
        "rotor disk outer radius in metres",
    ),
    ParameterSpec::new(
        "rotor_height_m",
        Some(Unit::Metre),
        1.0e-9,
        0.3e-9,
        20.0e-9,
        false,
        "rotor axial height in metres",
    ),
    ParameterSpec::new(
        "bearing_gap_m",
        Some(Unit::Metre),
        0.8e-9,
        0.1e-9,
        5.0e-9,
        false,
        "radial clearance between rotor and bearing in metres",
    ),
    ParameterSpec::new(
        "bearing_wall_m",
        Some(Unit::Metre),
        1.0e-9,
        0.3e-9,
        10.0e-9,
        false,
        "bearing sleeve wall thickness in metres",
    ),
];

static PUMP_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new(
        "pump_bore_radius_m",
        Some(Unit::Metre),
        2.5e-9,
        0.5e-9,
        50.0e-9,
        false,
        "pump cylinder bore (inner) radius in metres",
    ),
    ParameterSpec::new(
        "pump_wall_m",
        Some(Unit::Metre),
        0.8e-9,
        0.3e-9,
        10.0e-9,
        false,
        "pump cylinder wall thickness in metres",
    ),
    ParameterSpec::new(
        "pump_length_m",
        Some(Unit::Metre),
        2.0e-9,
        0.5e-9,
        30.0e-9,
        false,
        "pump cylinder axial length in metres",
    ),
    ParameterSpec::new(
        "piston_clearance_m",
        Some(Unit::Metre),
        0.3e-9,
        0.05e-9,
        2.0e-9,
        false,
        "radial clearance between piston and bore in metres",
    ),
    ParameterSpec::new(
        "piston_length_m",
        Some(Unit::Metre),
        1.0e-9,
        0.3e-9,
        20.0e-9,
        false,
        "piston axial length in metres",
    ),
    ParameterSpec::new(
        "stroke_m",
        Some(Unit::Metre),
        1.0e-9,
        0.1e-9,
        20.0e-9,
        false,
        "piston stroke length in metres",
    ),
    ParameterSpec::new(
        "seal_clearance_m",
        Some(Unit::Metre),
        0.2e-9,
        0.16e-9,
        1.0e-9,
        false,
        "radial running clearance between the piston seal ring and the bore in metres",
    ),
    ParameterSpec::new(
        "seal_land_m",
        Some(Unit::Metre),
        0.6e-9,
        0.2e-9,
        5.0e-9,
        false,
        "axial length of the piston seal ring land in metres",
    ),
];

static TANK_PARAMETERS: &[ParameterSpec] = &[
    ParameterSpec::new(
        "tank_inner_radius_m",
        Some(Unit::Metre),
        2.0e-9,
        0.5e-9,
        50.0e-9,
        false,
        "tank interior cavity radius in metres",
    ),
    ParameterSpec::new(
        "tank_wall_m",
        Some(Unit::Metre),
        1.0e-9,
        0.3e-9,
        10.0e-9,
        false,
        "tank shell wall thickness in metres",
    ),
];

/// An estimate of flow through a thin annular clearance gap.
///
/// The annulus is unrolled to a plane channel of width `2 * pi * R`. For a
/// Newtonian fluid of dynamic viscosity `mu`, the volumetric flow is the sum
/// of a pressure-driven Poiseuille term and a shear-driven Couette term:
///
/// ```text
/// Q_p = pi * R * h^3 / (6 * mu) * (dp / L)      (Poiseuille)
/// Q_c = pi * R * h * U                          (Couette)
/// ```
///
/// The pressure conductance is `dQ_p / d(dp) = pi * R * h^3 / (6 * mu * L)`.
/// Source: R. B. Bird, W. E. Stewart, E. N. Lightfoot, "Transport Phenomena",
/// 2nd ed. (2002), sections 2.3 and 2.4. This is a thin-gap estimate. It
/// ignores end effects, eccentricity, and surface slip. A positive pressure
/// difference drives flow from the high-pressure end to the low-pressure end.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AnnularGapFlow {
    /// The radial clearance between the two surfaces, in metres.
    pub gap_m: f64,
    /// The mean radius of the annular clearance, in metres.
    pub mean_radius_m: f64,
    /// The axial length of the leak path, in metres.
    pub length_m: f64,
    /// The dynamic viscosity of the fluid, in pascal seconds.
    pub viscosity_pa_s: f64,
    /// The pressure conductance `pi * R * h^3 / (6 * mu * L)`, in
    /// cubic metres per pascal second.
    pub pressure_conductance_m3_per_pa_s: f64,
}

impl AnnularGapFlow {
    /// Builds a thin-gap flow estimate. Every input must be finite and positive.
    pub fn new(
        gap_m: f64,
        mean_radius_m: f64,
        length_m: f64,
        viscosity_pa_s: f64,
    ) -> Result<Self, PartError> {
        for (label, value) in [
            ("gap", gap_m),
            ("mean radius", mean_radius_m),
            ("length", length_m),
            ("viscosity", viscosity_pa_s),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(PartError::InvalidGeometry(format!(
                    "annular gap {label} {value} must be finite and positive"
                )));
            }
        }
        let conductance = std::f64::consts::PI * mean_radius_m * gap_m * gap_m * gap_m
            / (6.0 * viscosity_pa_s * length_m);
        Ok(Self {
            gap_m,
            mean_radius_m,
            length_m,
            viscosity_pa_s,
            pressure_conductance_m3_per_pa_s: conductance,
        })
    }

    /// Returns the pressure conductance in cubic metres per pascal second.
    pub fn pressure_conductance_m3_per_pa_s(&self) -> f64 {
        self.pressure_conductance_m3_per_pa_s
    }

    /// Returns the Poiseuille flow for a pressure difference in pascals.
    ///
    /// The pressure difference must be finite. A negative difference reverses
    /// the flow direction.
    pub fn pressure_driven_flow_m3_per_s(
        &self,
        pressure_difference_pa: f64,
    ) -> Result<f64, PartError> {
        if !pressure_difference_pa.is_finite() {
            return Err(PartError::InvalidGeometry(format!(
                "pressure difference {pressure_difference_pa} Pa must be finite"
            )));
        }
        Ok(self.pressure_conductance_m3_per_pa_s * pressure_difference_pa)
    }

    /// Returns the Couette flow for a sliding speed in metres per second.
    ///
    /// The sliding speed must be finite. It is the relative axial speed of the
    /// two surfaces.
    pub fn couette_flow_m3_per_s(&self, sliding_speed_m_per_s: f64) -> Result<f64, PartError> {
        if !sliding_speed_m_per_s.is_finite() {
            return Err(PartError::InvalidGeometry(format!(
                "sliding speed {sliding_speed_m_per_s} m/s must be finite"
            )));
        }
        Ok(std::f64::consts::PI * self.mean_radius_m * self.gap_m * sliding_speed_m_per_s)
    }

    /// Returns the total flow, the Poiseuille term plus the Couette term.
    pub fn total_flow_m3_per_s(
        &self,
        pressure_difference_pa: f64,
        sliding_speed_m_per_s: f64,
    ) -> Result<f64, PartError> {
        Ok(self.pressure_driven_flow_m3_per_s(pressure_difference_pa)?
            + self.couette_flow_m3_per_s(sliding_speed_m_per_s)?)
    }
}

/// The rotor, bearing, and their design clearances.
///
/// The rotor is a solid diamondoid disk. The bearing is a concentric sleeve.
/// The design gap is `bearing_inner_radius_m - outer_radius_m`. The reported
/// rotational degree of freedom is the free spin of the rotor about the `z`
/// axis inside the sleeve.
#[derive(Clone, Debug, PartialEq)]
pub struct RespirocyteRotor {
    /// The rotor disk.
    pub rotor: Part,
    /// The concentric bearing sleeve.
    pub bearing: Part,
    /// The rotor disk outer radius, in metres.
    pub outer_radius_m: f64,
    /// The rotor disk axial height, in metres.
    pub height_m: f64,
    /// The bearing sleeve inner radius, in metres.
    pub bearing_inner_radius_m: f64,
    /// The bearing sleeve outer radius, in metres.
    pub bearing_outer_radius_m: f64,
    /// The design radial gap, in metres.
    pub gap_m: f64,
    /// The smallest atom-to-atom radial clearance measured from the generated
    /// geometry, in metres. It is never less than `gap_m`.
    pub actual_gap_m: f64,
    /// The number of rotational degrees of freedom. It is one.
    pub rotational_dof: usize,
    /// The rotation axis, a unit vector in the part frame.
    pub rotation_axis: [f64; 3],
}

/// Builds a respirocyte rotor and its concentric bearing sleeve.
#[derive(Clone, Copy, Debug, Default)]
pub struct RespirocyteRotorGenerator;

impl RespirocyteRotorGenerator {
    /// Resolves the inputs and builds the rotor, the bearing, and the metadata.
    pub fn build(&self, parameters: &ParameterSet) -> Result<RespirocyteRotor, PartError> {
        let resolved = self.resolve(parameters)?;
        let rotor_outer_radius_m = resolved.require("rotor_outer_radius_m")?;
        let rotor_height_m = resolved.require("rotor_height_m")?;
        let bearing_gap_m = resolved.require("bearing_gap_m")?;
        let bearing_wall_m = resolved.require("bearing_wall_m")?;

        if rotor_outer_radius_m < DIAMONDOID_BOND_M {
            return Err(PartError::InvalidGeometry(format!(
                "rotor outer radius {rotor_outer_radius_m} m is below one bond length"
            )));
        }
        if rotor_height_m < DIAMONDOID_BOND_M {
            return Err(PartError::InvalidGeometry(format!(
                "rotor height {rotor_height_m} m is below one bond length"
            )));
        }
        if bearing_wall_m < DIAMONDOID_BOND_M {
            return Err(PartError::InvalidGeometry(format!(
                "bearing wall {bearing_wall_m} m is below one bond length"
            )));
        }

        let bearing_inner_radius_m = rotor_outer_radius_m + bearing_gap_m;
        let bearing_outer_radius_m = bearing_inner_radius_m + bearing_wall_m;
        let half_height_m = rotor_height_m / 2.0;

        let mut rotor_topology = Topology::new();
        add_lattice_atoms(&mut rotor_topology, rotor_outer_radius_m, |position_m| {
            radial_xy_m(position_m) <= rotor_outer_radius_m && position_m[2].abs() <= half_height_m
        });
        add_first_shell_bonds(&mut rotor_topology, DIAMONDOID_BOND_M)?;
        if rotor_topology.atom_count() < MIN_BODY_ATOMS {
            return Err(PartError::InvalidGeometry(
                "rotor region cut no atoms from the diamond lattice".to_owned(),
            ));
        }

        let mut bearing_topology = Topology::new();
        add_lattice_atoms(
            &mut bearing_topology,
            bearing_outer_radius_m,
            |position_m| {
                let radial_m = radial_xy_m(position_m);
                radial_m >= bearing_inner_radius_m
                    && radial_m <= bearing_outer_radius_m
                    && position_m[2].abs() <= half_height_m
            },
        );
        add_first_shell_bonds(&mut bearing_topology, DIAMONDOID_BOND_M)?;
        if bearing_topology.atom_count() < MIN_BODY_ATOMS {
            return Err(PartError::InvalidGeometry(
                "bearing region cut no atoms from the diamond lattice".to_owned(),
            ));
        }

        let actual_gap_m = minimum_atom_gap_m(&rotor_topology, &bearing_topology);

        let mut rotor = Part::new("respirocyte-rotor", rotor_topology).with_material("diamondoid");
        rotor
            .metadata
            .insert("generator".to_owned(), self.id().to_owned());
        rotor.metadata.insert(
            "outer_radius_m".to_owned(),
            format!("{rotor_outer_radius_m:e}"),
        );
        rotor
            .metadata
            .insert("rotational_dof".to_owned(), "1".to_owned());

        let mut bearing =
            Part::new("respirocyte-bearing", bearing_topology).with_material("diamondoid");
        bearing
            .metadata
            .insert("generator".to_owned(), self.id().to_owned());
        bearing.metadata.insert(
            "inner_radius_m".to_owned(),
            format!("{bearing_inner_radius_m:e}"),
        );

        Ok(RespirocyteRotor {
            rotor,
            bearing,
            outer_radius_m: rotor_outer_radius_m,
            height_m: rotor_height_m,
            bearing_inner_radius_m,
            bearing_outer_radius_m,
            gap_m: bearing_gap_m,
            actual_gap_m,
            rotational_dof: 1,
            rotation_axis: [0.0, 0.0, 1.0],
        })
    }
}

impl PartGenerator for RespirocyteRotorGenerator {
    fn id(&self) -> &'static str {
        "respirocyte_rotor"
    }

    fn name(&self) -> &'static str {
        "Respirocyte rotor and bearing"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        ROTOR_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        Ok(self.build(parameters)?.rotor)
    }
}

impl RespirocyteRotor {
    /// Estimates the annular leak through the design bearing gap.
    ///
    /// The leak path is the design radial gap `gap_m` along the rotor height
    /// `height_m` at the mean radius. Pass the fluid dynamic viscosity in pascal
    /// seconds. The returned [`AnnularGapFlow`] gives the Poiseuille flow for a
    /// pressure difference and the Couette flow for a sliding speed.
    pub fn bearing_leakage_model(&self, viscosity_pa_s: f64) -> Result<AnnularGapFlow, PartError> {
        let mean_radius_m = self.outer_radius_m + 0.5 * self.gap_m;
        AnnularGapFlow::new(self.gap_m, mean_radius_m, self.height_m, viscosity_pa_s)
    }
}

/// The pump cylinder, piston, and the stroke.
#[derive(Clone, Debug, PartialEq)]
pub struct RespirocytePump {
    /// The hollow pump cylinder.
    pub cylinder: Part,
    /// The solid piston, parked at the bottom of its stroke.
    pub piston: Part,
    /// The seal ring land carried on the piston top.
    pub seal: Part,
    /// The cylinder bore (inner) radius, in metres.
    pub bore_radius_m: f64,
    /// The cylinder outer radius, in metres.
    pub outer_radius_m: f64,
    /// The cylinder wall thickness, in metres.
    pub wall_m: f64,
    /// The piston radius, in metres.
    pub piston_radius_m: f64,
    /// The cylinder axial length, in metres.
    pub cylinder_length_m: f64,
    /// The piston axial length, in metres.
    pub piston_length_m: f64,
    /// The piston stroke length, in metres.
    pub stroke_m: f64,
    /// The seal ring inner radius, in metres. It equals the piston radius.
    pub seal_inner_radius_m: f64,
    /// The seal ring outer radius, in metres.
    pub seal_outer_radius_m: f64,
    /// The radial running clearance between the seal ring and the bore, in
    /// metres.
    pub seal_clearance_m: f64,
    /// The axial length of the seal ring land, in metres.
    pub seal_land_m: f64,
}

/// Builds a respirocyte cylinder and piston pump.
#[derive(Clone, Copy, Debug, Default)]
pub struct RespirocytePumpGenerator;

impl RespirocytePumpGenerator {
    /// Resolves the inputs and builds the cylinder, the piston, and the metadata.
    pub fn build(&self, parameters: &ParameterSet) -> Result<RespirocytePump, PartError> {
        let resolved = self.resolve(parameters)?;
        let bore_radius_m = resolved.require("pump_bore_radius_m")?;
        let wall_m = resolved.require("pump_wall_m")?;
        let cylinder_length_m = resolved.require("pump_length_m")?;
        let piston_clearance_m = resolved.require("piston_clearance_m")?;
        let piston_length_m = resolved.require("piston_length_m")?;
        let stroke_m = resolved.require("stroke_m")?;
        let seal_clearance_m = resolved.require("seal_clearance_m")?;
        let seal_land_m = resolved.require("seal_land_m")?;

        if wall_m < DIAMONDOID_BOND_M {
            return Err(PartError::InvalidGeometry(format!(
                "pump wall {wall_m} m is below one bond length"
            )));
        }
        if piston_length_m < DIAMONDOID_BOND_M {
            return Err(PartError::InvalidGeometry(format!(
                "piston length {piston_length_m} m is below one bond length"
            )));
        }
        if piston_clearance_m >= bore_radius_m {
            return Err(PartError::InvalidGeometry(format!(
                "piston clearance {piston_clearance_m} m is not smaller than the bore \
                 radius {bore_radius_m} m"
            )));
        }
        if piston_length_m + stroke_m > cylinder_length_m {
            return Err(PartError::InvalidGeometry(format!(
                "piston length {piston_length_m} m plus stroke {stroke_m} m exceeds the \
                 cylinder length {cylinder_length_m} m"
            )));
        }
        if seal_clearance_m <= DIAMONDOID_BOND_M {
            return Err(PartError::InvalidGeometry(format!(
                "seal clearance {seal_clearance_m} m is not above one bond length, so the \
                 ring and the bore would bond across the running gap"
            )));
        }
        if seal_clearance_m >= piston_clearance_m {
            return Err(PartError::InvalidGeometry(format!(
                "seal clearance {seal_clearance_m} m is not smaller than the piston \
                 clearance {piston_clearance_m} m, so the seal does not tighten the gap"
            )));
        }
        if seal_land_m > piston_length_m {
            return Err(PartError::InvalidGeometry(format!(
                "seal land {seal_land_m} m is longer than the piston {piston_length_m} m"
            )));
        }

        let outer_radius_m = bore_radius_m + wall_m;
        let piston_radius_m = bore_radius_m - piston_clearance_m;
        let half_length_m = cylinder_length_m / 2.0;

        let mut cylinder_topology = Topology::new();
        add_lattice_atoms(&mut cylinder_topology, outer_radius_m, |position_m| {
            let radial_m = radial_xy_m(position_m);
            radial_m >= bore_radius_m
                && radial_m <= outer_radius_m
                && position_m[2].abs() <= half_length_m
        });
        add_first_shell_bonds(&mut cylinder_topology, DIAMONDOID_BOND_M)?;
        if cylinder_topology.atom_count() < MIN_BODY_ATOMS {
            return Err(PartError::InvalidGeometry(
                "pump cylinder region cut no atoms from the diamond lattice".to_owned(),
            ));
        }

        let piston_z_min_m = -half_length_m;
        let piston_z_max_m = piston_z_min_m + piston_length_m;
        let mut piston_topology = Topology::new();
        add_lattice_atoms(&mut piston_topology, piston_radius_m, |position_m| {
            radial_xy_m(position_m) <= piston_radius_m
                && position_m[2] >= piston_z_min_m
                && position_m[2] <= piston_z_max_m
        });
        add_first_shell_bonds(&mut piston_topology, DIAMONDOID_BOND_M)?;
        if piston_topology.atom_count() < MIN_BODY_ATOMS {
            return Err(PartError::InvalidGeometry(
                "piston region cut no atoms from the diamond lattice".to_owned(),
            ));
        }

        let seal_inner_radius_m = piston_radius_m;
        let seal_outer_radius_m = bore_radius_m - seal_clearance_m;
        let seal_z_max_m = piston_z_max_m;
        let seal_z_min_m = seal_z_max_m - seal_land_m;
        let mut seal_topology = Topology::new();
        add_lattice_atoms(&mut seal_topology, seal_outer_radius_m, |position_m| {
            let radial_m = radial_xy_m(position_m);
            radial_m >= seal_inner_radius_m
                && radial_m <= seal_outer_radius_m
                && position_m[2] >= seal_z_min_m
                && position_m[2] <= seal_z_max_m
        });
        add_first_shell_bonds(&mut seal_topology, DIAMONDOID_BOND_M)?;
        if seal_topology.atom_count() < MIN_BODY_ATOMS {
            return Err(PartError::InvalidGeometry(
                "seal ring region cut no atoms from the diamond lattice".to_owned(),
            ));
        }

        let mut cylinder =
            Part::new("respirocyte-pump-cylinder", cylinder_topology).with_material("diamondoid");
        cylinder
            .metadata
            .insert("generator".to_owned(), self.id().to_owned());
        cylinder
            .metadata
            .insert("stroke_m".to_owned(), format!("{stroke_m:e}"));

        let mut piston =
            Part::new("respirocyte-pump-piston", piston_topology).with_material("diamondoid");
        piston
            .metadata
            .insert("generator".to_owned(), self.id().to_owned());

        let mut seal =
            Part::new("respirocyte-pump-seal", seal_topology).with_material("diamondoid");
        seal.metadata
            .insert("generator".to_owned(), self.id().to_owned());
        seal.metadata.insert(
            "seal_clearance_m".to_owned(),
            format!("{seal_clearance_m:e}"),
        );

        Ok(RespirocytePump {
            cylinder,
            piston,
            seal,
            bore_radius_m,
            outer_radius_m,
            wall_m,
            piston_radius_m,
            cylinder_length_m,
            piston_length_m,
            stroke_m,
            seal_inner_radius_m,
            seal_outer_radius_m,
            seal_clearance_m,
            seal_land_m,
        })
    }
}

impl PartGenerator for RespirocytePumpGenerator {
    fn id(&self) -> &'static str {
        "respirocyte_pump"
    }

    fn name(&self) -> &'static str {
        "Respirocyte pump cylinder and piston"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        PUMP_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        Ok(self.build(parameters)?.cylinder)
    }
}

impl RespirocytePump {
    /// Estimates the annular leak past the piston seal ring.
    ///
    /// The leak path is the seal running clearance `seal_clearance_m` along the
    /// land `seal_land_m` at the seal mean radius. Pass the fluid dynamic
    /// viscosity in pascal seconds.
    pub fn seal_leakage_model(&self, viscosity_pa_s: f64) -> Result<AnnularGapFlow, PartError> {
        let mean_radius_m = self.seal_outer_radius_m + 0.5 * self.seal_clearance_m;
        AnnularGapFlow::new(
            self.seal_clearance_m,
            mean_radius_m,
            self.seal_land_m,
            viscosity_pa_s,
        )
    }
}

/// The gas tank shell and its interior cavity.
#[derive(Clone, Debug, PartialEq)]
pub struct RespirocyteTank {
    /// The hollow diamondoid shell.
    pub shell: Part,
    /// The interior cavity radius, in metres.
    pub inner_radius_m: f64,
    /// The shell outer radius, in metres.
    pub outer_radius_m: f64,
    /// The shell wall thickness, in metres.
    pub wall_m: f64,
}

/// Builds a hollow respirocyte gas tank.
#[derive(Clone, Copy, Debug, Default)]
pub struct RespirocyteTankGenerator;

impl RespirocyteTankGenerator {
    /// Resolves the inputs and builds the shell and its metadata.
    pub fn build(&self, parameters: &ParameterSet) -> Result<RespirocyteTank, PartError> {
        let resolved = self.resolve(parameters)?;
        let inner_radius_m = resolved.require("tank_inner_radius_m")?;
        let wall_m = resolved.require("tank_wall_m")?;

        if wall_m < DIAMONDOID_BOND_M {
            return Err(PartError::InvalidGeometry(format!(
                "tank wall {wall_m} m is below one bond length"
            )));
        }
        if inner_radius_m < DIAMONDOID_BOND_M {
            return Err(PartError::InvalidGeometry(format!(
                "tank inner radius {inner_radius_m} m is below one bond length"
            )));
        }

        let outer_radius_m = inner_radius_m + wall_m;
        let mut shell_topology = Topology::new();
        add_lattice_atoms(&mut shell_topology, outer_radius_m, |position_m| {
            let radial_m = radial_m(position_m);
            radial_m >= inner_radius_m && radial_m <= outer_radius_m
        });
        add_first_shell_bonds(&mut shell_topology, DIAMONDOID_BOND_M)?;
        if shell_topology.atom_count() < MIN_BODY_ATOMS {
            return Err(PartError::InvalidGeometry(
                "tank shell region cut no atoms from the diamond lattice".to_owned(),
            ));
        }

        let mut shell = Part::new("respirocyte-tank", shell_topology).with_material("diamondoid");
        shell
            .metadata
            .insert("generator".to_owned(), self.id().to_owned());
        shell
            .metadata
            .insert("inner_radius_m".to_owned(), format!("{inner_radius_m:e}"));
        shell
            .metadata
            .insert("outer_radius_m".to_owned(), format!("{outer_radius_m:e}"));

        Ok(RespirocyteTank {
            shell,
            inner_radius_m,
            outer_radius_m,
            wall_m,
        })
    }
}

impl PartGenerator for RespirocyteTankGenerator {
    fn id(&self) -> &'static str {
        "respirocyte_tank"
    }

    fn name(&self) -> &'static str {
        "Respirocyte gas tank shell"
    }

    fn parameters(&self) -> &'static [ParameterSpec] {
        TANK_PARAMETERS
    }

    fn generate(&self, parameters: &ParameterSet) -> Result<Part, PartError> {
        Ok(self.build(parameters)?.shell)
    }
}

/// Builds a rotor and bearing pair from a parameter set.
pub fn respirocyte_rotor(parameters: &ParameterSet) -> Result<RespirocyteRotor, PartError> {
    RespirocyteRotorGenerator.build(parameters)
}

/// Builds a pump cylinder and piston pair from a parameter set.
pub fn respirocyte_pump(parameters: &ParameterSet) -> Result<RespirocytePump, PartError> {
    RespirocytePumpGenerator.build(parameters)
}

/// Builds a gas tank shell from a parameter set.
pub fn respirocyte_tank(parameters: &ParameterSet) -> Result<RespirocyteTank, PartError> {
    RespirocyteTankGenerator.build(parameters)
}

/// Adds every diamond lattice atom inside an extent box that `keep` accepts.
///
/// The cut region is stated by `keep`. The function visits the conventional
/// cells that cover the extent box, so no lattice site inside the box is missed.
fn add_lattice_atoms(topology: &mut Topology, extent_m: f64, keep: impl Fn([f64; 3]) -> bool) {
    let a_m = DIAMOND_LATTICE_CONSTANT_M;
    let cells = (extent_m / a_m).ceil() as i64 + 1;
    for ix in -cells..=cells {
        for iy in -cells..=cells {
            for iz in -cells..=cells {
                for basis in DIAMOND_BASIS {
                    let position_m = [
                        (ix as f64 + basis[0]) * a_m,
                        (iy as f64 + basis[1]) * a_m,
                        (iz as f64 + basis[2]) * a_m,
                    ];
                    if keep(position_m) {
                        topology.add_atom(Atom::new(Element::CARBON, position_m, 0.0, "C"));
                    }
                }
            }
        }
    }
}

/// Adds one single bond for every atom pair at the diamond first-shell distance.
///
/// A grid of cells the width of the tolerance upper bound turns the search into
/// a scan of the 27 neighbouring cells. The search visits each unordered pair
/// once, so a bond cannot duplicate.
fn add_first_shell_bonds(topology: &mut Topology, expected_m: f64) -> Result<(), PartError> {
    let low_m = expected_m * (1.0 - BOND_TOLERANCE_RELATIVE);
    let high_m = expected_m * (1.0 + BOND_TOLERANCE_RELATIVE);
    let cell_m = high_m;

    let mut buckets: HashMap<(i64, i64, i64), Vec<u32>> = HashMap::new();
    for index in 0..topology.atom_count() {
        let Some(position_m) = topology.position_m(index) else {
            continue;
        };
        buckets
            .entry(cell_key(position_m, cell_m))
            .or_default()
            .push(index as u32);
    }

    for u in 0..topology.atom_count() {
        let Some(u_m) = topology.position_m(u) else {
            continue;
        };
        let base = cell_key(u_m, cell_m);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let Some(neighbors) = buckets.get(&(base.0 + dx, base.1 + dy, base.2 + dz))
                    else {
                        continue;
                    };
                    for &v in neighbors {
                        if v as usize <= u {
                            continue;
                        }
                        let Some(v_m) = topology.position_m(v as usize) else {
                            continue;
                        };
                        let length_m = distance_m(u_m, v_m);
                        if length_m >= low_m && length_m <= high_m {
                            topology.add_bond(Bond::new(u as u32, v, 1, BondType::Single))?;
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

fn cell_key(position_m: [f64; 3], cell_m: f64) -> (i64, i64, i64) {
    (
        (position_m[0] / cell_m).floor() as i64,
        (position_m[1] / cell_m).floor() as i64,
        (position_m[2] / cell_m).floor() as i64,
    )
}

/// Returns the smallest radial clearance between two atom sets, in metres.
///
/// The signed radial difference `bearing_radius - rotor_radius` separates, so
/// the minimum over all pairs is the smallest bearing radius minus the largest
/// rotor radius. A positive result means the bodies do not overlap radially.
fn minimum_atom_gap_m(rotor: &Topology, bearing: &Topology) -> f64 {
    let max_rotor_radius_m = rotor
        .atoms()
        .map(|atom| radial_xy_m(atom.position_m))
        .fold(0.0_f64, f64::max);
    let min_bearing_radius_m = bearing
        .atoms()
        .map(|atom| radial_xy_m(atom.position_m))
        .fold(f64::INFINITY, f64::min);
    min_bearing_radius_m - max_rotor_radius_m
}

/// Returns the smallest atom-to-atom distance between two atom sets, in metres.
///
/// The two sets are separate bodies, so a distance inside the bonding band
/// means a bond would form if the sets shared a topology.
#[cfg(test)]
fn minimum_pair_distance_m(a: &Topology, b: &Topology) -> f64 {
    let mut minimum_m = f64::INFINITY;
    for atom_a in a.atoms() {
        for atom_b in b.atoms() {
            minimum_m = minimum_m.min(distance_m(atom_a.position_m, atom_b.position_m));
        }
    }
    minimum_m
}

fn radial_xy_m(position_m: [f64; 3]) -> f64 {
    (position_m[0] * position_m[0] + position_m[1] * position_m[1]).sqrt()
}

fn radial_m(position_m: [f64; 3]) -> f64 {
    (position_m[0] * position_m[0] + position_m[1] * position_m[1] + position_m[2] * position_m[2])
        .sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_bonds_match_the_diamond_bond(part: &Part) {
        assert!(part.bond_count() > 0, "expected bonds");
        for bond in part.topology.bonds() {
            let u_m = part.topology.position_m(bond.u as usize).expect("atom u");
            let v_m = part.topology.position_m(bond.v as usize).expect("atom v");
            let length_m = distance_m(u_m, v_m);
            let error_m = (length_m - DIAMONDOID_BOND_M).abs();
            assert!(
                error_m <= DIAMONDOID_BOND_M * BOND_TOLERANCE_RELATIVE,
                "bond length {length_m:e} m differs from {DIAMONDOID_BOND_M:e} m"
            );
        }
    }

    #[test]
    fn atom_counts_are_nonzero_and_stable() {
        let first = RespirocyteRotorGenerator
            .build(&ParameterSet::new())
            .expect("rotor");
        let second = RespirocyteRotorGenerator
            .build(&ParameterSet::new())
            .expect("rotor");
        assert!(first.rotor.atom_count() >= MIN_BODY_ATOMS);
        assert!(first.bearing.atom_count() >= MIN_BODY_ATOMS);
        assert_eq!(first.rotor.atom_count(), second.rotor.atom_count());
        assert_eq!(first.bearing.atom_count(), second.bearing.atom_count());
        assert_eq!(first.rotor.material, "diamondoid");
        assert_eq!(first.bearing.name, "respirocyte-bearing");
        assert_eq!(first.rotational_dof, 1);
        assert_eq!(first.rotation_axis, [0.0, 0.0, 1.0]);
    }

    #[test]
    fn the_bearing_gap_is_positive_and_within_design_tolerance() {
        let build = RespirocyteRotorGenerator
            .build(&ParameterSet::new())
            .expect("rotor");
        assert!(build.gap_m > 0.0);
        assert!(build.actual_gap_m > 0.0, "gap {}", build.actual_gap_m);
        assert!(
            build.actual_gap_m >= build.gap_m - f64::EPSILON,
            "actual gap {:e} m is below the design gap {:e} m",
            build.actual_gap_m,
            build.gap_m
        );
        assert!(
            build.actual_gap_m <= build.gap_m + 2.0 * DIAMOND_LATTICE_CONSTANT_M,
            "actual gap {:e} m exceeds the design gap {:e} m by too much",
            build.actual_gap_m,
            build.gap_m
        );
        assert!(build.bearing_inner_radius_m > build.outer_radius_m);
    }

    #[test]
    fn rotor_and_bearing_bonds_match_the_diamond_bond() {
        let build = RespirocyteRotorGenerator
            .build(&ParameterSet::new())
            .expect("rotor");
        assert_bonds_match_the_diamond_bond(&build.rotor);
        assert_bonds_match_the_diamond_bond(&build.bearing);
    }

    #[test]
    fn the_tank_interior_is_empty_of_atoms() {
        let build = RespirocyteTankGenerator
            .build(&ParameterSet::new())
            .expect("tank");
        assert!(build.shell.atom_count() >= MIN_BODY_ATOMS);
        assert!(build.inner_radius_m > 0.0);
        assert!(build.outer_radius_m > build.inner_radius_m);
        for atom in build.shell.topology.atoms() {
            let radius_m = radial_m(atom.position_m);
            assert!(
                radius_m >= build.inner_radius_m,
                "atom at radius {radius_m:e} m is inside the tank cavity {:e} m",
                build.inner_radius_m
            );
            assert!(radius_m <= build.outer_radius_m + f64::EPSILON);
        }
        assert_bonds_match_the_diamond_bond(&build.shell);
    }

    #[test]
    fn the_pump_stroke_is_positive_and_the_piston_fits() {
        let build = RespirocytePumpGenerator
            .build(&ParameterSet::new())
            .expect("pump");
        assert!(build.stroke_m > 0.0);
        assert!(build.piston_radius_m > 0.0);
        assert!(build.piston_radius_m < build.bore_radius_m);
        assert!(build.outer_radius_m > build.bore_radius_m);
        assert!((build.outer_radius_m - build.bore_radius_m - build.wall_m).abs() <= 1.0e-24);
        assert!(build.cylinder.atom_count() >= MIN_BODY_ATOMS);
        assert!(build.piston.atom_count() >= MIN_BODY_ATOMS);
        assert_bonds_match_the_diamond_bond(&build.cylinder);
        assert_bonds_match_the_diamond_bond(&build.piston);
    }

    #[test]
    fn the_piston_sits_inside_the_bore() {
        let build = RespirocytePumpGenerator
            .build(&ParameterSet::new())
            .expect("pump");
        for atom in build.piston.topology.atoms() {
            assert!(radial_xy_m(atom.position_m) <= build.bore_radius_m);
        }
    }

    /// Hand calculation of the thin-gap annular Poiseuille flow.
    fn annulus_poiseuille_flow_m3_per_s(
        gap_m: f64,
        mean_radius_m: f64,
        length_m: f64,
        viscosity_pa_s: f64,
        pressure_difference_pa: f64,
    ) -> f64 {
        std::f64::consts::PI * mean_radius_m * gap_m.powi(3) * pressure_difference_pa
            / (6.0 * viscosity_pa_s * length_m)
    }

    #[test]
    fn the_pump_seal_tightens_the_running_gap() {
        let build = RespirocytePumpGenerator
            .build(&ParameterSet::new())
            .expect("pump");
        assert!(build.seal.atom_count() >= MIN_BODY_ATOMS);
        assert!(build.seal_clearance_m > 0.0);
        assert!(build.seal_clearance_m < build.bore_radius_m - build.piston_radius_m);
        assert!((build.seal_inner_radius_m - build.piston_radius_m).abs() <= f64::EPSILON);
        let expected_outer_m = build.bore_radius_m - build.seal_clearance_m;
        assert!((build.seal_outer_radius_m - expected_outer_m).abs() <= f64::EPSILON);
        assert!(build.seal_inner_radius_m < build.seal_outer_radius_m);
        assert_eq!(build.seal.name, "respirocyte-pump-seal");
        assert_bonds_match_the_diamond_bond(&build.seal);
    }

    #[test]
    fn no_bond_crosses_the_seal_running_gap() {
        let build = RespirocytePumpGenerator
            .build(&ParameterSet::new())
            .expect("pump");
        let bond_upper_m = DIAMONDOID_BOND_M * (1.0 + BOND_TOLERANCE_RELATIVE);
        let seal_to_bore_m =
            minimum_pair_distance_m(&build.seal.topology, &build.cylinder.topology);
        assert!(
            seal_to_bore_m > bond_upper_m,
            "seal-to-bore minimum distance {seal_to_bore_m:e} m is inside the bonding band"
        );
        assert!(
            build.seal_clearance_m <= seal_to_bore_m + f64::EPSILON,
            "measured gap {seal_to_bore_m:e} m is tighter than the design clearance {:e} m",
            build.seal_clearance_m
        );
        let seal_to_piston_m =
            minimum_pair_distance_m(&build.seal.topology, &build.piston.topology);
        assert!(
            seal_to_piston_m > 0.0,
            "the seal overlaps the piston by {seal_to_piston_m:e} m"
        );
    }

    #[test]
    fn the_rotor_bearing_leakage_matches_the_hand_calculation() {
        let build = RespirocyteRotorGenerator
            .build(&ParameterSet::new())
            .expect("rotor");
        let viscosity_pa_s = 1.0e-3;
        let pressure_difference_pa = 1.0e5;
        let model = build
            .bearing_leakage_model(viscosity_pa_s)
            .expect("valid model");
        let mean_radius_m = build.outer_radius_m + 0.5 * build.gap_m;
        let expected_m3_per_s = annulus_poiseuille_flow_m3_per_s(
            build.gap_m,
            mean_radius_m,
            build.height_m,
            viscosity_pa_s,
            pressure_difference_pa,
        );
        let actual_m3_per_s = model
            .pressure_driven_flow_m3_per_s(pressure_difference_pa)
            .expect("finite pressure");
        let relative = (actual_m3_per_s - expected_m3_per_s).abs() / expected_m3_per_s;
        assert!(relative < 1.0e-12, "leakage relative error {relative}");
        let rounded_m3_per_s = 9.1148e-20;
        assert!(
            (actual_m3_per_s - rounded_m3_per_s).abs() / rounded_m3_per_s < 1.0e-3,
            "leakage {actual_m3_per_s:e} m^3/s differs from the hand value {rounded_m3_per_s:e}"
        );
        assert_close_rel(
            model.pressure_conductance_m3_per_pa_s(),
            expected_m3_per_s / pressure_difference_pa,
            1.0e-12,
            "conductance",
        );
        println!(
            "bearing leakage: {actual_m3_per_s:.6e} m^3/s at {pressure_difference_pa:.1e} Pa \
             (hand {expected_m3_per_s:.6e}) conductance {:.6e} m^3/(Pa*s)",
            model.pressure_conductance_m3_per_pa_s()
        );
    }

    #[test]
    fn the_seal_leakage_matches_the_hand_calculation() {
        let build = RespirocytePumpGenerator
            .build(&ParameterSet::new())
            .expect("pump");
        let viscosity_pa_s = 1.0e-3;
        let pressure_difference_pa = 1.0e5;
        let model = build
            .seal_leakage_model(viscosity_pa_s)
            .expect("valid model");
        let mean_radius_m = build.seal_outer_radius_m + 0.5 * build.seal_clearance_m;
        let expected_m3_per_s = annulus_poiseuille_flow_m3_per_s(
            build.seal_clearance_m,
            mean_radius_m,
            build.seal_land_m,
            viscosity_pa_s,
            pressure_difference_pa,
        );
        let actual_m3_per_s = model
            .pressure_driven_flow_m3_per_s(pressure_difference_pa)
            .expect("finite pressure");
        let relative = (actual_m3_per_s - expected_m3_per_s).abs() / expected_m3_per_s;
        assert!(relative < 1.0e-12, "leakage relative error {relative}");
        println!(
            "seal leakage: {actual_m3_per_s:.6e} m^3/s at {pressure_difference_pa:.1e} Pa \
             (hand {expected_m3_per_s:.6e})"
        );
    }

    #[test]
    fn the_couette_term_matches_the_plane_channel_estimate() {
        let model = AnnularGapFlow::new(0.5e-9, 2.0e-9, 1.0e-9, 1.0e-3).expect("valid model");
        let speed_m_per_s = 1.0e-3;
        let expected_m3_per_s =
            std::f64::consts::PI * model.mean_radius_m * model.gap_m * speed_m_per_s;
        let actual_m3_per_s = model
            .couette_flow_m3_per_s(speed_m_per_s)
            .expect("finite speed");
        assert!((actual_m3_per_s - expected_m3_per_s).abs() <= 1.0e-30);
        let total_m3_per_s = model
            .total_flow_m3_per_s(1.0e5, speed_m_per_s)
            .expect("finite inputs");
        let pressure_only_m3_per_s = model
            .pressure_driven_flow_m3_per_s(1.0e5)
            .expect("finite pressure");
        assert!((total_m3_per_s - pressure_only_m3_per_s - actual_m3_per_s).abs() <= 1.0e-30);
        assert!(AnnularGapFlow::new(0.0, 2.0e-9, 1.0e-9, 1.0e-3).is_err());
        assert!(AnnularGapFlow::new(0.5e-9, 2.0e-9, 1.0e-9, -1.0).is_err());
    }

    fn assert_close_rel(actual: f64, expected: f64, tolerance: f64, label: &str) {
        let relative = (actual - expected).abs() / expected.abs();
        assert!(
            relative <= tolerance,
            "{label}: actual {actual} expected {expected} relative {relative}"
        );
    }

    #[test]
    fn parameters_out_of_range_are_errors_not_panics() {
        assert!(RespirocyteRotorGenerator
            .build(&ParameterSet::new().with("bearing_gap_m", -1.0))
            .is_err());
        assert!(RespirocyteRotorGenerator
            .build(&ParameterSet::new().with("rotor_outer_radius_m", 0.0))
            .is_err());
        assert!(RespirocyteRotorGenerator
            .build(&ParameterSet::new().with("nope", 1.0))
            .is_err());
        assert!(RespirocytePumpGenerator
            .build(&ParameterSet::new().with("stroke_m", 0.0))
            .is_err());
        assert!(RespirocytePumpGenerator
            .build(&ParameterSet::new().with("piston_clearance_m", 4.0e-9))
            .is_err());
        assert!(RespirocytePumpGenerator
            .build(&ParameterSet::new().with("seal_clearance_m", 1.0e-9))
            .is_err());
        assert!(RespirocyteTankGenerator
            .build(&ParameterSet::new().with("tank_wall_m", 0.0))
            .is_err());
        assert!(RespirocyteTankGenerator
            .build(&ParameterSet::new().with("tank_inner_radius_m", -1.0))
            .is_err());
    }

    #[test]
    fn every_generator_builds_with_defaults_through_the_trait() {
        assert!(RespirocyteRotorGenerator.generate_with_defaults().is_ok());
        assert!(RespirocytePumpGenerator.generate_with_defaults().is_ok());
        assert!(RespirocyteTankGenerator.generate_with_defaults().is_ok());
        assert!(respirocyte_rotor(&ParameterSet::new()).is_ok());
        assert!(respirocyte_pump(&ParameterSet::new()).is_ok());
        assert!(respirocyte_tank(&ParameterSet::new()).is_ok());
    }

    #[test]
    fn a_generator_reports_its_identity() {
        assert_eq!(RespirocyteRotorGenerator.id(), "respirocyte_rotor");
        assert_eq!(
            RespirocyteRotorGenerator.name(),
            "Respirocyte rotor and bearing"
        );
        assert_eq!(RespirocyteRotorGenerator.parameters().len(), 4);
        assert_eq!(RespirocytePumpGenerator.id(), "respirocyte_pump");
        assert_eq!(RespirocytePumpGenerator.parameters().len(), 8);
        assert_eq!(RespirocyteTankGenerator.id(), "respirocyte_tank");
        assert_eq!(RespirocyteTankGenerator.parameters().len(), 2);
    }
}
