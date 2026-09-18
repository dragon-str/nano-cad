//! Frozen model charges for the functional groups.
//!
//! Each group is modelled as that group on a methyl carbon. The geometry comes
//! from RDKit 2026.03.6 with ETKDGv3, the seed 20260917 and MMFF94, and each
//! charge is the Gasteiger-Marsili PEOE value with 12 iterations, multiplied
//! by the elementary charge to give coulombs. The method
//! matches the guest molecules, so the wall and the guest share one rule.
//!
//! The host charge is the charge of the methyl carbon that carries the group.
//! The lattice carbon that carries a group in a pocket takes that charge. The
//! remaining lattice carbons keep zero charge.
//!
//! The charges are a model. They are not computed for the pocket geometry and
//! they are not validated against experiment. See ADR-0063.

/// The stated model charges of one functional group.
pub(crate) struct GroupCharges {
    /// The charge of the lattice carbon that carries the group.
    pub(crate) host_charge_c: f64,
    /// The charge of the group heavy atom.
    pub(crate) heavy_charge_c: f64,
    /// The charge of each hydrogen on the group heavy atom.
    pub(crate) hydrogen_charge_c: f64,
}

/// The group charges, in the order of `FUNCTIONAL_GROUPS`.
pub(crate) const GROUP_CHARGES: [GroupCharges; 6] = [
    GroupCharges {
        host_charge_c: 5.117461713e-21,
        heavy_charge_c: -6.402782385e-20,
        hydrogen_charge_c: 3.358637501e-20,
    },
    GroupCharges {
        host_charge_c: -3.124625755e-21,
        heavy_charge_c: -5.339558924e-20,
        hydrogen_charge_c: 1.892131821e-20,
    },
    GroupCharges {
        host_charge_c: -1.093683957e-20,
        heavy_charge_c: -1.093683957e-20,
        hydrogen_charge_c: 3.645613190e-21,
    },
    GroupCharges {
        host_charge_c: 1.257836145e-20,
        heavy_charge_c: -4.081220350e-20,
        hydrogen_charge_c: 0.000000000e+00,
    },
    GroupCharges {
        host_charge_c: 1.736628986e-21,
        heavy_charge_c: -2.089729145e-20,
        hydrogen_charge_c: 0.000000000e+00,
    },
    GroupCharges {
        host_charge_c: -3.436801071e-21,
        heavy_charge_c: -2.928839119e-20,
        hydrogen_charge_c: 1.632850683e-20,
    },
];
