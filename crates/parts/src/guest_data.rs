//! Frozen guest molecule data.
//!
//! Each geometry comes from RDKit 2026.03.6: ETKDGv3 embedding with the
//! fixed seed 20260917, then MMFF94 optimization. Each partial charge is the
//! Gasteiger-Marsili PEOE value that RDKit computes for that geometry with 12
//! iterations. The charge of a molecule sums to zero.
//!
//! The charges are a model. They are not computed for the pocket, and they are
//! not validated against experiment. Treat every result that uses them as a
//! model result. See ADR-0063.

use nanocad_model::Element;

use crate::guest::{Guest, GuestAtom, GuestBond};

/// Every guest molecule in the table.
pub(crate) const GUEST_DATA: &[Guest] = &[
    Guest {
        id: "methanol",
        label: "Methanol",
        formula: "CH4O",
        molar_mass_kg_per_mol: 3.204200000e-02,
        atoms: &[
            GuestAtom {
                element: Element::CARBON,
                position_m: [-3.700288e-11, -2.627629e-12, -1.387555e-12],
                charge_c: 3.194068e-02,
                atom_type: "C_sp3",
            },
            GuestAtom {
                element: Element::OXYGEN,
                position_m: [9.634647e-11, -4.511955e-11, -2.272788e-11],
                charge_c: -3.996302e-01,
                atom_type: "O_hydroxyl",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [-5.444452e-11, 1.056440e-11, 1.056777e-10],
                charge_c: 5.268663e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [-1.051734e-10, -7.873551e-11, -4.007338e-11],
                charge_c: 5.268663e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [-5.420368e-11, 9.155090e-11, -5.407293e-11],
                charge_c: 5.268663e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [1.544780e-10, 2.436740e-11, 1.258408e-11],
                charge_c: 2.096297e-01,
                atom_type: "H_hydroxyl",
            },
        ],
        bonds: &[
            GuestBond {
                i: 0,
                j: 1,
                order: 1,
            },
            GuestBond {
                i: 0,
                j: 2,
                order: 1,
            },
            GuestBond {
                i: 0,
                j: 3,
                order: 1,
            },
            GuestBond {
                i: 0,
                j: 4,
                order: 1,
            },
            GuestBond {
                i: 1,
                j: 5,
                order: 1,
            },
        ],
    },
    Guest {
        id: "ethanol",
        label: "Ethanol",
        formula: "C2H6O",
        molar_mass_kg_per_mol: 4.606900000e-02,
        atoms: &[
            GuestAtom {
                element: Element::CARBON,
                position_m: [-1.015116e-10, -6.754788e-12, 1.239509e-11],
                charge_c: -4.183849e-02,
                atom_type: "C_sp3",
            },
            GuestAtom {
                element: Element::CARBON,
                position_m: [4.684111e-11, -2.344919e-11, -1.263506e-11],
                charge_c: 4.022058e-02,
                atom_type: "C_sp3",
            },
            GuestAtom {
                element: Element::OXYGEN,
                position_m: [1.119096e-10, 1.005818e-10, 1.057974e-11],
                charge_c: -3.966637e-01,
                atom_type: "O_hydroxyl",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [-1.433543e-10, 7.005809e-11, -5.348340e-11],
                charge_c: 2.537329e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [-1.198325e-10, 2.614166e-11, 1.152209e-10],
                charge_c: 2.537329e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [-1.548051e-10, -1.007152e-10, -4.633551e-12],
                charge_c: 2.537329e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [8.896257e-11, -9.852550e-11, 5.489220e-11],
                charge_c: 5.606972e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [6.512385e-11, -5.401820e-11, -1.160850e-10],
                charge_c: 5.606972e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [2.066663e-10, 8.668134e-11, -6.250931e-12],
                charge_c: 2.100223e-01,
                atom_type: "H_hydroxyl",
            },
        ],
        bonds: &[
            GuestBond {
                i: 0,
                j: 1,
                order: 1,
            },
            GuestBond {
                i: 1,
                j: 2,
                order: 1,
            },
            GuestBond {
                i: 0,
                j: 3,
                order: 1,
            },
            GuestBond {
                i: 0,
                j: 4,
                order: 1,
            },
            GuestBond {
                i: 0,
                j: 5,
                order: 1,
            },
            GuestBond {
                i: 1,
                j: 6,
                order: 1,
            },
            GuestBond {
                i: 1,
                j: 7,
                order: 1,
            },
            GuestBond {
                i: 2,
                j: 8,
                order: 1,
            },
        ],
    },
    Guest {
        id: "dimethyl_ether",
        label: "Dimethyl ether",
        formula: "C2H6O",
        molar_mass_kg_per_mol: 4.606900000e-02,
        atoms: &[
            GuestAtom {
                element: Element::CARBON,
                position_m: [-1.171631e-10, -9.717365e-12, -1.701194e-12],
                charge_c: 3.507837e-02,
                atom_type: "C_sp3",
            },
            GuestAtom {
                element: Element::OXYGEN,
                position_m: [-4.524330e-12, 3.898404e-12, 8.388013e-11],
                charge_c: -3.879283e-01,
                atom_type: "O_ether",
            },
            GuestAtom {
                element: Element::CARBON,
                position_m: [1.167159e-10, 1.010272e-11, 9.992625e-12],
                charge_c: 3.507837e-02,
                atom_type: "C_sp3",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [-2.071764e-10, -1.396475e-11, 6.020554e-11],
                charge_c: 5.296193e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [-1.249129e-10, 7.637142e-11, -6.881232e-11],
                charge_c: 5.296193e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [-1.102309e-10, -1.022437e-10, -5.971905e-11],
                charge_c: 5.296193e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [1.995847e-10, 2.050626e-11, 8.054336e-11],
                charge_c: 5.296193e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [1.311946e-10, -8.178406e-11, -4.764790e-11],
                charge_c: 5.296193e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [1.165125e-10, 9.683102e-11, -5.674119e-11],
                charge_c: 5.296193e-02,
                atom_type: "H_c",
            },
        ],
        bonds: &[
            GuestBond {
                i: 0,
                j: 1,
                order: 1,
            },
            GuestBond {
                i: 1,
                j: 2,
                order: 1,
            },
            GuestBond {
                i: 0,
                j: 3,
                order: 1,
            },
            GuestBond {
                i: 0,
                j: 4,
                order: 1,
            },
            GuestBond {
                i: 0,
                j: 5,
                order: 1,
            },
            GuestBond {
                i: 2,
                j: 6,
                order: 1,
            },
            GuestBond {
                i: 2,
                j: 7,
                order: 1,
            },
            GuestBond {
                i: 2,
                j: 8,
                order: 1,
            },
        ],
    },
    Guest {
        id: "benzene",
        label: "Benzene",
        formula: "C6H6",
        molar_mass_kg_per_mol: 7.811400000e-02,
        atoms: &[
            GuestAtom {
                element: Element::CARBON,
                position_m: [-1.640512e-11, 1.385102e-10, 1.090394e-12],
                charge_c: -6.226857e-02,
                atom_type: "C_ar",
            },
            GuestAtom {
                element: Element::CARBON,
                position_m: [-1.281596e-10, 5.504828e-11, 4.236854e-13],
                charge_c: -6.226857e-02,
                atom_type: "C_ar",
            },
            GuestAtom {
                element: Element::CARBON,
                position_m: [-1.117545e-10, -8.346192e-11, -6.678807e-13],
                charge_c: -6.226857e-02,
                atom_type: "C_ar",
            },
            GuestAtom {
                element: Element::CARBON,
                position_m: [1.640514e-11, -1.385102e-10, -1.092992e-12],
                charge_c: -6.226857e-02,
                atom_type: "C_ar",
            },
            GuestAtom {
                element: Element::CARBON,
                position_m: [1.281597e-10, -5.504833e-11, -4.261183e-13],
                charge_c: -6.226857e-02,
                atom_type: "C_ar",
            },
            GuestAtom {
                element: Element::CARBON,
                position_m: [1.117545e-10, 8.346189e-11, 6.654150e-13],
                charge_c: -6.226857e-02,
                atom_type: "C_ar",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [-2.918674e-11, 2.464255e-10, 1.943524e-12],
                charge_c: 6.226857e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [-2.280106e-10, 9.793746e-11, 7.571507e-13],
                charge_c: 6.226857e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [-1.988242e-10, -1.484883e-10, -1.184029e-12],
                charge_c: 6.226857e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [2.918669e-11, -2.464257e-10, -1.940973e-12],
                charge_c: 6.226857e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [2.280108e-10, -9.793731e-11, -7.549084e-13],
                charge_c: 6.226857e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [1.988240e-10, 1.484884e-10, 1.186732e-12],
                charge_c: 6.226857e-02,
                atom_type: "H_c",
            },
        ],
        bonds: &[
            GuestBond {
                i: 0,
                j: 1,
                order: 4,
            },
            GuestBond {
                i: 1,
                j: 2,
                order: 4,
            },
            GuestBond {
                i: 2,
                j: 3,
                order: 4,
            },
            GuestBond {
                i: 3,
                j: 4,
                order: 4,
            },
            GuestBond {
                i: 4,
                j: 5,
                order: 4,
            },
            GuestBond {
                i: 5,
                j: 0,
                order: 4,
            },
            GuestBond {
                i: 0,
                j: 6,
                order: 1,
            },
            GuestBond {
                i: 1,
                j: 7,
                order: 1,
            },
            GuestBond {
                i: 2,
                j: 8,
                order: 1,
            },
            GuestBond {
                i: 3,
                j: 9,
                order: 1,
            },
            GuestBond {
                i: 4,
                j: 10,
                order: 1,
            },
            GuestBond {
                i: 5,
                j: 11,
                order: 1,
            },
        ],
    },
    Guest {
        id: "cyclohexane",
        label: "Cyclohexane",
        formula: "C6H12",
        molar_mass_kg_per_mol: 8.416200000e-02,
        atoms: &[
            GuestAtom {
                element: Element::CARBON,
                position_m: [-8.594197e-11, -1.100980e-10, 4.754066e-11],
                charge_c: -5.330597e-02,
                atom_type: "C_sp3",
            },
            GuestAtom {
                element: Element::CARBON,
                position_m: [-1.426297e-10, 2.955888e-11, 2.346785e-11],
                charge_c: -5.330597e-02,
                atom_type: "C_sp3",
            },
            GuestAtom {
                element: Element::CARBON,
                position_m: [-3.585316e-11, 1.373145e-10, 4.033580e-11],
                charge_c: -5.330597e-02,
                atom_type: "C_sp3",
            },
            GuestAtom {
                element: Element::CARBON,
                position_m: [8.594204e-11, 1.100980e-10, -4.754041e-11],
                charge_c: -5.330597e-02,
                atom_type: "C_sp3",
            },
            GuestAtom {
                element: Element::CARBON,
                position_m: [1.426297e-10, -2.955896e-11, -2.346794e-11],
                charge_c: -5.330597e-02,
                atom_type: "C_sp3",
            },
            GuestAtom {
                element: Element::CARBON,
                position_m: [3.585306e-11, -1.373146e-10, -4.033579e-11],
                charge_c: -5.330597e-02,
                atom_type: "C_sp3",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [-5.773282e-11, -1.202721e-10, 1.530622e-10],
                charge_c: 2.665299e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [-1.632393e-10, -1.852061e-10, 2.775054e-11],
                charge_c: 2.665299e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [-1.842040e-10, 3.502364e-11, -7.790169e-11],
                charge_c: 2.665299e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [-2.252684e-10, 4.831245e-11, 9.294387e-11],
                charge_c: 2.665299e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [-4.715289e-12, 1.416074e-10, 1.454361e-10],
                charge_c: 2.665299e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [-7.810241e-11, 2.353255e-10, 1.550411e-11],
                charge_c: 2.665299e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [5.773300e-11, 1.202724e-10, -1.530620e-10],
                charge_c: 2.665299e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [1.632394e-10, 1.852061e-10, -2.775005e-11],
                charge_c: 2.665299e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [1.842045e-10, -3.502389e-11, 7.790138e-11],
                charge_c: 2.665299e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [2.252682e-10, -4.831224e-11, -9.294427e-11],
                charge_c: 2.665299e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [4.714976e-12, -1.416073e-10, -1.454360e-10],
                charge_c: 2.665299e-02,
                atom_type: "H_c",
            },
            GuestAtom {
                element: Element::HYDROGEN,
                position_m: [7.810218e-11, -2.353256e-10, -1.550437e-11],
                charge_c: 2.665299e-02,
                atom_type: "H_c",
            },
        ],
        bonds: &[
            GuestBond {
                i: 0,
                j: 1,
                order: 1,
            },
            GuestBond {
                i: 1,
                j: 2,
                order: 1,
            },
            GuestBond {
                i: 2,
                j: 3,
                order: 1,
            },
            GuestBond {
                i: 3,
                j: 4,
                order: 1,
            },
            GuestBond {
                i: 4,
                j: 5,
                order: 1,
            },
            GuestBond {
                i: 5,
                j: 0,
                order: 1,
            },
            GuestBond {
                i: 0,
                j: 6,
                order: 1,
            },
            GuestBond {
                i: 0,
                j: 7,
                order: 1,
            },
            GuestBond {
                i: 1,
                j: 8,
                order: 1,
            },
            GuestBond {
                i: 1,
                j: 9,
                order: 1,
            },
            GuestBond {
                i: 2,
                j: 10,
                order: 1,
            },
            GuestBond {
                i: 2,
                j: 11,
                order: 1,
            },
            GuestBond {
                i: 3,
                j: 12,
                order: 1,
            },
            GuestBond {
                i: 3,
                j: 13,
                order: 1,
            },
            GuestBond {
                i: 4,
                j: 14,
                order: 1,
            },
            GuestBond {
                i: 4,
                j: 15,
                order: 1,
            },
            GuestBond {
                i: 5,
                j: 16,
                order: 1,
            },
            GuestBond {
                i: 5,
                j: 17,
                order: 1,
            },
        ],
    },
];
