//! Writes the default three-scale planetary gearbox scene as JSON.
//!
//! The first argument is the output path. The default is `site/scene.json`.
//! The static viewer in `site/` loads this file. The scene is a schematic
//! snapshot of a simulated design, not a built device.
//!
//! The example also writes two companion files:
//!
//! - `<stem>.data.js` assigns the same scene to `window.NANOCAD_SCENE`. A
//!   browser blocks a `fetch` of a local file from a `file://` page, but it
//!   allows a `<script>` tag, so the viewer auto-loads when you double-click
//!   `site/index.html`.
//! - `<stem>.bonds.json` carries the real bond topology for the video's
//!   atomistic layer. It holds the diamond cubic lattice from the
//!   `nanocad-parts` [`DiamondGenerator`] with its first-shell bonds, and it
//!   holds the topology bonds of the generated [`PlanetarySet`]. The diamond
//!   lattice is the only repository generator whose geometry is true
//!   diamondoid carbon: the `PlanetaryGenerator` emits a skeletal gear
//!   outline, not a diamond network.
//!
//! Run it from the repository root:
//!
//! ```sh
//! cargo run -p nanocad-jigs --example scene_json -- site/scene.json
//! ```

use std::error::Error;
use std::path::{Path, PathBuf};

use serde_json::{json, Value};

use nanocad_jigs::{assemble_planetary, build_scene, scene_to_json_pretty};
use nanocad_model::{Part, Topology};
use nanocad_parts::{
    DiamondGenerator, ParameterSet, PartGenerator, PlanetaryGenerator, PlanetarySet,
};

/// The diamond conventional cubic lattice constant at 300 K, in metres.
///
/// Source: CRC Handbook of Chemistry and Physics; N. W. Ashcroft and N. D.
/// Mermin, Solid State Physics (1976), chapter 4. The first-shell bond length
/// is `a sqrt(3) / 4`, which equals `1.544e-10 m`.
const DIAMOND_LATTICE_CONSTANT_M: f64 = 3.567e-10;

/// The side of the rendered diamond block, in conventional cubic cells.
const DIAMOND_CELLS: f64 = 3.0;

fn write_file(path: &Path, contents: &str) -> Result<(), Box<dyn Error>> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    std::fs::write(path, contents)?;
    Ok(())
}

/// Serializes the bond list of a topology as `[u, v]` index pairs.
fn bond_pairs(topology: &Topology) -> Value {
    let pairs: Vec<Value> = (0..topology.bond_count())
        .filter_map(|index| {
            let u = topology.bond_u(index)?;
            let v = topology.bond_v(index)?;
            Some(json!([u, v]))
        })
        .collect();
    Value::Array(pairs)
}

/// Serializes the atoms of a topology as `{element, position_m}` records.
fn atom_records(topology: &Topology) -> Value {
    let atoms: Vec<Value> = (0..topology.atom_count())
        .filter_map(|index| topology.atom(index))
        .map(|atom| {
            let element = atom.element.symbol().unwrap_or("C").to_owned();
            json!({
                "element": element,
                "position_m": [atom.position_m[0], atom.position_m[1], atom.position_m[2]],
            })
        })
        .collect();
    Value::Array(atoms)
}

/// Builds the companion bond document.
///
/// The `diamond` block is the renderer's atomistic layer. The `planetary`
/// block is the topology of the generated gear set, for provenance. The scene
/// atom order equals the topology index order, so the planetary pairs index
/// the atoms in `<stem>.json` directly.
fn build_bonds_document(set: &PlanetarySet) -> Result<Value, Box<dyn Error>> {
    let diamond_parameters = ParameterSet::new()
        .with("cells_x", DIAMOND_CELLS)
        .with("cells_y", DIAMOND_CELLS)
        .with("cells_z", DIAMOND_CELLS);
    let diamond: Part = DiamondGenerator.generate(&diamond_parameters)?;
    let bond_length_m = DIAMOND_LATTICE_CONSTANT_M * 3.0_f64.sqrt() / 4.0;

    let planetary = &set.part.topology;

    Ok(json!({
        "schema": "nanocad.scene.bonds",
        "version": 1,
        "length_unit": "m",
        "angle_unit": "deg",
        "note": "Real generator output. The diamond block is the atomistic layer. \
                 The planetary block is the skeletal gear topology, for provenance.",
        "reference": {
            "c_c_bond_m": 1.544e-10,
            "tetrahedral_angle_deg": 109.4712,
        },
        "diamond": {
            "generator": DiamondGenerator.id(),
            "lattice_constant_m": DIAMOND_LATTICE_CONSTANT_M,
            "bond_length_m": bond_length_m,
            "atom_count": diamond.topology.atom_count(),
            "bond_count": diamond.topology.bond_count(),
            "atoms": atom_records(&diamond.topology),
            "bonds": bond_pairs(&diamond.topology),
        },
        "planetary": {
            "generator": "planetary",
            "atom_count": planetary.atom_count(),
            "bond_count": planetary.bond_count(),
            "bonds": bond_pairs(planetary),
        },
    }))
}

/// Checks that the scene atom order matches the planetary topology order.
///
/// The companion planetary bonds index the topology directly. The check makes
/// the assumption explicit, so a future reordering fails loudly.
fn check_scene_order(
    scene: &nanocad_jigs::Scene,
    topology: &Topology,
) -> Result<(), Box<dyn Error>> {
    if scene.atomistic.atoms.len() != topology.atom_count() {
        return Err(format!(
            "scene atom count {} does not match topology atom count {}",
            scene.atomistic.atoms.len(),
            topology.atom_count()
        )
        .into());
    }
    for (index, atom) in scene.atomistic.atoms.iter().enumerate() {
        let Some(position_m) = topology.position_m(index) else {
            return Err(format!("topology is missing atom {index}").into());
        };
        if atom.position_m != position_m {
            return Err(format!("scene atom {index} does not match the topology order").into());
        }
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let path = std::env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("site/scene.json"));

    let parameters = ParameterSet::new().with("planet_count", 3.0);
    let assembly = assemble_planetary(&parameters)?;
    let set = PlanetaryGenerator.build(&parameters)?;
    let scene = build_scene(&assembly, &set)?;
    check_scene_order(&scene, &set.part.topology)?;
    let json = scene_to_json_pretty(&scene)?;

    write_file(&path, &format!("{json}\n"))?;
    let companion = path.with_extension("data.js");
    write_file(&companion, &format!("window.NANOCAD_SCENE = {json};\n"))?;

    let bonds = build_bonds_document(&set)?;
    let bonds_path = path.with_extension("bonds.json");
    write_file(
        &bonds_path,
        &format!("{}\n", serde_json::to_string_pretty(&bonds)?),
    )?;

    let diamond = bonds
        .get("diamond")
        .and_then(|value| value.get("atom_count"))
        .and_then(Value::as_u64)
        .unwrap_or(0);
    let diamond_bonds = bonds
        .get("diamond")
        .and_then(|value| value.get("bond_count"))
        .and_then(Value::as_u64)
        .unwrap_or(0);

    println!(
        "wrote {} and {}: schema {} version {} bodies {} joints {} planets {} ratio {} atoms {}",
        path.display(),
        companion.display(),
        scene.schema,
        scene.version,
        scene.device.body_count,
        scene.device.joints.len(),
        scene.design.planet_count,
        scene.design.gear_ratio,
        scene.atomistic.atom_count,
    );
    println!(
        "wrote {}: diamond {} atoms {} bonds; planetary {} atoms {} bonds",
        bonds_path.display(),
        diamond,
        diamond_bonds,
        set.part.topology.atom_count(),
        set.part.topology.bond_count(),
    );
    Ok(())
}
