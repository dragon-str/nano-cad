# Three-scale scene export

This document describes the M8-01 scene export. A visualization tool reads one
JSON file and animates the planetary gearbox across three scales.

## Where the code lives

The export is in `crates/jigs/src/scene.rs`. It belongs to `nanocad-jigs`
because that crate owns `PlanetaryAssembly`, the rigid-body device, and the
URDF export.

## The three layers

The scene has three layers. A viewer shows one layer and cross-fades to the
next.

- `atomistic` is the L1 view. It holds every atom with its element and its
  position.
- `device` is the L2 view. It holds every rigid body with its pose and its
  velocities. It also holds the joints, the gear couplings, and the gear
  constraints.
- `coarse` is the handoff view. It holds a gear pitch circle for each gear and
  a bounding cylinder for each body.

## Coordinate system and units

- The frame is right-handed and Cartesian.
- The z axis is the gear axis. The gears lie in the x-y plane.
- Length is SI metres. Angle is SI radians.
- An orientation is a unit quaternion `[w, x, y, z]`. It maps a body-frame
  vector to the world frame.
- The atomistic positions are the world-frame positions at the zero
  configuration. The device poses are the current poses of the rigid-body
  system. The two agree at the zero configuration.

A viewer animates an atomistic body with the device transform of the body
named by `atomistic.atoms[].body`. That index points into `device.bodies`.

## Determinism

The export is deterministic. The body order, the joint order, and the atom
order are the construction order. The code does not iterate a map.

## Public API

```
pub fn build_scene(
    assembly: &PlanetaryAssembly,
    set: &PlanetarySet,
) -> Result<Scene, SceneError>

pub fn scene_to_json(scene: &Scene) -> Result<String, SceneError>
pub fn scene_to_json_pretty(scene: &Scene) -> Result<String, SceneError>
pub fn scene_from_json(json: &str) -> Result<Scene, SceneError>
pub fn write_scene_json(
    scene: &Scene,
    path: impl AsRef<Path>,
    pretty: bool,
) -> Result<(), SceneError>
```

`build_scene` requires the `PlanetarySet` that produced the assembly records.
It compares the embedded design and returns `SceneError::DesignMismatch` when
the two differ.

## Schema

The top level holds the schema name, the version, the units, the frame note,
the derived design, and the three layers.

```json
{
  "schema": "nanocad.scene",
  "version": 1,
  "length_unit": "m",
  "angle_unit": "rad",
  "frame": "right-handed, z is the gear axis, x-y is the gear plane, metres and radians, quaternion [w, x, y, z]",
  "design": {
    "module_m": 5e-10,
    "sun_teeth": 24,
    "planet_teeth": 18,
    "ring_teeth": 60,
    "planet_count": 3,
    "gear_ratio": 3.5,
    "sun_pitch_radius_m": 6e-9,
    "planet_pitch_radius_m": 4.5e-9,
    "ring_pitch_radius_m": 1.5e-8,
    "carrier_radius_m": 1.05e-8
  },
  "atomistic": {
    "atom_count": 1,
    "atoms": [
      {
        "element": "C",
        "atomic_number": 6,
        "position_m": [0.0, 0.0, 0.0],
        "body": 1,
        "name": "sun",
        "part_id": "planetary.sun"
      }
    ]
  },
  "device": {
    "body_count": 7,
    "bodies": [
      {
        "index": 0,
        "name": "ground",
        "role": "ground",
        "part_id": null,
        "position_m": [0.0, 0.0, 0.0],
        "orientation_wxyz": [1.0, 0.0, 0.0, 0.0],
        "mass_kg": 1.0,
        "inertia_diagonal_kg_m2": [1.0, 1.0, 1.0],
        "fixed": true,
        "linear_velocity_m_per_s": [0.0, 0.0, 0.0],
        "angular_velocity_rad_per_s": [0.0, 0.0, 0.0]
      }
    ],
    "joints": [
      {
        "name": "sun_joint",
        "kind": "revolute",
        "body_a": 0,
        "body_b": 1,
        "anchor_world_m": [0.0, 0.0, 0.0],
        "axis_world": [0.0, 0.0, 1.0],
        "axis_a_body": [0.0, 0.0, 1.0],
        "axis_b_body": [0.0, 0.0, 1.0]
      }
    ],
    "gear_couplings": [],
    "gear_constraints": [
      {
        "terms": [
          { "joint": 0, "coefficient_m": 6e-9 },
          { "joint": 1, "coefficient_m": -6e-9 },
          { "joint": 3, "coefficient_m": 4.5e-9 }
        ]
      }
    ]
  },
  "coarse": {
    "body_count": 7,
    "bodies": [
      {
        "index": 1,
        "name": "sun",
        "role": "sun",
        "pitch_circle": {
          "center_m": [0.0, 0.0, 0.0],
          "axis": [0.0, 0.0, 1.0],
          "radius_m": 6e-9
        },
        "bounding_cylinder": {
          "center_m": [0.0, 0.0, 0.0],
          "axis": [0.0, 0.0, 1.0],
          "radius_m": 7e-9,
          "half_length_m": 0.0
        }
      }
    ]
  }
}
```

The example is shortened. A real file holds every atom, every body, and every
joint.

Notes on the schema:

- `element` is the chemical symbol. An unsupported element becomes `Z<n>`,
  where `<n>` is the atomic number.
- `part_id` is `null` for the ground. The other bodies use the part ids
  `planetary.sun`, `planetary.planet.<i>`, `planetary.ring`, and
  `planetary.carrier`.
- `anchor_world_m` is `null` for a prismatic joint. A prismatic joint does not
  expose its anchor.
- `axis_b_body` is `null` for a prismatic joint. It does not expose the value.
- `pitch_circle` is `null` for the ground and for the carrier. The carrier is a
  plate, not a gear.
- `bounding_cylinder` is `null` for the ground. The ground has no atoms.

## Tests

`crates/jigs/src/scene.rs` holds eight tests. They check the layer counts for
the default design, finite positions, the atom count against the generated
parts, the atom positions against the topology, the JSON round trip, the
determinism, the pitch circles, and the file write.

## Limitations

- The atomistic layer is a static snapshot. It has no bonds and no atom type.
- The device layer has no joint limits and no joint names for a joint that the
  assembly does not name. Such a joint falls back to `revolute_<i>`.
- The coarse layer uses the device axis only. A tilted body does not tilt its
  bounding cylinder.
- The export does not read the parameter library. It reads geometry from the
  generated parts and the design.
