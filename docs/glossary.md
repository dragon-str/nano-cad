# Glossary

Terms used in this repository. The level names L0 to L4 follow
`ARCHITECTURE.md`.

## Data model

- **Atom** — The smallest component of a part. An atom has an element, a
  position, a charge, and a type. The model stores atoms in struct-of-arrays
  form.
- **Bond** — A connection between two atoms. A bond has an order and a type.
  Topology uses atom indices, not pointers.
- **Part** — A named collection of atoms and bonds, plus a material and
  metadata. A part is the unit that the generators build.
- **Assembly** — Parts placed by transforms, with constraints and mates. An
  assembly is the L2 arrangement of parts.
- **Jig** — A constraint that acts on a simulation. The jig types are anchor,
  motor, spring, thermostat, rigid body, and joint.
- **Constraint** — A rule that limits the motion of atoms or bodies. A jig
  applies a constraint.

## Jig types

- **Anchor** — A jig that fixes an atom or a group of atoms in place.
- **Motor** — A jig that drives motion. A motor applies a constant angular
  velocity or a torque.
- **Spring** — A jig that applies a spring force between atoms or bodies.
- **Thermostat** — A jig that holds the temperature at a set value. Langevin
  and Berendsen are the planned types.
- **Rigid body** — A body whose shape does not change. At L2 the engine moves a
  rigid body by its position and its orientation.
- **Joint** — A constraint between rigid bodies. The joint types are revolute,
  prismatic, fixed, and gear coupling.

## Model levels

- **L0** — The quantum adapter. It reaches PySCF or xTB through ASE.
- **L1** — The atomistic engine. This is our Rust engine for parts, forces, and
  jigs.
- **L2** — The device layer. It holds rigid bodies, joints, and the parameter
  handoff.
- **L3** — The continuum adapter. It reaches external FEA and flow solvers.
- **L4** — The system layer. It holds lumped ODE models with CellML or SBML
  interop.

## Engine terms

- **Pair list** — A flat list of 32-bit atom index pairs that are close enough
  to interact. The engine reads the pair list for every force evaluation.
- **Skin** — A Verlet skin. This is a buffer distance around the force cutoff.
  The pair list stays valid until an atom moves more than half the skin.

## Formats

- **NCZ** — The native document format. It is a versioned zip of JSON metadata
  and binary arrays, little endian, with a schema version in the header.
- **MMP** — The NanoEngineer file format. nano-cad imports and exports MMP for
  interop with the old tree.
