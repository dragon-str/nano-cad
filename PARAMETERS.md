# PARAMETERS.md — The Multiscale Parameter Schema

This schema is the central contract. Every level reads it and writes it. It is
the reason the levels stay decoupled. Version it carefully. A change needs a
version bump and an ADR.

## Purpose

A lower level computes or measures properties. A higher level consumes them.
The schema carries the value, the unit, the provenance, and the uncertainty.

## Rules

1. Every value has a unit. Internally SI.
2. Every value has a provenance record. Who measured or computed it, how, and
   when.
3. Every value has an uncertainty. If unknown, write `null`, not zero.
4. Every record has a schema version.
5. A record names the method that produced it: `md`, `dft`, `experiment`,
   `fit`, or `estimate`.
6. A record states its validation status: `unverified`, `cross-checked`, or
   `certified`.

## Record types

### PartRecord (L1 output, L2 input)

```
schema: "nanocad.param.part"
version: 1
part_id: string
geometry_ref: string        # hash of the NCZ part
material: string
atoms: integer
mass_kg: Quantity
inertia_kg_m2: [Quantity, Quantity, Quantity]
elastic_modulus_Pa: Quantity
shear_modulus_Pa: Quantity
poisson_ratio: Quantity
failure_stress_Pa: Quantity
friction_coefficient: Quantity
thermal_conductivity_W_m_K: Quantity
specific_heat_J_kg_K: Quantity
method: enum
validation: enum
provenance: Provenance
```

### JointRecord (L2 output, L3 input)

Revolute or prismatic joints carry stiffness, damping, friction, and backlash.
Each is a Quantity with provenance.

### BodyRecord (L2 output, L3 input)

Mass, inertia, center of mass, and the joints that attach the body.

### SystemRecord (L4)

Lumped parameters: capacity, rate constants, diffusion coefficients, and
volumes. Units explicit.

### Provenance

```
Provenance:
  source: string            # tool, person, or dataset
  method: enum              # md | dft | experiment | fit | estimate
  code_version: string
  force_field: string|null
  timestamp: ISO-8601
  uncertainty: Quantity|null
  validation: enum
  notes: string
```

## The extraction pipeline (L1 to L2)

1. Load the part. Relax it to a minimum.
2. Measure stiffness: apply small strains, record the stress response, fit the
   modulus. Report the fit uncertainty.
3. Measure friction: slide one part over another at a set load and speed.
   Record the mean and the spread.
4. Measure failure: increase the load until bonds break. Record the stress.
5. Measure thermal properties: run NVT, extract conductivity and specific heat.
6. Write a PartRecord. Attach the provenance and the validation status.
7. Store the record in the parameter library. Version the library.

## Example

```
schema: "nanocad.param.part"
version: 1
part_id: "gear.sun.diamond.2nm.12t"
geometry_ref: "sha256:..."
material: "diamond"
elastic_modulus_Pa: { value: 1.05e12, unit: "Pa", uncertainty: 0.08e12 }
friction_coefficient: { value: 0.05, unit: "1", uncertainty: 0.02 }
method: "md"
validation: "unverified"
provenance: { source: "nanocad-engine", code_version: "0.1.0",
              force_field: "MM4", timestamp: "2026-09-13T00:00:00Z" }
```

## Versioning

- The schema version is an integer. Bump it on any breaking change.
- Keep one migration function per version step.
- Golden files pin the current version.
