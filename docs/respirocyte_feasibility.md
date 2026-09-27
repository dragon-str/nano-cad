# Respirocyte oxygen chain: feasibility summary

Status: simulated estimates, 2026-09-28. Nothing here is built or validated.
Each line names the test, the method, the result and the main uncertainty.
The decision records (ADR) in `DECISIONS.md` hold the details.

## The design under test

1. A rotor with small pockets on its rim takes oxygen from blood.
2. Each pocket holds one heme-like iron site (Fe-porphine with an axial
   ligand). The site binds oxygen when the iron is flat.
3. At the outlet a cam pulls the axial ligand. The iron domes 0.3 A out of
   the ring plane, the bond weakens by about 50 kJ/mol, and the oxygen
   leaves into the tank.
4. Each pocket also carries about four water molecules. The tank wall has
   single-file drain channels lined with hydroxyl groups. The tank pressure
   pushes the water out through them; the oxygen stays in.
5. The tank is built gas-filled and never falls below a 10 atm cushion.

## Results

| Link | Test | Method | Result | Main uncertainty | ADR |
|---|---|---|---|---|---|
| Gear pump (rejected) | release at the mesh | geometry | the other gear never reaches the bound O2; each gap carries about 100 waters | none: the design fails | 0089, 0090 |
| Site selects O2 over water | competition at alveolar pressures | Fe(NH3)4 proxy, PBE; Langmuir | O2 0.995, water 7e-12 | proxy site, PBE | 0084, 0086, 0089 |
| Pocket size | water in the pocket | OpenMM, SPC/E, 1 ns | 20 waters (7.1 A bore), 3.8 beside bound O2 (5.6 A bore) | UFF walls, frozen solid | 0090, 0091 |
| Squeezing the water out | rod raised under the O2 | same MD | 3.6-3.7 waters remain | as above | 0096 |
| Gas chamber (option A, rejected) | slot wetting | same MD, 1-2 nm slots | every dry slot floods | UFF carbon-water well not fitted | 0095 |
| Drain: size selection | rigid scan | UFF, all orientations | no size window: 3.6 A blocks both, 3.8 A passes both | lattice gives discrete sizes | 0097 |
| Drain: water flow | pressure-driven MD | 2 ns at 0, 200, 500 MPa | about 9e9 water/s per channel at 1000 atm; one channel suffices | force-driven pressure model | 0098 |
| Drain: O2 loss while draining | Widom insertion | MD | about 0.5 % of pumped O2 | noisy insertion average | 0097 |
| Drain: idle leak | umbrella sampling, forward and reverse | 82 windows, WHAM | 1-4 % of a full tank per 60 s cycle | assumed diffusion coefficient; end mismatch 7-9 kJ/mol | 0099 |
| Release switch, proxy | strain scan | Fe(NH3)4, TPSSh and B3LYP | doming weakens the bond by about 100 kJ/mol | proxy, rigid | 0094 |
| Release switch, heme-like site | dome scan | Fe-porphine, PBE and B3LYP | switch 50 (PBE) and 57 (B3LYP) kJ/mol; at least 32 needed | absolute level differs by 63 kJ/mol between methods; anchored to myoglobin P50 | 0101, 0103 |
| Tank start-up | nucleation | classical nucleation theory | needs 850-1700 atm supersaturation; use a 10 atm gas cushion | rate prefactor | 0102 |
| Energy and power | budget against Freitas | arithmetic | cam work 50-57 kJ/mol (minimum 23); peak power 2x Freitas at a 10 s fill; mean power inside his figures | Freitas's own engine figure is below the compression minimum | 0100, 0101 |

## What passes, what fails, what is open

Passes as simulated estimates: selective loading, release by a strain
switch, water removal by a drain, idle leak, tank start-up, mean power.

Fails or exceeds: peak power at a 10 s fill is about 2x the Freitas maximum.
A 20 s fill, or a cam that returns more of its energy, would close it.

Not modelled:

- How the porphine attaches to diamond, and how a cam pulls its lower
  ligand without breaking the attachment.
- Relaxed geometries. Every quantum-chemistry number uses a rigid geometry.
- Autoxidation of Fe(II) to Fe(III), which stops binding (met-hemoglobin).
- Nitrogen build-up in the tank (a small vent is proposed, not tested).
- Fouling, immune response, manufacturing.

## How far to trust this

Every number stacks model chemistry: proxy or idealized sites, UFF walls,
rigid parts, two density functionals that disagree on absolute binding by
63 kJ/mol. Each link carries an uncertainty of tens of kJ/mol or a factor of
a few. The chain is internally consistent and no link fails outright, so it
earns a better model. It does not show that the device works.
