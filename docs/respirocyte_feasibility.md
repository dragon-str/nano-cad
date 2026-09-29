# Respirocyte oxygen chain: feasibility summary

Status: simulated estimates, 2026-09-28. Nothing here is built or validated.
Each line names the test, the method, the result and the main uncertainty.
The decision records (ADR) in `DECISIONS.md` hold the details.

## The design under test

1. A rotor with small pockets on its rim takes oxygen from blood.
2. Each pocket holds one heme-like iron site (Fe-porphine with an axial
   ligand). The site binds oxygen when the iron is flat.
3. At the outlet a cam withdraws the axial ligand about 2 A from the iron.
   The O2 bond weakens by 42-63 kJ/mol, and the oxygen leaves into the
   tank. Forcing the iron to dome was tested and fails (ADR-0104, 0105).
4. Each pocket also carries about four water molecules. The tank wall has
   single-file drain channels lined with hydroxyl groups. The tank pressure
   pushes the water out through them; the oxygen stays in.
5. The tank is built gas-filled and never falls below a 10 atm cushion.
6. A one-electron donor on the housing, on the arc where each site is empty,
   reduces any oxidized Fe(III) back to Fe(II) as the site passes.

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
| Cam lever through the axial ligand | pull scan | Fe-porphine, PBE | fails: the axial Fe-N bond gives way at about 1.3 nN; doming needs 3.2-6.6 nN | rigid ring | 0104 |
| Cam lever by squeezing the ring | core-size scan | Fe-porphine, PBE | fails: the loaded iron stays flat; the ring absorbs the squeeze | uniform ring scaling | 0105 |
| Release by withdrawing the axial ligand | empty and loaded pull scans | Fe-porphine, PBE and B3LYP | passes: switch 42 (PBE) and 63 (B3LYP) kJ/mol, 32 needed; cam work equals the switch; peak power 1.7-2.4x Freitas | rigid ring, NH3 ligand | 0106, 0107 |
| Oxidation of the iron site | Fe/Co redox scan; repair budget; tunnelling estimate | PBE + C-PCM; Moser-Dutton | keep iron (Co oxidizes more easily in this model); 3e-10 oxidations per O2; a donor station on the empty arc repairs far faster than sites oxidize, at about 1e-20 W | absolute redox energies unphysical; only the Fe-Co difference used | 0108 |
| O2 loading rate | k_on c budget; MD arrival count | myoglobin-like k_on; OpenMM | arrival into the pocket about 3e6/s (MD), so bond formation limits; design basis 2e3/s per site: about 9,700 rotors at about 300 rev/s | bond-formation step not modelled; large seed spread | 0109, 0110 |
| Tank start-up | nucleation | classical nucleation theory | needs 850-1700 atm supersaturation; use a 10 atm gas cushion | rate prefactor | 0102 |
| Energy and power | budget against Freitas | arithmetic | cam work 50-57 kJ/mol (minimum 23); peak power 2x Freitas at a 10 s fill; mean power inside his figures | Freitas's own engine figure is below the compression minimum | 0100, 0101 |

## What passes, what fails, what is open

Passes as simulated estimates: selective loading, release by a strain
switch, water removal by a drain, idle leak, tank start-up, mean power.

Fails or exceeds: peak power at a 10 s fill is 1.7-2.4x the Freitas maximum.
A 20-25 s fill brings it inside.
A 20 s fill, or a cam that returns more of its energy, would close it.

Fails: a cam cannot force the iron to dome, by pulling its axial ligand
(ADR-0104) or by squeezing the ring (ADR-0105). Withdrawing the axial
ligand works instead (ADR-0106).

Not modelled:

- How the porphine attaches to diamond.
- Relaxed geometries. Every quantum-chemistry number uses a rigid geometry.
- Irreversible damage of the porphyrin ring by superoxide and peroxide.
  This sets the device lifetime.
- The donor chemistry of the repair station and its potential.
- Nitrogen build-up in the tank (a small vent is proposed, not tested).
- Fouling, immune response, manufacturing.

## How far to trust this

Every number stacks model chemistry: proxy or idealized sites, UFF walls,
rigid parts, two density functionals that disagree on absolute binding by
63 kJ/mol. Each link carries an uncertainty of tens of kJ/mol or a factor of
a few. The chain is internally consistent and no link fails outright, so it
earns a better model. It does not show that the device works.
