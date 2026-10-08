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
| Porphine mount | ring-to-wall gap; tether | CAD build and valence relaxation | gap 0.55 nm; four 4-bond zigzag tethers, relaxed to 1.42-1.55 A, no overlap | flexible, not rigid; reaction not modelled | 0113 |
| Tank drain channel (CAD) | hydroxyl-lined radial channel | part build and probe test | 72 lining oxygens; free radius 2.3 A to centres, about 1.2 A accessible | flow and leak from ADR-0098; not re-run | 0114 |
| Module pocket water (as built) | water carried into the tank | OpenMM, SPC/E, 2 ns, as-built slot + bore | about 1-2 waters per loaded pocket; iron never wetted; carry is negligible | open cut at the arm back, frozen solid, UFF walls | 0112 |
| Squeezing the water out | rod raised under the O2 | same MD | 3.6-3.7 waters remain | as above | 0096 |
| Gas chamber (option A, rejected) | slot wetting | same MD, 1-2 nm slots | every dry slot floods | UFF carbon-water well not fitted | 0095 |
| Drain: size selection | rigid scan | UFF, all orientations | no size window: 3.6 A blocks both, 3.8 A passes both | lattice gives discrete sizes | 0097 |
| Drain: water flow | pressure-driven MD | 2 ns at 0, 200, 500 MPa | about 9e9 water/s per channel at 1000 atm; one channel suffices | force-driven pressure model; the box held a vapour bubble (corrected run: 12-13.5 waters/ns at 200 MPa); frozen lining only, a mobile lining is untested for this channel | 0098, 0119 |
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
| Hull wall and seals | relaxed CAD section; rigid-probe flood fill of the free space | valence relax; `scripts/seal_audit.py` | rotor in a hull-wall slab, seals at H-H contact; O2 probe (1.73 A) finds no path blood to tank, also at 0.9x radii; water (1.40 A) passes the rim seal, closed at 1.50 A | static rigid-sphere test; leak rates need MD; the as-built drain bore does not pass a rigid water probe | 0117, 0118 |
| Drain as built in the hull wall | MD water transit, frozen and mobile hydroxyl lining | OpenMM SPC/E, constant force | no net water at 100-200 MPa in 2 ns: frozen lining traps fixed waters, mobile lining fills the bore; the drain is unproven | the tested 3.8 A channel conducts only with a frozen lining; design scan open (M28-38) | 0119 |
| Whole device, real size | lumped system model from the atomistic parameters | `respirocyte_system.py`, backward Euler | central: 3 stages fill 1000 atm in 111 s, power-limited (duty 0.5 at 7.3e-12 W), 5.5e-19 J per delivered O2, 48 % of the hull; pessimistic chemistry stalls near 440 atm at any stage count | pocket-to-pocket scatter of the free-gas carry; many assumed values (gas binding, chamber volume, plunger work) | 0124 |
| Whole device, sized by the circulation | loop time, capillary windows, red-cell target, diameter-pressure sweep | `respirocyte_circulation.py` | optimum pressure is the 10-40 atm cushion (not 1000 atm), single stage, rotor speed about 1000 rev/s at every size; nanocrit 0.4-14 % at peak exercise, far below the 45 % hematocrit; the limit is energy, 650 W chemical at exercise with the cam (342 W reversible, about 75 W with a passive release) | the power density, the cam work and the device count at small diameters are assumed or unbounded | 0125 |
| Whole device in the real circulation | lung/tissue cycle, tank over repeated loops | `cycle_steady_state`, backward Euler | the 0.75 s lung transit charges the tank to only 3.7 atm, not 1000, and the discharge empties it each loop, so it never accumulates; the device is a flow-through pump near ambient; 1.28e7 O2 per loop; the 1 um device needs 5.5 % nanocrit at rest, 16 % at exercise; the tank drains below the 10 atm cushion; the fill is power-limited by the cascade banks | no tissue-demand feedback; design rotor speeds | 0132 |
| Top-stage carry, every pocket | released and loaded pocket MD at 600-1000 atm, 600 ps | `scripts/top_stage_sweep.py`, `pocket_carry_md.py` | the released pocket carries 0.30 O2 at the top stage, so a single stage would pump at 1000 atm (s_back 0.060 +- 0.013, s_fwd 0.065 +- 0.022); but the load is not converged (0.04 at 100 ps, 0.38 at 600 ps, 1.5-2.5 in the earlier runs), so the top-stage stall lies between about 500 and 1000 atm | under-equilibrated pocket at 600 ps; the far-density slab is depleted by the rotor | 0134 |
| Overall design: carrier vs pump | passive Hill carrier vs the active pump, on nanocrit and metabolic power | `respirocyte_design.py` | the carrier model reproduces the real red cell (free hemoglobin needs a 12.3 % volume fraction); a heme carrier with a red-cell-like site meets the flux at 1-3 % nanocrit and no engine, and beats the pump on both counts; the pump is only needed for storage or against a gradient | the carrier p50 is assumed tunable and the scaffold is assumed large enough to stay in blood | 0126 |
| Device specs, counts and lifetimes | device count, mass, O2 per volume, power, site lifetime | `respirocyte_specs.py` | the carrier is 4-100x denser in O2 than a red cell, so it carries the flux at 1-3 % by volume; the pump needs 1e14-1e15 machines, the carrier 1e20-1e21 molecules; the site lifetime spans 12 min to 800 days across the references | the autoxidation rate of a synthetic pocket is unmeasured (the M28-22 risk) | 0127 |

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

## CAD model

`respirocyte_module` in the viewer (Machines): the rotor with 12 heme
pockets, 12 [111] ligand arms and the housing, relaxed with the valence
field and audited (ADR-0111). Hide the housing and set the clip to 0.5 to
see the pockets in section. The attachment is four covalent tethers per site; the cam, repair donor
and tank are not in the model.

## How far to trust this

Every number stacks model chemistry: proxy or idealized sites, UFF walls,
rigid parts, two density functionals that disagree on absolute binding by
63 kJ/mol. Each link carries an uncertainty of tens of kJ/mol or a factor of
a few. The chain is internally consistent and no link fails outright, so it
earns a better model. It does not show that the device works.
