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
4. Each pocket also carries about four water molecules. A water-return
   rotor takes the water back to the blood; every pocket returns more water
   than it brings, so the chamber stays dry (M28-31). The drain channel was
   tested and failed (ADR-0119, 0120); the tank water path is not modelled.
5. The tank is built gas-filled. The 10 atm cushion is the initial state,
   not a margin the circulation keeps (ADR-0132).
6. A one-electron donor on the housing, on the arc where each site is empty,
   reduces any oxidized Fe(III) back to Fe(II) as the site passes.

The carrier is the baseline. It binds oxygen and releases it by the p50;
it has no tank, no engine and no 1000 atm state. The pump columns in this
table are the storage branch of the study, and every tank or fill number
applies to that branch only (ADR-0126, 0132).

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
| Hull wall and seals | relaxed CAD section; rigid-probe flood fill of the free space | valence relax; `scripts/seal_audit.py` | rotor in a hull-wall slab, seals at H-H contact; O2 probe (1.73 A) finds no path blood to tank, also at 0.9x radii; water (1.40 A) passes the rim seal, closed at 1.50 A | static rigid-sphere test; leak rates by a barrier model (ADR-0148), MD/PMF open; the as-built drain bore does not pass a rigid water probe | 0117, 0118 |
| Drain as built in the hull wall | MD water transit, frozen and mobile hydroxyl lining | OpenMM SPC/E, constant force | no net water at 100-200 MPa in 2 ns: frozen lining traps fixed waters, mobile lining fills the bore; the drain is unproven | the tested 3.8 A channel conducts only with a frozen lining; design scan open (M28-38) | 0119 |
| Whole device, real size | lumped system model from the atomistic parameters | `respirocyte_system.py`, backward Euler | central: 3 stages fill 1000 atm in 97 s, power-limited (duty 0.73 at 7.3e-12 W), 4.9e-19 J per delivered O2, 48.1 % of the hull area and 3.2 % of the free interior; stage stalls 829, 1912, 2824 atm; pessimistic chemistry stalls at 1439 atm at three stages, limited by the binding rate (500 /s), not by carry | pocket-to-pocket scatter of the free-gas carry; many assumed values (gas binding, chamber volume, plunger work); no tank-water bookkeeping unless a seal leak is set | 0124, 0149 |
| Whole device, sized by the circulation | loop time, capillary windows, red-cell target, diameter-pressure sweep | `respirocyte_circulation.py` | optimum pressure is the 10-40 atm cushion (not 1000 atm), single stage, rotor speed about 1000 rev/s at every size; nanocrit 0.4-14 % at peak exercise, far below the 45 % hematocrit; the limit is energy, 650 W chemical at exercise with the cam (342 W reversible, about 75 W with a passive release) | the power density, the cam work and the device count at small diameters are assumed or unbounded | 0125 |
| Whole device in the real circulation | lung/tissue cycle, tank over repeated loops | `cycle_steady_state`, backward Euler | the 0.75 s lung transit charges the tank to only 3.7 atm, not 1000, and the discharge empties it each loop, so it never accumulates; the device is a flow-through pump near ambient; 1.28e7 O2 per loop; the 1 um device needs 5.5 % nanocrit at rest, 16 % at exercise; the tank drains below the 10 atm cushion; the fill is power-limited by the cascade banks | no tissue-demand feedback; design rotor speeds | 0132 |
| Top-stage carry, every pocket | released and loaded pocket MD at 600-1000 atm, 400 ps, bound O2 removed, tall box, empty and prefilled starts | `scripts/top_stage_sweep.py`, `pocket_carry_md.py` | the corrected released carry is s_back 0.076 +- 0.009, about 1.06 O2 at 1000 atm, so a single stage from blood stalls near 920 atm; the top stage still pumps, net 0.57 per pass (forward 0.63, back 1.06), but the margin is thin and one pocket of eight has a negative net | the earlier ADR-0134 carry (0.30) measured a pocket that still held its bound O2; wide per-pocket scatter, one seed | 0138, 0134 |
| Overall design: carrier vs pump | passive Hill carrier vs the active pump, on nanocrit and metabolic power | `respirocyte_design.py` | the carrier model reproduces the real red cell (free hemoglobin needs a 12.3 % volume fraction); a heme carrier with a red-cell-like site meets the flux at 1-3 % nanocrit and no engine, and beats the pump on both counts; the pump is only needed for storage or against a gradient | the carrier p50 is assumed tunable and the scaffold is assumed large enough to stay in blood | 0126 |
| Device specs, counts and lifetimes | device count, mass, O2 per volume, power, site lifetime | `respirocyte_specs.py` | the carrier is 4-100x denser in O2 than a red cell, so it carries the flux at 1-3 % by volume; the pump needs 1e14-1e15 machines, the carrier 1e20-1e21 molecules; the site lifetime comes from a sourced reference table (PubMed ids) and a two-channel mechanistic estimate of the pocket: the central case is 31 h, so a carrier without repair needs about 30 g of new carrier each day | no measured heme carrier reaches the 8-day target; the pocket factors f_ox and f_bi are assumed, and the guard's effect is untested (the M28-22 risk) | 0136, 0127 |
| Porphine mount to the diamond wall | build audit and an independent force-field hold test | `scripts/mount_uff.py`, `scripts/tether_strain.py`, `crates/parts/src/heme_rotor.rs` | the old build bonded each tether to a wall carbon without removing its hydrogen, giving five-bonded carbons and 90 close contacts; the fix removes the anchor hydrogen and uses a rigid ethynyl link (meso-C#C-C). The corrected build has 36 close contacts and a stiff ring (tilt rms 0.56 deg at 310 K); with the tethers cut the ring is soft (1.5 deg), so the tethers hold it | the ethynyl link is a design, not a synthesis; UFF is crude for a metal-porphyrin | 0137 |
| Self-sustaining site repair and lifespan | H2O2 race between the heme and the catalyst; electron budget | `respirocyte_repair.py` | the site needs no new material: the iron is re-reduced by a donor (16,000-fold margin). The ring is damaged by H2O2 (the superoxide converts almost completely to it), so the rate is bounded by the measured heme + H2O2 reaction (3.6e3 /M/s). The optimized pocket (100-fold suppression) with a Mn-porphyrin mimic (1e5 /M/s) gives a 10-year ring lifetime, so the site outlives the 3.2-year rotor; with a catalase it gives about 1000 years | the heme + H2O2 rate is for myoglobin, not this pocket, and is an upper bound; the Mn-porphyrin catalase rate is an estimate; the 100-fold pocket suppression is a design goal | 0139, 0140 |
| Donor and catalyst in the site CAD | build audit of the self-sustaining groups | `crates/parts/src/heme_rotor.rs` | `CarrierSiteGenerator` builds a proximal pocket below the ring with a donor hydroxyl on its floor (Fe-O 3.61 A) and a manganese center on the pocket wall (Fe-Mn 6.37 A); every carbon stays four-bonded and the tether angles stay within 95-140 deg | the catalyst is a design, not a synthesis; the Mn center is placed by a build rule, not by energy | 0139, 0140 |
| Donor and catalyst in the oxygen rotor | build audit of the self-sustaining groups on each rim pocket | `crates/parts/src/heme_rotor.rs` | the device's `HemeRotorGenerator` builds the catalyst on the distal pocket wall (Fe-Mn 6.16-6.86 A, three bonds, one open site) and the donor on the nearest proximal scaffold carbon (Fe-O 7.14-7.63 A) on every one of the eight pockets; no carbon gains a fifth bond | the arm collar occupies the proximal axis, so the rotor donor sits at about 7.4 A, twice the free site; the distance is not a limit (row below), but the donor is static, so its potential must separate repair from back-oxidation | 0141 |
| Donor distance against the repair rate | Moser-Dutton tunnelling from the ADR-0108 points | `python/nanocad/respirocyte_repair.py` | `donor_repair_rate` (beta 1.38 /A) gives 1.1e10 /s at the free site (3.61 A) and 5.6e7 /s at the rotor (7.4 A), both far above the model's 22 A worst case of 0.1 /s; the margins over autoxidation are 1.7e15 and 9.0e12 | Moser-Dutton is distance only; the driving force is assumed; the donor potential window is uncomputed | 0142 |
| Five-coordinate domed iron in an unloaded pocket | build audit under the 0.3 A dome | `crates/parts/src/heme_rotor.rs` | with `iron_dome_m` set, the iron sits 0.30 A out of the N4 plane, the four Fe-N bonds are 2.02 A, the loaded pockets stay in-plane, and the four meso tethers are unchanged (no fifth bond) | 0.3 A is the deoxy-myoglobin value, not computed here; the ring is rigid, so only the iron moves | 0143 |
| Discharge rotor | cam-phase swap of the loading rotor | `crates/parts/src/heme_rotor.rs` | `build_cam` takes any cam phase pair; `DischargeRotorGenerator` (id `discharge_rotor`) swaps the loading phases, so the arm goes in on the tank side (load at tank pressure) and out on the blood side (release at the tissues). The loaded pockets are complementary: loading 5 of 8, discharge 3 of 8, every pocket flipping. The disc body and pocket azimuths are identical | geometric only; the discharge loading rate against the tank pressure and the seal leak at that pressure are unmodelled | 0144 |
| 62. Repair-donor redox window | QM couples from the redox scan, against a literature window | `scripts/heme_switch_scan.py --redox Fe`, `respirocyte_repair.py` | the model's empty-site couple E(Fe(III))-E(Fe(II)) is +5.95 eV (gas phase), the bare ionization energy, not an aqueous potential; a rigid C-PCM point cannot pay the solvent reorganization that takes the real heme couple to near zero, so the method cannot place a donor between the empty-site and bound-O2 couples. The real window is narrow, about 0.21 V (O2/O2- about -0.16 V to a heme couple about +0.05 V, literature estimates) | the absolute QM energies are unphysical; the 0.21 V window is literature, not computed; the empty-arc donor station is a design (M29-19) | 0145, 0108 |
| 63. Gated repair donor on the housing | build audit of the empty-arc station, with the tunnelling rate | `crates/parts/src/heme_rotor.rs`, `respirocyte_repair.py` | `HullWallSpec.donor_station` places one hydroxyl on the top end plate at the midpoint of the arm-out arc (150 deg for the loading rotor), 14.98 A from the iron; Moser-Dutton tunnelling there is 1.6e3 /s, a 2.6e8 margin over the autoxidation. The site meets the station only while it is empty, so the bound O2 never sees the donor. The module builds the rotor donor off and the station on | the tunnelling driving force is assumed (ADR-0108); the station is one group at one azimuth, so the repair window is not the whole empty arc; the Moser-Dutton beta 1.38 /A is a protein value, not a diamond or vacuum value | 0147 |
| 64. Seal leak rates | Hertz-Knudsen permeation from the static seal barriers, at 310 K | `python/nanocad/respirocyte_leak.py` | the O2 rotor O2 barrier (146 kJ/mol, 57 kT) gives 4e-13 /s per nm^2, so the product does not leak; the O2 rotor water barrier (57 kJ/mol, 22 kT) gives 1.3e3 /s, about 3e15x the O2 leak; the water rotor is tighter (O2 78 kT, water 58 kT). The water leak is the soft one and the water-return rotor must match it | the barrier is static (no entropy term) and bulk impingement over-states a narrow channel, so these are upper bounds; a PMF is the upgrade and the balance against the water rotor is open | 0148 |

## What passes, what fails, what is open

Passes as simulated estimates: selective loading, release by a strain
switch, water return by the plunger rotor, idle leak, mean power.

The baseline is the passive carrier (ADR-0126): it meets the flux at 1-3 %
by volume with no engine. The active pump is the storage branch of the
study. Its peak power at a 10 s fill is 1.7-2.4x the Freitas maximum, and a
cam that returns more of its energy would close that gap (ADR-0100, 0106).

The fill time cannot be shortened to fix a power gap. A lung transit lasts
about 0.75 s, so the tank reaches only a few atm per pass (ADR-0132); 97 s
is the time to 1000 atm in the lumped model, not a transit time.

Fails: a cam cannot force the iron to dome, by pulling its axial ligand
(ADR-0104) or by squeezing the ring (ADR-0105). Withdrawing the axial
ligand works instead (ADR-0106).

Not modelled:

- How the porphine attaches to diamond.
- Relaxed geometries. Every quantum-chemistry number uses a rigid geometry.
- Irreversible damage of the porphyrin ring by superoxide and peroxide.
  This sets the device lifetime.
- The donor chemistry of the repair station. The donor potential window is
  not computable with the current QM (ADR-0145), so the donor is gated by
  geometry: a station on the housing on the empty arc (ADR-0147). The
  tunnelling driving force is assumed.
- Nitrogen build-up in the tank (a small vent is proposed, not tested).
- Fouling, immune response, manufacturing.

## CAD model

`respirocyte_module` in the viewer (Machines): the rotor with 8 heme
pockets, 8 [111] ligand arms and the housing, relaxed with the valence
field and audited (ADR-0111). Hide the housing and set the clip to 0.5 to
see the pockets in section. The attachment is four covalent tethers per site; the repair donor
is a station on the housing (ADR-0147) and the tank is a separate part.

## How far to trust this

Every number stacks model chemistry: proxy or idealized sites, UFF walls,
rigid parts, two density functionals that disagree on absolute binding by
63 kJ/mol. Each link carries an uncertainty of tens of kJ/mol or a factor of
a few. The chain is internally consistent and no link fails outright, so it
earns a better model. It does not show that the device works.
