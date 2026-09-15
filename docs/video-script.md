# nano-cad gearbox — 90-second video script and storyboard

This document is the M8-02 deliverable. It is the script and the storyboard
for a 90-second video about the nano-cad planetary gearbox. The recording
itself is a human task. This document contains only the plan.

> **Claim rule.** The narration says "simulated". The narration never claims a
> physical device or a medical validation. Every number on screen cites its
> source in `docs/report.md` or `TASKS.md`.

## Format

| Item | Value |
|---|---|
| Total length | 90 seconds |
| Aspect ratio | 16:9 |
| Frame rate | 30 frames per second |
| Audio | Voice-over plus captions |
| Captions | Burned in, and also a sidecar `.srt` |
| On-screen text | Monospace for commands and output |

## Time budget

The shot list sums to 90 seconds.

| Shot | Start (s) | End (s) | Length (s) |
|---|---|---|---|
| 1. Cold open: the simulated gearbox | 0 | 8 | 8 |
| 2. What this is: the five levels | 8 | 19 | 11 |
| 3. Clean checkout and one command | 19 | 29 | 10 |
| 4. The Rust gates pass | 29 | 40 | 11 |
| 5. The agent generates the gear set | 40 | 52 | 12 |
| 6. The agent assembles the L2 device | 52 | 62 | 10 |
| 7. The agent measures the gear ratio | 62 | 72 | 10 |
| 8. URDF and the three-scale scene | 72 | 82 | 10 |
| 9. Benchmark timings | 82 | 88 | 6 |
| 10. Caveats and closing card | 88 | 90 | 2 |
| **Total** | | | **90** |

## Shot list and narration

### Shot 1 — Cold open: the simulated gearbox (0 to 8 s)

**Visual.** A slow orbit of the three-scale scene. The atomistic layer
cross-fades to the device layer. The gears turn. Overlay: `simulated`.

**Narration.** "This is a planetary gearbox. It is simulated at the atomic
scale and at the device scale. The sun gear turns. The planets turn. The
carrier turns."

**On-screen text.** `simulated`, `sun ratio 3.500000000`.

**Source for the number.** `TASKS.md` M6-06 and `docs/reproduce.md`.

### Shot 2 — What this is: the five levels (8 to 19 s)

**Visual.** An animated diagram of L0 to L1 to L2 to L3 to L4. L0, L3, and L4
are gray and marked "later". L1 and L2 light up.

**Narration.** "nano-cad is a multiscale design tool. It has five model
levels. Only L1, the atomistic engine, and L2, the device layer, are
implemented. L0, L3, and L4 are later work."

**On-screen text.** `L0 quantum — later`, `L1 atomistic — implemented`,
`L2 device — implemented`, `L3 continuum — later`, `L4 system — later`.

**Source.** `ARCHITECTURE.md` and `docs/glossary.md`.

### Shot 3 — Clean checkout and one command (19 to 29 s)

**Visual.** A terminal in a clean directory. The narrator types the clone and
then the reproduce command. The terminal prints the seven steps.

**Command on screen.**

```sh
git clone <repository-url>
cd nano-cad
just reproduce
```

**Narration.** "A clean checkout reproduces the whole demo with one command.
The first run needs network access. It downloads the build and test tools."

**On-screen text.** `7 steps`, `needs network on the first run`.

**Source.** `docs/reproduce.md` and `README.md`. The repository URL is unknown
in the sources, so the script shows a placeholder.

### Shot 4 — The Rust gates pass (29 to 40 s)

**Visual.** The terminal runs the test gates. The test summary lines scroll.
The finite-difference tests are highlighted.

**Commands on screen.**

```sh
cargo test --workspace
sh scripts/reproduce.sh
```

**Narration.** "Every force term has a finite-difference gradient test. Every
integrator has an energy-conservation test. The tests are a merge gate, not a
phase."

**On-screen text.** `bond 1.28e-17 N`, `angle 8.3e-18 N`, `total 5.7e-17 N`,
`NVE drift 1.278e-5`.

**Source for the numbers.** `TASKS.md` M2-02, M2-03, M2-08, and M2-10.

### Shot 5 — The agent generates the gear set (40 to 52 s)

**Visual.** A split screen. Left: the JSON-RPC request. Right: an animation of
the sun, planets, and ring appearing. The handle `part-1` appears.

**Command on screen.**

```json
{"jsonrpc":"2.0","id":2,"method":"tools/call","params":{"name":"generate_gear","arguments":{"name":"planetary","specs":{"planet_count":3}}}}
```

**Narration.** "The agent calls typed tools. It never edits coordinates. It
asks for a planetary set with three planets. The generator returns the atoms
and the bonds."

**On-screen text.** `sun 24`, `planet 18`, `ring 60`, `planets 3`.

**Source.** `TASKS.md` M4-06, `docs/three-scale.md`, and `mcp/README.md`.

### Shot 6 — The agent assembles the L2 device (52 to 62 s)

**Visual.** The JSON-RPC request, then the device view. Seven bodies and six
joints appear. Revolute joints are marked with axis arrows.

**Command on screen.**

```json
{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"assemble_planetary","arguments":{"specs":{"planet_count":3}}}}
```

**Narration.** "The agent assembles the L2 device. The device has seven bodies
and six joints. The constraints are position-level projections."

**On-screen text.** `7 bodies`, `6 joints`, `position-level projection`.

**Source.** `TASKS.md` M7-04, M8-01, and `mcp/tests/test_agent_demo.py`.

### Shot 7 — The agent measures the gear ratio (62 to 72 s)

**Visual.** The JSON-RPC request, then the running device with a live ratio
readout. The readout settles at 3.500000000.

**Command on screen.**

```json
{"jsonrpc":"2.0","id":4,"method":"tools/call","params":{"name":"measure_gear_ratio","arguments":{"assembly_id":"assembly-1"}}}
```

**Narration.** "The agent drives the sun and measures the sun-to-carrier
ratio. The analytic ratio is three point five. The simulated ratio is three
point five. The relative error is one point one four six times ten to the
minus ten."

**On-screen text.** `analytic 3.500000000`, `measured 3.500000000`,
`relative error 1.146e-10`, `tolerance 1.0e-3`.

**Source for the numbers.** `TASKS.md` M6-06 and `docs/reproduce.md`.

### Shot 8 — URDF and the three-scale scene (72 to 82 s)

**Visual.** The JSON-RPC request, then the URDF tree. Next, the three-scale
scene cross-fades from atoms to bodies to pitch circles.

**Command on screen.**

```json
{"jsonrpc":"2.0","id":5,"method":"tools/call","params":{"name":"export_urdf","arguments":{"assembly_id":"assembly-1"}}}
```

**Narration.** "The agent exports a URDF with seven links and six joints. The
scene export carries three layers: atomistic, device, and coarse. A viewer
animates across them."

**On-screen text.** `URDF: 7 links, 6 joints`, `atomistic`, `device`,
`coarse`.

**Source.** `TASKS.md` M8-01, `docs/three-scale.md`, and
`mcp/tests/test_agent_demo.py`.

### Shot 9 — Benchmark timings (82 to 88 s)

**Visual.** A terminal runs the benchmark. A bar chart draws the three means.

**Commands on screen.**

```sh
cargo bench -p nanocad-engine
./scripts/bench.sh
```

**Narration.** "The engine benchmark runs on a periodic water box with six
hundred atoms. The serial force pass takes six hundred eighty microseconds.
The four-worker pass takes three hundred seventy-one microseconds."

**On-screen text.** `energy 680.97 us`, `serial 679.83 us`,
`4-worker 370.97 us`, `Apple M4, 10 cores`.

**Source for the numbers.**
`benchmarks/results/2026-09-14-engine-timings.txt` and `docs/benchmarks.md`.

**Caveat on screen.** `one host; your numbers will differ`.

### Shot 10 — Caveats and closing card (88 to 90 s)

**Visual.** Closing card. The title. The repository. The word `simulated`.

**Narration.** "All results are simulated. The OpenMM cross-check agrees to
1e-14 relative, but it is a coding check, not a physics validation. This is a
design hypothesis, not a validated device."

**On-screen text.** `simulated`, `Apache-2.0`, `github.com/<owner>/nano-cad`.

**Source.** `PLAN.md`, `README.md`, and `TASKS.md` M2-11. The repository URL
is unknown in the sources, so the script uses a placeholder.

## Full command list

These commands appear on screen. They come from `docs/reproduce.md` and
`mcp/tools.json`.

| Shot | Command or message | Source |
|---|---|---|
| 3 | `git clone <repository-url>` | placeholder; no URL in the sources |
| 3 | `cd nano-cad` | placeholder path |
| 3 | `just reproduce` | `docs/reproduce.md`, `README.md` |
| 4 | `cargo test --workspace` | `docs/reproduce.md`, `justfile` |
| 4 | `sh scripts/reproduce.sh` | `docs/reproduce.md` |
| 4 | `python -m pytest mcp/tests/test_reproduce_script.py` | `docs/reproduce.md` |
| 4 | `python -m pytest mcp/tests` | `mcp/README.md` |
| 5 | JSON-RPC `generate_gear` | `mcp/tools.json`, `mcp/README.md` |
| 5 | `initialize` handshake | `mcp/README.md` |
| 6 | JSON-RPC `assemble_planetary` | `mcp/tools.json` |
| 7 | JSON-RPC `measure_gear_ratio` | `mcp/tools.json` |
| 8 | JSON-RPC `export_urdf` | `mcp/tools.json` |
| 9 | `cargo bench -p nanocad-engine` | `docs/benchmarks.md` |
| 9 | `./scripts/bench.sh` | `docs/benchmarks.md` |

The JSON-RPC messages follow the shape in `mcp/README.md`:

```json
{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2025-06-18","capabilities":{},"clientInfo":{"name":"demo","version":"1.0"}}}
```

The handle names are `part-1`, `document-1`, and `assembly-1`. Source:
`mcp/README.md`.

## Assets needed

### Rendered assets

- A three-scale scene JSON file from `build_scene`. The default design has 7
  bodies, 6 joints, and 3 planets. Source: `TASKS.md` M8-01.
- A rendered orbit of the device layer, 8 seconds, 16:9.
- A rendered cross-fade from the atomistic layer to the coarse layer, 4
  seconds.
- A rendered planetary set with 24 sun teeth, 18 planet teeth, and 60 ring
  teeth. Source: `docs/three-scale.md`.
- A rendered running device with a live ratio readout.
- A URDF tree graphic with 7 links and 6 joints. Source: `TASKS.md` M7-04.
- A bar chart of the three benchmark means. Source:
  `benchmarks/results/2026-09-14-engine-timings.txt`.

### Screen recordings

- A clean terminal session for shot 3.
- A terminal session of `cargo test --workspace` for shot 4.
- A terminal session of `cargo bench -p nanocad-engine` for shot 9.
- A terminal session that sends the JSON-RPC messages, or a screen-capture of
  an MCP client.

### Text and graphics

- The L0 to L4 animated diagram.
- The count-up text: `1.28e-17 N`, `8.3e-18 N`, `5.7e-17 N`, `1.278e-5`.
- The title card and the closing card.
- The word `simulated`, shown in shots 1 and 10.
- The caveat line `one host; your numbers will differ`.

### Audio

- Voice-over recording, 90 seconds.
- Background music, low volume, royalty-free.
- A caption file, `.srt`.

### Files and tools

- A terminal with a monospace font that shows `µs` and `e` notation.
- The repository at a known commit.
- A machine that can run `just reproduce`.

## Production caveats

- The video must not show a physical device. No physical device exists.
- The video must label every result as simulated.
- The video must show the repository URL as a placeholder until a real URL
  exists. No URL is in the sources.
- The benchmark shot must show the host and the caveat. The timings are from
  one Apple M4 host. Source: `docs/benchmarks.md`.
- Shot 4 may state the OpenMM cross-check, but must call it a coding check,
  not a physics validation. Source: `benchmarks/results/openmm-water-crosscheck.txt`.
- The disk estimate, the clone URL, and the repository path in shot 3 are
  placeholders or estimates. The sources do not give them.
