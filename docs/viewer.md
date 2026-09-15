# Three-scale viewer

This document tells you how to generate the scene file and how to open the
viewer. The viewer is a static page. It has no external library. The view is a
schematic projection, not a measurement. When WebGL2 is available the atom
layer is a lit sphere render with screen-space ambient occlusion; otherwise
the viewer falls back to flat dots. A local server adds a live parameter panel
and a chat interface; see "Interactive mode" below.

## What the viewer shows

The viewer draws the three layers of a `nanocad.scene` document.

- **Atomistic (L1).** One dot for each atom, in colour by element.
- **Device (L2).** One marker for each rigid body, and a line for each joint.
  A fixed body shows a small white centre dot.
- **Coarse (handoff).** A dashed gear pitch circle for each gear and a bounding
  cylinder for each body.

Only the **Atoms** layer is on by default, so the view shows the atom layer or
the schematic and nothing else. The **Device markers** and **Coarse** checkboxes
add the rigid-body markers and the handoff circles when you want them. The
**Hide** button at the top of the panel collapses the panel to a small tab, so
the render is not covered at high zoom.

Drag to rotate. Hold Shift and drag, or drag with the middle button, to pan.
Use the wheel to zoom. Double-click to reset. The arrow keys pan, `+` and `-`
zoom, `R` resets, `F` toggles fullscreen, `T` toggles a turntable, `P` toggles
motion, and Escape clears a measurement.

The atomistic layer has a level of detail. When the projected carbon-carbon
spacing is six pixels or more, the viewer draws the atoms with WebGL (or flat
depth-shaded dots when WebGL2 is absent). Below six pixels it draws the exact
involute schematic instead. Both use one profile function and one world fit, so
the switch does not move the geometry.

## Controls

The panel groups the controls.

- **View.** Reset, Iso, Top and Front set fixed views. Play animates the
  planetary kinematics: the sun spins, each planet spins and revolves on the
  carrier, the ring stays fixed. The Motion speed slider scales it. Turntable
  orbits the camera. Fullscreen expands the view. Save PNG writes a composite
  image.
- **Display.** Ambient occlusion turns the SSAO pass on and off. Atom size
  scales the sphere radius from space-filling (1.0) to ball (0.4). Section z
  clips the atoms above a plane, so you can see inside the solid gear.
- **Elements.** One checkbox for each element in the scene. Clear the
  Hydrogen box to see the carbon lattice alone.
- **Readout.** Move the pointer over an atom to read its element, index, body
  and position. Click one atom, then a second, to measure the distance between
  them.

The scale bar shows the length of a segment at the middle of the view. The
triad shows the world x, y and z axes. Each atom sphere has a dark silhouette
outline, so single atoms stay visible when they overlap.

The renderer is a display effect only. The occlusion, the lighting and the
sphere radius are not measurements and carry no physical result.

The scene is the default planetary design: 7 bodies, 6 joints, 3 planets, and
an analytic sun-to-carrier gear ratio of 3.5.

## Generate the scene file

Run the example from the repository root. The first argument is the output
path. The default path is `site/scene.json`.

```sh
cargo run -p nanocad-jigs --example scene_json -- site/scene.json
```

The command prints a one-line summary. The default output has 70070 atoms in
four axial layers, `2.67525e-10 m` thick. The file is a few MB. The example
creates the parent directory.

The example also accepts `key=value` arguments that override the generator
parameters. The keys are the names in `PLANETARY_PARAMETERS`.

```sh
cargo run -p nanocad-jigs --example scene_json -- site/scene.json layers=6
```

This writes a scene with six atomic layers. The interactive app uses the same
arguments.

The example also writes a companion file next to the scene, for example
`site/scene.data.js`. The companion assigns the same scene to
`window.NANOCAD_SCENE`. The viewer uses the companion to load the scene from
`file://`, because a browser blocks `fetch` of a local file.

## Open the viewer

Open `site/index.html` in a web browser by double-clicking it, or with:

```sh
open site/index.html
```

The viewer loads the scene automatically. It tries the companion
`site/scene.data.js` first. That works from `file://`. If the companion is
absent, the viewer tries `fetch("scene.json")`. A browser can block that
request from a `file://` page. When both paths fail, the panel shows a file
picker. Pick `site/scene.json` in the picker. You can also drop the file onto
the page.

To serve the folder with a local web server instead:

```sh
python3 -m http.server 8000 --directory site
```

Then open `http://localhost:8000/`. The server is local. It needs no network.

## Interactive mode: parameters and chat

`app/server.py` serves the same page and adds a live parameter panel and a chat
interface. It runs the Rust generator for each change, so the panel and the
chat always show engine output, never an invented number.

```sh
python3 app/server.py --port 8000
```

Then open `http://localhost:8000/`. The page shows numeric fields for the
module, the tooth counts, the planet count, and the number of atomic layers.
Change a field to rebuild the gears. The chat box accepts short commands:

- `one atomic layer thicker` or `make it 2 layers thinner`
- `6 atoms thick` or `layers to 8`
- `add more teeth to the gears` or `sun teeth to 30`
- `4 planets`
- `reset`

One atomic layer is one diamond (001) plane. The crystal fixes the spacing at
89.175 pm, so the layer count is the free control and the spacing is not.

The parser is rule-based and honest: it recognizes these commands and clamps
each value to a stated limit. It returns the new state and says why. A command
it does not understand changes nothing and returns help. `app/chat.py` holds
the parser; `app/tests/test_chat_parser.py` tests it.

The server binds to `127.0.0.1` only. It serves these routes:

| Route | Purpose |
|---|---|
| `GET /api/meta` | The parameters, their display units, and their limits. |
| `GET /api/build?<params>` | Build a scene with the given parameters. |
| `POST /api/chat` | Parse a message, apply the change, and build. |
| `GET /api/scene` | The last generated `site/scene.json`. |

When the engine rejects a value, for example a planet with too few teeth for
the module, the reply states the engine error and keeps the old parameters.


## Generate the documentation site

The docs generator reads every top-level markdown file and every `docs/*.md`
file. It writes one HTML page for each file under `site/docs/`, plus an index.
It uses the Python standard library only.

```sh
python3 scripts/build_docs_site.py
```

Open `site/docs/index.html` to read the result. The generator fails with a
non-zero status when a required input file is missing.

## Check the result

The check script verifies the scene counts, the viewer files, and the docs
outputs.

```sh
python3 site/check.py
```

Expected output:

```
scene ok: bodies=7 joints=6 planets=3 ratio=3.5 atoms=70070
viewer ok: index.html, viewer.js, style.css, render_atoms.js present
docs ok: 18 markdown inputs have HTML pages
check: all checks passed
```

## Honesty caveats

- The label **simulated / schematic** is always visible in the viewer.
- The scene is a snapshot of a simulated design. It is not a built device.
- The atomistic layer animates the rigid-body kinematics only. It shows the
  zero-stress geometry turning. It does not simulate forces, heat or strain.
- The projection is a schematic. The WebGL render is a display effect. The
  ambient occlusion, the lighting and the sphere radii are approximate. An atom
  sphere is drawn at the covalent radius, not at the true van der Waals radius.
  The scale bar is a projection aid, not a calibrated measurement.
- The viewer draws no bonds and no forces. The gears are solid
  hydrogen-capped diamond, so the atom spheres form a solid with the crystal
  thickness; the render is an illustration of that solid. The ball-and-stick
  `setBonds` path exists in `site/render_atoms.js` but the viewer does not load
  `site/scene.bonds.json`, so no bonds are drawn.
- The gear pitch circles are design values. They are not measured.
- `site/scene.json`, `site/docs/`, and the viewer are generated or static. They
  are not part of `just verify`. The app parser has its own tests under
  `app/tests/`.
- The chat parser is rule-based, not a language model. It accepts the listed
  commands only. It clamps every value and reports the engine error verbatim.
- No medical claim is made. The output is a design hypothesis.
