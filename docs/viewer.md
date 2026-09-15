# Three-scale viewer

This document tells you how to generate the scene file and how to open the
viewer. The viewer is a static page. It has no external library. The view is a
schematic projection, not a physical render. A local server adds a live
parameter panel and a chat interface; see "Interactive mode" below.

## What the viewer shows

The viewer draws the three layers of a `nanocad.scene` document.

- **Atomistic (L1).** One dot for each atom, in colour by element.
- **Device (L2).** One marker for each rigid body, and a line for each joint.
  A fixed body shows a small white centre dot.
- **Coarse (handoff).** A dashed gear pitch circle for each gear and a bounding
  cylinder for each body.

The scale slider blends between the layers. The three checkboxes hide or show
one layer. Drag to rotate, use the wheel to zoom, and double-click to reset.

The scene is the default planetary design: 7 bodies, 6 joints, 3 planets, and
an analytic sun-to-carrier gear ratio of 3.5.

## Generate the scene file

Run the example from the repository root. The first argument is the output
path. The default path is `site/scene.json`.

```sh
cargo run -p nanocad-jigs --example scene_json -- site/scene.json
```

The command prints a one-line summary. The default output has 8340 atoms in
four axial layers, `4.632e-10 m` thick. The file is a few MB. The example
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
module, the tooth counts, the planet count, the number of atomic layers, and
the layer spacing. Change a field to rebuild the gears. The chat box accepts
short commands:

- `one atomic layer thicker` or `make it 2 layers thinner`
- `6 atoms thick` or `layers to 8`
- `add more teeth to the gears` or `sun teeth to 30`
- `4 planets`
- `move the layers 20 pm apart`
- `reset`

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
scene ok: bodies=7 joints=6 planets=3 ratio=3.5 atoms=8340
viewer ok: index.html, viewer.js, style.css present
docs ok: 17 markdown inputs have HTML pages
check: all checks passed
```

## Honesty caveats

- The label **simulated / schematic** is always visible in the viewer.
- The scene is a snapshot of a simulated design. It is not a built device.
- The atomistic layer is the zero configuration. The viewer does not animate
  the device transform. At the zero configuration the two agree.
- The projection is a flat schematic. It is not a physical or space-filling
  render. Atom sizes are not to scale.
- The viewer draws no bonds, no atom types, and no forces. The gears have a
  stated axial thickness, so the atom dots form a thin solid; the render is
  still schematic.
- The gear pitch circles are design values. They are not measured.
- `site/scene.json`, `site/docs/`, and the viewer are generated or static. They
  are not part of `just verify`. The app parser has its own tests under
  `app/tests/`.
- The chat parser is rule-based, not a language model. It accepts the listed
  commands only. It clamps every value and reports the engine error verbatim.
- No medical claim is made. The output is a design hypothesis.
