# Three-scale viewer

This document tells you how to generate the scene file and how to open the
viewer. The viewer is a static page. It has no network access and no external
library. The view is a schematic projection, not a physical render.

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

The command prints a one-line summary. The default output has 2202 atoms. The
file is about 600 KB. The example creates the parent directory.

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
scene ok: bodies=7 joints=6 planets=3 ratio=3.5 atoms=2202
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
- The viewer draws no bonds, no atom types, and no forces.
- The gear pitch circles are design values. They are not measured.
- `site/scene.json`, `site/docs/`, and the viewer are generated or static. They
  are not part of `just verify`.
- No medical claim is made. The output is a design hypothesis.
