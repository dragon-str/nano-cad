# The 90-second gearbox video

`docs/media/nano-cad-gearbox.mp4` is the M8-02 artifact. It is a 90-second
video about the nano-cad planetary gearbox. This page tells you how to rebuild
it and what it is honest about.

## Files

| File | What it is |
|---|---|
| `docs/media/nano-cad-gearbox.mp4` | H.264, 1920x1080, 30 fps, ~90 s, with an AAC audio stream and burned-in captions. |
| `docs/media/nano-cad-gearbox.srt` | Sidecar captions. They match the narration and the shot timings. |
| `scripts/make_video.sh` | The POSIX `sh` entry point. It regenerates both files from scratch. |
| `scripts/render_video.py` | The frame renderer. It draws with Pillow and calls ffmpeg. |
| `docs/video-script.md` | The storyboard and the narration. |

## Regenerate

```sh
sh scripts/make_video.sh
```

The script is idempotent. It removes the old outputs first. It renders every
frame again and writes both files.

It needs these tools:

| Tool | Default | Purpose |
|---|---|---|
| Python 3 | `/opt/homebrew/bin/python3` | Frame rendering. Needs Pillow and numpy. |
| ffmpeg | `/opt/homebrew/bin/ffmpeg` | H.264 and AAC encode. |
| ffprobe | `/opt/homebrew/bin/ffprobe` | Duration check after the encode. |
| say | `/usr/bin/say` | Text-to-speech narration. |

Set `PYTHON`, `FFMPEG`, `FFPROBE`, `SAY`, or `NCAD_VIDEO_BUILD` in the
environment to override a path. The script exits nonzero and prints a clear
error when a required tool is missing.

The render takes several minutes. It is not a merge gate, because it is slow
and it needs ffmpeg. It is not part of `just verify`.

A silent variant is possible. `sh scripts/make_video.sh` with `SAY` missing
writes a silent video and still writes the `.srt`.

## What is generated and what is captured

**Every frame is generated.** No camera and no screen recording exist in the
pipeline. The terminal shots are drawn, not filmed. The gear scenes are
schematic 2D drawings:

- Gears are circles with teeth. The teeth counts and the pitch radii follow
  `docs/three-scale.md`: 24 sun teeth, 18 planet teeth, 60 ring teeth, module
  `5e-10 m`, sun pitch radius `6e-9 m`, planet `4.5e-9 m`, ring `1.5e-8 m`,
  carrier `1.05e-8 m`.
- Atoms are dots. The three-scale cross-fade is a drawing, not a render of
  computed positions.
- The device bodies and joints are marked with circles and arrows.

**The narration is synthetic.** It is machine text-to-speech, not a human
voice. The voice is the macOS voice `Samantha` (`en_US`) through `/usr/bin/say`
at 180 words per minute. A shot is time-stretched with `ffmpeg atempo` when its
voice line is longer than the shot window. The narration is therefore honest
about being a machine reading, not an actor.

**No number is invented.** Every figure on screen comes from the files that
`docs/video-script.md` names in its source column: `TASKS.md`, `docs/report.md`,
`docs/benchmarks.md`, `benchmarks/results/openmm-water-crosscheck.txt`, and
`benchmarks/results/2026-09-14-engine-timings.txt`.

## Honesty caveats

- **All results are simulated.** The word `SIMULATED` is on screen at all
  times. No physical device exists.
- **No medical claim is made.** The video calls the gearbox a design
  hypothesis. It does not imply a validated device.
- **The visuals are schematic.** They help a viewer see the structure. They are
  not photorealism and they are not computed atom positions.
- **The benchmark is from one host.** The timings are from a single Apple M4
  with 10 cores and cargo 1.98.1. The video shows the caveat "one host; your
  numbers will differ". `docs/benchmarks.md` says the same.
- **The OpenMM cross-check is a coding check.** The video calls it that. It is
  not a physics validation. One fixed configuration is checked.

## Shot timing

| Shot | Start (s) | End (s) | Length (s) |
|---|---|---|---|
| 1. Cold open: the simulated gearbox | 0 | 8 | 8 |
| 2. What this is: the five model levels | 8 | 19 | 11 |
| 3. Clean checkout and one command | 19 | 29 | 10 |
| 4. The Rust gates pass | 29 | 40 | 11 |
| 5. The agent generates the gear set | 40 | 52 | 12 |
| 6. The agent assembles the L2 device | 52 | 62 | 10 |
| 7. The agent measures the gear ratio | 62 | 72 | 10 |
| 8. URDF and the three-scale scene | 72 | 82 | 10 |
| 9. Benchmark timings | 82 | 88 | 6 |
| 10. Caveats and closing card | 88 | 90 | 2 |
| **Total** | | | **90** |

## Limitations

- The narration for shots 9 and 10 is shorter than the storyboard text, because
  the storyboard gives those shots 6 s and 2 s. The video keeps the storyboard
  wording for the other shots.
- The video has no background music. The storyboard lists it as an asset. The
  generated artifact uses narration only.
- The three-scale cross-fade is a schematic illustration, not a viewer reading
  a real `nanocad.scene` JSON file.
