#!/bin/sh
# Regenerate the nano-cad 90-second gearbox video and its caption sidecar.
#
# This script is POSIX sh. It calls /opt/homebrew/bin/python3 for the frame
# rendering and ffmpeg for the encode. It exits nonzero when a required tool is
# missing, or when the atom-geometry check fails, so a broken host or a bad
# atomistic layer fails loudly.
#
# Usage:
#   sh scripts/make_video.sh
#
# Outputs:
#   docs/media/nano-cad-gearbox.mp4   H.264, 1920x1080, 30 fps, ~90 s, audio
#   docs/media/nano-cad-gearbox.srt   sidecar captions
#
# Side effects:
#   site/scene.json, site/scene.data.js, site/scene.bonds.json are regenerated
#   from the real Rust generators before the render.
#
# The video is schematic and simulated. It is not a physical device.

set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"

PYTHON=${PYTHON:-/opt/homebrew/bin/python3}
CARGO=${CARGO:-"$HOME/.cargo/bin/cargo"}
FFMPEG=${FFMPEG:-/opt/homebrew/bin/ffmpeg}
FFPROBE=${FFPROBE:-/opt/homebrew/bin/ffprobe}
SAY=${SAY:-/usr/bin/say}

OUT="$ROOT/docs/media/nano-cad-gearbox.mp4"
SRT="$ROOT/docs/media/nano-cad-gearbox.srt"
SCENE="$ROOT/site/scene.json"
BONDS="$ROOT/site/scene.bonds.json"
BUILD=${NCAD_VIDEO_BUILD:-"$ROOT/target/video-build"}

die() {
    echo "make_video: error: $1" >&2
    exit 1
}

[ -x "$PYTHON" ] || die "python interpreter not found or not executable: $PYTHON (set PYTHON=...)"
[ -x "$CARGO" ] || die "cargo not found or not executable: $CARGO (set CARGO=...)"
[ -x "$FFMPEG" ] || die "ffmpeg not found or not executable: $FFMPEG (set FFMPEG=...)"
[ -x "$FFPROBE" ] || die "ffprobe not found or not executable: $FFPROBE (set FFPROBE=...)"

"$PYTHON" -c "import PIL, numpy" 2>/dev/null \
    || die "$PYTHON lacks Pillow or numpy (use /opt/homebrew/bin/python3)"

# Regenerate the coordinates and the bond topology from the real generators.
# The atomistic layer is only real when it comes from this step.
echo "make_video: generating scene and bond data..."
"$CARGO" run --quiet -p nanocad-jigs --example scene_json -- "$SCENE" \
    || die "the scene_json example failed"

# The accuracy gate. The atom layer must be real diamondoid carbon.
echo "make_video: checking atom geometry..."
"$PYTHON" "$ROOT/scripts/check_atom_geometry.py" --bonds "$BONDS" \
    || die "the atom-geometry check failed; the atomistic layer is not diamondoid carbon"

NO_SAY=""
if [ ! -x "$SAY" ]; then
    echo "make_video: warning: say not found; writing a silent video plus srt" >&2
    NO_SAY="--no-say"
fi

mkdir -p "$ROOT/docs/media"
rm -f "$OUT" "$SRT"

echo "make_video: rendering -> $OUT"
# shellcheck disable=SC2086
"$PYTHON" "$ROOT/scripts/render_video.py" \
    --out "$OUT" \
    --srt "$SRT" \
    --build "$BUILD" \
    --bonds "$BONDS" \
    --ffmpeg "$FFMPEG" \
    --ffprobe "$FFPROBE" \
    --say "$SAY" \
    $NO_SAY

[ -s "$OUT" ] || die "video was not written: $OUT"
[ -s "$SRT" ] || die "srt was not written: $SRT"

DUR=$("$FFPROBE" -v error -show_entries format=duration \
    -of default=noprint_wrappers=1:nokey=1 "$OUT")
echo "make_video: done: $OUT (${DUR}s)"
echo "make_video: done: $SRT"
