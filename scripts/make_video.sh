#!/bin/sh
# Regenerate the nano-cad 90-second gearbox video and its caption sidecar.
#
# This script is POSIX sh. It calls /opt/homebrew/bin/python3 for the frame
# rendering and ffmpeg for the encode. It exits nonzero when a required tool is
# missing, so a broken host fails loudly.
#
# Usage:
#   sh scripts/make_video.sh
#
# Outputs:
#   docs/media/nano-cad-gearbox.mp4   H.264, 1920x1080, 30 fps, ~90 s, audio
#   docs/media/nano-cad-gearbox.srt   sidecar captions
#
# The video is schematic and simulated. It is not a physical device.

set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"

PYTHON=${PYTHON:-/opt/homebrew/bin/python3}
FFMPEG=${FFMPEG:-/opt/homebrew/bin/ffmpeg}
FFPROBE=${FFPROBE:-/opt/homebrew/bin/ffprobe}
SAY=${SAY:-/usr/bin/say}

OUT="$ROOT/docs/media/nano-cad-gearbox.mp4"
SRT="$ROOT/docs/media/nano-cad-gearbox.srt"
BUILD=${NCAD_VIDEO_BUILD:-"$ROOT/target/video-build"}

die() {
    echo "make_video: error: $1" >&2
    exit 1
}

[ -x "$PYTHON" ] || die "python interpreter not found or not executable: $PYTHON (set PYTHON=...)"
[ -x "$FFMPEG" ] || die "ffmpeg not found or not executable: $FFMPEG (set FFMPEG=...)"
[ -x "$FFPROBE" ] || die "ffprobe not found or not executable: $FFPROBE (set FFPROBE=...)"

"$PYTHON" -c "import PIL, numpy" 2>/dev/null \
    || die "$PYTHON lacks Pillow or numpy (use /opt/homebrew/bin/python3)"

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
