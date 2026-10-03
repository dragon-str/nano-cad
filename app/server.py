#!/usr/bin/env python3
"""The nano-cad local web app.

This serves the static viewer in `site/` and adds a small JSON API. The API
runs the real Rust scene generator, so a parameter change in the browser goes
through the same engine as the tests and the CLI. There is no second geometry
implementation.

Run it from anywhere:

    python3 app/server.py --port 8000

Then open http://127.0.0.1:8000/ . The server binds to the loopback address, so
it is not reachable from another machine. Standard library only.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import threading
import urllib.parse
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chat  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SITE = os.path.join(ROOT, "site")
sys.path.insert(0, os.path.join(ROOT, "python"))


def system_report(edits: dict | None) -> dict:
    """The respirocyte system model (M29-01). Pure Python, about 1 s."""
    from nanocad import respirocyte_system  # noqa: PLC0415

    return respirocyte_system.page_report(edits or {})
# The respirocyte device relaxes in about 150 s in a release build (ADR-0116),
# and the first request after a code change also compiles the example.
BUILD_TIMEOUT_S = 600

# The engine runs as a release build. The default gear set has 142091 atoms,
# and the release binary is far faster than the debug one. Set `--debug-engine`
# to rebuild and run the debug profile instead.
CARGO_PROFILE = "--release"
# The score of one parameter set does not change, so keep the last results. The
# sweep is expensive, and a page reload repeats the same request.
_SCORE_CACHE = {}
_SCORE_LOCK = threading.Lock()
# The rotor scene is deterministic. Build it and pack it once, then serve the
# same bytes. A page reload repeats the request, and the build costs seconds.
_ROTOR_CACHE = None
_ROTOR_CACHE_GZIP = None
_ROTOR_LOCK = threading.Lock()
# A part scene is deterministic too. The respirocyte device costs about 100 s
# to build, so cache each part response by id and serve the same bytes.
_PART_CACHE = {}
_PART_LOCK = threading.Lock()

_CONTENT_TYPES = {
    ".html": "text/html; charset=utf-8",
    ".js": "text/javascript; charset=utf-8",
    ".css": "text/css; charset=utf-8",
    ".json": "application/json; charset=utf-8",
    ".svg": "image/svg+xml",
    ".png": "image/png",
    ".ico": "image/x-icon",
}


class BuildError(RuntimeError):
    """The generator rejected the parameters."""


def _cargo() -> str:
    return shutil.which("cargo") or os.path.expanduser("~/.cargo/bin/cargo")


def _run_example(package: str, example: str, args: list) -> subprocess.CompletedProcess:
    """Run a Rust example in the chosen profile and capture its output."""
    command = [
        _cargo(), "run", "--quiet", CARGO_PROFILE,
        "-p", package, "--example", example, "--",
    ]
    command.extend(args)
    return subprocess.run(
        command, cwd=ROOT, capture_output=True, text=True, timeout=BUILD_TIMEOUT_S,
    )


def _param_args(params: dict) -> list:
    """Turn a parameter map into the `key=value` arguments of an example."""
    args = []
    for key, value in params.items():
        args.append(f"{key}={value!r}" if isinstance(value, float) else f"{key}={value}")
    return args


def build_scene(params: dict) -> dict:
    """Run the Rust generator with `params` and return the scene and bonds."""
    workdir = tempfile.mkdtemp(prefix="ncad-app-")
    try:
        output = os.path.join(workdir, "scene.json")
        args = [output]
        args.extend(_param_args(params))
        result = _run_example("nanocad-jigs", "scene_json", args)
        if result.returncode != 0:
            detail = result.stderr.strip().splitlines()
            raise BuildError(detail[-1] if detail else "the generator failed")
        with open(output, encoding="utf-8") as handle:
            scene = json.load(handle)
        with open(os.path.join(workdir, "scene.bonds.json"), encoding="utf-8") as handle:
            bonds = json.load(handle)
        return {"scene": scene, "bonds": bonds, "log": result.stdout.strip()}
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def score_scene(params: dict) -> dict:
    """Run the Rust metric example with `params` and return the metric list."""
    key = tuple(sorted(params.items()))
    with _SCORE_LOCK:
        cached = _SCORE_CACHE.get(key)
    if cached is not None:
        return cached
    result = _run_example("nanocad-meter", "score_json", _param_args(params))
    if result.returncode != 0:
        detail = result.stderr.strip().splitlines()
        raise BuildError(detail[-1] if detail else "the metric run failed")
    lines = result.stdout.strip().splitlines()
    if not lines:
        raise BuildError("the metric run printed nothing")
    try:
        score = json.loads(lines[-1])
    except json.JSONDecodeError as error:
        raise BuildError(f"the metric output was not JSON: {error}") from None
    with _SCORE_LOCK:
        _SCORE_CACHE[key] = score
    return score


def optimize_scene(args: list) -> dict:
    """Run the CMA-ES parameter search and return its result."""
    result = _run_example("nanocad-opt", "optimize_planetary", args)
    if result.returncode != 0:
        detail = result.stderr.strip().splitlines()
        raise BuildError(detail[-1] if detail else "the search failed")
    lines = result.stdout.strip().splitlines()
    if not lines:
        raise BuildError("the search printed nothing")
    try:
        return json.loads(lines[-1])
    except json.JSONDecodeError as error:
        raise BuildError(f"the search output was not JSON: {error}") from None


def rotor_result() -> dict:
    """Run the sorting-rotor example and return its scene and its facts."""
    workdir = tempfile.mkdtemp(prefix="ncad-rotor-")
    try:
        scene_path = os.path.join(workdir, "scene.json")
        result = _run_example("nanocad-jigs", "rotor_json", [scene_path])
        if result.returncode != 0:
            detail = result.stderr.strip().splitlines()
            raise BuildError(detail[-1] if detail else "the rotor run failed")
        lines = result.stdout.strip().splitlines()
        if not lines:
            raise BuildError("the rotor run printed nothing")
        try:
            facts = json.loads(lines[-1])
        except json.JSONDecodeError as error:
            raise BuildError(f"the rotor output was not JSON: {error}") from None
        with open(scene_path, encoding="utf-8") as handle:
            scene = json.load(handle)
        return {"scene": scene, "facts": facts}
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def rotor_body() -> tuple[bytes, bytes]:
    """Return the rotor response body, as raw bytes and as gzip bytes.

    The payload is deterministic, so build it once and keep it. The raw bytes
    serve a client without gzip support. The gzip bytes cut the transfer for a
    gzip-capable browser.
    """
    global _ROTOR_CACHE, _ROTOR_CACHE_GZIP
    with _ROTOR_LOCK:
        if _ROTOR_CACHE is None:
            body = json.dumps(rotor_result()).encode("utf-8")
            _ROTOR_CACHE = body
            _ROTOR_CACHE_GZIP = gzip.compress(body, compresslevel=1)
    return _ROTOR_CACHE, _ROTOR_CACHE_GZIP


def warm_rotor_cache() -> None:
    """Build the rotor response off the request path, so the first load is quick."""
    try:
        rotor_body()
    except (BuildError, subprocess.TimeoutExpired) as error:
        sys.stderr.write(f"app: rotor warm-up failed: {error}\n")


def part_body(part_id: str) -> tuple[bytes, bytes]:
    """Return a part response body as raw and gzip bytes, cached by part id."""
    with _PART_LOCK:
        cached = _PART_CACHE.get(part_id)
        if cached is None:
            body = json.dumps(part_result(part_id)).encode("utf-8")
            cached = (body, gzip.compress(body, compresslevel=1))
            _PART_CACHE[part_id] = cached
    return cached


def warm_default_scene() -> None:
    """Build the first-load device scene off the request path."""
    try:
        part_body("respirocyte_module")
    except (BuildError, subprocess.TimeoutExpired) as error:
        sys.stderr.write(f"app: default scene warm-up failed: {error}\n")


def scenes_result() -> dict:
    """List the machines and the library parts that the viewer can show."""
    result = _run_example("nanocad-jigs", "scenes_json", [])
    if result.returncode != 0:
        detail = result.stderr.strip().splitlines()
        raise BuildError(detail[-1] if detail else "the scene list failed")
    lines = result.stdout.strip().splitlines()
    if not lines:
        raise BuildError("the scene list printed nothing")
    try:
        return json.loads(lines[-1])
    except json.JSONDecodeError as error:
        raise BuildError(f"the scene list was not JSON: {error}") from None


def part_result(part_id: str) -> dict:
    """Run one library part through the part-scene example and return its scene."""
    result = _run_example("nanocad-jigs", "part_scene_json", [part_id])
    if result.returncode != 0:
        detail = result.stderr.strip().splitlines()
        raise BuildError(detail[-1] if detail else "the part scene failed")
    lines = result.stdout.strip().splitlines()
    if not lines:
        raise BuildError("the part scene printed nothing")
    try:
        return json.loads(lines[-1])
    except json.JSONDecodeError as error:
        raise BuildError(f"the part scene was not JSON: {error}") from None


def selectivity_result(radius_m: float | None = None) -> dict:
    """Run the binding example and return the wall-group selectivity table."""
    args = [] if radius_m is None else [repr(radius_m)]
    result = _run_example("nanocad-meter", "binding_json", args)
    if result.returncode != 0:
        detail = result.stderr.strip().splitlines()
        raise BuildError(detail[-1] if detail else "the selectivity run failed")
    lines = result.stdout.strip().splitlines()
    if not lines:
        raise BuildError("the selectivity run printed nothing")
    try:
        return json.loads(lines[-1])
    except json.JSONDecodeError as error:
        raise BuildError(f"the selectivity output was not JSON: {error}") from None


def _params_from_query(query: str) -> dict:
    params = {}
    for key, values in urllib.parse.parse_qs(query).items():
        if key not in chat.LIMITS:
            continue
        try:
            params[key] = float(values[0])
        except (TypeError, ValueError):
            raise BuildError(f"bad value for {key}") from None
    if not params:
        params = dict(chat.DEFAULT_PARAMS)
    return params


class Handler(BaseHTTPRequestHandler):
    server_version = "nanocad-app/1.0"

    def _send_bytes(
        self, status: int, body: bytes, content_type: str, encoding: str = ""
    ) -> None:
        """Send a prepared body, optionally gzip-encoded."""
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        if encoding:
            self.send_header("Content-Encoding", encoding)
            self.send_header("Vary", "Accept-Encoding")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _send_json(self, status: int, payload: dict) -> None:
        body = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _send_file(self, path: str) -> None:
        if not os.path.isfile(path):
            self.send_error(404, "not found")
            return
        extension = os.path.splitext(path)[1].lower()
        with open(path, "rb") as handle:
            body = handle.read()
        self.send_response(200)
        self.send_header("Content-Type", _CONTENT_TYPES.get(extension, "application/octet-stream"))
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _site_path(self, relative: str) -> str:
        cleaned = os.path.normpath(relative).lstrip(os.sep)
        if cleaned.startswith(".."):
            return ""
        return os.path.join(SITE, cleaned)

    def do_GET(self) -> None:  # noqa: N802 (http.server API)
        parsed = urllib.parse.urlparse(self.path)
        route = parsed.path
        if route in ("/", ""):
            self._send_file(os.path.join(SITE, "index.html"))
            return
        if route in ("/system", "/system.html"):
            self._send_file(os.path.join(SITE, "system.html"))
            return
        if route == "/api/system":
            try:
                self._send_json(200, system_report(None))
            except Exception as error:  # noqa: BLE001
                self._send_json(500, {"error": str(error)})
            return
        if route == "/api/scene":
            try:
                with open(os.path.join(SITE, "scene.json"), encoding="utf-8") as handle:
                    scene = json.load(handle)
                self._send_json(200, {"scene": scene})
            except OSError as error:
                self._send_json(500, {"error": str(error)})
            return
        if route == "/api/meta":
            self._send_json(200, {
                "parameters": [
                    {"key": key, "label": label, "scale": scale, "unit": unit}
                    for key, label, scale, unit in chat.PARAM_META
                ],
                "limits": chat.LIMITS,
                "defaults": chat.DEFAULT_PARAMS,
            })
            return
        if route == "/api/build":
            try:
                params = _params_from_query(parsed.query)
                self._send_json(200, build_scene(params))
            except BuildError as error:
                self._send_json(400, {"error": str(error)})
            except subprocess.TimeoutExpired:
                self._send_json(504, {"error": "the generator timed out"})
            return
        if route == "/api/optimize":
            try:
                query = urllib.parse.parse_qs(parsed.query)
                args = []
                for key, limit in (
                    ("population", 16),
                    ("generations", 12),
                    ("steps", 240),
                    ("seed", 1000000),
                ):
                    values = query.get(key)
                    if values and values[0]:
                        args.append(f"{key}={min(int(values[0]), limit)}")
                self._send_json(200, optimize_scene(args))
            except ValueError:
                self._send_json(400, {"error": "the search budget must be a number"})
            except BuildError as error:
                self._send_json(400, {"error": str(error)})
            except subprocess.TimeoutExpired:
                self._send_json(504, {"error": "the search timed out"})
            return
        if route == "/api/rotor":
            try:
                raw, packed = rotor_body()
            except BuildError as error:
                self._send_json(400, {"error": str(error)})
                return
            except subprocess.TimeoutExpired:
                self._send_json(504, {"error": "the rotor run timed out"})
                return
            if "gzip" in self.headers.get("Accept-Encoding", ""):
                self._send_bytes(200, packed, "application/json; charset=utf-8", "gzip")
            else:
                self._send_bytes(200, raw, "application/json; charset=utf-8")
            return
        if route == "/api/scenes":
            try:
                self._send_json(200, scenes_result())
            except BuildError as error:
                self._send_json(400, {"error": str(error)})
            except subprocess.TimeoutExpired:
                self._send_json(504, {"error": "the scene list timed out"})
            return
        if route == "/api/part":
            try:
                query = urllib.parse.parse_qs(parsed.query)
                values = query.get("id")
                part_id = values[0] if values and values[0] else ""
                if not part_id:
                    raise BuildError("the part id is missing")
                raw, packed = part_body(part_id)
            except BuildError as error:
                self._send_json(400, {"error": str(error)})
                return
            except subprocess.TimeoutExpired:
                self._send_json(504, {"error": "the part scene timed out"})
                return
            if "gzip" in self.headers.get("Accept-Encoding", ""):
                self._send_bytes(200, packed, "application/json; charset=utf-8", "gzip")
            else:
                self._send_bytes(200, raw, "application/json; charset=utf-8")
            return
        if route == "/api/selectivity":
            try:
                query = urllib.parse.parse_qs(parsed.query)
                values = query.get("radius_m")
                radius_m = float(values[0]) if values and values[0] else None
                self._send_json(200, selectivity_result(radius_m))
            except ValueError:
                self._send_json(400, {"error": "the well radius must be a number"})
            except BuildError as error:
                self._send_json(400, {"error": str(error)})
            except subprocess.TimeoutExpired:
                self._send_json(504, {"error": "the selectivity run timed out"})
            return
        if route == "/api/score":
            try:
                params = _params_from_query(parsed.query)
                self._send_json(200, score_scene(params))
            except BuildError as error:
                self._send_json(400, {"error": str(error)})
            except subprocess.TimeoutExpired:
                self._send_json(504, {"error": "the metric run timed out"})
            return
        self._send_file(self._site_path(route.lstrip("/")))

    def do_POST(self) -> None:  # noqa: N802 (http.server API)
        parsed = urllib.parse.urlparse(self.path)
        if parsed.path == "/api/system":
            length = int(self.headers.get("Content-Length", "0"))
            raw = self.rfile.read(length) if length else b"{}"
            try:
                self._send_json(200, system_report(json.loads(raw or b"{}")))
            except json.JSONDecodeError:
                self._send_json(400, {"error": "bad JSON body"})
            except Exception as error:  # noqa: BLE001
                self._send_json(500, {"error": str(error)})
            return
        if parsed.path != "/api/chat":
            self.send_error(404, "not found")
            return
        length = int(self.headers.get("Content-Length", "0"))
        raw = self.rfile.read(length) if length else b"{}"
        try:
            payload = json.loads(raw or b"{}")
        except json.JSONDecodeError:
            self._send_json(400, {"error": "bad JSON body"})
            return
        message = str(payload.get("message", ""))
        current = payload.get("params") or dict(chat.DEFAULT_PARAMS)
        command = chat.parse(message, current)
        new_params = dict(current)
        new_params.update(command.changes)
        response = {
            "reply": command.reply,
            "intent": command.intent,
            "changed": bool(command.changes),
            "params": new_params,
        }
        if command.changes:
            try:
                response.update(build_scene(new_params))
            except BuildError as error:
                response["reply"] = (
                    command.reply + " The engine rejected it: " + str(error)
                )
                response["changed"] = False
                response["params"] = current
            except subprocess.TimeoutExpired:
                response["reply"] = command.reply + " The generator timed out."
                response["changed"] = False
                response["params"] = current
        self._send_json(200, response)

    def log_message(self, fmt: str, *args) -> None:
        sys.stderr.write("app: " + (fmt % args) + "\n")


def main() -> int:
    global CARGO_PROFILE
    parser = argparse.ArgumentParser(description="nano-cad local web app")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument(
        "--debug-engine", action="store_true",
        help="run the debug profile of the Rust engine instead of the release profile",
    )
    arguments = parser.parse_args()
    if arguments.debug_engine:
        CARGO_PROFILE = "--debug"
    server = ThreadingHTTPServer((arguments.host, arguments.port), Handler)
    print(f"nano-cad app on http://{arguments.host}:{arguments.port}/ (Ctrl-C stops)")
    print(f"engine profile: {CARGO_PROFILE.lstrip('-')}")
    threading.Thread(target=warm_rotor_cache, daemon=True).start()
    threading.Thread(target=warm_default_scene, daemon=True).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nstopping")
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
