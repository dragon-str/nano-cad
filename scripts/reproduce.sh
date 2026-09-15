#!/bin/sh
# Reproduce the nano-cad demo from a clean checkout.
#
# This script needs no prebuilt artifact. It builds the Python extension from
# source, runs the Rust test suite, then runs the MCP gearbox demo. It exits
# nonzero when any step fails.
#
# The first run needs network access to download maturin and pytest. Override
# the virtual environment location with NCAD_VENV=/path/to/venv.
#
# Usage: sh scripts/reproduce.sh

set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(CDPATH= cd -- "$script_dir/.." && pwd)

venv_dir=${NCAD_VENV:-"$repo_root/.venv"}

step() {
    printf '\n==> %s\n' "$1"
}

fail() {
    printf 'ERROR: %s\n' "$1" >&2
    exit 1
}

# Step 1: the toolchains. Stop loudly when one is absent.
step "Step 1/7: check for the Rust and Python toolchains"
if ! command -v cargo >/dev/null 2>&1; then
    printf 'ERROR: cargo is not on PATH.\n' >&2
    printf 'Install the Rust toolchain from https://rustup.rs and retry.\n' >&2
    exit 127
fi
if ! command -v rustc >/dev/null 2>&1; then
    printf 'ERROR: rustc is not on PATH.\n' >&2
    printf 'Install the Rust toolchain from https://rustup.rs and retry.\n' >&2
    exit 127
fi
if ! command -v python3 >/dev/null 2>&1; then
    printf 'ERROR: python3 is not on PATH.\n' >&2
    printf 'Install Python 3.10 or newer and retry.\n' >&2
    exit 127
fi
printf 'cargo:   %s\n' "$(cargo --version)"
printf 'rustc:   %s\n' "$(rustc --version)"
printf 'python3: %s\n' "$(python3 --version)"

# Step 2: the virtual environment. Create it when it is absent.
step "Step 2/7: create or reuse the Python virtual environment"
if [ -x "$venv_dir/bin/python" ]; then
    printf 'Reusing %s\n' "$venv_dir"
else
    printf 'Creating %s\n' "$venv_dir"
    python3 -m venv "$venv_dir" || fail "python3 -m venv failed"
fi
venv_python="$venv_dir/bin/python"
"$venv_python" --version || fail "the virtual environment has no working Python"

# Step 3: maturin and pytest. The download needs network access.
step "Step 3/7: install maturin and pytest (network access needed on the first run)"
if [ -x "$venv_dir/bin/maturin" ] && [ -x "$venv_dir/bin/pytest" ]; then
    printf 'maturin and pytest are already installed; skipping the download\n'
else
    "$venv_python" -m pip install --quiet --upgrade pip \
        || fail "pip install --upgrade pip failed (check network access)"
    "$venv_python" -m pip install --quiet "maturin>=1.7,<2.0" pytest \
        || fail "pip install maturin pytest failed (check network access)"
fi
"$venv_python" -m maturin --version || fail "maturin is not runnable in the venv"

# Step 4: build the extension. Run maturin from python/, because the pyproject
# there sets module-name and python-source. A build from the repo root with
# --manifest-path misnames the module.
step "Step 4/7: build the nanocad extension with maturin from python/"
export VIRTUAL_ENV="$venv_dir"
( cd "$repo_root/python" && "$venv_python" -m maturin develop --release ) \
    || fail "maturin develop failed (it must run from the python/ directory)"

# Step 5: the Rust test suite.
step "Step 5/7: run the Rust test suite"
( cd "$repo_root" && cargo test --workspace ) || fail "cargo test --workspace failed"

# Step 6: the MCP gearbox demo. Keep the log so step 7 can read the result.
step "Step 6/7: run the MCP gearbox demo test"
demo_log=$(mktemp "${TMPDIR:-/tmp}/ncad-reproduce.XXXXXX") || fail "mktemp failed"
trap 'rm -f "$demo_log"' EXIT HUP INT TERM
if ! ( cd "$repo_root" && "$venv_python" -m pytest -s -q mcp/tests/test_agent_demo.py ) \
    > "$demo_log" 2>&1; then
    cat "$demo_log" >&2
    fail "the MCP gearbox demo test failed"
fi
cat "$demo_log"

# Step 7: the key measured result.
step "Step 7/7: key measured result"
result_line=$(grep 'gearbox demo:' "$demo_log" | tail -n 1)
if [ -z "$result_line" ]; then
    fail "the demo test printed no result line"
fi
printf 'REPRODUCE OK: %s\n' "$result_line"
printf 'The agent built the planetary gearbox through the MCP tool surface only.\n'
