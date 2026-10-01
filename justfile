set shell := ["bash", "-uc"]

# Make cargo reachable even when the caller's shell did not source rustup.
export PATH := env_var("HOME") + "/.cargo/bin:" + env_var("PATH")

# List available recipes.
default:
    @just --list

# Format all Rust code.
fmt:
    cargo fmt --all

# Check formatting without changing files (gate).
fmt-check:
    cargo fmt --all -- --check

# Lint with clippy, warnings are errors (gate).
lint:
    cargo clippy --workspace --all-targets --all-features -- -D warnings

# Run the Rust test suite (gate).
test:
    cargo test --workspace

# Run the Python test suite (gate). uv supplies pytest, numpy and scipy, so
# the gate does not depend on the system Python having them (M28-37).
py-test:
    uv run --quiet --no-project --python 3.12 --with pytest --with "numpy>=2" --with scipy env PYTHONPATH=python:app python -m pytest -q python/tests app/tests

# Check the site and the generated scene (needs python3).
site-check:
    python3 site/check.py

# Run the viewer JavaScript checks and renderer tests (needs node). Not a gate.
js-test:
    node --check site/render_atoms.js
    node --check site/viewer.js
    node site/render_atoms.test.js
    node --check site/solid_mesh.js
    node site/solid_mesh.test.js

# Run timing benchmarks. Not a gate.
bench:
    cargo bench

# Reproduce the whole demo from a clean checkout. Needs network on first run.
reproduce:
    ./scripts/reproduce.sh

# Run every merge gate. Must return zero.
verify: fmt-check lint test py-test
    @echo "verify: all gates passed"
