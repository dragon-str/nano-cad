# nanocad Python package

This directory holds the `nanocad` Python package and its tests. The compiled
core lives in the Rust crate at `crates/python`. The build backend is
[maturin](https://www.maturin.rs/).

## Status

The Rust core is built as the extension module `nanocad._core`. It exposes the
units, model, format, part, parameter, engine, jig, and device tools from
`ARCHITECTURE.md`. The `python` Cargo feature turns on the optional `pyo3`
dependency. Without that feature the crate stays a plain library, so
`cargo test` needs no Python interpreter. This follows ADR-0012.

## Requirements

- Python 3.10 or newer.
- `pytest` to run the tests.
- `maturin` to build the extension.
- A Rust toolchain.

## Install

Create and use a virtual environment. Do not change the system Python.

```sh
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install pytest maturin
```

Build and install the package in editable mode. This step selects the `python`
feature through `[tool.maturin]` and installs `nanocad._core`.

```sh
python3 -m pip install -e python
```

## Test

Run the tests from the repository root. The `pythonpath` setting in
`pyproject.toml` makes `nanocad` importable without an install. The extension
must be built first.

```sh
python3 -m pytest python/tests
```

## Build the extension without maturin

For a quick local check, build the cdylib and copy it into the package. On
macOS the extension module needs the dynamic-lookup linker flags.

```sh
RUSTFLAGS="-C link-arg=-undefined -C link-arg=dynamic_lookup" \
  cargo build -p nanocad-python --features python
cp target/debug/libnanocad_python.dylib python/nanocad/_core.abi3.so
```

## Lint and format

Line length is 88 characters for both tools.

```sh
ruff check python
black --check python
```
