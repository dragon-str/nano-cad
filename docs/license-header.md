# License Header Policy

All source files are Apache-2.0. See `LICENSE` and `NOTICE`.

## Rules

1. Do not copy code from NanoEngineer or any other GPL-licensed project. Read
   the old C engine for algorithms only, then reimplement in clean code.
2. Every new Rust crate carries `license.workspace = true`, which resolves to
   `Apache-2.0`.
3. Do not add a dependency without checking its license. Record the choice in
   `ARCHITECTURE.md`.
4. Keep LGPL dependencies (ASE, xTB) dynamic-link only.
5. Vendored third-party code is not allowed. Depend on a crate instead.

Every crate in this workspace uses `#![forbid(unsafe_code)]` except one audited
FFI crate, if one is ever needed. That exception needs an ADR.
