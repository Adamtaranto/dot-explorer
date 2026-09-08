"""Pure-Python core logic for the dot-explorer browser app.

These modules contain everything the Shiny app does that is independent of
the UI layer: FASTA parsing, alignment-method dispatch, session caching, and
plot-configuration state.  They import :mod:`dot_explorer` lazily (inside
functions) so the app can install the wasm wheel at startup before any Rust
code is touched, and so the modules stay importable for native unit tests.
"""
