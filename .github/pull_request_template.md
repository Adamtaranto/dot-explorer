## What and why

<!-- What does this change, and what problem does it solve?  Link any issue:
     "Fixes #123" / "Part of #123". -->

## Type of change

<!-- Tick all that apply.  The PR title is squash-merged as the commit message,
     so please make it a Conventional Commit line, e.g.
     `feat(app): colour alignments by identity` or `fix(plot): ...`.
     Use `!` (`refactor!:`) for a breaking change. -->

- [ ] Bug fix
- [ ] New feature
- [ ] Performance
- [ ] Documentation
- [ ] Build / CI / dependencies
- [ ] Breaking change (API, CLI, or output format)

## How was this tested?

<!-- New/updated tests, and anything you checked by hand (e.g. rendered a plot,
     ran the app, opened an HTML report). -->

## Checklist

- [ ] `cargo fmt --all -- --check`, `cargo clippy -- -D warnings` and `cargo test` pass
- [ ] `ruff format`/`ruff check` clean and `pytest tests/ -v` passes
- [ ] Tests added or updated for the change
- [ ] `CHANGELOG.md` updated under `## [Unreleased]`
- [ ] Docs updated (`docs/`, module docstrings, or `app/README.md`) if behaviour changed

## Project-specific checks

<!-- Delete the sections that do not apply to this PR. -->

**Touched `src/*.rs`**

- [ ] Rebuilt with `maturin develop` and re-ran the Python suite
- [ ] New `needletail` / `rayon` use is `#[cfg(feature = ...)]`-gated, and Python
      callers degrade gracefully without it (the wasm build uses
      `--no-default-features`)
- [ ] `python/dot_explorer/_dot_explorer.pyi` stubs updated to match any changed
      signature
- [ ] Nothing added that holds GIL-bound cross-call state (statics,
      `GILOnceCell`, `thread_local`, `unsendable` pyclasses) — the module is
      declared `gil_used = false`

**Touched the matching or contig-ordering path**

- [ ] `tests/test_kmer_index_parity.py` still passes unchanged
- [ ] Gravity ordering changed in *both* engines (Rust
      `SequenceIndex.optimal_contig_order` and Python
      `paf_io.compute_gravity_contigs`), with
      `tests/test_gravity_parity.py` passing

**Performance-sensitive change**

- [ ] Benchmark numbers included below (`cargo bench`, `pytest python/benchmarks
      --codspeed`, and/or `python scripts/mem_profile_index.py`) — memory counts
      as much as speed, since the browser app runs in a ~4 GB wasm heap

**Touched `app/`**

- [ ] Logic went in `app/core/` (with tests in `tests/test_app_*.py`), not
      `app/app.py`
- [ ] Verified natively (`shiny run app/app.py`)
- [ ] Verified in a Shinylive export, if the change could behave differently
      under Pyodide (no threads, no native FASTA reader, sourmash 4.8 APIs only)
- [ ] Pinned wasm toolchain versions unchanged, or the change is explained in
      the PR body (see `app/README.md`)

**Touched clustering / similarity**

- [ ] Restricted to sourmash 4.8-era APIs (Pyodide bundles 4.8.11)
- [ ] Optional imports still `importorskip`-ed so the suite runs without the
      `cluster` extra

**Touched `docs/tutorials/`**

- [ ] Edited the `.ipynb` notebooks, not the generated `.md` build products
