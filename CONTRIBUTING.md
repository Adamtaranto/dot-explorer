# Contributing to dot-explorer

Thanks for your interest in dot-explorer! Bug reports, feature ideas,
documentation fixes and code contributions are all welcome.

By contributing you agree that your work is licensed under the project's
[GPL-3.0-or-later](LICENSE) licence.

## Reporting bugs and requesting features

Open an issue using one of the templates on the
[issue tracker](https://github.com/Adamtaranto/dot-explorer/issues):

- **Bug report** — please include the dot-explorer version, your Python and
  Rust versions, the OS, and a minimal reproducer. For browser-app bugs,
  include the browser and anything the JavaScript console printed.
- **Feature request** — describe the analysis you are trying to do, not only
  the API you have in mind; there may already be a way to do it.

If you are unsure whether something is a bug, open an issue and ask — that is
never wasted effort.

## Development setup

Full setup instructions — installing Rust, creating the conda environment,
`maturin develop`, the pre-commit hooks, running tests and benchmarks, and
building the docs — live in the
**[Development Guide](https://adamtaranto.github.io/dot-explorer/development/)**
(source: [`docs/development.md`](docs/development.md)).

The short version, once you have a Rust toolchain and Python ≥ 3.12:

```bash
git clone https://github.com/<your-username>/dot-explorer.git
cd dot-explorer
pip install maturin
maturin develop --extras dev,docs
pre-commit install
pre-commit install --hook-type pre-push
```

Re-run `maturin develop` after **any** change to `src/*.rs`. Pure-Python
changes under `python/dot_explorer/` take effect immediately.

Work on the browser app has its own, separately pinned toolchain — see
[the app README](python/dot_explorer/app/README.md) before touching anything wasm-related.

## Where code goes

| Area | Location | Notes |
|---|---|---|
| Rust core (k-mer index, matching, PAF, gravity ordering) | `src/` | Must stay buildable for `wasm32-unknown-emscripten` with `--no-default-features` |
| Python API (plotting, PAF I/O, annotation, clustering, HTML report) | `python/dot_explorer/` | |
| Type stubs for the Rust extension | `python/dot_explorer/_dot_explorer.pyi` | Update alongside any pyclass/pyfunction signature change |
| Browser app logic | `python/dot_explorer/app/core/` | Pure Python, unit-tested natively via `tests/test_app_*.py` |
| Browser app UI wiring | `python/dot_explorer/app/app.py` | Keep thin — logic belongs in `core/` |
| Tests | `tests/` | |
| Benchmarks | `benches/` (Rust), `python/benchmarks/` (Python) | |

A few constraints worth knowing before you start:

- **Anything new in `src/` that uses `needletail` or `rayon` must be feature
  gated** (`fasta` / `parallel`), and Python callers must degrade gracefully
  without it. The wasm build excludes both.
- **`similarity.py`, `tree.py` and `heatmap.py` are limited to sourmash 4.8-era
  APIs** — Pyodide bundles sourmash 4.8.11, and the app has to run the same code.
- **The gravity contig ordering is implemented twice**, in Rust
  (`SequenceIndex.optimal_contig_order`) and in Python
  (`paf_io.compute_gravity_contigs`). `tests/test_gravity_parity.py` requires
  them to agree; if you change one, change the other.
- **`tests/test_kmer_index_parity.py`** checks the ntHash engine against a naive
  brute-force reference. It should keep passing unchanged — if a change to the
  matching path requires editing it, say so explicitly in the PR.

## Making a change

1. Fork the repository and create a branch off `main`.
2. Make your change, with tests.
3. Run the checks below locally (pre-commit does most of this for you).
4. Add a `CHANGELOG.md` entry under `## [Unreleased]`, in the matching
   `### Added` / `### Changed` / `### Fixed` / `### Performance` / `### Removed`
   section.
5. Open a pull request against `main` describing what changed and why.

### Checks

These are the same gates CI enforces (`.github/workflows/ci.yml`):

```bash
cargo fmt --all -- --check
cargo clippy -- -D warnings
cargo test
ruff format python/ tests/
ruff check --fix python/ tests/
pytest tests/ -v
```

Both `pytest tests/ -x -q` and `cargo test --lib` also run automatically at
push time if you installed the pre-push hooks.

CI additionally builds the package on free-threaded Python 3.14t and asserts
the extension does not re-enable the GIL, builds a native wheel and a Pyodide
wasm wheel, and builds the docs site.

### Tests

- Python tests go in `tests/`, named `test_*.py`; Rust unit tests live
  alongside the code they cover in `src/`.
- Run a single test with
  `pytest tests/test_index.py::test_get_paf_all_returns_paf_lines -v`.
- Tests must not depend on real assemblies. Benchmarks and tests use
  deterministic seeded synthetic DNA; anything needing private local data goes
  behind the `private_data` marker (opt in with `DOT_EXPLORER_PRIVATE_DATA=1`)
  and the data is never committed.
- Optional dependencies (`sourmash`, `scipy`, `shiny`) should be
  `importorskip`-ed so the suite still runs without them.

### Style

Formatting and linting are automated, so mostly just run the tools:

- **Python** — ruff (`pyproject.toml`): single quotes, 88-column lines,
  numpy-convention docstrings. Docstrings are required on library code and
  `python/dot_explorer/app/core/`, not on tests.
- **Rust** — `cargo fmt`, and Clippy warnings are errors.
- Comments are most useful where they explain *why* a non-obvious choice was
  made (there is a lot of that in this codebase — pinned toolchains, memory
  layouts, wasm constraints). Please keep that habit up.

### Commit messages and PRs

Commits follow [Conventional Commits](https://www.conventionalcommits.org/):
`feat(app): ...`, `fix(plot): ...`, `perf(lib): ...`, `docs: ...`,
`build(deps): ...`, with `!` for breaking changes. Pull requests are squash
merged, so the PR title becomes the commit message — please make it a
well-formed conventional-commit line.

### Documentation

Docs are built with [Zensical](https://zensical.org) from `docs/`. Tutorials
are Jupyter notebooks: **edit `docs/tutorials/*.ipynb`, never the generated
`.md`** files, which are build products.

```bash
python scripts/notebooks_to_md.py
zensical serve
```

API pages are generated by mkdocstrings, which imports `dot_explorer` — run
`maturin develop` first, and make sure no stale `.so` from a different Python
version is sitting in `python/dot_explorer/`.

### Performance work

Performance is tracked in CI with [CodSpeed](https://codspeed.io/). If you are
changing the indexing or matching hot paths, please include benchmark numbers:

```bash
cargo bench --bench bench_compare
pytest python/benchmarks --codspeed
python scripts/mem_profile_index.py --mb 10 50 100   # memory, which CodSpeed does not measure
```

Memory matters as much as speed here: the browser app runs in a ~4 GB wasm
heap, and the index's per-base-pair footprint is what sets the app's input
size limit.
