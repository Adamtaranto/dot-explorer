# Development Guide

This page explains how to set up a local development environment for dot-explorer so you can edit both the Rust extension and the Python package and run the full test suite.

## Prerequisites

### Install Rust

dot-explorer requires a working Rust toolchain (stable channel).
Install it with [rustup](https://rustup.rs/):

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env
rustup update stable
```

Verify the installation:

```bash
rustc --version
cargo --version
```

### Create a conda environment

A minimal conda environment file is provided at `environment.yml`.
It pins Python to 3.14 (free-threaded build) and installs the key Python dependencies.

```bash
conda env create -f environment.yml
conda activate dot-explorer
```

!!! note
    This environment is for **native** development. Building the browser app's
    wasm wheel needs a separate environment (`environment-wasm.yml`) because
    the wasm build cannot use the `python-freethreading` build pinned here,
    and it needs a pinned Emscripten toolchain. See
    [the app README](https://github.com/Adamtaranto/dot-explorer/blob/main/python/dot_explorer/app/README.md).

!!! note "Released wheels"
    Releases publish binary wheels for CPython **3.12**, **3.13**, **3.14**
    and **3.14t** (free-threaded) on Linux, macOS and Windows, an sdist, and a
    PEP 783 `pyemscripten_*_wasm32` wheel — see
    `.github/workflows/publish.yml`. That covers the whole supported range
    (`requires-python = ">=3.12,<3.15"`).

Alternatively, use any Python ≥ 3.12 virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
```

## Fork and clone the repository

[Fork dot-explorer](https://github.com/Adamtaranto/dot-explorer/fork) on GitHub, then clone your fork:

```bash
git clone https://github.com/<your-username>/dot-explorer.git
cd dot-explorer
git remote add upstream https://github.com/Adamtaranto/dot-explorer.git
```

## Install in development mode

[maturin](https://www.maturin.rs/) is the build backend for the Rust extension.
Install it and then build the package in editable mode:

```bash
pip install maturin
# Build Rust extension + install package in editable mode
maturin develop --extras dev,docs
```

The `--extras dev,docs` flag installs all optional development and documentation
dependencies declared in `pyproject.toml`.

Verify the install:

```python
import dot_explorer
print(dot_explorer.__version__)
```

!!! note
    Always re-run `maturin develop` after modifying any Rust source files
    (`src/*.rs`). Pure Python changes in `python/dot_explorer/` take effect
    immediately without a rebuild.

## Installing pre-commit hooks

dot-explorer uses [pre-commit](https://pre-commit.com/) to enforce code quality
on every commit. Install the hooks once after cloning:

```bash
pip install pre-commit   # already included in the dev extras
pre-commit install                 # install commit-stage hooks
pre-commit install --hook-type pre-push   # install push-stage hooks
```

The hooks include:

| Hook | Stage | What it checks |
|------|-------|---------------|
| `check-ast` / `check-yaml` / `check-toml` | commit | File syntax |
| `trailing-whitespace` / `end-of-file-fixer` | commit | Whitespace |
| `cargo fmt --check` | commit | Rust formatting |
| `cargo check` | commit | Rust compilation |
| `cargo clippy -- -D warnings` | commit | Rust linting |
| `ruff format` / `ruff check` | commit | Python formatting & linting |
| `pydocstyle` (`.pyi` stubs only) | commit | Docstring style |
| `cargo test --lib` | push | Rust unit tests |
| `pytest tests/ -x -q` | push | Python test suite |

Run all hooks manually at any time:

```bash
pre-commit run --all-files
```

## Running Python tests

```bash
pytest tests/ -v
```

To run a specific test file or test function:

```bash
pytest tests/test_index.py -v
pytest tests/test_index.py::test_get_paf_all_returns_paf_lines -v
```

## Python code style — ruff

dot-explorer uses [ruff](https://docs.astral.sh/ruff/) for Python linting and
formatting (configured in `pyproject.toml`).

**Check for issues:**

```bash
ruff check python/ tests/
```

**Auto-fix fixable issues:**

```bash
ruff check --fix python/ tests/
```

**Format code:**

```bash
ruff format python/ tests/
```

## Rust code quality

### Format

```bash
cargo fmt
```

Check formatting without modifying files (also used by the pre-commit hook):

```bash
cargo fmt --all -- --check
```

### Lint (Clippy)

```bash
cargo clippy -- -D warnings
```

### Compile check

Quickly verify the crate compiles without producing a binary:

```bash
cargo check
```

### Rust unit tests

Run the Rust-side library tests:

```bash
cargo test --lib
```

Run all Rust tests (including integration tests, if any):

```bash
cargo test
```

## Benchmarks

dot-explorer ships both Rust and Python benchmarks, run in CI via
[CodSpeed](https://codspeed.io/) (`.github/workflows/codspeed.yml`). Both feed on
deterministic, seeded synthetic DNA so results are reproducible; the real
assemblies used during development are private and git-ignored.

### Rust benchmarks

The Rust benches live in `benches/` and use
[`codspeed-criterion-compat`](https://docs.rs/codspeed-criterion-compat) (a
drop-in criterion API), so the same targets run under plain `cargo bench` and
under `cargo codspeed`:

```bash
cargo bench                       # all benches
cargo bench --bench bench_compare # forward/reverse matching only
```

| Bench | Measures |
|-------|----------|
| `bench_index` | Building the rolling-hash k-mer index across sizes and `k` |
| `bench_compare` | Forward and reverse-strand matching for a homologous pair |
| `bench_kmer` | Index primitives (k-mer set, FM-index, k-mer index build) |

!!! note
    On macOS, linking a standalone `cargo bench`/`cargo test` binary against the
    pyo3 extension needs the undefined-symbol linker flags configured in
    `.cargo/config.toml`; these are applied automatically for Apple targets.

### Python benchmarks

The Python benches live in `python/benchmarks/` and use
[`pytest-codspeed`](https://github.com/CodSpeedHQ/pytest-codspeed) (declared in
the `test` extra). Run them with:

```bash
pip install '.[test]'
pytest python/benchmarks --codspeed
```

They cover FASTA byte parsing, the `CrossIndex` build + `compute_matches` path, contig reordering, dotplot rendering, and HTML report generation.

### K-mer memory profile

CodSpeed measures instructions, not memory. `scripts/mem_profile_index.py`
builds `CrossIndex` pairs over synthetic homologous sequences in clean
subprocesses and reports peak RSS plus the index's own accounting
(`SequenceIndex.approx_bytes()`: retained sequence bytes, the CSR k-mer
tables, and the pairwise coordinate cache):

```bash
python scripts/mem_profile_index.py --mb 10 50 100
```

Measured natively (macOS arm64, k=15, 2 % divergence, per-side sizes; the
audit that sized the web app's limits):

| Per side | Peak RSS | Index total | seq bytes | CSR tables | Records |
|---:|---:|---:|---:|---:|---:|
| 10 Mb | 1.6 GB | 0.31 GB | 19 MB | 303 MB | 2.3 M |
| 50 Mb | 4.9 GB | 1.53 GB | 95 MB | 1.47 GB | 15.1 M |
| 100 Mb | 6.4 GB | 2.97 GB | 191 MB | 2.85 GB | 39.2 M |

Take-aways: match records store **coordinates only** (no sequence copies);
the fixed cost is the CSR tables at ~15 B/bp, the retained raw sequence
copies are minor (~1 B/bp), and the rest of the peak is the materialised
Python `PafRecord` objects — `min_block_len` cuts the record count (39 M →
1.0 M at 25 bp here) without changing the index size. This matches the
in-browser measurement (~32 B/bp of combined input at the 2.9 GB peak for a
90 Mb run), which is why the web app's k-mer gate sits at 80 Mb combined:
the 4 GB wasm heap minus the Pyodide baseline supports ~110 Mb in theory,
and 80 Mb keeps a ~30 % margin.

## Building the documentation locally

The site is built with [Zensical](https://zensical.org) (configured in
`zensical.toml`), the successor to MkDocs + Material by the same authors.
Install the docs dependencies (included in `maturin develop --extras dev,docs`),
then render the tutorial notebooks and serve:

```bash
python scripts/notebooks_to_md.py   # docs/tutorials/*.ipynb -> *.md
zensical serve
```

Open <http://127.0.0.1:8000> in your browser. The site rebuilds automatically
when you save a documentation file.

To build a static site into `site/`:

```bash
zensical build --clean
```

!!! note
    Zensical has no notebook support, so `scripts/notebooks_to_md.py` converts
    `docs/tutorials/*.ipynb` to Markdown first (without executing them). The
    generated `.md` files and `*_files/` image directories are build products
    and are gitignored — edit the notebooks, never the Markdown.

    The API pages are generated by mkdocstrings, which imports `dot_explorer`.
    That includes the compiled `_dot_explorer` extension, so run
    `maturin develop` first, and make sure no stale `.so` built for another
    Python version is left in `python/dot_explorer/` — the build fails with
    `AliasResolutionError` if the extension cannot be imported.
