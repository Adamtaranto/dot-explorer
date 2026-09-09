# dot-explorer browser app

A fully client-side assembly-comparison app (Shiny for Python, deployed as a
static Shinylive/Pyodide site). Uploads never leave the browser.

## Layout

- `app.py` — UI + server wiring (thin; logic lives in `core/`)
- `core/` — pure-Python logic, unit-tested natively (`tests/test_app_*.py`)
- `www/` — custom CSS/JS
- `requirements.txt` — packages micropip-installs at startup (Pyodide builds)

This directory ships inside the wheel as package data, so a `pip install
"dot-explorer[app]"` can run it (see `dot_explorer/launch.py`). It is **not**
an importable subpackage: `app.py` uses flat `core.*` imports that resolve once
its own directory is on `sys.path`, which is what both `shiny run` and
`shinylive export` do. There is deliberately no `wheels/` directory here — a
`.whl` inside the package tree would end up nested inside the built wheel, so
`scripts/build_shinylive.py` stages a copy of this directory outside the
package and drops the wasm wheel in there instead.

## Run natively (development)

The fastest edit/reload loop: Python runs on your machine against a natively
built dot-explorer, so there is no wasm wheel and no export step. Use the
native dev environment (`environment.yml`, see
[docs/development.md](../../../docs/development.md)), not the wasm environment
below.

```bash
pip install ".[app]"         # shiny; plus a native dot-explorer:
maturin develop --release    # from the repo root
dot-explorer-app             # or: shiny run --launch-browser python/dot_explorer/app/app.py
```

`dot-explorer-app` is the console script installed with the package; it puts
this directory on `sys.path` and hands the `App` object to `shiny.run_app`.
Running `shiny run` against `app.py` directly works too — `app.py` uses
path-relative `core.*` imports that resolve once Shiny puts this directory on
`sys.path`.

Everything works except the k-mer method's provenance: it uses your local
dot-explorer instead of the wheel. `minimap2`/`nucmer` still run in the browser
via biowasm, so they behave identically — note this means the aligners need
network access to the biowasm CDN even when the app itself runs locally.
The k-mer upload-size gate is wasm-only; native runs are bounded by machine
RAM instead.

### Run on a remote server / HPC

The same command works headless — bind a port and tunnel to it:

```bash
dot-explorer-app --host 0.0.0.0 --port 8000 --no-browser   # on the server
ssh -L 8000:localhost:8000 user@server            # from your laptop
# then open http://localhost:8000
```

On shared systems, prefer `--host 127.0.0.1` plus the SSH tunnel so the app
is not exposed to other users on the network.

## Build the static Shinylive site

The exported site runs entirely in the browser, so dot-explorer must be
cross-compiled to a Pyodide-compatible wasm wheel. Four things have to agree
exactly, and all of them follow from the Pyodide release **shinylive bundles**
— not from what is newest:

| Component | Pinned to | Why |
|---|---|---|
| Pyodide | 0.27.7 (whatever shinylive bundles) | supplies the runtime |
| CPython | 3.12 | Pyodide 0.27 ships CPython 3.12.7 |
| Emscripten | 3.1.58 | the wheel's platform tag derives from it |
| Rust | `nightly-2025-05-01` | `-Z link-native-libraries=no`, and the pre-wasm-eh ABI |

The wheel filename encodes the middle two
(`…-cp312-cp312-emscripten_3_1_58_wasm32.whl`). If either is wrong, micropip
rejects the wheel at app startup and the app never finishes loading.

> **Do not bump Emscripten on its own.** 4.0.9 is not a drop-in upgrade: it is
> the Pyodide 0.28/0.29 ABI (CPython 3.13, WebAssembly exception handling, a
> different platform tag, and a Rust `std` that `rustup target add` does not
> provide). Its wheels cannot be installed by the runtime shinylive ships, so
> that move has to wait for shinylive itself. The PyPI wasm wheel built by
> `publish.yml` is a separate artifact and does not change this one.

To check what your shinylive actually bundles (this recipe was verified
against shinylive 0.10.14 → Pyodide 0.27.7, CPython 3.12.7,
`emscripten_3_1_58`):

```bash
python -c "import json;print(json.load(open('site/app/shinylive/pyodide/pyodide-lock.json'))['info'])"
```

### 1. Prerequisites

**Rust nightly, pinned.** `.cargo/config.toml` passes
`-Z link-native-libraries=no` to rustc for the wasm target, which is a
nightly-only flag, so a stable toolchain fails with
`error: the option 'Z' is only accepted on the nightly compiler`. The pin
matters in both directions: newer nightlies default to wasm-exceptions
(incompatible with Pyodide 0.27's non-EH ABI) and older ones predate stdlib
APIs used by our dependencies (`bio` 4.x needs `is_multiple_of`, stabilised in
Rust 1.87).

```bash
rustup toolchain install nightly-2025-05-01 --target wasm32-unknown-emscripten
```

This does **not** change your default toolchain — native builds and
`cargo test` keep using stable. The build command below opts in per-invocation
with `RUSTUP_TOOLCHAIN`. (The repo deliberately has no `rust-toolchain.toml`
for exactly this reason.)

**Python + Emscripten + shinylive**, via the dedicated conda environment
(`environment-wasm.yml` at the repo root). It is separate from
`environment.yml` because that env pins `python-freethreading`, which the wasm
build cannot use. `emscripten` comes from conda-forge at exactly 3.1.58 — no
`emsdk` clone or `emsdk_env.sh` sourcing needed:

```bash
conda env create -f environment-wasm.yml
conda activate dot-explorer-wasm
emcc --version   # must report 3.1.58
```

Check nothing else shadows it: `which emcc` must point inside
`dot-explorer-wasm`. Another activated env with its own emscripten is the most
common cause of a wheel that builds but will not install.

### 2. Build the wasm wheel

```bash
rm -rf target/wheels     # never let a stale wheel reach the staging copy
RUSTUP_TOOLCHAIN=nightly-2025-05-01 maturin build --release \
  --target wasm32-unknown-emscripten \
  --no-default-features \
  -i python3.12
ls target/wheels/*cp312-cp312-emscripten_3_1_58_wasm32.whl   # tag check
```

- `--no-default-features` drops `needletail` (native-only FASTA I/O) and
  `rayon` (the target is single-threaded).
- A cold build (empty `target/wasm32-unknown-emscripten/`, warm cargo
  registry) takes well under a minute; the first ever build also downloads
  crates.
- `-i python3.12` does **not** need a CPython 3.12 on your machine: for
  Emscripten targets maturin uses a built-in interpreter config, so this works
  from the Python 3.14 environment above.
- maturin ≥ 1.14 prints a warning about falling back to the "legacy"
  `emscripten_3_1_58_wasm32` platform tag and suggests
  `MATURIN_PYEMSCRIPTEN_PLATFORM_VERSION` / `MATURIN_PYODIDE_ABI_VERSION`
  (PEP 783). **Ignore it and do not set those variables** — Pyodide 0.27 wants
  the legacy tag; a PEP 783 tag will not install.
- There are **two different wasm wheels** in this repo, and they are not
  interchangeable:

  | | this one (the app) | the PyPI one |
  | --- | --- | --- |
  | tag | `cp312-cp312-emscripten_3_1_58_wasm32` | `cp314-cp314-pyemscripten_2026_0_wasm32` |
  | Pyodide | 0.27.7, what shinylive bundles | 314.x (CPython 3.14) |
  | built by | `ci.yml` / `docs.yml`, plain `maturin build` | `publish.yml`, via `pyodide build` |
  | on PyPI | rejected — the tag is not allowed | yes |

  Never let a `pyemscripten_*` wheel reach the staging copy: the app will fail to
  install it at startup.

Optional, mirrors CI (needs node — `environment-wasm.yml` installs it):

```bash
npm install pyodide@0.27.7
node scripts/wasm_smoke_test.mjs target/wheels/*cp312*emscripten*.whl
```

This loads the wheel in a real Pyodide 0.27.7 runtime and runs a comparison —
it catches a mis-tagged or mis-linked wheel in ~30s, before a full export.

### 3. Run the app from the wheel

```bash
python scripts/build_shinylive.py --out site/app
python scripts/add_loading_splash.py site/app/index.html   # optional, matches deploy
python -m http.server --directory site 8741
# open http://127.0.0.1:8741/app/
```

`build_shinylive.py` copies this directory to `build/shinylive/`, drops the one
matching wheel from `target/wheels/` into `build/shinylive/wheels/`, and
exports that staging copy. It refuses to run if `target/wheels/` holds no
matching wheel or more than one: `shinylive export` bundles **everything**
under the directory it is given, so a leftover wheel from an earlier build
would ship alongside the new one. `pick_wasm_wheel()` in `core/wheels.py`
selects by platform tag rather than sort order for exactly this reason, but two
wheels with the same tag and different versions is still ambiguous.

The export bundles the staging copy — including the wheel in its `wheels/`
directory — into the static site; `ensure_dot_explorer()` in `app.py` installs
it from the Pyodide virtual filesystem at startup. The site must be served over
http (the shinylive service worker does not run from `file://`).

A successful start looks like: loading splash → "Loading Python" → the
sidebar renders. If it hangs on the splash, open the browser console: a
micropip `ValueError` naming the platform tag means the wheel does not match
the runtime (see the version table above).

**Testing changes: clear the service worker first.** A previously loaded
shinylive site caches aggressively, so a rebuild can appear to do nothing.
In DevTools → Application, unregister service workers and clear storage for
`127.0.0.1:8741`, then hard-reload.

### Alternative: use the CI-built wheel

To run the app without a local Rust/Emscripten setup, download the
`wasm-wheel` artifact from a recent run of the "Wasm (Pyodide) Wheel Build"
job in [ci.yml](../../../.github/workflows/ci.yml), drop it into `target/wheels/`, and
start at step 3. The same wheel is built by `docs.yml` for the deployed site.

### Troubleshooting

| Symptom | Cause |
|---|---|
| `the option 'Z' is only accepted on the nightly compiler` | missing `RUSTUP_TOOLCHAIN=nightly-2025-05-01` |
| wheel tag is not `…emscripten_3_1_58_wasm32` | wrong `emcc` on PATH (check `emcc --version`, and that the conda env is activated) |
| app hangs on the loading splash; micropip error in the console | wheel tag mismatch, or a stale wheel left in `target/wheels/` |
| edits do not appear after re-export | stale service worker / browser cache (see above) |
| `maturin` picks the wrong interpreter | an unrelated env is activated; `CONDA_PREFIX` should point at `dot-explorer-wasm` |
| `use of unstable library feature 'unsigned_is_multiple_of'` while compiling `bio` | toolchain older than rustc 1.87 |
| `emcc: error: invalid export name: _ZN…` | Emscripten 4.x on PATH — that ABI needs Pyodide's patched emsdk and is not what this build targets |

## Usage notes

- **Input modes**: upload two assemblies (FASTA/FASTA.gz) *or* a
  precomputed PAF alignment. In PAF mode the aligner options are hidden; an
  optional query-assembly upload enables the reordered-FASTA download.
- **Self-alignment**: tick "Align assembly to itself" to compare one
  assembly against itself without uploading it twice.
- **Identity colouring**: for methods that report identity (minimap2,
  nucmer, PAF import) alignments can be coloured by % identity
  instead of forward/reverse strand.
- **Instant display options**: line width and min match length are applied
  client-side inside the interactive report (CSS + postMessage) — changing
  them never re-renders the plot. Contig orders are cached per ordering
  mode, so re-selecting "maximise colinearity" after trying another mode is
  instant.
- **Trees & clustering** (self-alignments): upload a newick / IQ-TREE
  `.treefile` to fix the contig order and draw the tree beside the matrix,
  or tick "Compute clustering tree" to sketch every contig with sourmash
  and cluster hierarchically. Adds a cluster-assignment table (rows
  highlight their clusters in the plot), bold cluster outlines, a
  similarity-heatmap tab, and CSV downloads. sourmash + scipy install on
  first use (micropip under Pyodide, `pip install "dot-explorer[cluster]"`
  natively); see `docs/clustering.md` for choosing a metric.
- **GFF annotations**: upload GFF3 files (`.gff`/`.gff3`/`.gz`) for the
  query and/or target assembly. Detected feature types get per-type
  visibility toggles and colour pickers; features shade self-vs-self
  panels behind the alignments, and the focused (drill-down) pair view
  adds side annotation tracks with strand arrows and multi-part CDS
  connectors. Clicking a feature in the interactive report shows its
  name, type, coordinates, strand and parent.

## Alignment methods

| Method | Runs | Status |
|---|---|---|
| k-mer matching | dot-explorer wasm wheel (Pyodide) | ✅ |
| PAF import | pure Python | ✅ |
| minimap2 / nucmer | biowasm (Aioli WebWorker; fetched from the biowasm CDN at runtime) | ✅ |
| BLAST | — | not possible: no production WASM build of BLAST+ exists |
