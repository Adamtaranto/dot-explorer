# Installation

## Requirements

- Python 3.12 – 3.14 (including free-threaded 3.14t)
- Nothing else for a normal install — binary wheels are published for Linux
  (x86_64/aarch64), macOS (arm64/x86_64) and Windows, so pip never has to
  compile the Rust extension.
- A Rust toolchain and [maturin](https://www.maturin.rs/) ≥ 1.0 only if you
  build from a checkout.

## Install from PyPI

```bash
pip install dot-explorer
```

That gives you the library and the interactive HTML reports. To also run the
browser app locally, add the `app` extra and use the console script:

```bash
pip install "dot-explorer[app]"
dot-explorer-app
```

`dot-explorer-app` starts the Shiny app on <http://127.0.0.1:8000> and opens a
browser. See [Web App](webapp.md) for the options.

## Optional extras

| Extra | Installs | Gives you |
|-------|----------|-----------|
| `app` | `shiny` | The local browser app (`dot-explorer-app`) |
| `cluster` | `sourmash ≥ 4.8`, `scipy ≥ 1.10` | Similarity metrics, clustering and trees |
| `docs` | `zensical`, `mkdocstrings-python`, `nbconvert` | Building these docs |

Combine them as usual:

```bash
pip install "dot-explorer[app,cluster]"
```

Both are optional rather than core dependencies for concrete reasons: `shiny`
pulls in `orjson`, which has no free-threaded wheel, so requiring it would make
the package uninstallable on 3.14t; and `scipy` + `sourmash` add roughly 110 MB,
more than doubling the install for users who only want dot plots.

### Core dependencies

Installed automatically:

| Package | Purpose |
|---------|---------|
| `matplotlib ≥ 3.5` | Dotplot visualisation |
| `numpy ≥ 1.21` | Array operations used by matplotlib |
| `pyfaidx` | Lazy, indexed FASTA access (pure Python, ~0.2 MB) |

## Install from a checkout

Use this if you want to modify the Rust core. It needs a Rust toolchain:

```bash
git clone https://github.com/Adamtaranto/dot-explorer.git
cd dot-explorer
pip install maturin
maturin develop --release
```

The `--release` flag enables full Rust compiler optimisations, which is
strongly recommended for any non-trivial dataset. To work on the app as well:

```bash
pip install ".[app]"
dot-explorer-app
```

See the [Development Guide](development.md) for the full contributor setup.

## Verify the installation

```python
import dot_explorer
print(dot_explorer.__version__)  # e.g. 0.1.0

from dot_explorer import SequenceIndex
idx = SequenceIndex(k=10)
idx.add_sequence("test", "ACGTACGTACGT")
print(idx)  # SequenceIndex(k=10, sequences=1)
```
