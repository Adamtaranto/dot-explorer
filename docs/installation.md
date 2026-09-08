# Installation

## Requirements

- Python 3.9 – 3.14
- A working Rust toolchain (for building from source)
- [maturin](https://www.maturin.rs/) ≥ 1.0

## Install from source

Clone the repository and build the Rust extension with [maturin](https://www.maturin.rs/):

```bash
git clone https://github.com/Adamtaranto/dot-explorer.git
cd dot-explorer
pip install maturin
maturin develop --release
```

The `--release` flag enables full Rust compiler optimisations, which is strongly recommended for any non-trivial dataset.

## Install Python dependencies

dot-explorer depends on:

| Package | Purpose |
|---------|---------|
| `matplotlib ≥ 3.5` | Dotplot visualisation |
| `numpy ≥ 1.21` | Array operations used by matplotlib |

These are declared as package dependencies and will be installed automatically by pip.

## Optional: documentation dependencies

To build the documentation locally:

```bash
pip install dot-explorer[docs]
python scripts/notebooks_to_md.py   # render the tutorial notebooks
zensical serve
```

## Verify the installation

```python
import dot_explorer
print(dot_explorer.__version__)  # 0.1.0

from dot_explorer import SequenceIndex
idx = SequenceIndex(k=10)
idx.add_sequence("test", "ACGTACGTACGT")
print(idx)  # SequenceIndex(k=10, sequences=1)
```
