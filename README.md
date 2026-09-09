[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# dot-explorer

A python library with Shiny Web App frontend for fast dot plot comparisons of DNA sequences. Written in Rust with PyO3 python bindings.

## Browser app

[![Web App](https://img.shields.io/badge/Web%20App-live-teal)](https://adamtaranto.github.io/dot-explorer/app/)

Try dot-explorer without installing anything:
**[Launch Dot-Explorer Live](https://adamtaranto.github.io/dot-explorer/app/)**.

The app runs entirely in your browser (WebAssembly + Pyodide) — uploaded
assemblies never leave your machine. Align with the k-mer engine, minimap2,
or nucmer (or import a PAF file), overlay GFF3 annotations, and download SVG/PDF
plots, PAF alignments, and a reordered/reoriented query FASTA.

In-browser computation is memory limited (the wasm heap caps at ~4 GB, so
the k-mer method is gated above ~80 Mb of combined input); larger genomes
should use minimap2 instead of k-mer based plotting. Alternatively, run the
same app locally with no size limits (`pip install "dot-explorer[app]"`, then
`dot-explorer-app`), or generate plots using the Python API locally or on
Google Colab.

See [docs/webapp.md](https://adamtaranto.github.io/dot-explorer/webapp/) for capabilities, limits, and the local/HPC run guide, and
[Python library tutorials](https://adamtaranto.github.io/dot-explorer/tutorials/quickstart/) to run analysis locally.

## Installation

Requires Python >=3.12 <3.15. Binary wheels are published for 3.12, 3.13, 3.14
and 3.14t on Linux, macOS and Windows, so pip never has to compile Rust.

```bash
# Python library
pip install dot-explorer

# Library plus the local browser app
pip install "dot-explorer[app]"

# Optional: similarity metrics, clustering and trees (sourmash + scipy)
pip install "dot-explorer[cluster]"
```

## Run the app locally

```bash
dot-explorer-app
```

Serves the app on <http://127.0.0.1:8000> and opens a browser. Use `--host`,
`--port` and `--no-browser` when running on a remote machine or HPC node — see
the [web app guide](https://adamtaranto.github.io/dot-explorer/webapp/) for the
SSH-tunnel recipe and for where the app writes temporary files.

## Build from source

Only needed to modify the Rust core.

Requirements:

- Rust: See [rust-lang.org](https://rust-lang.org/tools/install/)
- Python >=3.12 <3.15

```bash
# Clone this project repo
git clone https://github.com/Adamtaranto/dot-explorer.git && cd dot-explorer

# Install maturin build tool
pip install maturin

# Build and install the python package
maturin develop --release
```

See the [development guide](https://adamtaranto.github.io/dot-explorer/development/)
for the full contributor setup.
