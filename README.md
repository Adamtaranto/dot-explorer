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
`shiny run --launch-browser app/app.py` from a checkout), or generate plots
using the Python API locally or on Google Colab.

See [docs/webapp.md](https://adamtaranto.github.io/dot-explorer/webapp/) for capabilities, limits, and the local/HPC run guide, and
[Python library tutorials](https://adamtaranto.github.io/dot-explorer/tutorials/quickstart/) to run analysis locally.

## Installation

Requirements:

- Rust: See [rust-lang.org](https://rust-lang.org/tools/install/)
- Python >=3.12 <3.15 (binary wheels are published for 3.12, 3.13, 3.14 and 3.14t)

```bash
# Clone this project repo
git clone https://github.com/Adamtaranto/dot-explorer.git && cd dot-explorer

# Install maturin build tool
pip install maturin

# Build and install the python package
maturin develop --release
```
