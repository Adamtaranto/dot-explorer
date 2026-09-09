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

## Quick Start — single multi-FASTA index

Each sequence added to a `SequenceIndex` gets its **own independent k-mer
index**.

Calling `add_sequence` or `load_fasta` multiple times **accumulates** sequences
— it never merges or replaces the existing collection.

Re-using an existing sequence name emits a `UserWarning` and **overwrites** that
entry.

If a FASTA file contains duplicate sequence names, `load_fasta` raises a
`ValueError` before adding any sequences.

```python
from dot_explorer import SequenceIndex
from dot_explorer.dotplot import DotPlotter

# Build an index from a multi-sequence FASTA file
# Each sequence in the file gets its own independent k-mer index entry
idx = SequenceIndex(k=15)
names = idx.load_fasta('assembly.fasta')

# load_fasta accumulates: calling it again adds more sequences, keeps existing ones
# idx.load_fasta("more_sequences.fasta")   # would add to the same index

# List the sequences now held in the index
print(idx.sequence_names())  # ['contig1', 'contig2', 'contig3', ...]

# Print all pairwise PAF lines (every i ≠ j combination)
for line in idx.get_paf_all():
    print(line)

# Print PAF lines for one specific pair
for line in idx.get_paf('contig1', 'contig2'):
    print(line)

# All-vs-all dotplot
# Forward (+) hits are drawn in blue, reverse-complement (-) hits in red.
# Sequence names appear once per column (bottom) and once per row (left).
plotter = DotPlotter(idx)
plotter.plot(output_path='all_vs_all.png', title='All vs All')

# Save as an SVG vector image instead of PNG
plotter.plot(output_path='all_vs_all.svg', title='All vs All')

# Filter out short alignments (< 500 bp) before plotting
plotter.plot(output_path='filtered.png', min_length=500)

# Single pairwise dotplot
plotter.plot_single('contig1', 'contig2', output_path='pair.png')
```

## All-vs-All Dotplot Between Two Genomes

Compare sequences from two separate FASTA files (e.g. two genome assemblies) and
plot an all-vs-all grid with subpanels scaled by relative sequence length.

```python
from dot_explorer.dotplot import DotPlotter
from dot_explorer.paf_io import CrossIndex, PafAlignment, PafRecord

# --- Build a cross-index for two assemblies ---
cross = CrossIndex(k=15)
cross.load_fasta('genome_a.fasta', group='a')  # query sequences (rows)
cross.load_fasta('genome_b.fasta', group='b')  # target sequences (columns)

# Compute and cache the cross-group matches.  This is the primary
# computation step and is required before any reordering; pass
# min_block_len= to drop short match blocks at compute time.
cross.compute_matches()

# --- Sort contigs for maximum collinearity ---
# Option 1: via CrossIndex (the d-genies gravity algorithm; its result
# matches SequenceIndex.optimal_contig_order by construction).  Each query
# contig is assigned to its best-matching target chromosome and ordered by
# its gravity centre there; reverse-oriented contigs are detected and
# exposed via cross.reversed_contigs("a").
q_sorted, t_sorted = cross.reorder_contigs()
reversed_a = cross.reversed_contigs('a')  # names to render flipped

# Option 2: via PafAlignment gravity-centre algorithm
# Retrieve all cross-group PAF lines
paf_lines = cross.get_paf_all()

records = [PafRecord.from_line(line) for line in paf_lines]
aln = PafAlignment.from_records(records)
q_sorted, t_sorted = aln.reorder_contigs(
    query_names=cross.query_names,
    target_names=cross.target_names,
)
# Unmatched contigs are placed at the end, sorted by descending length.

# --- Plot with relative scaling ---
# For a CrossIndex, select each axis by group; contig_order='colinearity'
# applies the gravity ordering at plot time, and reverse-oriented query
# contigs are auto-detected (pass reverse_contigs= to override).

plotter = DotPlotter(cross)
plotter.plot(
    query_group='a',
    target_group='b',
    contig_order='colinearity',
    output_path='cross_dotplot.png',
    scale_sequences=True,  # subplot size proportional to sequence length
    title='Genome A vs Genome B',
    # Render reverse-oriented query contigs flipped so they read along the main
    # diagonal.  Pass an explicit set, or omit to auto-pull the detected set.
    reverse_contigs=reversed_a,
)

# Save as SVG vector image for publication-quality output
plotter.plot(
    query_group='a',
    target_group='b',
    contig_order='colinearity',
    output_path='cross_dotplot.svg',
    scale_sequences=True,
    title='Genome A vs Genome B',
)

# Suppress short alignments (e.g. < 500 bp) from the plot
plotter.plot(
    query_group='a',
    target_group='b',
    contig_order='colinearity',
    output_path='cross_dotplot_filtered.png',
    scale_sequences=True,
    min_length=500,
    title='Genome A vs Genome B (≥500 bp alignments)',
)

# --- Persist the collinearity layout as FASTA ---
# Reorder + reorient assembly B against a FIXED assembly A, then write both.
# Reverse-oriented B contigs are reverse-complemented on write; A is untouched.
# (compute_matches must have been called for the (query, target) pair first.)
cross.compute_matches(query_group='b', target_group='a')
cross.reorder_for_colinearity('b', 'a', reorder_target=False)
cross.write_fasta('assembly_a.sorted.fasta', 'a')  # forward reference, unchanged
cross.write_fasta('assembly_b.sorted.fasta', 'b')  # reordered + reoriented
```

## GFF Annotation Overlays

Overlay GFF3 features on dotplots. Feature types are auto-coloured from a
qualitative palette (override per type), and a colour legend is added
automatically.

```python
from dot_explorer import DotPlotter, SequenceIndex
from dot_explorer.annotation import GffAnnotation

# From a file (gzip detected automatically), text, or raw bytes.
ann = GffAnnotation.from_file('features.gff3.gz')
ann = ann.keep_feature_types(['gene', 'CDS', 'repeat_region'])
ann.set_colors({'gene': '#2c7fb8'})

plotter = DotPlotter(idx)

# Diagonal squares: features shade self-vs-self panels behind the alignments.
plotter.plot(annotation=ann, output_path='annotated_grid.png')

# Focused single-pair view with side annotation tracks: lane-packed feature
# shapes left of the y axis and below the x axis, direction arrows for
# stranded types (gene/mRNA/exon/CDS/ORF), and connector lines joining
# multi-part CDS groups.
plotter.plot(
    query_names=['chr1'],
    target_names=['chr2'],
    annotation_query=ann,
    annotation_target=ann,
    annotation_tracks=True,
    output_path='annotated_pair.png',
)

# Interactive HTML reports make diagonal features clickable (name, type,
# coordinates, strand, parent shown in the detail bar).
plotter.to_html('report.html', annotation=ann)
```

See the dotplot tutorial for rendered examples.

## Filtering PAF Alignments by Length

Use `PafAlignment.filter_by_min_length` to remove short alignment records after
loading a PAF file. This is particularly useful for cleaned-up visualisations
when alignments have been merged from k-mer runs (which can be longer than the
k-mer size) or when working with a pre-computed PAF file.

```python
from dot_explorer.paf_io import PafAlignment

aln = PafAlignment.from_file('alignments.paf')

# Keep only alignments of at least 500 bp on the query
aln_long = aln.filter_by_min_length(500)
print(f'Records before: {len(aln)}, after: {len(aln_long)}')
```

## Writing PAF Lines to a File

```python
# All pairwise alignments within a single index
paf_lines = idx.get_paf_all()

# Or one specific pair
paf_lines = idx.get_paf('contig1', 'contig2', merge=True)

with open('alignments.paf', 'w') as f:
    for line in paf_lines:
        f.write(line + '\n')
```

## Saving and Loading Indexes

```python
# Save the current index to a compact binary file
idx.save('my_index.bin')

# Load into a new index (k must match the saved index)
idx2 = SequenceIndex(k=15)
idx2.load('my_index.bin')
```
