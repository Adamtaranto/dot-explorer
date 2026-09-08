# PAF I/O

This module provides classes and helpers for reading, writing, and reordering
[PAF (Pairwise mApping Format)](https://github.com/lh3/miniasm/blob/master/PAF.md)
alignment records.

## PafAlignment — Alignment record collection

`PafAlignment` wraps a list of [`PafRecord`](#dot_explorer.paf_io.PafRecord)
objects and provides filtering, contig reordering, and sequence-length lookup
utilities.  It can be passed directly to
[`DotPlotter`](dotplot.md#dot_explorer.dotplot.DotPlotter) — no
[`SequenceIndex`](sequence_index.md#dot_explorer._dot_explorer.SequenceIndex) is
required:

```python
from dot_explorer.paf_io import PafAlignment
from dot_explorer.dotplot import DotPlotter

aln = PafAlignment.from_file("alignments.paf")
q_order, t_order = aln.reorder_contigs()

plotter = DotPlotter(aln)
plotter.plot(
    query_names=q_order,
    target_names=t_order,
    output_path="dotplot.png",
    scale_sequences=True,
)
```

::: dot_explorer.paf_io.PafRecord

::: dot_explorer.paf_io.PafAlignment

## Functions

::: dot_explorer.paf_io.parse_paf_file

::: dot_explorer.paf_io.compute_gravity_contigs

::: dot_explorer.paf_io.compute_reversed_contigs

::: dot_explorer.paf_io.reverse_complement
