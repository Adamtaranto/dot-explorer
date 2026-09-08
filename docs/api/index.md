# API Reference

dot-explorer exposes its functionality through the following classes and functions.

## Classes

| Class | Module | Description |
|-------|--------|-------------|
| [`SequenceIndex`](sequence_index.md) | `dot_explorer` | Rust-backed rolling-hash k-mer index for sequence comparison |
| [`DotPlotter`](dotplot.md) | `dot_explorer.dotplot` | All-vs-all dotplot visualisation |
| [`CrossIndex`](cross_index.md) | `dot_explorer.paf_io` | Multi-group cross-index for cross-group pairwise comparisons; DotPlotter-compatible |
| [`PafRecord`](paf_io.md#dot_explorer.paf_io.PafRecord) | `dot_explorer.paf_io` | Single PAF alignment record |
| [`PafAlignment`](paf_io.md#dot_explorer.paf_io.PafAlignment) | `dot_explorer.paf_io` | Collection of PAF records with reordering utilities; DotPlotter-compatible |

## Functions

| Function | Module | Description |
|----------|--------|-------------|
| [`py_read_fasta`](functions.md#dot_explorer._dot_explorer.py_read_fasta) | `dot_explorer` | Read a FASTA or gzipped FASTA file |
| [`py_build_kmer_set`](functions.md#dot_explorer._dot_explorer.py_build_kmer_set) | `dot_explorer` | Build the k-mer set for a sequence |
| [`py_find_kmer_coords`](functions.md#dot_explorer._dot_explorer.py_find_kmer_coords) | `dot_explorer` | Find k-mer positions in a sequence via FM-index |
| [`py_merge_runs`](functions.md#dot_explorer._dot_explorer.py_merge_runs) | `dot_explorer._dot_explorer` | Unified strand-aware merge: forward and both RC patterns |
| [`py_merge_kmer_runs`](functions.md#dot_explorer._dot_explorer.py_merge_kmer_runs) | `dot_explorer` | Merge forward-strand (`+`) co-linear k-mer hits into blocks |
| [`py_merge_rev_runs`](functions.md#dot_explorer._dot_explorer.py_merge_rev_runs) | `dot_explorer._dot_explorer` | Merge RC anti-diagonal k-mer hits (standard inverted repeat) |
| [`py_merge_rev_fwd_runs`](functions.md#dot_explorer._dot_explorer.py_merge_rev_fwd_runs) | `dot_explorer._dot_explorer` | Merge RC co-diagonal k-mer hits (both arms same direction) |
| [`py_coords_to_paf`](functions.md#dot_explorer._dot_explorer.py_coords_to_paf) | `dot_explorer` | Convert coordinate tuples to PAF lines |
| [`py_save_index`](functions.md#dot_explorer._dot_explorer.py_save_index) | `dot_explorer` | Serialise an index collection to disk |
| [`py_load_index`](functions.md#dot_explorer._dot_explorer.py_load_index) | `dot_explorer` | Load a serialised index from disk |
| [`parse_paf_file`](paf_io.md#dot_explorer.paf_io.parse_paf_file) | `dot_explorer.paf_io` | Yield PAF records from a file |
| [`compute_gravity_contigs`](paf_io.md#dot_explorer.paf_io.compute_gravity_contigs) | `dot_explorer.paf_io` | Sort contigs by best-chromosome gravity centre; report reverse-oriented contigs |
| [`compute_reversed_contigs`](paf_io.md#dot_explorer.paf_io.compute_reversed_contigs) | `dot_explorer.paf_io` | Detect reverse-oriented query contigs (d-genies orientation check) |
| [`reverse_complement`](paf_io.md#dot_explorer.paf_io.reverse_complement) | `dot_explorer.paf_io` | Reverse-complement a nucleotide string |
