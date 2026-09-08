# Low-Level Functions

These functions are implemented in Rust and exposed via PyO3.
They provide fine-grained access to the underlying k-mer machinery.  (Note that
`py_find_kmer_coords` performs a one-off FM-index search and is independent of the
rolling-hash index used by `SequenceIndex` for comparisons.)
For most use cases, the [`SequenceIndex`](sequence_index.md) class is more convenient.

## FASTA I/O

::: dot_explorer._dot_explorer.py_read_fasta

## K-mer Operations

::: dot_explorer._dot_explorer.py_build_kmer_set

::: dot_explorer._dot_explorer.py_find_kmer_coords

## Merging K-mer Runs

dot-explorer provides four merge functions covering all k-mer alignment orientations.
`py_merge_runs` is the recommended entry-point for new code; the strand-specific
functions are available for lower-level control.

### Unified entry-point

::: dot_explorer._dot_explorer.py_merge_runs

### Forward-strand merge

::: dot_explorer._dot_explorer.py_merge_kmer_runs

### Reverse-complement merges

Two complementary algorithms cover all reverse-complement alignment patterns:

| Pattern | When to use |
|---------|-------------|
| `py_merge_rev_runs` | RC target positions *decrease* as query advances (query +1, target −1 per step) — standard inverted-repeat alignment where the two arms face each other |
| `py_merge_rev_fwd_runs` | RC target positions *increase* as query advances (query +1, target +1 per step) — both repeat arms run in the same left-to-right direction |

`py_merge_runs(strand="-")` calls both and deduplicates the results automatically.

::: dot_explorer._dot_explorer.py_merge_rev_runs

::: dot_explorer._dot_explorer.py_merge_rev_fwd_runs

## PAF Formatting

::: dot_explorer._dot_explorer.py_coords_to_paf

## Index Serialization

::: dot_explorer._dot_explorer.py_save_index

::: dot_explorer._dot_explorer.py_load_index
