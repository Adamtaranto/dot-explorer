"""dot-explorer: fast dot plots of DNA via a canonical ntHash k-mer index.

This package provides:
- Rust-backed canonical-hash k-mer index for fast lookup (via PyO3 bindings)
- Python classes for building and querying sequence indexes
- Dotplot visualization tools
- PAF file I/O with CIGAR string parsing
- Gravity-based contig reordering for maximising dotplot collinearity

Examples
--------
Basic usage:

>>> from dot_explorer import SequenceIndex
>>> idx = SequenceIndex(k=10)
>>> idx.add_sequence("seq1", "ACGTACGTACGT")
>>> idx.add_sequence("seq2", "TACGTACGTACG")
>>> matches = idx.compare_sequences_stranded("seq1", "seq2")
"""

from dot_explorer._dot_explorer import (  # noqa: F401
    SequenceIndex,
    py_build_kmer_set,
    py_coords_to_paf,
    py_find_kmer_coords,
    py_load_index,
    py_merge_kmer_runs,
    py_save_index,
)

try:
    from dot_explorer._dot_explorer import py_read_fasta  # noqa: F401
except ImportError:  # pragma: no cover - wasm-only branch
    # wasm32-unknown-emscripten wheels are built --no-default-features, which
    # excludes the needletail-backed FASTA reader (the `fasta` cargo feature).
    # Everything else works; in-browser callers parse FASTA in Python and use
    # add_sequence() instead.
    py_read_fasta = None  # type: ignore[assignment]
from dot_explorer.annotation import GffAnnotation, GffFeature  # noqa: F401
from dot_explorer.dotplot import DotPlotter  # noqa: F401
from dot_explorer.heatmap import plot_similarity_heatmap  # noqa: F401
from dot_explorer.paf_io import (  # noqa: F401
    CrossIndex,
    PafAlignment,
    PafRecord,
    compute_gravity_contigs,
    compute_reversed_contigs,
    parse_paf_file,
    reverse_complement,
)
from dot_explorer.similarity import (  # noqa: F401
    ClusterResult,
    SimilarityMatrix,
    SketchParams,
    assign_clusters,
    assign_clusters_dual,
    compute_sketches,
    linkage_from_similarity,
    pairwise_similarity,
)
from dot_explorer.style import NATURE_RC, nature_style, use_nature_style  # noqa: F401
from dot_explorer.tree import Tree, TreeNode, draw_tree  # noqa: F401

__version__ = '0.1.0'
__all__ = [
    'SequenceIndex',
    'DotPlotter',
    'CrossIndex',
    'PafRecord',
    'PafAlignment',
    'parse_paf_file',
    'compute_gravity_contigs',
    'compute_reversed_contigs',
    'reverse_complement',
    'GffAnnotation',
    'GffFeature',
    'py_read_fasta',
    'py_build_kmer_set',
    'py_find_kmer_coords',
    'py_merge_kmer_runs',
    'py_coords_to_paf',
    'py_save_index',
    'py_load_index',
    'NATURE_RC',
    'nature_style',
    'use_nature_style',
    'Tree',
    'TreeNode',
    'draw_tree',
    'SketchParams',
    'SimilarityMatrix',
    'ClusterResult',
    'compute_sketches',
    'pairwise_similarity',
    'linkage_from_similarity',
    'assign_clusters',
    'assign_clusters_dual',
    'plot_similarity_heatmap',
]
