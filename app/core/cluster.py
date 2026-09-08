"""Clustering support for the app: provider adapters and table rows.

Pure helpers between the lazy :class:`core.seqs.SequenceProvider` world of
the app and the :mod:`dot_explorer.similarity` API, kept import-light so the
module loads even when the optional `cluster` dependencies (sourmash,
scipy) are absent.
"""

from __future__ import annotations

import importlib.util
from typing import TYPE_CHECKING, Callable, Iterable, Mapping

if TYPE_CHECKING:  # pragma: no cover - only for type checkers
    from dot_explorer.paf_io import PafRecord
    from dot_explorer.similarity import ClusterResult, SimilarityMatrix

    from .seqs import SequenceProvider


def cluster_deps_missing() -> list[str]:
    """Return which of the optional clustering dependencies are missing.

    Returns
    -------
    list of str
        Importable-module names (``sourmash``, ``scipy``) that cannot be
        found; empty when clustering is fully available.
    """
    return [
        name for name in ('sourmash', 'scipy') if importlib.util.find_spec(name) is None
    ]


class ProviderIndex:
    """Adapt a :class:`~core.seqs.SequenceProvider` to the sketching API.

    :func:`dot_explorer.similarity.compute_sketches` expects an index-like
    object with ``sequence_names()`` and ``get_sequence(name)``; providers
    expose lazy slices instead. The adapter fetches each sequence as one
    string only while it is being sketched, so peak residency stays at a
    single contig.
    """

    def __init__(self, provider: 'SequenceProvider'):
        self._provider = provider
        self._lengths = dict(provider.lengths())

    def sequence_names(self) -> list[str]:
        """Return the provider's sequence names in input order."""
        return list(self._provider.names)

    def get_sequence(self, name: str) -> str:
        """Fetch one full sequence as a plain string.

        Parameters
        ----------
        name : str
            Sequence name.

        Returns
        -------
        str
            The sequence.
        """
        return self._provider.get_slice(name, 0, self._lengths[name])


def cluster_table_rows(
    clusters: 'ClusterResult',
    lengths: Mapping[str, int],
    sim: 'SimilarityMatrix | None' = None,
) -> list[dict]:
    """Build the rows of the cluster-assignment table.

    Parameters
    ----------
    clusters : ClusterResult
        Cluster assignments keyed by contig name.
    lengths : mapping of str to int
        Contig name → length in bp.
    sim : SimilarityMatrix, optional
        When given, each row carries the contig's mean similarity to the
        other members of its cluster (blank for singletons).

    Returns
    -------
    list of dict
        One row per contig with keys ``cluster``, ``contig``, ``length``,
        ``members`` (cluster size) and ``mean_sim`` (float or None),
        grouped by cluster (clusters in first-seen order, members in
        assignment order).
    """
    rows: list[dict] = []
    for cluster_name, members in clusters.clusters.items():
        for contig in members:
            mean_sim = None
            if sim is not None and len(members) > 1:
                others = [m for m in members if m != contig]
                try:
                    mean_sim = sum(sim[(contig, m)] for m in others) / len(others)
                except ValueError:  # contig not in the matrix
                    mean_sim = None
            rows.append(
                {
                    'cluster': cluster_name,
                    'contig': contig,
                    'length': int(lengths.get(contig, 0)),
                    'members': len(members),
                    'mean_sim': mean_sim,
                }
            )
    return rows


def tree_layout_order(
    tree_leaves: list[str],
    plotted: list[str],
) -> list[str] | None:
    """Return *plotted* reordered to the tree's leaf order, if possible.

    The tree wins over every other ordering mode, but the plot may show a
    subset of the tree's tips (the min-contig-length filter): tips absent
    from *plotted* are skipped. When *plotted* contains names the tree
    does not know, the tree cannot define their position and ``None`` is
    returned (caller falls back to the non-tree order).

    Parameters
    ----------
    tree_leaves : list of str
        Tree tip names, top-to-bottom.
    plotted : list of str
        Contig names the layout wants to show.

    Returns
    -------
    list of str or None
        The reordered names, or ``None`` when the tree does not cover
        *plotted*.
    """
    plotted_set = set(plotted)
    ordered = [name for name in tree_leaves if name in plotted_set]
    if len(ordered) != len(plotted_set):
        return None
    return ordered


def merge_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    """Union possibly-overlapping half-open ``(start, end)`` intervals.

    Parameters
    ----------
    intervals : iterable of (int, int)
        Interval endpoints; ``start > end`` pairs are swapped.

    Returns
    -------
    list of (int, int)
        Sorted, non-overlapping intervals covering the same positions.
    """
    fixed = sorted(
        (min(s, e), max(s, e)) for s, e in intervals if min(s, e) != max(s, e)
    )
    merged: list[tuple[int, int]] = []
    for start, end in fixed:
        if merged and start <= merged[-1][1]:
            if end > merged[-1][1]:
                merged[-1] = (merged[-1][0], end)
            continue
        merged.append((start, end))
    return merged


def alignment_coverage_matrix(
    records: Iterable['PafRecord'],
    names: list[str],
    lengths: Mapping[str, int],
    *,
    normalize: Callable[[str], str] | None = None,
) -> 'SimilarityMatrix':
    """Pairwise coverage from alignment records (an ANI-robust containment).

    ``values[i, j]`` is the fraction of contig *i* covered by the union of
    its alignment blocks against contig *j* — from whichever aligner
    produced the current result (minimap2, nucmer, the k-mer engine, or an
    imported PAF). Unlike sourmash containment, block coverage does not
    collapse when relatives differ by scattered SNPs. Asymmetric, like
    containment: a nested fragment is fully covered by its parent, not
    vice versa.

    Parameters
    ----------
    records : iterable of PafRecord
        Alignment records for the self-comparison (both orientations of a
        pair contribute: the query side covers ``query_name``, the target
        side ``target_name``).
    names : list of str
        Display contig names, in matrix order.
    lengths : mapping of str to int
        Display name → contig length in bp.
    normalize : callable, optional
        Maps record names to display names (e.g. stripping the
        CrossIndex ``group:`` prefix). Default: identity.

    Returns
    -------
    SimilarityMatrix
        ``metric='aln_coverage'`` matrix with a unit diagonal.
    """
    import numpy as np

    from dot_explorer.similarity import SimilarityMatrix

    norm = normalize or (lambda name: name)
    per_pair: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for rec in records:
        q, t = norm(rec.query_name), norm(rec.target_name)
        per_pair.setdefault((q, t), []).append((rec.query_start, rec.query_end))
        per_pair.setdefault((t, q), []).append((rec.target_start, rec.target_end))

    n = len(names)
    values = np.eye(n)
    for i, a in enumerate(names):
        length = int(lengths.get(a, 0))
        for j, b in enumerate(names):
            if i == j:
                continue
            intervals = per_pair.get((a, b))
            if not intervals or length <= 0:
                continue
            covered = sum(end - start for start, end in merge_intervals(intervals))
            values[i, j] = min(1.0, covered / length)
    return SimilarityMatrix(names=list(names), values=values, metric='aln_coverage')
