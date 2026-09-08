"""Clustering support for the app: provider adapters and table rows.

Pure helpers between the lazy :class:`core.seqs.SequenceProvider` world of
the app and the :mod:`rusty_dot.similarity` API, kept import-light so the
module loads even when the optional `cluster` dependencies (sourmash,
scipy) are absent.
"""

from __future__ import annotations

import importlib.util
from typing import TYPE_CHECKING, Mapping

if TYPE_CHECKING:  # pragma: no cover - only for type checkers
    from rusty_dot.similarity import ClusterResult, SimilarityMatrix

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

    :func:`rusty_dot.similarity.compute_sketches` expects an index-like
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
