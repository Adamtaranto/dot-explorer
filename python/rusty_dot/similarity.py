"""Pairwise sequence similarity via sourmash sketches, plus clustering.

Builds FracMinHash sketches of indexed sequences and compares them with
sourmash's Jaccard, angular (abundance-weighted cosine), containment, and
ANI estimators; hierarchical clustering and cluster assignment sit on top
via scipy. Both dependencies are optional — install with
``pip install 'rusty-dot[cluster]'``.

The code restricts itself to sourmash APIs available in 4.8.11, the
version Pyodide bundles for the browser app; native installs may run any
``sourmash>=4.8``.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

import numpy as np

if TYPE_CHECKING:  # pragma: no cover - only for type checkers
    from sourmash import MinHash

_METRICS = (
    'jaccard',
    'angular',
    'ani',
    'containment',
    'max_containment',
    'avg_containment',
)


def _require_sourmash() -> Any:
    """Import sourmash, raising a helpful error when missing."""
    try:
        import sourmash
    except ImportError as exc:  # pragma: no cover - depends on environment
        raise ImportError(
            'sourmash is required for sequence similarity sketching; '
            "install it with: pip install 'rusty-dot[cluster]'"
        ) from exc
    return sourmash


def _require_scipy_hierarchy() -> Any:
    """Import scipy.cluster.hierarchy, raising a helpful error when missing."""
    try:
        from scipy.cluster import hierarchy
    except ImportError as exc:  # pragma: no cover - depends on environment
        raise ImportError(
            'scipy is required for hierarchical clustering; '
            "install it with: pip install 'rusty-dot[cluster]'"
        ) from exc
    return hierarchy


@dataclass(frozen=True)
class SketchParams:
    """Parameters for building sourmash FracMinHash sketches.

    Attributes
    ----------
    ksize : int
        K-mer size (sourmash default for DNA comparisons is 21).
    scaled : int
        FracMinHash scaling factor: roughly one hash is kept per *scaled*
        bp. Lower values keep more hashes, improving ANI estimates on
        similar sequences at the cost of memory.
    track_abundance : bool
        Record k-mer abundances (sourmash ``-p abund``). Required for the
        angular similarity metric; other metrics can ignore it.
    """

    ksize: int = 21
    scaled: int = 1000
    track_abundance: bool = True


def compute_sketches(
    index: Any,
    names: Sequence[str] | None = None,
    *,
    group: str | None = None,
    params: SketchParams | None = None,
) -> dict[str, 'MinHash']:
    """Sketch indexed sequences with sourmash FracMinHash.

    Parameters
    ----------
    index : SequenceIndex or CrossIndex
        Index holding the sequences (anything with ``sequence_names()``
        and ``get_sequence(name)``; ``CrossIndex`` also accepts *group*).
    names : sequence of str, optional
        Sequences to sketch (default: all names in the index/group).
    group : str, optional
        For a ``CrossIndex``, restrict to this group and fetch sequences
        from it.
    params : SketchParams, optional
        Sketch parameters (k, scaled, abundance tracking); defaults to
        ``SketchParams()``.

    Returns
    -------
    dict of str to sourmash.MinHash
        One sketch per sequence name, in input order.
    """
    sourmash = _require_sourmash()
    if params is None:
        params = SketchParams()
    if names is None:
        if group is not None:
            names = index.sequence_names(group)
        else:
            names = index.sequence_names()
    sketches: dict[str, MinHash] = {}
    for name in names:
        seq = (
            index.get_sequence(name, group)
            if group is not None
            else index.get_sequence(name)
        )
        mh = sourmash.MinHash(
            n=0,
            ksize=params.ksize,
            scaled=params.scaled,
            track_abundance=params.track_abundance,
        )
        mh.add_sequence(seq, force=True)  # force: skip non-ACGT k-mers
        sketches[name] = mh
    return sketches


@dataclass
class SimilarityMatrix:
    """A named pairwise similarity matrix in ``[0, 1]``.

    Attributes
    ----------
    names : list of str
        Sequence names; row/column order of *values*.
    values : numpy.ndarray
        ``(n, n)`` similarity values. Symmetric for every metric except
        ``'containment'``, where ``values[i, j]`` is the fraction of
        sequence *i* contained in sequence *j*.
    metric : str
        One of ``'jaccard'``, ``'angular'``, ``'ani'``, ``'containment'``,
        ``'max_containment'``, ``'avg_containment'``.
    params : SketchParams
        The sketch parameters the sketches were built with.
    ci_low, ci_high : numpy.ndarray or None
        95% confidence bounds; only populated for ``metric='ani'``.
    """

    names: list[str]
    values: np.ndarray
    metric: str
    params: SketchParams = field(default_factory=SketchParams)
    ci_low: np.ndarray | None = None
    ci_high: np.ndarray | None = None

    def __getitem__(self, key: tuple[str, str]) -> float:
        """Return the similarity for a ``(name_a, name_b)`` pair."""
        a, b = key
        i, j = self.names.index(a), self.names.index(b)
        return float(self.values[i, j])

    def reorder(self, names: Sequence[str]) -> 'SimilarityMatrix':
        """Return a copy with rows/columns permuted to *names*.

        Parameters
        ----------
        names : sequence of str
            Permutation of :attr:`names`.

        Returns
        -------
        SimilarityMatrix
            Reordered copy (CI arrays permuted alongside).
        """
        if sorted(names) != sorted(self.names):
            raise ValueError('reorder names must be a permutation of the matrix names')
        idx = [self.names.index(n) for n in names]
        take = np.ix_(idx, idx)
        return SimilarityMatrix(
            names=list(names),
            values=self.values[take],
            metric=self.metric,
            params=self.params,
            ci_low=None if self.ci_low is None else self.ci_low[take],
            ci_high=None if self.ci_high is None else self.ci_high[take],
        )

    def to_csv(self, path: str | Path) -> None:
        """Write the matrix as CSV with names as header row and column.

        Parameters
        ----------
        path : str or Path
            Output file path.
        """
        with open(path, 'w', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow([''] + self.names)
            for name, row in zip(self.names, self.values):
                writer.writerow([name] + [f'{v:.6g}' for v in row])


def pairwise_similarity(
    sketches: dict[str, 'MinHash'],
    *,
    metric: str = 'jaccard',
    ignore_abundance: bool = False,
) -> SimilarityMatrix:
    """Compare sketches all-vs-all with a sourmash similarity metric.

    Parameters
    ----------
    sketches : dict of str to sourmash.MinHash
        Sketches from :func:`compute_sketches`.
    metric : str
        ``'jaccard'`` — flat Jaccard similarity (abundance ignored);
        ``'angular'`` — abundance-weighted angular similarity (requires
        sketches built with ``track_abundance=True``);
        ``'ani'`` — ANI point estimate from max containment, with 95%
        confidence bounds in :attr:`SimilarityMatrix.ci_low`/``ci_high``
        (0 where the sketches share no hashes — see the docs);
        ``'containment'`` — asymmetric: row *i*, column *j* holds the
        fraction of *i*'s hashes found in *j* (a coverage analogue);
        ``'max_containment'`` / ``'avg_containment'`` — symmetric
        containment variants.
    ignore_abundance : bool
        For ``'jaccard'`` this is implied; for ``'angular'`` it is an
        error (angular similarity is abundance-weighted by definition).

    Returns
    -------
    SimilarityMatrix
        Matrix with the diagonal set to 1.

    Raises
    ------
    ValueError
        On an unknown metric, fewer than 2 sketches, mismatched sketch
        parameters, or ``'angular'`` without abundance data.
    """
    if metric not in _METRICS:
        raise ValueError(f'unknown metric {metric!r}; choose from {_METRICS}')
    if len(sketches) < 2:
        raise ValueError('need at least 2 sketches to compare')
    if metric == 'angular':
        if ignore_abundance:
            raise ValueError(
                'angular similarity is abundance-weighted by definition; '
                "use metric='jaccard' to ignore abundance"
            )
        no_abund = [n for n, mh in sketches.items() if not mh.track_abundance]
        if no_abund:
            raise ValueError(
                'angular similarity requires abundance-tracking sketches; '
                f'built without abundance: {no_abund} '
                '(rebuild with SketchParams(track_abundance=True))'
            )

    names = list(sketches)
    n = len(names)
    values = np.eye(n)
    ci_low = np.eye(n) if metric == 'ani' else None
    ci_high = np.eye(n) if metric == 'ani' else None

    # Flatten (drop abundance) once for the metrics that do not use it;
    # containment/ANI require flat sketches in sourmash 4.8.
    if metric == 'angular':
        flat = dict(sketches)
    else:
        flat = {name: mh.flatten() for name, mh in sketches.items()}

    for i in range(n):
        for j in range(i + 1, n):
            a, b = flat[names[i]], flat[names[j]]
            if metric == 'jaccard':
                sim = a.jaccard(b, downsample=True)
            elif metric == 'angular':
                sim = a.similarity(b, ignore_abundance=False, downsample=True)
            elif metric == 'containment':
                values[i, j] = a.contained_by(b, downsample=True)
                values[j, i] = b.contained_by(a, downsample=True)
                continue
            elif metric == 'max_containment':
                sim = a.max_containment(b, downsample=True)
            elif metric == 'avg_containment':
                sim = a.avg_containment(b, downsample=True)
            else:  # ani
                est = a.max_containment_ani(b, estimate_ci=True)
                sim = est.ani if est.ani is not None else 0.0
                low = est.ani_low if est.ani_low is not None else sim
                high = est.ani_high if est.ani_high is not None else sim
                ci_low[i, j] = ci_low[j, i] = low
                ci_high[i, j] = ci_high[j, i] = high
            values[i, j] = values[j, i] = sim

    params = _sketch_params_of(next(iter(sketches.values())))
    return SimilarityMatrix(
        names=names,
        values=values,
        metric=metric,
        params=params,
        ci_low=ci_low,
        ci_high=ci_high,
    )


def _sketch_params_of(mh: 'MinHash') -> SketchParams:
    return SketchParams(
        ksize=mh.ksize, scaled=mh.scaled, track_abundance=mh.track_abundance
    )


def linkage_from_similarity(
    sim: SimilarityMatrix, *, method: str = 'average'
) -> np.ndarray:
    """Hierarchically cluster a similarity matrix.

    Parameters
    ----------
    sim : SimilarityMatrix
        Pairwise similarities; converted to distances as ``1 - values``
        (clipped to ``[0, 1]`` and symmetrised, so the asymmetric
        ``'containment'`` metric averages its two directions).
    method : str
        Linkage method for :func:`scipy.cluster.hierarchy.linkage`
        (``'average'`` = UPGMA, ``'complete'``, ``'single'``, ...).

    Returns
    -------
    numpy.ndarray
        scipy linkage matrix; feed to :meth:`Tree.from_linkage` and
        :func:`assign_clusters`.
    """
    hierarchy = _require_scipy_hierarchy()
    from scipy.spatial.distance import squareform

    dist = np.clip(1.0 - sim.values, 0.0, 1.0)
    dist = (dist + dist.T) / 2.0
    np.fill_diagonal(dist, 0.0)
    return hierarchy.linkage(squareform(dist, checks=False), method=method)


@dataclass
class ClusterResult:
    """Cluster assignments for a set of sequences.

    Attributes
    ----------
    assignments : dict of str to str
        Maps each sequence name to its cluster name (``cluster_1``, ...);
        clusters are numbered by first-seen member.
    cutoff : float
        The similarity cutoff used (identity cutoff in dual mode).
    mode : str
        ``'similarity'`` or ``'identity_coverage'``.
    metric : str
        The similarity metric the clustering was based on.
    reciprocal : bool or None
        Dual mode only: whether coverage had to pass in both directions.
    """

    assignments: dict[str, str]
    cutoff: float
    mode: str
    metric: str
    reciprocal: bool | None = None

    @property
    def clusters(self) -> dict[str, list[str]]:
        """Return cluster name → member sequence names."""
        out: dict[str, list[str]] = {}
        for seq, cluster in self.assignments.items():
            out.setdefault(cluster, []).append(seq)
        return out

    def to_csv(self, path: str | Path) -> None:
        """Write ``contig,cluster`` rows (with header) to *path*.

        Parameters
        ----------
        path : str or Path
            Output file path.
        """
        with open(path, 'w', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow(['contig', 'cluster'])
            for seq, cluster in self.assignments.items():
                writer.writerow([seq, cluster])


def _name_clusters(names: Sequence[str], labels: Sequence[int]) -> dict[str, str]:
    """Map names to cluster names, numbering clusters by first-seen member."""
    label_to_cluster: dict[int, str] = {}
    assignments: dict[str, str] = {}
    for name, label in zip(names, labels):
        if label not in label_to_cluster:
            label_to_cluster[label] = f'cluster_{len(label_to_cluster) + 1}'
        assignments[name] = label_to_cluster[label]
    return assignments


def assign_clusters(
    sim: SimilarityMatrix,
    cutoff: float,
    *,
    linkage: np.ndarray | None = None,
    method: str = 'average',
) -> ClusterResult:
    """Assign sequences to clusters at a similarity *cutoff*.

    Cuts the hierarchical clustering tree at distance ``1 - cutoff``, so
    members of a cluster have (linkage-aggregated) similarity of at
    least *cutoff*.

    Parameters
    ----------
    sim : SimilarityMatrix
        Pairwise similarities.
    cutoff : float
        Similarity threshold in ``[0, 1]``.
    linkage : numpy.ndarray, optional
        Precomputed linkage from :func:`linkage_from_similarity` (rebuilt
        from *sim* when omitted). Pass the same linkage used to draw the
        tree so cluster blocks are contiguous along the tree order.
    method : str
        Linkage method when *linkage* is omitted.

    Returns
    -------
    ClusterResult
        Assignments keyed by sequence name.
    """
    hierarchy = _require_scipy_hierarchy()
    if linkage is None:
        linkage = linkage_from_similarity(sim, method=method)
    labels = hierarchy.fcluster(linkage, t=1.0 - cutoff, criterion='distance')
    return ClusterResult(
        assignments=_name_clusters(sim.names, labels),
        cutoff=cutoff,
        mode='similarity',
        metric=sim.metric,
    )


def assign_clusters_dual(
    identity: SimilarityMatrix,
    coverage: SimilarityMatrix,
    *,
    identity_cutoff: float = 0.80,
    coverage_cutoff: float = 0.80,
    reciprocal: bool = True,
) -> ClusterResult:
    """Cluster sequences that pass identity AND coverage thresholds.

    Two sequences are linked when their identity (e.g. ANI) is at least
    *identity_cutoff* and their containment passes *coverage_cutoff* —
    in both directions when *reciprocal*, in at least one otherwise.
    Clusters are the connected components of that link graph.

    Parameters
    ----------
    identity : SimilarityMatrix
        Identity estimates (typically ``metric='ani'``).
    coverage : SimilarityMatrix
        Coverage analogue (typically the asymmetric ``'containment'``
        metric, where both directions are available).
    identity_cutoff, coverage_cutoff : float
        Thresholds in ``[0, 1]`` (default 0.80 / 0.80).
    reciprocal : bool
        Require coverage in both directions (min) instead of either
        direction (max).

    Returns
    -------
    ClusterResult
        Assignments keyed by sequence name.
    """
    if identity.names != coverage.names:
        cov = coverage.reorder(identity.names)
    else:
        cov = coverage
    names = identity.names
    n = len(names)
    cov_pair = (
        np.minimum(cov.values, cov.values.T)
        if reciprocal
        else np.maximum(cov.values, cov.values.T)
    )
    ident = (identity.values + identity.values.T) / 2.0
    linked = (ident >= identity_cutoff) & (cov_pair >= coverage_cutoff)

    # Connected components via union-find (no scipy needed here).
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(n):
        for j in range(i + 1, n):
            if linked[i, j]:
                ri, rj = find(i), find(j)
                if ri != rj:
                    parent[rj] = ri
    labels = [find(i) for i in range(n)]
    return ClusterResult(
        assignments=_name_clusters(names, labels),
        cutoff=identity_cutoff,
        mode='identity_coverage',
        metric=identity.metric,
        reciprocal=reciprocal,
    )
