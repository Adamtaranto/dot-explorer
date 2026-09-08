"""Pairwise similarity heatmap with optional clustering tree and outlines.

Renders a :class:`rusty_dot.similarity.SimilarityMatrix` as a matrix
heatmap: tree (user-supplied or from linkage) drawn left of the y-axis,
cluster membership outlined on the cells, colorbar at the right.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Sequence

if TYPE_CHECKING:  # pragma: no cover - only for type checkers
    import matplotlib.figure

    from .similarity import ClusterResult, SimilarityMatrix
    from .tree import Tree

logger = logging.getLogger(__name__)


def plot_similarity_heatmap(
    sim: 'SimilarityMatrix',
    *,
    tree: 'Tree | None' = None,
    order: Sequence[str] | None = None,
    clusters: 'ClusterResult | None' = None,
    cluster_border_color: str = 'black',
    cluster_border_lw: float = 2.0,
    cmap: str = 'viridis',
    colorbar: bool = True,
    cutoff: float | None = None,
    tree_scalebar: bool = True,
    tree_width: float = 1.2,
    figsize: tuple[float, float] = (6.0, 6.0),
    annotate: bool = False,
    title: str | None = None,
    output_path: str | None = None,
    dpi: int = 150,
    format: str | None = None,
) -> 'matplotlib.figure.Figure':
    """Plot a pairwise similarity matrix as a heatmap.

    Parameters
    ----------
    sim : SimilarityMatrix
        Pairwise similarities (needs at least 2 sequences).
    tree : Tree, optional
        Drawn left of the y-axis; rows/columns are permuted to its leaf
        order. Tip labels must match ``sim.names`` exactly.
    order : sequence of str, optional
        Explicit row/column order; mutually exclusive with *tree*.
    clusters : ClusterResult, optional
        Draw a bold outline around each cluster's block of cells. A
        cluster whose members are not contiguous in the display order is
        outlined per contiguous run (with a log warning).
    cluster_border_color : str
        Outline colour.
    cluster_border_lw : float
        Outline line width.
    cmap : str
        Matplotlib colormap name.
    colorbar : bool
        Draw a colour scale bar to the right of the heatmap.
    cutoff : float, optional
        Similarity cutoff, drawn as a dashed line through the tree at
        distance ``cutoff`` from the tips (ultrametric trees).
    tree_scalebar : bool
        Show a branch-length scale bar under the tree.
    tree_width : float
        Width of the tree gutter in inches.
    figsize : tuple of float
        Size of the heatmap portion in inches (tree/colorbar gutters are
        added on top of this).
    annotate : bool
        Write the similarity value in each cell (readable for small
        matrices only).
    title : str, optional
        Figure title.
    output_path : str, optional
        When given, save the figure (SVG/PDF/PNG by extension). The HTML
        report format is not supported for heatmaps.
    dpi : int
        Raster resolution when saving.
    format : str, optional
        Explicit save format overriding the extension.

    Returns
    -------
    matplotlib.figure.Figure
        The rendered figure.

    Raises
    ------
    ValueError
        With fewer than 2 sequences, when both *tree* and *order* are
        given, or when tree tips do not match the matrix names.
    """
    import matplotlib.pyplot as plt

    from .tree import draw_tree

    if len(sim.names) < 2:
        raise ValueError('heatmap needs at least 2 sequences to compare')
    if tree is not None and order is not None:
        raise ValueError('tree fixes the display order; do not also pass order')
    if output_path is not None and str(output_path).lower().endswith(('.html', '.htm')):
        raise ValueError('HTML output is not supported for similarity heatmaps')

    if tree is not None:
        tree.validate_labels(sim.names)
        display = sim.reorder(tree.leaf_names())
    elif order is not None:
        display = sim.reorder(list(order))
    else:
        display = sim

    n = len(display.names)
    hm_w, hm_h = figsize
    widths = [hm_w]
    if tree is not None:
        widths.insert(0, tree_width)
    if colorbar:
        widths.append(0.25)
    fig_w = sum(widths) + 0.4 * (len(widths) - 1)
    fig = plt.figure(figsize=(fig_w, hm_h))
    gs = fig.add_gridspec(1, len(widths), width_ratios=widths, wspace=0.15)

    col = 0
    tree_ax = None
    if tree is not None:
        tree_ax = fig.add_subplot(gs[0, col])
        col += 1
    ax = fig.add_subplot(gs[0, col])
    col += 1
    cax = fig.add_subplot(gs[0, col]) if colorbar else None

    image = ax.imshow(
        display.values,
        cmap=cmap,
        vmin=0.0,
        vmax=1.0,
        aspect='auto',
        interpolation='nearest',
    )
    image.set_gid('rd-heatmap')
    ax.set_xticks(range(n), display.names, rotation=90, fontsize='small')
    ax.set_yticks(range(n), display.names, fontsize='small')
    if tree is not None:
        # The tree carries the row identity; avoid doubled labels.
        ax.set_yticks([])

    if annotate:
        for i in range(n):
            for j in range(n):
                value = display.values[i, j]
                ax.text(
                    j,
                    i,
                    f'{value:.2f}',
                    ha='center',
                    va='center',
                    fontsize='x-small',
                    color='white' if value < 0.5 else 'black',
                )

    if clusters is not None:
        _outline_clusters(
            ax,
            display.names,
            clusters,
            color=cluster_border_color,
            lw=cluster_border_lw,
        )

    if tree_ax is not None and tree is not None:
        tree_ax.set_ylim(ax.get_ylim())  # imshow: (n-0.5, -0.5), rows aligned
        draw_tree(
            tree_ax,
            tree,
            {name: i for i, name in enumerate(display.names)},
            cutoff=cutoff,
            scalebar=tree_scalebar,
            leaf_labels=True,
        )

    if cax is not None:
        cbar = fig.colorbar(image, cax=cax)
        cbar.set_label(_metric_label(display.metric))
        cbar.ax.set_gid('rd-heatmap-scale')

    if title:
        fig.suptitle(title)

    if output_path is not None:
        fig.savefig(output_path, dpi=dpi, format=format, bbox_inches='tight')
    return fig


_METRIC_LABELS = {
    'jaccard': 'Jaccard similarity',
    'angular': 'Angular similarity',
    'ani': 'ANI',
    'containment': 'Containment',
    'max_containment': 'Max containment',
    'avg_containment': 'Average containment',
}


def _metric_label(metric: str) -> str:
    return _METRIC_LABELS.get(metric, metric)


def _outline_clusters(
    ax,
    names: Sequence[str],
    clusters: 'ClusterResult',
    *,
    color: str,
    lw: float,
) -> None:
    """Outline each cluster's contiguous runs of rows/columns."""
    from matplotlib.patches import Rectangle

    index_of = {name: i for i, name in enumerate(names)}
    for cluster_name, members in clusters.clusters.items():
        rows = sorted(index_of[m] for m in members if m in index_of)
        if not rows:
            continue
        runs = _contiguous_runs(rows)
        if len(runs) > 1:
            logger.warning(
                'cluster %s is not contiguous in the display order; '
                'outlining %d separate blocks',
                cluster_name,
                len(runs),
            )
        for start, stop in runs:
            rect = Rectangle(
                (start - 0.5, start - 0.5),
                stop - start + 1,
                stop - start + 1,
                fill=False,
                edgecolor=color,
                linewidth=lw,
                zorder=5,
            )
            rect.set_gid(f'rd-hm-cluster-{cluster_name}')
            ax.add_patch(rect)


def _contiguous_runs(indices: list[int]) -> list[tuple[int, int]]:
    """Split sorted *indices* into inclusive (start, stop) runs."""
    runs: list[tuple[int, int]] = []
    start = prev = indices[0]
    for idx in indices[1:]:
        if idx == prev + 1:
            prev = idx
            continue
        runs.append((start, prev))
        start = prev = idx
    runs.append((start, prev))
    return runs
