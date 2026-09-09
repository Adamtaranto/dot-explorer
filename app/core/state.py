"""Plot configuration state for the browser app.

A single dataclass captures every display option that can change without
recomputing an alignment.  The app re-renders from cached alignment results
whenever this config changes, and remembers it across sub-plot drill-down.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

#: Contig-ordering choices exposed in the UI, mapped to display labels.
ORDER_CHOICES: dict[str, str] = {
    'input': 'Input order',
    'length': 'Sort by size',
    'colinearity': 'Maximise colinearity',
    'colinearity_ref': 'Colinearity vs fixed reference',
}

#: Ordering modes that derive orientation info (reversed contigs).
COLINEARITY_MODES: frozenset[str] = frozenset({'colinearity', 'colinearity_ref'})

#: Match line-cap choices exposed in the UI (matplotlib names), mapped to
#: display labels.  Square is the default: with flat caps a match shorter than
#: the line width draws wider across its diagonal than along it, so it reads as
#: a mark rotated 90 degrees.
CAP_STYLE_CHOICES: dict[str, str] = {
    'projecting': 'Square',
    'round': 'Round',
    'butt': 'Flat',
}

#: matplotlib cap name -> SVG ``stroke-linecap`` value.
_SVG_LINECAPS: dict[str, str] = {
    'butt': 'butt',
    'round': 'round',
    'projecting': 'square',
}


def svg_linecap(cap_style: str) -> str:
    """Map a matplotlib cap style to its SVG ``stroke-linecap`` value.

    The embedded report applies the cap as a CSS override, which speaks SVG's
    vocabulary: matplotlib's ``'projecting'`` is SVG's ``'square'``.

    Parameters
    ----------
    cap_style : str
        A key of :data:`CAP_STYLE_CHOICES`.

    Returns
    -------
    str
        The corresponding ``stroke-linecap`` value; unknown input falls back
        to ``'square'`` (the default cap).
    """
    return _SVG_LINECAPS.get(cap_style, 'square')


@dataclass
class PlotConfig:
    """Display options for the dotplot, independent of the alignment.

    Attributes
    ----------
    contig_order : str
        One of ``'input'``, ``'length'``, ``'colinearity'`` or
        ``'colinearity_ref'`` (keys of :data:`ORDER_CHOICES`).
        ``'colinearity_ref'`` keeps the target axis fixed and reorders only
        the query contigs against it.
    auto_reverse : bool
        Flip contigs detected as reverse-oriented so colinear alignments
        render as forward diagonals.
    hide_internal_axes : bool
        Remove internal ticks/spines so the grid reads as one continuous
        plot.
    dot_size : float
        Line width for match segments.
    cap_style : str
        Line cap for match segments — a key of :data:`CAP_STYLE_CHOICES`
        (matplotlib names: ``'butt'``, ``'round'`` or ``'projecting'``).
    min_length : int
        Minimum match length (bp) to draw.
    title : str or None
        Optional figure title.
    color_by_identity : bool
        Colour match segments by percent identity instead of the default
        forward/reverse colouring (meaningful for tool/PAF alignments only;
        k-mer matches are always 100% identity).
    identity_palette : str
        Matplotlib colormap name used when *color_by_identity* is on.
    """

    contig_order: str = 'input'
    auto_reverse: bool = False
    hide_internal_axes: bool = False
    dot_size: float = 0.5
    cap_style: str = 'projecting'
    min_length: int = 0
    title: str | None = field(default=None)
    color_by_identity: bool = False
    identity_palette: str = 'viridis'

    def plot_kwargs(self) -> dict[str, Any]:
        """Translate this config into ``DotPlotter.plot`` keyword arguments.

        ``'colinearity_ref'`` maps to ``contig_order=None`` because the
        target-fixed reordering is not a plot-level strategy: the app
        computes explicit name orders (see :func:`core.panels.resolve_orders`)
        and passes ``query_names=`` / ``target_names=`` /
        ``reverse_contigs=`` itself.

        Returns
        -------
        dict[str, Any]
            Keyword arguments accepted by
            :meth:`dot_explorer.dotplot.DotPlotter.plot`.
        """
        order = (
            self.contig_order
            if self.contig_order in ('length', 'colinearity')
            else None
        )
        return {
            'contig_order': order,
            'auto_reverse': self.auto_reverse,
            'hide_internal_axes': self.hide_internal_axes,
            'dot_size': self.dot_size,
            'cap_style': self.cap_style,
            'min_length': self.min_length,
            'title': self.title,
            'color_by_identity': self.color_by_identity,
            'identity_palette': self.identity_palette,
        }
