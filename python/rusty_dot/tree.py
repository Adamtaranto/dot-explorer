"""Phylogenetic / clustering trees for dot-plot and heatmap axes.

Provides a minimal, dependency-free newick parser (compatible with plain
newick files and IQ-TREE ``.treefile`` output), a converter from
:func:`scipy.cluster.hierarchy.linkage` matrices, and a matplotlib
renderer that draws a rectangular dendrogram alongside a plot axis.

Trees fix the row order of a dot-plot matrix or similarity heatmap:
``Tree.leaf_names()`` (left to right) maps to plot rows (top to bottom).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Iterator, Sequence

if TYPE_CHECKING:  # pragma: no cover - only for type checkers
    import matplotlib.axes
    import numpy as np

__all__ = ['Tree', 'TreeNode', 'draw_tree']


@dataclass
class TreeNode:
    """A node in a :class:`Tree`.

    Attributes
    ----------
    name : str or None
        Leaf label, or internal-node label when present (IQ-TREE writes
        support values here).
    length : float or None
        Branch length to the parent, or ``None`` when the newick source
        omitted it.
    support : float or None
        Numeric support value parsed from an internal-node label
        (``)95:`` or ``)label/95:`` styles), when recognisable.
    children : list of TreeNode
        Child nodes; empty for leaves.
    """

    name: str | None = None
    length: float | None = None
    support: float | None = None
    children: list['TreeNode'] = field(default_factory=list)

    @property
    def is_leaf(self) -> bool:
        """Return True when this node has no children."""
        return not self.children


class Tree:
    """A rooted tree with named leaves.

    Construct via :meth:`from_newick`, :meth:`read`, or
    :meth:`from_linkage`; the leaf order (:meth:`leaf_names`) defines the
    plot row order when the tree is passed to
    :meth:`rusty_dot.DotPlotter.plot` or
    :func:`rusty_dot.plot_similarity_heatmap`.
    """

    def __init__(self, root: TreeNode):
        self.root = root

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    @classmethod
    def from_newick(cls, text: str) -> 'Tree':
        """Parse a newick string (IQ-TREE ``.treefile`` compatible).

        Handles quoted labels (``'...'`` with ``''`` escapes), bracketed
        comments, missing branch lengths, and internal support labels in
        both the ``)95:0.1`` and ``)label/95:0.1`` styles.

        Parameters
        ----------
        text : str
            Newick source; everything after the first ``;`` is ignored.

        Returns
        -------
        Tree
            The parsed tree.

        Raises
        ------
        ValueError
            On malformed newick (unbalanced parentheses, unterminated
            quotes, empty input, or duplicate leaf names).
        """
        root = _parse_newick(text)
        tree = cls(root)
        names = tree.leaf_names()
        if not names:
            raise ValueError('newick tree has no labelled leaves')
        from collections import Counter

        dupes = sorted(n for n, c in Counter(names).items() if c > 1)
        if dupes:
            raise ValueError(f'duplicate leaf names in tree: {dupes}')
        return tree

    @classmethod
    def read(cls, path: str | Path) -> 'Tree':
        """Read a newick tree from *path* (``.nwk``/``.newick``/``.treefile``/``.txt``).

        Parameters
        ----------
        path : str or Path
            File containing a single newick tree.

        Returns
        -------
        Tree
            The parsed tree.
        """
        return cls.from_newick(Path(path).read_text())

    @classmethod
    def from_linkage(cls, Z: 'np.ndarray', labels: Sequence[str]) -> 'Tree':
        """Build a tree from a scipy ``linkage`` matrix.

        Branch lengths are derived from merge heights, so cutting the
        tree at distance ``d`` from the tips matches
        ``scipy.cluster.hierarchy.fcluster(Z, t=d, criterion='distance')``.
        The leaf order matches scipy's default (unsorted) ``dendrogram``.

        Parameters
        ----------
        Z : numpy.ndarray
            ``(n - 1, 4)`` linkage matrix from
            :func:`scipy.cluster.hierarchy.linkage`.
        labels : sequence of str
            Names of the ``n`` original observations, in the order they
            were given to ``linkage``.

        Returns
        -------
        Tree
            Ultrametric tree with one leaf per label.
        """
        n = len(labels)
        if len(Z) != n - 1:
            raise ValueError(
                f'linkage matrix has {len(Z)} rows; expected {n - 1} for {n} labels'
            )
        nodes: dict[int, TreeNode] = {
            i: TreeNode(name=str(labels[i])) for i in range(n)
        }
        heights: dict[int, float] = dict.fromkeys(range(n), 0.0)
        for k, (a, b, height, _count) in enumerate(Z):
            a, b = int(a), int(b)
            parent = TreeNode(children=[nodes[a], nodes[b]])
            for child in (a, b):
                nodes[child].length = float(height) - heights[child]
            nodes[n + k] = parent
            heights[n + k] = float(height)
        return cls(nodes[n + len(Z) - 1])

    # ------------------------------------------------------------------
    # Queries
    # ------------------------------------------------------------------

    def _postorder(self) -> Iterator[TreeNode]:
        """Yield nodes children-first (iterative; safe for deep trees)."""
        stack: list[tuple[TreeNode, bool]] = [(self.root, False)]
        while stack:
            node, expanded = stack.pop()
            if expanded:
                yield node
            else:
                stack.append((node, True))
                for child in reversed(node.children):
                    stack.append((child, False))

    def leaves(self) -> list[TreeNode]:
        """Return the leaf nodes in left-to-right (plot top-to-bottom) order."""
        out = []
        stack = [self.root]
        while stack:
            node = stack.pop()
            if node.is_leaf:
                out.append(node)
            else:
                stack.extend(reversed(node.children))
        return out

    def leaf_names(self) -> list[str]:
        """Return leaf labels in left-to-right (plot top-to-bottom) order."""
        return [leaf.name or '' for leaf in self.leaves()]

    def max_depth(self) -> float:
        """Return the largest root-to-leaf cumulative branch length.

        When the tree carries no branch lengths at all, every edge counts
        as 1.0; otherwise missing lengths count as 0.0.
        """
        default = 1.0 if not self._has_lengths() else 0.0
        deepest = 0.0
        stack: list[tuple[TreeNode, float]] = [(self.root, 0.0)]
        while stack:
            node, depth = stack.pop()
            if node.is_leaf:
                deepest = max(deepest, depth)
            for child in node.children:
                length = child.length if child.length is not None else default
                stack.append((child, depth + length))
        return deepest

    def _has_lengths(self) -> bool:
        return any(node.length is not None for node in self._postorder())

    def validate_labels(self, seq_names: Sequence[str]) -> None:
        """Check that tree tips and sequence names match exactly.

        Parameters
        ----------
        seq_names : sequence of str
            The sequence names the tree will be plotted against.

        Raises
        ------
        ValueError
            When any tip lacks a sequence or any sequence lacks a tip;
            the message lists both directions of the mismatch.
        """
        tips = set(self.leaf_names())
        seqs = set(seq_names)
        problems = []
        extra_tips = sorted(tips - seqs)
        if extra_tips:
            problems.append(f'tree tips with no matching sequence: {extra_tips}')
        missing_tips = sorted(seqs - tips)
        if missing_tips:
            problems.append(f'sequences with no tip in the tree: {missing_tips}')
        if problems:
            raise ValueError(
                'tree tip labels do not match sequence names — ' + '; '.join(problems)
            )


# ----------------------------------------------------------------------
# Newick parsing
# ----------------------------------------------------------------------


def _strip_comments(text: str) -> str:
    """Remove ``[...]`` comments, leaving quoted sections untouched."""
    out: list[str] = []
    i, n = 0, len(text)
    while i < n:
        ch = text[i]
        if ch == "'":
            j = i + 1
            while j < n:
                if text[j] == "'":
                    if j + 1 < n and text[j + 1] == "'":  # escaped quote
                        j += 2
                        continue
                    break
                j += 1
            else:
                raise ValueError('unterminated quoted label in newick')
            out.append(text[i : j + 1])
            i = j + 1
        elif ch == '[':
            depth = 1
            j = i + 1
            while j < n and depth:
                if text[j] == '[':
                    depth += 1
                elif text[j] == ']':
                    depth -= 1
                j += 1
            if depth:
                raise ValueError('unterminated [comment] in newick')
            i = j
        else:
            out.append(ch)
            i += 1
    return ''.join(out)


def _split_chunk(chunk: str) -> tuple[str | None, float | None]:
    """Split a ``label:length`` chunk, honouring quoted labels."""
    chunk = chunk.strip()
    if not chunk:
        return None, None
    if chunk.startswith("'"):
        j = 1
        n = len(chunk)
        label_parts: list[str] = []
        while j < n:
            if chunk[j] == "'":
                if j + 1 < n and chunk[j + 1] == "'":
                    label_parts.append("'")
                    j += 2
                    continue
                break
            label_parts.append(chunk[j])
            j += 1
        else:
            raise ValueError(f'unterminated quoted label: {chunk!r}')
        name = ''.join(label_parts)
        rest = chunk[j + 1 :].strip()
    else:
        name, sep, rest = chunk.partition(':')
        name = name.strip() or None
        if not sep:
            return name, None
        return name, _parse_length(rest)
    if not rest:
        return name, None
    if not rest.startswith(':'):
        raise ValueError(f'unexpected text after quoted label: {chunk!r}')
    return name, _parse_length(rest[1:])


def _parse_length(text: str) -> float | None:
    text = text.strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError as exc:
        raise ValueError(f'invalid branch length: {text!r}') from exc


def _parse_support(name: str | None) -> float | None:
    """Extract a numeric support value from an internal-node label."""
    if not name:
        return None
    candidate = name.rsplit('/', 1)[-1]  # ')lbl/95' IQ-TREE style
    try:
        return float(candidate)
    except ValueError:
        return None


def _parse_newick(text: str) -> TreeNode:
    text = _strip_comments(text).strip()
    end = _find_semicolon(text)
    if end is not None:
        text = text[:end]
    text = text.strip()
    if not text:
        raise ValueError('empty newick input')

    root = TreeNode()
    node = root
    ancestors: list[TreeNode] = []
    chunk_start = 0
    i = 0
    n = len(text)

    def apply_chunk(upto: int) -> None:
        name, length = _split_chunk(text[chunk_start:upto])
        if name is not None:
            node.name = name
            if node.children:
                node.support = _parse_support(name)
        if length is not None:
            node.length = length

    while i < n:
        ch = text[i]
        if ch == "'":  # skip over quoted section inside the current chunk
            j = i + 1
            while j < n:
                if text[j] == "'":
                    if j + 1 < n and text[j + 1] == "'":
                        j += 2
                        continue
                    break
                j += 1
            i = j + 1
            continue
        if ch == '(':
            child = TreeNode()
            node.children.append(child)
            ancestors.append(node)
            node = child
            chunk_start = i + 1
        elif ch == ',':
            apply_chunk(i)
            if not ancestors:
                raise ValueError('unbalanced parentheses in newick (stray comma)')
            node = TreeNode()
            ancestors[-1].children.append(node)
            chunk_start = i + 1
        elif ch == ')':
            apply_chunk(i)
            if not ancestors:
                raise ValueError('unbalanced parentheses in newick')
            node = ancestors.pop()
            chunk_start = i + 1
        i += 1
    apply_chunk(n)
    if ancestors:
        raise ValueError('unbalanced parentheses in newick (unclosed group)')
    return root


def _find_semicolon(text: str) -> int | None:
    """Index of the first ``;`` outside quotes, or None."""
    i, n = 0, len(text)
    while i < n:
        ch = text[i]
        if ch == "'":
            j = i + 1
            while j < n:
                if text[j] == "'":
                    if j + 1 < n and text[j + 1] == "'":
                        j += 2
                        continue
                    break
                j += 1
            i = j + 1
            continue
        if ch == ';':
            return i
        i += 1
    return None


# ----------------------------------------------------------------------
# Drawing
# ----------------------------------------------------------------------


def draw_tree(
    ax: 'matplotlib.axes.Axes',
    tree: Tree,
    leaf_pos: dict[str, float],
    *,
    orientation: str = 'left',
    cutoff: float | None = None,
    scalebar: bool = True,
    leaf_labels: bool = False,
    label_size: float | None = None,
    color: str = '0.2',
    lw: float = 1.0,
    gid_prefix: str = 'rd-tree',
) -> None:
    """Draw *tree* as a rectangular dendrogram on *ax*.

    The tree is drawn with the root at the left edge and tips pointing
    right (toward the plot the axis sits beside). The caller supplies
    the y position of every leaf via *leaf_pos*, so tips line up with
    plot rows of any height; the y limits of *ax* are left untouched.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to draw into; spines and ticks are hidden.
    tree : Tree
        Tree to draw.
    leaf_pos : dict of str to float
        Maps each leaf name to its y coordinate in *ax* data space.
    orientation : str
        Only ``'left'`` (tree left of the plot, tips pointing right) is
        supported.
    cutoff : float or None
        When given, draw a dashed vertical line at this distance from
        the tips (matching a clustering cutoff for ultrametric trees).
    scalebar : bool
        Draw a branch-length scale bar below the tree.
    leaf_labels : bool
        Write each leaf name at its tip, right-aligned against the axis
        edge nearest the plot.
    label_size : float or None
        Font size for leaf labels (default: matplotlib small).
    color : str
        Line colour.
    lw : float
        Line width.
    gid_prefix : str
        SVG gid prefix: segments get ``{gid_prefix}``, the cutoff line
        ``{gid_prefix}-cutoff``, the scale bar ``{gid_prefix}-scalebar``.

    Raises
    ------
    ValueError
        On unsupported *orientation* or a leaf missing from *leaf_pos*.
    """
    from matplotlib.collections import LineCollection
    from matplotlib.transforms import blended_transform_factory

    if orientation != 'left':
        raise ValueError(f"orientation must be 'left', got {orientation!r}")
    missing = [n for n in tree.leaf_names() if n not in leaf_pos]
    if missing:
        raise ValueError(f'leaf_pos missing positions for leaves: {missing}')

    default_len = 1.0 if not tree._has_lengths() else 0.0

    # Node x = cumulative distance from the root; leaf y from leaf_pos,
    # internal y = midpoint of the outermost children.
    depth: dict[int, float] = {id(tree.root): 0.0}
    stack = [tree.root]
    while stack:
        parent = stack.pop()
        for child in parent.children:
            length = child.length if child.length is not None else default_len
            depth[id(child)] = depth[id(parent)] + length
            stack.append(child)

    ypos: dict[int, float] = {}
    for node in tree._postorder():
        if node.is_leaf:
            ypos[id(node)] = leaf_pos[node.name or '']
        else:
            child_ys = [ypos[id(c)] for c in node.children]
            ypos[id(node)] = (min(child_ys) + max(child_ys)) / 2.0

    segments: list[list[tuple[float, float]]] = []
    for node in tree._postorder():
        if node.is_leaf:
            continue
        x = depth[id(node)]
        child_ys = [ypos[id(c)] for c in node.children]
        # Vertical connector spanning the children.
        segments.append([(x, min(child_ys)), (x, max(child_ys))])
        # Horizontal branch out to each child.
        for child in node.children:
            segments.append([(x, ypos[id(child)]), (depth[id(child)], ypos[id(child)])])
    lines = LineCollection(segments, colors=color, linewidths=lw, capstyle='projecting')
    lines.set_gid(gid_prefix)
    ax.add_collection(lines)

    max_depth = max(depth.values()) if depth else 1.0
    if max_depth <= 0:
        max_depth = 1.0
    ax.set_xlim(0, max_depth * 1.02)

    if cutoff is not None:
        line = ax.axvline(max_depth - cutoff, color='0.4', lw=0.9, ls='--', zorder=1)
        line.set_gid(f'{gid_prefix}-cutoff')

    if leaf_labels:
        label_tf = blended_transform_factory(ax.transAxes, ax.transData)
        for leaf in tree.leaves():
            name = leaf.name or ''
            ax.text(
                0.99,
                leaf_pos[name],
                name,
                transform=label_tf,
                ha='right',
                va='center',
                fontsize=label_size if label_size is not None else 'small',
                clip_on=False,
                gid=f'{gid_prefix}-label',
            )

    if scalebar:
        bar = _round_1sf(max_depth / 4.0)
        bar_tf = blended_transform_factory(ax.transData, ax.transAxes)
        (bar_line,) = ax.plot(
            [0, bar],
            [-0.015, -0.015],
            transform=bar_tf,
            color=color,
            lw=lw,
            clip_on=False,
            solid_capstyle='butt',
        )
        bar_line.set_gid(f'{gid_prefix}-scalebar')
        ax.text(
            bar / 2.0,
            -0.025,
            f'{bar:g}',
            transform=bar_tf,
            ha='center',
            va='top',
            fontsize='x-small',
            clip_on=False,
        )

    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def _round_1sf(value: float) -> float:
    """Round *value* to one significant figure (for the scale bar)."""
    if value <= 0:
        return 1.0
    from math import floor, log10

    exp = floor(log10(value))
    return round(value, -exp)
