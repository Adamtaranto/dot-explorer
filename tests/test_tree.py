"""Tests for dot_explorer.tree: newick parsing, linkage conversion, drawing."""

from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pytest

from dot_explorer import Tree, draw_tree

DATA = Path(__file__).parent / 'data'


class TestNewickParsing:
    def test_simple_topology_and_lengths(self):
        tree = Tree.from_newick('((A:0.1,B:0.2):0.05,C:0.3);')
        assert tree.leaf_names() == ['A', 'B', 'C']
        leaves = {leaf.name: leaf for leaf in tree.leaves()}
        assert leaves['A'].length == pytest.approx(0.1)
        assert leaves['C'].length == pytest.approx(0.3)

    def test_no_branch_lengths(self):
        tree = Tree.from_newick('((A,B),(C,D));')
        assert tree.leaf_names() == ['A', 'B', 'C', 'D']
        # Unit edges when no lengths anywhere: leaves sit two levels deep.
        assert tree.max_depth() == pytest.approx(2.0)

    def test_iqtree_support_plain(self):
        tree = Tree.from_newick('((A:0.1,B:0.2)95:0.05,C:0.3);')
        internal = tree.root.children[0]
        assert internal.support == pytest.approx(95.0)
        assert internal.name == '95'

    def test_iqtree_support_slash_style(self):
        tree = Tree.from_newick('((A:0.1,B:0.2)87/0.99:0.05,C:0.3);')
        internal = tree.root.children[0]
        assert internal.support == pytest.approx(0.99)
        assert internal.name == '87/0.99'

    def test_quoted_labels_with_escapes(self):
        tree = Tree.from_newick("(('contig 1':0.1,'it''s':0.2):0.05,C:0.3);")
        assert tree.leaf_names() == ['contig 1', "it's", 'C']

    def test_comments_stripped(self):
        tree = Tree.from_newick('((A[comment]:0.1,B:0.2):0.05,C:0.3);[trailing]')
        assert tree.leaf_names() == ['A', 'B', 'C']

    def test_text_after_semicolon_ignored(self):
        tree = Tree.from_newick('(A:1,B:2);\nnot a tree\n')
        assert tree.leaf_names() == ['A', 'B']

    def test_read_treefile_fixture(self):
        tree = Tree.read(DATA / 'example.treefile')
        assert tree.leaf_names() == ['chr1', 'chr2', 'chr3', 'chr4']
        assert tree.max_depth() == pytest.approx(0.0234567 + 0.02, abs=1e-6)

    def test_deep_ladder_tree_no_recursion_error(self):
        n = 5000
        text = '(' * n + 'L0'
        for i in range(1, n + 1):
            text += f',L{i})'
        tree = Tree.from_newick(text + ';')
        assert len(tree.leaf_names()) == n + 1

    def test_unbalanced_raises(self):
        with pytest.raises(ValueError, match='unbalanced'):
            Tree.from_newick('((A,B);')
        with pytest.raises(ValueError, match='unbalanced'):
            Tree.from_newick('(A,B));')

    def test_unterminated_quote_raises(self):
        with pytest.raises(ValueError, match='unterminated'):
            Tree.from_newick("('A,B);")

    def test_empty_raises(self):
        with pytest.raises(ValueError):
            Tree.from_newick('   ;')

    def test_duplicate_leaves_raise(self):
        with pytest.raises(ValueError, match='duplicate'):
            Tree.from_newick('((A,A),B);')


class TestValidateLabels:
    def test_matching_passes(self):
        tree = Tree.from_newick('((A,B),C);')
        tree.validate_labels(['C', 'A', 'B'])  # order-insensitive

    def test_mismatch_lists_both_directions(self):
        tree = Tree.from_newick('((A,B),C);')
        with pytest.raises(ValueError) as exc:
            tree.validate_labels(['A', 'B', 'D'])
        message = str(exc.value)
        assert "tree tips with no matching sequence: ['C']" in message
        assert "sequences with no tip in the tree: ['D']" in message


class TestFromLinkage:
    def test_leaf_order_matches_scipy_dendrogram(self):
        np = pytest.importorskip('numpy')
        hierarchy = pytest.importorskip('scipy.cluster.hierarchy')
        rng = np.random.default_rng(42)
        points = rng.random((7, 3))
        Z = hierarchy.linkage(points, method='average')
        labels = [f'seq{i}' for i in range(7)]
        tree = Tree.from_linkage(Z, labels)
        ivl = hierarchy.dendrogram(Z, no_plot=True, labels=labels)['ivl']
        assert tree.leaf_names() == ivl

    def test_ultrametric_depth_equals_root_height(self):
        np = pytest.importorskip('numpy')
        hierarchy = pytest.importorskip('scipy.cluster.hierarchy')
        rng = np.random.default_rng(7)
        Z = hierarchy.linkage(rng.random((5, 4)), method='complete')
        tree = Tree.from_linkage(Z, list('abcde'))
        assert tree.max_depth() == pytest.approx(Z[-1, 2])

    def test_row_count_mismatch_raises(self):
        np = pytest.importorskip('numpy')
        Z = np.zeros((2, 4))
        with pytest.raises(ValueError, match='expected'):
            Tree.from_linkage(Z, ['a', 'b'])


class TestDrawTree:
    def _tree(self):
        return Tree.from_newick('((A:0.1,B:0.2):0.05,C:0.3);')

    def test_draws_segments_with_gids(self):
        fig, ax = plt.subplots()
        try:
            tree = self._tree()
            pos = {'A': 2.5, 'B': 1.5, 'C': 0.5}
            ax.set_ylim(0, 3)
            draw_tree(ax, tree, pos, cutoff=0.1, leaf_labels=True)
            gids = {a.get_gid() for a in ax.get_children() if a.get_gid()}
            assert 'de-tree' in gids
            assert 'de-tree-cutoff' in gids
            assert 'de-tree-scalebar' in gids
            assert 'de-tree-label' in gids
        finally:
            plt.close(fig)

    def test_missing_leaf_pos_raises(self):
        fig, ax = plt.subplots()
        try:
            with pytest.raises(ValueError, match='leaf_pos missing'):
                draw_tree(ax, self._tree(), {'A': 0.0, 'B': 1.0})
        finally:
            plt.close(fig)

    def test_bad_orientation_raises(self):
        fig, ax = plt.subplots()
        try:
            with pytest.raises(ValueError, match='orientation'):
                draw_tree(
                    ax,
                    self._tree(),
                    {'A': 0, 'B': 1, 'C': 2},
                    orientation='top',
                )
        finally:
            plt.close(fig)


class TestBranchLengthFidelity:
    def test_segment_extents_proportional_to_branch_lengths(self):
        """User treefiles with branch lengths keep their relative depths."""
        from matplotlib.collections import LineCollection

        # A at depth 0.1, B at 0.3; the (A,B) clade node at depth... root
        # children: clade (length 0.05) and C (0.4).
        tree = Tree.from_newick('((A:0.05,B:0.25)90:0.05,C:0.4);')
        fig, ax = plt.subplots()
        try:
            ax.set_ylim(0, 3)
            draw_tree(ax, tree, {'A': 2.5, 'B': 1.5, 'C': 0.5}, scalebar=False)
            lines = next(a for a in ax.collections if isinstance(a, LineCollection))
            xs = {x for seg in lines.get_segments() for x, _y in seg}
            # Cumulative depths that must appear as segment endpoints:
            # root 0, clade node 0.05, A tip 0.10, B tip 0.30, C tip 0.40.
            for depth in (0.0, 0.05, 0.10, 0.30, 0.40):
                assert any(abs(x - depth) < 1e-9 for x in xs), depth
        finally:
            plt.close(fig)
