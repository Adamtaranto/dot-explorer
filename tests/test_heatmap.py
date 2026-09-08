"""Tests for dot_explorer.heatmap (similarity heatmap rendering)."""

import io

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pytest

from dot_explorer import (
    ClusterResult,
    SimilarityMatrix,
    Tree,
    plot_similarity_heatmap,
)


@pytest.fixture
def sim():
    names = ['A', 'B', 'C', 'D']
    values = np.array(
        [
            [1.0, 0.9, 0.2, 0.1],
            [0.9, 1.0, 0.25, 0.15],
            [0.2, 0.25, 1.0, 0.85],
            [0.1, 0.15, 0.85, 1.0],
        ]
    )
    return SimilarityMatrix(names=names, values=values, metric='jaccard')


@pytest.fixture
def tree():
    return Tree.from_newick('((A:0.1,B:0.1):0.6,(C:0.15,D:0.15):0.55);')


@pytest.fixture
def clusters():
    return ClusterResult(
        assignments={
            'A': 'cluster_1',
            'B': 'cluster_1',
            'C': 'cluster_2',
            'D': 'cluster_2',
        },
        cutoff=0.8,
        mode='similarity',
        metric='jaccard',
    )


def _svg_of(fig) -> str:
    buf = io.StringIO()
    fig.savefig(buf, format='svg')
    return buf.getvalue()


class TestHeatmap:
    def test_basic_axes(self, sim):
        fig = plot_similarity_heatmap(sim)
        try:
            # Colorbar hangs off the heatmap box as an inset (child) axis.
            assert len(fig.axes) == 1
            assert len(fig.axes[0].child_axes) == 1
        finally:
            plt.close(fig)

    def test_no_colorbar(self, sim):
        fig = plot_similarity_heatmap(sim, colorbar=False)
        try:
            assert len(fig.axes) == 1
        finally:
            plt.close(fig)

    def test_tree_adds_axis_and_reorders(self, sim, tree):
        fig = plot_similarity_heatmap(sim, tree=tree)
        try:
            # Tree + colorbar are inset (child) axes of the heatmap box.
            assert len(fig.axes[0].child_axes) == 2
            svg = _svg_of(fig)
            assert 'de-heatmap' in svg
            assert 'de-tree' in svg
        finally:
            plt.close(fig)

    def test_cluster_outline_gids_in_svg(self, sim, tree, clusters):
        fig = plot_similarity_heatmap(sim, tree=tree, clusters=clusters, cutoff=0.5)
        try:
            svg = _svg_of(fig)
            assert 'de-hm-cluster-cluster_1' in svg
            assert 'de-hm-cluster-cluster_2' in svg
            assert 'de-tree-cutoff' in svg
        finally:
            plt.close(fig)

    def test_noncontiguous_cluster_outlines_runs(self, sim, caplog):
        scattered = ClusterResult(
            assignments={'A': 'c1', 'C': 'c1', 'B': 'c2', 'D': 'c2'},
            cutoff=0.5,
            mode='similarity',
            metric='jaccard',
        )
        with caplog.at_level('WARNING', logger='dot_explorer.heatmap'):
            fig = plot_similarity_heatmap(sim, clusters=scattered)
        try:
            assert 'not contiguous' in caplog.text
            svg = _svg_of(fig)
            assert svg.count('de-hm-cluster-c1') == 2
        finally:
            plt.close(fig)

    def test_explicit_order(self, sim):
        fig = plot_similarity_heatmap(sim, order=['D', 'C', 'B', 'A'], colorbar=False)
        try:
            labels = [t.get_text() for t in fig.axes[0].get_yticklabels()]
            assert labels == ['D', 'C', 'B', 'A']
        finally:
            plt.close(fig)

    def test_save_and_annotate(self, sim, tmp_path):
        out = tmp_path / 'heatmap.svg'
        fig = plot_similarity_heatmap(sim, annotate=True, output_path=str(out))
        try:
            assert out.exists()
            assert 'de-heatmap-scale' in out.read_text()
        finally:
            plt.close(fig)


class TestHeatmapValidation:
    def test_single_sequence_raises(self):
        one = SimilarityMatrix(names=['A'], values=np.ones((1, 1)), metric='jaccard')
        with pytest.raises(ValueError, match='at least 2'):
            plot_similarity_heatmap(one)

    def test_tree_and_order_conflict(self, sim, tree):
        with pytest.raises(ValueError, match='order'):
            plot_similarity_heatmap(sim, tree=tree, order=list('ABCD'))

    def test_tree_label_mismatch(self, sim):
        bad = Tree.from_newick('((A,B),(C,X));')
        with pytest.raises(ValueError, match='do not match'):
            plot_similarity_heatmap(sim, tree=bad)

    def test_html_output_rejected(self, sim, tmp_path):
        with pytest.raises(ValueError, match='HTML'):
            plot_similarity_heatmap(sim, output_path=str(tmp_path / 'x.html'))


class TestHeatmapGeometry:
    def test_axes_box_is_square(self, sim, tree):
        fig = plot_similarity_heatmap(sim, tree=tree)
        try:
            hm_ax = next(
                a
                for a in fig.axes
                if a.get_images()  # the imshow axis
            )
            pos = hm_ax.get_position()
            w_in = pos.width * fig.get_figwidth()
            h_in = pos.height * fig.get_figheight()
            assert w_in == pytest.approx(h_in, rel=0.05)
        finally:
            plt.close(fig)

    def test_names_shown_beside_tree(self, sim, tree):
        # With a tree the names are the heatmap's y tick labels (in the
        # spacer between tree and cells), not tree-axis text.
        fig = plot_similarity_heatmap(sim, tree=tree)
        try:
            hm_ax = next(a for a in fig.axes if a.get_images())
            labels = [t.get_text() for t in hm_ax.get_yticklabels()]
            assert labels == tree.leaf_names()
        finally:
            plt.close(fig)

    def test_ani_annotate_shows_ci(self):
        names = ['A', 'B']
        values = np.array([[1.0, 0.98], [0.98, 1.0]])
        ci = np.array([[1.0, 0.96], [0.96, 1.0]])
        ani = SimilarityMatrix(
            names=names,
            values=values,
            metric='ani',
            ci_low=ci,
            ci_high=np.array([[1.0, 0.99], [0.99, 1.0]]),
        )
        fig = plot_similarity_heatmap(ani, annotate=True, colorbar=False)
        try:
            texts = [t.get_text() for t in fig.axes[0].texts]
            assert any('0.96' in t and '0.99' in t for t in texts)
        finally:
            plt.close(fig)
