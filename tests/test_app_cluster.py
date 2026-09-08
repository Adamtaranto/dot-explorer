"""Tests for the app's clustering helpers (app/core/cluster.py)."""

from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'app'))

from core.cluster import (  # noqa: E402
    ProviderIndex,
    cluster_deps_missing,
    cluster_table_rows,
    tree_layout_order,
)
from core.fasta import parse_fasta_bytes  # noqa: E402
from core.seqs import InMemoryProvider  # noqa: E402

from rusty_dot import ClusterResult, SimilarityMatrix  # noqa: E402


@pytest.fixture
def provider():
    fasta = parse_fasta_bytes(b'>c1\nACGTACGTAC\n>c2\nGGGGCCCC\n')
    return InMemoryProvider(fasta)


class TestProviderIndex:
    def test_names_and_sequences(self, provider):
        adapter = ProviderIndex(provider)
        assert adapter.sequence_names() == ['c1', 'c2']
        assert adapter.get_sequence('c1') == 'ACGTACGTAC'
        assert adapter.get_sequence('c2') == 'GGGGCCCC'

    def test_sketching_through_adapter(self, provider):
        pytest.importorskip('sourmash')
        from rusty_dot import SketchParams, compute_sketches

        sketches = compute_sketches(
            ProviderIndex(provider),
            params=SketchParams(ksize=5, scaled=1),
        )
        assert set(sketches) == {'c1', 'c2'}
        assert len(sketches['c1'].hashes) > 0


class TestClusterTableRows:
    def _clusters(self):
        return ClusterResult(
            assignments={
                'a': 'cluster_1',
                'b': 'cluster_1',
                'c': 'cluster_2',
            },
            cutoff=0.8,
            mode='similarity',
            metric='jaccard',
        )

    def test_rows_grouped_by_cluster(self):
        rows = cluster_table_rows(self._clusters(), {'a': 100, 'b': 50, 'c': 25})
        assert [(r['cluster'], r['contig']) for r in rows] == [
            ('cluster_1', 'a'),
            ('cluster_1', 'b'),
            ('cluster_2', 'c'),
        ]
        assert rows[0]['length'] == 100
        assert rows[0]['members'] == 2
        assert rows[2]['members'] == 1

    def test_mean_similarity(self):
        sim = SimilarityMatrix(
            names=['a', 'b', 'c'],
            values=np.array([[1.0, 0.9, 0.2], [0.9, 1.0, 0.3], [0.2, 0.3, 1.0]]),
            metric='jaccard',
        )
        rows = cluster_table_rows(self._clusters(), {}, sim=sim)
        by_contig = {r['contig']: r for r in rows}
        assert by_contig['a']['mean_sim'] == pytest.approx(0.9)
        assert by_contig['c']['mean_sim'] is None  # singleton

    def test_missing_length_defaults_to_zero(self):
        rows = cluster_table_rows(self._clusters(), {})
        assert all(r['length'] == 0 for r in rows)


class TestTreeLayoutOrder:
    def test_full_cover_reorders(self):
        assert tree_layout_order(['c', 'a', 'b'], ['a', 'b', 'c']) == [
            'c',
            'a',
            'b',
        ]

    def test_filtered_subset_keeps_tree_order(self):
        # min-length filter dropped 'b'; the tree still orders the rest.
        assert tree_layout_order(['c', 'a', 'b'], ['a', 'c']) == ['c', 'a']

    def test_unknown_plotted_name_returns_none(self):
        assert tree_layout_order(['a', 'b'], ['a', 'b', 'x']) is None


def test_cluster_deps_missing_reports_installed_state():
    missing = cluster_deps_missing()
    # In the test env both are installed; the contract is the type/shape.
    assert isinstance(missing, list)
    assert all(m in ('sourmash', 'scipy') for m in missing)
