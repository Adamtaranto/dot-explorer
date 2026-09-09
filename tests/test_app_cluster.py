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

from dot_explorer import ClusterResult, SimilarityMatrix  # noqa: E402


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
        from dot_explorer import SketchParams, compute_sketches

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


class TestMergeIntervals:
    def test_overlapping_and_adjacent(self):
        from core.cluster import merge_intervals

        assert merge_intervals([(0, 10), (5, 15), (15, 20), (30, 40)]) == [
            (0, 20),
            (30, 40),
        ]

    def test_reversed_and_empty_intervals(self):
        from core.cluster import merge_intervals

        assert merge_intervals([(10, 0), (5, 5)]) == [(0, 10)]
        assert merge_intervals([]) == []


class TestAlignmentCoverageMatrix:
    def _record(self, q, qs, qe, t, ts, te):
        from dot_explorer import PafRecord

        return PafRecord(
            query_name=q,
            query_len=100,
            query_start=qs,
            query_end=qe,
            strand='+',
            target_name=t,
            target_len=100,
            target_start=ts,
            target_end=te,
            residue_matches=qe - qs,
            alignment_block_len=qe - qs,
            mapping_quality=60,
        )

    def test_asymmetric_nested_fragment(self):
        from core.cluster import alignment_coverage_matrix

        # frag (100 bp) fully aligns into a 400 bp region of chrom.
        records = [self._record('frag', 0, 100, 'chrom', 100, 200)]
        sim = alignment_coverage_matrix(
            records, ['frag', 'chrom'], {'frag': 100, 'chrom': 400}
        )
        assert sim.metric == 'aln_coverage'
        assert sim[('frag', 'chrom')] == pytest.approx(1.0)
        assert sim[('chrom', 'frag')] == pytest.approx(0.25)

    def test_overlapping_blocks_not_double_counted(self):
        from core.cluster import alignment_coverage_matrix

        records = [
            self._record('a', 0, 60, 'b', 0, 60),
            self._record('a', 40, 100, 'b', 40, 100),
        ]
        sim = alignment_coverage_matrix(records, ['a', 'b'], {'a': 100, 'b': 200})
        assert sim[('a', 'b')] == pytest.approx(1.0)
        assert sim[('b', 'a')] == pytest.approx(0.5)

    def test_missing_pair_and_normalize(self):
        from core.cluster import alignment_coverage_matrix

        records = [self._record('query:a', 0, 50, 'target:b', 0, 50)]
        sim = alignment_coverage_matrix(
            records,
            ['a', 'b', 'c'],
            {'a': 100, 'b': 100, 'c': 100},
            normalize=lambda n: n.split(':', 1)[1] if ':' in n else n,
        )
        assert sim[('a', 'b')] == pytest.approx(0.5)
        assert sim[('a', 'c')] == 0.0
        assert sim[('c', 'a')] == 0.0

    def test_secondary_alignments_excluded(self):
        from core.cluster import alignment_coverage_matrix

        primary = self._record('a', 0, 50, 'b', 0, 50)
        primary.tags['tp'] = 'P'
        secondary = self._record('a', 50, 100, 'b', 50, 100)
        secondary.tags['tp'] = 'S'
        sim = alignment_coverage_matrix(
            [primary, secondary], ['a', 'b'], {'a': 100, 'b': 100}
        )
        assert sim[('a', 'b')] == pytest.approx(0.5)
        assert sim[('b', 'a')] == pytest.approx(0.5)


class TestCleanMinimap2:
    def test_truth_table(self):
        from core.cluster import is_clean_minimap2

        assert is_clean_minimap2('minimap2', {'P': False})
        assert is_clean_minimap2('minimap2', {})
        assert is_clean_minimap2('minimap2', None)
        assert not is_clean_minimap2('minimap2', {'P': True})
        assert not is_clean_minimap2('nucmer', {})
        assert not is_clean_minimap2('kmer', {})
        assert not is_clean_minimap2('paf_upload', {})
        assert not is_clean_minimap2(None, None)

    def test_coverage_params_are_clean_and_deterministic(self):
        from core.align import build_tool_args
        from core.cluster import coverage_align_params, is_clean_minimap2

        params = coverage_align_params()
        assert params == coverage_align_params()
        assert params['P'] is False
        assert is_clean_minimap2('minimap2', params)
        args = build_tool_args('minimap2', params)
        assert '-P' not in args
