"""Tests for dot_explorer.similarity (sourmash sketching + clustering).

Skipped wholesale when the `cluster` extra (sourmash/scipy) is absent.
"""

import random

import numpy as np
import pytest

pytest.importorskip('sourmash')

from dot_explorer import (
    ClusterResult,
    SequenceIndex,
    SimilarityMatrix,
    SketchParams,
    assign_clusters,
    assign_clusters_dual,
    compute_sketches,
    linkage_from_similarity,
    pairwise_similarity,
)

PARAMS = SketchParams(ksize=21, scaled=10, track_abundance=True)


def _random_seq(rng: random.Random, length: int) -> str:
    return ''.join(rng.choice('ACGT') for _ in range(length))


@pytest.fixture(scope='module')
def indexed_family():
    """Index of related and unrelated sequences.

    seqA and seqB share their first 15 kb; seqC is seqA's first 8 kb
    (nested, for containment asymmetry); seqD is unrelated.
    """
    rng = random.Random(42)
    a = _random_seq(rng, 20_000)
    b = a[:15_000] + _random_seq(rng, 5_000)
    c = a[:8_000]
    d = _random_seq(rng, 20_000)
    idx = SequenceIndex(k=21)
    for name, seq in [('seqA', a), ('seqB', b), ('seqC', c), ('seqD', d)]:
        idx.add_sequence(name, seq)
    return idx


NAMES = ['seqA', 'seqB', 'seqC', 'seqD']


@pytest.fixture(scope='module')
def sketches(indexed_family):
    return compute_sketches(indexed_family, NAMES, params=PARAMS)


class TestComputeSketches:
    def test_explicit_names_keep_order(self, sketches):
        assert list(sketches) == NAMES

    def test_default_names_cover_index(self, indexed_family):
        got = compute_sketches(indexed_family, params=PARAMS)
        assert sorted(got) == NAMES

    def test_deterministic(self, indexed_family, sketches):
        again = compute_sketches(indexed_family, params=PARAMS)
        assert sketches['seqA'].hashes == again['seqA'].hashes

    def test_subset_and_params(self, indexed_family):
        got = compute_sketches(
            indexed_family,
            ['seqA'],
            params=SketchParams(ksize=31, scaled=100, track_abundance=False),
        )
        assert list(got) == ['seqA']
        assert got['seqA'].ksize == 31
        assert not got['seqA'].track_abundance


class TestPairwiseSimilarity:
    def test_jaccard_symmetric_unit_diagonal(self, sketches):
        sim = pairwise_similarity(sketches, metric='jaccard')
        assert np.allclose(sim.values, sim.values.T)
        assert np.allclose(np.diag(sim.values), 1.0)
        assert sim[('seqA', 'seqB')] > 0.5
        assert sim[('seqA', 'seqD')] < 0.05

    def test_angular_differs_from_jaccard(self, sketches):
        jac = pairwise_similarity(sketches, metric='jaccard')
        ang = pairwise_similarity(sketches, metric='angular')
        assert ang.values[0, 1] != pytest.approx(jac.values[0, 1])

    def test_angular_without_abundance_raises(self, indexed_family):
        flat = compute_sketches(
            indexed_family,
            params=SketchParams(scaled=10, track_abundance=False),
        )
        with pytest.raises(ValueError, match='abundance'):
            pairwise_similarity(flat, metric='angular')

    def test_angular_ignore_abundance_raises(self, sketches):
        with pytest.raises(ValueError, match='abundance'):
            pairwise_similarity(sketches, metric='angular', ignore_abundance=True)

    def test_ani_with_confidence_intervals(self, sketches):
        sim = pairwise_similarity(sketches, metric='ani')
        i, j = 0, 1  # seqA vs seqB
        assert sim.ci_low is not None and sim.ci_high is not None
        assert sim.ci_low[i, j] <= sim.values[i, j] <= sim.ci_high[i, j]
        assert sim.values[i, j] > 0.9  # near-identical shared prefix

    def test_containment_asymmetric_for_nested(self, sketches):
        sim = pairwise_similarity(sketches, metric='containment')
        # seqC is a subsequence of seqA: C-in-A ~1, A-in-C ~0.4.
        c_in_a = sim[('seqC', 'seqA')]
        a_in_c = sim[('seqA', 'seqC')]
        assert c_in_a > 0.95
        assert a_in_c < 0.6
        assert c_in_a != pytest.approx(a_in_c)

    def test_max_containment_symmetric(self, sketches):
        sim = pairwise_similarity(sketches, metric='max_containment')
        assert np.allclose(sim.values, sim.values.T)
        assert sim[('seqC', 'seqA')] > 0.95

    def test_unknown_metric_raises(self, sketches):
        with pytest.raises(ValueError, match='unknown metric'):
            pairwise_similarity(sketches, metric='euclidean')

    def test_single_sketch_raises(self, sketches):
        with pytest.raises(ValueError, match='at least 2'):
            pairwise_similarity({'seqA': sketches['seqA']})


class TestSimilarityMatrix:
    def test_reorder_permutes_values(self, sketches):
        sim = pairwise_similarity(sketches, metric='jaccard')
        new_order = ['seqD', 'seqB', 'seqA', 'seqC']
        re = sim.reorder(new_order)
        assert re.names == new_order
        assert re[('seqA', 'seqB')] == pytest.approx(sim[('seqA', 'seqB')])

    def test_reorder_bad_names_raises(self, sketches):
        sim = pairwise_similarity(sketches, metric='jaccard')
        with pytest.raises(ValueError, match='permutation'):
            sim.reorder(['seqA', 'seqB'])

    def test_csv_round_trip(self, sketches, tmp_path):
        sim = pairwise_similarity(sketches, metric='jaccard')
        out = tmp_path / 'sim.csv'
        sim.to_csv(out)
        lines = out.read_text().strip().splitlines()
        assert lines[0] == ',seqA,seqB,seqC,seqD'
        cells = lines[1].split(',')
        assert cells[0] == 'seqA'
        assert float(cells[1]) == pytest.approx(1.0)
        assert float(cells[2]) == pytest.approx(sim[('seqA', 'seqB')], rel=1e-4)


class TestClustering:
    def test_linkage_shape(self, sketches):
        pytest.importorskip('scipy')
        sim = pairwise_similarity(sketches, metric='jaccard')
        Z = linkage_from_similarity(sim)
        assert Z.shape == (3, 4)

    def test_cutoff_extremes(self, sketches):
        pytest.importorskip('scipy')
        sim = pairwise_similarity(sketches, metric='jaccard')
        everyone = assign_clusters(sim, cutoff=0.0)
        assert len(set(everyone.assignments.values())) == 1
        singletons = assign_clusters(sim, cutoff=1.0)
        assert len(set(singletons.assignments.values())) == 4

    def test_moderate_cutoff_groups_family(self, sketches):
        pytest.importorskip('scipy')
        sim = pairwise_similarity(sketches, metric='jaccard')
        result = assign_clusters(sim, cutoff=0.3)
        assert result.assignments['seqA'] == result.assignments['seqB']
        assert result.assignments['seqD'] != result.assignments['seqA']
        assert result.mode == 'similarity'

    def test_clusters_property_and_csv(self, sketches, tmp_path):
        pytest.importorskip('scipy')
        sim = pairwise_similarity(sketches, metric='jaccard')
        result = assign_clusters(sim, cutoff=0.3)
        members = result.clusters[result.assignments['seqA']]
        assert 'seqA' in members and 'seqB' in members
        out = tmp_path / 'clusters.csv'
        result.to_csv(out)
        lines = out.read_text().strip().splitlines()
        assert lines[0] == 'contig,cluster'
        assert len(lines) == 5


class TestDualClustering:
    def _matrices(self):
        names = ['a', 'b', 'c']
        ident = SimilarityMatrix(
            names=names,
            values=np.array(
                [
                    [1.0, 0.95, 0.95],
                    [0.95, 1.0, 0.5],
                    [0.95, 0.5, 1.0],
                ]
            ),
            metric='ani',
        )
        # a<->b reciprocal coverage; a->c passes but c->a fails.
        cov = SimilarityMatrix(
            names=names,
            values=np.array(
                [
                    [1.0, 0.9, 0.9],
                    [0.9, 1.0, 0.1],
                    [0.3, 0.1, 1.0],
                ]
            ),
            metric='containment',
        )
        return ident, cov

    def test_reciprocal_excludes_one_way_pair(self):
        ident, cov = self._matrices()
        result = assign_clusters_dual(ident, cov, reciprocal=True)
        assert result.assignments['a'] == result.assignments['b']
        assert result.assignments['c'] != result.assignments['a']
        assert result.reciprocal is True

    def test_non_reciprocal_includes_one_way_pair(self):
        ident, cov = self._matrices()
        result = assign_clusters_dual(ident, cov, reciprocal=False)
        assert result.assignments['c'] == result.assignments['a']

    def test_identity_threshold_applies(self):
        ident, cov = self._matrices()
        result = assign_clusters_dual(
            ident, cov, identity_cutoff=0.99, reciprocal=False
        )
        assert len(set(result.assignments.values())) == 3

    def test_result_type(self):
        ident, cov = self._matrices()
        result = assign_clusters_dual(ident, cov)
        assert isinstance(result, ClusterResult)
        assert result.mode == 'identity_coverage'
