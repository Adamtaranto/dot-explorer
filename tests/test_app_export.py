"""Tests for core.export: selection FASTA and per-cluster zip bundles."""

import io
from pathlib import Path
import sys
import zipfile

APP_DIR = Path(__file__).resolve().parent.parent / 'python' / 'dot_explorer' / 'app'
sys.path.insert(0, str(APP_DIR))

from core.export import (  # noqa: E402
    cluster_fasta_zip,
    reordered_fasta_text,
    selected_regions_fasta,
)

SEQS = {'c1': 'ACGTACGTAA', 'c2': 'GGGGCCCCTT', 'c3': 'TTTTAAAACC'}


def _get(contig, start, end):
    seq = SEQS.get(contig)
    return None if seq is None else seq[start:end]


def _records(text):
    out = {}
    name = None
    for line in text.splitlines():
        if line.startswith('>'):
            name = line[1:]
            out[name] = ''
        elif name is not None:
            out[name] += line
    return out


def test_selected_regions_dedupe_and_revcomp():
    regions = [
        ('c1', 0, 4, '+', 'query'),
        ('c1', 0, 4, '+', 'target'),  # same coords: written once
        ('c1', 0, 4, '-', 'query'),  # other strand: a distinct record
        ('c2', 2, 6, '-', 'target'),
        ('nope', 0, 3, '+', 'query'),  # no sequence: skipped
        ('c3', 5, 5, '+', 'query'),  # empty: skipped
    ]
    recs = _records(selected_regions_fasta(regions, _get))
    assert list(recs) == [
        'c1:1-4(+) side=query',
        'c1:1-4(-) side=query',
        'c2:3-6(-) side=target',
    ]
    assert recs['c1:1-4(+) side=query'] == 'ACGT'
    assert recs['c1:1-4(-) side=query'] == 'ACGT'  # ACGT is its own revcomp
    assert recs['c2:3-6(-) side=target'] == 'GGCC'  # revcomp(GGCC)


def test_selected_regions_wraps_lines():
    text = selected_regions_fasta([('c1', 0, 10, '+', 'q')], _get, line_width=4)
    assert text.splitlines()[1:] == ['ACGT', 'ACGT', 'AA']


def test_cluster_zip_one_fasta_per_cluster_with_orientation():
    records = list(SEQS.items())
    assignments = {'c1': 'cluster_2', 'c2': 'cluster_10', 'c3': 'cluster_2'}
    data = cluster_fasta_zip(records, assignments, reverse={'c3'}, order=['c3', 'c1'])
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        names = zf.namelist()
        assert names == ['cluster_2.fasta', 'cluster_10.fasta']
        two = zf.read('cluster_2.fasta').decode()
        ten = zf.read('cluster_10.fasta').decode()
    # Display order inside the cluster, and the flipped contig marked + revcomp'd.
    assert two == reordered_fasta_text(
        [('c3', SEQS['c3']), ('c1', SEQS['c1'])], ['c3', 'c1'], {'c3'}
    )
    assert two.startswith('>c3 reverse_complement\nGGTTTTAAAA\n>c1\n')
    assert ten == '>c2\nGGGGCCCCTT\n'


def test_cluster_zip_unassigned_bucket_last():
    records = list(SEQS.items())
    data = cluster_fasta_zip(records, {'c1': 'cluster_1'}, reverse=set())
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        assert zf.namelist() == ['cluster_1.fasta', 'unassigned.fasta']
        assert 'c2' in zf.read('unassigned.fasta').decode()
        assert 'c3' in zf.read('unassigned.fasta').decode()
