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
    reoriented_paf_records,
    selected_regions_fasta,
)

from dot_explorer.paf_io import PafRecord  # noqa: E402

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


def _rec(line):
    return PafRecord.from_line(line)


def test_reoriented_paf_untouched_without_flips():
    rec = _rec('q\t100\t10\t30\t+\tt\t200\t40\t70\t20\t30\t60\tcg:Z:10M2I8M')
    out = reoriented_paf_records([rec], set())
    assert out == [rec]
    assert out[0] is rec


def test_reoriented_paf_query_flip_mirrors_coords_strand_and_cigar():
    rec = _rec('q\t100\t10\t30\t+\tt\t200\t40\t70\t20\t30\t60\tcg:Z:10M2I8M')
    (out,) = reoriented_paf_records([rec], {'q'})
    assert (out.query_start, out.query_end) == (70, 90)
    assert (out.target_start, out.target_end) == (40, 70)
    assert out.strand == '-'
    assert out.cigar == '8M2I10M'
    assert out.to_line().endswith('cg:Z:8M2I10M')
    # Lengths, counts and names survive; the input is untouched.
    assert (out.query_len, out.residue_matches, out.mapping_quality) == (100, 20, 60)
    assert rec.strand == '+' and rec.cigar == '10M2I8M'


def test_reoriented_paf_target_flip_and_double_flip():
    rec = _rec('q\t100\t10\t30\t-\tt\t200\t40\t70\t20\t30\t60\tcg:Z:10M2I8M')
    (t_only,) = reoriented_paf_records([rec], set(), {'t'})
    assert (t_only.target_start, t_only.target_end) == (130, 160)
    assert (t_only.query_start, t_only.query_end) == (10, 30)
    assert t_only.strand == '+'
    assert t_only.cigar == '8M2I10M'
    # Self mode: the same contig flipped on both axes keeps strand and CIGAR.
    both = _rec('c\t100\t10\t30\t-\tc\t100\t40\t70\t20\t30\t60\tcg:Z:10M2I8M')
    (b,) = reoriented_paf_records([both], {'c'}, {'c'})
    assert (b.query_start, b.query_end, b.target_start, b.target_end) == (
        70,
        90,
        30,
        60,
    )
    assert b.strand == '-'
    assert b.cigar == '10M2I8M'


def test_reoriented_paf_is_an_involution():
    rec = _rec('q\t100\t10\t30\t+\tt\t200\t40\t70\t20\t30\t60\tcg:Z:10M2I8M')
    twice = reoriented_paf_records(reoriented_paf_records([rec], {'q'}), {'q'})
    assert twice[0].to_line() == rec.to_line()
