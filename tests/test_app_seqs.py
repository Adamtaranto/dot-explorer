"""Tests for the app's lazy sequence providers (core/seqs.py)."""

import gzip
from pathlib import Path
import sys

import pytest

sys.path.insert(
    0, str(Path(__file__).resolve().parents[1] / 'python' / 'dot_explorer' / 'app')
)

from core.fasta import content_digest, parse_fasta_bytes  # noqa: E402
from core.seqs import (  # noqa: E402
    FaidxProvider,
    InMemoryProvider,
    SequenceProvider,
    file_digest,
    provider_from_fasta_input,
    provider_from_path,
)

FASTA = b'>ctg1 desc\nACGTACGTAC\nGGGG\n>ctg2\nTTTTCCCC\n'
EXPECTED = {'ctg1': 'ACGTACGTACGGGG', 'ctg2': 'TTTTCCCC'}


@pytest.fixture
def fasta_path(tmp_path):
    p = tmp_path / 'asm.fasta'
    p.write_bytes(FASTA)
    return p


@pytest.fixture
def gz_path(tmp_path):
    p = tmp_path / 'asm.fasta.gz'
    p.write_bytes(gzip.compress(FASTA))
    return p


def _providers(tmp_path):
    plain = tmp_path / 'p.fasta'
    plain.write_bytes(FASTA)
    return [
        InMemoryProvider(parse_fasta_bytes(FASTA)),
        FaidxProvider(plain, content_digest(FASTA)),
    ]


def test_providers_agree_on_names_and_lengths(tmp_path):
    for prov in _providers(tmp_path):
        assert prov.names == ['ctg1', 'ctg2']
        assert prov.lengths() == {'ctg1': 14, 'ctg2': 8}
        assert prov.total_length == 22
        assert prov.digest == content_digest(FASTA)


def test_providers_agree_on_sequences(tmp_path):
    for prov in _providers(tmp_path):
        assert dict(prov.iter_records()) == EXPECTED
        for name, seq in EXPECTED.items():
            lazy = prov.get_lazy(name)
            assert len(lazy) == len(seq)
            assert lazy[2:6] == seq[2:6]
            assert prov.get_slice(name, 0, 4) == seq[:4]
        assert prov.get_lazy('missing') is None


def test_providers_satisfy_protocol(tmp_path):
    for prov in _providers(tmp_path):
        assert isinstance(prov, SequenceProvider)
    # A bare FastaInput is not a provider (the app's isinstance gate).
    assert not isinstance(parse_fasta_bytes(FASTA), SequenceProvider)


def test_lazy_slices_are_plain_strings(tmp_path):
    """Handlers revcomp/concatenate slices, so they must be real str."""
    for prov in _providers(tmp_path):
        window = prov.get_lazy('ctg1')[0:5]
        assert type(window) is str


def test_provider_from_path_plain(fasta_path):
    prov = provider_from_path(fasta_path)
    assert isinstance(prov, FaidxProvider)
    assert dict(prov.iter_records()) == EXPECTED
    assert prov.digest == content_digest(FASTA)


def test_provider_from_path_gzip(gz_path):
    prov = provider_from_path(gz_path)
    assert dict(prov.iter_records()) == EXPECTED
    # Digest is over the raw (compressed) upload, matching the old contract.
    assert prov.digest == content_digest(gz_path.read_bytes())


def test_file_digest_matches_content_digest(fasta_path):
    assert file_digest(fasta_path) == content_digest(FASTA)


def test_provider_from_path_invalid_raises_valueerror(tmp_path):
    """Invalid input surfaces the pure-Python parser's messages."""
    bad = tmp_path / 'bad.fasta'
    bad.write_bytes(b'ACGT\nno header here\n')
    with pytest.raises(ValueError, match='header'):
        provider_from_path(bad)


def test_provider_from_path_duplicate_names(tmp_path):
    dup = tmp_path / 'dup.fasta'
    dup.write_bytes(b'>a\nACGT\n>a\nGGGG\n')
    with pytest.raises(ValueError, match='[Dd]uplicate'):
        provider_from_path(dup)


def test_provider_from_fasta_input(tmp_path):
    fi = parse_fasta_bytes(FASTA)
    prov = provider_from_fasta_input(fi, tmp_path)
    assert dict(prov.iter_records()) == EXPECTED
    assert prov.digest == fi.digest


def test_provider_from_fasta_input_unwritable_falls_back(tmp_path):
    fi = parse_fasta_bytes(FASTA)
    prov = provider_from_fasta_input(fi, tmp_path / 'does-not-exist')
    assert isinstance(prov, InMemoryProvider)
    assert dict(prov.iter_records()) == EXPECTED
