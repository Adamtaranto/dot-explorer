"""CodSpeed benchmarks for the browser app's FASTA byte parser.

Run with ``pytest python/benchmarks --codspeed``.  The app parses uploads in
pure Python under Pyodide, so parser throughput directly bounds how quickly a
run can start on real multi-megabyte assemblies.
"""

from __future__ import annotations

import gzip
import sys

from _synth import random_dna

from dot_explorer.launch import APP_DIR

# The parser ships inside the wheel as package data rather than as an
# importable subpackage: app.py and its siblings use flat `core.*` imports
# that resolve once the app directory is on sys.path (what `shiny run` and
# `shinylive export` both do).  Take the location from the package so this
# keeps working whether the benchmarks run against a checkout or an install.
sys.path.insert(0, str(APP_DIR))

from core.fasta import parse_fasta_bytes  # noqa: E402

# ~2 Mb across a handful of contigs, wrapped at 70 columns like real
# assembler output.  Kept modest so CodSpeed's instruction-count simulation
# stays quick while still exercising the per-record block path.
_N_CONTIGS = 5
_CONTIG_LEN = 400_000
_WRAP = 70


def _fasta_bytes() -> bytes:
    """Build a deterministic multi-contig, line-wrapped FASTA payload."""
    chunks: list[str] = []
    for i in range(_N_CONTIGS):
        seq = random_dna(_CONTIG_LEN, seed=100 + i)
        chunks.append(f'>contig_{i} length={_CONTIG_LEN}\n')
        chunks.extend(seq[j : j + _WRAP] + '\n' for j in range(0, len(seq), _WRAP))
    return ''.join(chunks).encode()


_DATA = _fasta_bytes()
_DATA_GZ = gzip.compress(_DATA, compresslevel=6)


def test_bench_parse_fasta_bytes(benchmark):
    """Parse a multi-contig line-wrapped FASTA from raw bytes."""
    result = benchmark(parse_fasta_bytes, _DATA)
    assert len(result.records) == _N_CONTIGS


def test_bench_parse_fasta_bytes_gzip(benchmark):
    """Parse the same payload gzip-compressed (decompress + parse)."""
    result = benchmark(parse_fasta_bytes, _DATA_GZ)
    assert result.total_length == _N_CONTIGS * _CONTIG_LEN
