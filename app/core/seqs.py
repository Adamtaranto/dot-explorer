"""Lazy sequence access for uploaded assemblies.

Uploads used to be parsed into :class:`~core.fasta.FastaInput`, holding every
contig as a resident Python string for the whole session.  This module wraps
sequence access behind a small provider interface so the app can instead
index the uploaded file with :mod:`pyfaidx` and fetch only the windows it
needs (a 1,000 bp preview, a copy-request region, one contig at a time while
building the k-mer index).

Two implementations:

- :class:`FaidxProvider` — a ``pyfaidx.Fasta`` over the uploaded file's
  on-disk path (under Pyodide the "disk" is the RAM-backed MEMFS, which
  still saves the decoded-string copy of every contig).
- :class:`InMemoryProvider` — a thin adapter over a parsed ``FastaInput``,
  used as the fallback when ``pyfaidx`` is unavailable or cannot index the
  file (the pure-Python parser then also produces the app's user-facing
  validation errors).

Sequence values returned by :meth:`SequenceProvider.get_lazy` are "lazy
sequences": objects supporting ``len()`` and ``[start:stop]`` slicing that
returns a plain ``str``.  A Python string satisfies this, and so does a
``pyfaidx`` record with ``as_raw=True`` — callers slice first and only the
sliced window is ever materialised.
"""

from __future__ import annotations

import gzip
import hashlib
import logging
from pathlib import Path
import shutil
from typing import Iterator, Protocol, runtime_checkable

from .fasta import FastaInput, parse_fasta_bytes

logger = logging.getLogger(__name__)

_GZIP_MAGIC = b'\x1f\x8b'


class LazySeq(Protocol):
    """A sequence sliceable to ``str`` without materialising the whole."""

    def __len__(self) -> int:
        """Return the sequence length."""
        ...

    def __getitem__(self, item: slice) -> str:
        """Return the bases selected by *item* as a plain string."""
        ...


@runtime_checkable
class SequenceProvider(Protocol):
    """Read access to an uploaded assembly's sequences.

    ``isinstance`` checks against this protocol replace the app's previous
    ``isinstance(x, FastaInput)`` checks.
    """

    @property
    def digest(self) -> str:
        """Stable content digest of the raw uploaded bytes (cache key)."""
        ...

    @property
    def names(self) -> list[str]:
        """Contig names in file order."""
        ...

    @property
    def total_length(self) -> int:
        """Total number of residues across all contigs."""
        ...

    def lengths(self) -> dict[str, int]:
        """Return ``{contig name: length}`` without reading sequences."""
        ...

    def get_lazy(self, name: str) -> LazySeq | None:
        """Return a lazily sliceable sequence for *name*, or ``None``."""
        ...

    def get_slice(self, name: str, start: int, end: int) -> str:
        """Return ``sequence[start:end]`` for contig *name*."""
        ...

    def iter_records(self) -> Iterator[tuple[str, str]]:
        """Yield ``(name, sequence)`` pairs one contig at a time."""
        ...


class InMemoryProvider:
    """Sequence provider over a parsed :class:`~core.fasta.FastaInput`.

    Parameters
    ----------
    fasta : FastaInput
        The parsed upload to serve sequences from.
    """

    def __init__(self, fasta: FastaInput) -> None:
        self._fasta = fasta
        self._by_name = dict(fasta.records)

    @property
    def digest(self) -> str:
        """Stable content digest of the raw uploaded bytes.

        Returns
        -------
        str
            First 16 hex characters of the SHA-256 digest.
        """
        return self._fasta.digest

    @property
    def names(self) -> list[str]:
        """Contig names in file order.

        Returns
        -------
        list[str]
            Names in the order they appear in the file.
        """
        return self._fasta.names

    @property
    def total_length(self) -> int:
        """Total number of residues across all contigs.

        Returns
        -------
        int
            Sum of sequence lengths.
        """
        return self._fasta.total_length

    def lengths(self) -> dict[str, int]:
        """Return per-contig sequence lengths.

        Returns
        -------
        dict[str, int]
            ``{contig name: length}``.
        """
        return {name: len(seq) for name, seq in self._fasta.records}

    def get_lazy(self, name: str) -> LazySeq | None:
        """Return the sequence string for *name*, or ``None`` if unknown.

        Parameters
        ----------
        name : str
            Contig name.

        Returns
        -------
        str or None
            The full sequence (strings satisfy the lazy-slice protocol).
        """
        return self._by_name.get(name)

    def get_slice(self, name: str, start: int, end: int) -> str:
        """Return ``sequence[start:end]`` for contig *name*.

        Parameters
        ----------
        name : str
            Contig name (must exist).
        start : int
            0-based inclusive start.
        end : int
            0-based exclusive end.

        Returns
        -------
        str
            The requested window.
        """
        return self._by_name[name][start:end]

    def iter_records(self) -> Iterator[tuple[str, str]]:
        """Yield ``(name, sequence)`` pairs in file order.

        Yields
        ------
        tuple[str, str]
            One parsed record per contig.
        """
        yield from self._fasta.records


class FaidxProvider:
    """Sequence provider backed by a ``pyfaidx`` index over an on-disk file.

    Sequences are fetched from the file per request; nothing but the ``.fai``
    offsets (built on first open) is held in memory.

    Parameters
    ----------
    path : pathlib.Path
        Plain (uncompressed) FASTA file to index.  A ``.fai`` is written
        next to it.
    digest : str
        Content digest of the raw upload, preserved verbatim so cache keys
        match the previous ``FastaInput.digest`` contract.

    Raises
    ------
    Exception
        Whatever ``pyfaidx`` raises for unindexable input; callers are
        expected to fall back to :class:`InMemoryProvider` (see
        :func:`provider_from_path`).
    """

    def __init__(self, path: Path, digest: str) -> None:
        import pyfaidx

        self._path = path
        self._digest = digest
        # as_raw=True makes record slices return plain str; keep the file's
        # own case (the report shows sequences as uploaded).
        self._fasta = pyfaidx.Fasta(str(path), as_raw=True, sequence_always_upper=False)
        if not len(self._fasta.keys()):
            raise ValueError('No FASTA records found in input')

    @property
    def digest(self) -> str:
        """Stable content digest of the raw uploaded bytes.

        Returns
        -------
        str
            First 16 hex characters of the SHA-256 digest.
        """
        return self._digest

    @property
    def names(self) -> list[str]:
        """Contig names in file order.

        Returns
        -------
        list[str]
            Names in the order they appear in the file.
        """
        return list(self._fasta.keys())

    @property
    def total_length(self) -> int:
        """Total number of residues across all contigs.

        Returns
        -------
        int
            Sum of sequence lengths (from the ``.fai``, no sequence reads).
        """
        return sum(len(self._fasta[name]) for name in self._fasta.keys())

    def lengths(self) -> dict[str, int]:
        """Return per-contig sequence lengths from the ``.fai`` index.

        Returns
        -------
        dict[str, int]
            ``{contig name: length}``.
        """
        return {name: len(self._fasta[name]) for name in self._fasta.keys()}

    def get_lazy(self, name: str) -> LazySeq | None:
        """Return a lazily sliceable record for *name*, or ``None``.

        Parameters
        ----------
        name : str
            Contig name.

        Returns
        -------
        LazySeq or None
            A ``pyfaidx`` record: ``len()`` comes from the ``.fai`` and
            ``[start:stop]`` reads only that window from the file.
        """
        try:
            return self._fasta[name]
        except KeyError:
            return None

    def get_slice(self, name: str, start: int, end: int) -> str:
        """Return ``sequence[start:end]`` for contig *name*.

        Parameters
        ----------
        name : str
            Contig name (must exist).
        start : int
            0-based inclusive start.
        end : int
            0-based exclusive end.

        Returns
        -------
        str
            The requested window, read from the file.
        """
        return self._fasta[name][start:end]

    def iter_records(self) -> Iterator[tuple[str, str]]:
        """Yield ``(name, sequence)`` pairs one contig at a time.

        Only one contig's sequence is resident per iteration step, so peak
        memory while e.g. feeding the k-mer index stays at one contig.

        Yields
        ------
        tuple[str, str]
            One record per contig, in file order.
        """
        for name in self._fasta.keys():
            yield name, self._fasta[name][:]


def file_digest(path: Path) -> str:
    """Compute the content digest of a file without loading it whole.

    Streams the file through SHA-256 and matches
    :func:`core.fasta.content_digest` output for the same bytes.

    Parameters
    ----------
    path : pathlib.Path
        File to hash.

    Returns
    -------
    str
        First 16 hex characters of the SHA-256 digest.
    """
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()[:16]


def _ensure_plain_fasta(path: Path) -> Path:
    """Return a plain-FASTA path for *path*, decompressing gzip if needed.

    ``pyfaidx`` cannot index plain gzip (only bgzip), so a gzip upload is
    decompressed once, streaming, to a sibling file in the same (per-upload
    temporary) directory.

    Parameters
    ----------
    path : pathlib.Path
        Uploaded file; may be gzip-compressed.

    Returns
    -------
    pathlib.Path
        *path* itself, or the decompressed sibling for gzip input.
    """
    with open(path, 'rb') as fh:
        magic = fh.read(2)
    if magic != _GZIP_MAGIC:
        return path
    plain = path.parent / (path.name + '.plain.fasta')
    with gzip.open(path, 'rb') as src, open(plain, 'wb') as dst:
        shutil.copyfileobj(src, dst)
    return plain


def provider_from_path(path: Path) -> SequenceProvider:
    """Build a sequence provider for an uploaded FASTA file.

    Prefers a :class:`FaidxProvider` over the file on disk; any failure
    (``pyfaidx`` missing, malformed FASTA, duplicate names, …) falls back to
    the pure-Python parser, which either succeeds in memory or raises the
    app's usual user-facing :class:`ValueError` messages.

    Parameters
    ----------
    path : pathlib.Path
        The upload's ``datapath``.

    Returns
    -------
    SequenceProvider
        Faidx-backed when possible, in-memory otherwise.

    Raises
    ------
    ValueError
        If the file is not valid FASTA (raised by the fallback parser).
    """
    digest = file_digest(path)
    try:
        provider = FaidxProvider(_ensure_plain_fasta(path), digest)
        logger.info(
            'Indexed %d contig(s) with pyfaidx (digest %s)',
            len(provider.names),
            digest,
        )
        return provider
    except Exception as exc:  # noqa: BLE001 - any faidx failure falls back
        logger.info('pyfaidx unavailable or failed (%s); parsing in memory', exc)
    return InMemoryProvider(parse_fasta_bytes(path.read_bytes()))


def provider_from_fasta_input(fasta: FastaInput, scratch_dir: Path) -> SequenceProvider:
    """Build a provider for already-parsed records (e.g. a GenBank upload).

    Writes the records to a FASTA file in *scratch_dir* and indexes it with
    ``pyfaidx`` so the parsed strings can be released; falls back to serving
    them from memory when indexing fails.

    Parameters
    ----------
    fasta : FastaInput
        Parsed records (typically converted from GenBank ORIGIN blocks).
    scratch_dir : pathlib.Path
        Writable directory (the upload's own temporary directory).

    Returns
    -------
    SequenceProvider
        Faidx-backed when possible, in-memory otherwise.
    """
    try:
        path = scratch_dir / f'rd_{fasta.digest}.fasta'
        if not path.exists():
            with open(path, 'w', encoding='utf-8') as fh:
                for name, seq in fasta.records:
                    fh.write(f'>{name}\n{seq}\n')
        return FaidxProvider(path, fasta.digest)
    except Exception as exc:  # noqa: BLE001 - any faidx failure falls back
        logger.info('pyfaidx unavailable or failed (%s); serving from memory', exc)
        return InMemoryProvider(fasta)
