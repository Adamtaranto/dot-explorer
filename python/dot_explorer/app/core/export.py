"""Reordered / reoriented FASTA export from in-memory sequences.

The k-mer method stores sequences inside the ``CrossIndex`` (which has its
own ``write_fasta``), but PAF-import and external-tool methods only carry
coordinates — for those the app keeps the parsed upload
(:class:`core.fasta.FastaInput`) and exports directly from its records
using this module.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
import dataclasses
import io
import re
from typing import TYPE_CHECKING, Iterable, Sequence
import zipfile

if TYPE_CHECKING:  # pragma: no cover - only for type checkers
    from dot_explorer.paf_io import PafRecord

#: One sequence region to export: ``(contig, start, end, strand, side)`` with
#: 0-based half-open coordinates.  ``side`` is a free label (``'query'`` /
#: ``'target'``) carried into the FASTA header; it plays no part in
#: de-duplication.
Region = tuple[str, int, int, str, str]


def reordered_fasta_text(
    records: Sequence[tuple[str, str]],
    order: Iterable[str],
    reverse: set[str],
    line_width: int = 60,
) -> str:
    """Serialise sequence records to FASTA in a given order and orientation.

    Contigs named in *order* are written first, in that order; any remaining
    contigs from *records* follow in their original order (so an ordering
    derived from alignments — which may not cover every contig — still
    exports the complete assembly).  Contigs named in *reverse* are written
    reverse-complemented with a ``reverse_complement`` note in the header,
    matching :meth:`dot_explorer.paf_io.CrossIndex.write_fasta` output.

    Parameters
    ----------
    records : sequence of (str, str)
        ``(name, sequence)`` pairs, e.g. :attr:`core.fasta.FastaInput.records`.
    order : iterable of str
        Contig names in the desired output order.  Names not present in
        *records* are ignored.
    reverse : set[str]
        Contig names to reverse-complement.
    line_width : int, optional
        Wrap sequence lines at this many bases; ``0`` or negative writes
        each sequence on one line.  Default is ``60``.

    Returns
    -------
    str
        The FASTA text.
    """
    from dot_explorer.paf_io import reverse_complement  # noqa: PLC0415 - lazy

    seq_map = dict(records)
    names = [n for n in order if n in seq_map]
    seen = set(names)
    names.extend(n for n, _ in records if n not in seen)

    chunks: list[str] = []
    for name in names:
        seq = seq_map[name]
        if name in reverse:
            seq = reverse_complement(seq)
            chunks.append(f'>{name} reverse_complement\n')
        else:
            chunks.append(f'>{name}\n')
        if line_width and line_width > 0:
            for i in range(0, len(seq), line_width):
                chunks.append(seq[i : i + line_width] + '\n')
        else:
            chunks.append(seq + '\n')
    return ''.join(chunks)


def selected_regions_fasta(
    regions: Iterable[Region],
    get_sequence: Callable[[str, int, int], str | None],
    line_width: int = 60,
) -> str:
    """Serialise the sequences of selected alignment regions as multi-FASTA.

    Regions sharing a contig, coordinates and strand are written once, no
    matter how many selected alignments (or which sides) produced them.
    Minus-strand regions are reverse-complemented so each record reads in
    alignment orientation, matching the report's copy buttons.

    Parameters
    ----------
    regions : Iterable[Region]
        ``(contig, start, end, strand, side)`` entries, in selection order.
    get_sequence : callable
        ``(contig, start, end) -> str | None`` slicing the forward-strand
        sequence; ``None`` skips the region (sequence not available).
    line_width : int, optional
        Wrap sequence lines at this many bases; ``0`` writes one line.
        Default is ``60``.

    Returns
    -------
    str
        FASTA text with headers ``>contig:start-end(strand) side=...``
        (1-based inclusive coordinates, as in the report's detail bar).
    """
    from dot_explorer.paf_io import reverse_complement  # noqa: PLC0415 - lazy

    seen: set[tuple[str, int, int, str]] = set()
    chunks: list[str] = []
    for contig, start, end, strand, side in regions:
        strand = '-' if strand == '-' else '+'
        key = (contig, int(start), int(end), strand)
        if key in seen or end <= start:
            continue
        seq = get_sequence(contig, int(start), int(end))
        if seq is None:
            continue
        seen.add(key)
        if strand == '-':
            seq = reverse_complement(seq)
        chunks.append(f'>{contig}:{start + 1}-{end}({strand}) side={side}\n')
        if line_width and line_width > 0:
            for i in range(0, len(seq), line_width):
                chunks.append(seq[i : i + line_width] + '\n')
        else:
            chunks.append(seq + '\n')
    return ''.join(chunks)


#: Name of the archive member holding contigs without a cluster assignment.
UNASSIGNED_CLUSTER = 'unassigned'


def cluster_fasta_zip(
    records: Sequence[tuple[str, str]],
    assignments: Mapping[str, str],
    reverse: set[str],
    order: Iterable[str] = (),
    line_width: int = 60,
) -> bytes:
    """Bundle one multi-FASTA per cluster into a zip archive.

    Parameters
    ----------
    records : sequence of (str, str)
        ``(name, sequence)`` pairs for every contig.
    assignments : Mapping[str, str]
        Contig name -> cluster name (e.g. from
        :attr:`dot_explorer.similarity.ClusterResult.assignments`).  Contigs
        absent from the mapping go to ``unassigned.fasta``.
    reverse : set[str]
        Contigs to write reverse-complemented (the plot's displayed
        orientation, manual flips included).
    order : Iterable[str], optional
        Contig names in display order; members of each cluster are written
        in this order, unlisted contigs afterwards.  Default: record order.
    line_width : int, optional
        FASTA line wrap.  Default is ``60``.

    Returns
    -------
    bytes
        The zip archive (deflate-compressed), members named
        ``<cluster>.fasta`` in sorted cluster order.
    """
    by_name = dict(records)
    rank = {name: i for i, name in enumerate(order)}
    members: dict[str, list[str]] = {}
    for name, _seq in records:
        cluster = assignments.get(name, UNASSIGNED_CLUSTER)
        members.setdefault(cluster, []).append(name)
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, 'w', compression=zipfile.ZIP_DEFLATED) as zf:
        for cluster in sorted(members, key=_cluster_sort_key):
            names = sorted(members[cluster], key=lambda n: rank.get(n, len(rank)))
            text = reordered_fasta_text(
                [(n, by_name[n]) for n in names], names, reverse, line_width
            )
            zf.writestr(f'{_safe_member(cluster)}.fasta', text)
    return buf.getvalue()


def _cluster_sort_key(name: str) -> tuple[int, int | str]:
    """Sort ``cluster_2`` before ``cluster_10``; the unassigned bucket last."""
    if name == UNASSIGNED_CLUSTER:
        return (2, 0)
    head, _sep, tail = name.rpartition('_')
    if head and tail.isdigit():
        return (0, int(tail))
    return (1, name)


def _safe_member(name: str) -> str:
    """Make a cluster name safe as a zip member name."""
    return ''.join(c if c.isalnum() or c in '-_.' else '_' for c in name) or 'cluster'


_CIGAR_OP = re.compile(r'(\d+)([MIDNSHP=X])')


def _reverse_cigar(cigar: str) -> str:
    """Reverse the operation order of a CIGAR string (``10M2I5M`` → ``5M2I10M``)."""
    ops = _CIGAR_OP.findall(cigar)
    if ''.join(n + op for n, op in ops) != cigar:
        return cigar  # not a plain CIGAR: leave it untouched
    return ''.join(n + op for n, op in reversed(ops))


def reoriented_paf_records(
    records: Iterable['PafRecord'],
    reverse_query: set[str],
    reverse_target: set[str] = frozenset(),
) -> list['PafRecord']:
    """Rewrite PAF records against reverse-complemented contigs.

    The reordered-FASTA export writes contigs in *reverse_query* (and, for
    a self-comparison, *reverse_target*) reverse-complemented; a PAF
    exported alongside must describe the same sequences.  For each
    mirrored axis the coordinates become ``len - end, len - start`` and
    the strand flips once, so a record whose query and target are both
    flipped keeps its strand.  A CIGAR walks the alignment from the
    query's start, so its operations reverse when exactly one axis flips
    (two flips restore the original order).  Everything else — names,
    lengths, match counts, tags — is unchanged, and records on unflipped
    contigs are returned as-is.

    Parameters
    ----------
    records : Iterable[PafRecord]
        Genomic-orientation records (as parsed or computed).
    reverse_query : set[str]
        Query contig names exported reverse-complemented.
    reverse_target : set[str], optional
        Target contig names exported reverse-complemented (self mode).

    Returns
    -------
    list[PafRecord]
        New records; the inputs are not modified.
    """
    out = []
    for rec in records:
        rq = rec.query_name in reverse_query
        rt = rec.target_name in reverse_target
        if not (rq or rt):
            out.append(rec)
            continue
        changes: dict = {}
        if rq:
            changes['query_start'] = rec.query_len - rec.query_end
            changes['query_end'] = rec.query_len - rec.query_start
        if rt:
            changes['target_start'] = rec.target_len - rec.target_end
            changes['target_end'] = rec.target_len - rec.target_start
        if rq != rt:
            changes['strand'] = '-' if rec.strand == '+' else '+'
            if rec.cigar:
                flipped = _reverse_cigar(rec.cigar)
                changes['cigar'] = flipped
                if 'cg' in rec.tags:
                    tags = dict(rec.tags)
                    tags['cg'] = flipped
                    changes['tags'] = tags
        out.append(dataclasses.replace(rec, **changes))
    return out
