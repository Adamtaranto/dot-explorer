"""Cross-checks between an uploaded assembly and a pre-computed PAF.

When a user pairs a PAF alignment with supplementary query / target
assemblies (for sequence fetching and the reordered-FASTA export), the
contig names must line up with the matching PAF column — a mismatch
silently produces empty or wrong exports.  These checks turn the common
mistakes (wrong file, swapped query/target, partial assemblies) into
explicit errors and warnings.
"""

from __future__ import annotations

from collections.abc import Iterable


def _preview(names: set[str], limit: int = 5) -> str:
    shown = ', '.join(sorted(names)[:limit])
    return shown + (', …' if len(names) > limit else '')


def validate_paf_names(
    assembly_names: Iterable[str],
    paf_names: Iterable[str],
    other_paf_names: Iterable[str],
    role: str = 'query',
) -> tuple[list[str], list[str]]:
    """Check an uploaded assembly's contig names against one PAF column.

    Parameters
    ----------
    assembly_names : Iterable[str]
        Contig names from the uploaded assembly for *role*.
    paf_names : Iterable[str]
        Names from the PAF column this assembly must cover (column 1 for
        the query, column 6 for the target).
    other_paf_names : Iterable[str]
        Names from the PAF's other column, for the swapped-inputs hint.
    role : str, optional
        ``'query'`` or ``'target'``, used in the messages.  Default
        ``'query'``.

    Returns
    -------
    tuple[list[str], list[str]]
        ``(errors, warnings)``.  Errors are conditions under which
        sequences the PAF refers to cannot be found — no assembly name in
        the PAF column (with a hint when they match the *other* column),
        or PAF names missing from the assembly.  Warnings cover assembly
        contigs with no alignments and names that appear in both PAF
        columns (ambiguous role).
    """
    assembly = set(assembly_names)
    paf = set(paf_names)
    other = set(other_paf_names)
    errors: list[str] = []
    warnings: list[str] = []
    other_role = 'target' if role == 'query' else 'query'

    matched = assembly & paf
    if assembly and paf and not matched:
        msg = (
            f'None of the uploaded {role} assembly contig names appear in '
            f'the PAF {role} column.'
        )
        if assembly & other:
            msg += (
                f' They DO match the PAF {other_role} column: the file looks '
                f'like the {other_role} assembly, or the PAF was generated '
                'with query and target swapped.'
            )
        return [msg], warnings

    missing_from_assembly = paf - assembly
    if missing_from_assembly:
        errors.append(
            f'{len(missing_from_assembly)} PAF {role} name(s) are not in the '
            f'uploaded {role} assembly ({_preview(missing_from_assembly)}); '
            'their sequences cannot be found. Upload the assembly the PAF '
            'was generated from.'
        )

    missing_from_paf = assembly - paf
    if missing_from_paf:
        warnings.append(
            f'{len(missing_from_paf)} {role} assembly contig(s) have no '
            f'alignments in the PAF ({_preview(missing_from_paf)}); they '
            'will appear unplaced in the reordered FASTA.'
        )

    ambiguous = paf & other
    if ambiguous and role == 'query':
        warnings.append(
            f'{len(ambiguous)} name(s) appear in BOTH the PAF query and '
            f'target columns ({_preview(ambiguous)}); their role is '
            'ambiguous and reordering/orientation may be wrong.'
        )
    return errors, warnings


def validate_query_names(
    assembly_names: Iterable[str],
    paf_query_names: Iterable[str],
    paf_target_names: Iterable[str],
) -> list[str]:
    """Check uploaded query-assembly contig names against a PAF's name sets.

    Thin wrapper over :func:`validate_paf_names` returning every finding as
    a single list of messages (errors first).

    Parameters
    ----------
    assembly_names : Iterable[str]
        Contig names from the uploaded query assembly FASTA.
    paf_query_names : Iterable[str]
        Names from the PAF query-name column (column 1).
    paf_target_names : Iterable[str]
        Names from the PAF target-name column (column 6).

    Returns
    -------
    list[str]
        Human-readable messages, empty when everything lines up.
    """
    errors, warnings = validate_paf_names(
        assembly_names, paf_query_names, paf_target_names, 'query'
    )
    return errors + warnings


def validate_annotation_names(
    assembly_names: Iterable[str],
    annotation_names: Iterable[str],
    role: str,
) -> list[str]:
    """Check a GFF's sequence names against the assembly it annotates.

    A GFF whose ``seqname`` column does not match the FASTA headers is the
    single most common annotation mistake — a different assembly version,
    or names rewritten by an intermediate tool.  Nothing errors: the
    features simply never appear, which is easy to misread as "the plot is
    broken".

    Parameters
    ----------
    assembly_names : Iterable[str]
        Contig names from the assembly FASTA for this role.
    annotation_names : Iterable[str]
        Distinct ``seqname`` values from the parsed GFF.
    role : str
        ``'query'`` or ``'target'``, used in the message.

    Returns
    -------
    list[str]
        Human-readable warnings, empty when every GFF sequence is present.
    """
    assembly = set(assembly_names)
    annotated = set(annotation_names)
    if not assembly or not annotated:
        return []

    def _preview(names: set[str], limit: int = 5) -> str:
        shown = ', '.join(sorted(names)[:limit])
        return shown + (', …' if len(names) > limit else '')

    missing = annotated - assembly
    if not missing:
        return []
    if missing == annotated:
        # Nothing lines up at all: almost always the wrong file or the
        # wrong role, so say so rather than listing every name.
        return [
            f'None of the {role} GFF sequence names match the {role} '
            f'assembly ({_preview(annotated)}); no features will be drawn. '
            'Check the GFF belongs to this assembly and is assigned to the '
            'right role.'
        ]
    return [
        f'{len(missing)} sequence name(s) in the {role} GFF are not in the '
        f'{role} assembly ({_preview(missing)}); features on them will not '
        'be drawn.'
    ]
