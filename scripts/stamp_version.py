#!/usr/bin/env python3
"""Stamp the git-tag version into every file that carries one.

The most recent tag matching ``v[0-9]*`` is the single source of truth.
This script derives the version from ``git describe`` and writes it to:

- ``Cargo.toml`` ``[package] version`` (plain ``X.Y.Z`` — maturin reads it
  via ``dynamic = ["version"]`` in ``pyproject.toml``)
- ``Cargo.lock`` (the ``dot-explorer`` package block, keeping the committed
  lock consistent with the manifest)
- ``python/dot_explorer/_version.py`` (full PEP 440 string, e.g.
  ``X.Y.Z.postN`` for commits after a tag)
- ``CITATION.cff`` ``version:`` and ``date-released:`` — only when HEAD is
  exactly at a release tag, so between releases the last released values
  stand (matching what the tag-triggered citation-cff workflow enforces)

Run with no arguments to stamp files in place (idempotent: exits 0 with no
writes when everything is current). Run with ``--check`` to fail (exit 1)
when any file is out of date, without modifying it — used in CI.
"""

from __future__ import annotations

from pathlib import Path
import re
import subprocess
import sys

REPO_ROOT = Path(__file__).resolve().parent.parent

# Used when the repository has no v* tags yet (fresh clone of an untagged
# repo). Matches the last manually released version.
FALLBACK_VERSION = '0.1.0'

DESCRIBE_RE = re.compile(
    r'^v(?P<base>[0-9]+\.[0-9]+\.[0-9]+[0-9A-Za-z.\-]*?)'
    r'(?:-(?P<distance>[0-9]+)-g(?P<hash>[0-9a-f]+))?$'
)


def git_describe(repo: Path = REPO_ROOT) -> str | None:
    """Return ``git describe`` output for the latest ``v*`` tag, or None.

    Parameters
    ----------
    repo : Path
        Repository to describe.

    Returns
    -------
    str or None
        The raw ``git describe --tags --match 'v[0-9]*'`` output, or None
        when no matching tag exists (or git is unavailable).
    """
    try:
        out = subprocess.run(
            ['git', 'describe', '--tags', '--match', 'v[0-9]*'],
            cwd=repo,
            capture_output=True,
            text=True,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return out.stdout.strip() or None


def parse_describe(text: str | None) -> tuple[str, str, bool]:
    """Map ``git describe`` output to Python and Cargo version strings.

    Parameters
    ----------
    text : str or None
        Output of :func:`git_describe`; None when the repo has no tags.

    Returns
    -------
    tuple of (str, str, bool)
        ``(py_version, cargo_version, exact)`` where ``py_version`` is the
        full PEP 440 string (``X.Y.Z`` at a tag, ``X.Y.Z.postN`` after it),
        ``cargo_version`` is the plain tag version (SemVer, no PEP 440
        suffix — avoids Cargo.lock churn between tags), and ``exact`` is
        True when HEAD sits exactly on the tag.
    """
    if not text:
        return FALLBACK_VERSION, FALLBACK_VERSION, False
    m = DESCRIBE_RE.match(text)
    if not m:
        return FALLBACK_VERSION, FALLBACK_VERSION, False
    base = m.group('base')
    distance = m.group('distance')
    if distance is None or int(distance) == 0:
        return base, base, True
    return f'{base}.post{distance}', base, False


def tag_date(version: str, repo: Path = REPO_ROOT) -> str | None:
    """Return the commit date (YYYY-MM-DD) of the tag ``v<version>``.

    Parameters
    ----------
    version : str
        Version without the ``v`` prefix.
    repo : Path
        Repository to query.

    Returns
    -------
    str or None
        ISO date of the tagged commit, or None when the tag is missing.
    """
    try:
        out = subprocess.run(
            ['git', 'log', '-1', '--format=%cs', f'v{version}'],
            cwd=repo,
            capture_output=True,
            text=True,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return out.stdout.strip() or None


def _sub_once(text: str, pattern: str, replacement: str, path: Path) -> str:
    new, n = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
    if n != 1:
        raise SystemExit(f'stamp_version: could not find version line in {path}')
    return new


def stamp_cargo_toml(text: str, version: str, path: Path) -> str:
    """Return ``Cargo.toml`` content with the ``[package]`` version updated.

    Parameters
    ----------
    text : str
        Current file content.
    version : str
        Plain SemVer version to write.
    path : Path
        File path, for error messages only.

    Returns
    -------
    str
        Updated content.
    """
    # Anchor on the [package] section so dependency `version =` keys are
    # never touched.
    pattern = r'(\[package\][^\[]*?^version = ")[^"]+(")'
    return _sub_once(text, pattern, rf'\g<1>{version}\g<2>', path)


def stamp_cargo_lock(text: str, version: str, path: Path) -> str:
    """Return ``Cargo.lock`` content with the dot-explorer block updated.

    Parameters
    ----------
    text : str
        Current file content.
    version : str
        Plain SemVer version to write.
    path : Path
        File path, for error messages only.

    Returns
    -------
    str
        Updated content.
    """
    pattern = r'(name = "dot-explorer"\nversion = ")[^"]+(")'
    return _sub_once(text, pattern, rf'\g<1>{version}\g<2>', path)


def stamp_version_py(version: str) -> str:
    """Return the full content of ``_version.py`` for a version.

    Parameters
    ----------
    version : str
        PEP 440 version string.

    Returns
    -------
    str
        Complete module source.
    """
    return (
        '"""Version stamped from the latest git tag by scripts/stamp_version.py.\n'
        '\n'
        'Do not edit by hand: tag a release (``git tag vX.Y.Z``) and the\n'
        'pre-commit hook / CI stamp step rewrites this file.\n'
        '"""\n'
        '\n'
        f"__version__ = '{version}'\n"
    )


def stamp_citation(text: str, version: str, date: str | None, path: Path) -> str:
    """Return ``CITATION.cff`` content with version and release date updated.

    Parameters
    ----------
    text : str
        Current file content.
    version : str
        Released version (HEAD is at the matching tag).
    date : str or None
        ISO release date; the date line is left alone when None.
    path : Path
        File path, for error messages only.

    Returns
    -------
    str
        Updated content.
    """
    text = _sub_once(text, r'^(version: ).+$', rf'\g<1>{version}', path)
    if date:
        text = _sub_once(text, r'^(date-released: ).+$', rf'\g<1>{date}', path)
    return text


def compute_updates(repo: Path = REPO_ROOT) -> dict[Path, str]:
    """Compute the target content for every stamped file.

    Parameters
    ----------
    repo : Path
        Repository root.

    Returns
    -------
    dict of Path to str
        Mapping of file path to the content it should have. Files already
        current are included; the caller diffs against disk.
    """
    py_version, cargo_version, exact = parse_describe(git_describe(repo))
    updates: dict[Path, str] = {}

    cargo_toml = repo / 'Cargo.toml'
    updates[cargo_toml] = stamp_cargo_toml(
        cargo_toml.read_text(), cargo_version, cargo_toml
    )
    cargo_lock = repo / 'Cargo.lock'
    if cargo_lock.exists():
        updates[cargo_lock] = stamp_cargo_lock(
            cargo_lock.read_text(), cargo_version, cargo_lock
        )
    updates[repo / 'python' / 'dot_explorer' / '_version.py'] = stamp_version_py(
        py_version
    )
    if exact:
        citation = repo / 'CITATION.cff'
        updates[citation] = stamp_citation(
            citation.read_text(), cargo_version, tag_date(cargo_version, repo), citation
        )
    return updates


def main(argv: list[str]) -> int:
    """Stamp (or with ``--check`` verify) all version-bearing files.

    Parameters
    ----------
    argv : list of str
        Command-line arguments (excluding the program name).

    Returns
    -------
    int
        Process exit code: 0 when current/stamped, 1 when ``--check`` finds
        stale files.
    """
    check = '--check' in argv
    stale = []
    for path, content in compute_updates().items():
        current = path.read_text() if path.exists() else None
        if current == content:
            continue
        stale.append(path)
        if not check:
            path.write_text(content)
            print(f'stamp_version: updated {path.relative_to(REPO_ROOT)}')
    if check and stale:
        names = ', '.join(str(p.relative_to(REPO_ROOT)) for p in stale)
        print(f'stamp_version: stale version in: {names}', file=sys.stderr)
        print('Run `python scripts/stamp_version.py` and commit.', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
