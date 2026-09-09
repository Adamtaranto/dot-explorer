#!/usr/bin/env python3
"""Stamp the git-tag version into every file that carries one.

The most recent tag matching ``v[0-9]*`` is the single source of truth.
This script derives the version from ``git describe`` and writes it to:

- ``Cargo.toml`` ``[package] version`` (plain ``X.Y.Z`` — maturin reads it
  via ``dynamic = ["version"]`` in ``pyproject.toml``)
- ``Cargo.lock`` (the ``dot-explorer`` package block, keeping the committed
  lock consistent with the manifest)
- ``python/dot_explorer/_version.py`` and ``CITATION.cff``
  (``version:`` / ``date-released:``) — **only when HEAD is exactly at a
  release tag**, so between releases the last released values stand. That
  keeps ``__version__`` in step with the wheel metadata built from
  ``Cargo.toml``, stops an ordinary commit from dirtying the tree, and
  matches what the tag-triggered citation-cff workflow enforces.

No ``.postN`` suffix is written to a tracked file: the wheel version comes
from ``Cargo.toml``, so a suffix in ``_version.py`` alone would just
disagree with ``importlib.metadata.version('dot-explorer')``.

Development builds are identified instead through an **untracked**
``python/dot_explorer/_version_local.py``, rewritten on every run between
tags and deleted at a tag. It carries a PEP 440 local version —
``X.Y.Z+<short hash>``, plus ``.dirty`` when the tree has uncommitted
changes — and ``_version.py`` imports it when present. Local versions are
exactly the mechanism PEP 440 provides for "this release, built from a
specific local state", so a dev build says where it came from while every
tracked file stays put.

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


def parse_describe(text: str | None) -> tuple[str, bool]:
    """Map ``git describe`` output to a version and a tagged-ness flag.

    The version is always the plain tag version: no ``.postN`` development
    suffix. A dev suffix in ``_version.py`` would disagree with the wheel's
    own metadata (which comes from ``Cargo.toml``, stamped without one), and
    would rewrite the file on every commit made after a tag.

    Parameters
    ----------
    text : str or None
        Output of :func:`git_describe`; None when the repo has no tags.

    Returns
    -------
    tuple of (str, bool)
        ``(version, exact)`` — the plain tag version, and whether HEAD sits
        exactly on that tag.
    """
    if not text:
        return FALLBACK_VERSION, False
    m = DESCRIBE_RE.match(text)
    if not m:
        return FALLBACK_VERSION, False
    base = m.group('base')
    distance = m.group('distance')
    return base, distance is None or int(distance) == 0


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
        '\n'
        'Builds made between releases additionally write an untracked\n'
        '``_version_local.py`` carrying a PEP 440 local version\n'
        '(``X.Y.Z+<short hash>``); when present it wins, so a development\n'
        'build identifies the commit it came from without this file — or\n'
        'Cargo.toml — changing on every commit.\n'
        '"""\n'
        '\n'
        f"__version__ = '{version}'\n"
        '\n'
        'try:  # pragma: no cover - only present in development builds\n'
        '    from dot_explorer._version_local import __version__  # noqa: F401,F811\n'
        'except ImportError:\n'
        '    pass\n'
    )


def stamp_local_version_py(version: str) -> str:
    """Return the content of the untracked ``_version_local.py``.

    Parameters
    ----------
    version : str
        Full local version, e.g. ``0.1.0+1a2b3c4``.

    Returns
    -------
    str
        Complete module source.
    """
    return (
        '"""Development-build version (untracked, written at build time).\n'
        '\n'
        'Generated by scripts/stamp_version.py for builds made between release\n'
        'tags, and deleted again at a tag. Never commit this file.\n'
        '"""\n'
        '\n'
        f"__version__ = '{version}'\n"
    )


def short_hash(repo: Path = REPO_ROOT) -> str | None:
    """Return the abbreviated HEAD commit hash, or None outside a repo.

    Parameters
    ----------
    repo : Path
        Repository to query.

    Returns
    -------
    str or None
        Short hash, or None when git is unavailable (e.g. an sdist build).
    """
    try:
        out = subprocess.run(
            ['git', 'rev-parse', '--short', 'HEAD'],
            cwd=repo,
            capture_output=True,
            text=True,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return out.stdout.strip() or None


def is_dirty(repo: Path = REPO_ROOT) -> bool:
    """Whether the working tree has uncommitted changes.

    Parameters
    ----------
    repo : Path
        Repository to query.

    Returns
    -------
    bool
        True when ``git status --porcelain`` reports anything.
    """
    try:
        out = subprocess.run(
            ['git', 'status', '--porcelain'],
            cwd=repo,
            capture_output=True,
            text=True,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return False
    return bool(out.stdout.strip())


def local_version(version: str, repo: Path = REPO_ROOT) -> str | None:
    """Compose the development local version for the current checkout.

    Parameters
    ----------
    version : str
        The last released version.
    repo : Path
        Repository to query.

    Returns
    -------
    str or None
        ``X.Y.Z+<hash>`` (with ``.dirty`` appended when the tree has
        uncommitted changes), or None when the hash cannot be determined.
    """
    sha = short_hash(repo)
    if sha is None:
        return None
    return f'{version}+{sha}.dirty' if is_dirty(repo) else f'{version}+{sha}'


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

    Notes
    -----
    ``_version.py`` and ``CITATION.cff`` are only rewritten when HEAD sits
    exactly on a release tag. Between tags they keep the last released
    values, which is also what ``Cargo.toml`` (and therefore the wheel's
    metadata) carries — so ``dot_explorer.__version__`` and
    ``importlib.metadata.version('dot-explorer')`` always agree, and an
    ordinary commit never dirties them.
    """
    version, exact = parse_describe(git_describe(repo))
    updates: dict[Path, str] = {}

    cargo_toml = repo / 'Cargo.toml'
    updates[cargo_toml] = stamp_cargo_toml(cargo_toml.read_text(), version, cargo_toml)
    cargo_lock = repo / 'Cargo.lock'
    if cargo_lock.exists():
        updates[cargo_lock] = stamp_cargo_lock(
            cargo_lock.read_text(), version, cargo_lock
        )
    if exact:
        updates[repo / 'python' / 'dot_explorer' / '_version.py'] = stamp_version_py(
            version
        )
        citation = repo / 'CITATION.cff'
        updates[citation] = stamp_citation(
            citation.read_text(), version, tag_date(version, repo), citation
        )
    return updates


def refresh_local_version(repo: Path = REPO_ROOT) -> str | None:
    """Write (or remove) the untracked development-version override.

    Deliberately separate from :func:`compute_updates`: the file is
    untracked and changes with every commit, so ``--check`` must not see it.

    Parameters
    ----------
    repo : Path
        Repository root.

    Returns
    -------
    str or None
        The local version written, or None when the file was removed (at a
        release tag) or could not be determined (no git).
    """
    path = repo / 'python' / 'dot_explorer' / '_version_local.py'
    version, exact = parse_describe(git_describe(repo))
    if exact:
        # At a tag the released version is the truth; never let a stale dev
        # override ride along into a release wheel.
        if path.exists():
            path.unlink()
            print(f'stamp_version: removed {path.relative_to(repo)} (at a tag)')
        return None
    dev = local_version(version, repo)
    if dev is None:
        return None
    content = stamp_local_version_py(dev)
    if not path.exists() or path.read_text() != content:
        path.write_text(content)
        print(f'stamp_version: {path.relative_to(repo)} -> {dev}')
    return dev


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
    if not check:
        refresh_local_version()
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
