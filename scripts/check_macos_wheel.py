#!/usr/bin/env python3
"""Fail if a macOS wheel links a library outside the OS.

A Mach-O extension records an absolute path for every dynamic dependency.
Build machines carry libraries users do not: GitHub's macOS runners ship
Homebrew, so a wheel that picks up ``/opt/homebrew/opt/xz/lib/liblzma.5.dylib``
imports fine in CI and dies with ``Library not loaded`` on any Mac without
that exact path — which is what shipped in 0.1.0 and 0.2.0.

Only the OS-provided prefixes are safe to depend on. Anything under a package
manager's prefix must be statically linked or bundled into the wheel instead.

Usage
-----
    python scripts/check_macos_wheel.py dist/*.whl
"""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile

#: Load-command prefixes that resolve on a stock macOS install. ``/usr/lib``
#: and ``/System`` come from the dyld shared cache; the ``@`` prefixes are
#: relative to the loading binary, i.e. bundled inside the wheel.
SAFE_PREFIXES = ('/usr/lib/', '/System/', '@loader_path', '@rpath', '@executable_path')


def dependencies(binary: Path) -> list[str]:
    """Return the dynamic-library paths a Mach-O binary loads.

    Parameters
    ----------
    binary : Path
        Path to a ``.so``/``.dylib``.

    Returns
    -------
    list of str
        One path per ``otool -L`` entry, excluding the binary's own id line.
    """
    out = subprocess.run(
        ['otool', '-L', str(binary)], capture_output=True, text=True, check=True
    ).stdout
    # First line is the file name; each following line is "<path> (compat ...)".
    return [line.split(' (')[0].strip() for line in out.splitlines()[1:] if line.strip()]


def unsafe_dependencies(binary: Path) -> list[str]:
    """Return the dependencies that will not resolve on a stock macOS.

    Parameters
    ----------
    binary : Path
        Path to a ``.so``/``.dylib``.

    Returns
    -------
    list of str
        Paths outside :data:`SAFE_PREFIXES`.
    """
    return [
        dep
        for dep in dependencies(binary)
        if not dep.startswith(SAFE_PREFIXES) and not dep.startswith('/usr/lib')
    ]


def check_wheel(wheel: Path) -> list[tuple[str, str]]:
    """Check every Mach-O object inside a wheel.

    Parameters
    ----------
    wheel : Path
        The ``.whl`` to inspect.

    Returns
    -------
    list of (str, str)
        ``(member, dependency)`` pairs that would fail to load.
    """
    bad: list[tuple[str, str]] = []
    with tempfile.TemporaryDirectory() as tmp:
        with zipfile.ZipFile(wheel) as zf:
            members = [n for n in zf.namelist() if n.endswith(('.so', '.dylib'))]
            for name in members:
                extracted = Path(zf.extract(name, tmp))
                bad += [(name, dep) for dep in unsafe_dependencies(extracted)]
    return bad


def main(argv: list[str]) -> int:
    """Check each wheel named on the command line.

    Parameters
    ----------
    argv : list of str
        Wheel paths.

    Returns
    -------
    int
        0 when every wheel is clean, 1 otherwise.
    """
    if not argv:
        print('usage: check_macos_wheel.py <wheel> [<wheel> ...]', file=sys.stderr)
        return 1
    failed = False
    for arg in argv:
        wheel = Path(arg)
        bad = check_wheel(wheel)
        if bad:
            failed = True
            print(f'FAIL {wheel.name}', file=sys.stderr)
            for member, dep in bad:
                print(f'  {member} -> {dep}', file=sys.stderr)
        else:
            print(f'ok   {wheel.name}')
    if failed:
        print(
            '\nThese paths exist on the build machine but not on a stock macOS '
            'install, so the wheel would fail to import. Link the library '
            'statically or bundle it into the wheel.',
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
