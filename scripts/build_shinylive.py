#!/usr/bin/env python3
"""Assemble and export the Shinylive (browser) build of the app.

``shinylive export`` bundles *everything* under the directory it is given, so
it cannot run against ``python/dot_explorer/app/`` directly: the wasm wheel has
to sit next to ``app.py`` at export time, and a ``.whl`` inside the package
tree would end up nested inside the built wheel.

So the app is copied to a staging directory outside the package, the wasm wheel
is dropped in beside it, and the staging copy is exported.  ``app.py`` finds the
wheel through ``APP_DIR.glob('wheels/*.whl')`` and micropip-installs it from the
Pyodide virtual filesystem at startup, exactly as before.

Usage
-----
    python scripts/build_shinylive.py [--out site/app]

Expects the Pyodide-compatible wheel in ``target/wheels/`` — build it with the
pinned toolchain first, as described in ``python/dot_explorer/app/README.md``.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import subprocess
import sys

REPO_ROOT = Path(__file__).resolve().parent.parent
APP_SRC = REPO_ROOT / 'python' / 'dot_explorer' / 'app'
STAGING = REPO_ROOT / 'build' / 'shinylive'
WHEEL_DIR = REPO_ROOT / 'target' / 'wheels'

#: The app's wheel is the Pyodide 0.27 / CPython 3.12 / Emscripten 3.1.58 one.
#: The PEP 783 `pyemscripten_*` wheel published to PyPI is a different ABI and
#: must never be used here (see the app README's "two wasm wheels" warning).
WHEEL_GLOB = '*cp312-cp312-emscripten_3_1_58_wasm32.whl'


def find_wheel(wheel_dir: Path = WHEEL_DIR) -> Path:
    """Locate exactly one Pyodide-compatible wasm wheel.

    Parameters
    ----------
    wheel_dir : Path
        Directory to search (maturin's output directory).

    Returns
    -------
    Path
        The wheel to bundle.

    Raises
    ------
    SystemExit
        If no wheel matches, or more than one does — a stale wheel alongside a
        fresh one is the classic cause of an app that hangs on the splash
        screen with a micropip error.
    """
    matches = sorted(wheel_dir.glob(WHEEL_GLOB))
    if not matches:
        raise SystemExit(
            f'No wasm wheel matching {WHEEL_GLOB} in {wheel_dir}.\n'
            'Build it first with the pinned toolchain — see '
            'python/dot_explorer/app/README.md.'
        )
    if len(matches) > 1:
        names = '\n  '.join(w.name for w in matches)
        raise SystemExit(
            f'Multiple wasm wheels in {wheel_dir}:\n  {names}\n'
            'Remove the stale ones (rm -rf target/wheels) and rebuild.'
        )
    return matches[0]


def stage_app(wheel: Path, staging: Path = STAGING) -> Path:
    """Copy the app next to the wasm wheel, ready for export.

    Parameters
    ----------
    wheel : Path
        The wasm wheel to bundle.
    staging : Path
        Staging directory; replaced if it already exists.

    Returns
    -------
    Path
        The staging directory.
    """
    if staging.exists():
        shutil.rmtree(staging)
    staging.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(
        APP_SRC, staging, ignore=shutil.ignore_patterns('__pycache__', '*.pyc')
    )
    wheels = staging / 'wheels'
    wheels.mkdir()
    shutil.copy2(wheel, wheels / wheel.name)
    return staging


def main(argv: list[str] | None = None) -> int:
    """Stage the app and run ``shinylive export``.

    Parameters
    ----------
    argv : list of str or None
        Argument list; ``None`` reads ``sys.argv[1:]``.

    Returns
    -------
    int
        Process exit code.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        '--out',
        default=str(REPO_ROOT / 'site' / 'app'),
        help='Export destination (default: %(default)s).',
    )
    parser.add_argument(
        '--stage-only',
        action='store_true',
        help='Assemble the staging directory but skip the export.',
    )
    args = parser.parse_args(argv)

    wheel = find_wheel()
    staging = stage_app(wheel)
    print(f'staged {APP_SRC.relative_to(REPO_ROOT)} + {wheel.name} -> {staging}')
    if args.stage_only:
        return 0

    # shinylive only mkdir()s the leaf, so the parent has to exist. In CI
    # `zensical build` has already created site/; locally it may not have.
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)

    result = subprocess.run(
        ['shinylive', 'export', str(staging), args.out], check=False
    )
    if result.returncode != 0:
        return result.returncode
    print(f'exported to {args.out}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
