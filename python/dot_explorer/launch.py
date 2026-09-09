"""Console-script launcher for the bundled Shiny app.

The app ships inside the wheel as package data at ``dot_explorer/app/`` rather
than as an importable subpackage: ``app.py`` uses flat ``core.*`` imports that
resolve once its own directory is on ``sys.path``.  That is what ``shiny run``
does with a file path, and it is also what ``shinylive export`` does in the
browser build, so keeping the layout flat lets one copy of the app serve both.

``shiny`` is an optional dependency (it pulls ``orjson``, which has no
free-threaded wheel, so it cannot be required without breaking the 3.14t
build), hence the import lives inside :func:`main` behind a helpful error.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

#: The packaged app directory (contains ``app.py``, ``core/`` and ``www/``).
APP_DIR = Path(__file__).parent / 'app'

_MISSING_SHINY = (
    'The dot-explorer app needs Shiny, which is an optional dependency.\n'
    'Install it with:\n\n'
    '    pip install "dot-explorer[app]"\n'
)


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    """Parse the launcher's command-line arguments.

    Parameters
    ----------
    argv : list of str or None
        Argument list; ``None`` reads ``sys.argv[1:]``.

    Returns
    -------
    argparse.Namespace
        Parsed options.
    """
    parser = argparse.ArgumentParser(
        prog='dot-explorer-app',
        description='Launch the dot-explorer Shiny app in a local browser.',
    )
    parser.add_argument(
        '--host',
        default='127.0.0.1',
        help='Interface to bind (default: %(default)s). Use 0.0.0.0 to expose '
        'the app on the network, e.g. when running on an HPC node.',
    )
    parser.add_argument(
        '--port',
        type=int,
        default=8000,
        help='Port to listen on (default: %(default)s).',
    )
    parser.add_argument(
        '--no-browser',
        action='store_true',
        help='Do not open a browser window (useful over SSH).',
    )
    parser.add_argument(
        '--version',
        action='store_true',
        help='Print the dot-explorer version and exit.',
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Run the packaged Shiny app.

    Parameters
    ----------
    argv : list of str or None
        Argument list; ``None`` reads ``sys.argv[1:]``.

    Returns
    -------
    int
        Process exit code: ``0`` on a clean shutdown, ``1`` when Shiny is not
        installed or the packaged app is missing.
    """
    args = _parse_args(argv)

    if args.version:
        from dot_explorer import __version__  # noqa: PLC0415 - cheap, on demand

        print(__version__)
        return 0

    try:
        from shiny import run_app  # noqa: PLC0415 - optional dependency
    except ImportError:
        print(_MISSING_SHINY, file=sys.stderr)
        return 1

    app_file = APP_DIR / 'app.py'
    if not app_file.is_file():
        print(
            f'The packaged app is missing (looked for {app_file}). This '
            'usually means an incomplete install — try reinstalling '
            'dot-explorer.',
            file=sys.stderr,
        )
        return 1

    # app.py imports its siblings as top-level `core.*`, so its own directory
    # has to come first on the path -- the same thing `shiny run app.py` and
    # `shinylive export` do.  Prepending also means the app's `core` wins over
    # any unrelated module of that name.
    sys.path.insert(0, str(APP_DIR))

    # Import the App object and hand it to run_app directly rather than passing
    # a path or a "module:attr" string: run_app accepts `str | shiny.App`, and
    # a str would be re-resolved against the working directory.
    import app as app_module  # noqa: PLC0415 - needs the sys.path entry above

    run_app(
        app_module.app,
        host=args.host,
        port=args.port,
        launch_browser=not args.no_browser,
    )
    return 0


if __name__ == '__main__':
    sys.exit(main())
