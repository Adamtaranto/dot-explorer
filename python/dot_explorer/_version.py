"""Version stamped from the latest git tag by scripts/stamp_version.py.

Do not edit by hand: tag a release (``git tag vX.Y.Z``) and the
pre-commit hook / CI stamp step rewrites this file.

Builds made between releases additionally write an untracked
``_version_local.py`` carrying a PEP 440 local version
(``X.Y.Z+<short hash>``); when present it wins, so a development
build identifies the commit it came from without this file — or
Cargo.toml — changing on every commit.
"""

__version__ = '0.1.0'

try:  # pragma: no cover - only present in development builds
    from dot_explorer._version_local import __version__  # noqa: F401,F811
except ImportError:
    pass
