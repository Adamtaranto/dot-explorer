"""Tests for the dot-explorer-app console-script launcher."""

import builtins
from pathlib import Path

import pytest

from dot_explorer import launch


def test_app_dir_points_at_packaged_app():
    """The launcher resolves the app shipped alongside the package."""
    assert launch.APP_DIR.name == 'app'
    assert launch.APP_DIR.parent.name == 'dot_explorer'
    assert (launch.APP_DIR / 'app.py').is_file()
    assert (launch.APP_DIR / 'core' / 'align.py').is_file()
    assert (launch.APP_DIR / 'www' / 'bridge.js').is_file()


def test_packaged_app_carries_no_wheels_dir():
    """A .whl inside the package would end up nested in the built wheel.

    The Shinylive staging copy gets the wasm wheel instead
    (scripts/build_shinylive.py).
    """
    assert not (launch.APP_DIR / 'wheels').exists()


def test_defaults():
    args = launch._parse_args([])
    assert args.host == '127.0.0.1'
    assert args.port == 8000
    assert args.no_browser is False


def test_arguments_parsed():
    args = launch._parse_args(['--host', '0.0.0.0', '--port', '9999', '--no-browser'])
    assert args.host == '0.0.0.0'
    assert args.port == 9999
    assert args.no_browser is True


def test_missing_shiny_exits_with_install_hint(monkeypatch, capsys):
    """Shiny is optional, so the failure must name the extra to install."""
    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == 'shiny':
            raise ImportError('No module named shiny')
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, '__import__', fake_import)
    assert launch.main([]) == 1
    err = capsys.readouterr().err
    assert 'dot-explorer[app]' in err


def test_version_flag(capsys):
    import dot_explorer

    assert launch.main(['--version']) == 0
    assert capsys.readouterr().out.strip() == dot_explorer.__version__


def test_missing_app_file_reports_broken_install(monkeypatch, capsys, tmp_path):
    pytest.importorskip('shiny')
    monkeypatch.setattr(launch, 'APP_DIR', Path(tmp_path) / 'absent')
    assert launch.main([]) == 1
    assert 'reinstalling' in capsys.readouterr().err
