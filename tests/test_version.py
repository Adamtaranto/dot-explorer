"""Version plumbing: __version__ export and the stamp script's tag mapping."""

import importlib.util
from pathlib import Path
import re
import sys

import dot_explorer

REPO_ROOT = Path(__file__).resolve().parent.parent


def _load_stamp_module():
    spec = importlib.util.spec_from_file_location(
        'stamp_version', REPO_ROOT / 'scripts' / 'stamp_version.py'
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules['stamp_version'] = mod
    spec.loader.exec_module(mod)
    return mod


def test_dunder_version_is_pep440_like():
    """A release reports X.Y.Z; a dev build appends +<hash>[.dirty]."""
    assert re.match(
        r'^\d+\.\d+\.\d+(\+[0-9a-f]+(\.dirty)?)?$', dot_explorer.__version__
    )


def test_version_module_importable():
    from dot_explorer import _version

    assert _version.__version__ == dot_explorer.__version__


def test_parse_describe_exact_tag():
    stamp = _load_stamp_module()
    assert stamp.parse_describe('v1.2.3') == ('1.2.3', True)


def test_parse_describe_after_tag_keeps_the_released_version():
    """No .postN suffix: it would disagree with the wheel's own metadata."""
    stamp = _load_stamp_module()
    assert stamp.parse_describe('v1.2.3-5-gabc1234') == ('1.2.3', False)


def test_parse_describe_without_a_tag_returns_none():
    """No tag means no known version — never a guessed default.

    Shallow CI checkouts fetch no tags; guessing there would stamp the
    committed manifest backwards to whatever the default happened to be.
    """
    stamp = _load_stamp_module()
    assert stamp.parse_describe(None) == (None, False)
    assert stamp.parse_describe('not-a-tag') == (None, False)


def test_compute_updates_is_a_noop_without_tags(monkeypatch):
    stamp = _load_stamp_module()
    monkeypatch.setattr(stamp, 'git_describe', lambda repo=None: None)
    assert stamp.compute_updates() == {}


def test_refresh_local_version_is_a_noop_without_tags(tmp_path, monkeypatch):
    stamp = _load_stamp_module()
    monkeypatch.setattr(stamp, 'git_describe', lambda repo=None: None)
    assert stamp.refresh_local_version(tmp_path) is None


def test_version_files_untouched_between_tags(tmp_path, monkeypatch):
    """An ordinary commit must not dirty _version.py or CITATION.cff."""
    stamp = _load_stamp_module()
    monkeypatch.setattr(stamp, 'git_describe', lambda repo=None: 'v1.2.3-5-gabc1234')
    targets = stamp.compute_updates()
    names = {p.name for p in targets}
    assert '_version.py' not in names
    assert 'CITATION.cff' not in names
    assert 'Cargo.toml' in names


def test_version_files_stamped_at_a_tag(monkeypatch):
    stamp = _load_stamp_module()
    monkeypatch.setattr(stamp, 'git_describe', lambda repo=None: 'v1.2.3')
    monkeypatch.setattr(stamp, 'tag_date', lambda version, repo=None: '2026-01-02')
    names = {p.name for p in stamp.compute_updates()}
    assert {'_version.py', 'CITATION.cff', 'Cargo.toml'} <= names


def test_stamp_cargo_toml_only_touches_package_section():
    stamp = _load_stamp_module()
    text = (
        '[package]\nname = "x"\nversion = "0.0.1"\n\n'
        '[dependencies]\npyo3 = { version = "0.29" }\n'
    )
    out = stamp.stamp_cargo_toml(text, '9.9.9', Path('Cargo.toml'))
    assert 'version = "9.9.9"' in out
    assert 'pyo3 = { version = "0.29" }' in out


def test_stamp_cargo_lock_targets_own_package():
    stamp = _load_stamp_module()
    text = (
        '[[package]]\nname = "ahash"\nversion = "0.8.1"\n\n'
        '[[package]]\nname = "dot-explorer"\nversion = "0.0.1"\n'
    )
    out = stamp.stamp_cargo_lock(text, '9.9.9', Path('Cargo.lock'))
    assert 'name = "ahash"\nversion = "0.8.1"' in out
    assert 'name = "dot-explorer"\nversion = "9.9.9"' in out


def test_compute_updates_is_idempotent_on_repo():
    stamp = _load_stamp_module()
    for path, content in stamp.compute_updates().items():
        assert path.read_text() == content, f'{path} is stale — run stamp_version.py'


def test_local_version_composition(monkeypatch):
    """Dev builds carry a PEP 440 local version naming the commit."""
    stamp = _load_stamp_module()
    monkeypatch.setattr(stamp, 'short_hash', lambda repo=None: 'abc1234')
    monkeypatch.setattr(stamp, 'is_dirty', lambda repo=None: False)
    assert stamp.local_version('0.1.0') == '0.1.0+abc1234'
    monkeypatch.setattr(stamp, 'is_dirty', lambda repo=None: True)
    assert stamp.local_version('0.1.0') == '0.1.0+abc1234.dirty'


def test_local_version_without_git(monkeypatch):
    stamp = _load_stamp_module()
    monkeypatch.setattr(stamp, 'short_hash', lambda repo=None: None)
    assert stamp.local_version('0.1.0') is None


def test_local_versions_are_pep440(monkeypatch):
    stamp = _load_stamp_module()
    monkeypatch.setattr(stamp, 'short_hash', lambda repo=None: 'abc1234')
    monkeypatch.setattr(stamp, 'is_dirty', lambda repo=None: True)
    pep440 = re.compile(
        r'^([0-9]+!)?[0-9]+(\.[0-9]+)*((a|b|rc)[0-9]+)?(\.post[0-9]+)?'
        r'(\.dev[0-9]+)?(\+[a-zA-Z0-9]+([._-][a-zA-Z0-9]+)*)?$'
    )
    assert pep440.match(stamp.local_version('0.1.0'))


def test_refresh_local_version_writes_between_tags(tmp_path, monkeypatch):
    stamp = _load_stamp_module()
    pkg = tmp_path / 'python' / 'dot_explorer'
    pkg.mkdir(parents=True)
    monkeypatch.setattr(stamp, 'git_describe', lambda repo=None: 'v1.2.3-4-gabc1234')
    monkeypatch.setattr(stamp, 'short_hash', lambda repo=None: 'abc1234')
    monkeypatch.setattr(stamp, 'is_dirty', lambda repo=None: False)
    assert stamp.refresh_local_version(tmp_path) == '1.2.3+abc1234'
    assert "__version__ = '1.2.3+abc1234'" in (pkg / '_version_local.py').read_text()


def test_refresh_local_version_removes_the_override_at_a_tag(tmp_path, monkeypatch):
    """A stale dev override must never ride along into a release wheel."""
    stamp = _load_stamp_module()
    pkg = tmp_path / 'python' / 'dot_explorer'
    pkg.mkdir(parents=True)
    stale = pkg / '_version_local.py'
    stale.write_text("__version__ = '1.2.2+deadbee'\n")
    monkeypatch.setattr(stamp, 'git_describe', lambda repo=None: 'v1.2.3')
    assert stamp.refresh_local_version(tmp_path) is None
    assert not stale.exists()


def test_version_module_prefers_the_local_override():
    """_version.py must import the untracked override when it exists."""
    src = (REPO_ROOT / 'python' / 'dot_explorer' / '_version.py').read_text()
    assert 'from dot_explorer._version_local import __version__' in src
    assert 'except ImportError' in src
