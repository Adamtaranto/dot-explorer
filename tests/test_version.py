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
    assert re.match(r'^\d+\.\d+\.\d+(\.post\d+)?$', dot_explorer.__version__)


def test_version_module_importable():
    from dot_explorer import _version

    assert _version.__version__ == dot_explorer.__version__


def test_parse_describe_exact_tag():
    stamp = _load_stamp_module()
    assert stamp.parse_describe('v1.2.3') == ('1.2.3', '1.2.3', True)


def test_parse_describe_post_tag():
    stamp = _load_stamp_module()
    assert stamp.parse_describe('v1.2.3-5-gabc1234') == ('1.2.3.post5', '1.2.3', False)


def test_parse_describe_no_tags_falls_back():
    stamp = _load_stamp_module()
    py, cargo, exact = stamp.parse_describe(None)
    assert py == cargo == stamp.FALLBACK_VERSION
    assert exact is False


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
