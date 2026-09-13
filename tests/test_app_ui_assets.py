"""Contract tests for the app's static assets (www/*.js, www/*.css).

The JS side has no test harness, so these pin behaviour at the source level:
a regression that re-introduces a removed hook fails here rather than in a
browser nobody is watching.
"""

from pathlib import Path

APP_DIR = Path(__file__).parent.parent / 'python' / 'dot_explorer' / 'app'
FULLSCREEN_JS = (APP_DIR / 'www' / 'fullscreen.js').read_text()
APP_CSS = (APP_DIR / 'www' / 'app.css').read_text()
ALIGNERS_JS = (APP_DIR / 'www' / 'aligners.js').read_text()


def test_fullscreen_keeps_the_active_tab():
    """Entering fullscreen from the Heatmap tab used to snap back to Plot.

    The tab strip is no longer hidden in fullscreen, so there is nothing to
    strand the user on and no reason to force the first tab.
    """
    assert 'activatePlotTab' not in FULLSCREEN_JS
    hidden = APP_CSS.split('html.de-fullscreen .de-nav-hint')[1].split('}')[0]
    assert '.nav-tabs' not in hidden
    assert 'html.de-fullscreen .de-plot-area .de-heatmap-wrap' in APP_CSS


def test_aligner_inputs_are_unlinked_before_every_remount():
    """Aioli keeps files already present at a mounted name.

    Running Q vs T and then Q vs Q on the same cached CLI instance mounted
    the new pair under the same fixed names, so the stale target.fa won and
    the "self" run was still Q vs T.  Both inputs are removed before each
    mount so the new data always lands.
    """
    body = ALIGNERS_JS.split('async function runAligner(')[1]
    before_mount = body.split('await CLI.mount(')[0]
    assert 'var stale = [QUERY_FILENAME, TARGET_FILENAME];' in before_mount
    assert 'await CLI.fs.unlink(stale[i]);' in before_mount
    assert 'mountedOn[msg.tool] = null;' in before_mount
