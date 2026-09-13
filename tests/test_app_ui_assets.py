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
BRIDGE_JS = (APP_DIR / 'www' / 'bridge.js').read_text()
APP_PY = (APP_DIR / 'app.py').read_text()
HTML_DIR = APP_DIR.parent / '_html'
REPORT_JS = (HTML_DIR / 'report.js').read_text()
TEMPLATE = (HTML_DIR / 'template.html').read_text()
REPORT_CSS = (HTML_DIR / 'report.css').read_text()


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


def test_report_context_menu_is_wired():
    """Right-click offers shadow selection, bulk download and flips."""
    assert '<div id="de-ctx" hidden></div>' in TEMPLATE
    assert '#de-ctx {' in REPORT_CSS
    assert "svg.addEventListener('contextmenu'" in REPORT_JS
    for label in (
        'Select all alignments in shadow of selection',
        'Download selected alignments (FASTA)',
        'Select alignments in shadow of feature',
    ):
        assert label in REPORT_JS
    for msg in ('de-flip-query', 'de-download-selected'):
        assert msg in REPORT_JS
        assert msg in BRIDGE_JS
    # Shadow tests run in data space with a binary search per axis.
    assert 'function inRanges(' in REPORT_JS
    assert 'function selectShadow(' in REPORT_JS
    # Feature shadows undo the per-axis mirroring before comparing.
    body = REPORT_JS.split('function selectFeatureShadow(')[1].split('\n  }')[0]
    assert 'panel.reverse_query' in body and 'panel.reverse_target' in body


def test_bridge_forwards_flips_selection_and_download():
    assert "setInputValue(\n          'flip_query'" in BRIDGE_JS
    assert "setInputValue(\n          'download_selected'" in BRIDGE_JS
    assert "setInputValue('match_selection', { matches: lastMatchSel })" in BRIDGE_JS
    assert "addCustomMessageHandler('rd_click_download'" in BRIDGE_JS
    assert 'input.flip_query' in APP_PY
    assert 'input.download_selected' in APP_PY
    assert 'input.match_selection' in APP_PY
    assert "'rd_click_download', {'id': 'dl_selected_fasta'}" in APP_PY


def test_per_side_copy_flags_reach_the_report():
    """A PAF run may carry one assembly: the copy buttons are per side."""
    assert "reply['copy_query'] = have_q" in APP_PY
    assert "reply['copy_target'] = have_t" in APP_PY
    handler = REPORT_JS.split("msg.type === 'de-seq-response'")[1].split(
        'de-copy-response'
    )[0]
    assert 'msg.copy_query' in handler and 'msg.copy_target' in handler


def test_sequence_lookup_never_crosses_roles():
    """Two assemblies can share a contig name: no cross-role fallback."""
    body = APP_PY.split('def _sequence_for(')[1].split('\n    def ')[0]
    assert 'provider = meta.get(side)' in body
    assert "('query', 'target')" not in body


def test_paf_mode_offers_supplementary_uploads():
    for input_id in (
        'paf_seq_format',
        'paf_query_fasta',
        'paf_target_fasta',
        'paf_query_gbk',
        'paf_target_gbk',
    ):
        assert f"'{input_id}'" in APP_PY
    body = APP_PY.split('def _parse_paf_sequences(')[1].split('\n    @')[0]
    assert 'validate_paf_names(' in body
    assert 'raise ValueError' in body


def test_layout_mirrors_targets_only_in_self_mode_and_applies_flips():
    body = APP_PY.split('def layout():')[1].split('\n    # ---')[0]
    assert 'apply_manual_flips(reverse, manual_flips())' in body
    assert "'reverse_targets': set(reverse) if self_mode() else set()" in body
    figure = APP_PY.split('def make_figure(')[1].split('\n    def ')[0]
    assert "kwargs['reverse_targets']" in figure
    # Every sequence export reads the same orientation as the plot.
    for fn in ('def dl_fasta(', 'def dl_cluster_fasta('):
        export = APP_PY.split(fn)[1].split('\n    @')[0]
        assert "lay['reverse']" in export
