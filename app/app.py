"""dot-explorer assembly comparison — browser app (Shiny for Python / Shinylive).

Runs entirely client-side: uploads never leave the browser.  Under Pyodide
the dot-explorer wasm wheel bundled in ``wheels/`` is installed at startup;
run natively (``shiny run app/app.py``) it uses the installed dot-explorer.
"""

from __future__ import annotations

import importlib.util
import io
import logging
import os
from pathlib import Path
import re
import sys
import tempfile
import time
import uuid

from core.align import (
    AVAILABLE_METHODS,
    BIOWASM_TOOLS,
    METHOD_LABELS,
    MINIMAP2_PRESET_DEFAULTS,
    MINIMAP2_PRESETS,
    alignment_from_tool_output,
    build_tool_args,
    fasta_text,
    paf_alignment_from_text,
    paf_text_from_alignment,
)
from core.annotation_colors import (
    assign_shared_colors,
    color_map_for,
    display_name,
    normalise_type,
)
from core.annotation_state import (
    ANNOTATION_ROLES,
    apply_annotation_config,
    apply_feature_overrides,
    build_feature_rows,
    count_pending_overrides,
    count_pending_type_changes,
    dedupe_feature_rows,
    merge_annotations,
    replace_source,
    type_slug_map,
)
from core.cache import QUERY_GROUP, TARGET_GROUP, SessionCache
from core.cluster import (
    ProviderIndex,
    alignment_coverage_matrix,
    cluster_deps_missing,
    cluster_table_rows,
    coverage_align_params,
    is_clean_minimap2,
    tree_layout_order,
)
from core.export import reordered_fasta_text
from core.fasta import content_digest
from core.genbank import parse_genbank_bytes
from core.panels import (
    filter_by_min_length,
    has_self_pair,
    nav_tips,
    panel_pair,
    resolve_orders,
)
from core.seqs import (
    SequenceProvider,
    provider_from_fasta_input,
    provider_from_path,
)
from core.state import CAP_STYLE_CHOICES, ORDER_CHOICES, PlotConfig, svg_linecap
from core.validate import validate_annotation_names, validate_query_names
from core.wheels import pick_wasm_wheel, runtime_platform_tag
import matplotlib  # noqa: F401  (ensures shinylive bundles the pyodide package)
import numpy  # noqa: F401
from shiny import App, reactive, render, req, ui

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('dot_explorer_app')

APP_DIR = Path(__file__).parent

_boot_done = False


async def ensure_dot_explorer() -> None:
    """Make :mod:`dot_explorer` importable, installing the wasm wheel if needed.

    Under Pyodide (Shinylive) the wheel bundled in the app's ``wheels/``
    directory is installed with micropip on first call.  Natively this is a
    no-op when dot-explorer is already installed.

    Raises
    ------
    RuntimeError
        If dot-explorer is unavailable and cannot be installed (no bundled
        wheel under Pyodide, or not pip-installed natively).
    """
    global _boot_done
    if _boot_done or importlib.util.find_spec('dot_explorer') is not None:
        _boot_done = True
        return
    if sys.platform == 'emscripten':
        import micropip  # noqa: PLC0415 - pyodide-only module

        wheel = pick_wasm_wheel(
            list(APP_DIR.glob('wheels/*.whl')), runtime_platform_tag()
        )
        logger.info('Installing bundled wheel %s', wheel.name)
        await micropip.install(f'emfs:{wheel}')
        # find_spec() above already cached a FileFinder for site-packages
        # that predates the install, so without this the very next
        # `import dot_explorer` can still raise ModuleNotFoundError.
        importlib.invalidate_caches()
        _boot_done = True
        return
    raise RuntimeError(
        'dot-explorer is not installed — run `pip install dot-explorer` '
        '(or `maturin develop` from the repo).'
    )


def _method_choices() -> dict[str, str]:
    """Return UI choices for the method selector, marking unimplemented ones.

    Returns
    -------
    dict[str, str]
        Method key -> label, with unavailable methods suffixed.
    """
    return {
        key: (label if key in AVAILABLE_METHODS else f'{label} — coming soon')
        for key, label in METHOD_LABELS.items()
        # PAF import is an input *mode* (the input_mode radio), not an
        # alignment method — with a PAF there is nothing left to align.
        if key != 'paf'
    }


def _lbl(text: str, tip: str):
    """Input label with a hoverable info icon explaining the setting.

    Parameters
    ----------
    text : str
        The visible label text.
    tip : str
        One-line explanation shown in a tooltip on the icon.

    Returns
    -------
    htmltools.Tag
        Label span suitable as any input's ``label`` argument.
    """
    return ui.span(
        text,
        ' ',
        ui.tooltip(ui.tags.span('ⓘ', class_='de-info'), tip, placement='right'),
    )


# --- W2: interactive plot ---------------------------------------------------
# Script injected into the generated HTML report before embedding: posts a
# message to the parent app window when a dotplot panel is double-clicked.
# app/www/bridge.js listens for it and forwards to the 'panel_dblclick' input.
_PANEL_DBLCLICK_JS = """
<script>
(function () {
  'use strict';
  // Tell report.js that double-click drills down here, so it defers the
  // single-click focus zoom briefly and a double-click cancels it —
  // without this, double-clicking a panel zooms in and then swaps views.
  window.RD_DBLCLICK_DRILLDOWN = true;
  document.addEventListener('dblclick', function (ev) {
    // Walk ancestors testing the exact panel id.  closest() with a prefix
    // selector would stop at the axes background group, whose id merely
    // starts the same way, and the drill-down would silently never fire.
    var node = ev.target;
    var m = null;
    while (node && node.nodeType === 1) {
      m = /^de-panel-(\\d+)-(\\d+)$/.exec(node.id || '');
      if (m) { break; }
      node = node.parentNode;
    }
    if (!m) { return; }
    window.parent.postMessage(
      {type: 'de-panel-dblclick',
       row: parseInt(m[1], 10),
       col: parseInt(m[2], 10)},
      '*'
    );
  });
})();
</script>
"""


#: Columns of the drill-down annotations table, with how each one sorts.
#: 'check' reads the checkbox, 'color' the picker's value, 'num' strips the
#: thousands separators the cells are formatted with, 'text' compares
#: case-insensitively.
_FEATURE_COLUMNS: tuple[tuple[str, str], ...] = (
    ('Show', 'check'),
    ('Colour', 'color'),
    ('Axis', 'text'),
    ('Contig', 'text'),
    ('Type', 'text'),
    ('Name', 'text'),
    ('Start', 'num'),
    ('End', 'num'),
    ('Length', 'num'),
    ('Strand', 'text'),
    ('Source', 'text'),
    ('Attributes', 'text'),
)


def inject_panel_bridge(html: str) -> str:
    """Insert the panel double-click bridge script into a report document.

    Parameters
    ----------
    html : str
        Full HTML report text produced by ``DotPlotter.to_html``.

    Returns
    -------
    str
        The report with the bridge script spliced in before ``</body>``
        (appended at the end if no closing body tag is found).
    """
    idx = html.rfind('</body>')
    if idx == -1:
        return html + _PANEL_DBLCLICK_JS
    return html[:idx] + _PANEL_DBLCLICK_JS + html[idx:]


def debounce(delay_secs: float):
    """Debounce a reactive calculation (standard Shiny recipe).

    Wraps a ``reactive.calc``-decorated function so downstream consumers
    only see a new value once the source has been stable for *delay_secs*.
    Used for numeric inputs whose spinner arrows would otherwise trigger a
    full plot re-render for every 0.1 step.

    Parameters
    ----------
    delay_secs : float
        Quiet period before the new value propagates.

    Returns
    -------
    Callable
        Decorator for a ``reactive.calc`` function.
    """

    def wrapper(f):
        when = reactive.value(None)
        trigger = reactive.value(0)

        @reactive.effect(priority=102)
        def _primer():
            try:
                f()
            except Exception:  # noqa: BLE001 - value read only to register deps
                pass
            with reactive.isolate():
                when.set(time.monotonic() + delay_secs)

        @reactive.effect(priority=101)
        def _timer():
            deadline = when()
            if deadline is None:
                return
            now = time.monotonic()
            if now >= deadline:
                when.set(None)
                with reactive.isolate():
                    trigger.set(trigger() + 1)
            else:
                reactive.invalidate_later(deadline - now)

        @reactive.calc
        @reactive.event(trigger, ignore_none=False)
        def debounced():
            with reactive.isolate():
                return f()

        return debounced

    return wrapper


# Injected into embedded reports: hide the report's own header, and fit
# the SVG to the iframe viewport.  Without the fit, a figure taller than
# the frame keeps its intrinsic size and the page scrolls — a drag-zoomed
# region then fills the whole (mostly off-screen) element and appears
# shifted up.  Sizing the SVG to the viewport (as fullscreen already did)
# letterboxes via preserveAspectRatio, so every zoom lands exactly where
# it was drawn.  Standalone to_html exports are unaffected.
_HIDE_REPORT_HEADER_CSS = (
    '<style>'
    '#de-header{display:none}'
    'body{padding-bottom:0}'
    '#de-figure{padding:0.5rem}'
    '#de-figure svg{width:100%;height:calc(100vh - 1rem)}'
    '</style>'
)

# Fullscreen-toggle icons (expand / collapse corners); www/app.css shows
# exactly one of the pair depending on the html.de-fullscreen class.
_FS_EXPAND_SVG = (
    '<svg class="de-fs-expand" width="15" height="15" viewBox="0 0 24 24"'
    ' fill="none" stroke="currentColor" stroke-width="2.4"'
    ' stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">'
    '<path d="M8 3H3v5M16 3h5v5M8 21H3v-5M16 21h5v-5"/></svg>'
)
_FS_COLLAPSE_SVG = (
    '<svg class="de-fs-collapse" width="15" height="15" viewBox="0 0 24 24"'
    ' fill="none" stroke="currentColor" stroke-width="2.4"'
    ' stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">'
    '<path d="M3 8h5V3M21 8h-5V3M3 16h5v5M21 16h-5v5"/></svg>'
)


def strip_report_header(html: str) -> str:
    """Hide the report's built-in title/navigation header for embedding.

    The standalone HTML report carries its own header with navigation hints
    (``#de-header`` in ``_html/template.html``); inside the app that
    duplicates the app-level hint bar, so the embed hides it with CSS.
    Standalone ``to_html`` exports are unaffected.

    Parameters
    ----------
    html : str
        Full HTML report text produced by ``DotPlotter.to_html``.

    Returns
    -------
    str
        The report with a header-hiding style spliced in before ``</head>``
        (prepended if no closing head tag is found).
    """
    idx = html.find('</head>')
    if idx == -1:
        return _HIDE_REPORT_HEADER_CSS + html
    return html[:idx] + _HIDE_REPORT_HEADER_CSS + html[idx:]


# --- end W2 ------------------------------------------------------------------


app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.h5('Input'),
        ui.input_radio_buttons(
            'input_mode',
            None,
            {
                'fasta': 'Assemblies (FASTA)',
                'genbank': 'Assemblies (GenBank)',
                'paf': 'Alignment (PAF)',
            },
            selected='fasta',
            inline=True,
        ),
        # Everything except the PAF branch shares the self-align checkbox,
        # the method selector and the per-method parameter panels; only the
        # upload widgets differ between FASTA and GenBank.
        ui.panel_conditional(
            "input.input_mode !== 'paf'",
            ui.panel_conditional(
                "input.input_mode === 'fasta'",
                ui.input_file(
                    'query_fasta',
                    'Query assembly (FASTA / .gz)',
                    accept=['.fa', '.fasta', '.fna', '.gz'],
                ),
            ),
            ui.panel_conditional(
                "input.input_mode === 'genbank'",
                ui.input_file(
                    'query_gbk',
                    _lbl(
                        'Query assembly (GenBank / .gz)',
                        'Sequences and annotations are read from the same '
                        'file; any GFF uploaded below is merged with them.',
                    ),
                    accept=['.gb', '.gbk', '.gbff', '.genbank', '.gz'],
                ),
            ),
            ui.input_checkbox(
                'self_align',
                _lbl(
                    'Align assembly to itself',
                    'Compare one assembly against itself to reveal repeats '
                    'and segmental duplications — no second upload needed. '
                    'Only the query assembly is used.',
                ),
                False,
            ),
            ui.panel_conditional(
                "!input.self_align && input.input_mode === 'fasta'",
                ui.input_file(
                    'target_fasta',
                    'Target / reference assembly (FASTA / .gz)',
                    accept=['.fa', '.fasta', '.fna', '.gz'],
                ),
            ),
            ui.panel_conditional(
                "!input.self_align && input.input_mode === 'genbank'",
                ui.input_file(
                    'target_gbk',
                    'Target / reference assembly (GenBank / .gz)',
                    accept=['.gb', '.gbk', '.gbff', '.genbank', '.gz'],
                ),
            ),
            ui.input_select(
                'method',
                _lbl(
                    'Alignment method',
                    'k-mer: fast exact matching, fully offline; minimap2 / '
                    'nucmer: full aligners compiled to WebAssembly (fetched '
                    'from the biowasm CDN at runtime).',
                ),
                choices=_method_choices(),
            ),
            ui.panel_conditional(
                "input.method === 'kmer'",
                ui.input_slider(
                    'k',
                    _lbl(
                        'k-mer size',
                        'Exact-match seed length; larger k gives fewer, more '
                        'specific matches.',
                    ),
                    min=8,
                    max=64,
                    value=21,
                    step=1,
                ),
                ui.input_checkbox(
                    'merge',
                    _lbl(
                        'Merge adjacent matches',
                        'Join runs of adjacent k-mer hits into longer '
                        'diagonal segments.',
                    ),
                    True,
                ),
                # Compute-time filter: repeat-rich assembly pairs can produce
                # millions of short match blocks that exhaust browser memory;
                # dropping them natively keeps real genomes workable.
                ui.input_numeric(
                    'kmer_min_block',
                    _lbl(
                        'Min match block (bp, 0 = keep all)',
                        'Drop merged match blocks shorter than this before '
                        'they reach the plot; keeps repeat-rich genome pairs '
                        'from producing millions of records.',
                    ),
                    50,
                    min=0,
                ),
            ),
            # --- W1: biowasm aligner options ---
            ui.panel_conditional(
                "input.method === 'minimap2'",
                ui.input_select(
                    'mm2_preset',
                    _lbl(
                        'Preset (-x)',
                        'Divergence preset: asm5 ≲1% sequence divergence, '
                        'asm10 ≲5%, asm20 ≲10% (safe default across '
                        'strains/isolates).',
                    ),
                    choices={p: p for p in MINIMAP2_PRESETS},
                    selected='asm20',
                ),
                # k/w/m are pre-filled with the selected preset's actual
                # values (asm20 defaults) and refreshed by _sync_mm2_defaults
                # whenever the preset changes.
                ui.input_numeric(
                    'mm2_k',
                    _lbl(
                        'k-mer size (-k)',
                        'Minimizer k-mer length; smaller values increase '
                        'sensitivity for diverged repeats.',
                    ),
                    MINIMAP2_PRESET_DEFAULTS['asm20']['k'],
                    min=1,
                    max=28,
                ),
                ui.input_numeric(
                    'mm2_w',
                    _lbl(
                        'Minimizer window (-w)',
                        'Minimizer window size; smaller windows increase seed density.',
                    ),
                    MINIMAP2_PRESET_DEFAULTS['asm20']['w'],
                    min=1,
                ),
                ui.input_numeric(
                    'mm2_m',
                    _lbl(
                        'Min chaining score (-m)',
                        'Discard chains scoring below this; raise (e.g. 200) '
                        'to filter weak or noisy repeat matches.',
                    ),
                    MINIMAP2_PRESET_DEFAULTS['asm20']['m'],
                    min=1,
                ),
                ui.input_checkbox(
                    'mm2_c',
                    _lbl(
                        'Base-level alignment (-c)',
                        'Compute exact CIGARs and identity tags (cg/NM/de): '
                        'enables gap-compressed identity colouring and '
                        'aligned-sequence display, but is slower on large '
                        'assemblies.',
                    ),
                    False,
                ),
                ui.input_checkbox(
                    'mm2_p',
                    _lbl(
                        'Retain all chains (-P)',
                        'Keep every chain instead of dropping secondary '
                        'matches, so each repeat copy maps — larger output.',
                    ),
                    False,
                ),
                ui.panel_conditional(
                    'input.self_align',
                    ui.input_checkbox(
                        'mm2_d',
                        _lbl(
                            'Skip self-diagonal matches (-D)',
                            'Drop the trivial full-length match of each '
                            'contig against itself on the main diagonal.',
                        ),
                        False,
                    ),
                ),
            ),
            ui.panel_conditional(
                "input.method === 'nucmer'",
                # mummer's -l 20 -c 65 defaults are tuned for small regions;
                # for assembly-scale dotplots they mostly add noise clusters
                # and minutes of runtime.  -l 100 -c 200 is the conventional
                # whole-genome comparison setting.
                ui.input_numeric(
                    'nucmer_l',
                    _lbl(
                        'Min match length (-l)',
                        'Minimum exact-match anchor length; larger values '
                        'run faster with fewer spurious matches.',
                    ),
                    100,
                    min=1,
                ),
                ui.input_numeric(
                    'nucmer_c',
                    _lbl(
                        'Min cluster length (-c)',
                        'Minimum combined anchor length for a cluster to be '
                        'reported as an alignment.',
                    ),
                    200,
                    min=1,
                ),
                ui.input_checkbox(
                    'nucmer_maxmatch',
                    _lbl(
                        'Use all matches (--maxmatch)',
                        'Use every anchor match regardless of uniqueness — '
                        'required to see all copies of repetitive regions.',
                    ),
                    False,
                ),
                ui.input_checkbox(
                    'nucmer_nosimplify',
                    _lbl(
                        'Keep repeat-induced alignments (--nosimplify)',
                        'Keep shadowed clusters instead of simplifying them '
                        'away — needed to find inexact repeats in '
                        'self-alignments.',
                    ),
                    False,
                ),
            ),
            # --- end W1 ---
        ),
        ui.panel_conditional(
            "input.input_mode === 'paf'",
            ui.input_file('paf_file', 'PAF file', accept=['.paf', '.txt']),
            ui.input_file(
                'paf_query_fasta',
                'Query assembly (optional — enables the reordered-FASTA download)',
                accept=['.fa', '.fasta', '.fna', '.gz'],
            ),
        ),
        ui.input_action_button('run', 'Run comparison', class_='btn-primary'),
        ui.hr(),
        ui.h5('Plot options'),
        # --- W2: interactive plot ---
        ui.input_checkbox(
            'interactive',
            _lbl(
                'Interactive plot (zoom & drill-down)',
                'Render a zoomable HTML report with clickable matches '
                'instead of a static image.',
            ),
            True,
        ),
        ui.input_select(
            'contig_order',
            _lbl(
                'Contig order',
                'How contigs are arranged along the axes: upload order, '
                'longest first, or reordered to maximise colinearity with '
                'the other assembly.',
            ),
            choices=ORDER_CHOICES,
        ),
        ui.input_numeric(
            'min_contig_len',
            _lbl(
                'Min contig length (bp, 0 = keep all)',
                'Leave contigs shorter than this out of the plot. One panel '
                'per contig means a few chromosomes can be buried under '
                'hundreds of short scaffolds. Excluded contigs are still '
                'written to the reordered-FASTA download.',
            ),
            value=0,
            min=0,
            step=1000,
        ),
        ui.input_checkbox(
            'auto_reverse',
            _lbl(
                'Auto-flip reversed contigs',
                'Display contigs that align mostly in reverse on their '
                'reverse strand so synteny reads as a forward diagonal.',
            ),
            False,
        ),
        ui.input_checkbox(
            'hide_internal_axes',
            _lbl(
                'Hide internal axes',
                'Remove the internal panel borders and ticks so the '
                'all-vs-all grid reads as one continuous plot.',
            ),
            False,
        ),
        # Identity colouring only makes sense for tool/PAF alignments (k-mer
        # matches are always 100% identity); output.result_kind is a hidden
        # text output that tracks the current result.
        ui.panel_conditional(
            "output.result_kind === 'paf'",
            ui.input_checkbox(
                'color_by_identity',
                _lbl(
                    'Colour by % identity',
                    'Colour each match by its percent identity instead of '
                    'forward/reverse strand (aligner and PAF results only).',
                ),
                False,
            ),
            ui.panel_conditional(
                'input.color_by_identity',
                ui.input_select(
                    'identity_palette',
                    _lbl(
                        'Identity palette',
                        'Colour map used for the identity scale.',
                    ),
                    choices=['viridis', 'plasma', 'cividis', 'coolwarm'],
                ),
            ),
        ),
        ui.input_numeric(
            'min_length',
            _lbl(
                'Min match length (bp)',
                'Hide matches shorter than this many base pairs — applied '
                'instantly, without re-rendering.',
            ),
            0,
            min=0,
        ),
        ui.input_numeric(
            'dot_size',
            _lbl(
                'Line width',
                'Stroke width of the match segments — applied instantly, '
                'without re-rendering.',
            ),
            0.5,
            min=0.1,
            max=5,
            step=0.1,
        ),
        ui.input_select(
            'cap_style',
            _lbl(
                'Line cap',
                'Shape of the match-segment ends. Square and round keep a '
                'match shorter than the line width sitting on its own '
                'diagonal; flat draws it square-on, so it looks rotated. '
                'Applied instantly, without re-rendering.',
            ),
            choices=CAP_STYLE_CHOICES,
            selected='projecting',
        ),
        ui.hr(),
        # --- Trees & clustering ----------------------------------------------
        # Self-alignment only: one tree can only order one shared axis, so a
        # cross-assembly comparison (different row and column sets) has no
        # coherent tree order.  The whole section is hidden otherwise; the
        # inputs stay bound (panel_conditional toggles CSS display only).
        ui.panel_conditional(
            "input.input_mode !== 'paf' && input.self_align",
            ui.h5('Trees & clustering'),
            ui.input_file(
                'tree_file',
                _lbl(
                    'Tree (newick / IQ-TREE .treefile)',
                    'Order the contigs by this tree and draw it beside the '
                    'plot. Tip labels must match the sequence names exactly; '
                    'while a tree is active it overrides the contig-order '
                    'and auto-flip options.',
                ),
                accept=['.nwk', '.newick', '.treefile', '.txt'],
            ),
            ui.input_checkbox(
                'cluster_enabled',
                _lbl(
                    'Compute clustering tree',
                    'Sketch every contig with sourmash, compare all pairs, '
                    'and build a hierarchical clustering tree that orders '
                    'the plot (an uploaded tree takes precedence). Adds '
                    'cluster assignment, a heatmap tab, and CSV exports.',
                ),
                False,
            ),
            ui.panel_conditional(
                'input.cluster_enabled',
                ui.input_select(
                    'cluster_metric',
                    _lbl(
                        'Similarity metric',
                        'Jaccard: fraction of shared k-mers. Angular: '
                        'k-mer sharing weighted by copy number (needs '
                        'abundance). ANI: estimated nucleotide identity. '
                        'See the Similarity & Clustering docs page for '
                        'guidance.',
                    ),
                    choices={
                        'jaccard': 'Jaccard similarity',
                        'angular': 'Angular similarity (abundance)',
                        'ani': 'ANI approximation',
                    },
                ),
                ui.input_numeric(
                    'sketch_k',
                    _lbl(
                        'Sketch k-mer size',
                        'K-mer length for the sourmash sketches. 21 (the '
                        'sourmash DNA default) suits most genomes; 31 for '
                        'strain-level comparisons.',
                    ),
                    21,
                    min=4,
                    max=51,
                    step=1,
                ),
                ui.input_numeric(
                    'sketch_scaled',
                    _lbl(
                        'Sketch scaled factor',
                        'Keep roughly one hash per this many bp. Lower '
                        'values keep more hashes — better ANI estimates '
                        '(tighter confidence intervals) and short-sequence '
                        'resolution, at higher memory cost. Aim for at '
                        'least a few hundred hashes per contig.',
                    ),
                    1000,
                    min=1,
                    step=100,
                ),
                ui.input_checkbox(
                    'sketch_abund',
                    _lbl(
                        'Track k-mer abundance (-p abund)',
                        'Record how often each kept k-mer occurs. Required '
                        'by the angular metric; ignored by the others.',
                    ),
                    True,
                ),
                ui.input_radio_buttons(
                    'cluster_mode',
                    _lbl(
                        'Cluster assignment',
                        'Similarity cutoff: cut the clustering tree at one '
                        'similarity threshold. Identity + coverage: link '
                        'contigs passing BOTH an ANI threshold and a '
                        'containment (coverage) threshold, then take '
                        'connected components.',
                    ),
                    choices={
                        'similarity': 'Similarity cutoff',
                        'identity_coverage': 'Identity + coverage thresholds',
                    },
                ),
                ui.panel_conditional(
                    "input.cluster_mode === 'similarity'",
                    ui.input_slider(
                        'cluster_cutoff',
                        _lbl(
                            'Similarity cutoff',
                            'Contigs whose clustered similarity is at least '
                            'this value share a cluster.',
                        ),
                        min=0.0,
                        max=1.0,
                        value=0.8,
                        step=0.01,
                    ),
                    ui.input_checkbox(
                        'cluster_cutoff_line',
                        _lbl(
                            'Show cutoff line on tree',
                            'Draw the cutoff as a dashed line through the '
                            'dendrogram, at 1 − cutoff from the tips.',
                        ),
                        True,
                    ),
                ),
                ui.panel_conditional(
                    "input.cluster_mode === 'identity_coverage'",
                    ui.input_select(
                        'coverage_source',
                        _lbl(
                            'Coverage source',
                            'Containment: fraction of shared sourmash '
                            'k-mers (collapses when relatives differ by '
                            'scattered SNPs — at 98% identity only ~65% '
                            'of 21-mers survive). Alignment: fraction of '
                            'the contig covered by minimap2 alignment '
                            'blocks — robust to SNPs. Always computed '
                            'from a clean minimap2 run without -P (run '
                            'in the background when the current result '
                            'is nucmer, k-mer, an uploaded PAF, or '
                            'minimap2 with -P).',
                        ),
                        choices={
                            'containment': 'sourmash containment',
                            'alignment': 'Alignment block coverage',
                        },
                    ),
                    ui.input_slider(
                        'identity_cutoff',
                        _lbl(
                            'Min identity (ANI)',
                            'Estimated average nucleotide identity both '
                            'contigs must share. Lower the sketch scaled '
                            'factor for tighter ANI confidence intervals.',
                        ),
                        min=0.0,
                        max=1.0,
                        value=0.8,
                        step=0.01,
                    ),
                    ui.input_slider(
                        'coverage_cutoff',
                        _lbl(
                            'Min coverage (containment)',
                            'Fraction of a contig’s k-mers found in the '
                            'other. Note: SNPs depress exact-k-mer '
                            'containment steeply (≈ coverage × ANI^k) — '
                            'see the docs before relying on high values.',
                        ),
                        min=0.0,
                        max=1.0,
                        value=0.8,
                        step=0.01,
                    ),
                    ui.input_checkbox(
                        'cov_reciprocal',
                        _lbl(
                            'Reciprocal coverage',
                            'Require the coverage threshold in both '
                            'directions — a short fragment nested in a '
                            'longer contig then stays out of its cluster.',
                        ),
                        True,
                    ),
                ),
                ui.input_checkbox(
                    'cluster_borders_on',
                    _lbl(
                        'Outline clusters in the plot',
                        'Draw a bold border around each cluster’s block of '
                        'sub-plots (and cells in the heatmap).',
                    ),
                    True,
                ),
                ui.input_select(
                    'heatmap_cmap',
                    _lbl(
                        'Heatmap palette',
                        'Colour map for the similarity heatmap tab.',
                    ),
                    choices=[
                        'viridis',
                        'magma',
                        'plasma',
                        'cividis',
                        'coolwarm',
                        'YlGnBu',
                    ],
                ),
                ui.input_checkbox(
                    'heatmap_values',
                    _lbl(
                        'Show values in heatmap cells',
                        'Write each pairwise score inside its heatmap '
                        'cell (with the 95% CI for ANI). Readable for '
                        'small numbers of contigs only.',
                    ),
                    False,
                ),
                ui.div(
                    ui.input_action_button(
                        'apply_cluster',
                        'Apply changes',
                        class_='btn-primary btn-sm',
                    ),
                    ui.span(
                        'Setting edits are held until you apply them.',
                        class_='de-ft-apply-note',
                    ),
                    class_='de-ft-apply',
                ),
            ),
        ),
        ui.hr(),
        # --- GFF annotations -------------------------------------------------
        ui.h5('Annotations (GFF3)'),
        # Same query-then-target order as the assembly uploads above.
        ui.input_file(
            'query_gff',
            'Query annotations (.gff / .gff3 / .gz)',
            accept=['.gff', '.gff3', '.gz'],
        ),
        # Hidden alongside the target assembly when self-aligning: both axes
        # are then the query assembly, so there is no target to annotate.
        # PAF input has no self-align notion but does have both roles, so it
        # keeps the upload regardless of a stale checkbox value.
        ui.panel_conditional(
            "input.input_mode === 'paf' || !input.self_align",
            ui.input_file(
                'target_gff',
                'Target annotations (.gff / .gff3 / .gz)',
                accept=['.gff', '.gff3', '.gz'],
            ),
        ),
        ui.panel_conditional(
            "output.gff_mode === 'plain' || output.gff_mode === 'self'",
            ui.input_action_button(
                'clear_gff',
                'Clear annotations',
                class_='btn-outline-secondary btn-sm',
            ),
        ),
        # Static, so pressing Run never destroys them.  Only their visibility
        # is reactive (output.gff_mode); panel_conditional toggles CSS display
        # and leaves the inputs bound, so they can be read unconditionally.
        ui.panel_conditional(
            "output.gff_mode === 'plain' || output.gff_mode === 'self'",
            ui.input_checkbox('gff_tracks', 'Side tracks in focused pair view', True),
        ),
        ui.panel_conditional(
            "output.gff_mode === 'self'",
            ui.tooltip(
                ui.input_checkbox('gff_diagonal', 'Shade features on diagonal', True),
                'Draws each feature as a square on the diagonal of '
                'self-comparison panels.',
            ),
        ),
        ui.output_ui('gff_controls'),
        ui.hr(),
        ui.h5('Downloads'),
        ui.output_ui('downloads'),
        width=320,
    ),
    # Fixed memory note (bottom-right; hidden while the readout is empty,
    # e.g. on native runs where the wasm heap does not exist).
    ui.div(ui.output_text('app_memory'), class_='de-mem-fixed'),
    ui.div(ui.output_text('result_kind'), class_='de-hidden'),
    ui.div(ui.output_text('gff_mode'), class_='de-hidden'),
    ui.output_ui('status'),
    # --- W2: interactive plot ---
    ui.output_ui('plot_area'),
    ui.output_ui('aligner_log_ui'),
    ui.head_content(
        # Space Grotesk for the wordmark and headings, embedded as a data
        # URI so it loads offline under shinylive too.
        ui.include_css(APP_DIR / 'www' / 'font.css'),
        ui.include_css(APP_DIR / 'www' / 'app.css'),
        ui.include_js(APP_DIR / 'www' / 'bridge.js'),
        # Custom Shiny binding for native <input type="color"> pickers.
        ui.include_js(APP_DIR / 'www' / 'color-input.js', method='inline'),
        # Delegated events for the drill-down Annotations table (one
        # listener instead of ~1200 Shiny-bound inputs).
        ui.include_js(APP_DIR / 'www' / 'feature-table.js', method='inline'),
        # Row-click selection on the cluster-assignment table.
        ui.include_js(APP_DIR / 'www' / 'cluster-table.js', method='inline'),
        # Pan/zoom controls on the heatmap image.
        ui.include_js(APP_DIR / 'www' / 'heatmap-zoom.js', method='inline'),
        # Hold the sidebar's scroll position across dynamic-UI re-renders.
        ui.include_js(APP_DIR / 'www' / 'sidebar-scroll.js', method='inline'),
        # Mirror ui.Progress messages into the header's task-status slot
        # (the popup card itself is hidden in app.css).
        ui.include_js(APP_DIR / 'www' / 'task-status.js', method='inline'),
        # Plot-area fullscreen toggle (delegated; state on <html>).
        ui.include_js(APP_DIR / 'www' / 'fullscreen.js', method='inline'),
    ),
    # The busy pill and the header task status already cover "something is
    # running"; Shiny's own sliding banner on top of them reads as noise.
    ui.busy_indicators.use(pulse=False),
    # W1: biowasm aligner bridge (Aioli loaded lazily from the CDN on use).
    ui.head_content(ui.include_js(APP_DIR / 'www' / 'aligners.js', method='inline')),
    title=ui.div(
        ui.span(
            'dot',
            ui.span('·explorer', class_='de-wordmark-dot'),
            class_='de-wordmark',
        ),
        ui.span('live assembly comparison', class_='de-wordmark-sub'),
        # Filled by www/task-status.js with the active ui.Progress message.
        ui.span(id='de-task-status', class_='de-task-status'),
        # Manual light/dark toggle; defaults to the system preference and
        # stamps data-bs-theme on <html>, which the --de-* palette keys off.
        ui.span(ui.input_dark_mode(id='dark_mode'), class_='de-theme-toggle'),
        class_='de-header-flex',
    ),
    window_title='dot-explorer · assembly comparison',
    fillable=True,
)


def server(input, output, session) -> None:  # noqa: A002, D103
    cache = SessionCache()
    ready = reactive.value(False)
    boot_error = reactive.value('')
    # (kind, alignment-object, {'query': SequenceProvider|None, 'target': ...,
    #  'method': 'kmer'|'paf_upload'|'minimap2'|'nucmer', 'params': dict})
    result = reactive.value(None)

    @reactive.effect
    async def _boot():
        try:
            await ensure_dot_explorer()
            ready.set(True)
            logger.info('dot-explorer ready (platform=%s)', sys.platform)
        except RuntimeError as exc:
            boot_error.set(str(exc))

    def _parse_upload(file_input, label: str) -> SequenceProvider:
        files = file_input()
        if not files:
            raise ValueError(f'Please upload a {label} assembly.')
        return provider_from_path(Path(files[0]['datapath']))

    def _parse_seq_upload(role: str) -> SequenceProvider:
        """Parse one assembly upload for the current input mode.

        GenBank carries its annotations in the same file, so parsing also
        registers them as an annotation source for *role* — merged with
        any GFF the user uploads separately.
        """
        if input.input_mode() != 'genbank':
            fasta_input = input.query_fasta if role == 'query' else input.target_fasta
            return _parse_upload(fasta_input, role)

        gbk_input = input.query_gbk if role == 'query' else input.target_gbk
        files = gbk_input()
        if not files:
            raise ValueError(f'Please upload a {role} assembly.')
        datapath = Path(files[0]['datapath'])
        parsed = parse_genbank_bytes(datapath.read_bytes())
        # digest is over the raw upload, so an unchanged file re-run keeps the
        # existing annotation source (and the user's annotation choices).
        _set_ann_source(
            role,
            'genbank',
            files[0]['name'],
            parsed.gff_text,
            key=(parsed.fasta.digest, files[0]['name']),
        )
        # Re-serve the ORIGIN sequences through a faidx index so the parsed
        # strings can be released.
        return provider_from_fasta_input(parsed.fasta, datapath.parent)

    def _parse_inputs(progress=None) -> tuple[SequenceProvider, SequenceProvider]:
        """Parse the query (and target, or reuse query when self-aligning)."""
        if progress is not None:
            progress.set(0, message='Parsing query assembly…')
        query = _parse_seq_upload('query')
        if input.self_align():
            return query, query
        if progress is not None:
            progress.set(1, message='Parsing target assembly…')
        return query, _parse_seq_upload('target')

    # The canonical CSR k-mer index costs ~13-16 bytes/bp (both strands).
    # Measured in Chrome with Pyodide 0.27 (synthetic pairs, 1% SNPs):
    # 90 Mb combined completes with the wasm heap peaking at ~2.9 GB;
    # 100 Mb aborts the interpreter, and the abort is silent — the worker
    # restarts, all session state is lost, and no error is ever shown.
    # The guard sits at 80 Mb for headroom over that measured ceiling.
    _KMER_HARD_LIMIT = 80 * 1024 * 1024
    _KMER_WARN_LIMIT = 40 * 1024 * 1024

    def _check_kmer_memory(query: SequenceProvider, target: SequenceProvider) -> None:
        if sys.platform != 'emscripten':
            return  # native runs are bounded by system RAM, not the wasm heap
        total = query.total_length + (0 if target is query else target.total_length)
        if total > _KMER_HARD_LIMIT:
            raise ValueError(
                f'Combined assemblies are {total / 1e6:.0f} Mb — the k-mer '
                'index needs more memory than the browser allows above '
                '~80 Mb and would crash the app. Use minimap2 (or nucmer) '
                'for assemblies this size.'
            )
        if total > _KMER_WARN_LIMIT:
            ui.notification_show(
                'Large input for the k-mer method — expect a few minutes of '
                'processing; minimap2 is much faster at this scale.',
                type='warning',
                duration=10,
            )

    _TUTORIALS_URL = 'https://adamtaranto.github.io/dot-explorer/tutorials/quickstart/'

    def _combined_upload_size() -> int | None:
        """Return the combined uploaded-assembly size from file metadata.

        Raw file bytes approximate sequence length well for plain FASTA
        (gzip uploads under-count; the parse-time check in
        ``_check_kmer_memory`` stays as the exact backstop).  ``None``
        until both required uploads are present, and in PAF mode.
        """
        mode = input.input_mode()
        if mode == 'paf':
            return None
        q_in = input.query_gbk if mode == 'genbank' else input.query_fasta
        t_in = input.target_gbk if mode == 'genbank' else input.target_fasta
        q = q_in()
        if not q:
            return None
        total = int(q[0].get('size') or 0)
        if not input.self_align():
            t = t_in()
            if not t:
                return None
            total += int(t[0].get('size') or 0)
        return total

    # Upload-size gating: beyond the measured k-mer ceiling the method is
    # removed from the selector entirely (running it would silently kill
    # the Python runtime), and past the biowasm comfort zone the user is
    # pointed at precomputed-PAF input.  Both react to uploads, so the
    # guidance appears before Run is ever pressed.
    kmer_gated = reactive.value(False)
    paf_hint_shown = reactive.value(False)

    @reactive.effect
    def _gate_kmer_on_size():
        # RD_FORCE_GATE exercises the gate on native runs (for testing);
        # otherwise native is bounded by system RAM, not the wasm heap.
        if sys.platform != 'emscripten' and not os.environ.get('RD_FORCE_GATE'):
            return
        total = _combined_upload_size()
        too_big = total is not None and total > _KMER_HARD_LIMIT
        with reactive.isolate():
            changed = too_big != kmer_gated()
        if changed:
            kmer_gated.set(too_big)
            if too_big:
                choices = {k: v for k, v in _method_choices().items() if k != 'kmer'}
                with reactive.isolate():
                    selected = input.method()
                if selected == 'kmer':
                    selected = 'minimap2'
                ui.update_select('method', choices=choices, selected=selected)
                ui.notification_show(
                    ui.HTML(
                        f'Combined upload is ~{total / 1e6:.0f} MB — beyond '
                        'the ~80 Mb the in-browser k-mer index can handle, '
                        'so that method is disabled here. Use minimap2 or '
                        'nucmer, or run the dot-explorer Python library '
                        'locally — see the '
                        f'<a href="{_TUTORIALS_URL}" target="_blank" '
                        'rel="noopener">tutorial notebooks</a>.'
                    ),
                    type='warning',
                    duration=15,
                )
            else:
                with reactive.isolate():
                    selected = input.method()
                ui.update_select('method', choices=_method_choices(), selected=selected)
        if total is not None and total > _BIOWASM_SIZE_WARN:
            with reactive.isolate():
                already = paf_hint_shown()
            if not already:
                paf_hint_shown.set(True)
                ui.notification_show(
                    f'Inputs this large (~{total / 1e6:.0f} MB) can crash '
                    'the browser tab even with minimap2 / nucmer. If you '
                    'can align locally, switch the input mode to '
                    '"Alignment (PAF)" and upload the precomputed PAF '
                    'instead — plotting handles large alignments far '
                    'better than in-browser aligning does.',
                    type='warning',
                    duration=15,
                )
        elif total is not None:
            with reactive.isolate():
                if paf_hint_shown():
                    paf_hint_shown.set(False)

    @reactive.effect
    @reactive.event(input.run)
    async def _run():
        req(ready())
        mode = input.input_mode()
        try:
            if mode == 'paf':
                with ui.Progress(min=0, max=3) as progress:
                    progress.set(1, message='Parsing PAF…')
                    files = input.paf_file()
                    if not files:
                        raise ValueError('Please upload a PAF file.')
                    text = Path(files[0]['datapath']).read_text()
                    alignment = paf_alignment_from_text(text)
                    query = None
                    if input.paf_query_fasta():
                        progress.set(2, message='Parsing query assembly…')
                        query = _parse_upload(input.paf_query_fasta, 'query')
                        for warning in validate_query_names(
                            query.names,
                            alignment.query_names,
                            alignment.target_names,
                        ):
                            ui.notification_show(warning, type='warning', duration=12)
                    progress.set(3, message=f'{len(alignment)} alignment(s) loaded')
                    result.set(
                        (
                            'paf',
                            alignment,
                            {
                                'query': query,
                                'target': None,
                                'method': 'paf_upload',
                                'params': {},
                            },
                        )
                    )
                return
            method = input.method()
            if method not in AVAILABLE_METHODS:
                ui.notification_show(
                    f'{METHOD_LABELS[method]} is not available yet.',
                    type='warning',
                )
                return
            if method != 'kmer':
                return  # biowasm tools are handled by _run_biowasm
            with ui.Progress(min=0, max=4) as progress:
                query, target = _parse_inputs(progress)
                _check_kmer_memory(query, target)
                progress.set(
                    2,
                    message=(
                        f'Building k-mer index (k={input.k()}, '
                        f'{len(query.names)}×{len(target.names)} contigs)…'
                    ),
                )
                index = cache.kmer_index(
                    input.k(),
                    query,
                    target,
                    merge=input.merge(),
                    min_block_len=int(input.kmer_min_block() or 0),
                )
                progress.set(3, message='Rendering dotplot…')
                result.set(
                    (
                        'kmer',
                        index,
                        {
                            'query': query,
                            'target': target,
                            'method': 'kmer',
                            'params': {},
                        },
                    )
                )
                progress.set(4, message='Done')
        except ValueError as exc:
            ui.notification_show(str(exc), type='error', duration=8)

    # --- W1: biowasm aligners ---
    # Runs minimap2 / nucmer in a browser WebWorker via biowasm
    # (Aioli).  Python builds the CLI args and plain-FASTA payload, sends a
    # 'rd_run_aligner' custom message to www/aligners.js, and receives the
    # tool's text output back through the 'aligner_result' input.

    # biowasm tools run in their own workers with separate memory.
    # Measured in Chrome: minimap2 completes at 200 MB combined and nucmer
    # at 250 MB (Pyodide heap under 1 GB); at 300 MB the entire browser
    # tab crashed, for both tools.  Warn above 200 MB.
    _BIOWASM_SIZE_WARN = 200 * 1024 * 1024
    _NOTIF_ID = 'rd_aligner_progress'
    # Extra budget granted each time the user chooses to keep waiting.
    _TIMEOUT_EXTEND_MS = 300_000
    # In-flight request: {'request_id', 'method', 'params', 'query', 'target'}
    aligner_pending = reactive.value(None)
    # Background minimap2 run feeding alignment-based coverage clustering.
    # Coverage must come from a clean minimap2 run (no -P: secondary chains
    # inflate covered span), so when the displayed result is anything else
    # a dedicated run is dispatched. Separate from aligner_pending: it
    # never touches `result` and is not cancelled by a user Run (the JS
    # side serialises runs, so the two simply queue).
    # In-flight: {'request_id', 'query'}
    coverage_pending = reactive.value(None)
    # Completed: (query digest, PafAlignment)
    coverage_alignment = reactive.value(None)
    # Query digest whose background run failed (fall back to containment
    # instead of retrying in a loop).
    coverage_failed = reactive.value(None)
    _COV_NOTIF_ID = 'rd_coverage_progress'
    # Rolling log of completed tool runs: [{'tool', 'cmd', 'stderr', 'error'}]
    aligner_log = reactive.value([])
    _ALIGNER_LOG_MAX = 10
    # Dataset digests already shipped to aligners.js this session.  FASTA
    # payloads are sent once per upload and referenced by digest afterwards,
    # so re-runs and tool switches never re-copy whole genomes across the
    # Pyodide/JS boundary.
    sent_datasets: set[str] = set()

    def _on_coverage_result(res: dict, cov_info: dict) -> None:
        """Land the background minimap2 coverage run (never touches `result`)."""
        coverage_pending.set(None)
        ui.notification_remove(_COV_NOTIF_ID)
        _log_run(
            'minimap2',
            res.get('cmd') or '',
            res.get('stderr') or '',
            res.get('error') or None,
        )
        digest = cov_info['query'].digest
        if res.get('error'):
            coverage_failed.set(digest)
            ui.notification_show(
                'Background minimap2 coverage run failed: '
                f'{res["error"]} — clustering will use sourmash '
                'containment for coverage instead.',
                type='warning',
                duration=12,
            )
            return
        try:
            alignment = alignment_from_tool_output('minimap2', res.get('output') or '')
        except ValueError as exc:
            coverage_failed.set(digest)
            ui.notification_show(
                f'Could not parse the minimap2 coverage run output: {exc} '
                '— clustering will use sourmash containment for coverage '
                'instead.',
                type='warning',
                duration=12,
            )
            return
        cache.put_paf('minimap2', coverage_align_params(), alignment, digest, digest)
        coverage_alignment.set((digest, alignment))

    def _log_run(tool: str, cmd: str, stderr: str, error: str | None) -> None:
        entries = list(aligner_log())
        entries.append({'tool': tool, 'cmd': cmd, 'stderr': stderr, 'error': error})
        aligner_log.set(entries[-_ALIGNER_LOG_MAX:])

    async def _cancel_pending(reason: str) -> None:
        """Cancel the in-flight aligner run (JS drops/terminates it)."""
        info = aligner_pending()
        if info is None:
            return
        aligner_pending.set(None)
        ui.modal_remove()
        ui.notification_remove(_NOTIF_ID)
        await session.send_custom_message(
            'rd_cancel_aligner', {'request_id': info['request_id']}
        )
        ui.notification_show(
            f'{METHOD_LABELS[info["method"]]} run cancelled ({reason}).',
            type='message',
            duration=5,
        )

    @reactive.effect
    @reactive.event(input.method)
    async def _cancel_on_method_change():
        info = aligner_pending()
        if info is not None and info['method'] != input.method():
            await _cancel_pending('method changed')

    @reactive.effect
    @reactive.event(input.mm2_preset, ignore_init=True)
    def _sync_mm2_defaults():
        # Selecting a preset means adopting its parameters: refresh the
        # k/w/m inputs to the preset's actual values so what is shown is
        # always what runs.
        defaults = MINIMAP2_PRESET_DEFAULTS.get(input.mm2_preset())
        if defaults is None:
            return
        ui.update_numeric('mm2_k', value=defaults['k'])
        ui.update_numeric('mm2_w', value=defaults['w'])
        ui.update_numeric('mm2_m', value=defaults['m'])

    def _tool_params(method: str) -> dict:
        if method == 'minimap2':
            return {
                'preset': input.mm2_preset(),
                'k': int(input.mm2_k() or 0),
                'w': int(input.mm2_w() or 0),
                'm': int(input.mm2_m() or 0),
                'c': bool(input.mm2_c()),
                'P': bool(input.mm2_p()),
                # -D only applies when a sequence can meet itself; storing
                # the effective value keeps the cache key honest when the
                # checkbox stays ticked but self-align is turned off.
                'D': bool(input.mm2_d()) and bool(input.self_align()),
            }
        return {
            'l': int(input.nucmer_l() or 100),
            'c': int(input.nucmer_c() or 200),
            'maxmatch': bool(input.nucmer_maxmatch()),
            'nosimplify': bool(input.nucmer_nosimplify()),
        }

    async def _send_dataset(data: SequenceProvider) -> None:
        """Ship a parsed assembly to aligners.js once, keyed by digest."""
        if data.digest in sent_datasets:
            return
        await session.send_custom_message(
            'rd_mount_fasta',
            {'dataset_id': data.digest, 'text': fasta_text(data.iter_records())},
        )
        sent_datasets.add(data.digest)

    @reactive.effect
    @reactive.event(input.run)
    async def _run_biowasm():
        req(ready())
        # Every Run click supersedes any in-flight aligner request: cancel
        # it JS-side too so the superseded wasm computation actually stops
        # instead of burning CPU ahead of the run the user wants.
        if aligner_pending() is not None:
            await _cancel_pending('superseded by a new run')
        if input.input_mode() not in ('fasta', 'genbank'):
            return
        method = input.method()
        if method not in BIOWASM_TOOLS:
            return
        try:
            query, target = _parse_inputs()
            params = _tool_params(method)
            args = build_tool_args(method, params)
        except ValueError as exc:
            ui.notification_show(str(exc), type='error', duration=8)
            return
        cached = cache.get_paf(method, params, query.digest, target.digest)
        if cached is not None:
            logger.info('%s alignment cache hit', method)
            result.set(
                (
                    'paf',
                    cached,
                    {
                        'query': query,
                        'target': target,
                        'method': method,
                        'params': params,
                    },
                )
            )
            return
        if query.total_length + target.total_length > _BIOWASM_SIZE_WARN:
            ui.notification_show(
                'Combined input exceeds ~200 MB — the in-app aligners run '
                'inside the browser tab (biowasm), and inputs this large '
                'can crash the tab outright (300 MB did in testing). '
                'Consider aligning outside the app with native minimap2 and '
                'uploading the PAF via the "Alignment (PAF)" input mode '
                'instead.',
                type='warning',
                duration=12,
            )
        request_id = uuid.uuid4().hex
        aligner_pending.set(
            {
                'request_id': request_id,
                'method': method,
                'params': params,
                'query': query,
                'target': target,
            }
        )
        ui.notification_show(
            f'Running {METHOD_LABELS[method]} in your browser (the first '
            'run downloads the tool from the biowasm CDN)…',
            id=_NOTIF_ID,
            duration=None,
        )
        try:
            await _send_dataset(query)
            await _send_dataset(target)
        except ValueError as exc:
            aligner_pending.set(None)
            ui.notification_remove(_NOTIF_ID)
            ui.notification_show(str(exc), type='error', duration=8)
            return
        await session.send_custom_message(
            'rd_run_aligner',
            {
                'tool': method,
                'args': args,
                'query_id': query.digest,
                'target_id': target.digest,
                'request_id': request_id,
            },
        )

    @reactive.effect
    @reactive.event(input.aligner_result)
    def _on_aligner_result():
        res = input.aligner_result()
        info = aligner_pending()
        if not res or res.get('cancelled'):
            return  # cancelled runs were already reported when cancelled
        cov_info = coverage_pending()
        if cov_info is not None and res.get('request_id') == cov_info['request_id']:
            _on_coverage_result(res, cov_info)
            return
        if not info or res.get('request_id') != info['request_id']:
            return  # stale or unsolicited result
        aligner_pending.set(None)
        # A result can land while the "still running?" modal is open — the
        # run finished on its own.  Take the modal down with it.
        ui.modal_remove()
        ui.notification_remove(_NOTIF_ID)
        method = info['method']
        _log_run(
            method,
            res.get('cmd') or '',
            res.get('stderr') or '',
            res.get('error') or None,
        )
        if res.get('error'):
            ui.notification_show(
                f'{METHOD_LABELS[method]} failed: {res["error"]}',
                type='error',
                duration=12,
            )
            return
        try:
            alignment = alignment_from_tool_output(method, res.get('output') or '')
        except ValueError as exc:
            ui.notification_show(
                f'Could not parse {METHOD_LABELS[method]} output: {exc}',
                type='error',
                duration=12,
            )
            return
        cache.put_paf(
            method,
            info['params'],
            alignment,
            info['query'].digest,
            info['target'].digest,
        )
        ui.notification_show(
            f'{METHOD_LABELS[method]}: {len(alignment)} alignment(s).',
            type='message',
            duration=5,
        )
        result.set(
            (
                'paf',
                alignment,
                {
                    'query': info['query'],
                    'target': info['target'],
                    'method': method,
                    'params': info['params'],
                },
            )
        )

    @reactive.effect
    @reactive.event(input.aligner_progress)
    def _on_aligner_progress():
        prog = input.aligner_progress()
        info = aligner_pending()
        if not prog or not info or prog.get('request_id') != info['request_id']:
            return  # stale or pre-warm activity, not the run in flight
        stage_msgs = {
            'loading-aioli': 'Loading the biowasm runtime…',
            'initialising-tool': (
                f'Downloading {METHOD_LABELS[info["method"]]} from the '
                'biowasm CDN (first use only)…'
            ),
            'mounting-data': 'Mounting assemblies into the tool sandbox…',
            'aligning': f'Running {METHOD_LABELS[info["method"]]}…',
            'reading-output': 'Reading alignment output…',
        }
        msg = stage_msgs.get(prog.get('stage'))
        if msg:
            ui.notification_show(msg, id=_NOTIF_ID, duration=None)

    @reactive.effect
    @reactive.event(input.aligner_timeout)
    def _on_aligner_timeout():
        """Offer to keep waiting when a run outlives its watchdog budget.

        The tool has not failed and has not been stopped — aligners.js only
        notifies — so the choice really is "wait longer" vs "give up", and
        it can be offered again each time the extended budget expires.
        """
        ev = input.aligner_timeout()
        info = aligner_pending()
        if not ev or not info or ev.get('request_id') != info['request_id']:
            return  # stale: this run already finished or was superseded
        minutes = max(1, round((ev.get('total_ms') or 0) / 60000))
        label = METHOD_LABELS[info['method']]
        ui.modal_show(
            ui.modal(
                ui.p(
                    f'{label} has been running for about {minutes} minute'
                    f'{"s" if minutes != 1 else ""} and has not finished yet.'
                ),
                ui.p(
                    'It is still working in the background — nothing has been '
                    'lost. Large genomes (especially minimap2 with base-level '
                    'alignment) can take a while in the browser.',
                    class_='text-muted',
                ),
                title='Still aligning',
                easy_close=False,
                footer=ui.TagList(
                    ui.input_action_button(
                        'aligner_give_up', 'Cancel run', class_='btn-outline-secondary'
                    ),
                    ui.input_action_button(
                        'aligner_wait',
                        f'Wait another {_TIMEOUT_EXTEND_MS // 60000} minutes',
                        class_='btn-primary',
                    ),
                ),
            )
        )

    @reactive.effect
    @reactive.event(input.aligner_wait, ignore_init=True)
    async def _on_aligner_wait():
        info = aligner_pending()
        ui.modal_remove()
        if info is None:
            return  # finished while the modal was up
        await session.send_custom_message(
            'rd_extend_aligner',
            {'request_id': info['request_id'], 'ms': _TIMEOUT_EXTEND_MS},
        )
        ui.notification_show(
            f'Still running {METHOD_LABELS[info["method"]]} — waiting another '
            f'{_TIMEOUT_EXTEND_MS // 60000} minutes…',
            id=_NOTIF_ID,
            duration=None,
        )

    @reactive.effect
    @reactive.event(input.aligner_give_up, ignore_init=True)
    async def _on_aligner_give_up():
        await _cancel_pending('timed out')

    # --- end W1: biowasm aligners ---

    # Numeric spinner arrows fire one input event per 0.1 step; without a
    # quiet period every step triggers a full plot re-render.  Only the
    # settled value should reach config().
    @debounce(0.7)
    @reactive.calc
    def dot_size_settled() -> float:
        return input.dot_size() or 0.5

    @debounce(0.7)
    @reactive.calc
    def min_length_settled() -> int:
        return int(input.min_length() or 0)

    # --- GFF annotations -----------------------------------------------------
    # Annotation sources per role (upload-derived state, deliberately outside
    # PlotConfig).  A role can hold both a GenBank-derived set and a GFF
    # upload; they are merged, tagged by the file they came from.  Each entry
    # is {'kind': 'gff'|'genbank', 'filename': str, 'annotation': GffAnnotation}.
    ann_sources = {role: reactive.value(()) for role in ANNOTATION_ROLES}

    def _set_ann_source(role: str, kind: str, filename: str, ann, key=None) -> None:
        """Add, replace or clear one annotation source for *role*.

        *ann* may be a parsed GffAnnotation or GFF3 text.  Every record is
        tagged with *filename* so the drill-down can say where a feature
        came from once several sources are merged.

        *key* identifies the upload — normally ``(content_digest, filename)``.
        When it matches what is already stored this is a no-op: re-running an
        alignment re-parses the same GenBank file every time, and blindly
        re-setting the reactive value would rebuild the feature-type controls
        and wipe the user's toggles, colours and per-feature overrides.
        """
        from dot_explorer.annotation import GffAnnotation  # noqa: PLC0415

        entries = ann_sources[role]()
        # Cheap identity check first, so an unchanged GenBank upload does not
        # even pay for re-parsing its GFF text.
        if key is not None:
            current = next((e for e in entries if e['kind'] == kind), None)
            if current is not None and current.get('key') == key:
                return
        if isinstance(ann, str):
            ann = GffAnnotation.from_text(ann) if ann.strip() else None
        updated = replace_source(entries, kind, filename, ann, key)
        if updated is None:
            return  # nothing changed; leave the controls and overrides alone
        if ann is not None and len(ann):
            for rec in ann.records:
                rec.source_file = filename
        ann_sources[role].set(updated)
        _reset_feature_overrides()

    def _parse_gff_upload(file_input, role: str) -> None:
        from dot_explorer.annotation import GffAnnotation  # noqa: PLC0415

        files = file_input()
        if not files:
            _set_ann_source(role, 'gff', '', None)
            return
        raw = Path(files[0]['datapath']).read_bytes()
        try:
            ann = GffAnnotation.from_bytes(raw)
        except ValueError as exc:
            _set_ann_source(role, 'gff', '', None)
            ui.notification_show(
                f'Could not parse {role} GFF: {exc}', type='error', duration=10
            )
            return
        if len(ann) == 0:
            _set_ann_source(role, 'gff', '', None)
            ui.notification_show(
                f'No features found in the {role} GFF file.',
                type='warning',
                duration=8,
            )
            return
        _set_ann_source(
            role,
            'gff',
            files[0]['name'],
            ann,
            key=(content_digest(raw), files[0]['name']),
        )
        kinds = {e['kind'] for e in ann_sources[role]()}
        if 'genbank' in kinds:
            # Overlapping duplicates are likely, but silently deduping would
            # discard features the user deliberately supplied.  Say so and
            # let them toggle sources off instead.
            ui.notification_show(
                f'{role.capitalize()} now has annotations from both a GenBank '
                'file and a GFF; features present in both will be drawn twice.',
                type='warning',
                duration=10,
            )
        ui.notification_show(
            f'{role.capitalize()} annotations: {len(ann)} feature(s), '
            f'{len(ann.feature_types())} type(s).',
            type='message',
            duration=5,
        )

    def gff_raw_for(role: str):
        """Return the merged annotation for *role* across all its sources."""
        return merge_annotations([e['annotation'] for e in ann_sources[role]()])

    # Per-feature state from the drill-down Annotations tab, keyed by the
    # positional uids build_feature_rows hands out.  Cleared whenever a
    # role's sources change, since re-uploading renumbers those uids.
    # Edits accumulate in the *pending* values and only reach the applied
    # ones when the user presses "Apply changes".  Rebuilding the figure on
    # every checkbox would make working through a long feature list
    # unusable, and a debounce cannot help someone who wants to review a
    # set of changes before committing them.
    feature_hidden_pending = reactive.value(frozenset())
    feature_colors_pending = reactive.value({})
    feature_hidden = reactive.value(frozenset())
    feature_colors = reactive.value({})
    # Per-type choices the plot was last drawn with.  The sidebar controls
    # are the pending side; they reach this value automatically while no
    # plot exists yet, and only through "Apply changes" afterwards.  An
    # empty dict means every type is at its default (visible, type colour).
    gff_type_applied = reactive.value({})

    def _reset_feature_overrides() -> None:
        feature_hidden_pending.set(frozenset())
        feature_colors_pending.set({})
        feature_hidden.set(frozenset())
        feature_colors.set({})
        # New sources rebuild the type controls at their defaults, so the
        # applied side resets with them -- anything else would show the
        # fresh controls as phantom pending changes.
        gff_type_applied.set({})

    @reactive.effect
    @reactive.event(input.query_gff)
    def _on_query_gff():
        req(ready())
        _parse_gff_upload(input.query_gff, 'query')

    @reactive.effect
    @reactive.event(input.target_gff)
    def _on_target_gff():
        req(ready())
        _parse_gff_upload(input.target_gff, 'target')

    @reactive.effect
    async def _drop_target_annotations_when_self_aligning():
        """Empty the target annotations when self-aligning.

        A self-alignment puts the query assembly on both axes, so a target
        annotation set has nothing of its own to annotate: left loaded it
        either mismatches every contig or draws the query's own features a
        second time.  Its upload is hidden in this mode, so hiding alone
        would leave it in the merged set with no visible control to remove
        it.

        PAF input is exempt: it has both roles and no self-alignment, and
        ``self_align`` keeps its last value while its panel is hidden.
        """
        if input.input_mode() == 'paf' or not input.self_align():
            return
        if not ann_sources['target']():
            return  # nothing loaded; also stops this re-firing on its own write
        for kind in ('gff', 'genbank'):
            _set_ann_source('target', kind, '', None)
        await session.send_custom_message(
            'rd_clear_file_inputs', {'ids': ['target_gff']}
        )
        ui.notification_show(
            'Self-alignment uses one assembly, so the target annotations were cleared.',
            type='message',
            duration=6,
        )

    # Features the interactive report is currently banding.  Held so a
    # saved figure can carry the same highlights: the report draws its
    # bands from on-screen geometry, but an export is redrawn from
    # coordinates and would otherwise come out plain.
    highlight_bands = reactive.value(())

    _HEX_COLOR_RE = re.compile(r'^#(?:[0-9a-fA-F]{3}|[0-9a-fA-F]{6})$')

    @reactive.effect
    @reactive.event(input.track_bands)
    def _on_track_bands():
        ev = input.track_bands() or {}
        bands = []
        for band in ev.get('bands') or []:
            try:
                bands.append(
                    {
                        'axis': band['axis'],
                        'seqname': band['seqname'],
                        'start': int(band['start']),
                        'end': int(band['end']),
                        # Straight from a sandboxed frame into
                        # matplotlib, which raises on anything it cannot
                        # parse -- accept only plain hex.
                        'color': (
                            band.get('color')
                            if _HEX_COLOR_RE.match(str(band.get('color') or ''))
                            else '#888888'
                        ),
                    }
                )
            except (KeyError, TypeError, ValueError):
                continue  # malformed entry from the frame; ignore it
        highlight_bands.set(tuple(bands))

    _MIN_LEN_NOTIF_ID = 'rd_min_contig_len'

    @reactive.effect
    def _report_excluded_contigs():
        """Say how many contigs the length filter is hiding.

        Silently dropping panels looks like missing data, so keep a
        standing note while the filter is active.  Uses a fixed id so
        adjusting the threshold replaces the note instead of stacking up.
        """
        if result() is None:
            return
        lay = layout()
        dropped = len(lay['excluded_query']) + len(lay['excluded_target'])
        min_len = max(0, int(input.min_contig_len() or 0))
        if not min_len:
            ui.notification_remove(_MIN_LEN_NOTIF_ID)
            return
        if not dropped:
            ui.notification_show(
                f'No contigs are shorter than {min_len:,} bp; nothing excluded.',
                id=_MIN_LEN_NOTIF_ID,
                type='message',
                duration=6,
            )
            return
        ui.notification_show(
            f'Hiding {dropped} contig(s) shorter than {min_len:,} bp. '
            'They are still included in the reordered-FASTA download.',
            id=_MIN_LEN_NOTIF_ID,
            type='message',
            duration=8,
        )

    @reactive.effect
    def _warn_on_annotation_name_mismatch():
        """Warn when a GFF annotates contigs the assembly does not have.

        Keyed on both the result and the annotation sources, so it fires
        whichever order the user supplies them in — a GFF can be uploaded
        long before the first run, and a new assembly can be run against an
        already-loaded GFF.  The plotter only logs about this, which the
        browser user never sees.
        """
        res = result()
        if res is None:
            return
        _kind, _obj, meta = res
        for role in ANNOTATION_ROLES:
            ann = gff_raw_for(role)
            fasta = meta.get(role)
            if ann is None or not isinstance(fasta, SequenceProvider):
                continue
            for warning in validate_annotation_names(
                fasta.names, ann.sequence_names(), role
            ):
                ui.notification_show(warning, type='warning', duration=12)

    @reactive.effect
    @reactive.event(input.clear_gff, ignore_init=True)
    async def _on_clear_gff():
        """Drop every uploaded GFF and reset the file inputs.

        Shiny's file-input binding has a no-op ``setValue``, so the widget
        cannot be cleared from the server: without the custom message the
        filename would linger, and re-selecting the same file would fire no
        change event and so never reload.
        """
        for role in ANNOTATION_ROLES:
            _set_ann_source(role, 'gff', '', None)
        await session.send_custom_message(
            'rd_clear_file_inputs', {'ids': ['query_gff', 'target_gff']}
        )
        ui.notification_show(
            'Cleared uploaded annotations.', type='message', duration=4
        )

    def _read_dynamic(input_id: str, default):
        """Read a dynamically rendered input, tolerating its absence."""
        try:
            return input[input_id]()
        except Exception:  # noqa: BLE001 - silent until the control renders
            # Expected before gff_controls() has rendered; a genuine id typo
            # looks identical, so leave a breadcrumb rather than a silent
            # fallback to the default for every type.
            logger.debug('dynamic input %r not available; using default', input_id)
            return default

    @reactive.calc
    def gff_type_index():
        """Union of feature types across roles, with shared colours.

        Returns ``(rows, slugs, shared)`` where *rows* is one entry per
        normalised type -- ``{'key', 'label', 'roles', 'color'}`` -- and
        *slugs* maps normalised key to input-id fragment.  Colours are
        assigned over the union so the same type never gets two colours on
        the two axes (see core/annotation_colors.py).
        """
        types_by_role = {
            role: list(ann.feature_types())
            for role in ANNOTATION_ROLES
            if (ann := gff_raw_for(role)) is not None
        }
        if not types_by_role:
            return [], {}, {}
        shared = assign_shared_colors(types_by_role)
        spellings: dict[str, list[str]] = {}
        roles_with: dict[str, list[str]] = {}
        for role, fts in types_by_role.items():
            for ft in fts:
                key = normalise_type(ft)
                spellings.setdefault(key, []).append(ft)
                if role not in roles_with.setdefault(key, []):
                    roles_with[key].append(role)
        slugs = type_slug_map(sorted(spellings))
        rows = [
            {
                'key': key,
                'label': display_name(spellings[key]),
                'roles': roles_with[key],
                'color': shared[key],
            }
            for key in sorted(spellings)
        ]
        return rows, slugs, shared

    @render.ui
    def gff_controls():
        rows, slugs, _shared = gff_type_index()
        if not rows:
            return None
        controls = []
        for row in rows:
            badges = ''.join(
                f'<span class="de-gff-role" title="{r} annotations">'
                f'{r[0].upper()}</span>'
                for r in row['roles']
            )
            controls.append(
                ui.div(
                    ui.input_checkbox(f'gtyp_{slugs[row["key"]]}', row['label'], True),
                    ui.HTML(badges),
                    ui.HTML(
                        f'<input type="color" class="de-color-input" '
                        f'id="gcol_{slugs[row["key"]]}" value="{row["color"]}">'
                    ),
                    class_='de-gff-type-row',
                )
            )
        # Deliberately depends on gff_type_index() alone.  The gff_diagonal /
        # gff_tracks toggles live in the static sidebar because anything that
        # made this output re-render — self_panels() used to, via result() —
        # rebuilds every checkbox and colour picker below at its default,
        # silently discarding the user's choices on each run.
        return ui.div(
            ui.div(
                ui.h6('Feature types'),
                *controls,
                # Static within this output: its label and disabled state
                # are driven by ui.update_action_button from
                # _sync_apply_gff_types_button, so the render still
                # depends on gff_type_index() alone.
                ui.input_action_button(
                    'apply_gff_types',
                    'Apply changes',
                    class_='btn-primary btn-sm de-gff-apply',
                    disabled=True,
                ),
                class_='de-gff-section',
            ),
        )

    @reactive.calc
    def gff_type_choices_live() -> dict:
        """Per-type visibility and colour as the sidebar controls stand.

        One control per *normalised* type, shared by both roles.  This is
        the pending side: the figure draws from ``gff_type_applied``,
        which these values reach automatically before the first plot and
        through the "Apply changes" button afterwards.

        A colour picker always holds a real value -- it is initialised to
        the type's assigned colour -- so that default reads back as ``''``
        (no override) here; otherwise every type would count as a pending
        colour change from the moment the controls render.
        """
        rows, slugs, _shared = gff_type_index()
        out = {}
        for row in rows:
            color = str(_read_dynamic(f'gcol_{slugs[row["key"]]}', '') or '')
            if color.lower() == str(row['color']).lower():
                color = ''
            out[row['key']] = (
                bool(_read_dynamic(f'gtyp_{slugs[row["key"]]}', True)),
                color,
            )
        return out

    @reactive.effect
    def _sync_types_before_first_plot():
        """Mirror the type controls straight into the applied state.

        Only while no plot exists yet: choices made before the first run
        should simply be there when the plot appears, without asking for
        an explicit apply.  Once a result exists this backs off and the
        button takes over.
        """
        if result() is not None:
            return
        live = gff_type_choices_live()
        with reactive.isolate():
            if gff_type_applied() != live:
                gff_type_applied.set(live)

    @reactive.effect
    @reactive.event(input.apply_gff_types, ignore_init=True)
    def _apply_gff_types():
        """Commit the pending type toggles/colours, redrawing once."""
        # reactive.event isolates the body, so these reads add no
        # dependencies -- this runs on the button press and nothing else.
        live = gff_type_choices_live()
        if gff_type_applied() != live:
            gff_type_applied.set(live)

    @reactive.effect
    def _sync_apply_gff_types_button():
        """Say whether there is anything to apply, and how much."""
        if result() is None:
            # Auto-apply mode: nothing can ever be pending.
            ui.update_action_button(
                'apply_gff_types', label='Apply changes', disabled=True
            )
            return
        pending = count_pending_type_changes(
            gff_type_choices_live(), gff_type_applied()
        )
        ui.update_action_button(
            'apply_gff_types',
            label=f'Apply changes ({pending})' if pending else 'Apply changes',
            disabled=not pending,
        )

    @reactive.calc
    def annotations():
        """Return (query_ann, target_ann) with the user's type/colour choices."""
        _rows, _slugs, shared = gff_type_index()
        chosen = gff_type_applied()
        hidden = feature_hidden()
        feat_colors = feature_colors()
        result = {}
        for role in ANNOTATION_ROLES:
            ann = gff_raw_for(role)
            if ann is None:
                result[role] = None
                continue
            # Per-feature choices first: their uids are positions in the
            # *unfiltered* record list, so filtering by type first would
            # renumber them.  The type filter runs second and therefore
            # wins — a type switched off hides its features regardless.
            ann = apply_feature_overrides(ann, hidden, feat_colors, role)
            if ann is None:
                result[role] = None
                continue
            fts = list(ann.feature_types())
            enabled = {ft: chosen.get(normalise_type(ft), (True, ''))[0] for ft in fts}
            picked = {key: color for key, (_on, color) in chosen.items() if color}
            colors = color_map_for(fts, {**shared, **picked})
            result[role] = apply_annotation_config(ann, enabled, colors)
        return result['query'], result['target']

    # --- end GFF annotations -------------------------------------------------

    @reactive.calc
    def config() -> PlotConfig:
        return PlotConfig(
            contig_order=input.contig_order(),
            auto_reverse=input.auto_reverse(),
            hide_internal_axes=input.hide_internal_axes(),
            dot_size=dot_size_settled(),
            cap_style=input.cap_style() or 'projecting',
            min_length=min_length_settled(),
            color_by_identity=bool(input.color_by_identity()),
            identity_palette=input.identity_palette() or 'viridis',
        )

    @reactive.calc
    def render_config() -> PlotConfig:
        """Structural plot options — the ones that genuinely need a re-render.

        Display-only options (line width, line cap, min match length) are
        applied client-side inside the embedded report via ``rd_display_opts``
        messages, so they are deliberately absent here: the report HTML is
        rendered with ``min_length=0``, the default line width and the default
        (square) cap, and the client owns them from then on.  The static plot
        and the SVG/PDF downloads keep using the full :func:`config`
        (server-side semantics unchanged there).
        """
        return PlotConfig(
            contig_order=input.contig_order(),
            auto_reverse=input.auto_reverse(),
            hide_internal_axes=input.hide_internal_axes(),
            color_by_identity=bool(input.color_by_identity()),
            identity_palette=input.identity_palette() or 'viridis',
        )

    @reactive.effect
    async def _send_display_opts():
        # Forward the debounced display options to the embedded report
        # (bridge.js relays them into the iframe), where they are applied
        # instantly — no matplotlib re-render in Pyodide.
        await session.send_custom_message(
            'rd_display_opts',
            {
                'dot_size': float(dot_size_settled()),
                # The wire carries the CSS value: the report applies the cap
                # as a stroke-linecap override, not through matplotlib.
                'cap_style': svg_linecap(input.cap_style() or 'projecting'),
                'min_length': int(min_length_settled()),
            },
        )

    # --- W2: interactive plot ------------------------------------------------
    # Focused (query, target) contig pair for the drill-down view, or None
    # for the full overview grid.
    focus = reactive.value(None)
    # Per-result memos, keyed by the alignment object's identity: computed
    # contig orders per ordering mode (gravity is seconds-long on real
    # assemblies) and the DotPlotter/PafAlignment pair (whose construction
    # copies the full record list).  Both are invalidated with the result.
    order_cache: dict = {}
    figure_ctx_cache: dict = {}
    # (query_name, target_name) -> [PafRecord], filled lazily from the
    # current result's records for aligned-sequence lookups.  The grouping
    # mirrors DotPlotter._records_for_pair, so a payload segment index maps
    # straight onto the per-pair list.
    paf_pair_index: dict = {}

    @reactive.effect
    def _reset_focus_on_new_result():
        result()
        focus.set(None)
        order_cache.clear()
        figure_ctx_cache.clear()
        # focus.set(None) is a no-op when the user re-runs from the overview,
        # so the focus-event clear below never fires and match lookups would
        # keep resolving against the previous run's records.
        paf_pair_index.clear()

    @reactive.effect
    @reactive.event(focus)
    def _clear_bands_on_view_change():
        # Bands are gids from one rendered report, so a different focused
        # pair -- or returning to the overview -- leaves them meaningless.
        # A new result resets focus, so this covers that too.
        if highlight_bands():
            highlight_bands.set(())
        paf_pair_index.clear()

    @reactive.calc
    def ordering_config() -> tuple[str, bool]:
        """Return only the config fields that affect contig ordering.

        ``layout()`` depends on this instead of ``config()`` so display-only
        edits (line width, min length, palettes, styling) never re-run the
        gravity ordering — on real assemblies that reorder is the expensive
        step.  Reads the raw inputs directly: going through ``config()``
        would re-invalidate ``layout()`` on every config recompute.
        """
        return (
            input.contig_order(),
            input.auto_reverse(),
            max(0, int(input.min_contig_len() or 0)),
        )

    def _axis_inputs(res) -> tuple[list[str], list[str], dict[str, int]]:
        """Return ``(query_names, target_names, lengths)`` before ordering.

        Shared by :func:`layout` and :func:`grid_panel_count` so the two
        cannot disagree about which contigs the grid contains.
        """
        kind, obj, meta = res
        if kind == 'kmer':
            q_in = list(meta['query'].names)
            t_in = list(meta['target'].names)
            lengths = dict(meta['target'].lengths())
            lengths.update(meta['query'].lengths())
        else:
            q_in = list(obj.query_names)
            t_in = list(obj.target_names)
            lengths = {n: obj.get_sequence_length(n) for n in (*q_in, *t_in)}
        return q_in, t_in, lengths

    @reactive.calc
    def grid_panel_count() -> int:
        """Count the panels the grid will draw, from the result alone.

        Deliberately does NOT go through ``layout()``: that would pull
        ``input.contig_order`` / ``input.auto_reverse`` into ``plot_area``,
        and re-rendering ``plot_area`` rebuilds the report ``<iframe>``,
        throwing away its client-side state (zoom, highlight bands, pushed
        display options, the selected match).  Reading ``result()`` — which
        ``plot_area`` already depends on — costs no extra invalidation.

        This is an upper bound: the ``colinearity`` ordering modes delegate
        to ``compute_gravity_contigs``, which may return fewer contigs.  The
        worst case is one over-advertised navigation tip.
        """
        res = result()
        if res is None:
            return 0
        q_in, t_in, lengths = _axis_inputs(res)
        min_len = max(0, int(input.min_contig_len() or 0))
        q_keep, _ = filter_by_min_length(q_in, lengths, min_len)
        t_keep, _ = filter_by_min_length(t_in, lengths, min_len)
        # Same empty-axis fallback as layout(); the two must agree or the
        # hint advertises a click-to-focus the report has disabled.
        return (len(q_keep) or len(q_in)) * (len(t_keep) or len(t_in))

    @reactive.calc
    def self_panels() -> bool:
        """Whether the panel grid contains a self-comparison panel.

        Gates the diagonal-shading control, which can only ever draw on
        such a panel.  Judged on the whole grid rather than the focused
        pair, so the control does not appear and vanish as the user drills
        in and out.  Returns ``False`` rather than blocking before the
        first run, so the feature-type list stays visible while a plot is
        being set up.
        """
        if input.input_mode() != 'paf' and input.self_align():
            return True
        if result() is None:
            return False
        lay = layout()
        return has_self_pair(lay['query_names'], lay['target_names'])

    @reactive.calc
    def layout():
        """Explicit plotted axis orders for the current result + config.

        Single source of truth for row/column order: the plot call, the
        panel double-click mapping and the FASTA download all use it.
        """
        res = result()
        req(res)
        contig_order, auto_reverse, min_len = ordering_config()
        kind, obj, _meta = res
        # Memoise per ordering mode: re-selecting e.g. 'maximise colinearity'
        # after trying another mode must not recompute the gravity sort.
        # The length filter is part of the key because it changes which
        # contigs the ordering runs over.  auto_reverse stays outside it —
        # it only selects whether the cached reversed set is applied.
        cache_key = (id(obj), contig_order, min_len)
        cached = order_cache.get(cache_key)
        if cached is None:
            q_in, t_in, lengths = _axis_inputs(res)
            records = (
                obj.get_records_for_pair(QUERY_GROUP, TARGET_GROUP)
                if kind == 'kmer'
                else obj.records
            )
            q_keep, q_drop = filter_by_min_length(q_in, lengths, min_len)
            t_keep, t_drop = filter_by_min_length(t_in, lengths, min_len)
            # A threshold above every contig would leave nothing to draw;
            # keep that axis whole rather than rendering an empty grid.
            if not q_keep:
                q_keep, q_drop = q_in, []
            if not t_keep:
                t_keep, t_drop = t_in, []
            cached = (
                *resolve_orders(contig_order, records, q_keep, t_keep, lengths),
                q_drop,
                t_drop,
            )
            order_cache[cache_key] = cached
        else:
            logger.info('contig-order cache hit for mode %r', contig_order)
        q_order, t_order, reversed_q, q_drop, t_drop = cached
        reverse = reversed_q if auto_reverse else set()
        # An active tree overrides every ordering mode (and orientation
        # flipping — the tree's leaf order is authoritative).  Applied on
        # top of the memoised base order so switching the tree on and off
        # never invalidates the gravity-sort cache.
        tree = active_tree()
        if tree is not None:
            ordered = tree_layout_order(tree.leaf_names(), list(q_order))
            if ordered is not None:
                q_order = ordered
                if set(t_order) == set(q_order):
                    t_order = list(ordered)
                reverse = set()
        # Return copies so downstream mutation cannot poison the memo.
        return {
            'query_names': list(q_order),
            'target_names': list(t_order),
            'reverse': set(reverse),
            # Excluded by the length filter.  Kept so the FASTA export can
            # stay complete: CrossIndex.write_fasta writes exactly the names
            # it is given, so anything omitted here is silently lost.
            'excluded_query': list(q_drop),
            'excluded_target': list(t_drop),
        }

    # --- Trees & clustering --------------------------------------------------
    # Whether the optional deps (sourmash + scipy) are importable.  Native
    # runs know at startup; under Pyodide the first 'compute clustering'
    # tick downloads them from the Pyodide channel via micropip.
    cluster_deps_ready = reactive.value(not cluster_deps_missing())

    @reactive.calc
    def self_mode() -> bool:
        """Self-alignment input mode (one assembly on both axes)."""
        return input.input_mode() != 'paf' and bool(input.self_align())

    @reactive.calc
    def query_provider():
        """Return the result's query sequences, when lazily accessible."""
        res = result()
        if res is None:
            return None
        prov = res[2].get('query')
        return prov if isinstance(prov, SequenceProvider) else None

    @reactive.calc
    def clustering_possible() -> bool:
        """Self-comparison with at least 2 contigs and sequences on hand."""
        if not self_mode():
            return False
        prov = query_provider()
        return prov is not None and len(prov.names) >= 2

    @reactive.calc
    def parsed_tree():
        """Parse the uploaded tree (None when absent; errors notify)."""
        files = input.tree_file()
        if not files:
            return None
        from dot_explorer import Tree  # noqa: PLC0415 - after ensure_dot_explorer

        try:
            return Tree.read(files[0]['datapath'])
        except (ValueError, OSError) as exc:
            ui.notification_show(
                f'Could not parse the tree file: {exc}', type='error', duration=12
            )
            return None

    @reactive.calc
    def user_tree():
        """Return the uploaded tree, validated against the query contigs."""
        tree = parsed_tree()
        if tree is None or not clustering_possible():
            return None
        prov = query_provider()
        try:
            tree.validate_labels(list(prov.names))
        except ValueError as exc:
            ui.notification_show(str(exc), type='error', duration=15)
            return None
        return tree

    async def _install_wasm_cluster_packages() -> None:
        """Fetch sourmash + scipy (and their deps) from the Pyodide CDN.

        The shinylive export bundles only the Pyodide packages the app
        needs at startup, so the local distribution has no sourmash/scipy
        wheels and a plain ``micropip.install`` 404s against it.  The
        bundled lockfile still describes every package in the release, so
        resolve the dependency closure there and fetch each missing file
        from the versioned CDN mirror of the same Pyodide release.
        """
        from pyodide.ffi import to_js  # noqa: PLC0415
        import pyodide_js  # noqa: PLC0415 - pyodide-only module

        # The live lockfile object loadPackage resolves from.  Pointing
        # each missing entry's file_name at the CDN (absolute URLs are
        # used verbatim) reroutes the fetch without touching bundled
        # packages — loadPackage-by-URL cannot do that: it accepts *.whl
        # URLs only, and scipy's openblas dependency ships as a .zip.
        lock = pyodide_js._api.lockfile_packages.as_object_map()
        loaded = set(pyodide_js.loadedPackages.as_object_map().keys())
        cdn = f'https://cdn.jsdelivr.net/pyodide/v{pyodide_js.version}/full/'
        names: list[str] = []
        seen: set[str] = set()

        def visit(name: str) -> None:
            if name in seen:
                return
            seen.add(name)
            info = lock.get(name)
            if info is None:
                return
            for dep in info.depends:
                visit(str(dep))
            if name not in loaded:
                names.append(name)

        visit('sourmash')
        visit('scipy')
        for name in names:
            info = lock[name]
            file_name = str(info.file_name)
            if not file_name.startswith(('http://', 'https://')):
                info.file_name = cdn + file_name
        if names:
            await pyodide_js.loadPackage(to_js(names))

    async def ensure_cluster_deps() -> bool:
        """Make sourmash + scipy importable, installing them under Pyodide.

        Deliberately NOT in app/requirements.txt: that file installs before
        the app starts, so listing them would put the ~25 MB download in
        every visitor's first load.  Instead they fetch on first use, like
        the biowasm aligner binaries.
        """
        if not cluster_deps_missing():
            cluster_deps_ready.set(True)
            return True
        if sys.platform == 'emscripten':
            with ui.Progress(min=0, max=1) as progress:
                progress.set(
                    0,
                    message='Downloading clustering packages '
                    '(sourmash + scipy, ~25 MB, one-time)…',
                )
                try:
                    await _install_wasm_cluster_packages()
                except Exception:  # noqa: BLE001 - fall back to micropip
                    logger.exception(
                        'CDN package load failed; falling back to micropip'
                    )
                    import micropip  # noqa: PLC0415 - pyodide-only module

                    await micropip.install(['sourmash', 'scipy'])
                importlib.invalidate_caches()
                progress.set(1, message='Clustering packages installed')
            if cluster_deps_missing():
                ui.notification_show(
                    'Could not download the clustering packages — check the '
                    'network connection and try again.',
                    type='error',
                    duration=12,
                )
                return False
            cluster_deps_ready.set(True)
            return True
        ui.notification_show(
            'Clustering needs sourmash and scipy — install them with: '
            'pip install "dot-explorer[cluster]"',
            type='error',
            duration=12,
        )
        return False

    @reactive.effect
    @reactive.event(input.cluster_enabled)
    async def _install_cluster_deps():
        if input.cluster_enabled():
            if await ensure_cluster_deps() and cluster_settings() is None:
                # Seed the applied settings so first activation works
                # without an extra Apply click; later edits wait for it.
                cluster_settings.set(_snapshot_cluster_inputs())

    def _snapshot_cluster_inputs() -> dict:
        """Read every Trees & clustering input into a plain settings dict."""
        return {
            'metric': input.cluster_metric() or 'jaccard',
            'ksize': max(4, int(input.sketch_k() or 21)),
            'scaled': max(1, int(input.sketch_scaled() or 1000)),
            'abund': bool(input.sketch_abund()),
            'mode': input.cluster_mode() or 'similarity',
            'cutoff': float(input.cluster_cutoff()),
            'cutoff_line': bool(input.cluster_cutoff_line()),
            'identity_cutoff': float(input.identity_cutoff()),
            'coverage_cutoff': float(input.coverage_cutoff()),
            'reciprocal': bool(input.cov_reciprocal()),
            'coverage_source': input.coverage_source() or 'containment',
            'borders_on': bool(input.cluster_borders_on()),
            'cmap': input.heatmap_cmap() or 'viridis',
            'cell_values': bool(input.heatmap_values()),
        }

    # Applied Trees & clustering settings.  Sidebar edits do nothing until
    # the Apply button snapshots them here (same held-until-applied flow as
    # the annotations table), so slider drags never trigger a re-sketch or
    # figure redraw.  None until clustering is first enabled.
    cluster_settings = reactive.value(None)

    @reactive.effect
    @reactive.event(input.apply_cluster)
    def _apply_cluster_settings():
        snap = _snapshot_cluster_inputs()
        if snap != cluster_settings():
            cluster_settings.set(snap)

    @reactive.calc
    def clustering_on() -> bool:
        """Whether the sourmash-based clustering pipeline should run."""
        return (
            clustering_possible()
            and bool(input.cluster_enabled())
            and cluster_deps_ready()
            and cluster_settings() is not None
        )

    @reactive.calc
    def sketch_params():
        from dot_explorer import SketchParams  # noqa: PLC0415

        settings = cluster_settings() or {}
        return SketchParams(
            ksize=settings.get('ksize', 21),
            scaled=settings.get('scaled', 1000),
            track_abundance=settings.get('abund', True),
        )

    @reactive.calc
    def sketches():
        """One sourmash sketch per query contig (cached until k/scaled change)."""
        req(clustering_on())
        prov = query_provider()
        req(prov)
        from dot_explorer import compute_sketches  # noqa: PLC0415

        with ui.Progress(min=0, max=1) as progress:
            progress.set(0, message='Sketching contigs (sourmash)…')
            sk = compute_sketches(ProviderIndex(prov), params=sketch_params())
            progress.set(1, message='Sketches ready')
        return sk

    @reactive.calc
    def sim_matrix():
        """All-vs-all similarity in the applied metric, or None on error."""
        if not clustering_on():
            return None
        from dot_explorer import pairwise_similarity  # noqa: PLC0415

        metric = cluster_settings()['metric']
        if metric == 'ani':
            return ani_matrix()  # shared with dual mode; carries the CIs
        try:
            return pairwise_similarity(sketches(), metric=metric)
        except ValueError as exc:
            # e.g. angular metric on abundance-free sketches.
            ui.notification_show(str(exc), type='error', duration=10)
            return None

    @reactive.calc
    def containment_matrix():
        """Asymmetric sourmash containment matrix, or None."""
        if not clustering_on():
            return None
        from dot_explorer import pairwise_similarity  # noqa: PLC0415

        try:
            return pairwise_similarity(sketches(), metric='containment')
        except ValueError as exc:
            ui.notification_show(str(exc), type='error', duration=10)
            return None

    def _display_name(name: str) -> str:
        """Strip a CrossIndex ``group:`` prefix from an internal name."""
        if name.startswith(('query:', 'target:')):
            return name.split(':', 1)[1]
        return name

    @reactive.calc
    def coverage_records():
        """Alignment records fit to drive coverage clustering, or None.

        Only clean minimap2 output (no ``-P``) qualifies: the displayed
        result when it is one, else a cached or background-computed
        dedicated run (``_ensure_coverage_alignment`` dispatches it).
        nucmer / k-mer / uploaded-PAF results never feed coverage.
        """
        res = result()
        prov = query_provider()
        if res is None or prov is None:
            return None
        kind, obj, meta = res
        if kind == 'paf' and is_clean_minimap2(meta.get('method'), meta.get('params')):
            return obj.records
        cov = coverage_alignment()
        if cov is not None and cov[0] == prov.digest:
            return cov[1].records
        cached = cache.get_paf(
            'minimap2', coverage_align_params(), prov.digest, prov.digest
        )
        if cached is not None:
            return cached.records
        return None

    @reactive.effect
    async def _ensure_coverage_alignment():
        """Dispatch the background minimap2 run when coverage needs one."""
        if not clustering_on():
            return
        settings = cluster_settings()
        if (
            settings['mode'] != 'identity_coverage'
            or settings['coverage_source'] != 'alignment'
        ):
            return
        if coverage_records() is not None:
            return
        prov = query_provider()
        if prov is None or coverage_failed() == prov.digest:
            return
        pending = coverage_pending()
        if pending is not None and pending['query'].digest == prov.digest:
            return
        request_id = uuid.uuid4().hex
        coverage_pending.set({'request_id': request_id, 'query': prov})
        ui.notification_show(
            'Alignment coverage needs a clean minimap2 run (no -P) — '
            'computing one in the background; the dot plot keeps showing '
            'the current alignment.',
            id=_COV_NOTIF_ID,
            duration=None,
        )
        try:
            await _send_dataset(prov)
        except ValueError as exc:
            coverage_pending.set(None)
            coverage_failed.set(prov.digest)
            ui.notification_remove(_COV_NOTIF_ID)
            ui.notification_show(str(exc), type='error', duration=8)
            return
        await session.send_custom_message(
            'rd_run_aligner',
            {
                'tool': 'minimap2',
                'args': build_tool_args('minimap2', coverage_align_params()),
                'query_id': prov.digest,
                'target_id': prov.digest,
                'request_id': request_id,
            },
        )

    @reactive.calc
    def alignment_coverage():
        """Asymmetric coverage from clean minimap2 records, or None."""
        prov = query_provider()
        records = coverage_records()
        if records is None or prov is None:
            return None
        return alignment_coverage_matrix(
            records,
            list(prov.names),
            dict(prov.lengths()),
            normalize=_display_name,
        )

    @reactive.calc
    def linkage_and_tree():
        """(scipy linkage, Tree) for the computed similarity, or None."""
        sim = sim_matrix()
        if sim is None:
            return None
        from dot_explorer import Tree, linkage_from_similarity  # noqa: PLC0415

        Z = linkage_from_similarity(sim)
        return Z, Tree.from_linkage(Z, sim.names)

    @reactive.calc
    def active_tree():
        """Return the tree ordering the plot (user upload wins)."""
        tree = user_tree()
        if tree is not None:
            return tree
        if not clustering_on():
            return None
        lt = linkage_and_tree()
        return None if lt is None else lt[1]

    @reactive.calc
    def ani_matrix():
        """ANI matrix with confidence bounds, or None."""
        if not clustering_on():
            return None
        from dot_explorer import pairwise_similarity  # noqa: PLC0415

        with ui.Progress(min=0, max=1) as progress:
            progress.set(0, message='Computing ANI…')
            ani = pairwise_similarity(sketches(), metric='ani')
            progress.set(1, message='Done')
        if ani.ci_low is not None and bool(numpy.all(ani.ci_low == ani.values)):
            # Every pair fell back to the raw distance (sketches keep too
            # few hashes for sourmash to trust its ANI model).
            ui.notification_show(
                'ANI confidence intervals unavailable: the sketches keep '
                'too few hashes at scaled='
                f'{cluster_settings()["scaled"]} for these sequence '
                'lengths — lower the sketch scaled factor.',
                type='warning',
                duration=12,
            )
        return ani

    @reactive.calc
    def dual_matrices():
        """(ANI, coverage) matrices for identity+coverage clustering.

        Coverage comes from the applied source: sourmash containment, or
        block coverage from a clean minimap2 alignment (SNP-robust). When
        the background minimap2 run is still in flight this returns None
        (clustering appears once it lands); when it failed, sourmash
        containment stands in.
        """
        if not clustering_on():
            return None
        ani = ani_matrix()
        if ani is None:
            return None
        if cluster_settings()['coverage_source'] == 'alignment':
            cov = alignment_coverage()
            if cov is None:
                prov = query_provider()
                if prov is None or coverage_failed() != prov.digest:
                    return None  # background minimap2 run pending
                cov = containment_matrix()
        else:
            cov = containment_matrix()
        if cov is None:
            return None
        return ani, cov

    @reactive.calc
    def cluster_result():
        """Cluster assignments for the applied mode, or None."""
        if not clustering_on():
            return None
        from dot_explorer import assign_clusters, assign_clusters_dual  # noqa: PLC0415

        settings = cluster_settings()
        if settings['mode'] == 'identity_coverage':
            matrices = dual_matrices()
            if matrices is None:
                return None
            ani, cov = matrices
            return assign_clusters_dual(
                ani,
                cov,
                identity_cutoff=settings['identity_cutoff'],
                coverage_cutoff=settings['coverage_cutoff'],
                reciprocal=settings['reciprocal'],
            )
        sim = sim_matrix()
        lt = linkage_and_tree()
        if sim is None or lt is None:
            return None
        return assign_clusters(sim, settings['cutoff'], linkage=lt[0])

    @reactive.calc
    def tree_cutoff_distance():
        """Distance-from-tips for the dendrogram cutoff line, or None."""
        if not clustering_on() or user_tree() is not None:
            # User trees have arbitrary branch-length units.
            return None
        settings = cluster_settings()
        if settings['mode'] == 'similarity' and settings['cutoff_line']:
            return 1.0 - settings['cutoff']
        return None

    @reactive.effect
    @reactive.event(input.cluster_select)
    async def _on_cluster_select():
        """Forward the cluster-table row selection to the report iframe."""
        sel = input.cluster_select() or {}
        names = sel.get('clusters') if isinstance(sel, dict) else None
        # Shiny deserialises JS arrays as tuples; JSON needs a list back.
        names = list(names) if isinstance(names, (list, tuple)) else []
        await session.send_custom_message(
            'rd_highlight_clusters',
            {'clusters': [str(n) for n in names]},
        )

    @reactive.effect
    async def _clear_cluster_selection_on_change():
        """Reset any highlight when the assignments themselves change."""
        cluster_result()
        await session.send_custom_message('rd_highlight_clusters', {'clusters': []})

    # --- end Trees & clustering ----------------------------------------------

    def make_figure(res, cfg: PlotConfig, lay, pair=None, output_path=None):
        from dot_explorer import DotPlotter
        from dot_explorer.paf_io import PafAlignment

        kind, obj, _meta = res
        q_names = [pair[0]] if pair is not None else lay['query_names']
        t_names = [pair[1]] if pair is not None else lay['target_names']
        kwargs = cfg.plot_kwargs()
        # Ordering is fully explicit (see layout()); never reorder in-plot.
        kwargs.update(contig_order=None, auto_reverse=False)
        if kind == 'kmer':
            # k-mer matches are exact: identity is uniformly 100%, so
            # identity colouring would be a misleading single-colour plot.
            kwargs['color_by_identity'] = False
        # Identity colouring is unreadable without a key.
        kwargs['identity_colorbar'] = bool(kwargs.get('color_by_identity'))
        kwargs['reverse_contigs'] = set(lay['reverse'])
        if pair is not None and cfg.title is None:
            kwargs['title'] = f'{pair[0]} vs {pair[1]}'
        if pair is None:
            # Tree gutter + cluster borders on the overview grid only (the
            # focused single-pair view has no rows to order).  The tree is
            # attached only when its tips exactly cover the plotted rows —
            # a min-contig-length filter can hide tips, in which case the
            # leaf *order* still applies (see layout()) but the dendrogram
            # itself would be wrong to draw.
            tree = active_tree()
            if (
                tree is not None
                and len(q_names) >= 2
                and set(tree.leaf_names()) == set(q_names)
            ):
                kwargs['tree'] = tree
                kwargs['tree_cutoff'] = tree_cutoff_distance()
            clusters = cluster_result()
            settings = cluster_settings() or {}
            if clusters is not None and settings.get('borders_on', True):
                kwargs['cluster_borders'] = clusters
        # GFF annotations: diagonal shading on self panels plus side tracks
        # in the focused (1x1) drill-down view.  Reading annotations() here
        # keeps every figure consumer reactive to toggle/colour changes
        # without touching layout()'s dependencies.
        ann_q, ann_t = annotations()
        if ann_q is not None or ann_t is not None:
            # Static inputs now, so read them directly.  The self_panels()
            # conjunct still matters: a hidden checkbox keeps its last value,
            # so a box ticked during a self-comparison must not shade the
            # diagonal of a later cross-assembly run.
            if self_panels() and bool(input.gff_diagonal()):
                # On a self-comparison both roles describe the same
                # sequences, so merge them: a user who uploaded the
                # annotation under either role — or split it across two
                # files — gets everything shaded.
                kwargs['annotation'] = merge_annotations([ann_q, ann_t])
            kwargs['annotation_query'] = ann_q
            kwargs['annotation_target'] = ann_t
            kwargs['annotation_tracks'] = bool(input.gff_tracks())
        if output_path is None and highlight_bands():
            # Only for figures the user takes away.  The HTML report draws
            # its own bands live, and would end up with two.
            kwargs['highlight_regions'] = [dict(b) for b in highlight_bands()]
        # Memoise the plotter per result: for the k-mer path, building the
        # PafAlignment copies the full record list — pointless to repeat on
        # every re-render of the same result.
        plotter = figure_ctx_cache.get(id(obj))
        if plotter is None:
            if kind == 'kmer':
                # Internal 'group:name' identifiers keep name resolution
                # unambiguous; the cached cross-group records are handed to
                # the plotter directly so nothing is recomputed.
                paf = PafAlignment(obj.get_records_for_pair(QUERY_GROUP, TARGET_GROUP))
                plotter = DotPlotter(obj, paf_alignment=paf)
            else:
                plotter = DotPlotter(obj)
            figure_ctx_cache[id(obj)] = plotter
        if kind == 'kmer':
            kwargs['query_names'] = [
                obj.make_internal_name(QUERY_GROUP, n) for n in q_names
            ]
            kwargs['target_names'] = [
                obj.make_internal_name(TARGET_GROUP, n) for n in t_names
            ]
        else:
            kwargs['query_names'] = list(q_names)
            kwargs['target_names'] = list(t_names)
        if output_path is not None:
            # embed_sequences stays at its False default: the in-app report
            # fetches previews/copies lazily over the bridge, so embedding
            # would only bloat the payload.
            return plotter.to_html(output_path, **kwargs)
        return plotter.plot(**kwargs)

    def _report_html(pair=None) -> str:
        """Render the interactive HTML report and return it as a string.

        Uses :func:`render_config` (structural options only): line width
        and min match length are applied client-side in the embedded
        report, so changing them never re-runs this seconds-long render.
        """
        import matplotlib.pyplot as plt

        res = result()
        req(res)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'report.html'
            fig = make_figure(
                res, render_config(), layout(), pair=pair, output_path=path
            )
            html = path.read_text(encoding='utf-8')
        plt.close(fig)
        return strip_report_header(inject_panel_bridge(html))

    @reactive.calc
    def overview_html() -> str:
        # Cached by reactive.calc: entering/leaving the drill-down view does
        # not invalidate it, so 'Back to overview' re-embeds the same string.
        return _report_html(None)

    @reactive.calc
    def focus_html() -> str:
        pair = focus()
        req(pair)
        return _report_html(pair)

    @render.ui
    def aligner_log_ui():
        entries = aligner_log()
        if not entries:
            return None
        blocks = []
        for e in entries:
            text = f'$ {e["cmd"]}\n{e["stderr"] or "(no tool output)"}'
            if e.get('error'):
                text += f'\nERROR: {e["error"]}'
            blocks.append(ui.tags.pre(text, class_='de-log'))
        return ui.accordion(
            ui.accordion_panel(
                f'Aligner log ({len(entries)} run(s))', *blocks, value='log'
            ),
            open=False,
            class_='de-log-accordion',
        )

    @output(suspend_when_hidden=False)
    @render.ui
    def report_frame():
        # suspend_when_hidden=False: switching to the drill-down's
        # Annotations tab hides this output, and Shiny would re-render it
        # on the way back — a seconds-long matplotlib pass under Pyodide
        # for a view that has not changed.
        #
        # Rendering the report re-runs layout (contig ordering) and the
        # matplotlib figure — seconds-long on real assemblies, so show what
        # is happening rather than freezing silently.
        with ui.Progress(min=0, max=2) as progress:
            progress.set(1, message='Computing layout & rendering report…')
            html = focus_html() if focus() is not None else overview_html()
            progress.set(2, message='Done')
        return ui.tags.iframe(
            srcdoc=html,
            class_='de-report-frame',
            sandbox='allow-scripts',
            title='Interactive dotplot report',
        )

    def _nav_hint(focused: bool, multi_panel: bool):
        """Build the navigation-tips box shown under the interactive report.

        Parameters
        ----------
        focused : bool
            Whether the focused single-pair view is active; panel-click tips
            do not apply there (click-to-focus is disabled on single-panel
            reports).
        multi_panel : bool
            Whether the grid has more than one panel.  See
            :func:`core.panels.nav_tips`.

        Returns
        -------
        htmltools.Tag
            The hint ``div`` with each action term in bold.
        """
        tips = nav_tips(focused, multi_panel)
        parts: list = [ui.tags.b('Navigate: ')]
        for i, (action, effect) in enumerate(tips):
            if i:
                parts.append(' · ')
            parts += [ui.tags.b(action), f' = {effect}']
        return ui.div(*parts, class_='de-nav-hint')

    @render.ui
    def plot_area():
        if result() is None:
            return None
        pair = focus()
        # The fullscreen button is a plain client-side toggle (state lives
        # as a class on <html>, wired in www/fullscreen.js, so it survives
        # this output's frequent re-renders); no Shiny input involved.
        toolbar = [
            ui.tags.button(
                ui.HTML(_FS_EXPAND_SVG),
                ui.HTML(_FS_COLLAPSE_SVG),
                class_='de-fs-btn',
                type='button',
                title='Toggle fullscreen (Esc exits)',
                aria_label='Toggle fullscreen plot',
            )
        ]
        if pair is not None:
            toolbar.append(
                ui.input_action_button(
                    'back_overview',
                    'Back to overview',
                    class_='btn-primary btn-sm de-back-btn',
                )
            )
            toolbar.append(ui.span(f'{pair[0]} vs {pair[1]}', class_='de-focus-label'))
        if input.interactive():
            body = ui.output_ui('report_frame')
        else:
            body = ui.output_plot('dotplot', height='72vh')
        hint = (
            _nav_hint(pair is not None, grid_panel_count() > 1)
            if input.interactive()
            else None
        )
        if pair is None and cluster_result() is not None:
            # Clustering active: the overview gets its own tab strip with
            # the assignment table and similarity heatmap beside the plot.
            return ui.div(
                ui.div(*toolbar, class_='de-plot-toolbar'),
                ui.navset_tab(
                    ui.nav_panel('Plot', body, hint),
                    ui.nav_panel(
                        'Clusters',
                        ui.div(
                            ui.download_button(
                                'dl_clusters_csv',
                                'Cluster table (CSV)',
                                class_='btn-sm',
                            ),
                            class_='de-cluster-actions',
                        ),
                        ui.output_ui('cluster_table'),
                    ),
                    ui.nav_panel(
                        'Matrix',
                        ui.div(
                            ui.input_radio_buttons(
                                'matrix_view',
                                None,
                                choices={
                                    'similarity': 'Similarity (applied metric)',
                                    'containment': 'Containment (sourmash)',
                                    'coverage': 'Coverage (alignment)',
                                },
                                inline=True,
                            ),
                            class_='de-cluster-actions',
                        ),
                        ui.output_ui('matrix_table'),
                    ),
                    ui.nav_panel(
                        'Heatmap',
                        ui.div(
                            ui.output_image('heatmap_plot', inline=True),
                            class_='de-heatmap-wrap',
                        ),
                        ui.div(
                            ui.tags.b('Navigate: '),
                            ui.tags.b('scroll'),
                            ' = pan up/down · ',
                            ui.tags.b('Shift+scroll'),
                            ' = pan left/right · ',
                            ui.tags.b('Cmd/Ctrl+scroll'),
                            ' = zoom · ',
                            ui.tags.b('drag'),
                            ' = zoom to region · ',
                            ui.tags.b('double-click / Esc'),
                            ' = reset',
                            class_='de-nav-hint',
                        ),
                    ),
                    id='overview_tabs',
                ),
                class_='de-plot-area',
            )
        if pair is None or not feature_rows():
            return ui.div(
                ui.div(*toolbar, class_='de-plot-toolbar'),
                body,
                hint,
                class_='de-plot-area',
            )
        # Tabs only in the drill-down, so the overview path is untouched.
        return ui.div(
            ui.div(*toolbar, class_='de-plot-toolbar'),
            ui.navset_tab(
                ui.nav_panel('Plot', body, hint),
                ui.nav_panel(
                    'Annotations',
                    ui.div(
                        ui.input_action_button(
                            'apply_features',
                            'Apply changes',
                            class_='btn-primary btn-sm',
                            disabled=True,
                        ),
                        ui.span(
                            'Show/hide and colour edits are held until you apply them.',
                            class_='de-ft-apply-note',
                        ),
                        class_='de-ft-apply',
                    ),
                    ui.output_ui('annotation_table'),
                ),
                id='drill_tabs',
            ),
            class_='de-plot-area',
        )

    @reactive.calc
    def feature_rows() -> list[dict]:
        """Every feature on the focused pair's sequences, as table rows."""
        pair = focus()
        if pair is None:
            return []
        # Both roles are always listed: on a self panel the two axes carry
        # the same sequence but may still hold different uploads.
        q_ann = gff_raw_for('query')
        t_ann = gff_raw_for('target')
        rows = build_feature_rows(q_ann, pair[0], 'query')
        if t_ann is not None:
            rows += build_feature_rows(t_ann, pair[1], 'target')
        elif pair[1] != pair[0]:
            # No target sources (self-alignment clears them): the target
            # axis still carries the query assembly, so its features live
            # in the query annotation set.  uids stay query-role so
            # show/hide/colour edits hit the records that are drawn; only
            # the displayed axis says where the feature sits in this pair.
            fallback = build_feature_rows(q_ann, pair[1], 'query')
            for r in fallback:
                r['axis'] = 'target'
            rows += fallback
        if pair[0] == pair[1]:
            # Same sequence on both axes: the same upload under both roles
            # would list every feature twice.
            rows = dedupe_feature_rows(rows)
        return rows

    @reactive.effect
    @reactive.event(input.feature_table_change)
    def _on_feature_table_change():
        """Fold one delta from the annotations table into the override state.

        The table ships raw HTML inputs and posts deltas through
        www/feature-table.js rather than binding one Shiny input per
        feature: the real GFFs here run to ~600 features per sequence, and
        ~1200 bound inputs visibly freezes Pyodide.
        """
        ev = input.feature_table_change()
        if not ev:
            return
        kind = ev.get('kind')
        uids = ev.get('uids') or ([ev['uid']] if ev.get('uid') else [])
        if not uids:
            return
        if kind in ('vis', 'bulk'):
            visible = bool(ev.get('value'))
            hidden = set(feature_hidden_pending())
            hidden.difference_update(uids) if visible else hidden.update(uids)
            feature_hidden_pending.set(frozenset(hidden))
        elif kind == 'color':
            colors = dict(feature_colors_pending())
            value = ev.get('value') or ''
            for uid in uids:
                if value:
                    colors[uid] = value
                else:
                    colors.pop(uid, None)  # reset to the type colour
            feature_colors_pending.set(colors)

    @reactive.effect
    @reactive.event(input.apply_features, ignore_init=True)
    def _apply_feature_overrides_now():
        """Commit the pending per-feature edits, redrawing once."""
        # reactive.event isolates the body, so these reads add no
        # dependencies -- this runs on the button press and nothing else.
        if feature_hidden_pending() == feature_hidden() and (
            feature_colors_pending() == feature_colors()
        ):
            return  # equal-but-new objects would still invalidate
        feature_hidden.set(feature_hidden_pending())
        feature_colors.set(dict(feature_colors_pending()))

    @reactive.effect
    def _sync_apply_features_button():
        """Say whether there is anything to apply, and how much."""
        pending = count_pending_overrides(
            feature_hidden_pending(),
            feature_colors_pending(),
            feature_hidden(),
            feature_colors(),
        )
        ui.update_action_button(
            'apply_features',
            label=f'Apply changes ({pending})' if pending else 'Apply changes',
            disabled=not pending,
        )

    def _feature_table_payload() -> dict | None:
        """Everything feature-table.js needs to build the table body."""
        rows = feature_rows()
        pair = focus()
        if not rows or pair is None:
            return None
        # Isolated: the table shows the pending edits whenever it *is*
        # rebuilt, but must not be re-sent on every toggle -- the client
        # already painted the change it posted.
        with reactive.isolate():
            hidden = feature_hidden_pending()
            overrides = feature_colors_pending()
        _type_rows, _slugs, shared = gff_type_index()

        def type_color(ft: str) -> str:
            return shared.get(normalise_type(ft), '#888888')

        return {
            'pair': [pair[0], pair[1]],
            'rows': [
                {
                    'uid': r['uid'],
                    # Which axis of the focused pair the feature sits on;
                    # usually the role, but self-alignments list the query
                    # set on both axes.
                    'role': r.get('axis', r['role']),
                    'contig': r['seqname'],
                    'type': r['type'],
                    'name': r['name'] or '',
                    'start': r['start'],
                    'end': r['end'],
                    'length': r['length'],
                    'strand': r['strand'],
                    'source': r['source_file'] or r['source'] or '',
                    'attrs': ' · '.join(
                        f'{k}={v}' for k, v in list(r['attributes'].items())[:6]
                    ),
                    'visible': r['uid'] not in hidden,
                    'color': overrides.get(r['uid']) or type_color(r['type']),
                    'type_color': type_color(r['type']),
                }
                for r in rows
            ],
        }

    @reactive.effect
    async def _push_feature_rows():
        """Ship the focused pair's feature rows to the browser.

        The table body is built client-side from this message: rendering
        ~600 rows per sequence as raw HTML through a ``@render.ui`` was
        what made the Plot <-> Annotations tab switch stall, since Shiny
        re-inserted and re-scanned the whole blob on every re-show.
        """
        payload = _feature_table_payload()
        if payload is None:
            return
        await session.send_custom_message('rd_feature_rows', payload)

    @output(suspend_when_hidden=False)
    @render.ui
    def annotation_table():
        # Static shell only -- caption and body rows arrive through the
        # 'rd_feature_rows' custom message and are built by
        # www/feature-table.js.  suspend_when_hidden=False so returning to
        # the tab re-shows the already-built DOM instead of re-rendering.
        if not feature_rows() or focus() is None:
            return None
        # Each header carries how its column compares, so the client can
        # sort without a round trip -- the table deliberately does not
        # re-render, and re-rendering would drop the pending edits.
        header = (
            '<thead><tr>'
            + ''.join(
                f'<th data-sort="{kind}" role="columnheader" '
                f'aria-sort="none" tabindex="0" '
                f'title="Sort by {label.lower()}">{label}</th>'
                for label, kind in _FEATURE_COLUMNS
            )
            + '</tr></thead>'
        )
        return ui.div(
            ui.div('', class_='de-ft-caption'),
            ui.div(
                ui.HTML(
                    '<span class="de-ft-filters">'
                    '<select id="de-ft-filter-col" '
                    'title="Column the filter applies to"></select>'
                    '<input type="search" id="de-ft-filter" '
                    'placeholder="Filter text…">'
                    '<button type="button" id="de-ft-filter-add">'
                    'Add filter</button>'
                    '<button type="button" id="de-ft-filter-apply">'
                    'Apply filters</button>'
                    '</span>'
                    '<button type="button" data-bulk="show">Show all</button>'
                    '<button type="button" data-bulk="hide">Hide all</button>'
                    '<button type="button" data-bulk="reset">Reset colours</button>'
                ),
                class_='de-ft-tools',
            ),
            ui.div(ui.HTML(''), id='de-ft-chips', class_='de-ft-chips'),
            ui.HTML(
                f'<div class="de-ft-scroll"><table id="de-feature-table">'
                f'{header}<tbody></tbody></table></div>'
            ),
            ui.div(
                'Click a row to select it (⌘/Ctrl-click for several); '
                'click a selected row again or press Esc to clear. '
                'Selected features are banded on the Plot tab.',
                class_='de-ft-hint',
            ),
            class_='de-ft-panel',
        )

    @reactive.effect
    @reactive.event(input.panel_dblclick)
    def _drill_down():
        if result() is None or focus() is not None:
            # No result yet, or already in the standalone single-pair view
            # (whose only panel is (0, 0) and must not remap).
            return
        info = input.panel_dblclick()
        lay = layout()
        try:
            q, t = panel_pair(
                lay['query_names'],
                lay['target_names'],
                int(info['row']),
                int(info['col']),
            )
        except (KeyError, TypeError, ValueError, IndexError):
            logger.warning('Ignoring unmappable panel_dblclick payload: %r', info)
            return
        focus.set((q, t))

    @reactive.effect
    @reactive.event(input.back_overview)
    def _back_to_overview():
        focus.set(None)

    def _pair_records(obj, q: str, t: str) -> list:
        """Return the PAF records for one pair, building the index lazily."""
        if not paf_pair_index:
            for rec in obj.records:
                paf_pair_index.setdefault((rec.query_name, rec.target_name), []).append(
                    rec
                )
        return paf_pair_index.get((q, t), [])

    def _record_for_segment(obj, q: str, t: str, idx: int):
        """Return the PAF record behind identity-layer payload row *idx*."""
        recs = _pair_records(obj, q, t)
        return recs[idx] if 0 <= idx < len(recs) else None

    def _genomic_coords(info: dict, q: str, lay, qlen: int):
        """Map a clicked segment's payload coords to genomic values.

        The payload's query side is mirrored on reverse-oriented contigs;
        the target side is always genomic.

        Returns
        -------
        tuple or None
            ``(q_start, q_end, t_start, t_end, strand)`` or ``None`` when
            the message lacks usable coordinates.
        """
        try:
            gqs, gqe = int(info['qs']), int(info['qe'])
            ts_, te_ = int(info['ts']), int(info['te'])
        except (KeyError, TypeError, ValueError):
            return None
        strand = info.get('strand', '+')
        if strand not in ('+', '-'):
            strand = '+'
        if q in lay['reverse']:
            gqs, gqe = qlen - gqe, qlen - gqs
            strand = '-' if strand == '+' else '+'
        return gqs, gqe, ts_, te_, strand

    def _preview_slice(seq: str, start: int, end: int, minus: bool) -> str:
        """Return an alignment-oriented preview of ``seq[start:end]``.

        Clipped to 1,000 bases *before* any copy or revcomp so megabase
        matches never materialise a full slice for the preview.  On the
        minus strand the alignment-oriented sequence begins at the genomic
        end, so the window is taken from there.
        """
        from dot_explorer.alignment_view import revcomp

        cap = 1_000
        n = end - start
        if n <= cap:
            s = seq[start:end]
            return revcomp(s) if minus else s
        if minus:
            s = revcomp(seq[end - cap : end])
        else:
            s = seq[start : start + cap]
        return s + f'… [truncated at {cap:,} bases]'

    def _sequence_for(meta: dict, name: str, side: str):
        """Look up a lazily sliceable sequence by name, preferring *side*.

        Falls back to the other assembly so self-align mode (one file) and
        PAF uploads with a single FASTA still resolve.  The returned object
        supports ``len()`` and ``[start:stop]`` slicing to ``str``; with a
        faidx-backed provider only the sliced windows are ever read.
        """
        order = ('query', 'target') if side == 'query' else ('target', 'query')
        for key in order:
            provider = meta.get(key)
            if provider is not None:
                seq = provider.get_lazy(name)
                if seq is not None:
                    return seq
        return None

    def _match_context(info):
        """Resolve a clicked segment to sequences and genomic coordinates.

        Shared by the match-detail and copy-request handlers.

        Returns
        -------
        dict or None
            ``{'q', 't', 'rec', 'qseq', 'tseq', 'gqs', 'gqe', 'ts', 'te',
            'strand'}`` with genomic (unmirrored) coordinates, or ``None``
            when the panel or sequences cannot be resolved.  ``rec`` is
            the exact PAF record when one matches (its CIGAR, if any,
            enables the gapped view).
        """
        res = result()
        if not info or res is None:
            return None
        kind, obj, meta = res
        pair = focus()
        if pair is not None:
            q, t = pair
            lay = layout()
        else:
            try:
                lay = layout()
                q, t = panel_pair(
                    lay['query_names'],
                    lay['target_names'],
                    int(info['row']),
                    int(info['col']),
                )
            except (KeyError, TypeError, ValueError, IndexError):
                return None
        qseq = _sequence_for(meta, q, 'query')
        tseq = _sequence_for(meta, t, 'target')
        if qseq is None or tseq is None:
            return None
        rec = None
        if kind == 'paf' and info.get('layer') == 'identity':
            try:
                rec = _record_for_segment(obj, q, t, int(info['idx']))
            except (TypeError, ValueError):
                rec = None
        # Strand-coloured (fwd/rev) layers have no stable index->record
        # mapping (blocks may be chained), but an unchained block's coords
        # match its record exactly — recover the CIGAR that way so runs
        # with -c get the gapped view without identity colouring on.
        if rec is None and kind == 'paf':
            g = _genomic_coords(info, q, lay, len(qseq))
            if g is not None:
                gqs0, gqe0, ts0, te0, gstrand0 = g
                rec = next(
                    (
                        r
                        for r in _pair_records(obj, q, t)
                        if r.query_start == gqs0
                        and r.query_end == gqe0
                        and r.target_start == ts0
                        and r.target_end == te0
                        and r.strand == gstrand0
                    ),
                    None,
                )
        if rec is not None:
            gqs, gqe = rec.query_start, rec.query_end
            ts_, te_ = rec.target_start, rec.target_end
            strand = rec.strand
        else:
            g = _genomic_coords(info, q, lay, len(qseq))
            if g is None:
                return None
            gqs, gqe, ts_, te_, strand = g
        return {
            'q': q,
            't': t,
            'rec': rec,
            'qseq': qseq,
            'tseq': tseq,
            'gqs': gqs,
            'gqe': gqe,
            'ts': ts_,
            'te': te_,
            'strand': strand,
        }

    @reactive.effect
    @reactive.event(input.match_select)
    async def _on_match_select():
        """Serve the sequence preview for a clicked match.

        The report posts 'de-match-select'; the reply travels back through
        bridge.js as an 'de-seq-response' echoing row/col/layer/idx so the
        report can discard stale responses.  When the record has a CIGAR
        (minimap2 ``-c``) the reply is a gapped alignment view; otherwise
        the raw query and target slices are shown unaligned.  Only the
        clipped preview is sent — full sequences are fetched separately on
        copy (`_on_copy_request`), so click latency stays flat regardless
        of match size.
        """
        from dot_explorer.alignment_view import aligned_text

        info = input.match_select()
        ctx = _match_context(info)
        reply = {k: (info or {}).get(k) for k in ('row', 'col', 'layer', 'idx')}
        if ctx is None:
            reply['error'] = 'no-sequences'
            await session.send_custom_message('rd_match_seq', reply)
            return
        rec = ctx['rec']
        # Full sequences are available for on-demand copy in either branch.
        reply['copy'] = True
        if rec is not None and rec.cigar is not None:
            view = aligned_text(rec, ctx['qseq'], ctx['tseq'])
            reply['aligned'] = True
            reply['text'] = (
                f'{rec.query_name}:{rec.query_start:,}-{rec.query_end:,} '
                f'({rec.strand}) vs '
                f'{rec.target_name}:{rec.target_start:,}-{rec.target_end:,}'
                f'\n\n{view["text"]}'
            )
        else:
            reply['aligned'] = False
            reply['text'] = (
                f'{ctx["q"]}:{ctx["gqs"]:,}-{ctx["gqe"]:,} ({ctx["strand"]}) '
                f'vs {ctx["t"]}:{ctx["ts"]:,}-{ctx["te"]:,}\n'
                'No CIGAR for this match — sequences shown unaligned; run '
                'minimap2 with base-level alignment (-c) for a gapped '
                'alignment view.\n\n'
                f'>query {ctx["q"]}:{ctx["gqs"]:,}-{ctx["gqe"]:,} '
                f'({ctx["strand"]})\n'
                f'{_preview_slice(ctx["qseq"], ctx["gqs"], ctx["gqe"], ctx["strand"] == "-")}\n'
                f'>target {ctx["t"]}:{ctx["ts"]:,}-{ctx["te"]:,}\n'
                f'{_preview_slice(ctx["tseq"], ctx["ts"], ctx["te"], False)}'
            )
        await session.send_custom_message('rd_match_seq', reply)

    @reactive.effect
    @reactive.event(input.copy_request)
    async def _on_copy_request():
        """Serve one full match sequence for a copy-button press.

        Shipping megabases of sequence with every click made match details
        slow; instead the detail reply carries only the preview and the
        full sequence crosses the boundary once, here, when actually
        requested.  The reply echoes the segment key plus ``side`` so the
        report caches it client-side (repeat copies are then instant).
        """
        from dot_explorer.alignment_view import revcomp

        info = input.copy_request()
        side = (info or {}).get('side')
        reply = {k: (info or {}).get(k) for k in ('row', 'col', 'layer', 'idx')}
        reply['side'] = side
        ctx = _match_context(info)
        if ctx is None or side not in ('query', 'target'):
            reply['error'] = 'no-sequences'
        elif side == 'query':
            seq = ctx['qseq'][ctx['gqs'] : ctx['gqe']]
            reply['seq'] = revcomp(seq) if ctx['strand'] == '-' else seq
        else:
            reply['seq'] = ctx['tseq'][ctx['ts'] : ctx['te']]
        await session.send_custom_message('rd_copy_seq', reply)

    # --- end W2 --------------------------------------------------------------

    @render.ui
    def status():
        if boot_error():
            return ui.div(boot_error(), class_='de-status de-status-error')
        if not ready():
            msg = (
                'Loading dot-explorer (first visit compiles the WASM runtime — '
                'this can take a few seconds)…'
                if sys.platform == 'emscripten'
                else 'Loading dot-explorer…'
            )
            return ui.div(msg, class_='de-status')
        if result() is None:
            if sys.platform == 'emscripten':
                msg = (
                    'Upload two assemblies (or a PAF file) and press '
                    '"Run comparison". Everything runs in your browser — '
                    'files are never uploaded to a server. Large genomes are '
                    'limited by browser memory (the Python heap can grow to '
                    '≈4 GB); bacterial/fungal-scale assemblies work best.'
                )
            else:
                msg = (
                    'Upload two assemblies (or a PAF file) and press '
                    '"Run comparison". This local app has no upload size '
                    'limits — memory is bounded by your machine. The '
                    'minimap2/nucmer aligners still run in your browser tab '
                    '(fetched from the biowasm CDN), so they need network '
                    'access and large inputs can still strain the tab.'
                )
            return ui.div(msg, class_='de-status')
        return None

    @output(suspend_when_hidden=False)
    @render.text
    def app_memory():
        # suspend_when_hidden=False: the fixed note is CSS-hidden while its
        # text is empty, and Shiny would otherwise never render into a
        # hidden output — leaving it permanently empty (and hidden).
        # Refresh every 5 s.  The wasm linear memory only ever grows, so
        # this reports the high-water mark of the Python runtime's heap —
        # the constrained resource (Pyodide 0.27 can grow it to ~4 GB;
        # a 90 Mb k-mer run was observed at 2.9 GB before the ceiling).
        # Native runs show resident set size instead (no fixed ceiling).
        reactive.invalidate_later(5)
        if sys.platform != 'emscripten':
            try:
                import resource  # noqa: PLC0415 - unavailable on some hosts

                rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                # ru_maxrss is bytes on macOS, KiB on Linux.
                if sys.platform != 'darwin':
                    rss *= 1024
                return f'App memory: {rss / 1048576:.0f} MB peak RSS'
            except Exception:  # noqa: BLE001 - readout is best-effort only
                return ''
        try:
            import pyodide_js  # noqa: PLC0415 - pyodide-only module

            used = int(pyodide_js._module.HEAPU8.length)
        except Exception:  # noqa: BLE001 - internals unavailable -> hide
            return ''
        return f'App memory: {used / 1048576:.0f} MB of ~4096 MB wasm heap'

    @output(suspend_when_hidden=False)
    @render.text
    def result_kind():
        # Hidden output driving panel_conditional visibility (e.g. the
        # identity-colouring controls, shown only for tool/PAF results).
        # suspend_when_hidden=False: the output lives in a display:none div,
        # which Shiny would otherwise never update.
        res = result()
        return res[0] if res else ''

    @output(suspend_when_hidden=False)
    @render.text
    def gff_mode():
        """Visibility state for the static annotation toggles.

        ``''`` no annotations · ``'plain'`` cross-comparison ·
        ``'self'`` at least one self-comparison panel.  Conditions on this
        output test *positively* (``=== 'self'``), because before the first
        report it is ``undefined`` client-side and a negative test would
        flash the controls on a fresh page.

        Checks ``ann_sources`` rather than ``gff_type_index()`` — same truth
        value, without re-running the shared-colour assignment on each run.
        """
        if not any(gff_raw_for(role) is not None for role in ANNOTATION_ROLES):
            return ''
        return 'self' if self_panels() else 'plain'

    def _has_sequences(res) -> bool:
        """Whether a reordered-FASTA export is possible for this result."""
        kind, _obj, meta = res
        return (
            kind == 'kmer'
            or isinstance(meta.get('query'), SequenceProvider)
            or bool(input.paf_query_fasta())
        )

    @render.ui
    def downloads():
        res = result()
        if res is None:
            return ui.div(
                ui.tags.button(
                    'Plot (SVG)', class_='btn de-dl-disabled', disabled=True
                ),
                ui.tags.button(
                    'Plot (PDF)', class_='btn de-dl-disabled', disabled=True
                ),
                ui.tags.button(
                    'Alignment (PAF)', class_='btn de-dl-disabled', disabled=True
                ),
                ui.tags.button(
                    'Reordered query (FASTA)',
                    class_='btn de-dl-disabled',
                    disabled=True,
                ),
                ui.div('Run a comparison first.', class_='de-dl-note'),
            )
        parts = [
            ui.download_button('dl_svg', 'Plot (SVG)'),
            ui.download_button('dl_pdf', 'Plot (PDF)'),
            ui.download_button('dl_paf', 'Alignment (PAF)'),
        ]
        if _has_sequences(res):
            parts.append(ui.download_button('dl_fasta', 'Reordered query (FASTA)'))
        else:
            parts += [
                ui.tags.button(
                    'Reordered query (FASTA)',
                    class_='btn de-dl-disabled',
                    disabled=True,
                ),
                ui.div(
                    'FASTA export needs sequences — upload the query '
                    'assembly in the sidebar.',
                    class_='de-dl-note',
                ),
            ]
        if sim_matrix() is not None:
            # Clustering outputs: the matrix CSV exports whichever view the
            # Matrix tab currently shows (similarity/containment/coverage).
            parts += [
                ui.download_button('dl_heatmap_svg', 'Heatmap (SVG)'),
                ui.download_button('dl_heatmap_png', 'Heatmap (PNG)'),
                ui.download_button('dl_matrix_csv', 'Pairwise matrix (CSV)'),
            ]
        return ui.div(*parts)

    @render.plot
    def dotplot():
        # --- W2: interactive plot --- (static fallback; honours drill-down)
        res = result()
        req(res)
        with ui.Progress(min=0, max=2) as progress:
            progress.set(1, message='Computing layout & rendering plot…')
            fig = make_figure(res, config(), layout(), pair=focus())
            progress.set(2, message='Done')
        return fig

    def _figure_bytes(fmt: str) -> bytes:
        # --- W2: interactive plot --- (exports the currently shown view)
        res = result()
        req(res)
        import matplotlib.pyplot as plt

        fig = make_figure(res, config(), layout(), pair=focus())
        buf = io.BytesIO()
        fig.savefig(buf, format=fmt, bbox_inches='tight')
        plt.close(fig)
        return buf.getvalue()

    @render.download_button(filename='dotplot.svg')
    def dl_svg():
        yield _figure_bytes('svg')

    @render.download_button(filename='dotplot.pdf')
    def dl_pdf():
        yield _figure_bytes('pdf')

    @render.download_button(filename='alignment.paf')
    def dl_paf():
        res = result()
        req(res)
        kind, obj, _meta = res
        if kind == 'kmer':
            lines = obj.get_paf(group_pairs=[(QUERY_GROUP, TARGET_GROUP)], merge=True)
            yield '\n'.join(lines) + '\n'
        else:
            yield paf_text_from_alignment(obj)

    @render.ui
    def cluster_table():
        clusters = cluster_result()
        if clusters is None:
            return None
        prov = query_provider()
        lengths = dict(prov.lengths()) if prov is not None else {}
        rows = cluster_table_rows(clusters, lengths, sim=sim_matrix())
        body = [
            ui.tags.tr(
                ui.tags.td(row['cluster']),
                ui.tags.td(row['contig']),
                ui.tags.td(f'{row["length"]:,}'),
                ui.tags.td(str(row['members'])),
                ui.tags.td('' if row['mean_sim'] is None else f'{row["mean_sim"]:.3f}'),
                class_='de-cluster-row',
                data_cluster=row['cluster'],
            )
            for row in rows
        ]
        headers = ('Cluster', 'Contig', 'Length (bp)', 'Members', 'Mean similarity')
        settings = cluster_settings() or {}
        cov_label = (
            'alignment coverage'
            if settings.get('coverage_source') == 'alignment'
            else 'containment'
        )
        mode = (
            f'similarity ≥ {settings.get("cutoff", 0.8):.2f}'
            if clusters.mode == 'similarity'
            else (
                f'ANI ≥ {settings.get("identity_cutoff", 0.8):.2f} and '
                f'{cov_label} ≥ {settings.get("coverage_cutoff", 0.8):.2f}'
                + (' (reciprocal)' if clusters.reciprocal else '')
            )
        )
        return ui.div(
            ui.tags.table(
                ui.tags.thead(ui.tags.tr(*[ui.tags.th(h) for h in headers])),
                ui.tags.tbody(*body),
                class_='de-cluster-table',
            ),
            ui.div(
                f'{len(clusters.clusters)} cluster(s) at {mode}. '
                'Click rows to highlight their clusters in the plot; '
                'click again to clear.',
                class_='de-dl-note',
            ),
        )

    def _heatmap_figure():
        """Build the similarity heatmap from the applied settings."""
        sim = sim_matrix()
        req(sim)
        from dot_explorer import plot_similarity_heatmap  # noqa: PLC0415

        settings = cluster_settings() or {}
        clusters = cluster_result() if settings.get('borders_on', True) else None
        return plot_similarity_heatmap(
            sim,
            tree=active_tree(),
            clusters=clusters,
            cutoff=tree_cutoff_distance(),
            cmap=settings.get('cmap', 'viridis'),
            annotate=settings.get('cell_values', False),
        )

    @render.image(delete_file=True)
    def heatmap_plot():
        # Rendered as an image at the figure's own geometry: render.plot
        # would stretch the figure to the container, crushing the name
        # gutter and the square cells the layout reserves.
        import matplotlib.pyplot as plt  # noqa: PLC0415

        fig = _heatmap_figure()
        width_in = fig.get_size_inches()[0]
        handle = tempfile.NamedTemporaryFile(suffix='.png', delete=False)
        fig.savefig(handle.name, dpi=192)
        plt.close(fig)
        return {
            'src': handle.name,
            'width': f'{int(width_in * 96)}px',
            'alt': 'Pairwise similarity heatmap',
        }

    def _heatmap_bytes(fmt: str) -> bytes:
        import matplotlib.pyplot as plt  # noqa: PLC0415

        fig = _heatmap_figure()
        buf = io.BytesIO()
        fig.savefig(buf, format=fmt, bbox_inches='tight', dpi=200)
        plt.close(fig)
        return buf.getvalue()

    @render.download_button(filename='similarity_heatmap.svg')
    def dl_heatmap_svg():
        yield _heatmap_bytes('svg')

    @render.download_button(filename='similarity_heatmap.png')
    def dl_heatmap_png():
        yield _heatmap_bytes('png')

    @render.download_button(filename='cluster_assignments.csv')
    def dl_clusters_csv():
        clusters = cluster_result()
        req(clusters)
        import csv  # noqa: PLC0415

        buf = io.StringIO()
        writer = csv.writer(buf)
        writer.writerow(['contig', 'cluster'])
        for contig, name in clusters.assignments.items():
            writer.writerow([contig, name])
        yield buf.getvalue()

    @reactive.calc
    def displayed_matrix():
        """Return the matrix the Matrix tab currently shows, or None.

        The view radio is tab-local display state (not staged behind
        Apply): switching views only reads cached calcs.
        """
        view = input.matrix_view() if 'matrix_view' in input else 'similarity'
        if view == 'containment':
            return containment_matrix()
        if view == 'coverage':
            return alignment_coverage()
        return sim_matrix()

    _MATRIX_BLURBS = {
        'jaccard': (
            'Jaccard similarity: the fraction of distinct k-mers the two '
            'contigs share (intersection over union of their sketches). '
            'Symmetric; 1 = identical k-mer content.'
        ),
        'angular': (
            'Angular similarity: cosine-style similarity of the two '
            'sketches weighted by k-mer abundance, so repeat copy-number '
            'differences lower the score. Symmetric.'
        ),
        'ani': (
            'ANI: average nucleotide identity estimated from max '
            'containment under a Poisson mutation model. Cells show the '
            'point estimate with its 95% confidence interval; lower the '
            'sketch scaled factor to tighten the intervals. Pairs sharing '
            'no hashes show 0.'
        ),
        'containment': (
            'sourmash containment: the fraction of the ROW contig’s '
            'k-mers found in the COLUMN contig. Asymmetric — a fragment '
            'is fully contained in its parent, not vice versa — and '
            'depressed by SNPs (≈ coverage × ANI^k).'
        ),
        'aln_coverage': (
            'Alignment block coverage: the fraction of the ROW contig '
            'covered by the union of its alignment blocks against the '
            'COLUMN contig, from a clean minimap2 run (no -P; secondary '
            'alignments excluded). Asymmetric; robust to SNPs.'
        ),
    }

    @render.ui
    def matrix_table():
        matrix = displayed_matrix()
        if matrix is None:
            return ui.div(
                'No matrix for this view yet — alignment coverage needs a '
                'completed clean minimap2 run (computed in the background '
                'when coverage clustering is applied); containment needs '
                'clustering enabled.',
                class_='de-dl-note',
            )
        settings = cluster_settings() or {}
        show_ci = matrix.metric == 'ani' and matrix.ci_low is not None
        blurb = _MATRIX_BLURBS.get(matrix.metric, matrix.metric)
        sketch_note = ''
        if matrix.metric != 'aln_coverage':
            sketch_note = (
                f' Sketches: k={settings.get("ksize", 21)}, '
                f'scaled={settings.get("scaled", 1000)}, abundance '
                f'{"on" if settings.get("abund", True) else "off"}.'
            )
        asym_note = ''
        if matrix.metric in ('containment', 'aln_coverage'):
            asym_note = (
                ' Read row-wise: each value is the fraction of the row '
                'contig accounted for by the column contig, so the matrix '
                'is not symmetric.'
            )

        def cell(i: int, j: int) -> str:
            value = matrix.values[i, j]
            low, high = None, None
            if show_ci and i != j:
                low, high = matrix.ci_low[i, j], matrix.ci_high[i, j]
            # Collapsed bounds mean no CI was available for the pair
            # (size-inaccurate sketch): show just the point estimate.
            if low is not None and not (low == value == high):
                return f'{value:.4f} ({low:.4f}–{high:.4f})'
            return f'{value:.4f}'

        # Cells carry the heatmap palette, so the table and Heatmap tab
        # read the same way; text flips light/dark on cell luminance.

        colormap = matplotlib.colormaps[settings.get('cmap', 'viridis')]

        def cell_style(value: float) -> str:
            r, g, b, _a = colormap(min(1.0, max(0.0, float(value))))
            luminance = 0.299 * r + 0.587 * g + 0.114 * b
            fg = '#000000' if luminance > 0.5 else '#ffffff'
            return (
                f'background-color: rgb({int(r * 255)}, {int(g * 255)}, '
                f'{int(b * 255)}); color: {fg};'
            )

        n = len(matrix.names)
        header = ui.tags.tr(
            ui.tags.th(''), *[ui.tags.th(name) for name in matrix.names]
        )
        body = [
            ui.tags.tr(
                ui.tags.th(matrix.names[i]),
                *[
                    ui.tags.td(cell(i, j), style=cell_style(matrix.values[i, j]))
                    for j in range(n)
                ],
            )
            for i in range(n)
        ]
        return ui.div(
            ui.div(blurb + sketch_note + asym_note, class_='de-matrix-blurb'),
            ui.div(
                ui.tags.table(
                    ui.tags.thead(header),
                    ui.tags.tbody(*body),
                    class_='de-cluster-table de-matrix-table',
                ),
                class_='de-matrix-scroll',
            ),
        )

    @render.download_button(filename='pairwise_matrix.csv')
    def dl_matrix_csv():
        matrix = displayed_matrix()
        req(matrix)
        import csv  # noqa: PLC0415

        buf = io.StringIO()
        buf.write(f'# metric: {matrix.metric}\n')
        writer = csv.writer(buf)
        writer.writerow([''] + matrix.names)
        for name, row in zip(matrix.names, matrix.values):
            writer.writerow([name] + [f'{value:.6g}' for value in row])
        if matrix.metric == 'ani' and matrix.ci_low is not None:
            for label, bounds in (
                ('ani_ci_low', matrix.ci_low),
                ('ani_ci_high', matrix.ci_high),
            ):
                buf.write(f'# {label}\n')
                writer.writerow([''] + matrix.names)
                for name, row in zip(matrix.names, bounds):
                    writer.writerow([name] + [f'{value:.6g}' for value in row])
        yield buf.getvalue()

    @render.download_button(filename='query_reordered.fasta')
    def dl_fasta():
        # --- W2: interactive plot --- (uses the same explicit layout as
        # the plot, so the exported FASTA matches what is displayed)
        res = result()
        req(res)
        kind, obj, meta = res
        lay = layout()
        # Contigs the length filter left out of the plot still belong in the
        # export — they are simply unordered, appended after the plotted
        # ones.  CrossIndex.write_fasta writes exactly the names it is
        # given, so omitting them here would quietly drop sequence.
        order = lay['query_names'] + lay['excluded_query']
        reverse = lay['reverse']
        if kind == 'kmer':
            with tempfile.TemporaryDirectory() as tmp:
                out = Path(tmp) / 'query_reordered.fasta'
                obj.write_fasta(out, group=QUERY_GROUP, order=order, reverse=reverse)
                yield out.read_bytes()
        else:
            # Coordinate-only alignment (PAF import / external tool): export
            # from the query assembly sequences when they are available —
            # either attached to the result or from the sidebar upload.
            query = meta.get('query')
            if not isinstance(query, SequenceProvider):
                try:
                    query = _parse_upload(input.paf_query_fasta, 'query')
                except ValueError:
                    query = None
            if query is None:
                ui.notification_show(
                    'Reordered FASTA export needs sequences — upload the '
                    'query assembly in the sidebar (PAF files carry '
                    'coordinates only).',
                    type='warning',
                    duration=8,
                )
                req(False)
            yield reordered_fasta_text(list(query.iter_records()), order, reverse)


app = App(app_ui, server, static_assets=APP_DIR / 'www')
