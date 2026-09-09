# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-09-09

### Added

- The browser app now ships inside the wheel, so `pip install
  "dot-explorer[app]"` gives you a runnable app with no repository clone and no
  Rust toolchain. Launch it with the new **`dot-explorer-app`** console script
  (`--host`, `--port`, `--no-browser`); it serves on
  <http://127.0.0.1:8000> by default. Previously the `app` extra installed the
  app's dependencies but not the app itself.

### Changed

- `pyfaidx` moved from the `app` extra into the core dependencies — it is pure
  Python (~0.2 MB), free-threading safe, and backs the lazy sequence provider.
  `shiny` and the `cluster` extra stay optional: `shiny` pulls `orjson`, which
  has no free-threaded wheel, so requiring it would make the package
  uninstallable on the 3.14t interpreter we publish wheels for.
- **Contributor-facing move:** the app source now lives at
  `python/dot_explorer/app/` instead of `app/` at the repository root. It ships
  as package data rather than an importable subpackage, so `app.py` keeps its
  flat `core.*` imports (which is what `shinylive export` requires). Run it
  with `dot-explorer-app`, or `shiny run python/dot_explorer/app/app.py`.
- The Shinylive export runs through the new `scripts/build_shinylive.py`, which
  stages the app plus the wasm wheel under `build/` and exports that — keeping
  `.whl` files out of the package tree, where they would be nested inside the
  built wheel. It also fails loudly on a missing or ambiguous wasm wheel.
- Documented where the app writes temporary files (the system temp directory,
  never the install or launch directory) and how to redirect it with `TMPDIR`
  on clusters with a small `/tmp`.
- `scripts/stamp_version.py` now writes `python/dot_explorer/_version.py` only
  when HEAD sits exactly on a release tag, as it already did for
  `CITATION.cff`, and no longer emits a `.postN` development suffix. Between
  tags every source keeps the last released version, so
  `dot_explorer.__version__` and `importlib.metadata.version('dot-explorer')`
  agree (they previously diverged on dev commits, since the wheel version comes
  from `Cargo.toml`, which never carried the suffix) and an ordinary commit no
  longer dirties the working tree.
- Development builds identify the commit they came from: between release tags
  the stamper writes an untracked `python/dot_explorer/_version_local.py`
  holding a PEP 440 local version — `X.Y.Z+<short hash>`, with `.dirty`
  appended when the working tree has uncommitted changes — and `_version.py`
  imports it when present. `dot_explorer.__version__` therefore reads e.g.
  `0.1.0+1a2b3c4` in a development build while the distribution metadata stays
  `0.1.0`, which is exactly what a PEP 440 local version means. The file is
  gitignored and deleted at a release tag, so it can never reach a release
  wheel; no tracked file changes between tags.
- The Colab setup cells in the quickstart and minimap2 tutorials install
  `dot-explorer` from PyPI instead of building from the GitHub repository.

## [0.1.0] - 2026-09-09

First published release.

### Added — library

- Core engine: rolling-hash ntHash k-mer indexing, both-strand k-mer
  matching, strand-aware run merging, PAF output, `CrossIndex`
  cross-assembly comparisons, collinearity contig ordering with
  reordered/reoriented FASTA export, and matplotlib dotplot visualisation
  via `DotPlotter`.
- Similarity, clustering and trees (new optional `cluster` extra —
  `pip install "dot-explorer[cluster]"` for sourmash + scipy):
  - `dot_explorer.tree`: dependency-free newick / IQ-TREE `.treefile` parser
    (`Tree.from_newick`/`Tree.read`, quoted labels, support values,
    comments), `Tree.from_linkage` for scipy linkage matrices, bidirectional
    tip-vs-sequence-name validation, and `draw_tree` (rectangular
    dendrogram with cutoff line and scale bar).
  - `dot_explorer.similarity`: sourmash FracMinHash sketching
    (`compute_sketches`, `SketchParams` — k=21, scaled=1000, abundance on
    by default) and all-vs-all `pairwise_similarity` with `jaccard`,
    abundance-weighted `angular`, `ani` (with 95% confidence bounds),
    asymmetric `containment` and `max`/`avg_containment` metrics;
    `linkage_from_similarity`, similarity-cutoff `assign_clusters`, and
    dual identity+coverage `assign_clusters_dual` (reciprocal or one-way
    coverage); `SimilarityMatrix.to_csv` / `ClusterResult.to_csv` exports.
  - `plot_similarity_heatmap`: pairwise similarity heatmap with the tree on
    the y-axis, cluster block outlines, selectable colormap and a colour
    scale bar.
  - `DotPlotter.plot` gained `tree=` (dendrogram drawn left of the matrix,
    contig order fixed to the tree's leaf order, names moved onto the tree
    tips; incompatible with `contig_order`/`auto_reverse`), `tree_cutoff=`,
    `tree_scalebar=`, and `cluster_borders=` (bold outline around each
    cluster's block of panels); the HTML report payload carries the cluster
    blocks under a `clusters` key.
  - Docs: a "Similarity & Clustering" metric-selection guide, a
    "Clustering & Trees" tutorial notebook, and API reference pages.
- `DotPlotter.plot` / `plot_single` gained `cap_style`
  (`'butt'`/`'round'`/`'projecting'`, default `'projecting'`), setting the
  line cap on every match layer. Butt caps draw a match shorter than
  *dot_size* wider across its diagonal than along it, so it reads as a mark
  rotated 90°; square and round caps extend the stroke half a line width
  past each endpoint and keep it on its own diagonal.
- HTML report: a **FASTA header** toggle in the match detail bar prepends a
  header carrying the match's names, coordinates, strand, length and (on the
  identity layer) percent identity to copied query/target sequences. The
  setting persists across match popups — via `localStorage` in standalone
  reports, and via the embedding app across plot re-renders.

### Added — app

- Trees & clustering (self-alignment mode): upload a newick / IQ-TREE
  `.treefile` to fix the contig order and draw the tree left of the matrix
  (tip/sequence-name mismatches error with both directions listed), or
  compute a hierarchical clustering tree from sourmash similarity
  (Jaccard / angular / ANI; user-set k, scaled, abundance) — the tree
  overrides the contig-order and auto-flip options while active. Cluster
  assignment by similarity cutoff (with an optional dashed cutoff line
  through the dendrogram) or combined ANI + containment thresholds with a
  reciprocal-coverage toggle; a Clusters tab whose rows highlight their
  cluster blocks in the plot (non-members dim); bold cluster outlines; a
  similarity Heatmap tab with a selectable palette; and CSV downloads of
  the matrix and assignments. sourmash + scipy stay out of the first-load
  bundle — under Pyodide they micropip-install (~25 MB, one-time) when
  clustering is first enabled; natively install `dot-explorer[cluster]`.
- Clustering refinements: an **Apply changes** button gates every Trees &
  clustering setting (nothing recomputes or redraws until clicked); a
  **Matrix** tab shows the pairwise matrix with row/column names, a
  metric explanation, ANI 95% CIs per cell, asymmetric containment and
  alignment-coverage views, and the matrix CSV download (moved from the
  Clusters tab); the Heatmap tab gained in-cell values and SVG/PNG
  downloads; a **coverage source** option lets the identity+coverage
  mode use alignment block coverage from the current run (SNP-robust)
  instead of sourmash containment.
- Heatmap figure polish: square cells, sequence names placed between the
  tree and the heatmap (long-name safe, no tree-tip overlap), cluster
  outlines now white, dashed and heavier; the Heatmap tab gained the same
  navigation controls as the main plot (scroll to pan, Shift+scroll
  sideways, Cmd/Ctrl+scroll to zoom, drag to zoom to a region,
  double-click or Esc to reset; the figure opens fitted to the pane).
  The Apply-changes button moved to the bottom of the Trees & clustering
  section; the heatmap and pairwise-matrix downloads live in the sidebar
  Downloads section; the Matrix tab's cells are coloured with the
  selected heatmap palette. Dot-plot tree gutter widened (floored
  at 18% of the grid width) and the figure-level 'Position' label moved
  clear of the dendrogram.

- `SequenceIndex.approx_bytes()` (and `CrossIndex.approx_bytes()`): exact
  per-component heap accounting for the k-mer index (sequence bytes, CSR
  tables, pair cache), plus `scripts/mem_profile_index.py` to measure peak
  RSS and index size for synthetic pairs — the audit behind the web app's
  size limits is recorded in docs/development.md.

- Interactive single-file HTML dotplot reports (`DotPlotter.to_html()`, or an
  `.html` output path): click a panel to focus it, scroll/drag to pan and
  zoom, drag a box to zoom to a region, click a match for coordinates,
  identity and sequence details.
- Match selection in reports: a selected match stays highlighted after its
  detail bar is closed, and its query/target ranges are drawn as translucent
  bands through every panel in its grid row and column. `Shift`+drag
  box-selects every match intersecting the box (testing the drawn line, not
  its bounding box) without opening details; `Esc` clears selections first
  and only then resets/exits.
- Annotation support: GFF3 parsing upgrades (percent-decoded attributes,
  `feature_id`/`parent`/`name`, `from_text`/`from_bytes` with gzip
  detection, embedded `##FASTA`, multi-part grouping); diagonal feature
  squares drawn behind the matches with an automatic type legend; lane-packed
  side annotation tracks on focused plots (strand arrows, multi-part CDS
  connectors); clickable features in reports, with track features banding
  their row/column — bands survive into saved SVG/PNG/PDF figures.
  `GffFeature` gained `color` and `source_file`; the report payload gained a
  `tracks` key and per-segment `uid`s.
- Identity and aligned sequences: `PafRecord.identity` (gap-compressed `de`
  tag, else CIGAR-derived, else the BLAST-style estimate);
  `alignment_view.aligned_text()` for gapped query/match/target views;
  `to_line()` preserves optional SAM-style tags; identity colouring uses the
  best available metric and `plot(identity_colorbar=True)` appends a colour
  key.
- Plotting: `contig_order='length'|'colinearity'` with `auto_reverse=True`,
  `hide_internal_axes=True`, opt-in `nature_style()` styling, and a
  reordering/reorientation tutorial notebook.
- `SequenceIndex.compare_pairs_stranded`: batched both-strand comparison in
  one native call (GIL released) with a `min_block_len` filter;
  `CrossIndex.compute_matches` gains the same option and computes the whole
  grid in one batched call.
- A **Line cap** selector beside the line-width control (Square / Round /
  Flat), applied client-side inside the embedded report like the other
  display-only options — switching caps never re-renders matplotlib — and
  carried into the SVG/PDF downloads.

### Added — browser app

- Fully client-side assembly comparison (`app/`) built with Shiny for Python
  and deployed as a static Shinylive/Pyodide site alongside the docs. Upload
  FASTA/FASTA.gz assemblies, align with dot-explorer's k-mer engine or in-browser
  biowasm aligners (minimap2 2.22, nucmer/MUMmer4), reconfigure the dotplot
  without recomputing, and download SVG/PDF plots, PAF, and a reordered query
  FASTA. Includes a loading splash with staged progress, background aligner
  pre-download, and a wasm-heap memory readout.
- Input modes: PAF as a top-level mode (with contig-name validation against
  an optional query assembly), self-alignment ("Align assembly to itself",
  which clears/hides target annotations), and GenBank assemblies parsed by a
  pure-Python flat-file parser (gzip, multi-record, full location grammar) so
  sequences and annotations come from one upload.
- Annotations: GFF3/GenBank uploads per role with merged multi-source
  handling, per-feature-type visibility toggles and colour pickers gated
  behind an **Apply changes** button (pre-first-plot edits apply on render),
  self-panel shading, side tracks in the drill-down, a **Clear annotations**
  reset, and warnings when uploaded sequence names match nothing.
- Drill-down **Annotations** tab: lists every feature on both focused
  sequences (contig column included, self-panels deduplicated) with
  per-feature show/hide and colour overrides held until applied, header
  sorting, multiple per-column text filters with an Apply button, and bulk
  actions. Rows are built client-side, so switching tabs is instant. Click a
  row to select it (⌘/Ctrl-click for several); selected features are banded
  on the Plot tab, and the bands survive re-renders and exports.
- Interactive drill-down: double-click a panel for the focused single-pair
  view with contig-name axis labels and bp/Kbp/Mbp ticks; clicking an
  identity-coloured match shows the gapped CIGAR alignment (minimap2 gained
  a **Base-level alignment (`-c`)** option), other matches show their raw
  sequences; copy buttons fetch and cache full sequences on demand.
- Fullscreen plot mode: an expand button at the top-left of the plot panel
  requests true browser fullscreen (falling back to a window-filling
  overlay), scales the figure to the viewport, hides the hints, memory
  readout and tab strip, and keeps **Back to overview** working without
  leaving fullscreen. With a selection active, `Esc` deselects first; a
  second `Esc` exits.
- Manual dark-mode toggle in the header (defaults to the system preference).
- Aligner UX: per-stage progress messages surfaced at the right end of the
  header bar, a collapsible aligner log with exact command lines, in-flight
  runs cancelled on method switch or re-run, per-tool timeouts with a
  **Still aligning** wait/cancel prompt, repeat-genome options (minimap2
  `-m`/`-P`/`-D`, nucmer `--nosimplify`), presets pre-filling their real
  `-k`/`-w`/`-m` values, hoverable ⓘ help on every option, and a size guard
  that steers >80 Mb k-mer runs to minimap2 instead of exhausting the wasm
  heap.
- A **Min contig length** setting hides short contigs from the plot (the
  reordered-FASTA download keeps them).
- Cmd/Ctrl+click adds matches to (or removes them from) the selection one
  by one, in the overview and the drill-down. Selections — matches and the
  selected annotation feature — persist across overview <-> drill-down
  swaps: the app holds them in view-independent terms (names + data
  coordinates) and re-applies them to each rebuilt report, keeping
  selections on pairs a narrower view does not show.
- Upload-size guidance, from empirically measured ceilings (synthetic
  pairs in Chrome): beyond ~80 Mb combined the k-mer method is removed
  from the selector (90 Mb completes at a 2.9 GB heap peak; 100 Mb
  silently kills the Python runtime), with a pointer to the tutorial
  notebooks for running the library locally; above ~200 MB the app
  suggests uploading a precomputed PAF instead (minimap2 completed at
  200 MB, nucmer at 250 MB; 300 MB crashed the browser tab for both).

### Changed

- Match segments are drawn with projecting (square) line caps by default,
  matching the dendrogram branches. Segments longer than the line width are
  visually unchanged (each end grows by half a line width — 0.25 pt at the
  default `dot_size`); sub-linewidth matches stop reading as marks rotated
  90°. Pass `cap_style='butt'` to restore the previous rendering.
- Alignment-based coverage for identity+coverage clustering is now always
  computed from a minimap2 run **without** `-P`, instead of from whichever
  alignment is on screen. Retaining every chain re-covers repeat features
  and inflates the covered span, and nucmer / k-mer / imported-PAF records
  are no longer used for coverage at all. When the displayed result is not
  a clean minimap2 run, the app computes one in the background (cached by
  input digest) while the plot keeps showing the user's chosen alignment;
  if that run cannot complete, clustering falls back to sourmash
  containment and says so. Records tagged as secondary (`tp:A:S`) are
  excluded from the coverage union.
- Versioning is now derived from the latest `v*` git tag, replacing four
  hand-maintained version strings. `scripts/stamp_version.py` (run by a
  `stamp-version` pre-commit hook and by CI before every build) stamps
  `Cargo.toml`, `Cargo.lock` and the new `python/dot_explorer/_version.py`,
  and — only at an exact release tag — `CITATION.cff`; `pyproject.toml`
  declares `dynamic = ["version"]`, which maturin resolves from
  `Cargo.toml`. Wheel-building jobs check out full history with tags, and
  CI fails on drift (`stamp_version.py --check`).
- HTML reports no longer embed match sequences by default: `to_html()` /
  `plot()` gained `embed_sequences` (default `False`), keeping exported
  files small with large alignment sets. Pass `embed_sequences=True` to
  restore the standalone sequence preview/copy buttons (subject to the
  existing ~2 Mb residue cap); a standalone report without embedded
  sequences now shows a notice instead of an empty detail pane. The match
  detail preview is clipped at 1,000 bases/columns (was 20,000).
- App: status and warning copy is now platform-aware — the local Shiny app
  no longer claims "everything runs in your browser" or a 4 GB heap limit,
  and states that it has no upload size limits (the k-mer size gate and
  wasm-heap messaging remain browser-only; biowasm aligner warnings apply on
  both platforms since those tools run in the browser tab either way).
- App: uploaded assemblies are now served through a lazy sequence provider
  (`app/core/seqs.py`) backed by a pyfaidx index over the uploaded file
  instead of resident Python strings — previews, copies and the k-mer index
  build fetch only the windows/contigs they need. Falls back to in-memory
  parsing when pyfaidx cannot index the input. New `pip install
  "dot-explorer[app]"` extra (shiny + pyfaidx) for running the app locally or
  on an HPC (`shiny run --launch-browser app/app.py`).
- Releases publish real binary wheels: maturin builds for Linux
  (x86_64/aarch64), macOS (arm64/x86_64) and Windows (x64) across CPython
  3.12–3.14 and 3.14t, plus an sdist. The free-threaded wheel declares
  `gil_used = false` and CI fails if importing it re-enables the GIL. A
  PEP 783 `pyemscripten_2026_0_wasm32` wheel is published to PyPI for
  browser runtimes; the Shinylive app keeps its separate legacy
  `emscripten_3_1_58_wasm32` build (Pyodide 0.27.7).
- App identity: restyled around a CSS custom-property theme layer (teal
  accent, light/dark parity) with a bundled display font, a
  "dot·explorer — live assembly comparison" wordmark, a dot-grid page
  background, a translucent sidebar, and a matching loading splash. The
  progress popup, title-bar pulse and corner "Processing…" pill were
  replaced by the single header task status plus spinner.
- Feature-type colours are assigned once across both uploads
  (case-insensitive with aliases, conventional colours for common types),
  shown in a single **Feature types** list with `Q`/`T` badges;
  **Shade features on diagonal** only appears when a self-comparison panel
  exists.
- Navigation tips moved below the plot with bold action terms, tailored to
  the view; report copy buttons are labelled **Fetch …**/**Copy …** by what
  the next press does.
- Assembly-scale aligner defaults (nucmer `-l 100 -c 200`); docs build moved
  from MkDocs + Material to Zensical with pre-rendered notebooks and a CI
  docs-build job; the docs Web App page opens the app in its own window.

### Fixed

- Drag-zoom in the embedded overview lands exactly where the box was
  drawn. A figure taller than the iframe kept its intrinsic size and
  scrolled, so the zoomed viewBox filled a mostly off-screen element and
  appeared shifted up; embedded reports now fit the SVG to the iframe
  viewport (as fullscreen already did).
- The memory readout and docs no longer claim a 2 GB wasm cap: the
  Pyodide heap measurably grows to 4096 MB (~4.0 GB allocatable before a
  clean MemoryError), and the first-load download is 26 MB compressed
  (37 MB of assets), replacing the stale 30/45 MB claims.
- Sequence previews no longer go stale after re-running the aligner with
  different settings: the per-pair record cache is cleared on every new
  result, so a run without `-c` stops showing the previous run's base-level
  alignments.
- The focused-view title is centred over the dot-plot panel rather than the
  whole figure (the identity colorbar and legend shift the panel
  off-centre).
- The aligner log no longer covers the report's bottom detail bar in the
  drill-down (tab panes now form the same flex column as the overview).
- Annotations survive app startup in the browser: the wasm-wheel install and
  the annotation priming raced, caching an import error for the whole
  session.
- Panel-id collisions no longer break click-to-focus, dimming or
  double-click drill-down (axes backgrounds use their own `rd-plotbg-`
  prefix and panel ids are matched exactly).
- Feature-type toggles and colours survive "Run comparison"; annotation
  edits no longer rebuild the figure once per click; toggling six types
  costs one redraw.
- Unstranded track features no longer render as full-height blocks (the
  corner radius is now computed per axis).
- Extreme-ratio focused views keep the title, y label and ticks legible
  without distorting bp-per-inch parity; short contig row labels overhang
  instead of collapsing to an ellipsis.
- The sidebar keeps its scroll position through file-picker interactions;
  file inputs reset properly after **Clear annotations**.
- The deploy pipeline no longer bundles stale wasm wheels (build dirs are
  cleared and the app selects the wheel by the running Pyodide's platform
  tag).

### Performance

- The k-mer index was restructured into a canonical-hash CSR table
  (~13–16 bytes/bp vs ~100) with a sorted two-pointer match walk: on a real
  37 Mb × 37 Mb fungal pair, index build 6.4 s → 3.4 s, matching
  20 s → 3.5 s, peak memory roughly quartered — and the k-mer method now
  completes in the browser instead of exhausting the 2 GB wasm heap.
  `reorder_for_colinearity` reuses cached match records (72 s → 11 s).
- Browser app: `parse_fasta_bytes` rewritten single-pass (~7.5× faster);
  line width and min match length apply instantly inside the report over
  postMessage instead of re-rendering matplotlib; contig orders and plotter
  construction are memoised per result; assemblies ship to the biowasm
  aligners once and are referenced by digest; session caches are bounded
  LRUs.
- Match details ship only the clipped preview (full sequences are fetched on
  demand and cached), so megabase matches never materialise full slices for
  display.

### Removed

- LASTZ as an in-browser alignment method (its wasm build was impractically
  slow at assembly scale); precomputed LASTZ alignments still import as PAF.
