"""Interactive HTML report generation for dot-explorer dotplots.

This private package renders a :class:`~dot_explorer.dotplot.DotPlotter` figure
as a single self-contained HTML file: the matplotlib figure is embedded as
inline SVG, match coordinates (and optionally sequences) are embedded as a
JSON payload, and vanilla CSS/JS assets bundled with the package provide
panel selection, scroll zoom and click-to-inspect behaviour.
"""

from dot_explorer._html.render import render_html_report
from dot_explorer._html.serialize import build_panel_payload

__all__ = ['build_panel_payload', 'render_html_report']
