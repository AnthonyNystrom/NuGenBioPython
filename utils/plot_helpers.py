"""Shared matplotlib plot styling + SVG output helpers.

All diagram-producing routes should go through these helpers so that:
- Colors/fonts stay consistent across the app (matches the genomediagram
  renderer in routes/genomediagram_routes.py).
- Output is vector SVG by default — browsers scale it cleanly and the
  response stays small for text-heavy diagrams.
- Figure sizing is done in one place per diagram type, not copy-pasted.
"""
import base64
import math
from io import StringIO, BytesIO

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# Curated palette — muted, print-safe, colorblind-friendlier.
PALETTE = {
    'blue':       '#2563eb',
    'red':        '#dc2626',
    'green':      '#059669',
    'orange':     '#ea580c',
    'purple':     '#7c3aed',
    'brown':      '#92400e',
    'pink':       '#db2777',
    'cyan':       '#0891b2',
    'gray':       '#64748b',
    'grey':       '#64748b',
    'darkblue':   '#1e3a8a',
    'teal':       '#0d9488',
    'black':      '#0f172a',
    'yellow':     '#ca8a04',
    'lightblue':  '#93c5fd',
    'lightgreen': '#86efac',
    'lightcoral': '#fca5a5',
    'lightgray':  '#cbd5e1',
}

# Role colors for network diagrams
ROLE_COLORS = {
    'source':       '#10b981',  # green
    'sink':         '#ef4444',  # red
    'intermediate': '#3b82f6',  # blue
    'catalyst':     '#a855f7',  # purple
}

EDGE_COLOR = '#1e293b'
TITLE_COLOR = '#0f172a'
LABEL_COLOR = '#334155'
MUTED_COLOR = '#64748b'
AXIS_COLOR = '#94a3b8'
GRID_COLOR = '#e2e8f0'


def resolve_color(c, default='#2563eb'):
    if not c:
        return default
    c = str(c).strip()
    if c.startswith('#'):
        return c
    return PALETTE.get(c.lower(), default)


def figure_to_svg_data_url(fig, pad_inches=0.2):
    """Serialize a matplotlib Figure to an svg+xml base64 data URL.

    Callers should not call plt.close() themselves — this helper does.
    """
    buf = StringIO()
    fig.savefig(buf, format='svg', bbox_inches='tight', facecolor='white',
                pad_inches=pad_inches)
    plt.close(fig)
    svg = buf.getvalue().encode('utf-8')
    b64 = base64.b64encode(svg).decode()
    return f'data:image/svg+xml;base64,{b64}'


def figure_to_png_data_url(fig, dpi=150, pad_inches=0.2):
    """PNG fallback for the rare cases where SVG is not appropriate
    (e.g. when the backend draws pixels directly)."""
    buf = BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight', facecolor='white',
                dpi=dpi, pad_inches=pad_inches)
    plt.close(fig)
    buf.seek(0)
    b64 = base64.b64encode(buf.getvalue()).decode()
    return f'data:image/png;base64,{b64}'


def style_axes(ax, hide=('top', 'right'), spine_color=AXIS_COLOR,
               tick_color=MUTED_COLOR, grid=False):
    """Apply consistent axis styling: hide spines, mute tick labels."""
    for s in hide:
        ax.spines[s].set_visible(False)
    for s in ax.spines:
        if s not in hide:
            ax.spines[s].set_color(spine_color)
            ax.spines[s].set_linewidth(0.6)
    ax.tick_params(axis='both', labelsize=8, colors=tick_color,
                   length=3, pad=2)
    if grid:
        ax.grid(True, color=GRID_COLOR, linewidth=0.5, zorder=0)


def set_title(ax_or_fig, text, **kw):
    """Set a title with the app's consistent styling."""
    defaults = dict(fontsize=12, fontweight='600', color=TITLE_COLOR, pad=10)
    if hasattr(ax_or_fig, 'suptitle') and not hasattr(ax_or_fig, 'set_title'):
        defaults = dict(fontsize=13, fontweight='600', color=TITLE_COLOR, y=0.98)
    defaults.update(kw)
    if hasattr(ax_or_fig, 'set_title'):
        ax_or_fig.set_title(text, **defaults)
    else:
        ax_or_fig.suptitle(text, **defaults)


def fmt_bp(value):
    """Format a base-pair position as kb/Mb when large."""
    v = int(round(value))
    if abs(v) >= 1_000_000:
        s = f'{v/1_000_000:.2f}'.rstrip('0').rstrip('.')
        return f'{s} Mb'
    if abs(v) >= 10_000:
        if v % 1000 == 0:
            return f'{v//1000} kb'
        return f'{v/1000:.1f} kb'
    return f'{v:,}'


def nice_ticks(extent, target_count=8):
    """Choose tick positions at a round step covering [0, extent]."""
    if extent <= 0:
        return [0]
    raw = extent / target_count
    magnitude = 10 ** math.floor(math.log10(raw)) if raw > 0 else 1
    step = magnitude
    for m in (1, 2, 2.5, 5, 10):
        step = m * magnitude
        if extent / step <= target_count * 1.5:
            break
    ticks = []
    v = 0
    while v <= extent + 1e-6:
        ticks.append(int(round(v)))
        v += step
    return ticks
