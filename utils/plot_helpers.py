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
import matplotlib.colors as mcolors
import matplotlib.text as mtext
import matplotlib.patches as mpatches
import matplotlib.collections as mcoll


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


# --------------------------------------------------------------------------
# Theme support
# --------------------------------------------------------------------------
# The app has a real dark mode (html[data-theme="dark"], see static/css/
# tokens.css), but server-rendered figures were always drawn on a hardcoded
# white ground with dark slate text — a glaring white slab in every dark-mode
# card. Figures are delivered as opaque <img src="data:..."> URLs, so CSS
# cannot fix this client-side; the theme has to reach the renderer.
#
# Rather than thread a `theme` argument through every render_* signature, the
# theme is resolved from the current request at serialization time and the
# figure's *chrome* colors (text, spines, ticks, grid, background) are remapped
# just before saving. Data colors from PALETTE are deliberately left alone:
# they were chosen to read on both grounds, so a chart keeps its meaning and
# its legend across themes.

DARK_BG = '#141c2e'          # --color-surface
DARK_TITLE = '#f1f5f9'       # --color-text
DARK_LABEL = '#e8ecf1'       # --slate-700 (flipped)
DARK_MUTED = '#b8c2cf'       # --color-text-muted
DARK_AXIS = '#6b7a8f'        # --slate-400 (flipped)
DARK_GRID = '#2a3548'        # --color-border
DARK_EDGE = '#cbd5e1'

# Light chrome color -> dark equivalent. Keys cover every constant this module
# exports plus the plain black/white matplotlib defaults.
_DARK_CHROME = {
    TITLE_COLOR.lower(): DARK_TITLE,
    LABEL_COLOR.lower(): DARK_LABEL,
    MUTED_COLOR.lower(): DARK_MUTED,
    AXIS_COLOR.lower():  DARK_AXIS,
    GRID_COLOR.lower():  DARK_GRID,
    EDGE_COLOR.lower():  DARK_EDGE,
    '#cbd5e1': '#3e4b64',    # baseline / light rules
    '#000000': DARK_TITLE,
    'black':   DARK_TITLE,
    '#ffffff': DARK_BG,
    'white':   DARK_BG,
}


def resolve_theme(explicit=None):
    """Return 'dark' or 'light' for the figure about to be rendered.

    Order of preference: an explicit argument, then the current request (an
    ``X-Theme`` header, a ``theme`` query arg, or a ``theme`` key in a JSON
    body), then 'light'. Safe to call outside a request context — the offline
    scripts in tests/visual call the renderers directly.
    """
    if explicit in ('dark', 'light'):
        return explicit
    try:
        from flask import request, has_request_context
        if has_request_context():
            candidate = (
                request.headers.get('X-Theme')
                or request.args.get('theme')
                or (request.get_json(silent=True) or {}).get('theme')
                or (request.form.get('theme') if request.form else None)
            )
            if candidate and str(candidate).lower() == 'dark':
                return 'dark'
    except Exception:
        # Never let theme detection break a render.
        pass
    return 'light'


def _remap(color):
    """Map one light-theme chrome color to its dark equivalent, or None.

    Normalises through matplotlib's color converter first: artist getters
    return RGBA tuples rather than the hex string the caller passed in, so a
    plain string lookup silently misses every patch and line.
    """
    if color is None:
        return None
    # A fully transparent artist reports RGBA (0, 0, 0, 0), and to_hex drops
    # the alpha channel — so an invisible fill would look like pure black and
    # get remapped to near-white, turning 'none' facecolors into solid light
    # blobs in dark mode. Leave anything transparent alone.
    try:
        rgba = mcolors.to_rgba(color)
        if rgba[3] == 0:
            return None
    except (ValueError, TypeError):
        pass
    try:
        key = mcolors.to_hex(color).lower()
    except (ValueError, TypeError):
        try:
            key = str(color).lower().strip()
        except Exception:
            return None
    return _DARK_CHROME.get(key)


def apply_dark_theme(fig):
    """Recolor a finished figure's chrome for dark backgrounds, in place.

    Walks text, spines, ticks, grid lines and patch edges. Data colors are
    untouched, so series stay matched to their legends.
    """
    fig.patch.set_facecolor(DARK_BG)

    for ax in fig.get_axes():
        ax.set_facecolor(DARK_BG)

        for spine in ax.spines.values():
            new = _remap(spine.get_edgecolor()) or DARK_AXIS
            if spine.get_visible():
                spine.set_color(new)

        ax.tick_params(axis='both', which='both',
                       colors=DARK_MUTED, labelcolor=DARK_MUTED)

        for gridline in list(ax.get_xgridlines()) + list(ax.get_ygridlines()):
            gridline.set_color(DARK_GRID)

        for line in ax.get_lines():
            new = _remap(line.get_color())
            if new:
                line.set_color(new)

        for txt in ax.texts + [ax.title, ax.xaxis.label, ax.yaxis.label]:
            new = _remap(txt.get_color())
            txt.set_color(new or DARK_LABEL)

        for lbl in ax.get_xticklabels() + ax.get_yticklabels():
            lbl.set_color(DARK_MUTED)

        legend = ax.get_legend()
        if legend is not None:
            legend.get_frame().set_facecolor(DARK_BG)
            legend.get_frame().set_edgecolor(DARK_GRID)
            for txt in legend.get_texts():
                txt.set_color(DARK_LABEL)

    # Patches and collections are swept figure-wide rather than per-axes:
    # artists added via add_artist (and the wedges/collections the circular
    # genome diagram builds) never appear in ax.patches, so a per-axes loop
    # silently misses them. _remap only rewrites known chrome colors, so
    # data-colored features pass through untouched.
    for patch in fig.findobj(mpatches.Patch):
        new_edge = _remap(patch.get_edgecolor())
        if new_edge:
            patch.set_edgecolor(new_edge)
        new_face = _remap(patch.get_facecolor())
        if new_face:
            patch.set_facecolor(new_face)

    for coll in fig.findobj(mcoll.Collection):
        try:
            edges = coll.get_edgecolor()
            if len(edges):
                remapped = [_remap(e) or mcolors.to_hex(e, keep_alpha=True)
                            for e in edges]
                coll.set_edgecolor(remapped)
        except (AttributeError, ValueError, TypeError):
            pass

    # Final sweep for any text the per-axes pass cannot reach — most notably
    # fig.suptitle(), which matplotlib stores as fig._suptitle rather than in
    # fig.texts. findobj walks the artist tree, so this catches titles,
    # annotations and legend labels wherever they live. Tick labels were
    # already set to DARK_MUTED above, and DARK_MUTED is not a key in
    # _DARK_CHROME, so they are left alone here.
    for txt in fig.findobj(mtext.Text):
        new = _remap(txt.get_color())
        if new:
            txt.set_color(new)

    return fig


def _theme_facecolor(theme):
    return DARK_BG if theme == 'dark' else 'white'


def resolve_color(c, default='#2563eb'):
    if not c:
        return default
    c = str(c).strip()
    if c.startswith('#'):
        return c
    return PALETTE.get(c.lower(), default)


def figure_to_svg_data_url(fig, pad_inches=0.2, theme=None, transparent=False):
    """Serialize a matplotlib Figure to an svg+xml base64 data URL.

    The figure is recolored for dark mode when the current request asks for
    it (see resolve_theme). Pass transparent=True for figures that should sit
    directly on the card background rather than carrying their own ground —
    the chrome is still themed, only the backdrop is dropped. Callers should
    not call plt.close() themselves — this helper does.
    """
    theme = resolve_theme(theme)
    if theme == 'dark':
        apply_dark_theme(fig)
    buf = StringIO()
    if transparent:
        fig.savefig(buf, format='svg', bbox_inches='tight',
                    transparent=True, pad_inches=pad_inches)
    else:
        fig.savefig(buf, format='svg', bbox_inches='tight',
                    facecolor=_theme_facecolor(theme), pad_inches=pad_inches)
    plt.close(fig)
    svg = buf.getvalue().encode('utf-8')
    b64 = base64.b64encode(svg).decode()
    return f'data:image/svg+xml;base64,{b64}'


def figure_to_png_data_url(fig, dpi=150, pad_inches=0.2, theme=None,
                           transparent=False):
    """PNG fallback for the rare cases where SVG is not appropriate
    (e.g. when the backend draws pixels directly)."""
    theme = resolve_theme(theme)
    if theme == 'dark':
        apply_dark_theme(fig)
    buf = BytesIO()
    if transparent:
        fig.savefig(buf, format='png', bbox_inches='tight',
                    transparent=True, dpi=dpi, pad_inches=pad_inches)
    else:
        fig.savefig(buf, format='png', bbox_inches='tight',
                    facecolor=_theme_facecolor(theme), dpi=dpi,
                    pad_inches=pad_inches)
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
