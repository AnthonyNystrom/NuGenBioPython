"""Server-rendered figures must follow the app's dark theme.

The app has a real dark mode (html[data-theme="dark"]), but figures are
delivered as opaque `<img src="data:...">` URLs that CSS cannot reach, so
the theme has to travel to the renderer. utils.plot_helpers.resolve_theme
reads it from the request; these tests pin both the resolution rules and
the resulting colors.
"""
import base64
import re

import pytest

from utils.plot_helpers import (
    apply_dark_theme, figure_to_svg_data_url, resolve_theme,
    TITLE_COLOR, LABEL_COLOR, MUTED_COLOR, AXIS_COLOR, GRID_COLOR,
    DARK_BG, DARK_LABEL, PALETTE,
)

# Light chrome colors that must never survive into a dark render.
LIGHT_CHROME = {TITLE_COLOR.lower(), LABEL_COLOR.lower(), MUTED_COLOR.lower(),
                AXIS_COLOR.lower(), GRID_COLOR.lower(), '#ffffff'}


def _svg_colors(data_url):
    svg = (data_url if data_url.lstrip().startswith("<svg")
           else base64.b64decode(data_url.split(',', 1)[1]).decode())
    return set(c.lower() for c in re.findall(r'#[0-9a-fA-F]{6}', svg))


def _demo_figure():
    import matplotlib.pyplot as plt
    from utils.plot_helpers import style_axes, set_title
    fig, ax = plt.subplots(figsize=(4, 3))
    ax.plot([0, 1, 2], [1, 3, 2], color=PALETTE['blue'])
    ax.bar([0, 1, 2], [1, 2, 1], color=PALETTE['red'], edgecolor=TITLE_COLOR)
    ax.set_xlabel('position')
    style_axes(ax, grid=True)
    set_title(ax, 'Demo')
    fig.suptitle('Suptitle')
    return fig


def test_default_theme_is_light_outside_a_request():
    """Offline callers (tests/visual scripts) must keep the light rendering."""
    assert resolve_theme() == 'light'


def test_explicit_theme_wins():
    assert resolve_theme('dark') == 'dark'
    assert resolve_theme('light') == 'light'


def test_light_render_keeps_white_ground():
    colors = _svg_colors(figure_to_svg_data_url(_demo_figure(), theme='light'))
    assert '#ffffff' in colors


def test_dark_render_drops_all_light_chrome():
    colors = _svg_colors(figure_to_svg_data_url(_demo_figure(), theme='dark'))
    leaked = LIGHT_CHROME & colors
    assert not leaked, f"light chrome leaked into dark render: {sorted(leaked)}"
    assert DARK_BG.lower() in colors


def test_dark_render_preserves_data_colors():
    """Series colors must survive so charts still match their legends."""
    colors = _svg_colors(figure_to_svg_data_url(_demo_figure(), theme='dark'))
    assert PALETTE['blue'].lower() in colors
    assert PALETTE['red'].lower() in colors


def test_suptitle_is_recolored():
    """fig.suptitle lives outside fig.texts; it must still be swept."""
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    sup = fig.suptitle('Title', color=TITLE_COLOR)
    apply_dark_theme(fig)
    assert sup.get_color().lower() != TITLE_COLOR.lower()
    plt.close(fig)


@pytest.mark.parametrize("header,expect_white", [("dark", False), ("light", True)])
def test_theme_header_reaches_the_renderer(client, header, expect_white):
    """End-to-end: X-Theme on an API call changes the returned figure."""
    resp = client.post(
        "/api/phylo/parse",
        data={"tree_string": "((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);",
              "format": "newick"},
        headers={"X-Theme": header},
    )
    assert resp.status_code == 200
    data = resp.get_json()
    assert data.get("success") is True, data.get("error")
    colors = _svg_colors(data["tree_image"])
    assert ('#ffffff' in colors) is expect_white
    if not expect_white:
        assert DARK_BG.lower() in colors


def test_theme_query_arg_also_works(client):
    resp = client.post(
        "/api/phylo/parse?theme=dark",
        data={"tree_string": "((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);",
              "format": "newick"},
    )
    data = resp.get_json()
    assert data.get("success") is True
    assert DARK_BG.lower() in _svg_colors(data["tree_image"])


# --------------------------------------------------------------------------
# Transparency
# --------------------------------------------------------------------------

def test_transparent_fills_stay_invisible_in_dark_mode():
    """Regression: matplotlib.colors.to_hex discards alpha, so a fully
    transparent artist reported as (0,0,0,0) looked like pure black and was
    remapped to near-white — turning `facecolor='none'` shapes into solid
    light blobs over the whole plot."""
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse

    fig, ax = plt.subplots()
    patch = Ellipse((0, 0), 1, 1, facecolor='none', edgecolor='red')
    ax.add_patch(patch)
    apply_dark_theme(fig)
    assert mcolors.to_rgba(patch.get_facecolor())[3] == 0
    plt.close(fig)


def test_opaque_chrome_is_still_remapped():
    """The transparency guard must not stop real chrome from being themed."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    fig, ax = plt.subplots()
    patch = Rectangle((0, 0), 1, 1, facecolor=TITLE_COLOR)
    ax.add_patch(patch)
    apply_dark_theme(fig)
    from matplotlib.colors import to_hex
    assert to_hex(patch.get_facecolor()).lower() != TITLE_COLOR.lower()
    plt.close(fig)


def test_ramachandran_render_keeps_transparent_regions_transparent():
    """End-to-end on the plot that exposed the transparency bug.

    The glycine-only contours are drawn with facecolor='none'. If the guard
    in _remap regresses they become filled shapes. Figures are now delivered
    as inline SVG with presentation attributes, so the check is on fill="..."
    rather than a CSS declaration.
    """
    import re
    from utils.structure_plots import classify_all, render_ramachandran
    from utils.plot_helpers import DARK_TITLE

    _counts, annotated = classify_all([
        {'phi': -63, 'psi': -43, 'resname': 'GLY', 'residue': 'GLY1'},
        {'phi': -139, 'psi': 135, 'resname': 'ALA', 'residue': 'ALA2'},
    ])
    svg = render_ramachandran(annotated)
    assert svg.lstrip().startswith("<svg")
    fills = re.findall(r'fill="(#[0-9a-fA-F]{6})"', svg)
    assert fills, "expected presentation-attribute fills"
    assert DARK_TITLE.lower() not in [f.lower() for f in fills]
