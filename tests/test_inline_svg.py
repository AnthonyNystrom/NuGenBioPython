"""Inline SVG figure delivery.

Figures used to ship as `<img src="data:image/svg+xml;base64,...">` — opaque
to the page: no hover, no zoom, no selectable text, and CSS could not reach
inside them.

Injecting the SVG instead makes all of that possible, but matplotlib emits
~30 `style="..."` attributes per figure plus a `<style>` block, and this app's
CSP has no 'unsafe-inline' in style-src. Inlined unchanged, every one is
refused and the figure renders unstyled (measured in Chromium: 30 violations,
computed stroke "none"). svg_markup() rewrites them as SVG presentation
attributes, which are XML attributes rather than CSS and so sit outside CSP
entirely — and which CSS outranks, so a stylesheet can still restyle a figure.
"""
import re

import pytest

from utils.plot_helpers import (
    subplots, svg_markup, figure_to_svg_data_url, PALETTE, TITLE_COLOR,
)

STYLE_ATTR = re.compile(r'(?<![\w-])style\s*=\s*["\']')


def _figure():
    fig, ax = subplots(figsize=(4, 3))
    ax.plot([0, 1, 2], [1, 3, 2], color=PALETTE['blue'], linewidth=2)
    ax.set_title("Demo title")
    ax.set_xlabel("x label")
    return fig


# --------------------------------------------------------------------------
# The converter
# --------------------------------------------------------------------------

def test_output_is_injectable_markup():
    svg = svg_markup(_figure())
    assert svg.lstrip().startswith("<svg")
    assert "<?xml" not in svg          # invalid inside an HTML document
    assert "<!DOCTYPE" not in svg


def test_no_style_attributes_survive():
    """The whole reason this converter exists."""
    svg = svg_markup(_figure())
    assert not STYLE_ATTR.search(svg), "a style attribute would be refused by CSP"


def test_no_style_element_survives():
    """matplotlib emits `*{stroke-linejoin: round; stroke-linecap: butt}`."""
    svg = svg_markup(_figure())
    assert "<style" not in svg


def test_colors_become_presentation_attributes():
    svg = svg_markup(_figure())
    assert re.search(r'fill="#[0-9a-fA-F]{6}"', svg)
    assert re.search(r'stroke="#[0-9a-fA-F]{6}"', svg)
    assert PALETTE['blue'] in svg, "the series colour must survive"


def test_labels_are_real_text_not_glyph_outlines():
    """svg.fonttype='none' is what makes labels selectable and searchable."""
    svg = svg_markup(_figure())
    assert svg.count("<text") > 0
    assert "Demo title" in svg
    assert "x label" in svg


def test_unknown_css_properties_are_dropped_not_emitted():
    """Only real SVG presentation attributes may be written out."""
    import matplotlib
    fig, ax = subplots(figsize=(2, 2))
    ax.plot([0, 1], [0, 1])
    svg = svg_markup(fig)
    # a CSS-only property must never appear as an attribute
    assert not re.search(r'\s(?:box-shadow|display|position)=', svg)


def test_figure_is_marked_and_labelled():
    svg = svg_markup(_figure(), title="Phylogenetic tree")
    assert 'class="figure-svg"' in svg
    assert 'role="img"' in svg
    assert 'aria-label="Phylogenetic tree"' in svg


def test_quotes_in_a_title_cannot_break_the_attribute():
    svg = svg_markup(_figure(), title='a "quoted" title')
    assert 'aria-label="a &quot;quoted&quot; title"' in svg


def test_the_export_path_is_unchanged():
    """Downloads stay self-contained: glyph outlines, no font dependency."""
    url = figure_to_svg_data_url(_figure())
    assert url.startswith("data:image/svg+xml;base64,")
    import base64
    raw = base64.b64decode(url.split(",", 1)[1]).decode()
    assert "<?xml" in raw
    assert raw.count("<text") == 0, "exports keep glyph paths for portability"


def test_transparent_and_themed_variants_still_convert():
    for kwargs in ({'transparent': True}, {'theme': 'dark'}):
        svg = svg_markup(_figure(), **kwargs)
        assert svg.lstrip().startswith("<svg")
        assert not STYLE_ATTR.search(svg)


# --------------------------------------------------------------------------
# Endpoints deliver inline markup
# --------------------------------------------------------------------------

@pytest.mark.parametrize("path,payload,key", [
    ("/api/sequence/dotplot", {"sequence1": "ATGCGTACGT" * 20}, "plot"),
    ("/api/sequence/gc_analysis",
     {"sequence": "ATGCGTACGT" * 200, "window": 100}, "skew_plot"),
])
def test_json_endpoints_return_inline_svg(client, path, payload, key):
    data = client.post(path, json=payload).get_json()
    assert data["success"] is True
    figure = data[key]
    assert figure.lstrip().startswith("<svg")
    assert not STYLE_ATTR.search(figure)


def test_phylo_returns_both_inline_and_export_forms(client):
    data = client.post("/api/phylo/parse",
                       data={"tree_string": "((A:0.1,B:0.2):0.3,C:0.4);",
                             "format": "newick"}).get_json()
    assert data["success"] is True
    assert data["tree_svg"].lstrip().startswith("<svg")
    assert data["tree_image"].startswith("data:image/svg+xml;base64,")


def test_no_renderer_still_returns_a_data_url_for_display(repo_root):
    """Guard: a figure meant for the page must be inline markup."""
    import pathlib
    offenders = []
    for folder in ("routes", "utils"):
        for path in pathlib.Path(repo_root, folder).rglob("*.py"):
            if path.name in ("plot_helpers.py", "phylo_helpers.py"):
                continue     # both intentionally offer the export form too
            for i, line in enumerate(path.read_text().splitlines(), 1):
                if line.strip().startswith("#"):
                    continue
                if "figure_to_svg_data_url(" in line:
                    offenders.append(f"{path.name}:{i}")
    assert not offenders, (
        "figures for display should use svg_markup(); "
        f"figure_to_svg_data_url is for exports: {offenders}")
