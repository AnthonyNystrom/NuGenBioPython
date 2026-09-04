"""Font weights must be named, not numeric.

matplotlib's SVG text writer resolves weights through
font_manager.weight_dict, which is keyed by NAME ('semibold', 'bold', ...).
A numeric string like '600' is tolerated by newer matplotlib but raises
KeyError: 600 on 3.9 — and only on the svg.fonttype='none' path, i.e. exactly
the inline figures this app now serves. CI on Python 3.9 failed on 12 tests
for this reason while 3.12 passed.
"""
import pathlib
import re

import pytest
from matplotlib import font_manager


NUMERIC_WEIGHT = re.compile(r"""fontweight\s*=\s*['"](\d+)['"]""")


def _python_sources(repo_root):
    for folder in ("routes", "utils"):
        yield from pathlib.Path(repo_root, folder).rglob("*.py")


def test_weight_dict_is_keyed_by_name_not_number():
    """The premise: this is why numeric weights break."""
    assert "semibold" in font_manager.weight_dict
    assert "bold" in font_manager.weight_dict
    assert 600 not in font_manager.weight_dict
    assert "600" not in font_manager.weight_dict


def test_named_weights_resolve_to_the_intended_numbers():
    assert font_manager.weight_dict["semibold"] == 600
    assert font_manager.weight_dict["bold"] == 700


def test_no_numeric_font_weights_in_source(repo_root):
    offenders = []
    for path in _python_sources(repo_root):
        for i, line in enumerate(path.read_text().splitlines(), 1):
            if line.strip().startswith("#"):
                continue
            match = NUMERIC_WEIGHT.search(line)
            if match:
                offenders.append(f"{path.name}:{i} fontweight='{match.group(1)}'")
    assert not offenders, (
        "matplotlib resolves weights by name; a numeric string raises "
        f"KeyError on older matplotlib when writing inline SVG: {offenders}")


@pytest.mark.parametrize("weight", ["semibold", "bold", "normal"])
def test_inline_svg_renders_with_named_weights(weight):
    """End-to-end on the path that actually broke: svg.fonttype='none'."""
    from utils.plot_helpers import subplots, svg_markup
    fig, ax = subplots(figsize=(3, 2))
    ax.set_title("Title", fontweight=weight)
    svg = svg_markup(fig)
    assert svg.lstrip().startswith("<svg")
    assert "Title" in svg


def test_every_renderer_survives_inline_serialization():
    """A smoke pass over the real renderers, which is what CI exercised."""
    import utils.phylo_helpers as ph
    from utils.sequence_plots import render_gc_skew, render_dot_plot
    import routes.genomediagram_routes as gd

    tree = ph.parse_tree_from_string("((A:0.1,B:0.2):0.3,C:0.4);", "newick")
    tracks = [{"name": "G", "features": [
        {"start": 100, "end": 1200, "strand": 1, "name": "a", "color": "blue"}]}]

    for produced in (
        ph.tree_to_inline_svg(tree),
        render_gc_skew("ATGCGTACGT" * 200, window=200),
        render_dot_plot("ACGT" * 100, "ACGT" * 100, word_size=6),
        gd.render_linear_feature_diagram(1_000_000, tracks, "x"),
        gd.render_circular_feature_diagram(1_000_000, tracks, "x"),
    ):
        assert produced.lstrip().startswith("<svg")
