"""Sequence-logo rendering: theming and text contrast.

The logo used to serialize itself with a bespoke fig.savefig() and hardcoded
chrome colors, so it bypassed the app's design system entirely. Its tick
labels sat at #94a3b8 on white — 2.56:1, well under WCAG AA's 4.5:1 for small
text — in the *default* light theme.
"""
import base64
import re

import pytest
from Bio import motifs
from Bio.Seq import Seq

import utils.plot_helpers as ph
from utils.plot_helpers import DARK_BG, MUTED_COLOR

WCAG_AA_SMALL_TEXT = 4.5
LIGHT_CARD = "#ffffff"


def _relative_luminance(hex_color):
    h = hex_color.lstrip("#")
    channels = [int(h[i:i + 2], 16) / 255 for i in (0, 2, 4)]
    channels = [c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4
                for c in channels]
    return 0.2126 * channels[0] + 0.7152 * channels[1] + 0.0722 * channels[2]


def contrast_ratio(fg, bg):
    a, b = _relative_luminance(fg), _relative_luminance(bg)
    hi, lo = max(a, b), min(a, b)
    return (hi + 0.05) / (lo + 0.05)


@pytest.fixture
def motif():
    return motifs.create([Seq(s) for s in
                          ["TATAAA", "TATAAT", "TATATA",
                           "TACAAA", "TATAAA", "TAAAAT"]])


def _render(motif, theme, monkeypatch):
    import matplotlib.pyplot as plt
    from routes.motif_routes import generate_sequence_logo
    monkeypatch.setattr(ph, "resolve_theme", lambda explicit=None: theme)
    url = generate_sequence_logo(motif)
    plt.close("all")
    assert url, "logo generation returned None"
    svg = base64.b64decode(url.split(",", 1)[1]).decode()
    return svg, set(c.lower() for c in re.findall(r"#[0-9a-fA-F]{6}", svg))


def test_contrast_helper_matches_known_values():
    """Sanity-check the helper against hand-computed ratios."""
    assert contrast_ratio("#ffffff", "#000000") == pytest.approx(21.0, abs=0.01)
    assert contrast_ratio("#94a3b8", "#ffffff") == pytest.approx(2.56, abs=0.01)


def test_old_low_contrast_tick_color_is_gone(motif, monkeypatch):
    """#94a3b8 on white was 2.56:1 — it must not come back."""
    _, colors = _render(motif, "light", monkeypatch)
    assert "#94a3b8" not in colors


def test_light_chrome_meets_wcag_aa(motif, monkeypatch):
    _, colors = _render(motif, "light", monkeypatch)
    assert MUTED_COLOR.lower() in colors
    assert contrast_ratio(MUTED_COLOR, LIGHT_CARD) >= WCAG_AA_SMALL_TEXT


def test_dark_chrome_meets_wcag_aa(motif, monkeypatch):
    _, colors = _render(motif, "dark", monkeypatch)
    themed = [c for c in colors
              if contrast_ratio(c, DARK_BG) >= WCAG_AA_SMALL_TEXT]
    assert themed, "no chrome color in the dark render reaches AA contrast"
    assert MUTED_COLOR.lower() not in colors, "light chrome leaked into dark"


def test_logo_stays_transparent_in_both_themes(motif, monkeypatch):
    """The logo sits on the card; it must not paint its own ground."""
    light, _ = _render(motif, "light", monkeypatch)
    dark, _ = _render(motif, "dark", monkeypatch)
    assert "fill: #ffffff" not in light
    assert DARK_BG.lower() not in dark.lower()


def test_base_glyph_colors_are_preserved(motif, monkeypatch):
    """A/C/G/T keep their identity colors so the logo stays readable."""
    _, colors = _render(motif, "dark", monkeypatch)
    for base_color in ("#059669", "#2563eb", "#ca8a04", "#dc2626"):
        assert base_color in colors, f"{base_color} lost in dark render"
