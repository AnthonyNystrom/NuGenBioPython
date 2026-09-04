"""GC skew and dot plots.

Both are standard sequence diagnostics that only mean anything as pictures.
GC skew was computed and then shown to the user as "First 10 values:
0.1234, 0.2345, ..."; the dot plot did not exist.
"""
import random

import pytest

from utils.sequence_plots import (
    gc_skew_series, dot_plot_points, render_gc_skew, render_dot_plot,
    reverse_complement, axis_bp_labels, MAX_DOTPLOT_LEN,
)


def _rand_seq(n, seed=0):
    rng = random.Random(seed)
    return ''.join(rng.choice('ACGT') for _ in range(n))


def _synthetic_replicon(n=60000, ori_frac=0.20, ter_frac=0.70, seed=7):
    """Leading strand G-rich, lagging strand C-rich — the real biology."""
    rng = random.Random(seed)
    ori, ter = int(ori_frac * n), int(ter_frac * n)
    out = []
    for i in range(n):
        leading = ori <= i < ter
        r = rng.random()
        if leading:
            out.append('G' if r < 0.34 else 'C' if r < 0.50 else rng.choice('AT'))
        else:
            out.append('C' if r < 0.34 else 'G' if r < 0.50 else rng.choice('AT'))
    return ''.join(out), ori, ter


# --------------------------------------------------------------------------
# GC skew
# --------------------------------------------------------------------------

def test_gc_skew_sign_follows_strand_composition():
    assert gc_skew_series("G" * 100, window=50)[1][0] == pytest.approx(1.0)
    assert gc_skew_series("C" * 100, window=50)[1][0] == pytest.approx(-1.0)


def test_gc_skew_is_zero_without_g_or_c():
    _, skews, _ = gc_skew_series("ATATATAT" * 10, window=20)
    assert all(s == 0.0 for s in skews)


def test_gc_skew_handles_empty_sequence():
    assert gc_skew_series("", window=10) == ([], [], [])


def test_cumulative_skew_locates_origin_and_terminus():
    """The property that makes this plot worth drawing at all."""
    seq, ori, ter = _synthetic_replicon()
    positions, _skews, cumulative = gc_skew_series(seq, window=2000)
    detected_ori = positions[cumulative.index(min(cumulative))]
    detected_ter = positions[cumulative.index(max(cumulative))]
    tolerance = 0.05 * len(seq)
    assert abs(detected_ori - ori) < tolerance, \
        f"origin off by {abs(detected_ori - ori)} bp"
    assert abs(detected_ter - ter) < tolerance, \
        f"terminus off by {abs(detected_ter - ter)} bp"


def test_gc_skew_renders_svg():
    seq, _, _ = _synthetic_replicon(n=8000)
    url = render_gc_skew(seq, window=500)
    assert url.lstrip().startswith("<svg")


def test_gc_skew_render_returns_none_for_empty():
    assert render_gc_skew("", window=10) is None


# --------------------------------------------------------------------------
# Axis labels
# --------------------------------------------------------------------------

def test_axis_labels_use_one_consistent_unit():
    """fmt_bp switches unit per value, giving axes like 5,000 / 10 kb."""
    labels = axis_bp_labels([0, 5000, 10000, 15000, 20000])
    assert labels == ['0', '5 kb', '10 kb', '15 kb', '20 kb']


def test_axis_labels_scale_to_megabases():
    labels = axis_bp_labels([0, 500_000, 1_000_000, 1_500_000])
    assert labels == ['0', '0.5 Mb', '1 Mb', '1.5 Mb']


def test_axis_labels_stay_in_bp_when_small():
    assert axis_bp_labels([0, 100, 200]) == ['0', '100', '200']


def test_axis_labels_handle_empty():
    assert axis_bp_labels([]) == []


# --------------------------------------------------------------------------
# Dot plot
# --------------------------------------------------------------------------

def test_reverse_complement():
    assert reverse_complement("ATGC") == "GCAT"
    assert reverse_complement("AAAA") == "TTTT"


def test_self_comparison_fills_the_main_diagonal():
    seq = _rand_seq(300, seed=1)
    (fx, fy), _rev, _t = dot_plot_points(seq, seq, word_size=8)
    on_diagonal = sum(1 for x, y in zip(fx, fy) if x == y)
    assert on_diagonal == len(seq) - 8 + 1


def test_duplication_produces_off_diagonal_matches():
    block = _rand_seq(200, seed=2)
    seq = _rand_seq(100, seed=3) + block + _rand_seq(100, seed=4) + block
    (fx, fy), _rev, _t = dot_plot_points(seq, seq, word_size=10)
    off = sum(1 for x, y in zip(fx, fy) if abs(x - y) > 5)
    assert off > 0, "a duplicated block must show off the main diagonal"


def test_inversion_shows_as_reverse_complement_matches():
    """The regression: forward-only indexing made inversions invisible."""
    inv = _rand_seq(300, seed=5)
    seq1 = _rand_seq(200, seed=6) + inv
    seq2 = _rand_seq(200, seed=6) + reverse_complement(inv)
    _fwd, (rx, _ry), _t = dot_plot_points(seq1, seq2, word_size=10)
    assert len(rx) > 100, "the inversion must appear as reverse matches"


def test_reverse_matching_can_be_disabled():
    inv = _rand_seq(200, seed=8)
    _fwd, (rx, _ry), _t = dot_plot_points(
        inv, reverse_complement(inv), word_size=10, include_reverse=False)
    assert rx == []


def test_unrelated_sequences_have_almost_no_matches():
    a, b = _rand_seq(500, seed=11), _rand_seq(500, seed=12)
    (fx, _fy), _rev, _t = dot_plot_points(a, b, word_size=12)
    assert len(fx) < 5


def test_sequences_shorter_than_the_word_are_handled():
    assert dot_plot_points("AC", "AC", word_size=8) == (([], []), ([], []), False)


def test_dot_plot_renders_svg():
    seq = _rand_seq(400, seed=13)
    assert render_dot_plot(seq, seq, word_size=8).lstrip().startswith("<svg")


def test_dot_plot_renders_with_no_matches():
    """Must not blow up when nothing matches."""
    url = render_dot_plot(_rand_seq(200, seed=14), _rand_seq(200, seed=15),
                          word_size=20)
    assert url.lstrip().startswith("<svg")


# --------------------------------------------------------------------------
# Endpoint
# --------------------------------------------------------------------------

def test_dotplot_endpoint_reports_both_orientations(client):
    inv = _rand_seq(300, seed=16)
    resp = client.post("/api/sequence/dotplot", json={
        "sequence1": _rand_seq(200, seed=17) + inv,
        "sequence2": _rand_seq(200, seed=17) + reverse_complement(inv),
        "word_size": 10,
    })
    data = resp.get_json()
    assert data["success"] is True
    assert data["reverse_matches"] > 100
    assert data["plot"].lstrip().startswith("<svg")


def test_dotplot_endpoint_rejects_oversized_input(client):
    resp = client.post("/api/sequence/dotplot",
                       json={"sequence1": "A" * (MAX_DOTPLOT_LEN + 1)})
    data = resp.get_json()
    assert data["success"] is False
    assert "limited to" in data["error"]


def test_dotplot_endpoint_requires_a_sequence(client):
    data = client.post("/api/sequence/dotplot", json={}).get_json()
    assert data["success"] is False


def test_dotplot_defaults_to_self_comparison(client):
    data = client.post("/api/sequence/dotplot",
                       json={"sequence1": _rand_seq(200, seed=18)}).get_json()
    assert data["success"] is True
    assert data["self_comparison"] is True


def test_gc_analysis_returns_a_plot(client):
    seq, _, _ = _synthetic_replicon(n=6000)
    data = client.post("/api/sequence/gc_analysis",
                       json={"sequence": seq, "window": 500}).get_json()
    assert data["success"] is True
    assert data["skew_plot"].lstrip().startswith("<svg")
