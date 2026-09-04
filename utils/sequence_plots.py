"""Sequence-level diagnostic plots: GC skew and dot plots.

Both of these were either computed and then thrown away as a list of numbers
(GC skew was rendered as "First 10 values: 0.1234, 0.2345, ...") or absent
entirely. They are the two most standard sequence plots in the field and both
only mean anything as a picture.
"""
import logging

import numpy as np

log = logging.getLogger(__name__)

# Guard rails. A dot plot is O(n*m) in the worst case and a browser cannot
# usefully display more points than this anyway.
MAX_DOTPLOT_LEN = 6000
MAX_DOTPLOT_POINTS = 200_000



def axis_bp_labels(ticks):
    """Format a whole tick series in one consistent unit.

    plot_helpers.fmt_bp switches unit per value, which is right for a single
    figure caption but produces axes reading "5,000 / 10 kb / 15 kb". Pick the
    unit once from the largest tick and apply it to every label.
    """
    if not ticks:
        return []
    top = max(abs(t) for t in ticks)
    if top >= 1_000_000:
        div, unit, places = 1_000_000.0, 'Mb', 2
    elif top >= 10_000:
        div, unit, places = 1000.0, 'kb', 1
    else:
        return [f'{int(t):,}' for t in ticks]
    out = []
    for t in ticks:
        v = t / div
        text = f'{v:.{places}f}'.rstrip('0').rstrip('.')
        out.append(f'{text} {unit}' if t else '0')
    return out

def gc_skew_series(sequence, window=1000, step=None):
    """Return (positions, windowed_skew, cumulative_skew).

    GC skew is (G - C) / (G + C) over a sliding window. Its *cumulative* sum
    is the biologically useful curve: in most bacterial genomes the global
    minimum marks the replication origin and the maximum marks the terminus,
    because leading- and lagging-strand mutational bias flips there.
    """
    seq = str(sequence).upper()
    n = len(seq)
    if n == 0:
        return [], [], []
    window = max(1, min(int(window), n))
    step = step or max(1, window // 2)

    positions, skews = [], []
    for start in range(0, n, step):
        chunk = seq[start:start + window]
        if not chunk:
            break
        g = chunk.count('G')
        c = chunk.count('C')
        total = g + c
        skews.append((g - c) / total if total else 0.0)
        positions.append(start + len(chunk) // 2)

    cumulative = np.cumsum(skews).tolist() if skews else []
    return positions, skews, cumulative


def render_gc_skew(sequence, window=1000, title='GC skew'):
    """Render windowed + cumulative GC skew as an SVG data URL."""

    from utils.plot_helpers import (
        svg_markup,
        subplots as oo_subplots,
        PALETTE, MUTED_COLOR, GRID_COLOR,
        set_title, style_axes, fmt_bp, nice_ticks,
    )

    positions, skews, cumulative = gc_skew_series(sequence, window)
    if not positions:
        return None

    fig, (ax_top, ax_bot) = oo_subplots(
        2, 1, figsize=(11, 5.6), sharex=True,
        gridspec_kw={'height_ratios': [1, 1.3], 'hspace': 0.18})

    # --- windowed skew ---------------------------------------------------
    ax_top.axhline(0, color=GRID_COLOR, linewidth=0.8, zorder=1)
    ax_top.fill_between(positions, skews, 0,
                        where=[s >= 0 for s in skews],
                        color=PALETTE['blue'], alpha=0.55, linewidth=0,
                        interpolate=True, label='G-rich (leading)')
    ax_top.fill_between(positions, skews, 0,
                        where=[s < 0 for s in skews],
                        color=PALETTE['orange'], alpha=0.55, linewidth=0,
                        interpolate=True, label='C-rich (lagging)')
    ax_top.set_ylabel('GC skew', fontsize=9, color=MUTED_COLOR)
    style_axes(ax_top)
    ax_top.legend(loc='upper right', fontsize=7, ncol=2, framealpha=0.85)

    # --- cumulative skew -------------------------------------------------
    ax_bot.plot(positions, cumulative, color=PALETTE['darkblue'],
                linewidth=1.4, zorder=3)
    ax_bot.axhline(0, color=GRID_COLOR, linewidth=0.8, zorder=1)

    imin = int(np.argmin(cumulative))
    imax = int(np.argmax(cumulative))
    for idx, color, label in ((imin, PALETTE['green'], 'Origin (min)'),
                              (imax, PALETTE['red'], 'Terminus (max)')):
        ax_bot.axvline(positions[idx], color=color, linewidth=1.0,
                       linestyle='--', alpha=0.9, zorder=2)
        ax_bot.scatter([positions[idx]], [cumulative[idx]], s=42, color=color,
                       zorder=4, edgecolors='white', linewidths=0.8,
                       label=f'{label} ~{fmt_bp(positions[idx])}')

    ax_bot.set_ylabel('Cumulative skew', fontsize=9, color=MUTED_COLOR)
    ax_bot.set_xlabel('Position', fontsize=9, color=MUTED_COLOR)
    style_axes(ax_bot)
    ax_bot.legend(loc='best', fontsize=7.5, framealpha=0.85)

    ticks = nice_ticks(len(str(sequence)))
    ax_bot.set_xticks(ticks)
    ax_bot.set_xticklabels(axis_bp_labels(ticks))
    ax_bot.set_xlim(0, len(str(sequence)))

    set_title(fig, f'{title} (window {fmt_bp(window)})')
    return svg_markup(fig, title='GC skew')


_COMPLEMENT = str.maketrans('ACGTacgt', 'TGCAtgca')
_UNAMBIGUOUS_DNA = set('ACGT')


def reverse_complement(seq):
    return str(seq).translate(_COMPLEMENT)[::-1]


def looks_like_nucleotide(seq, threshold=0.85):
    """True when the sequence is predominantly A/C/G/T/U/N.

    Reverse-complement matching is only meaningful for nucleotides; run on a
    protein it silently produces garbage, because the complement table leaves
    every non-ACGT letter untouched and just reverses the string.
    """
    text = str(seq).upper()
    if not text:
        return False
    hits = sum(1 for c in text if c in 'ACGTUN')
    return hits / len(text) >= threshold


def dot_plot_points(seq1, seq2, word_size=8, include_reverse=True):
    """Return (forward, reverse, truncated) word matches between two sequences.

    `forward` and `reverse` are each (xs, ys). Reverse-complement matches are
    tracked separately because they are the whole point of looking at a dot
    plot: a run perpendicular to the main diagonal is an inversion, and it is
    invisible if you only index the forward strand.

    Words containing anything but A/C/G/T are skipped. An assembly gap is a
    long run of N, and a run of N matches every other run of N — in testing, a
    single 40 bp gap produced 93% of all the "matches" in a self-comparison,
    burying the real signal under a solid block.

    Indexes seq1's k-mers once and walks seq2, so this is linear in the
    sequences rather than an O(n*m) all-pairs comparison.
    """
    a = str(seq1).upper()
    b = str(seq2).upper()
    word_size = max(2, int(word_size))
    if len(a) < word_size or len(b) < word_size:
        return ([], []), ([], []), False

    # Only complement when both sides actually look like nucleotides.
    if include_reverse and not (looks_like_nucleotide(a)
                                and looks_like_nucleotide(b)):
        include_reverse = False

    nucleotide = looks_like_nucleotide(a) and looks_like_nucleotide(b)

    def informative(word):
        """Reject ambiguity runs, which match each other and nothing real."""
        if not nucleotide:
            return 'X' not in word
        return set(word) <= _UNAMBIGUOUS_DNA

    index = {}
    for i in range(len(a) - word_size + 1):
        word = a[i:i + word_size]
        if informative(word):
            index.setdefault(word, []).append(i)

    fx, fy, rx, ry = [], [], [], []
    total = 0
    for j in range(len(b) - word_size + 1):
        word = b[j:j + word_size]
        if not informative(word):
            continue
        pairs = [(index.get(word), fx, fy)]
        if include_reverse:
            pairs.append((index.get(reverse_complement(word)), rx, ry))
        for hits, xs, ys in pairs:
            if not hits:
                continue
            for i in hits:
                xs.append(i)
                ys.append(j)
                total += 1
                if total >= MAX_DOTPLOT_POINTS:
                    return (fx, fy), (rx, ry), True
    return (fx, fy), (rx, ry), False


def render_dot_plot(seq1, seq2, word_size=8, label1='Sequence 1',
                    label2='Sequence 2', title='Dot plot',
                    include_reverse=True):
    """Render a dot plot as an SVG data URL.

    Diagonal runs are conserved segments; a diagonal parallel to the main one
    is a repeat or duplication, and a run perpendicular to it is an inversion.
    None of that is visible in a pairwise alignment score.
    """

    from utils.plot_helpers import (
        svg_markup,
        subplots as oo_subplots,
        PALETTE, MUTED_COLOR, GRID_COLOR,
        set_title, style_axes,
    )

    (fx, fy), (rx, ry), truncated = dot_plot_points(
        seq1, seq2, word_size, include_reverse)

    fig, ax = oo_subplots(figsize=(7.2, 7.0))
    if fx:
        ax.scatter(fx, fy, s=1.2, c=PALETTE['darkblue'], alpha=0.55,
                   linewidths=0, marker='s',
                   label=f'Forward ({len(fx):,})')
    if rx:
        ax.scatter(rx, ry, s=1.2, c=PALETTE['red'], alpha=0.55,
                   linewidths=0, marker='s',
                   label=f'Reverse complement ({len(rx):,})')
    if not fx and not rx:
        ax.text(0.5, 0.5, f'No exact matches of {word_size} bp',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=10, color=MUTED_COLOR)

    ax.set_xlim(0, max(1, len(str(seq1))))
    ax.set_ylim(0, max(1, len(str(seq2))))
    ax.set_xlabel(label1, fontsize=9, color=MUTED_COLOR)
    ax.set_ylabel(label2, fontsize=9, color=MUTED_COLOR)
    ax.grid(True, color=GRID_COLOR, linewidth=0.4, zorder=0)
    style_axes(ax, hide=())
    if fx or rx:
        leg = ax.legend(loc='upper left', fontsize=8, framealpha=0.9,
                        markerscale=6)
        leg.get_frame().set_edgecolor(GRID_COLOR)

    subtitle = f'word size {word_size}'
    if truncated:
        subtitle += f' — first {MAX_DOTPLOT_POINTS:,} matches'
    set_title(ax, f'{title} ({subtitle})')
    return svg_markup(fig, title='Dot plot')
