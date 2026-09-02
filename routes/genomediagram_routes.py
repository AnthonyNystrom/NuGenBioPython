"""
Routes for genome diagram visualization.

The rendering is matplotlib-based (not BioPython's ReportLab pipeline) because
the ReportLab output is print-oriented: it pins to a page size, wastes most of
the vertical space, and — with multi-fragment linear mode — duplicates track
names across fragments. The custom renderer below targets web display:
a single wide track band, proper position axis, arrow-shaped features,
and smart label placement with leader lines when features crowd.

File upload (GenBank/EMBL/FASTA) still uses BioPython for parsing only.
"""
from flask import Blueprint, request, jsonify

from utils.request_helpers import error_response, safe_error_message
from io import StringIO
import math

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle, Polygon, FancyArrowPatch, PathPatch
from matplotlib.path import Path
import numpy as np

from dependencies import SeqIO, Seq
from utils.plot_helpers import (
    PALETTE, EDGE_COLOR, TITLE_COLOR, LABEL_COLOR, AXIS_COLOR, GRID_COLOR,
    resolve_color, figure_to_svg_data_url, fmt_bp as _fmt,
)

bp = Blueprint('genomediagram', __name__, url_prefix='/api')


# --------------------------------------------------------------------------
# Shared rendering helpers
# --------------------------------------------------------------------------

# Palette, role constants, color resolution, and the bp axis formatter are
# shared with every other diagram route via utils.plot_helpers — this module
# used to carry byte-identical copies of all of them, which meant a styling
# change had to be made twice to keep the app's diagrams consistent.
BASELINE_COLOR = '#cbd5e1'


def _draw_arrow_feature(ax, start, end, y, height, color, strand,
                        edge=EDGE_COLOR, lw=0.6):
    """Draw a feature as a rectangle with an arrow head pointing in the
    strand direction. For strand==0 (unstranded), draws a plain rectangle.
    """
    width = end - start
    if width <= 0:
        return
    # Arrow head occupies up to 30% of the feature but is capped in absolute
    # coords so short features don't become all-arrow.
    head = min(width * 0.30, max(width * 0.05, 1))
    body_half = height / 2.0

    if strand > 0:
        body_x = start
        body_w = max(width - head, 0)
        # Body
        ax.add_patch(Rectangle((body_x, y - body_half), body_w, height,
                               facecolor=color, edgecolor=edge, linewidth=lw,
                               joinstyle='miter', zorder=3))
        # Arrow head polygon
        ax.add_patch(Polygon([
            (body_x + body_w, y - body_half),
            (body_x + body_w, y + body_half),
            (end,             y),
        ], facecolor=color, edgecolor=edge, linewidth=lw,
           joinstyle='miter', zorder=3))
    elif strand < 0:
        body_w = max(width - head, 0)
        body_x = start + head
        ax.add_patch(Rectangle((body_x, y - body_half), body_w, height,
                               facecolor=color, edgecolor=edge, linewidth=lw,
                               joinstyle='miter', zorder=3))
        ax.add_patch(Polygon([
            (body_x, y - body_half),
            (body_x, y + body_half),
            (start,  y),
        ], facecolor=color, edgecolor=edge, linewidth=lw,
           joinstyle='miter', zorder=3))
    else:
        ax.add_patch(Rectangle((start, y - body_half), width, height,
                               facecolor=color, edgecolor=edge, linewidth=lw,
                               joinstyle='miter', zorder=3))


def _place_feature_labels(ax, features, y, height, genome_length,
                          font_size=8.0, above=True):
    """Place feature labels above or below the track, alternating into two
    rows if any would collide. Returns list of (x, label_y) pairs used so
    the caller can adjust ylim.
    """
    if not features:
        return 0
    feats = sorted(features, key=lambda f: (f['start'] + f['end']) / 2)
    rows_used = 1
    occupancy = [[]]  # list of rows; each row = list of (x_min, x_max)
    # Conservative character-width estimate, tuned for a ~13in wide figure
    # at 100dpi with default sans-serif at ~8pt. 0.009 * genome_length per
    # character includes inter-label padding so neighbors don't touch.
    char_w = max(genome_length * 0.009, 1)
    min_gap = genome_length * 0.004  # always leave some breathing room

    positions = []
    for f in feats:
        name = f.get('name', '')
        if not name:
            positions.append(None)
            continue
        mid = (f['start'] + f['end']) / 2
        est_w = len(name) * char_w + min_gap
        x_min = mid - est_w / 2
        x_max = mid + est_w / 2
        row = 0
        while row < len(occupancy):
            if all(x_max < s or x_min > e for (s, e) in occupancy[row]):
                break
            row += 1
        if row == len(occupancy):
            occupancy.append([])
            rows_used = len(occupancy)
        occupancy[row].append((x_min, x_max))
        positions.append((mid, row))

    # Draw labels. Row 0 sits closest to the track; higher rows stack outward.
    row_gap = height * 1.25
    for f, pos in zip(feats, positions):
        if pos is None:
            continue
        mid, row = pos
        if above:
            ly = y + height / 2 + height * 0.45 + row * row_gap
            va = 'bottom'
        else:
            ly = y - height / 2 - height * 0.45 - row * row_gap
            va = 'top'
        # Leader line when label is in an outer row (improves readability)
        if row > 0:
            if above:
                ax.plot([mid, mid], [y + height / 2, ly - height * 0.1],
                        color=AXIS_COLOR, linewidth=0.4, zorder=2)
            else:
                ax.plot([mid, mid], [y - height / 2, ly + height * 0.1],
                        color=AXIS_COLOR, linewidth=0.4, zorder=2)
        ax.text(mid, ly, f.get('name', ''),
                ha='center', va=va, fontsize=font_size,
                color=LABEL_COLOR, zorder=4, clip_on=False)

    return rows_used


def _setup_linear_axes(ax, genome_length, n_tracks, label_rows_above,
                       label_rows_below, track_height, track_gap):
    """Configure the linear-diagram axes with a clean position scale."""
    left_pad = genome_length * 0.02
    right_pad = genome_length * 0.02
    ax.set_xlim(-left_pad, genome_length + right_pad)

    # Vertical layout: tracks stack from top to bottom starting at y=top.
    top_label_space = max(label_rows_above, 1) * track_height * 1.25 + track_height
    bottom_label_space = max(label_rows_below, 1) * track_height * 1.25 + track_height
    total_tracks_h = n_tracks * track_height + (n_tracks - 1) * track_gap
    total = total_tracks_h + top_label_space + bottom_label_space
    ax.set_ylim(-bottom_label_space, total_tracks_h + top_label_space)

    # X ticks at sensible intervals
    ticks = _nice_ticks(genome_length)
    ax.set_xticks(ticks)
    ax.set_xticklabels([_fmt(t) for t in ticks])
    ax.tick_params(axis='x', labelsize=8, colors=LABEL_COLOR,
                   length=3, pad=2)
    ax.set_yticks([])
    for spine in ('top', 'right', 'left'):
        ax.spines[spine].set_visible(False)
    ax.spines['bottom'].set_color(AXIS_COLOR)
    ax.spines['bottom'].set_linewidth(0.6)
    ax.xaxis.grid(True, color=GRID_COLOR, linewidth=0.5, zorder=0)


def _nice_ticks(genome_length, target_count=8):
    """Choose tick positions at a round step covering the genome."""
    if genome_length <= 0:
        return [0]
    raw = genome_length / target_count
    magnitude = 10 ** math.floor(math.log10(raw))
    for m in (1, 2, 2.5, 5, 10):
        step = m * magnitude
        if genome_length / step <= target_count * 1.5:
            break
    ticks = []
    v = 0
    while v <= genome_length + 1e-6:
        ticks.append(int(round(v)))
        v += step
    return ticks


# --------------------------------------------------------------------------
# Linear + circular renderers (feature tracks)
# --------------------------------------------------------------------------

def render_linear_feature_diagram(genome_length, tracks, title='',
                                  cross_links=None, show_track_labels=True):
    """Render a linear genome diagram with one or more feature tracks.

    `tracks` is a list of dicts: {name, features, height?}
    Each feature is {name, start, end, color?, strand?, sigil?}.
    `cross_links` is a list of {track1, start1, end1, track2, start2, end2, color?}.
    """
    n = max(1, len(tracks))
    track_height = 0.55
    track_gap = 0.65

    # Pre-compute label row counts so we can size the figure.
    label_rows_above = 1
    label_rows_below = 1
    # Recompute with placeholder later, but we need some estimate for height.
    fig_w = 13
    base_h = 1.4 + n * 1.4
    fig, ax = plt.subplots(figsize=(fig_w, base_h), dpi=100)

    # Track positions: track 0 at top, track n-1 at bottom.
    # y_center[i] = top_y - i * (track_height + track_gap)
    top_y = (n - 1) * (track_height + track_gap)
    y_centers = [top_y - i * (track_height + track_gap) for i in range(n)]

    # Track axis lines
    for i, track in enumerate(tracks):
        y = y_centers[i]
        # Baseline through the middle of the track
        ax.plot([0, genome_length], [y, y],
                color=BASELINE_COLOR, linewidth=0.8, zorder=1)
        if show_track_labels and track.get('name'):
            ax.text(-genome_length * 0.01, y, track['name'],
                    ha='right', va='center',
                    fontsize=9, fontweight='600',
                    color=LABEL_COLOR, clip_on=False)

    # Draw cross-links behind features so feature patches sit on top.
    if cross_links:
        for link in cross_links:
            t1 = int(link.get('track1', 1)) - 1
            t2 = int(link.get('track2', 2)) - 1
            if not (0 <= t1 < n) or not (0 <= t2 < n):
                continue
            y1 = y_centers[t1]
            y2 = y_centers[t2]
            # If track1 is above track2, link goes from bottom of t1 to top of t2
            sign = 1 if y1 > y2 else -1
            y1_edge = y1 - sign * track_height / 2
            y2_edge = y2 + sign * track_height / 2
            color = resolve_color(link.get('color'), '#93c5fd')
            verts = [
                (link['start1'], y1_edge),
                (link['end1'],   y1_edge),
                (link['end2'],   y2_edge),
                (link['start2'], y2_edge),
            ]
            ax.add_patch(Polygon(verts, facecolor=color, edgecolor='none',
                                 alpha=0.4, zorder=1.5))

    # Feature patches + labels
    max_above = 1
    max_below = 1
    for i, track in enumerate(tracks):
        y = y_centers[i]
        feats = track.get('features', [])
        # Resolve colors up-front
        prepared = []
        for f in feats:
            start = int(f.get('start', 0))
            end = int(f.get('end', start + 1))
            if end < start:
                start, end = end, start
            prepared.append({
                'name': f.get('name', ''),
                'start': max(0, min(start, genome_length)),
                'end': max(0, min(end, genome_length)),
                'color': resolve_color(f.get('color')),
                'strand': int(f.get('strand', 1) or 0),
            })
        # Draw features
        for f in prepared:
            _draw_arrow_feature(ax, f['start'], f['end'], y,
                                track_height, f['color'], f['strand'])
        # Labels: above for forward/unstranded, below for reverse — but only
        # when a track contains both strands. Otherwise put everything above.
        has_fwd = any(f['strand'] >= 0 for f in prepared)
        has_rev = any(f['strand'] < 0 for f in prepared)
        if has_fwd and has_rev:
            fwd = [f for f in prepared if f['strand'] >= 0]
            rev = [f for f in prepared if f['strand'] < 0]
            above = _place_feature_labels(ax, fwd, y, track_height,
                                          genome_length, above=True)
            below = _place_feature_labels(ax, rev, y, track_height,
                                          genome_length, above=False)
        else:
            above = _place_feature_labels(ax, prepared, y, track_height,
                                          genome_length, above=True)
            below = 0
        if i == 0:
            max_above = max(max_above, above)
        if i == n - 1:
            max_below = max(max_below, below)

    _setup_linear_axes(ax, genome_length, n, max_above, max_below,
                       track_height, track_gap)

    ax.set_xlabel('Position (bp)', fontsize=9, color=LABEL_COLOR, labelpad=6)
    if title:
        ax.set_title(title, fontsize=12, fontweight='600',
                     color=TITLE_COLOR, pad=14, loc='left')

    return figure_to_svg_data_url(fig, pad_inches=0.25)


def render_circular_feature_diagram(genome_length, tracks, title='',
                                    cross_links=None):
    """Render a circular genome diagram. Tracks stack concentrically."""
    n = max(1, len(tracks))
    fig, ax = plt.subplots(figsize=(9, 9), dpi=100,
                           subplot_kw={'projection': 'polar'})

    # Outer tracks drawn at larger radii. Reserve inner 35% for a clean hub
    # and leave a small gap between tracks.
    r_outer = 1.0
    r_inner_bound = 0.35
    track_span = (r_outer - r_inner_bound) / n
    track_gap = track_span * 0.18
    track_width = track_span - track_gap

    # Baseline rings
    for i in range(n):
        r_mid = r_outer - track_span * i - track_width / 2
        thetas = np.linspace(0, 2 * np.pi, 500)
        ax.plot(thetas, np.full_like(thetas, r_mid),
                color=BASELINE_COLOR, linewidth=0.6, zorder=1)

    # Features as polar wedges
    for i, track in enumerate(tracks):
        r_top = r_outer - track_span * i
        r_bot = r_top - track_width
        for f in track.get('features', []):
            start = int(f.get('start', 0))
            end = int(f.get('end', start + 1))
            if end < start:
                start, end = end, start
            color = resolve_color(f.get('color'))
            strand = int(f.get('strand', 1) or 0)
            t_start = 2 * np.pi * start / genome_length
            t_end = 2 * np.pi * end / genome_length
            if t_end < t_start:
                t_end += 2 * np.pi

            # Arrow head: 15% of the arc, capped
            head_arc = min((t_end - t_start) * 0.25, 0.04)
            if strand < 0:
                head_arc = -head_arc
                body_t0 = t_start - head_arc
                body_t1 = t_end
            else:
                body_t0 = t_start
                body_t1 = t_end - head_arc if strand > 0 else t_end

            # Body polygon (trapezoid in polar)
            thetas_body = np.linspace(body_t0, body_t1, 40)
            r_outer_arr = np.full_like(thetas_body, r_top)
            r_inner_arr = np.full_like(thetas_body, r_bot)
            theta_ring = np.concatenate([thetas_body, thetas_body[::-1]])
            r_ring = np.concatenate([r_outer_arr, r_inner_arr[::-1]])
            ax.fill(theta_ring, r_ring, facecolor=color,
                    edgecolor=EDGE_COLOR, linewidth=0.5, zorder=3)

            # Arrow head (small polygon)
            if strand != 0:
                if strand > 0:
                    tip_theta = t_end
                    base_theta = body_t1
                else:
                    tip_theta = t_start
                    base_theta = body_t0
                head_r = (r_top + r_bot) / 2
                ax.fill(
                    [base_theta, base_theta, tip_theta],
                    [r_top,      r_bot,      head_r],
                    facecolor=color, edgecolor=EDGE_COLOR, linewidth=0.5,
                    zorder=3)

            # Label outside the ring for the outermost track only.
            # Labels are drawn radially (pointing outward from the centre);
            # labels whose natural rotation would make them read upside-
            # down get flipped by 180° so every label reads left-to-right
            # from the viewer's POV.
            if i == 0 and f.get('name'):
                mid_theta = (t_start + t_end) / 2
                r_label = r_outer + 0.05
                deg = math.degrees(mid_theta) % 360
                # With theta_zero='N' and theta_direction=-1, the radial-
                # outward screen angle (CCW from horizontal) is (90 - deg).
                rotation = (90 - deg) % 360
                # Rotations in (90°, 270°] put text upside-down or
                # reading-backwards. Flip by 180° and anchor at the
                # opposite side so the text grows outward and reads
                # normally.
                if 90 < rotation <= 270:
                    rotation = (rotation + 180) % 360
                    ha = 'right'
                else:
                    ha = 'left'
                ax.text(mid_theta, r_label, f['name'],
                        ha=ha, va='center',
                        fontsize=7.5, color=LABEL_COLOR, rotation=rotation,
                        rotation_mode='anchor', zorder=4)

    # Position markers on the inner hub
    for frac, pos in [(0, 0), (0.25, genome_length // 4),
                      (0.5, genome_length // 2), (0.75, 3 * genome_length // 4)]:
        theta = 2 * np.pi * frac
        ax.plot([theta, theta], [r_inner_bound - 0.02, r_inner_bound + 0.01],
                color=AXIS_COLOR, linewidth=0.6, zorder=2)
        ax.text(theta, r_inner_bound - 0.06, _fmt(pos),
                ha='center', va='center',
                fontsize=7.5, color=LABEL_COLOR, zorder=4)

    # Center text — total length
    ax.text(0, 0, f'{_fmt(genome_length)}\n{genome_length:,} bp',
            ha='center', va='center', fontsize=10, color=LABEL_COLOR,
            fontweight='600')

    # Track legend at top-left (shows track names)
    if any(t.get('name') for t in tracks):
        legend_handles = []
        for t in tracks:
            if t.get('name'):
                # Use first feature's color as indicator; fallback to gray.
                col = resolve_color(
                    (t.get('features') or [{}])[0].get('color'),
                    '#64748b')
                legend_handles.append(
                    mpatches.Patch(facecolor=col, edgecolor=EDGE_COLOR,
                                   linewidth=0.5, label=t['name']))
        if legend_handles:
            ax.legend(handles=legend_handles, loc='upper left',
                      bbox_to_anchor=(-0.12, 1.08), fontsize=8,
                      frameon=False, labelcolor=LABEL_COLOR)

    ax.set_ylim(0, r_outer + 0.22)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_rticks([])
    ax.set_xticks([])
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    if title:
        fig.suptitle(title, fontsize=13, fontweight='600',
                     color=TITLE_COLOR, y=0.97)

    return figure_to_svg_data_url(fig, pad_inches=0.25)


# --------------------------------------------------------------------------
# Data-track renderer (GC content, GC skew, custom)
# --------------------------------------------------------------------------

def render_data_tracks(genome_length, graphs, title='', diagram_type='linear'):
    """Render one or more data graphs (GC content, GC skew, custom)."""
    n = max(1, len(graphs))
    if diagram_type == 'circular':
        fig, ax = plt.subplots(figsize=(9, 9), dpi=100,
                               subplot_kw={'projection': 'polar'})
        r_outer = 1.0
        r_inner = 0.4
        span = (r_outer - r_inner) / n
        for i, g in enumerate(graphs):
            data = g['data']
            if not data:
                continue
            color = resolve_color(g.get('color'), '#2563eb')
            style = g.get('style', 'line')
            thetas = np.linspace(0, 2 * np.pi, len(data), endpoint=False)
            arr = np.array(data, dtype=float)
            # Normalize to [0, 1] relative to this track's band
            lo, hi = float(np.nanmin(arr)), float(np.nanmax(arr))
            rng = hi - lo if hi > lo else 1
            norm = (arr - lo) / rng
            band_outer = r_outer - span * i
            band_inner = band_outer - span * 0.9
            radii = band_inner + norm * (band_outer - band_inner)
            if style == 'bar':
                widths = (2 * np.pi / len(data)) * 0.9
                ax.bar(thetas, radii - band_inner, width=widths,
                       bottom=band_inner, color=color, edgecolor='none',
                       alpha=0.8, zorder=3)
            else:
                ax.plot(np.append(thetas, thetas[0]),
                        np.append(radii, radii[0]),
                        color=color, linewidth=1.2, zorder=3)
                ax.fill_between(np.append(thetas, thetas[0]),
                                np.append(radii, radii[0]),
                                band_inner, color=color, alpha=0.15, zorder=2)
            # Label
            ax.text(np.pi, band_outer + 0.02, g.get('label', 'data'),
                    ha='center', va='bottom', fontsize=8, color=LABEL_COLOR)
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.set_rticks([])
        ax.set_xticks([])
        ax.grid(False)
        for s in ax.spines.values():
            s.set_visible(False)
    else:
        fig, axes = plt.subplots(n, 1, figsize=(13, 1.5 + n * 1.6),
                                 dpi=100, sharex=True)
        if n == 1:
            axes = [axes]
        for ax, g in zip(axes, graphs):
            data = g['data']
            color = resolve_color(g.get('color'), '#2563eb')
            style = g.get('style', 'line')
            label = g.get('label', 'data')
            if not data:
                ax.text(0.5, 0.5, 'No data', transform=ax.transAxes,
                        ha='center', va='center', color=LABEL_COLOR)
            else:
                positions = np.linspace(0, genome_length, len(data))
                arr = np.array(data, dtype=float)
                if style == 'bar':
                    width = max(genome_length / len(data) * 0.9, 1)
                    ax.bar(positions, arr, width=width, color=color,
                           edgecolor='none', alpha=0.85, zorder=3)
                elif style == 'heat':
                    heat = arr.reshape(1, -1)
                    ax.imshow(heat, aspect='auto',
                              extent=(0, genome_length, -0.5, 0.5),
                              cmap='viridis', vmin=float(np.nanmin(arr)),
                              vmax=float(np.nanmax(arr)), zorder=3)
                    ax.set_ylim(-0.5, 0.5)
                else:
                    ax.plot(positions, arr, color=color, linewidth=1.4,
                            zorder=3)
                    ax.fill_between(positions, arr, float(np.nanmin(arr)),
                                    color=color, alpha=0.12, zorder=2)
            ax.set_xlim(0, genome_length)
            ax.set_ylabel(label, fontsize=9, color=LABEL_COLOR)
            ax.tick_params(axis='both', labelsize=8, colors=LABEL_COLOR)
            for spine in ('top', 'right'):
                ax.spines[spine].set_visible(False)
            ax.spines['left'].set_color(AXIS_COLOR)
            ax.spines['left'].set_linewidth(0.6)
            ax.spines['bottom'].set_color(AXIS_COLOR)
            ax.spines['bottom'].set_linewidth(0.6)
            ax.grid(True, color=GRID_COLOR, linewidth=0.5, zorder=0)
        axes[-1].set_xlabel('Position (bp)', fontsize=9, color=LABEL_COLOR)
        ticks = _nice_ticks(genome_length)
        axes[-1].set_xticks(ticks)
        axes[-1].set_xticklabels([_fmt(t) for t in ticks])

    if title:
        fig.suptitle(title, fontsize=12, fontweight='600',
                     color=TITLE_COLOR, y=0.98)
    fig.tight_layout(rect=(0, 0, 1, 0.96) if title else None)

    return figure_to_svg_data_url(fig, pad_inches=0.25)


# ============================================================================
# FILE UPLOAD - Parse GenBank/EMBL/FASTA files
# ============================================================================

@bp.route('/genomediagram/upload_file', methods=['POST'])
def upload_file():
    """Parse uploaded sequence file and extract features."""
    try:
        if 'file' not in request.files:
            return jsonify({'success': False, 'error': 'No file provided'})

        file = request.files['file']
        if file.filename == '':
            return jsonify({'success': False, 'error': 'No file selected'})

        filename = file.filename.lower()
        content = file.read().decode('utf-8')

        features = []
        seq_info = {}

        if filename.endswith(('.gb', '.gbk', '.genbank')):
            handle = StringIO(content)
            record = SeqIO.read(handle, 'genbank')
            seq_info = {
                'format': 'GenBank',
                'id': record.id,
                'name': record.name,
                'description': record.description,
                'length': len(record.seq),
                'features': len(record.features),
            }

            color_map = {
                'CDS': 'blue', 'gene': 'blue',
                'promoter': 'green', 'regulatory': 'green',
                'repeat_region': 'orange', 'repeat': 'orange',
                'terminator': 'red',
                'tRNA': 'purple', 'rRNA': 'purple',
            }

            for feature in record.features:
                if hasattr(feature.location, 'start') and hasattr(feature.location, 'end'):
                    feature_type = feature.type
                    feature_name = feature_type
                    if 'gene' in feature.qualifiers:
                        feature_name = feature.qualifiers['gene'][0]
                    elif 'label' in feature.qualifiers:
                        feature_name = feature.qualifiers['label'][0]
                    elif 'product' in feature.qualifiers:
                        feature_name = feature.qualifiers['product'][0]

                    features.append({
                        'name': feature_name,
                        'start': int(feature.location.start),
                        'end': int(feature.location.end),
                        'type': feature_type,
                        'strand': feature.location.strand if hasattr(feature.location, 'strand') else 1,
                        'color': color_map.get(feature_type, 'blue'),
                    })

        elif filename.endswith('.embl'):
            handle = StringIO(content)
            record = SeqIO.read(handle, 'embl')
            seq_info = {
                'format': 'EMBL',
                'id': record.id,
                'name': record.name,
                'description': record.description,
                'length': len(record.seq),
                'features': len(record.features),
            }
            for feature in record.features:
                if hasattr(feature.location, 'start') and hasattr(feature.location, 'end'):
                    feature_type = feature.type
                    feature_name = feature_type
                    if 'gene' in feature.qualifiers:
                        feature_name = feature.qualifiers['gene'][0]
                    features.append({
                        'name': feature_name,
                        'start': int(feature.location.start),
                        'end': int(feature.location.end),
                        'type': feature_type,
                        'strand': feature.location.strand if hasattr(feature.location, 'strand') else 1,
                        'color': 'blue',
                    })

        elif filename.endswith(('.fasta', '.fa', '.fna')):
            handle = StringIO(content)
            record = SeqIO.read(handle, 'fasta')
            seq_info = {
                'format': 'FASTA',
                'id': record.id,
                'description': record.description,
                'length': len(record.seq),
                'features': 0,
            }

        else:
            return jsonify({'success': False, 'error': 'Unsupported file format'})

        return jsonify({
            'success': True,
            'filename': file.filename,
            'info': seq_info,
            'features': features,
        })

    except Exception as e:
        return jsonify({'success': False, 'error': f'File parsing error: {safe_error_message(e)}'})


# ============================================================================
# BASIC DIAGRAM
# ============================================================================

@bp.route('/genomediagram/create', methods=['POST'])
def create_genome_diagram():
    try:
        data = request.get_json(silent=True) or {}
        genome_length = int(data.get('genome_length', 10000))
        features = data.get('features', [])
        diagram_type = data.get('diagram_type', 'linear')
        title = data.get('title', 'Genome Diagram')

        tracks = [{'name': 'Features', 'features': features}]
        if diagram_type == 'circular':
            url = render_circular_feature_diagram(genome_length, tracks, title)
        else:
            url = render_linear_feature_diagram(genome_length, tracks, title,
                                                show_track_labels=False)

        return jsonify({'success': True, 'diagram': url})

    except Exception as e:
        return error_response(e, context='genomediagram_routes.create_genome_diagram')


# ============================================================================
# MULTI-TRACK
# ============================================================================

@bp.route('/genomediagram/create_multitrack', methods=['POST'])
def create_multitrack_diagram():
    try:
        data = request.get_json(silent=True) or {}
        genome_length = int(data.get('genome_length', 10000))
        tracks_data = data.get('tracks', [])
        diagram_type = data.get('diagram_type', 'linear')
        title = data.get('title', 'Multi-Track Genome')

        # Normalize to list of {name, features}. Preserve input order; caller
        # chooses the stacking order.
        tracks = []
        for t in tracks_data:
            tracks.append({
                'name': t.get('name') or f"Track {t.get('track_number', len(tracks)+1)}",
                'features': t.get('features', []),
            })

        if not tracks:
            return jsonify({'success': False, 'error': 'No tracks provided'})

        if diagram_type == 'circular':
            url = render_circular_feature_diagram(genome_length, tracks, title)
        else:
            url = render_linear_feature_diagram(genome_length, tracks, title)

        return jsonify({'success': True, 'diagram': url})

    except Exception as e:
        return error_response(e, context='genomediagram_routes.create_multitrack_diagram')


# ============================================================================
# DATA TRACKS
# ============================================================================

@bp.route('/genomediagram/create_data_tracks', methods=['POST'])
def create_data_tracks():
    try:
        data = request.get_json(silent=True) or {}
        genome_length = int(data.get('genome_length', 10000))
        sequence = data.get('sequence', '')
        graphs = data.get('graphs', [])
        title = data.get('title', 'Genome Data Visualization')
        diagram_type = data.get('diagram_type', 'linear')

        seq = sequence.upper() if sequence else ''

        prepared = []
        for g in graphs:
            gtype = g.get('type', 'gc_content')
            window = max(int(g.get('window', 1000)), 10)
            color = g.get('color', 'blue')
            style = g.get('style', 'line')

            data_points = []
            label = 'Data'
            if gtype == 'gc_content' and seq:
                label = 'GC Content'
                for i in range(0, len(seq), window):
                    w = seq[i:i + window]
                    if len(w) > 0:
                        gc = (w.count('G') + w.count('C')) / len(w)
                        data_points.append(gc)
            elif gtype == 'gc_skew' and seq:
                label = 'GC Skew'
                for i in range(0, len(seq), window):
                    w = seq[i:i + window]
                    g_count = w.count('G')
                    c_count = w.count('C')
                    total = g_count + c_count
                    data_points.append((g_count - c_count) / total if total else 0)
            elif gtype == 'custom':
                label = 'Custom Data'
                data_points = list(g.get('custom_data') or [])

            prepared.append({
                'label': label, 'data': data_points,
                'color': color, 'style': style,
            })

        if not prepared:
            return jsonify({'success': False, 'error': 'No graphs provided'})

        url = render_data_tracks(genome_length, prepared, title, diagram_type)
        return jsonify({'success': True, 'diagram': url})

    except Exception as e:
        return error_response(e, context='genomediagram_routes.create_data_tracks')


# ============================================================================
# ADVANCED FEATURES
# ============================================================================

@bp.route('/genomediagram/create_advanced', methods=['POST'])
def create_advanced_diagram():
    try:
        data = request.get_json(silent=True) or {}
        genome_length = int(data.get('genome_length', 10000))
        features = data.get('features', [])
        cross_links = data.get('cross_links', [])
        title = data.get('title', 'Advanced Feature Styling')
        diagram_type = data.get('diagram_type', 'linear')

        # Split by strand into Forward / Reverse tracks so cross-strand
        # features read cleanly.
        fwd = [f for f in features if int(f.get('strand', 1) or 0) >= 0]
        rev = [f for f in features if int(f.get('strand', 1) or 0) < 0]
        tracks = []
        if fwd:
            tracks.append({'name': 'Forward Strand', 'features': fwd})
        if rev:
            tracks.append({'name': 'Reverse Strand', 'features': rev})
        if not tracks:
            tracks = [{'name': 'Features', 'features': features}]

        if diagram_type == 'circular':
            url = render_circular_feature_diagram(genome_length, tracks, title,
                                                  cross_links=cross_links)
        else:
            url = render_linear_feature_diagram(genome_length, tracks, title,
                                                cross_links=cross_links)

        return jsonify({'success': True, 'diagram': url})

    except Exception as e:
        return error_response(e, context='genomediagram_routes.create_advanced_diagram')


# ============================================================================
# EXPORT
# ============================================================================

@bp.route('/genomediagram/export', methods=['POST'])
def export_diagram():
    """Pass-through exporter: the diagram is already a data URL, so we just
    hand it back. Clients that want PDF/EPS can print the SVG from the
    browser; inline conversion would require more heavy deps.
    """
    try:
        data = request.get_json(silent=True) or {}
        diagram_data = data.get('diagram_data', '')
        export_format = data.get('format', 'svg').lower()

        if not diagram_data:
            return jsonify({'success': False, 'error': 'No diagram data provided'})

        # For SVG we can always hand the data URL back. For PNG we convert.
        if export_format == 'png' and diagram_data.startswith('data:image/svg'):
            # Decode SVG, rasterize via matplotlib using svg2paths is heavy; instead,
            # emit the SVG as-is with a .svg filename so browsers/tools can open it.
            return jsonify({
                'success': True,
                'file_data': diagram_data,
                'note': 'PNG export requires browser-side rasterization; SVG returned.'
            })

        return jsonify({'success': True, 'file_data': diagram_data})

    except Exception as e:
        return error_response(e, context='genomediagram_routes.export_diagram')
