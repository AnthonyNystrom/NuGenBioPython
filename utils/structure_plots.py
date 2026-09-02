"""Ramachandran classification and plotting.

The previous in-route classifier could never report an outlier: its final
branch tested ``-180 <= phi <= 180 and -180 <= psi <= 180``, which is true for
every possible pair of angles, so anything not caught by the two preceding
rectangles fell into "allowed". Every structure scored zero outliers — the
opposite of what a validation plot is for.

Regions here are elliptical basins around the three canonical conformations.
They approximate the *general-case* favored and allowed contours; they are not
the residue-type-specific MolProbity contours, and are not a substitute for a
full validation package. They are, however, sufficient to place a residue in
the right basin and to flag genuinely strained geometry.
"""
import math

# (name, phi_centre, psi_centre, favored_radii, allowed_radii)
# Angles in degrees. psi comparisons wrap at +-180.
_BASINS = (
    ('beta',   -120.0,  135.0, (58.0, 50.0), (78.0, 72.0)),
    ('alphaR',  -63.0,  -43.0, (38.0, 32.0), (58.0, 52.0)),
    ('alphaL',   57.0,   40.0, (26.0, 26.0), (42.0, 42.0)),
    # Extended/bridge region between beta and alphaR, populated in real
    # structures and otherwise misread as an outlier.
    ('bridge', -100.0,   10.0, (35.0, 28.0), (55.0, 56.0)),
)


def _delta(a, b):
    """Shortest signed separation between two angles, in degrees."""
    return (a - b + 180.0) % 360.0 - 180.0


def _in_ellipse(phi, psi, phi_c, psi_c, radii):
    dphi = _delta(phi, phi_c)
    dpsi = _delta(psi, psi_c)
    rphi, rpsi = radii
    return (dphi / rphi) ** 2 + (dpsi / rpsi) ** 2 <= 1.0


def classify_phi_psi(phi, psi, resname=None):
    """Return 'favored', 'allowed' or 'outlier' for one residue.

    phi/psi are in degrees. Pass `resname` so glycine is handled correctly:
    with no side chain its accessible region is far broader and roughly
    symmetric through the origin, so positive-phi conformations that would be
    strained for any other residue are perfectly normal for Gly. Without this,
    a well-refined structure reports its glycines as outliers.
    """
    candidates = [(phi, psi)]
    if resname and str(resname).upper() == 'GLY':
        # Gly's distribution is approximately symmetric under
        # (phi, psi) -> (-phi, -psi); test the mirrored point too.
        candidates.append((-phi, -psi))

    for test_phi, test_psi in candidates:
        for _name, phi_c, psi_c, favored, _allowed in _BASINS:
            if _in_ellipse(test_phi, test_psi, phi_c, psi_c, favored):
                return 'favored'
    for test_phi, test_psi in candidates:
        for _name, phi_c, psi_c, _favored, allowed in _BASINS:
            if _in_ellipse(test_phi, test_psi, phi_c, psi_c, allowed):
                return 'allowed'
    return 'outlier'


def classify_all(phi_psi_data):
    """Classify a list of {'phi': deg, 'psi': deg} dicts.

    Returns (counts, annotated) where annotated is the same list with a
    'region' key added to each entry.
    """
    counts = {'favored': 0, 'allowed': 0, 'outlier': 0}
    annotated = []
    for entry in phi_psi_data:
        region = classify_phi_psi(entry['phi'], entry['psi'],
                                  entry.get('resname'))
        counts[region] += 1
        item = dict(entry)
        item['region'] = region
        annotated.append(item)
    return counts, annotated


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def render_ramachandran(annotated, title='Ramachandran plot'):
    """Render a Ramachandran scatter as an SVG data URL.

    `annotated` is the output of classify_all(): dicts carrying phi, psi and
    region. Basin contours are drawn behind the points so a reader can see
    *why* a residue was classified as it was, rather than having to trust a
    bare count.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse

    from utils.plot_helpers import (
        PALETTE, MUTED_COLOR, GRID_COLOR, figure_to_svg_data_url,
        set_title, style_axes,
    )

    fig, ax = plt.subplots(figsize=(6.4, 6.0))

    # Basin contours: allowed first (lighter), then favored on top.
    for _name, phi_c, psi_c, favored, allowed in _BASINS:
        ax.add_patch(Ellipse((phi_c, psi_c), allowed[0] * 2, allowed[1] * 2,
                             facecolor=GRID_COLOR, edgecolor='none',
                             alpha=0.55, zorder=1))
    for _name, phi_c, psi_c, favored, allowed in _BASINS:
        ax.add_patch(Ellipse((phi_c, psi_c), favored[0] * 2, favored[1] * 2,
                             facecolor=PALETTE['lightblue'], edgecolor='none',
                             alpha=0.45, zorder=2))

    # Glycine is scored against the mirrored basins as well (see
    # classify_phi_psi). Without drawing them, a Gly sitting at positive phi
    # reads as a "favored" point floating in empty space, which makes the
    # plot look wrong. Outline them so the classification is legible.
    if any(str(e.get('resname', '')).upper() == 'GLY' for e in annotated):
        for _name, phi_c, psi_c, _favored, allowed in _BASINS:
            ax.add_patch(Ellipse((-phi_c, -psi_c), allowed[0] * 2,
                                 allowed[1] * 2, facecolor='none',
                                 edgecolor=MUTED_COLOR, linewidth=0.6,
                                 linestyle=(0, (4, 3)), alpha=0.5, zorder=2))
        ax.plot([], [], color=MUTED_COLOR, linewidth=0.6,
                linestyle=(0, (4, 3)), alpha=0.7,
                label='Glycine-only region')

    ax.axhline(0, color=GRID_COLOR, linewidth=0.6, zorder=3)
    ax.axvline(0, color=GRID_COLOR, linewidth=0.6, zorder=3)

    styles = {
        'favored': (PALETTE['blue'], 18, 'Favored'),
        'allowed': (PALETTE['orange'], 22, 'Allowed'),
        'outlier': (PALETTE['red'], 46, 'Outlier'),
    }
    for region in ('favored', 'allowed', 'outlier'):
        pts = [(e['phi'], e['psi']) for e in annotated
               if e.get('region') == region]
        if not pts:
            continue
        color, size, label = styles[region]
        ax.scatter([p[0] for p in pts], [p[1] for p in pts],
                   s=size, c=color, label=f'{label} ({len(pts)})',
                   edgecolors='white' if region == 'outlier' else 'none',
                   linewidths=0.6, zorder=5 if region == 'outlier' else 4,
                   alpha=0.9)

    # Label outliers — they are the reason to look at this plot.
    for e in annotated:
        if e.get('region') == 'outlier' and e.get('residue'):
            ax.annotate(e['residue'], (e['phi'], e['psi']),
                        textcoords='offset points', xytext=(6, 4),
                        fontsize=7, color=PALETTE['red'], zorder=6)

    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xticks(range(-180, 181, 60))
    ax.set_yticks(range(-180, 181, 60))
    ax.set_xlabel('φ (degrees)', fontsize=9, color=MUTED_COLOR)
    ax.set_ylabel('ψ (degrees)', fontsize=9, color=MUTED_COLOR)
    ax.set_aspect('equal', adjustable='box')
    style_axes(ax, hide=())
    set_title(ax, title)
    leg = ax.legend(loc='lower right', fontsize=8, framealpha=0.9,
                    borderpad=0.6)
    leg.get_frame().set_edgecolor(GRID_COLOR)

    return figure_to_svg_data_url(fig)
