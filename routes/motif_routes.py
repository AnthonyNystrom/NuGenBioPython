"""
Routes for comprehensive motif analysis
"""
import os
from collections import OrderedDict

from flask import Blueprint, request, jsonify, session
import base64
from io import BytesIO
import numpy as np
import uuid

from Bio import motifs
from Bio.Seq import Seq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from utils.motif_helpers import (
    create_motif_from_sequences, parse_motif_from_file, get_motif_info,
    search_motif_advanced, compare_motifs, export_motif, calculate_motif_statistics
)

bp = Blueprint('motifs', __name__, url_prefix='/api')

# Bounded in-memory motif cache. Evicts oldest entries at the cap.
_MOTIF_CACHE_MAX = int(os.environ.get('MOTIF_CACHE_MAX', '64'))
motif_cache = OrderedDict()


def _store_motif(motif_id, motif):
    if motif_id in motif_cache:
        motif_cache.move_to_end(motif_id)
    motif_cache[motif_id] = motif
    while len(motif_cache) > _MOTIF_CACHE_MAX:
        motif_cache.popitem(last=False)


@bp.route('/motifs/create', methods=['POST'])
def create_motif():
    """Create motif from sequences"""
    try:
        data = request.get_json(silent=True) or {}
        sequences = data.get('sequences', [])

        if not sequences or len(sequences) < 2:
            return jsonify({'success': False, 'error': 'At least 2 sequences required'})

        # Create motif
        m = create_motif_from_sequences(sequences)

        # Store in cache
        if 'motif_id' not in session:
            session['motif_id'] = str(uuid.uuid4())

        motif_id = session['motif_id']
        _store_motif(motif_id, m)

        # Get motif info
        info = get_motif_info(m)

        # Generate sequence logo
        logo_base64 = generate_sequence_logo(m)
        if logo_base64:
            info['logo'] = logo_base64

        return jsonify({
            'success': True,
            'motif': info
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/motifs/upload', methods=['POST'])
def upload_motif():
    """Upload motif from file"""
    try:
        uploaded_file = request.files.get('motifFile')
        file_format = request.form.get('format', 'jaspar')

        if not uploaded_file or not uploaded_file.filename:
            return jsonify({'success': False, 'error': 'No file uploaded'})

        file_content = uploaded_file.read().decode('utf-8')

        # Parse motif
        result = parse_motif_from_file(file_content, file_format)

        # Handle single motif or list
        if isinstance(result, list):
            if len(result) == 0:
                return jsonify({'success': False, 'error': 'No motifs found in file'})
            m = result[0]  # Use first motif
            multiple = len(result) > 1
        else:
            m = result
            multiple = False

        # Store in cache
        if 'motif_id' not in session:
            session['motif_id'] = str(uuid.uuid4())

        motif_id = session['motif_id']
        _store_motif(motif_id, m)

        # Get motif info
        info = get_motif_info(m)

        # Generate logo
        logo_base64 = generate_sequence_logo(m)
        if logo_base64:
            info['logo'] = logo_base64

        return jsonify({
            'success': True,
            'motif': info,
            'multiple_motifs': multiple
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/motifs/search', methods=['POST'])
def search_motif():
    """Search for motif in sequence using PSSM scoring"""
    try:
        data = request.get_json(silent=True) or {}
        sequence = data.get('sequence', '')
        threshold_type = data.get('threshold_type', 'rel')  # 'abs' or 'rel'
        threshold_value = float(data.get('threshold', 0.7))
        search_rc = data.get('search_rc', True)

        if not sequence:
            return jsonify({'success': False, 'error': 'No sequence provided'})

        # Get motif from cache
        motif_id = session.get('motif_id')
        if not motif_id or motif_id not in motif_cache:
            return jsonify({'success': False, 'error': 'No motif created. Create a motif first.'})

        m = motif_cache[motif_id]

        # Search
        matches = search_motif_advanced(m, sequence, threshold_type, threshold_value)

        # Calculate statistics
        stats = calculate_motif_statistics(m, sequence)

        return jsonify({
            'success': True,
            'matches': matches,
            'total_matches': len(matches),
            'statistics': stats
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/motifs/compare', methods=['POST'])
def compare_motifs_route():
    """Compare current motif with another"""
    try:
        data = request.get_json(silent=True) or {}
        sequences2 = data.get('sequences', [])

        if not sequences2 or len(sequences2) < 2:
            return jsonify({'success': False, 'error': 'At least 2 sequences required for second motif'})

        # Get first motif from cache
        motif_id = session.get('motif_id')
        if not motif_id or motif_id not in motif_cache:
            return jsonify({'success': False, 'error': 'No motif created. Create a motif first.'})

        m1 = motif_cache[motif_id]

        # Create second motif
        m2 = create_motif_from_sequences(sequences2)

        # Compare
        comparison = compare_motifs(m1, m2)

        # Get info for both motifs
        info1 = get_motif_info(m1)
        info2 = get_motif_info(m2)

        # Generate logos
        logo1 = generate_sequence_logo(m1)
        logo2 = generate_sequence_logo(m2)

        if logo1:
            info1['logo'] = logo1
        if logo2:
            info2['logo'] = logo2

        return jsonify({
            'success': True,
            'comparison': comparison,
            'motif1': info1,
            'motif2': info2
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/motifs/export', methods=['GET'])
def export_motif_route():
    """Export current motif"""
    try:
        format_type = request.args.get('format', 'jaspar')

        # Get motif from cache
        motif_id = session.get('motif_id')
        if not motif_id or motif_id not in motif_cache:
            return jsonify({'success': False, 'error': 'No motif to export'})

        m = motif_cache[motif_id]

        # Export
        exported = export_motif(m, format_type)

        return jsonify({
            'success': True,
            'content': exported,
            'format': format_type
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/motifs/info', methods=['GET'])
def get_motif_info_route():
    """Get detailed information about current motif"""
    try:
        # Get motif from cache
        motif_id = session.get('motif_id')
        if not motif_id or motif_id not in motif_cache:
            return jsonify({'success': False, 'error': 'No motif created'})

        m = motif_cache[motif_id]

        # Get comprehensive info
        info = get_motif_info(m)

        # Generate logo
        logo_base64 = generate_sequence_logo(m)
        if logo_base64:
            info['logo'] = logo_base64

        return jsonify({
            'success': True,
            'motif': info
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


def generate_sequence_logo(m, use_information_content=True):
    """Generate a proper scaled-glyph sequence logo.

    Uses the `logomaker` library when available so each letter glyph is
    scaled by its information-content contribution (the standard Schneider
    & Stephens 1990 convention). Falls back to a stacked bar chart with
    lettered segments if logomaker isn't installed.

    Returns: SVG data URL string (or PNG fallback on the bar path).
    """
    try:
        pwm = m.counts.normalize(pseudocounts=0.5)
        alphabet_list = ['A', 'C', 'G', 'T']

        # IC per position (bits, uniform 0.25 background)
        ic_per_pos = []
        for i in range(len(m)):
            ic = 0.0
            for base in alphabet_list:
                p = pwm[base][i]
                if p > 0:
                    ic += p * np.log2(p / 0.25)
            ic_per_pos.append(max(0.0, ic))

        # Figure width: ~0.55 in per position but bounded. Height fixed at
        # 2.6 in — enough for 2-bit y-axis without wasting space.
        width = max(6, min(22, 0.55 * len(m) + 2))
        height = 2.8

        try:
            import logomaker
            # logomaker consumes a pandas DataFrame with positions as rows
            # and letters as columns. Values are the per-letter information
            # content (freq * IC).
            import pandas as pd
            rows = []
            for i in range(len(m)):
                row = {}
                for base in alphabet_list:
                    freq = pwm[base][i]
                    row[base] = freq * ic_per_pos[i] if use_information_content else freq
                rows.append(row)
            df = pd.DataFrame(rows)

            color_scheme = {
                'A': '#059669',  # green
                'C': '#2563eb',  # blue
                'G': '#ca8a04',  # gold
                'T': '#dc2626',  # red
            }

            fig, ax = plt.subplots(figsize=(width, height), dpi=100)
            fig.patch.set_alpha(0.0)
            ax.set_facecolor('none')
            logomaker.Logo(
                df, ax=ax,
                color_scheme=color_scheme,
                baseline_width=0.0,
                show_spines=False,
            )
            ax.set_xlabel('Position', fontsize=9, color='#64748b')
            ax.set_ylabel('Bits' if use_information_content else 'Frequency',
                          fontsize=9, color='#64748b')
            ax.set_xticks(range(len(m)))
            ax.set_xticklabels([str(i + 1) for i in range(len(m))])
            ax.tick_params(axis='both', labelsize=8, colors='#94a3b8')
            ax.set_ylim(0, 2.05 if use_information_content else 1.02)
            ax.spines['left'].set_color('#cbd5e1')
            ax.spines['bottom'].set_color('#cbd5e1')
            ax.spines['left'].set_linewidth(0.6)
            ax.spines['bottom'].set_linewidth(0.6)

            from io import StringIO as _SIO
            buf = _SIO()
            fig.savefig(buf, format='svg', bbox_inches='tight',
                        facecolor='none', transparent=True, pad_inches=0.2)
            plt.close(fig)
            return 'data:image/svg+xml;base64,' + \
                base64.b64encode(buf.getvalue().encode('utf-8')).decode()

        except ImportError:
            # Fallback: stacked bar chart with lettered segments — better
            # contrast than the previous version (dark letter on white
            # boxes instead of white letter on saturated color).
            colors = {'A': '#059669', 'C': '#2563eb',
                      'G': '#ca8a04', 'T': '#dc2626'}
            fig, ax = plt.subplots(figsize=(width, height), dpi=100)
            fig.patch.set_alpha(0.0)
            ax.set_facecolor('none')
            for i in range(len(m)):
                base_freqs = [(b, pwm[b][i]) for b in alphabet_list]
                base_freqs.sort(key=lambda x: x[1])
                bottom = 0
                for base, freq in base_freqs:
                    if freq <= 0.01:
                        continue
                    h = freq * ic_per_pos[i] if use_information_content else freq
                    ax.bar(i, h, bottom=bottom, color=colors[base],
                           width=0.9, edgecolor='none', linewidth=0)
                    if h > 0.18:
                        ax.text(i, bottom + h / 2, base,
                                ha='center', va='center', fontsize=13,
                                fontweight='700', color='white')
                    bottom += h
            ax.set_xlabel('Position', fontsize=9, color='#64748b')
            ax.set_ylabel('Bits' if use_information_content else 'Frequency',
                          fontsize=9, color='#64748b')
            ax.set_xticks(range(len(m)))
            ax.set_xticklabels([str(i + 1) for i in range(len(m))])
            ax.tick_params(axis='both', labelsize=8, colors='#94a3b8')
            ax.set_ylim(0, 2.05 if use_information_content else 1.02)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_color('#cbd5e1')
            ax.spines['bottom'].set_color('#cbd5e1')
            buf = BytesIO()
            fig.savefig(buf, format='png', bbox_inches='tight',
                        facecolor='none', transparent=True, dpi=150, pad_inches=0.2)
            plt.close(fig)
            return 'data:image/png;base64,' + \
                base64.b64encode(buf.getvalue()).decode()

    except Exception:
        return None
