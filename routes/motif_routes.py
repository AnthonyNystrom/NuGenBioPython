"""
Routes for comprehensive motif analysis
"""
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

# In-memory storage for motifs (keyed by session-based ID)
motif_cache = {}


@bp.route('/motifs/create', methods=['POST'])
def create_motif():
    """Create motif from sequences"""
    try:
        data = request.json
        sequences = data.get('sequences', [])

        if not sequences or len(sequences) < 2:
            return jsonify({'success': False, 'error': 'At least 2 sequences required'})

        # Create motif
        m = create_motif_from_sequences(sequences)

        # Store in cache
        if 'motif_id' not in session:
            session['motif_id'] = str(uuid.uuid4())

        motif_id = session['motif_id']
        motif_cache[motif_id] = m

        # Get motif info
        info = get_motif_info(m)

        # Generate sequence logo
        logo_base64 = generate_sequence_logo(m)
        if logo_base64:
            info['logo'] = f'data:image/png;base64,{logo_base64}'

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
        motif_cache[motif_id] = m

        # Get motif info
        info = get_motif_info(m)

        # Generate logo
        logo_base64 = generate_sequence_logo(m)
        if logo_base64:
            info['logo'] = f'data:image/png;base64,{logo_base64}'

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
        data = request.json
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
        data = request.json
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
            info1['logo'] = f'data:image/png;base64,{logo1}'
        if logo2:
            info2['logo'] = f'data:image/png;base64,{logo2}'

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
            info['logo'] = f'data:image/png;base64,{logo_base64}'

        return jsonify({
            'success': True,
            'motif': info
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


def generate_sequence_logo(m, use_information_content=True):
    """
    Generate sequence logo with information content

    Args:
        m: Motif object
        use_information_content: Use IC for letter heights

    Returns:
        Base64 encoded PNG image
    """
    try:
        pwm = m.counts.normalize(pseudocounts=0.5)
        alphabet_list = ['A', 'C', 'G', 'T']
        colors = {'A': '#228B22', 'C': '#4169E1', 'G': '#FFA500', 'T': '#DC143C'}

        # Calculate information content per position
        ic_per_pos = []
        for i in range(len(m)):
            ic = 0
            for base in alphabet_list:
                p = pwm[base][i]
                if p > 0:
                    ic += p * np.log2(p / 0.25)  # Assuming uniform background
            ic_per_pos.append(ic)

        # Set figure size based on motif length
        width = max(8, min(20, len(m) * 0.6))
        fig, ax = plt.subplots(figsize=(width, 4))

        # Create stacked bar chart
        positions = np.arange(len(m))

        for i in range(len(m)):
            # Sort bases by frequency for this position
            base_freqs = [(base, pwm[base][i]) for base in alphabet_list]
            base_freqs.sort(key=lambda x: x[1])

            bottom = 0
            for base, freq in base_freqs:
                if freq > 0.01:  # Only show significant frequencies
                    if use_information_content:
                        height = freq * ic_per_pos[i]
                    else:
                        height = freq

                    ax.bar(i, height, bottom=bottom, color=colors[base],
                          width=0.9, edgecolor='white', linewidth=0.5)

                    # Add letter in center of bar if tall enough
                    if height > 0.15:
                        ax.text(i, bottom + height/2, base,
                               ha='center', va='center', fontsize=14,
                               fontweight='bold', color='white')

                    bottom += height

        # Style the plot
        ax.set_xlabel('Position', fontsize=12, fontweight='bold')
        if use_information_content:
            ax.set_ylabel('Information Content (bits)', fontsize=12, fontweight='bold')
            ax.set_ylim(0, 2.1)
        else:
            ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
            ax.set_ylim(0, 1.1)

        ax.set_title('Sequence Logo', fontsize=14, fontweight='bold', pad=15)
        ax.set_xticks(positions)
        ax.set_xticklabels([str(i+1) for i in positions])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', alpha=0.3, linestyle='--')

        # Save to base64
        img_buffer = BytesIO()
        plt.savefig(img_buffer, format='png', bbox_inches='tight', dpi=150, facecolor='white')
        img_buffer.seek(0)
        logo_base64 = base64.b64encode(img_buffer.getvalue()).decode()
        plt.close()

        return logo_base64
    except Exception as e:
        print(f"Logo generation error: {e}")
        return None
