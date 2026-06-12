"""
Routes for restriction enzyme analysis
"""
import logging
import math

from flask import Blueprint, request, jsonify, send_file
from io import StringIO, BytesIO
import csv
import json

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from dependencies import Seq, Restriction, SeqIO
from utils.plot_helpers import (
    PALETTE, LABEL_COLOR, MUTED_COLOR, AXIS_COLOR, GRID_COLOR,
    figure_to_svg_data_url, set_title, fmt_bp, nice_ticks,
)

log = logging.getLogger(__name__)
bp = Blueprint('restriction', __name__, url_prefix='/api')


# Distinct, accessible colors for up to ~10 enzymes. Cycles afterward.
_ENZYME_COLORS = [
    '#2563eb', '#dc2626', '#059669', '#ea580c', '#7c3aed',
    '#0891b2', '#ca8a04', '#db2777', '#0d9488', '#92400e',
]


def render_restriction_map(analysis, sequence_length):
    """Render a linear restriction map as an SVG data URL.

    `analysis` is the {enzyme: {...}} dict returned by restriction_analyze.
    Draws the sequence as a baseline with color-coded tick marks per enzyme.
    Labels are staggered across vertical rows to avoid collisions on crowded
    maps — this is the main advantage over the client-side HTML version.
    """
    enzymes = [(name, r) for name, r in analysis.items()
               if isinstance(r, dict) and not r.get('error')
               and r.get('cut_positions')]
    if not enzymes or sequence_length <= 0:
        return None

    # All cuts, with their source enzyme + color
    cuts = []
    for i, (name, r) in enumerate(enzymes):
        color = _ENZYME_COLORS[i % len(_ENZYME_COLORS)]
        for pos in r['cut_positions']:
            cuts.append({'pos': int(pos), 'enzyme': name, 'color': color})
    cuts.sort(key=lambda c: c['pos'])

    # Figure dimensions — width fixed, height depends on how many label rows
    # we need to avoid collisions.
    fig_w = 14
    # Estimate characters per label ~ 8. Use sequence_length as the x-extent
    # to decide whether labels collide.
    min_gap = sequence_length * 0.05  # ~5% of the span per label slot
    label_rows = 1
    occupancy = [[]]
    for c in cuts:
        name_w = max(len(c['enzyme']) * sequence_length * 0.008, min_gap)
        x_lo = c['pos'] - name_w / 2
        x_hi = c['pos'] + name_w / 2
        row = 0
        while row < len(occupancy):
            if all(x_hi < s or x_lo > e for s, e in occupancy[row]):
                break
            row += 1
        if row == len(occupancy):
            occupancy.append([])
            label_rows = len(occupancy)
        occupancy[row].append((x_lo, x_hi))
        c['label_row'] = row

    # Height: base + extra row per stacked label (~0.25in per row)
    fig_h = max(2.6, 2.0 + label_rows * 0.28)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=100)

    # Baseline — a thick rounded bar for the sequence
    bar_y = 0
    bar_h = 0.18
    ax.add_patch(Rectangle((0, bar_y - bar_h / 2), sequence_length, bar_h,
                           facecolor='#f1f5f9', edgecolor=AXIS_COLOR,
                           linewidth=0.6, zorder=2))

    # Tick marks per cut
    tick_half = 0.35
    for c in cuts:
        ax.plot([c['pos'], c['pos']],
                [bar_y - tick_half, bar_y + tick_half],
                color=c['color'], linewidth=1.4, zorder=3)
        # Label above the tick
        row_offset = 0.55 + c['label_row'] * 0.35
        label_y = bar_y + tick_half + row_offset
        ax.plot([c['pos'], c['pos']], [bar_y + tick_half, label_y - 0.12],
                color=c['color'], linewidth=0.5, alpha=0.5, zorder=2)
        ax.text(c['pos'], label_y,
                f"{c['enzyme']}\n{c['pos']:,}",
                ha='center', va='bottom', fontsize=7.5,
                color=c['color'], zorder=4)

    # Axis — position scale below the bar
    ticks = nice_ticks(sequence_length)
    ax.set_xticks(ticks)
    ax.set_xticklabels([fmt_bp(t) for t in ticks])
    ax.set_xlim(-sequence_length * 0.02, sequence_length * 1.02)
    ax.set_ylim(bar_y - tick_half - 0.6, bar_y + tick_half + 0.4 + label_rows * 0.35)
    ax.set_yticks([])
    for spine in ('top', 'right', 'left'):
        ax.spines[spine].set_visible(False)
    ax.spines['bottom'].set_color(AXIS_COLOR)
    ax.spines['bottom'].set_linewidth(0.6)
    ax.tick_params(axis='x', labelsize=8, colors=MUTED_COLOR, length=3, pad=2)
    ax.xaxis.grid(True, color=GRID_COLOR, linewidth=0.5, zorder=0)
    ax.set_xlabel('Position (bp)', fontsize=9, color=LABEL_COLOR, labelpad=4)
    set_title(ax, f'Restriction Map · {sequence_length:,} bp', pad=10)

    return figure_to_svg_data_url(fig)


@bp.route('/restriction/analyze', methods=['POST'])
def restriction_analyze():
    """Basic restriction enzyme analysis"""
    try:
        data = request.get_json(silent=True) or {}
        sequence = Seq.Seq(data.get('sequence', ''))
        enzymes = data.get('enzymes', ['EcoRI', 'BamHI', 'HindIII'])
        include_sequences = data.get('include_sequences', False)

        results = {}

        # Analyze with each enzyme
        for enzyme_name in enzymes:
            try:
                enzyme = getattr(Restriction, enzyme_name)
                cuts = enzyme.search(sequence)

                # Get enzyme properties
                overhang_type = 'blunt' if enzyme.is_blunt() else ('5-overhang' if enzyme.is_5overhang() else '3-overhang')

                results[enzyme_name] = {
                    'recognition_site': str(enzyme.site),
                    'cut_positions': cuts,
                    'number_of_cuts': len(cuts),
                    'overhang': enzyme.ovhg,
                    'overhang_seq': str(enzyme.ovhgseq) if hasattr(enzyme, 'ovhgseq') else '',
                    'overhang_type': overhang_type,
                    'is_blunt': enzyme.is_blunt(),
                    'is_5overhang': enzyme.is_5overhang(),
                    'is_3overhang': enzyme.is_3overhang(),
                    'fragments': []
                }

                # Calculate fragment sizes
                if cuts:
                    fragments = []
                    sorted_cuts = sorted(cuts)

                    # Add fragment from start to first cut
                    if sorted_cuts[0] > 0:
                        fragments.append(sorted_cuts[0])

                    # Add fragments between consecutive cuts
                    for i in range(len(sorted_cuts) - 1):
                        fragment_size = sorted_cuts[i + 1] - sorted_cuts[i]
                        if fragment_size > 0:
                            fragments.append(fragment_size)

                    # Add fragment from last cut to end
                    if sorted_cuts[-1] < len(sequence):
                        fragments.append(len(sequence) - sorted_cuts[-1])

                    results[enzyme_name]['fragments'] = fragments

                    # Add actual fragment sequences if requested
                    if include_sequences:
                        try:
                            fragment_seqs = enzyme.catalyze(sequence)
                            results[enzyme_name]['fragment_sequences'] = [str(s) for s in fragment_seqs]
                        except Exception:
                            results[enzyme_name]['fragment_sequences'] = []
                else:
                    # No cuts means one fragment of the entire sequence
                    results[enzyme_name]['fragments'] = [len(sequence)]
                    if include_sequences:
                        results[enzyme_name]['fragment_sequences'] = [str(sequence)]

            except AttributeError:
                results[enzyme_name] = {'error': f'Enzyme {enzyme_name} not found'}

        map_svg = render_restriction_map(results, len(sequence))
        return jsonify({'success': True, 'analysis': results, 'map': map_svg})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/restriction/advanced_analysis', methods=['POST'])
def advanced_analysis():
    """Advanced restriction enzyme analysis using Bio.Restriction.Analysis"""
    try:
        from Bio.Restriction import Analysis, RestrictionBatch, AllEnzymes

        data = request.get_json(silent=True) or {}
        sequence = Seq.Seq(data.get('sequence', ''))
        filter_type = data.get('filter', 'all')  # all, unique, non_cutters, blunt, 5overhang, 3overhang
        min_cuts = data.get('min_cuts', 0)
        max_cuts = data.get('max_cuts', None)
        enzyme_list = data.get('enzymes', None)

        # Create restriction batch
        if enzyme_list and len(enzyme_list) > 0:
            # Use specified enzymes
            rb = RestrictionBatch([e for e in enzyme_list if hasattr(Restriction, e)])
        else:
            # Use all enzymes
            rb = AllEnzymes

        # Create analysis
        analyzer = Analysis(rb, sequence)

        result = {
            'sequence_length': len(sequence),
            'total_enzymes_tested': len(rb),
            'filters_applied': []
        }

        # Apply filters
        if filter_type == 'unique' or (min_cuts == 1 and max_cuts == 1):
            enzymes_dict = analyzer.with_N_sites(1)
            result['filters_applied'].append('Enzymes cutting exactly once')
        elif filter_type == 'non_cutters':
            enzymes_dict = analyzer.without_site()
            result['filters_applied'].append('Enzymes that do not cut')
        elif filter_type == 'blunt':
            enzymes_dict = analyzer.blunt()
            result['filters_applied'].append('Enzymes producing blunt ends')
        elif filter_type == '5overhang':
            enzymes_dict = analyzer.with_sites()
            # Filter for 5' overhangs
            enzymes_dict = {e: sites for e, sites in enzymes_dict.items() if e.is_5overhang()}
            result['filters_applied'].append("Enzymes with 5' overhang")
        elif filter_type == '3overhang':
            enzymes_dict = analyzer.with_sites()
            # Filter for 3' overhangs
            enzymes_dict = {e: sites for e, sites in enzymes_dict.items() if e.is_3overhang()}
            result['filters_applied'].append("Enzymes with 3' overhang")
        else:
            # Get all cutting enzymes
            enzymes_dict = analyzer.with_sites()
            result['filters_applied'].append('All enzymes that cut')

        # Apply cut count filter
        if min_cuts > 0 or max_cuts is not None:
            filtered = {}
            for enzyme, sites in enzymes_dict.items():
                num_cuts = len(sites)
                if num_cuts >= min_cuts and (max_cuts is None or num_cuts <= max_cuts):
                    filtered[enzyme] = sites
            enzymes_dict = filtered
            result['filters_applied'].append(f'Cut count: {min_cuts}-{max_cuts if max_cuts else "∞"}')

        # Format results
        enzymes_result = []
        for enzyme, sites in sorted(enzymes_dict.items(), key=lambda x: len(x[1]), reverse=True):
            enzyme_info = {
                'name': str(enzyme),
                'site': str(enzyme.site),
                'cuts': len(sites),
                'positions': sorted(sites),
                'overhang': enzyme.ovhg,
                'overhang_type': 'blunt' if enzyme.is_blunt() else ('5-overhang' if enzyme.is_5overhang() else '3-overhang'),
                'is_blunt': enzyme.is_blunt()
            }
            enzymes_result.append(enzyme_info)

        result['enzymes'] = enzymes_result
        result['enzyme_count'] = len(enzymes_result)

        return jsonify({'success': True, 'result': result})

    except Exception as e:
        log.exception('restriction/advanced_analysis failed')
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/restriction/list_enzymes', methods=['GET'])
def list_restriction_enzymes():
    """List restriction enzymes with filters"""
    try:
        from Bio.Restriction import AllEnzymes, CommOnly

        filter_type = request.args.get('filter', 'common')  # common, commercial, all
        search = request.args.get('search', '').upper()
        site_size = request.args.get('site_size', None)
        overhang_filter = request.args.get('overhang', None)  # blunt, 5, 3

        # Get enzyme list based on filter
        if filter_type == 'all':
            enzyme_list = list(AllEnzymes)
        elif filter_type == 'commercial':
            enzyme_list = list(CommOnly)
        else:  # common
            common_names = [
                'EcoRI', 'BamHI', 'HindIII', 'XbaI', 'SalI', 'PstI', 'SmaI', 'KpnI', 'SacI', 'XhoI',
                'BglII', 'NcoI', 'NotI', 'SpeI', 'ApaI', 'ClaI', 'DraI', 'EcoRV', 'HaeIII', 'HhaI',
                'PvuII', 'NdeI', 'SphI', 'BstEII', 'MluI', 'NsiI', 'AvrII', 'SbfI', 'XmaI', 'BspEI'
            ]
            enzyme_list = [getattr(Restriction, name) for name in common_names if hasattr(Restriction, name)]

        enzyme_info = []
        for enzyme in enzyme_list:
            try:
                enzyme_name = str(enzyme)

                # Apply search filter
                if search and search not in enzyme_name.upper():
                    continue

                # Get site size
                site = str(enzyme.site)
                current_site_size = len(site)

                # Apply site size filter
                if site_size and int(site_size) != current_site_size:
                    continue

                # Apply overhang filter
                if overhang_filter == 'blunt' and not enzyme.is_blunt():
                    continue
                if overhang_filter == '5' and not enzyme.is_5overhang():
                    continue
                if overhang_filter == '3' and not enzyme.is_3overhang():
                    continue

                # Get overhang info
                overhang_type = 'blunt' if enzyme.is_blunt() else ('5-overhang' if enzyme.is_5overhang() else '3-overhang')

                # Get suppliers
                suppliers = []
                try:
                    if hasattr(enzyme, 'supplier_list'):
                        suppliers = enzyme.supplier_list()
                    elif hasattr(enzyme, 'suppl'):
                        suppliers = list(enzyme.suppl)
                except Exception:
                    pass

                enzyme_info.append({
                    'name': enzyme_name,
                    'site': site,
                    'site_size': current_site_size,
                    'overhang': enzyme.ovhg,
                    'overhang_seq': str(enzyme.ovhgseq) if hasattr(enzyme, 'ovhgseq') else '',
                    'overhang_type': overhang_type,
                    'is_blunt': enzyme.is_blunt(),
                    'is_5overhang': enzyme.is_5overhang(),
                    'is_3overhang': enzyme.is_3overhang(),
                    'suppliers': suppliers,
                    'commercial': len(suppliers) > 0
                })
            except Exception as e:
                continue

        return jsonify({
            'success': True,
            'enzymes': enzyme_info,
            'count': len(enzyme_info),
            'filter': filter_type
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/restriction/enzyme_details/<enzyme_name>', methods=['GET'])
def enzyme_details(enzyme_name):
    """Get detailed information about a specific enzyme"""
    try:
        enzyme = getattr(Restriction, enzyme_name)

        # Get suppliers
        suppliers = []
        try:
            if hasattr(enzyme, 'supplier_list'):
                suppliers = enzyme.supplier_list()
            elif hasattr(enzyme, 'suppl'):
                suppliers = list(enzyme.suppl)
        except Exception:
            pass

        # Get isoschizomers (enzymes with same recognition site)
        isoschizomers = []
        try:
            if hasattr(enzyme, 'isoschizomers'):
                isoschizomers = [str(e) for e in enzyme.isoschizomers()]
        except Exception:
            pass

        details = {
            'name': enzyme_name,
            'site': str(enzyme.site),
            'site_size': len(str(enzyme.site)),
            'overhang': enzyme.ovhg,
            'overhang_seq': str(enzyme.ovhgseq) if hasattr(enzyme, 'ovhgseq') else '',
            'overhang_type': 'blunt' if enzyme.is_blunt() else ('5-overhang' if enzyme.is_5overhang() else '3-overhang'),
            'is_blunt': enzyme.is_blunt(),
            'is_5overhang': enzyme.is_5overhang(),
            'is_3overhang': enzyme.is_3overhang(),
            'cut_twice': enzyme.cut_twice(),
            'suppliers': suppliers,
            'supplier_count': len(suppliers),
            'isoschizomers': isoschizomers,
            'commercial': len(suppliers) > 0
        }

        return jsonify({'success': True, 'enzyme': details})
    except AttributeError:
        return jsonify({'success': False, 'error': f'Enzyme {enzyme_name} not found'})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/restriction/export', methods=['POST'])
def export_results():
    """Export restriction analysis results"""
    try:
        data = request.get_json(silent=True) or {}
        results = data.get('results', {})
        format_type = data.get('format', 'csv')  # csv, json, tsv

        if format_type == 'json':
            output = BytesIO()
            output.write(json.dumps(results, indent=2).encode())
            output.seek(0)
            return send_file(output, mimetype='application/json', as_attachment=True, download_name='restriction_analysis.json')

        elif format_type in ['csv', 'tsv']:
            delimiter = '\t' if format_type == 'tsv' else ','
            output = StringIO()
            writer = csv.writer(output, delimiter=delimiter)

            # Write header
            writer.writerow(['Enzyme', 'Recognition Site', 'Cuts', 'Cut Positions', 'Fragment Sizes', 'Overhang Type', 'Overhang Seq'])

            # Write data
            analysis = results.get('analysis', {})
            for enzyme_name, enzyme_data in analysis.items():
                if 'error' not in enzyme_data:
                    writer.writerow([
                        enzyme_name,
                        enzyme_data.get('recognition_site', ''),
                        enzyme_data.get('number_of_cuts', 0),
                        ';'.join(map(str, enzyme_data.get('cut_positions', []))),
                        ';'.join(map(str, enzyme_data.get('fragments', []))),
                        enzyme_data.get('overhang_type', ''),
                        enzyme_data.get('overhang_seq', '')
                    ])

            output.seek(0)
            mimetype = 'text/tab-separated-values' if format_type == 'tsv' else 'text/csv'
            filename = f'restriction_analysis.{format_type}'

            return send_file(
                BytesIO(output.getvalue().encode()),
                mimetype=mimetype,
                as_attachment=True,
                download_name=filename
            )

        return jsonify({'success': False, 'error': 'Unsupported format'})

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/restriction/compatible_ends', methods=['POST'])
def compatible_ends():
    """Find enzymes with compatible cohesive ends"""
    try:
        data = request.get_json(silent=True) or {}
        enzyme_names = data.get('enzymes', [])

        if len(enzyme_names) < 2:
            return jsonify({'success': False, 'error': 'Please provide at least 2 enzymes'})

        enzymes = []
        for name in enzyme_names:
            try:
                enzymes.append(getattr(Restriction, name))
            except AttributeError:
                pass

        compatible_pairs = []

        # Check all pairs for compatible ends
        for i, enz1 in enumerate(enzymes):
            for enz2 in enzymes[i+1:]:
                # Compatible if same overhang sequence
                if not enz1.is_blunt() and not enz2.is_blunt():
                    ovhg1 = str(enz1.ovhgseq) if hasattr(enz1, 'ovhgseq') else ''
                    ovhg2 = str(enz2.ovhgseq) if hasattr(enz2, 'ovhgseq') else ''

                    if ovhg1 and ovhg2 and ovhg1 == ovhg2:
                        compatible_pairs.append({
                            'enzyme1': str(enz1),
                            'enzyme2': str(enz2),
                            'overhang_seq': ovhg1,
                            'type': 'identical'
                        })

        return jsonify({
            'success': True,
            'compatible_pairs': compatible_pairs,
            'count': len(compatible_pairs)
        })

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/restriction/upload_sequence', methods=['POST'])
def upload_sequence():
    """Parse uploaded sequence file and return DNA sequence"""
    try:
        if 'file' not in request.files:
            return jsonify({'success': False, 'error': 'No file provided'})

        file = request.files['file']

        if file.filename == '':
            return jsonify({'success': False, 'error': 'No file selected'})

        # Get file extension
        filename = file.filename.lower()

        # Read file content
        content = file.read().decode('utf-8')

        # Parse based on file type
        sequence = None
        seq_info = {}

        if filename.endswith(('.fasta', '.fa', '.fna')):
            # Parse FASTA
            try:
                from io import StringIO
                handle = StringIO(content)
                record = SeqIO.read(handle, 'fasta')
                sequence = str(record.seq).upper()
                seq_info = {
                    'format': 'FASTA',
                    'id': record.id,
                    'description': record.description,
                    'length': len(sequence)
                }
            except Exception as e:
                return jsonify({'success': False, 'error': f'Invalid FASTA format: {str(e)}'})

        elif filename.endswith(('.gb', '.gbk', '.genbank')):
            # Parse GenBank
            try:
                from io import StringIO
                handle = StringIO(content)
                record = SeqIO.read(handle, 'genbank')
                sequence = str(record.seq).upper()
                seq_info = {
                    'format': 'GenBank',
                    'id': record.id,
                    'name': record.name,
                    'description': record.description,
                    'length': len(sequence),
                    'features': len(record.features) if hasattr(record, 'features') else 0
                }
            except Exception as e:
                return jsonify({'success': False, 'error': f'Invalid GenBank format: {str(e)}'})

        elif filename.endswith('.embl'):
            # Parse EMBL
            try:
                from io import StringIO
                handle = StringIO(content)
                record = SeqIO.read(handle, 'embl')
                sequence = str(record.seq).upper()
                seq_info = {
                    'format': 'EMBL',
                    'id': record.id,
                    'description': record.description,
                    'length': len(sequence)
                }
            except Exception as e:
                return jsonify({'success': False, 'error': f'Invalid EMBL format: {str(e)}'})

        elif filename.endswith('.txt'):
            # Plain text - just extract DNA sequence
            # Remove whitespace and newlines
            sequence = ''.join(content.split()).upper()
            seq_info = {
                'format': 'Plain Text',
                'length': len(sequence)
            }
        else:
            return jsonify({'success': False, 'error': 'Unsupported file format. Please use FASTA, GenBank, EMBL, or plain text (.txt)'})

        # Validate it's a DNA sequence
        if sequence:
            # Remove any non-ATGC characters and check if reasonable
            clean_seq = ''.join(c for c in sequence if c in 'ATGC')

            if len(clean_seq) == 0:
                return jsonify({'success': False, 'error': 'No valid DNA sequence found in file'})

            if len(clean_seq) / len(sequence) < 0.9:
                return jsonify({'success': False, 'error': 'File contains too many non-DNA characters. Expected only A, T, G, C.'})

            sequence = clean_seq

        return jsonify({
            'success': True,
            'sequence': sequence,
            'info': seq_info,
            'filename': file.filename
        })

    except Exception as e:
        log.exception('restriction/upload_sequence failed')
        return jsonify({'success': False, 'error': str(e)})
