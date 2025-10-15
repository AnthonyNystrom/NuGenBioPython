"""
Routes for BLAST operations
"""
from flask import Blueprint, request, jsonify, session
import json
import uuid
from utils.blast_helpers import (
    run_ncbi_blast, parse_blast_xml, format_blast_results,
    get_alignment_view, validate_sequence, get_blast_program_info,
    filter_results_by_evalue, filter_results_by_identity,
    get_result_statistics, extract_sequence_from_file
)

bp = Blueprint('blast', __name__, url_prefix='/api')

# In-memory storage for BLAST results (keyed by session-based ID)
# In production, use Redis or database
blast_results_cache = {}


@bp.route('/blast/search', methods=['POST'])
def blast_search():
    """Run BLAST search against NCBI"""
    try:
        # Check if file was uploaded
        uploaded_file = request.files.get('sequenceFile')
        sequence = request.form.get('sequence', '').strip()

        if uploaded_file and uploaded_file.filename:
            # Extract sequence from file
            file_content = uploaded_file.read().decode('utf-8')
            file_ext = uploaded_file.filename.rsplit('.', 1)[-1].lower()

            # Map file extension to format
            format_map = {
                'fasta': 'fasta', 'fa': 'fasta', 'fna': 'fasta', 'faa': 'fasta',
                'gb': 'genbank', 'gbk': 'genbank', 'genbank': 'genbank',
                'txt': 'txt'
            }
            file_format = format_map.get(file_ext, 'fasta')

            try:
                sequence = extract_sequence_from_file(file_content, file_format)
            except ValueError as e:
                return jsonify({'success': False, 'error': str(e)})

        program = request.form.get('program', 'blastn')
        database = request.form.get('database', 'nt')
        expect = float(request.form.get('expect', 10))
        hitlist_size = int(request.form.get('hitlist_size', 50))
        word_size = request.form.get('word_size')
        matrix_name = request.form.get('matrix_name')
        gap_costs = request.form.get('gap_costs')
        filter_low_complexity = request.form.get('filter', 'F')
        genetic_code = request.form.get('genetic_code')
        nucl_penalty = request.form.get('nucl_penalty')
        nucl_reward = request.form.get('nucl_reward')

        if not sequence:
            return jsonify({'success': False, 'error': 'No sequence provided'})

        # Get program info
        prog_info = get_blast_program_info(program)
        if not prog_info:
            return jsonify({'success': False, 'error': 'Invalid BLAST program'})

        # Validate sequence
        is_valid, error_msg = validate_sequence(sequence, prog_info['query_type'])
        if not is_valid:
            return jsonify({'success': False, 'error': error_msg})

        # Build kwargs for optional parameters
        kwargs = {
            'expect': expect,
            'hitlist_size': hitlist_size,
            'filter': filter_low_complexity
        }

        if word_size:
            kwargs['word_size'] = int(word_size)
        if matrix_name:
            kwargs['matrix_name'] = matrix_name
        if gap_costs:
            kwargs['gapcosts'] = gap_costs
        if genetic_code:
            kwargs['genetic_code'] = int(genetic_code)
        if nucl_penalty:
            kwargs['nucl_penalty'] = int(nucl_penalty)
        if nucl_reward:
            kwargs['nucl_reward'] = int(nucl_reward)

        # Run BLAST
        blast_xml = run_ncbi_blast(sequence, program, database, **kwargs)

        # Generate unique result ID
        if 'blast_result_id' not in session:
            session['blast_result_id'] = str(uuid.uuid4())

        result_id = session['blast_result_id']

        # Store results in memory cache (not in session cookie)
        blast_results_cache[result_id] = {
            'xml': blast_xml,
            'params': {
                'program': program,
                'database': database,
                'sequence_length': len(sequence.replace(' ', '').replace('\n', ''))
            }
        }

        # Parse and format results
        blast_records = parse_blast_xml(blast_xml)
        results = format_blast_results(blast_records, hitlist_size)

        # Get statistics
        stats = get_result_statistics(results)

        return jsonify({
            'success': True,
            'results': results,
            'statistics': stats,
            'num_hits': len(results)
        })

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/blast/alignment/<int:hit_index>', methods=['GET'])
def get_alignment(hit_index):
    """Get detailed alignment view for a specific hit"""
    try:
        result_id = session.get('blast_result_id')
        if not result_id or result_id not in blast_results_cache:
            return jsonify({'success': False, 'error': 'No BLAST results in session'})

        blast_xml = blast_results_cache[result_id]['xml']
        blast_records = parse_blast_xml(blast_xml)
        results = format_blast_results(blast_records)

        if hit_index >= len(results):
            return jsonify({'success': False, 'error': 'Invalid hit index'})

        hit = results[hit_index]
        alignment_text = get_alignment_view(hit)

        return jsonify({
            'success': True,
            'hit': hit,
            'alignment': alignment_text
        })

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/blast/filter', methods=['POST'])
def filter_results():
    """Filter BLAST results by criteria"""
    try:
        result_id = session.get('blast_result_id')
        if not result_id or result_id not in blast_results_cache:
            return jsonify({'success': False, 'error': 'No BLAST results in session'})

        blast_xml = blast_results_cache[result_id]['xml']
        max_evalue = float(request.form.get('max_evalue', 10))
        min_identity = float(request.form.get('min_identity', 0))

        blast_records = parse_blast_xml(blast_xml)
        results = format_blast_results(blast_records)

        # Apply filters
        if max_evalue < 10:
            results = filter_results_by_evalue(results, max_evalue)
        if min_identity > 0:
            results = filter_results_by_identity(results, min_identity)

        stats = get_result_statistics(results)

        return jsonify({
            'success': True,
            'results': results,
            'statistics': stats,
            'num_hits': len(results)
        })

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/blast/export', methods=['GET'])
def export_results():
    """Export BLAST results in various formats"""
    try:
        format_type = request.args.get('format', 'json')
        result_id = session.get('blast_result_id')

        if not result_id or result_id not in blast_results_cache:
            return jsonify({'success': False, 'error': 'No BLAST results in session'})

        blast_xml = blast_results_cache[result_id]['xml']
        params = blast_results_cache[result_id]['params']

        if format_type == 'xml':
            return blast_xml, 200, {'Content-Type': 'application/xml'}

        blast_records = parse_blast_xml(blast_xml)
        results = format_blast_results(blast_records)

        if format_type == 'json':
            return jsonify({
                'success': True,
                'results': results,
                'parameters': params
            })

        elif format_type == 'text':
            text_output = []

            text_output.append("BLAST Search Results")
            text_output.append("=" * 80)
            text_output.append(f"Program: {params.get('program', 'N/A')}")
            text_output.append(f"Database: {params.get('database', 'N/A')}")
            text_output.append(f"Query Length: {params.get('sequence_length', 'N/A')}")
            text_output.append(f"Total Hits: {len(results)}")
            text_output.append("=" * 80)
            text_output.append("")

            for i, hit in enumerate(results, 1):
                text_output.append(f"Hit {i}: {hit['accession']}")
                text_output.append(f"  Title: {hit['title']}")
                text_output.append(f"  Score: {hit['score']} bits, E-value: {hit['e_value']}")
                text_output.append(f"  Identity: {hit['identity_percent']}%, Coverage: {hit['coverage_percent']}%")
                text_output.append("")

            return '\n'.join(text_output), 200, {'Content-Type': 'text/plain'}

        else:
            return jsonify({'success': False, 'error': 'Unsupported format'})

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/blast/validate', methods=['POST'])
def validate_blast_sequence():
    """Validate sequence for BLAST"""
    try:
        sequence = request.form.get('sequence', '').strip()
        program = request.form.get('program', 'blastn')

        prog_info = get_blast_program_info(program)
        if not prog_info:
            return jsonify({'success': False, 'error': 'Invalid BLAST program'})

        is_valid, error_msg = validate_sequence(sequence, prog_info['query_type'])

        if is_valid:
            clean_seq = sequence.replace(' ', '').replace('\n', '').upper()
            return jsonify({
                'success': True,
                'valid': True,
                'length': len(clean_seq),
                'type': prog_info['query_type']
            })
        else:
            return jsonify({
                'success': True,
                'valid': False,
                'error': error_msg
            })

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/blast/programs', methods=['GET'])
def get_blast_programs():
    """Get list of available BLAST programs with descriptions"""
    programs = ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']
    program_info = [get_blast_program_info(p) for p in programs]

    return jsonify({
        'success': True,
        'programs': program_info
    })
