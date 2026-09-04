"""
Routes for BLAST operations
"""
import logging
import os
import threading
import time

from flask import Blueprint, request, jsonify, session
import json
import uuid
from utils.blast_helpers import (
    run_ncbi_blast, parse_blast_xml, format_blast_results,
    get_alignment_view, validate_sequence, get_blast_program_info,
    filter_results_by_evalue, filter_results_by_identity,
    get_result_statistics, extract_sequence_from_file
)
from utils.request_helpers import (
    clamp_int, clamp_float, error_response, safe_error_message,
)

from utils.shared_store import SharedStore

log = logging.getLogger(__name__)
bp = Blueprint('blast', __name__, url_prefix='/api')

# Bounded in-memory BLAST result cache. Evicts oldest entries at the cap.
# For multi-worker deployments, replace with Redis or shared storage.
_BLAST_CACHE_MAX = int(os.environ.get('BLAST_CACHE_MAX', '64'))
blast_results_cache = SharedStore('blast', max_entries=_BLAST_CACHE_MAX)


def _store_blast_result(result_id, payload):
    blast_results_cache.set(result_id, payload)


def _collect_blast_request():
    """Parse and validate a BLAST request.

    Returns (params, error_message). Shared by the synchronous endpoint and
    the async submitter so the two cannot drift apart in what they accept.
    """
    uploaded_file = request.files.get('sequenceFile')
    sequence = request.form.get('sequence', '').strip()

    if uploaded_file and uploaded_file.filename:
        file_content = uploaded_file.read().decode('utf-8')
        file_ext = uploaded_file.filename.rsplit('.', 1)[-1].lower()
        format_map = {
            'fasta': 'fasta', 'fa': 'fasta', 'fna': 'fasta', 'faa': 'fasta',
            'gb': 'genbank', 'gbk': 'genbank', 'genbank': 'genbank',
            'txt': 'txt',
        }
        try:
            sequence = extract_sequence_from_file(
                file_content, format_map.get(file_ext, 'fasta'))
        except ValueError as exc:
            return None, str(exc)

    program = request.form.get('program', 'blastn')
    database = request.form.get('database', 'nt')
    hitlist_size = clamp_int(request.form.get('hitlist_size'), 50, lo=1, hi=500)

    if not sequence:
        return None, 'No sequence provided'

    prog_info = get_blast_program_info(program)
    if not prog_info:
        return None, 'Invalid BLAST program'

    is_valid, error_msg = validate_sequence(sequence, prog_info['query_type'])
    if not is_valid:
        return None, error_msg

    kwargs = {
        'expect': clamp_float(request.form.get('expect'), 10.0, lo=0.0, hi=1e6),
        'hitlist_size': hitlist_size,
        'filter': request.form.get('filter', 'F'),
    }
    for form_key, kwarg, caster in (
        ('word_size', 'word_size', int),
        ('matrix_name', 'matrix_name', str),
        ('gap_costs', 'gapcosts', str),
        ('genetic_code', 'genetic_code', int),
        ('nucl_penalty', 'nucl_penalty', int),
        ('nucl_reward', 'nucl_reward', int),
    ):
        value = request.form.get(form_key)
        if value:
            kwargs[kwarg] = caster(value)

    return {
        'sequence': sequence,
        'program': program,
        'database': database,
        'hitlist_size': hitlist_size,
        'kwargs': kwargs,
    }, None


def _run_blast_and_store(params, result_id):
    """Execute a BLAST search and return the formatted payload.

    Shared by the blocking endpoint and the background worker.
    """
    blast_xml = run_ncbi_blast(params['sequence'], params['program'],
                               params['database'], **params['kwargs'])
    _store_blast_result(result_id, {
        'xml': blast_xml,
        'params': {
            'program': params['program'],
            'database': params['database'],
            'sequence_length': len(params['sequence'].replace(' ', '')
                                   .replace('\n', '')),
        },
    })
    results = format_blast_results(parse_blast_xml(blast_xml),
                                   params['hitlist_size'])
    return {
        'results': results,
        'statistics': get_result_statistics(results),
        'num_hits': len(results),
    }


@bp.route('/blast/search', methods=['POST'])
def blast_search():
    """Run BLAST against NCBI, blocking until it completes.

    Kept for callers that want a single round trip. A BLAST can take minutes,
    so prefer /blast/submit + /blast/job/<id>, which does not hold a worker
    open for the whole search.
    """
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
        expect = clamp_float(request.form.get('expect'), 10.0, lo=0.0, hi=1e6)
        hitlist_size = clamp_int(request.form.get('hitlist_size'), 50, lo=1, hi=500)
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

        # Store results in bounded memory cache (not in session cookie)
        _store_blast_result(result_id, {
            'xml': blast_xml,
            'params': {
                'program': program,
                'database': database,
                'sequence_length': len(sequence.replace(' ', '').replace('\n', ''))
            }
        })

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
        return error_response(e, context='blast_routes.blast_search')



# ---------------------------------------------------------------------------
# Asynchronous submission
# ---------------------------------------------------------------------------
# A BLAST against NCBI routinely takes 30-180 seconds. Running it inline holds
# a worker open for the whole search, and any proxy or browser in between is
# free to time out first and leave the user with nothing. Submitting returns
# immediately with a job id the client polls.
#
# Jobs live in the shared store, so with REDIS_URL set the poll can be served
# by a different worker than the one that ran the search.

_JOB_TTL = int(os.environ.get('BLAST_JOB_TTL', '3600'))
blast_jobs = SharedStore('blast_jobs', max_entries=256, ttl=_JOB_TTL,
                         serializer='json')


def _set_job(job_id, **fields):
    job = blast_jobs.get(job_id) or {'job_id': job_id}
    job.update(fields)
    blast_jobs.set(job_id, job)
    return job


def _blast_worker(job_id, params, result_id):
    """Background thread body. Never raises into the thread runner."""
    _set_job(job_id, status='running', started_at=time.time())
    try:
        payload = _run_blast_and_store(params, result_id)
        _set_job(job_id, status='done', finished_at=time.time(), **payload)
    except Exception as exc:                     # noqa: BLE001 - reported below
        log.exception('blast job %s failed', job_id)
        _set_job(job_id, status='error', finished_at=time.time(),
                 error=safe_error_message(exc))


@bp.route('/blast/submit', methods=['POST'])
def blast_submit():
    """Start a BLAST search and return a job id immediately."""
    try:
        params, error = _collect_blast_request()
        if error:
            return jsonify({'success': False, 'error': error})

        job_id = str(uuid.uuid4())

        # Allocate the result id here, in the request context, and put it in
        # the session now — the worker has no request context, and
        # /blast/alignment and /blast/export read this key.
        result_id = session.get('blast_result_id') or str(uuid.uuid4())
        session['blast_result_id'] = result_id

        _set_job(job_id, status='queued', submitted_at=time.time(),
                 program=params['program'], database=params['database'],
                 sequence_length=len(params['sequence']))

        threading.Thread(target=_blast_worker,
                         args=(job_id, params, result_id),
                         name=f'blast-{job_id[:8]}',
                         daemon=True).start()

        return jsonify({'success': True, 'job_id': job_id, 'status': 'queued'})
    except Exception as e:
        return error_response(e, context='blast_routes.blast_submit')


@bp.route('/blast/job/<job_id>', methods=['GET'])
def blast_job_status(job_id):
    """Poll a submitted BLAST job.

    status is one of queued / running / done / error. Results are included
    once status is 'done'.
    """
    try:
        job = blast_jobs.get(job_id)
        if job is None:
            return jsonify({
                'success': False,
                'error': ('That BLAST job is unknown or has expired. Jobs are '
                          f'kept for {_JOB_TTL // 60} minutes.'),
            }), 404
        return jsonify({'success': True, **job})
    except Exception as e:
        return error_response(e, context='blast_routes.blast_job_status')

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
        return error_response(e, context='blast_routes.get_alignment')


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
        return error_response(e, context='blast_routes.filter_results')


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
        return error_response(e, context='blast_routes.export_results')


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
        return error_response(e, context='blast_routes.validate_blast_sequence')


@bp.route('/blast/programs', methods=['GET'])
def get_blast_programs():
    """Get list of available BLAST programs with descriptions"""
    programs = ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']
    program_info = [get_blast_program_info(p) for p in programs]

    return jsonify({
        'success': True,
        'programs': program_info
    })
