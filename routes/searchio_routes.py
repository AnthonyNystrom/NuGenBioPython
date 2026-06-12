"""
Routes for Bio.SearchIO operations
Handles parsing, reading, indexing, converting, filtering, and writing search results
"""
from flask import Blueprint, request, jsonify
from utils.upload_helpers import saved_upload, UploadError
from utils.searchio_helpers import (
    searchio_parse, searchio_read, searchio_index, searchio_convert,
    searchio_filter, searchio_write, get_supported_formats
)

bp = Blueprint('searchio', __name__, url_prefix='/api/searchio')


@bp.route('/parse', methods=['POST'])
def parse():
    """Parse search results file with multiple queries"""
    try:
        file = request.files.get('file')
        format_name = request.form.get('format', 'blast-xml')

        with saved_upload(file) as temp_path:
            results = searchio_parse(temp_path, format_name)

        return jsonify({
            'success': True,
            'results': results,
            'count': len(results)
        })

    except UploadError as e:
        return jsonify({'success': False, 'error': str(e)})
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        })


@bp.route('/read', methods=['POST'])
def read():
    """Read search results file with exactly one query"""
    try:
        file = request.files.get('file')
        format_name = request.form.get('format', 'blast-xml')

        with saved_upload(file) as temp_path:
            result = searchio_read(temp_path, format_name)

        return jsonify({
            'success': True,
            'result': result
        })

    except UploadError as e:
        return jsonify({'success': False, 'error': str(e)})
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        })


@bp.route('/index', methods=['POST'])
def index():
    """Create indexed access to search results"""
    try:
        file = request.files.get('file')
        format_name = request.form.get('format', 'blast-xml')
        limit = int(request.form.get('limit', 100))

        with saved_upload(file) as temp_path:
            results = searchio_index(temp_path, format_name, limit)

        return jsonify({
            'success': True,
            'results': results
        })

    except UploadError as e:
        return jsonify({'success': False, 'error': str(e)})
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        })


@bp.route('/convert', methods=['POST'])
def convert():
    """Convert search results from one format to another"""
    try:
        file = request.files.get('file')
        input_format = request.form.get('input_format', 'blast-xml')
        output_format = request.form.get('output_format', 'blast-tab')

        with saved_upload(file) as input_path:
            result = searchio_convert(input_path, input_format, output_format)

        return jsonify({
            'success': True,
            'result': result
        })

    except UploadError as e:
        return jsonify({'success': False, 'error': str(e)})
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        })


@bp.route('/filter', methods=['POST'])
def filter_results():
    """Filter search results by E-value, bit score, or identity"""
    try:
        file = request.files.get('file')
        format_name = request.form.get('format', 'blast-xml')

        # Get filter parameters
        evalue_threshold = request.form.get('evalue_threshold')
        bitscore_threshold = request.form.get('bitscore_threshold')
        min_identity = request.form.get('min_identity')

        # Convert to appropriate types
        evalue_threshold = float(evalue_threshold) if evalue_threshold else None
        bitscore_threshold = float(bitscore_threshold) if bitscore_threshold else None
        min_identity = float(min_identity) if min_identity else None

        with saved_upload(file) as temp_path:
            results = searchio_filter(temp_path, format_name, evalue_threshold, bitscore_threshold, min_identity)

        return jsonify({
            'success': True,
            'results': results,
            'count': len(results)
        })

    except UploadError as e:
        return jsonify({'success': False, 'error': str(e)})
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        })


@bp.route('/write', methods=['POST'])
def write():
    """Write search results to a different format"""
    try:
        file = request.files.get('file')
        input_format = request.form.get('input_format', 'blast-xml')
        output_format = request.form.get('output_format', 'blast-tab')
        max_queries = int(request.form.get('max_queries', 10))

        with saved_upload(file) as input_path:
            result = searchio_write(input_path, input_format, output_format, max_queries)

        return jsonify({
            'success': True,
            'result': result
        })

    except UploadError as e:
        return jsonify({'success': False, 'error': str(e)})
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        })


@bp.route('/formats', methods=['GET'])
def formats():
    """Get list of supported SearchIO formats"""
    try:
        formats = get_supported_formats()

        return jsonify({
            'success': True,
            'formats': formats
        })

    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        })
