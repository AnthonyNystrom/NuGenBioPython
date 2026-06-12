"""
UniGene Analysis Routes
Handles Bio.UniGene parsing and reading operations
"""
from flask import Blueprint, request, jsonify
from utils.upload_helpers import saved_upload, UploadError
from utils.request_helpers import clamp_int
from utils.unigene_helpers import unigene_parse, unigene_read

bp = Blueprint('unigene', __name__, url_prefix='/api/unigene')


@bp.route('/parse', methods=['POST'])
def parse():
    """Parse multiple UniGene records from file"""
    try:
        file = request.files.get('file')
        max_records = clamp_int(request.form.get('max_records'), 10, lo=1, hi=1000)

        with saved_upload(file) as temp_path:
            records = unigene_parse(temp_path, max_records)

        return jsonify({'success': True, 'records': records, 'count': len(records)})
    except UploadError as e:
        return jsonify({'success': False, 'error': str(e)})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/read', methods=['POST'])
def read():
    """Read single UniGene record with full details"""
    try:
        file = request.files.get('file')

        with saved_upload(file) as temp_path:
            record = unigene_read(temp_path)

        return jsonify({'success': True, 'record': record})
    except UploadError as e:
        return jsonify({'success': False, 'error': str(e)})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})
