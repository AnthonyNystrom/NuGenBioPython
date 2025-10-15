"""
UniGene Analysis Routes
Handles Bio.UniGene parsing and reading operations
"""
from flask import Blueprint, request, jsonify, current_app
import os
from werkzeug.utils import secure_filename
from utils.unigene_helpers import unigene_parse, unigene_read

bp = Blueprint('unigene', __name__, url_prefix='/api/unigene')


@bp.route('/parse', methods=['POST'])
def parse():
    """Parse multiple UniGene records from file"""
    try:
        file = request.files.get('file')
        max_records = int(request.form.get('max_records', 10))

        if not file:
            return jsonify({'success': False, 'error': 'No file provided'})

        filename = secure_filename(file.filename)
        temp_path = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
        file.save(temp_path)

        records = unigene_parse(temp_path, max_records)

        if os.path.exists(temp_path):
            os.remove(temp_path)

        return jsonify({'success': True, 'records': records, 'count': len(records)})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/read', methods=['POST'])
def read():
    """Read single UniGene record with full details"""
    try:
        file = request.files.get('file')

        if not file:
            return jsonify({'success': False, 'error': 'No file provided'})

        filename = secure_filename(file.filename)
        temp_path = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
        file.save(temp_path)

        record = unigene_read(temp_path)

        if os.path.exists(temp_path):
            os.remove(temp_path)

        return jsonify({'success': True, 'record': record})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})
