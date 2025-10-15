"""
Routes for Bio.SwissProt operations
Handles parsing and reading SwissProt/UniProt protein database records
"""
from flask import Blueprint, request, jsonify, current_app
import os
from werkzeug.utils import secure_filename
from utils.swissprot_helpers import swissprot_parse, swissprot_read

bp = Blueprint('swissprot', __name__, url_prefix='/api/swissprot')


@bp.route('/parse', methods=['POST'])
def parse():
    """Parse multiple SwissProt records from file"""
    try:
        file = request.files.get('file')
        max_records = int(request.form.get('max_records', 10))

        if not file:
            return jsonify({'success': False, 'error': 'No file provided'})

        # Save uploaded file
        filename = secure_filename(file.filename)
        temp_path = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
        file.save(temp_path)

        # Parse records
        records = swissprot_parse(temp_path, max_records)

        # Clean up
        if os.path.exists(temp_path):
            os.remove(temp_path)

        return jsonify({
            'success': True,
            'records': records,
            'count': len(records)
        })

    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        })


@bp.route('/read', methods=['POST'])
def read():
    """Read single SwissProt record with full details"""
    try:
        file = request.files.get('file')

        if not file:
            return jsonify({'success': False, 'error': 'No file provided'})

        # Save uploaded file
        filename = secure_filename(file.filename)
        temp_path = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
        file.save(temp_path)

        # Read single record
        record = swissprot_read(temp_path)

        # Clean up
        if os.path.exists(temp_path):
            os.remove(temp_path)

        return jsonify({
            'success': True,
            'record': record
        })

    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        })
