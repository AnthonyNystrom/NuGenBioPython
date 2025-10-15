"""
Complete routes for KEGG database operations using Bio.KEGG.REST
"""
from flask import Blueprint, request, jsonify
from utils.kegg_helpers import (
    kegg_find_search, kegg_list_entries, kegg_link_entries,
    kegg_convert_ids, kegg_get_info, kegg_get_entry, parse_kegg_entry
)

bp = Blueprint('kegg', __name__, url_prefix='/api')


@bp.route('/kegg/search', methods=['POST'])
def search():
    """Search KEGG database using kegg_find"""
    try:
        data = request.json
        database = data.get('database', 'pathway')
        query = data.get('query', '')
        organism = data.get('organism', None)

        if not query:
            return jsonify({'success': False, 'error': 'Search query is required'})

        results = kegg_find_search(database, query, organism)

        return jsonify({
            'success': True,
            'results': results[:50],  # Limit to 50 results
            'count': len(results)
        })

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/kegg/list', methods=['POST'])
def list_database():
    """List entries in KEGG database using kegg_list"""
    try:
        data = request.json
        database = data.get('database', 'pathway')
        organism = data.get('organism', None)
        limit = data.get('limit', 100)

        results = kegg_list_entries(database, organism)

        return jsonify({
            'success': True,
            'results': results[:limit],
            'total': len(results),
            'displayed': min(limit, len(results))
        })

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/kegg/link', methods=['POST'])
def link():
    """Find related entries using kegg_link"""
    try:
        data = request.json
        target_db = data.get('target_db', 'pathway')
        source_db = data.get('source_db', 'genes')
        source_id = data.get('source_id', None)

        if not target_db or not source_db:
            return jsonify({'success': False, 'error': 'Target and source databases are required'})

        results = kegg_link_entries(target_db, source_db, source_id)

        return jsonify({
            'success': True,
            'results': results[:100],  # Limit to 100 links
            'count': len(results)
        })

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/kegg/convert', methods=['POST'])
def convert():
    """Convert identifiers using kegg_conv"""
    try:
        data = request.json
        target_db = data.get('target_db', 'ncbi-geneid')
        source_db = data.get('source_db', 'hsa')
        ids = data.get('ids', None)

        if not target_db or not source_db:
            return jsonify({'success': False, 'error': 'Target and source databases are required'})

        # Parse IDs if provided as string
        if ids and isinstance(ids, str):
            ids = [i.strip() for i in ids.split(',') if i.strip()]

        results = kegg_convert_ids(target_db, source_db, ids)

        return jsonify({
            'success': True,
            'results': results[:100],  # Limit to 100 conversions
            'count': len(results)
        })

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/kegg/info', methods=['POST'])
def info():
    """Get database information using kegg_info"""
    try:
        data = request.json
        database = data.get('database', None)

        info_data = kegg_get_info(database)

        return jsonify({
            'success': True,
            'database': database or 'kegg',
            'info': info_data
        })

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/kegg/get/<entry_id>', methods=['GET'])
def get_entry(entry_id):
    """Get detailed entry information using kegg_get"""
    try:
        entry_data = kegg_get_entry(entry_id)
        parsed_data = parse_kegg_entry(entry_data)

        # Get pathway image URL if it's a pathway
        image_url = None
        if 'path:' in entry_id or entry_id.startswith('hsa') or entry_id.startswith('map'):
            pathway_code = entry_id.lower().replace('path:', '')
            if pathway_code.startswith('map'):
                image_url = f"https://www.kegg.jp/kegg/pathway/map/{pathway_code}.png"
            else:
                image_url = f"https://www.kegg.jp/kegg/pathway/{pathway_code}.png"

        return jsonify({
            'success': True,
            'entry_id': entry_id,
            'raw_data': entry_data,
            'parsed_data': parsed_data,
            'image_url': image_url
        })

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})
