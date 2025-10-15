"""
Complete routes for NCBI Entrez database access
"""
from flask import Blueprint, request, jsonify
from utils.entrez_helpers import (
    search_entrez, fetch_summaries, global_query,
    fetch_records, find_related_records, get_database_info,
    format_search_results
)

bp = Blueprint('database', __name__, url_prefix='/api')


@bp.route('/database/search', methods=['POST'])
def search():
    """Search NCBI database with advanced options"""
    try:
        data = request.json
        database = data.get('database', 'pubmed')
        term = data.get('term', '')
        email = data.get('email', 'user@example.com')
        retmax = int(data.get('retmax', 20))
        sort = data.get('sort', 'relevance')
        date_from = data.get('date_from')
        date_to = data.get('date_to')

        # Search
        results = search_entrez(database, term, email, retmax=retmax, sort=sort,
                               date_from=date_from, date_to=date_to)

        if results['IdList']:
            # Fetch summaries
            summaries = fetch_summaries(database, results['IdList'], email)
            formatted = format_search_results(database, summaries)

            return jsonify({
                'success': True,
                'count': int(results['Count']),
                'results': formatted
            })
        else:
            return jsonify({'success': True, 'count': 0, 'results': []})

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/database/global', methods=['POST'])
def global_search():
    """Search across all NCBI databases"""
    try:
        data = request.json
        term = data.get('term', '')
        email = data.get('email', 'user@example.com')

        results = global_query(term, email)

        return jsonify({
            'success': True,
            'results': results
        })

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/database/fetch', methods=['POST'])
def fetch():
    """Fetch full records from NCBI"""
    try:
        data = request.json
        database = data.get('database', 'nucleotide')
        ids_input = data.get('ids', '')

        # Handle both string and list inputs
        if isinstance(ids_input, list):
            ids = [str(i).strip() for i in ids_input if str(i).strip()]
        else:
            ids = [i.strip() for i in ids_input.split(',') if i.strip()]

        email = data.get('email', 'user@example.com')
        format_type = data.get('rettype', data.get('format', 'fasta'))

        if not ids:
            return jsonify({'success': False, 'error': 'No IDs provided'})

        # Map format names
        rettype_map = {
            'fasta': 'fasta',
            'gb': 'gb',
            'gp': 'gp',
            'xml': 'xml',
            'text': 'text'
        }
        rettype = rettype_map.get(format_type, 'fasta')

        records = fetch_records(database, ids, email, rettype=rettype)

        return jsonify({
            'success': True,
            'data': records,
            'format': format_type,
            'count': len(ids)
        })

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/database/link', methods=['POST'])
def link():
    """Find related records using ELink"""
    try:
        data = request.json
        record_id = data.get('id', '')
        from_db = data.get('from_db', 'nucleotide')
        to_db = data.get('to_db', 'protein')
        email = data.get('email', 'user@example.com')

        if not record_id:
            return jsonify({'success': False, 'error': 'No record ID provided'})

        linked_ids = find_related_records(record_id, from_db, to_db, email)

        # If we found links, get summaries
        results = []
        if linked_ids:
            summaries = fetch_summaries(to_db, linked_ids[:50], email)  # Limit to 50
            results = format_search_results(to_db, summaries)

        return jsonify({
            'success': True,
            'count': len(linked_ids),
            'results': results
        })

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/database/info', methods=['POST'])
def info():
    """Get database information using EInfo"""
    try:
        data = request.json
        database = data.get('database', '')
        email = data.get('email', 'user@example.com')

        results = get_database_info(database, email)

        # Format results
        if database:
            # Single database info
            db_info = results['DbInfo']
            formatted = {
                'name': db_info.get('DbName', ''),
                'menu_name': db_info.get('MenuName', ''),
                'description': db_info.get('Description', ''),
                'count': db_info.get('Count', '0'),
                'last_update': db_info.get('LastUpdate', ''),
                'fields': [f['Name'] for f in db_info.get('FieldList', [])[:20]],
                'links': [l['Name'] for l in db_info.get('LinkList', [])[:20]]
            }
            return jsonify({'success': True, 'database': formatted})
        else:
            # List of all databases
            db_list = [{'name': db} for db in results['DbList']]
            return jsonify({'success': True, 'databases': db_list})

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})
