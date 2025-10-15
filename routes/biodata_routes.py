"""
Routes for Bio.Data operations (CodonTable, IUPACData, PDBData)
"""
from flask import Blueprint, request, jsonify
from utils.biodata_helpers import (
    get_codon_tables, translate_sequence, get_iupac_codes,
    convert_protein_letters, calculate_molecular_weight,
    get_pdb_conversions, get_atom_weights, lookup_iupac_code
)

bp = Blueprint('biodata', __name__, url_prefix='/api/biodata')


# CodonTable Routes
@bp.route('/codon_tables', methods=['GET'])
def codon_tables():
    """Get all available codon tables"""
    try:
        tables = get_codon_tables()
        return jsonify({'success': True, 'tables': tables})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/translate', methods=['POST'])
def translate():
    """Translate DNA sequence using specified codon table"""
    try:
        data = request.get_json()
        sequence = data.get('sequence', '').upper().strip()
        table_id = int(data.get('table_id', 1))

        if not sequence:
            return jsonify({'success': False, 'error': 'No sequence provided'})

        result = translate_sequence(sequence, table_id)
        return jsonify(result)
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


# IUPACData Routes
@bp.route('/iupac_codes', methods=['GET'])
def iupac_codes():
    """Get IUPAC ambiguity codes"""
    try:
        code_type = request.args.get('type', 'dna')
        codes = get_iupac_codes(code_type)
        return jsonify({'success': True, 'codes': codes})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/iupac_lookup', methods=['POST'])
def iupac_lookup():
    """Look up what a specific IUPAC code represents"""
    try:
        data = request.get_json()
        code = data.get('code', '')
        code_type = data.get('type', 'dna')

        if not code:
            return jsonify({'success': False, 'error': 'No code provided'})

        result = lookup_iupac_code(code, code_type)
        return jsonify(result)
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/convert_protein', methods=['POST'])
def convert_protein():
    """Convert protein letters between 1-letter and 3-letter codes"""
    try:
        data = request.get_json()
        input_str = data.get('input', '')
        conversion_type = data.get('conversion_type', '1to3')

        if not input_str:
            return jsonify({'success': False, 'error': 'No input provided'})

        result = convert_protein_letters(input_str, conversion_type)
        return jsonify(result)
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/molecular_weight', methods=['POST'])
def molecular_weight():
    """Calculate molecular weight of a sequence"""
    try:
        data = request.get_json()
        sequence = data.get('sequence', '')
        seq_type = data.get('seq_type', 'protein')
        weight_type = data.get('weight_type', 'average')

        if not sequence:
            return jsonify({'success': False, 'error': 'No sequence provided'})

        result = calculate_molecular_weight(sequence, seq_type, weight_type)
        return jsonify(result)
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/atom_weights', methods=['GET'])
def atom_weights():
    """Get atomic weights"""
    try:
        weights = get_atom_weights()
        return jsonify({'success': True, 'weights': weights})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


# PDBData Routes
@bp.route('/pdb_conversions', methods=['GET'])
def pdb_conversions():
    """Get PDB residue name conversions"""
    try:
        conversion_type = request.args.get('type', 'protein_1to3')
        conversions = get_pdb_conversions(conversion_type)
        return jsonify({'success': True, 'conversions': conversions})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})
