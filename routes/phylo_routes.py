"""
Routes for phylogenetic tree operations
"""
from flask import Blueprint, request, jsonify, current_app, session
import os
from werkzeug.utils import secure_filename

from utils.phylo_helpers import (
    parse_tree_from_string, parse_tree_from_file, tree_to_string,
    tree_to_image_base64, get_tree_info, build_tree_from_alignment,
    calculate_distance_matrix, get_all_terminals, get_all_clades,
    get_path_to_root, get_common_ancestor, prune_tree, collapse_branch,
    ladderize_tree, root_tree, calculate_tree_distance
)

bp = Blueprint('phylo', __name__, url_prefix='/api')


@bp.route('/phylo/parse', methods=['POST'])
def parse_tree():
    try:
        file = request.files.get('file')
        tree_format = request.form.get('format', 'newick')
        show_confidence = request.form.get('show_confidence', 'false').lower() == 'true'

        if file:
            # Handle file upload
            filename = secure_filename(file.filename)
            filepath = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
            file.save(filepath)

            tree = parse_tree_from_file(filepath, tree_format)

            # Store tree in session
            session['current_tree'] = tree_to_string(tree, 'newick')
            session['current_tree_format'] = tree_format

            # Clean up file
            os.remove(filepath)
        else:
            # Handle tree string
            tree_string = request.form.get('tree_string', '').strip()
            if not tree_string:
                return jsonify({'success': False, 'error': 'No tree data provided'})

            tree = parse_tree_from_string(tree_string, tree_format)

            # Store tree in session
            session['current_tree'] = tree_string
            session['current_tree_format'] = tree_format

        # Get tree information
        tree_info = get_tree_info(tree)

        # Create tree visualization
        tree_image = tree_to_image_base64(tree, show_confidence=show_confidence)

        return jsonify({
            'success': True,
            'tree_image': tree_image,
            **tree_info
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/phylo/build', methods=['POST'])
def build_tree():
    """Build a phylogenetic tree from alignment"""
    try:
        file = request.files.get('file')
        alignment_format = request.form.get('alignment_format', 'fasta')
        method = request.form.get('method', 'nj')
        model = request.form.get('model', 'identity')

        if file:
            alignment_string = file.read().decode('utf-8')
        else:
            alignment_string = request.form.get('alignment_string', '').strip()
            if not alignment_string:
                return jsonify({'success': False, 'error': 'No alignment data provided'})

        # Build tree
        tree, distance_matrix = build_tree_from_alignment(
            alignment_string, alignment_format, method, model
        )

        # Store tree in session
        session['current_tree'] = tree_to_string(tree, 'newick')
        session['current_tree_format'] = 'newick'

        # Get tree information
        tree_info = get_tree_info(tree)

        # Create visualization
        tree_image = tree_to_image_base64(tree)

        # Format distance matrix
        dm_names = distance_matrix.names
        dm_matrix = [[distance_matrix[i, j] for j in range(len(dm_names))]
                     for i in range(len(dm_names))]

        return jsonify({
            'success': True,
            'tree_image': tree_image,
            'distance_matrix': {
                'names': dm_names,
                'matrix': dm_matrix
            },
            **tree_info
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/phylo/distance-matrix', methods=['POST'])
def calculate_distance():
    """Calculate distance matrix from alignment"""
    try:
        file = request.files.get('file')
        alignment_format = request.form.get('alignment_format', 'fasta')
        model = request.form.get('model', 'identity')

        if file:
            alignment_string = file.read().decode('utf-8')
        else:
            alignment_string = request.form.get('alignment_string', '').strip()
            if not alignment_string:
                return jsonify({'success': False, 'error': 'No alignment data provided'})

        distance_matrix = calculate_distance_matrix(alignment_string, alignment_format, model)

        dm_names = distance_matrix.names
        dm_matrix = [[distance_matrix[i, j] for j in range(len(dm_names))]
                     for i in range(len(dm_names))]

        return jsonify({
            'success': True,
            'distance_matrix': {
                'names': dm_names,
                'matrix': dm_matrix
            }
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/phylo/terminals', methods=['GET'])
def get_terminals():
    """Get all terminal nodes from current tree"""
    try:
        tree_string = session.get('current_tree')
        tree_format = session.get('current_tree_format', 'newick')

        if not tree_string:
            return jsonify({'success': False, 'error': 'No tree in session'})

        tree = parse_tree_from_string(tree_string, tree_format)
        terminals = get_all_terminals(tree)

        return jsonify({
            'success': True,
            'terminals': terminals
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/phylo/clades', methods=['GET'])
def get_clades():
    """Get all clades from current tree"""
    try:
        tree_string = session.get('current_tree')
        tree_format = session.get('current_tree_format', 'newick')

        if not tree_string:
            return jsonify({'success': False, 'error': 'No tree in session'})

        tree = parse_tree_from_string(tree_string, tree_format)
        clades = get_all_clades(tree)

        return jsonify({
            'success': True,
            'clades': clades
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/phylo/path', methods=['POST'])
def get_path():
    """Get path from terminal to root"""
    try:
        tree_string = session.get('current_tree')
        tree_format = session.get('current_tree_format', 'newick')
        terminal_name = request.json.get('terminal_name')

        if not tree_string:
            return jsonify({'success': False, 'error': 'No tree in session'})

        if not terminal_name:
            return jsonify({'success': False, 'error': 'No terminal name provided'})

        tree = parse_tree_from_string(tree_string, tree_format)
        path = get_path_to_root(tree, terminal_name)

        if path is None:
            return jsonify({'success': False, 'error': f'Terminal {terminal_name} not found'})

        return jsonify({
            'success': True,
            'path': path
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/phylo/common-ancestor', methods=['POST'])
def find_common_ancestor():
    """Find common ancestor of two terminals"""
    try:
        tree_string = session.get('current_tree')
        tree_format = session.get('current_tree_format', 'newick')
        name1 = request.json.get('name1')
        name2 = request.json.get('name2')

        if not tree_string:
            return jsonify({'success': False, 'error': 'No tree in session'})

        if not name1 or not name2:
            return jsonify({'success': False, 'error': 'Two terminal names required'})

        tree = parse_tree_from_string(tree_string, tree_format)
        ancestor = get_common_ancestor(tree, name1, name2)

        if ancestor is None:
            return jsonify({'success': False, 'error': 'Terminals not found'})

        return jsonify({
            'success': True,
            'common_ancestor': ancestor
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/phylo/prune', methods=['POST'])
def prune():
    """Prune a terminal from the tree"""
    try:
        tree_string = session.get('current_tree')
        tree_format = session.get('current_tree_format', 'newick')
        terminal_name = request.json.get('terminal_name')

        if not tree_string:
            return jsonify({'success': False, 'error': 'No tree in session'})

        if not terminal_name:
            return jsonify({'success': False, 'error': 'No terminal name provided'})

        tree = parse_tree_from_string(tree_string, tree_format)
        pruned_tree = prune_tree(tree, terminal_name)

        # Update session
        session['current_tree'] = tree_to_string(pruned_tree, 'newick')

        # Get tree info and visualization
        tree_info = get_tree_info(pruned_tree)
        tree_image = tree_to_image_base64(pruned_tree)

        return jsonify({
            'success': True,
            'tree_image': tree_image,
            **tree_info
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/phylo/collapse', methods=['POST'])
def collapse():
    """Collapse branches in the tree"""
    try:
        tree_string = session.get('current_tree')
        tree_format = session.get('current_tree_format', 'newick')
        threshold = request.json.get('threshold')
        clade_name = request.json.get('clade_name')

        if not tree_string:
            return jsonify({'success': False, 'error': 'No tree in session'})

        tree = parse_tree_from_string(tree_string, tree_format)

        if threshold is not None:
            threshold = float(threshold)

        collapsed_tree = collapse_branch(tree, clade_name, threshold)

        # Update session
        session['current_tree'] = tree_to_string(collapsed_tree, 'newick')

        # Get tree info and visualization
        tree_info = get_tree_info(collapsed_tree)
        tree_image = tree_to_image_base64(collapsed_tree)

        return jsonify({
            'success': True,
            'tree_image': tree_image,
            **tree_info
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/phylo/ladderize', methods=['POST'])
def ladderize():
    """Ladderize (sort) the tree"""
    try:
        tree_string = session.get('current_tree')
        tree_format = session.get('current_tree_format', 'newick')
        reverse = request.json.get('reverse', False)

        if not tree_string:
            return jsonify({'success': False, 'error': 'No tree in session'})

        tree = parse_tree_from_string(tree_string, tree_format)
        ladderized_tree = ladderize_tree(tree, reverse)

        # Update session
        session['current_tree'] = tree_to_string(ladderized_tree, 'newick')

        # Get tree info and visualization
        tree_info = get_tree_info(ladderized_tree)
        tree_image = tree_to_image_base64(ladderized_tree)

        return jsonify({
            'success': True,
            'tree_image': tree_image,
            **tree_info
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/phylo/reroot', methods=['POST'])
def reroot():
    """Re-root the tree"""
    try:
        tree_string = session.get('current_tree')
        tree_format = session.get('current_tree_format', 'newick')
        outgroup_name = request.json.get('outgroup_name')
        at_midpoint = request.json.get('at_midpoint', False)

        if not tree_string:
            return jsonify({'success': False, 'error': 'No tree in session'})

        tree = parse_tree_from_string(tree_string, tree_format)
        rerooted_tree = root_tree(tree, outgroup_name, at_midpoint)

        # Update session
        session['current_tree'] = tree_to_string(rerooted_tree, 'newick')

        # Get tree info and visualization
        tree_info = get_tree_info(rerooted_tree)
        tree_image = tree_to_image_base64(rerooted_tree)

        return jsonify({
            'success': True,
            'tree_image': tree_image,
            **tree_info
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/phylo/compare', methods=['POST'])
def compare_trees():
    """Compare two trees"""
    try:
        tree1_string = session.get('current_tree')
        tree1_format = session.get('current_tree_format', 'newick')

        tree2_string = request.json.get('tree2_string')
        tree2_format = request.json.get('tree2_format', 'newick')

        if not tree1_string:
            return jsonify({'success': False, 'error': 'No tree in session'})

        if not tree2_string:
            return jsonify({'success': False, 'error': 'No second tree provided'})

        tree1 = parse_tree_from_string(tree1_string, tree1_format)
        tree2 = parse_tree_from_string(tree2_string, tree2_format)

        distance_info = calculate_tree_distance(tree1, tree2)

        return jsonify({
            'success': True,
            **distance_info
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/phylo/export', methods=['POST'])
def export_tree():
    """Export tree in different format"""
    try:
        tree_string = session.get('current_tree')
        tree_format = session.get('current_tree_format', 'newick')
        export_format = request.json.get('export_format', 'newick')

        if not tree_string:
            return jsonify({'success': False, 'error': 'No tree in session'})

        tree = parse_tree_from_string(tree_string, tree_format)
        exported_string = tree_to_string(tree, export_format)

        return jsonify({
            'success': True,
            'tree_string': exported_string,
            'format': export_format
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/phylo/visualize', methods=['POST'])
def visualize():
    """Visualize tree with custom options"""
    try:
        tree_string = session.get('current_tree')
        tree_format = session.get('current_tree_format', 'newick')
        show_confidence = request.json.get('show_confidence', False)

        if not tree_string:
            return jsonify({'success': False, 'error': 'No tree in session'})

        tree = parse_tree_from_string(tree_string, tree_format)
        tree_image = tree_to_image_base64(tree, show_confidence=show_confidence)

        return jsonify({
            'success': True,
            'tree_image': tree_image
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})
