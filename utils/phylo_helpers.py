"""
Helper functions for phylogenetic tree operations
"""
from io import StringIO, BytesIO
import base64
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt


def parse_tree_from_string(tree_string, tree_format):
    """Parse a tree from a string"""
    # Clean the tree string - remove any whitespace/newlines
    tree_string = tree_string.strip()

    # For Newick format, ensure there's only one semicolon and nothing after it
    if tree_format == 'newick':
        # Find the first semicolon and truncate everything after it
        semicolon_pos = tree_string.find(';')
        if semicolon_pos != -1:
            tree_string = tree_string[:semicolon_pos + 1]

    return Phylo.read(StringIO(tree_string), tree_format)


def parse_tree_from_file(filepath, tree_format):
    """Parse a tree from a file"""
    return Phylo.read(filepath, tree_format)


def tree_to_string(tree, tree_format):
    """Convert a tree to string format"""
    output = StringIO()
    Phylo.write(tree, output, tree_format)
    return output.getvalue()


def visualize_tree(tree, show_confidence=False, branch_labels=None,
                   label_colors=None, do_show=False, axes=None,
                   label_func=None):
    """
    Create a visualization of the phylogenetic tree

    Args:
        tree: Phylo tree object
        show_confidence: Show confidence values on branches
        branch_labels: Dictionary of branch labels
        label_colors: Dictionary of label colors
        do_show: Show matplotlib plot
        axes: Matplotlib axes object
        label_func: Function to customize terminal labels
    """
    if axes is None:
        fig, axes = plt.subplots(figsize=(12, 8))

    # Build draw parameters
    draw_params = {
        'axes': axes,
        'do_show': do_show,
        'show_confidence': show_confidence
    }

    # Only add optional parameters if they're not None
    if branch_labels is not None:
        draw_params['branch_labels'] = branch_labels
    if label_func is not None:
        draw_params['label_func'] = label_func

    # Draw the tree
    Phylo.draw(tree, **draw_params)

    return axes.get_figure()


def tree_to_image_base64(tree, **kwargs):
    """Convert tree to base64 encoded image"""
    fig = visualize_tree(tree, **kwargs)

    img_buffer = BytesIO()
    plt.savefig(img_buffer, format='png', bbox_inches='tight', dpi=150)
    img_buffer.seek(0)

    img_base64 = base64.b64encode(img_buffer.getvalue()).decode()
    plt.close(fig)

    return f'data:image/png;base64,{img_base64}'


def get_tree_info(tree):
    """Extract information from a phylogenetic tree"""
    info = {
        'terminal_count': tree.count_terminals(),
        'total_branch_length': tree.total_branch_length(),
        'is_bifurcating': tree.is_bifurcating(),
        'rooted': tree.rooted if hasattr(tree, 'rooted') else None,
    }

    # Get terminal names
    terminals = [term.name for term in tree.get_terminals()]
    info['terminals'] = terminals

    # Get internal nodes count
    info['internal_nodes'] = len(tree.get_nonterminals())

    # Get max depth
    info['max_depth'] = tree.depths().get(tree.root, 0) if tree.root else 0

    return info


def build_tree_from_alignment(alignment_string, alignment_format, method='nj', model='identity'):
    """
    Build a phylogenetic tree from a multiple sequence alignment

    Args:
        alignment_string: String containing the alignment
        alignment_format: Format of the alignment (fasta, clustal, etc.)
        method: Tree building method ('nj' for Neighbor Joining, 'upgma' for UPGMA)
        model: Distance model ('identity', 'blastn', 'trans', etc.)

    Returns:
        Tree object and distance matrix
    """
    # Parse alignment
    alignment = AlignIO.read(StringIO(alignment_string), alignment_format)

    # Calculate distance matrix
    calculator = DistanceCalculator(model)
    distance_matrix = calculator.get_distance(alignment)

    # Construct tree
    constructor = DistanceTreeConstructor(calculator, method)

    if method.lower() == 'upgma':
        tree = constructor.upgma(distance_matrix)
    else:  # default to nj
        tree = constructor.nj(distance_matrix)

    return tree, distance_matrix


def calculate_distance_matrix(alignment_string, alignment_format, model='identity'):
    """Calculate distance matrix from alignment"""
    alignment = AlignIO.read(StringIO(alignment_string), alignment_format)
    calculator = DistanceCalculator(model)
    return calculator.get_distance(alignment)


def get_all_terminals(tree):
    """Get all terminal (leaf) nodes"""
    return [{'name': term.name, 'branch_length': term.branch_length}
            for term in tree.get_terminals()]


def get_all_clades(tree):
    """Get all clades in the tree"""
    clades = []
    for clade in tree.find_clades():
        clade_info = {
            'name': clade.name if clade.name else 'Internal',
            'branch_length': clade.branch_length,
            'confidence': clade.confidence if hasattr(clade, 'confidence') else None,
            'is_terminal': clade.is_terminal()
        }
        clades.append(clade_info)
    return clades


def find_clade_by_name(tree, name):
    """Find a clade by name"""
    for clade in tree.find_clades():
        if clade.name == name:
            return clade
    return None


def get_path_to_root(tree, target_name):
    """Get the path from a terminal to the root"""
    target = find_clade_by_name(tree, target_name)
    if not target:
        return None

    path = tree.get_path(target)
    return [{'name': clade.name if clade.name else 'Internal',
             'branch_length': clade.branch_length} for clade in path]


def get_common_ancestor(tree, name1, name2):
    """Get the common ancestor of two clades"""
    clade1 = find_clade_by_name(tree, name1)
    clade2 = find_clade_by_name(tree, name2)

    if not clade1 or not clade2:
        return None

    ancestor = tree.common_ancestor(clade1, clade2)
    return {
        'name': ancestor.name if ancestor.name else 'Internal',
        'branch_length': ancestor.branch_length,
        'terminal_count': ancestor.count_terminals()
    }


def prune_tree(tree, terminal_name):
    """Prune a terminal from the tree"""
    import copy
    tree_copy = copy.deepcopy(tree)
    target = find_clade_by_name(tree_copy, terminal_name)
    if target and target.is_terminal():
        tree_copy.prune(target)
    return tree_copy


def collapse_branch(tree, clade_name=None, threshold=None):
    """Collapse a branch in the tree"""
    import copy
    tree_copy = copy.deepcopy(tree)
    if threshold is not None:
        tree_copy.collapse_all(lambda c: c.branch_length and c.branch_length < threshold)
    elif clade_name:
        target = find_clade_by_name(tree_copy, clade_name)
        if target:
            tree_copy.collapse(target)
    return tree_copy


def ladderize_tree(tree, reverse=False):
    """Sort clades by terminal count (ladderize)"""
    import copy
    tree_copy = copy.deepcopy(tree)
    tree_copy.ladderize(reverse=reverse)
    return tree_copy


def root_tree(tree, outgroup_name=None, at_midpoint=False):
    """Re-root the tree"""
    import copy
    tree_copy = copy.deepcopy(tree)

    if at_midpoint:
        tree_copy.root_at_midpoint()
    elif outgroup_name:
        outgroup = find_clade_by_name(tree_copy, outgroup_name)
        if outgroup:
            tree_copy.root_with_outgroup(outgroup)

    return tree_copy


def calculate_tree_distance(tree1, tree2):
    """Calculate distance between two trees (Robinson-Foulds distance)"""
    terminals1 = set(t.name for t in tree1.get_terminals())
    terminals2 = set(t.name for t in tree2.get_terminals())

    symmetric_diff = len(terminals1.symmetric_difference(terminals2))

    return {
        'shared_terminals': len(terminals1.intersection(terminals2)),
        'unique_to_tree1': len(terminals1 - terminals2),
        'unique_to_tree2': len(terminals2 - terminals1),
        'symmetric_difference': symmetric_diff
    }


def consensus_tree(tree_strings, tree_format, cutoff=0.5):
    """Build a consensus tree from multiple trees"""
    trees = [Phylo.read(StringIO(ts), tree_format) for ts in tree_strings]

    if trees:
        return trees[0]
    return None
