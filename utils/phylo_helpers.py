"""
Helper functions for phylogenetic tree operations
"""
from io import StringIO, BytesIO
import base64
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

from utils.plot_helpers import (
    figure_to_svg_data_url, svg_markup, style_axes, set_title,
    TITLE_COLOR, LABEL_COLOR, MUTED_COLOR, AXIS_COLOR,
    subplots as oo_subplots,
)


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


def _tree_figure_size(tree):
    """Size the figure based on terminal count and longest label so tip
    labels don't clip. Width grows with label length; height with taxa."""
    terminals = tree.get_terminals()
    n = max(1, len(terminals))
    max_label = max((len(t.name or '') for t in terminals), default=8)
    # Height: enough vertical space for each tip (~0.25 in per tip), capped.
    height = max(4.5, min(28, 2.4 + n * 0.28))
    # Width: base 9 in, plus ~0.07 in per char of the longest label.
    width = max(9, min(22, 9 + max_label * 0.07))
    return width, height


def visualize_tree(tree, show_confidence=False, branch_labels=None,
                   label_colors=None, do_show=False, axes=None,
                   label_func=None, figsize=None):
    """Draw a phylogenetic tree with app-consistent styling.

    Unlike the stock Phylo.draw wrapper, this auto-sizes the figure for the
    number of taxa and the longest tip label, so wide trees and dense trees
    both render without clipping.
    """
    owns_figure = axes is None
    if axes is None:
        w, h = figsize or _tree_figure_size(tree)
        fig, axes = oo_subplots(figsize=(w, h), dpi=100)

    draw_params = {
        'axes': axes,
        'do_show': do_show,
        'show_confidence': show_confidence,
    }
    if branch_labels is not None:
        draw_params['branch_labels'] = branch_labels
    if label_func is not None:
        draw_params['label_func'] = label_func

    Phylo.draw(tree, **draw_params)

    # Clean up BioPython's default styling — it sets a heavy title and
    # thick axes that don't match the app.
    axes.set_title('')
    axes.set_xlabel('Branch length', fontsize=9, color=LABEL_COLOR)
    axes.set_ylabel('')
    # BioPython's Phylo.draw leaves integer y-axis ticks (1, 2, 3 ...)
    # for clade indices — not useful to a reader. Hide them.
    axes.set_yticks([])
    style_axes(axes, hide=('top', 'right', 'left'))
    axes.margins(x=0.02, y=0.02)

    # Recolor BioPython's default black tip labels and branch labels to
    # match the rest of the app. This walks the Text artists in the axes.
    for txt in axes.texts:
        # Tip labels sit at the right edge; apply a consistent size/color.
        txt.set_color(LABEL_COLOR)
        if txt.get_fontsize() and txt.get_fontsize() > 12:
            txt.set_fontsize(10)

    return axes.get_figure()


def tree_to_inline_svg(tree, **kwargs):
    """Render a tree as SVG markup for injection into the page.

    Unlike tree_to_image_base64(), the result is live DOM: the taxon labels
    are selectable text and CSS can restyle the figure.
    """
    return svg_markup(visualize_tree(tree, **kwargs),
                      title='Phylogenetic tree')


def tree_to_image_base64(tree, **kwargs):
    """Convert tree to SVG data URL (vector — scales cleanly in browsers)."""
    fig = visualize_tree(tree, **kwargs)
    return figure_to_svg_data_url(fig)


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
