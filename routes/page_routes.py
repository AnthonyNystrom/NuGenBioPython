"""
Routes for rendering HTML pages
"""
from flask import Blueprint, render_template

bp = Blueprint('pages', __name__)


@bp.route('/')
def index():
    return render_template('index.html')


# --- Phase-3 landing hubs -------------------------------------------------
# Sequences and External-data couldn't be fully consolidated (10 and 30 DOM
# id collisions respectively between source pages) without a JS-side
# migration. Instead, these hubs surface the sub-tools as linked cards so
# users still get the sidebar grouping without risking broken JS.

@bp.route('/sequences')
def sequences_hub():
    return render_template('sequences.html')


@bp.route('/data')
def data_hub():
    return render_template('data.html')


@bp.route('/sequence')
def sequence_tools():
    return render_template('sequence.html')


@bp.route('/features')
def features_tools():
    return render_template('features.html')


@bp.route('/seqio')
def seqio_tools():
    return render_template('seqio.html')


@bp.route('/alignment')
def alignment_tools():
    return render_template('compare.html', active_tab='alignment')


@bp.route('/compare')
def compare_hub():
    return render_template('compare.html', active_tab='alignment')


@bp.route('/phylo')
def phylo_tools():
    return render_template('phylogeny.html', active_tab='phylo')


@bp.route('/phylogeny')
def phylogeny_hub():
    return render_template('phylogeny.html', active_tab='phylo')


@bp.route('/structure')
def structure_tools():
    # Structure is already a self-contained hub (all analyses tabbed in one
    # page). No consolidation needed; we just add the URL alias.
    return render_template('structure.html')


@bp.route('/database')
def database_tools():
    return render_template('database.html')


@bp.route('/motifs')
def motifs_tools():
    # Phase-3 hub: /motifs renders the Patterns hub with Motifs tab active.
    # The old motifs.html template is preserved on disk for rollback.
    return render_template('patterns.html', active_tab='motifs')


@bp.route('/restriction')
def restriction_tools():
    return render_template('patterns.html', active_tab='restriction')


@bp.route('/patterns')
def patterns_hub():
    return render_template('patterns.html', active_tab='motifs')


@bp.route('/clustering')
def clustering_tools():
    return render_template('phylogeny.html', active_tab='clustering')


@bp.route('/blast')
def blast_tools():
    return render_template('compare.html', active_tab='blast')


@bp.route('/kegg')
def kegg_tools():
    return render_template('kegg.html')


@bp.route('/genomediagram')
def genomediagram_tools():
    return render_template('genomediagram.html')


@bp.route('/popgen')
def popgen_tools():
    return render_template('phylogeny.html', active_tab='popgen')


@bp.route('/pathway')
def pathway_tools():
    return render_template('pathway.html')


@bp.route('/unigene')
def unigene_tools():
    return render_template('unigene.html')


@bp.route('/hmm')
def hmm_tools():
    return render_template('hmm.html')


@bp.route('/searchio')
def searchio_tools():
    return render_template('compare.html', active_tab='searchio')


@bp.route('/swissprot')
def swissprot_tools():
    return render_template('swissprot.html')


@bp.route('/biodata')
def biodata_tools():
    return render_template('biodata.html')
