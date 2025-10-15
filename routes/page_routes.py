"""
Routes for rendering HTML pages
"""
from flask import Blueprint, render_template

bp = Blueprint('pages', __name__)


@bp.route('/')
def index():
    return render_template('index.html')


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
    return render_template('alignment.html')


@bp.route('/phylo')
def phylo_tools():
    return render_template('phylo.html')


@bp.route('/structure')
def structure_tools():
    return render_template('structure.html')


@bp.route('/database')
def database_tools():
    return render_template('database.html')


@bp.route('/motifs')
def motifs_tools():
    return render_template('motifs.html')


@bp.route('/restriction')
def restriction_tools():
    return render_template('restriction.html')


@bp.route('/clustering')
def clustering_tools():
    return render_template('clustering.html')


@bp.route('/blast')
def blast_tools():
    return render_template('blast.html')


@bp.route('/kegg')
def kegg_tools():
    return render_template('kegg.html')


@bp.route('/genomediagram')
def genomediagram_tools():
    return render_template('genomediagram.html')


@bp.route('/popgen')
def popgen_tools():
    return render_template('popgen.html')


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
    return render_template('searchio.html')


@bp.route('/swissprot')
def swissprot_tools():
    return render_template('swissprot.html')


@bp.route('/biodata')
def biodata_tools():
    return render_template('biodata.html')
