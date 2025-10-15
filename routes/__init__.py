"""
Blueprint registration for all route modules
"""
from routes import page_routes
from routes import sequence_routes
from routes import seqio_routes
from routes import features_routes
from routes import alignment_routes
from routes import phylo_routes
from routes import structure_routes
from routes import database_routes
from routes import motif_routes
from routes import restriction_routes
from routes import clustering_routes
from routes import blast_routes
from routes import kegg_routes
from routes import genomediagram_routes
from routes import searchio_routes
from routes import swissprot_routes
from routes import unigene_routes
from routes import biodata_routes
from routes import advanced_routes
from routes import specialty_routes
from routes import pathway_routes


def register_blueprints(app):
    """Register all blueprints with the Flask application"""
    # Register page blueprint first
    app.register_blueprint(page_routes.bp)

    # Create endpoint aliases for backward compatibility with templates
    # Templates use url_for('index') but blueprint creates 'pages.index'
    # So we create aliases to make both work
    endpoint_names = ['index', 'sequence_tools', 'features_tools', 'seqio_tools', 'alignment_tools', 'phylo_tools',
                      'structure_tools', 'database_tools', 'motifs_tools', 'restriction_tools',
                      'clustering_tools', 'blast_tools', 'kegg_tools', 'genomediagram_tools',
                      'popgen_tools', 'pathway_tools', 'unigene_tools', 'hmm_tools',
                      'searchio_tools', 'swissprot_tools', 'biodata_tools']

    for endpoint_name in endpoint_names:
        blueprint_endpoint = f'pages.{endpoint_name}'
        if blueprint_endpoint in app.view_functions:
            # Add URL rule with the original endpoint name
            # Find the URL rule for this blueprint endpoint
            for rule in app.url_map.iter_rules():
                if rule.endpoint == blueprint_endpoint:
                    app.add_url_rule(rule.rule, endpoint_name, app.view_functions[blueprint_endpoint])
                    break

    # API routes (all have /api prefix)
    app.register_blueprint(sequence_routes.bp)
    app.register_blueprint(seqio_routes.bp)
    app.register_blueprint(features_routes.bp)
    app.register_blueprint(alignment_routes.bp)
    app.register_blueprint(phylo_routes.bp)
    app.register_blueprint(structure_routes.bp)
    app.register_blueprint(database_routes.bp)
    app.register_blueprint(motif_routes.bp)
    app.register_blueprint(restriction_routes.bp)
    app.register_blueprint(clustering_routes.bp)
    app.register_blueprint(blast_routes.bp)
    app.register_blueprint(kegg_routes.bp)
    app.register_blueprint(genomediagram_routes.bp)
    app.register_blueprint(searchio_routes.bp)
    app.register_blueprint(swissprot_routes.bp)
    app.register_blueprint(unigene_routes.bp)
    app.register_blueprint(biodata_routes.bp)
    app.register_blueprint(advanced_routes.bp)
    app.register_blueprint(specialty_routes.bp)
    app.register_blueprint(pathway_routes.bp)
