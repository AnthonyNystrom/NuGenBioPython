"""
Routes for pathway analysis - Complete Bio.Pathway implementation
Supports: Reaction, System, Network, MultiGraph
"""
from flask import Blueprint, request, jsonify

from utils.request_helpers import error_response
import base64
import math
from io import BytesIO
import json

import matplotlib.patches as mpatches

from dependencies import Pathway, plt, np
from utils.plot_helpers import (
    ROLE_COLORS, EDGE_COLOR, TITLE_COLOR, LABEL_COLOR,
    figure_to_svg_data_url, set_title,
)

bp = Blueprint('pathway', __name__, url_prefix='/api')


# ============================================================================
# SYSTEM ANALYSIS - Bio.Pathway.System
# ============================================================================

@bp.route('/pathway/analyze_system', methods=['POST'])
def analyze_system():
    """Analyze complete reaction system using Bio.Pathway.System"""
    try:
        from Bio.Pathway import System, Reaction

        data = request.get_json(silent=True) or {}
        reactions_data = data.get('reactions', [])

        if not reactions_data:
            return jsonify({'success': False, 'error': 'No reactions provided'})

        # Create system
        system = System()
        reversible_reactions = []
        irreversible_reactions = []

        # Add reactions to system
        for rxn_data in reactions_data:
            # Create reaction with proper stoichiometry
            species_dict = rxn_data['species']  # Already combined reactants + products
            catalysts = tuple(rxn_data['catalysts']) if rxn_data['catalysts'] else ()
            reversible = 1 if rxn_data['reversible'] else 0

            reaction = Reaction(
                reactants=species_dict,
                catalysts=catalysts,
                reversible=reversible,
                data={'name': rxn_data['name']}
            )

            system.add_reaction(reaction)

            if rxn_data['reversible']:
                reversible_reactions.append(rxn_data['name'])
            else:
                irreversible_reactions.append(rxn_data['name'])

        # Get system information
        all_species = list(system.species())
        all_reactions = list(system.reactions())

        return jsonify({
            'success': True,
            'analysis': {
                'reaction_count': len(all_reactions),
                'species_count': len(all_species),
                'species': all_species,
                'reversible_count': len(reversible_reactions),
                'irreversible_count': len(irreversible_reactions),
                'reversible_reactions': reversible_reactions,
                'irreversible_reactions': irreversible_reactions
            }
        })

    except Exception as e:
        return error_response(e, context='pathway_routes.analyze_system')


# ============================================================================
# NETWORK ANALYSIS - Bio.Pathway.Network
# ============================================================================

@bp.route('/pathway/analyze_network', methods=['POST'])
def analyze_network():
    """Analyze pathway network using Bio.Pathway.Network"""
    try:
        from Bio.Pathway import Network

        data = request.get_json(silent=True) or {}
        reactions_data = data.get('reactions', [])

        if not reactions_data:
            return jsonify({'success': False, 'error': 'No reactions provided'})

        # Create network
        network = Network()
        all_species = set()

        # Add species
        for rxn_data in reactions_data:
            for species in rxn_data['reactants'].keys():
                all_species.add(species)
            for species in rxn_data['products'].keys():
                all_species.add(species)

        for species in all_species:
            network.add_species(species)

        # Add interactions
        interactions_list = []
        for rxn_data in reactions_data:
            for reactant in rxn_data['reactants'].keys():
                for product in rxn_data['products'].keys():
                    network.add_interaction(reactant, product, rxn_data['name'])
                    interactions_list.append({
                        'source': reactant,
                        'sink': product,
                        'reaction': rxn_data['name']
                    })

        # Analyze sources and sinks
        sources = []
        sinks = []
        intermediates = []

        for species in all_species:
            upstream = network.source(species)  # Species that produce this species
            downstream = network.sink(species)  # Species that this species produces

            if len(upstream) == 0 and len(downstream) > 0:
                sources.append(species)
            elif len(downstream) == 0 and len(upstream) > 0:
                sinks.append(species)
            else:
                intermediates.append(species)

        # Get species connections
        species_connections = {}
        for species in all_species:
            species_connections[species] = {
                'upstream': list(network.source(species)),
                'downstream': list(network.sink(species))
            }

        return jsonify({
            'success': True,
            'analysis': {
                'sources': sources,
                'sinks': sinks,
                'intermediates': intermediates,
                'interactions': interactions_list,
                'species_connections': species_connections,
                'total_species': len(all_species),
                'total_interactions': len(interactions_list)
            }
        })

    except Exception as e:
        return error_response(e, context='pathway_routes.analyze_network')


# ============================================================================
# MULTIGRAPH VISUALIZATION - Bio.Pathway.MultiGraph
# ============================================================================

@bp.route('/pathway/visualize', methods=['POST'])
def visualize_pathway():
    """Generate pathway visualization using MultiGraph"""
    try:
        data = request.get_json(silent=True) or {}
        reactions_data = data.get('reactions', [])

        if not reactions_data:
            return jsonify({'success': False, 'error': 'No reactions provided'})

        import networkx as nx

        G = nx.DiGraph()
        for rxn_data in reactions_data:
            for reactant in rxn_data['reactants'].keys():
                G.add_node(reactant)
            for product in rxn_data['products'].keys():
                G.add_node(product)
            for reactant in rxn_data['reactants'].keys():
                for product in rxn_data['products'].keys():
                    G.add_edge(reactant, product,
                               reaction=rxn_data['name'],
                               reversible=rxn_data.get('reversible', False))

        if G.number_of_nodes() == 0:
            return jsonify({'success': False,
                            'error': 'No species extracted from reactions'})

        # Classify nodes
        node_roles = {}
        for node in G.nodes():
            in_deg = G.in_degree(node)
            out_deg = G.out_degree(node)
            if in_deg == 0 and out_deg > 0:
                node_roles[node] = 'source'
            elif out_deg == 0 and in_deg > 0:
                node_roles[node] = 'sink'
            else:
                node_roles[node] = 'intermediate'

        # Pick a layout that matches graph size. Kamada-Kawai gives cleaner
        # results for small-to-medium graphs; spring scales better past ~40
        # nodes. Graphviz "dot" would be ideal for flow diagrams but needs a
        # native binary, so we avoid that dependency.
        n = G.number_of_nodes()
        try:
            if n <= 40:
                pos = nx.kamada_kawai_layout(G)
            else:
                pos = nx.spring_layout(G, k=1.5 / math.sqrt(n),
                                       iterations=80, seed=42)
        except Exception:
            pos = nx.spring_layout(G, seed=42)

        # Figure size scales with node count; wider for more nodes, taller
        # for tall paths.
        fig_w = max(9, min(20, 6 + math.sqrt(n) * 1.5))
        fig_h = max(6, min(14, 4 + math.sqrt(n) * 1.1))
        fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=100)

        # Node sizes grow with label length so text fits.
        max_label = max((len(str(node)) for node in G.nodes()), default=8)
        node_size = max(1400, min(4200, 380 * max_label))

        # Draw edges first so nodes sit on top. Bent edges help readability
        # on dense graphs.
        node_colors = [ROLE_COLORS[node_roles[n]] for n in G.nodes()]
        # Curve bidirectional pairs in opposite directions to avoid overlap.
        bidirectional = set()
        for u, v in G.edges():
            if G.has_edge(v, u) and (v, u) not in bidirectional:
                bidirectional.add((u, v))

        for u, v, attrs in G.edges(data=True):
            rad = 0.15 if (u, v) in bidirectional or (v, u) in bidirectional else 0.05
            nx.draw_networkx_edges(
                G, pos, edgelist=[(u, v)], ax=ax,
                edge_color='#94a3b8',
                arrows=True, arrowsize=14, arrowstyle='-|>',
                width=1.3, alpha=0.75,
                connectionstyle=f'arc3,rad={rad}',
                node_size=node_size,
            )

        nx.draw_networkx_nodes(
            G, pos, ax=ax,
            node_color=node_colors,
            node_size=node_size,
            edgecolors=EDGE_COLOR, linewidths=0.8, alpha=0.95,
        )
        nx.draw_networkx_labels(
            G, pos, ax=ax,
            font_size=9, font_color='white', font_weight='600',
        )

        # Edge labels: only when graph is small enough for them to be legible.
        if G.number_of_edges() <= 20:
            edge_labels = {(u, v): d.get('reaction', '') for u, v, d in G.edges(data=True)}
            nx.draw_networkx_edge_labels(
                G, pos, ax=ax, edge_labels=edge_labels,
                font_size=7, font_color=LABEL_COLOR,
                bbox=dict(boxstyle='round,pad=0.15',
                          fc='white', ec='none', alpha=0.85),
            )

        # Legend
        legend_items = []
        for role, color in [('source', ROLE_COLORS['source']),
                            ('intermediate', ROLE_COLORS['intermediate']),
                            ('sink', ROLE_COLORS['sink'])]:
            if any(r == role for r in node_roles.values()):
                legend_items.append(
                    mpatches.Patch(facecolor=color, edgecolor=EDGE_COLOR,
                                   linewidth=0.5, label=role.title()))
        if legend_items:
            ax.legend(handles=legend_items, loc='upper left',
                      frameon=False, fontsize=9, labelcolor=LABEL_COLOR)

        set_title(ax, 'Pathway Network', pad=14)
        # Extra margin so node labels at the extremes don't clip.
        ax.margins(0.18)
        ax.set_axis_off()

        graph_image = figure_to_svg_data_url(fig)

        return jsonify({
            'success': True,
            'visualization': {
                'graph_image': graph_image,
                'node_count': G.number_of_nodes(),
                'edge_count': G.number_of_edges(),
                'graph_type': 'Directed Graph',
                'sources': [n for n, r in node_roles.items() if r == 'source'],
                'sinks': [n for n, r in node_roles.items() if r == 'sink'],
            }
        })

    except Exception as e:
        return error_response(e, context='pathway_routes.visualize_pathway')


# ============================================================================
# EXPORT FORMATS
# ============================================================================

@bp.route('/pathway/export', methods=['POST'])
def export_pathway():
    """Export pathway in various formats"""
    try:
        data = request.get_json(silent=True) or {}
        reactions_data = data.get('reactions', [])
        export_format = data.get('format', 'json')

        if not reactions_data:
            return jsonify({'success': False, 'error': 'No reactions provided'})

        if export_format == 'json':
            # Export as JSON
            content = json.dumps(reactions_data, indent=2)

        elif export_format == 'txt':
            # Export as human-readable text
            lines = ['PATHWAY ANALYSIS EXPORT', '=' * 50, '']
            for i, rxn in enumerate(reactions_data, 1):
                lines.append(f'Reaction {i}: {rxn["name"]}')

                reactants_str = ' + '.join([f'{abs(c)} {s}' for s, c in rxn['reactants'].items()])
                products_str = ' + '.join([f'{c} {s}' for s, c in rxn['products'].items()])
                arrow = ' <=> ' if rxn['reversible'] else ' => '

                lines.append(f'  {reactants_str}{arrow}{products_str}')

                if rxn['catalysts']:
                    lines.append(f'  Catalysts: {", ".join(rxn["catalysts"])}')

                lines.append('')

            content = '\n'.join(lines)

        elif export_format == 'dot':
            # Export as GraphViz DOT format
            lines = ['digraph pathway {']
            lines.append('  rankdir=LR;')
            lines.append('  node [shape=box, style=filled];')
            lines.append('')

            # Add nodes
            all_species = set()
            for rxn in reactions_data:
                all_species.update(rxn['reactants'].keys())
                all_species.update(rxn['products'].keys())

            for species in all_species:
                lines.append(f'  "{species}" [fillcolor=lightblue];')

            lines.append('')

            # Add edges
            for rxn in reactions_data:
                for reactant in rxn['reactants'].keys():
                    for product in rxn['products'].keys():
                        label = rxn['name']
                        if rxn['catalysts']:
                            label += f'\\n({", ".join(rxn["catalysts"])})'
                        lines.append(f'  "{reactant}" -> "{product}" [label="{label}"];')

            lines.append('}')
            content = '\n'.join(lines)

        elif export_format == 'graphml':
            # Export as GraphML format
            lines = ['<?xml version="1.0" encoding="UTF-8"?>']
            lines.append('<graphml xmlns="http://graphml.graphdrawing.org/xmlns">')
            lines.append('  <key id="name" for="node" attr.name="name" attr.type="string"/>')
            lines.append('  <key id="reaction" for="edge" attr.name="reaction" attr.type="string"/>')
            lines.append('  <graph id="pathway" edgedefault="directed">')

            # Add nodes
            all_species = set()
            for rxn in reactions_data:
                all_species.update(rxn['reactants'].keys())
                all_species.update(rxn['products'].keys())

            for species in all_species:
                lines.append(f'    <node id="{species}">')
                lines.append(f'      <data key="name">{species}</data>')
                lines.append('    </node>')

            # Add edges
            edge_id = 0
            for rxn in reactions_data:
                for reactant in rxn['reactants'].keys():
                    for product in rxn['products'].keys():
                        lines.append(f'    <edge id="e{edge_id}" source="{reactant}" target="{product}">')
                        lines.append(f'      <data key="reaction">{rxn["name"]}</data>')
                        lines.append('    </edge>')
                        edge_id += 1

            lines.append('  </graph>')
            lines.append('</graphml>')
            content = '\n'.join(lines)

        else:
            return jsonify({'success': False, 'error': f'Unsupported format: {export_format}'})

        return jsonify({
            'success': True,
            'content': content,
            'format': export_format
        })

    except Exception as e:
        return error_response(e, context='pathway_routes.export_pathway')
