"""
Routes for pathway analysis - Complete Bio.Pathway implementation
Supports: Reaction, System, Network, MultiGraph
"""
from flask import Blueprint, request, jsonify
import base64
from io import BytesIO
import json

from dependencies import Pathway, plt, np

bp = Blueprint('pathway', __name__, url_prefix='/api')


# ============================================================================
# SYSTEM ANALYSIS - Bio.Pathway.System
# ============================================================================

@bp.route('/pathway/analyze_system', methods=['POST'])
def analyze_system():
    """Analyze complete reaction system using Bio.Pathway.System"""
    try:
        from Bio.Pathway import System, Reaction

        data = request.json
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
        return jsonify({'success': False, 'error': str(e)})


# ============================================================================
# NETWORK ANALYSIS - Bio.Pathway.Network
# ============================================================================

@bp.route('/pathway/analyze_network', methods=['POST'])
def analyze_network():
    """Analyze pathway network using Bio.Pathway.Network"""
    try:
        from Bio.Pathway import Network

        data = request.json
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
        return jsonify({'success': False, 'error': str(e)})


# ============================================================================
# MULTIGRAPH VISUALIZATION - Bio.Pathway.MultiGraph
# ============================================================================

@bp.route('/pathway/visualize', methods=['POST'])
def visualize_pathway():
    """Generate pathway visualization using MultiGraph"""
    try:
        data = request.json
        reactions_data = data.get('reactions', [])

        if not reactions_data:
            return jsonify({'success': False, 'error': 'No reactions provided'})

        # Build network graph for visualization
        import networkx as nx

        G = nx.DiGraph()

        # Add nodes and edges
        for rxn_data in reactions_data:
            for reactant in rxn_data['reactants'].keys():
                G.add_node(reactant, node_type='species')

            for product in rxn_data['products'].keys():
                G.add_node(product, node_type='species')

            # Add edges
            for reactant in rxn_data['reactants'].keys():
                for product in rxn_data['products'].keys():
                    G.add_edge(reactant, product,
                             reaction=rxn_data['name'],
                             reversible=rxn_data['reversible'])

        # Create visualization
        plt.figure(figsize=(12, 8))
        pos = nx.spring_layout(G, k=2, iterations=50)

        # Draw nodes
        node_colors = []
        for node in G.nodes():
            in_degree = G.in_degree(node)
            out_degree = G.out_degree(node)

            if in_degree == 0 and out_degree > 0:
                node_colors.append('lightgreen')  # Source
            elif out_degree == 0 and in_degree > 0:
                node_colors.append('lightcoral')  # Sink
            else:
                node_colors.append('lightblue')  # Intermediate

        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=3000, alpha=0.9)
        nx.draw_networkx_labels(G, pos, font_size=10, font_weight='bold')

        # Draw edges
        nx.draw_networkx_edges(G, pos, edge_color='gray', arrows=True,
                              arrowsize=20, arrowstyle='->', width=2, alpha=0.6)

        plt.title('Pathway Network Visualization', fontsize=16, fontweight='bold')
        plt.axis('off')
        plt.tight_layout()

        # Save to buffer
        img_buffer = BytesIO()
        plt.savefig(img_buffer, format='png', bbox_inches='tight', dpi=150, facecolor='white')
        img_buffer.seek(0)
        img_base64 = base64.b64encode(img_buffer.getvalue()).decode()
        plt.close()

        return jsonify({
            'success': True,
            'visualization': {
                'graph_image': f'data:image/png;base64,{img_base64}',
                'node_count': G.number_of_nodes(),
                'edge_count': G.number_of_edges(),
                'graph_type': 'Directed Graph'
            }
        })

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


# ============================================================================
# EXPORT FORMATS
# ============================================================================

@bp.route('/pathway/export', methods=['POST'])
def export_pathway():
    """Export pathway in various formats"""
    try:
        data = request.json
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
        return jsonify({'success': False, 'error': str(e)})
