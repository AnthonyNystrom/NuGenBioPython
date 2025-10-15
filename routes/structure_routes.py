"""
Routes for protein structure analysis
"""
from flask import Blueprint, request, jsonify, send_file, current_app
import os
from werkzeug.utils import secure_filename
import numpy as np
from io import StringIO

from dependencies import PDBParser, DSSP, PPBuilder, is_aa, NeighborSearch
from Bio.PDB import PDBIO, Superimposer, Select, Vector, calc_dihedral, calc_angle

bp = Blueprint('structure', __name__, url_prefix='/api')


@bp.route('/structure/sample')
def download_sample_pdb():
    """Download a sample PDB file for testing"""
    try:
        sample_file_path = os.path.join(current_app.static_folder, 'sample_protein.pdb')
        if os.path.exists(sample_file_path):
            return send_file(sample_file_path, as_attachment=True, download_name='sample_protein.pdb')
        else:
            return jsonify({'success': False, 'error': 'Sample PDB file not found'})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/structure/parse', methods=['POST'])
def parse_structure():
    try:
        file = request.files['file']

        if file:
            filename = secure_filename(file.filename)
            filepath = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
            file.save(filepath)

            parser = PDBParser(QUIET=True)
            structure = parser.get_structure('structure', filepath)

            models_list = list(structure)
            structure_info = {
                'models': len(models_list),
                'chains': [],
                'residue_count': 0,
                'atom_count': 0,
                'polypeptides': [],
                'secondary_structure': None,
                'neighbor_analysis': None
            }

            # Enhanced chain analysis with polypeptides
            ppb = PPBuilder()
            polypeptides = ppb.build_peptides(structure)

            for model in structure:
                # Neighbor search analysis
                try:
                    atom_list = list(model.get_atoms())
                    if len(atom_list) > 0:
                        ns = NeighborSearch(atom_list)
                        # Find atoms within 5Å of first atom
                        if len(atom_list) > 1:
                            nearby = ns.search(atom_list[0].coord, 5.0, level='A')
                            structure_info['neighbor_analysis'] = {
                                'total_atoms': len(atom_list),
                                'nearby_atoms_5A': len(nearby),
                                'density': len(nearby) / len(atom_list) if len(atom_list) > 0 else 0
                            }
                except Exception as e:
                    structure_info['neighbor_analysis'] = {'error': str(e)}

                for chain in model:
                    # Count amino acid residues
                    aa_residues = [res for res in chain if is_aa(res)]

                    chain_info = {
                        'id': chain.id,
                        'residue_count': len(list(chain)),
                        'aa_residues': len(aa_residues),
                        'atoms': sum(len(list(residue)) for residue in chain),
                        'polypeptide_count': 0
                    }
                    structure_info['chains'].append(chain_info)
                    structure_info['residue_count'] += chain_info['residue_count']
                    structure_info['atom_count'] += chain_info['atoms']

            # Analyze polypeptides
            for i, pp in enumerate(polypeptides):
                try:
                    sequence = pp.get_sequence()
                    phi_psi_list = pp.get_phi_psi_list()
                    ca_list = pp.get_ca_list()

                    pp_info = {
                        'id': i + 1,
                        'length': len(sequence),
                        'sequence': str(sequence)[:50] + ('...' if len(sequence) > 50 else ''),
                        'phi_psi_angles': len([x for x in phi_psi_list if x[0] is not None and x[1] is not None]),
                        'ca_atoms': len(ca_list)
                    }
                    structure_info['polypeptides'].append(pp_info)

                    # Update chain info with polypeptide count
                    for chain_info in structure_info['chains']:
                        if any(residue.get_parent().id == chain_info['id'] for residue in pp):
                            chain_info['polypeptide_count'] += 1
                except Exception as e:
                    structure_info['polypeptides'].append({'error': str(e)})

            # Try DSSP analysis if available
            if DSSP:
                try:
                    model = structure[0]
                    dssp = DSSP(model, filepath)
                    ss_summary = {}
                    for key in dssp.keys():
                        ss = dssp[key][2]  # Secondary structure
                        ss_summary[ss] = ss_summary.get(ss, 0) + 1
                    structure_info['secondary_structure'] = ss_summary
                except Exception as e:
                    structure_info['secondary_structure'] = {'error': str(e)}

            os.remove(filepath)

            return jsonify({
                'success': True,
                'structure_info': structure_info
            })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/structure/advanced_analysis', methods=['POST'])
def advanced_structure_analysis():
    try:
        if 'file' not in request.files:
            return jsonify({'success': False, 'error': 'No file uploaded'})

        file = request.files['file']
        filename = secure_filename(file.filename)
        filepath = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
        file.save(filepath)

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('structure', filepath)
        model = structure[0]  # First model

        analysis_results = {
            'basic_info': {},
            'polypeptide_analysis': [],
            'geometric_analysis': {},
            'secondary_structure': None
        }

        # Basic structure info
        models_list = list(structure)
        chains_list = list(model)
        residues_list = list(model.get_residues())
        atoms_list = list(model.get_atoms())

        analysis_results['basic_info'] = {
            'models': len(models_list),
            'chains': len(chains_list),
            'residues': len(residues_list),
            'atoms': len(atoms_list)
        }

        # Polypeptide analysis
        ppb = PPBuilder()
        for i, pp in enumerate(ppb.build_peptides(structure)):
            try:
                sequence = pp.get_sequence()
                phi_psi_list = pp.get_phi_psi_list()

                # Calculate phi/psi statistics
                valid_angles = [x for x in phi_psi_list if x[0] is not None and x[1] is not None]
                avg_phi = sum(x[0] for x in valid_angles) / len(valid_angles) if valid_angles else None
                avg_psi = sum(x[1] for x in valid_angles) / len(valid_angles) if valid_angles else None

                pp_analysis = {
                    'id': i + 1,
                    'length': len(sequence),
                    'sequence': str(sequence)[:30] + ('...' if len(sequence) > 30 else ''),
                    'phi_psi_count': len(valid_angles),
                    'avg_phi': round(avg_phi, 2) if avg_phi else None,
                    'avg_psi': round(avg_psi, 2) if avg_psi else None
                }
                analysis_results['polypeptide_analysis'].append(pp_analysis)
            except Exception as e:
                analysis_results['polypeptide_analysis'].append({'error': str(e)})

        # Geometric analysis with neighbor search
        try:
            atom_list = list(model.get_atoms())
            if len(atom_list) > 10:
                ns = NeighborSearch(atom_list)
                center_atom = atom_list[len(atom_list)//2]  # Middle atom
                nearby_5A = ns.search(center_atom.coord, 5.0, level='A')
                nearby_10A = ns.search(center_atom.coord, 10.0, level='A')

                analysis_results['geometric_analysis'] = {
                    'center_atom': str(center_atom),
                    'neighbors_5A': len(nearby_5A),
                    'neighbors_10A': len(nearby_10A),
                    'packing_density_5A': len(nearby_5A) / len(atom_list),
                    'packing_density_10A': len(nearby_10A) / len(atom_list)
                }
        except Exception as e:
            analysis_results['geometric_analysis'] = {'error': str(e)}

        # DSSP secondary structure analysis
        if DSSP:
            try:
                dssp = DSSP(model, filepath)
                ss_counts = {}
                accessibility_stats = []

                for key in dssp.keys():
                    ss = dssp[key][2]  # Secondary structure
                    accessibility = dssp[key][3]  # Relative accessibility

                    ss_counts[ss] = ss_counts.get(ss, 0) + 1
                    accessibility_stats.append(accessibility)

                analysis_results['secondary_structure'] = {
                    'structure_counts': ss_counts,
                    'avg_accessibility': sum(accessibility_stats) / len(accessibility_stats) if accessibility_stats else 0,
                    'total_residues_analyzed': len(accessibility_stats)
                }
            except Exception as e:
                analysis_results['secondary_structure'] = {'error': str(e)}

        os.remove(filepath)

        return jsonify({
            'success': True,
            'analysis': analysis_results
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/structure/export', methods=['POST'])
def export_structure():
    """Export structure in PDB format using PDBIO"""
    try:
        file = request.files['file']
        filename = secure_filename(file.filename)
        filepath = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
        file.save(filepath)

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('structure', filepath)

        # Write to string
        io = PDBIO()
        io.set_structure(structure)

        output = StringIO()
        io.save(output)
        pdb_content = output.getvalue()

        os.remove(filepath)

        return jsonify({
            'success': True,
            'pdb_data': pdb_content
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/structure/superimpose', methods=['POST'])
def superimpose_structures():
    """Superimpose two structures and calculate RMSD"""
    try:
        if 'file1' not in request.files or 'file2' not in request.files:
            return jsonify({'success': False, 'error': 'Two structure files required'})

        file1 = request.files['file1']
        file2 = request.files['file2']

        filename1 = secure_filename(file1.filename)
        filename2 = secure_filename(file2.filename)
        filepath1 = os.path.join(current_app.config['UPLOAD_FOLDER'], filename1)
        filepath2 = os.path.join(current_app.config['UPLOAD_FOLDER'], filename2)

        file1.save(filepath1)
        file2.save(filepath2)

        parser = PDBParser(QUIET=True)
        structure1 = parser.get_structure('structure1', filepath1)
        structure2 = parser.get_structure('structure2', filepath2)

        # Get CA atoms for superposition
        model1 = structure1[0]
        model2 = structure2[0]

        # Get all CA atoms
        ca_atoms1 = []
        ca_atoms2 = []

        for chain1 in model1:
            for residue1 in chain1:
                if is_aa(residue1) and 'CA' in residue1:
                    ca_atoms1.append(residue1['CA'])

        for chain2 in model2:
            for residue2 in chain2:
                if is_aa(residue2) and 'CA' in residue2:
                    ca_atoms2.append(residue2['CA'])

        # Use minimum number of atoms
        min_len = min(len(ca_atoms1), len(ca_atoms2))
        if min_len == 0:
            return jsonify({'success': False, 'error': 'No CA atoms found for alignment'})

        ca_atoms1 = ca_atoms1[:min_len]
        ca_atoms2 = ca_atoms2[:min_len]

        # Superimpose
        super_imposer = Superimposer()
        super_imposer.set_atoms(ca_atoms1, ca_atoms2)

        rmsd = super_imposer.rms
        rotation_matrix = super_imposer.rotran[0].tolist()
        translation_vector = super_imposer.rotran[1].tolist()

        os.remove(filepath1)
        os.remove(filepath2)

        return jsonify({
            'success': True,
            'rmsd': float(rmsd),
            'atoms_aligned': min_len,
            'rotation_matrix': rotation_matrix,
            'translation_vector': translation_vector
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/structure/geometry', methods=['POST'])
def geometric_analysis():
    """Calculate distances, angles, and dihedrals"""
    try:
        file = request.files['file']
        chain_id = request.form.get('chain_id', 'A')

        filename = secure_filename(file.filename)
        filepath = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
        file.save(filepath)

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('structure', filepath)
        model = structure[0]

        if chain_id not in model:
            return jsonify({'success': False, 'error': f'Chain {chain_id} not found'})

        chain = model[chain_id]
        residues = [res for res in chain if is_aa(res)]

        geometry_data = {
            'distances': [],
            'angles': [],
            'dihedrals': []
        }

        # Calculate consecutive CA-CA distances
        for i in range(len(residues) - 1):
            if 'CA' in residues[i] and 'CA' in residues[i+1]:
                ca1 = residues[i]['CA']
                ca2 = residues[i+1]['CA']
                distance = ca1 - ca2  # Distance in Angstroms
                geometry_data['distances'].append({
                    'residue1': f"{residues[i].get_resname()}{residues[i].id[1]}",
                    'residue2': f"{residues[i+1].get_resname()}{residues[i+1].id[1]}",
                    'distance': round(float(distance), 3)
                })

        # Calculate CA-CA-CA angles
        for i in range(len(residues) - 2):
            if all('CA' in residues[j] for j in [i, i+1, i+2]):
                v1 = residues[i]['CA'].get_vector()
                v2 = residues[i+1]['CA'].get_vector()
                v3 = residues[i+2]['CA'].get_vector()
                angle = calc_angle(v1, v2, v3)
                geometry_data['angles'].append({
                    'residues': f"{residues[i].id[1]}-{residues[i+1].id[1]}-{residues[i+2].id[1]}",
                    'angle_degrees': round(float(np.degrees(angle)), 2)
                })

        # Calculate dihedrals for first 10 residues
        for i in range(min(10, len(residues) - 3)):
            if all('CA' in residues[j] for j in [i, i+1, i+2, i+3]):
                v1 = residues[i]['CA'].get_vector()
                v2 = residues[i+1]['CA'].get_vector()
                v3 = residues[i+2]['CA'].get_vector()
                v4 = residues[i+3]['CA'].get_vector()
                dihedral = calc_dihedral(v1, v2, v3, v4)
                geometry_data['dihedrals'].append({
                    'residues': f"{residues[i].id[1]}-{residues[i+1].id[1]}-{residues[i+2].id[1]}-{residues[i+3].id[1]}",
                    'dihedral_degrees': round(float(np.degrees(dihedral)), 2)
                })

        os.remove(filepath)

        return jsonify({
            'success': True,
            'geometry': geometry_data,
            'chain': chain_id,
            'residue_count': len(residues)
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/structure/quality', methods=['POST'])
def quality_metrics():
    """Analyze structure quality metrics"""
    try:
        file = request.files['file']
        filename = secure_filename(file.filename)
        filepath = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
        file.save(filepath)

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('structure', filepath)
        model = structure[0]

        quality_data = {
            'bfactor_stats': {},
            'occupancy_stats': {},
            'atoms_by_element': {},
            'residue_completeness': []
        }

        # Collect B-factors and occupancy
        bfactors = []
        occupancies = []
        elements = {}

        for atom in model.get_atoms():
            bfactors.append(atom.get_bfactor())
            occupancies.append(atom.get_occupancy())
            element = atom.element
            elements[element] = elements.get(element, 0) + 1

        if bfactors:
            quality_data['bfactor_stats'] = {
                'mean': round(np.mean(bfactors), 2),
                'median': round(np.median(bfactors), 2),
                'min': round(np.min(bfactors), 2),
                'max': round(np.max(bfactors), 2),
                'std': round(np.std(bfactors), 2)
            }

        if occupancies:
            quality_data['occupancy_stats'] = {
                'mean': round(np.mean(occupancies), 3),
                'min': round(np.min(occupancies), 3),
                'max': round(np.max(occupancies), 3)
            }

        quality_data['atoms_by_element'] = elements

        # Check residue completeness
        for chain in model:
            for residue in chain:
                if is_aa(residue):
                    expected_atoms = ['N', 'CA', 'C', 'O']
                    present_atoms = [atom.id for atom in residue]
                    missing = [a for a in expected_atoms if a not in present_atoms]
                    if missing:
                        quality_data['residue_completeness'].append({
                            'chain': chain.id,
                            'residue': f"{residue.get_resname()}{residue.id[1]}",
                            'missing_atoms': missing
                        })

        os.remove(filepath)

        return jsonify({
            'success': True,
            'quality': quality_data
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/structure/contacts', methods=['POST'])
def contact_analysis():
    """Calculate residue-residue contacts and distance matrix"""
    try:
        file = request.files['file']
        chain_id = request.form.get('chain_id', 'A')
        cutoff = float(request.form.get('cutoff', 8.0))

        filename = secure_filename(file.filename)
        filepath = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
        file.save(filepath)

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('structure', filepath)
        model = structure[0]

        if chain_id not in model:
            return jsonify({'success': False, 'error': f'Chain {chain_id} not found'})

        chain = model[chain_id]
        residues = [res for res in chain if is_aa(res) and 'CA' in res]

        # Calculate distance matrix (CA-CA distances)
        n = len(residues)
        distance_matrix = []
        contacts = []

        for i in range(min(n, 20)):  # Limit to first 20 residues for performance
            row = []
            for j in range(min(n, 20)):
                if 'CA' in residues[i] and 'CA' in residues[j]:
                    distance = residues[i]['CA'] - residues[j]['CA']
                    row.append(round(float(distance), 2))

                    # Record contacts (non-adjacent residues within cutoff)
                    if abs(i - j) > 1 and distance < cutoff:
                        contacts.append({
                            'residue1': f"{residues[i].get_resname()}{residues[i].id[1]}",
                            'residue2': f"{residues[j].get_resname()}{residues[j].id[1]}",
                            'distance': round(float(distance), 2)
                        })
                else:
                    row.append(None)
            distance_matrix.append(row)

        os.remove(filepath)

        return jsonify({
            'success': True,
            'distance_matrix': distance_matrix,
            'contacts': contacts[:50],  # Limit to 50 contacts
            'contact_count': len(contacts),
            'cutoff': cutoff,
            'residues_analyzed': min(n, 20)
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/structure/hbonds', methods=['POST'])
def hydrogen_bonds():
    """Detect potential hydrogen bonds and salt bridges"""
    try:
        file = request.files['file']
        filename = secure_filename(file.filename)
        filepath = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
        file.save(filepath)

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('structure', filepath)
        model = structure[0]

        hbonds = []
        salt_bridges = []

        # Get all atoms
        atoms = list(model.get_atoms())
        ns = NeighborSearch(atoms)

        # Simplified H-bond detection (N-O distance < 3.5Å)
        for atom in atoms:
            if atom.element == 'N':
                # Find nearby oxygen atoms
                nearby = ns.search(atom.coord, 3.5, level='A')
                for nearby_atom in nearby:
                    if nearby_atom.element == 'O' and nearby_atom.get_parent() != atom.get_parent():
                        distance = atom - nearby_atom
                        hbonds.append({
                            'donor': f"{atom.get_parent().get_resname()}{atom.get_parent().id[1]}:{atom.id}",
                            'acceptor': f"{nearby_atom.get_parent().get_resname()}{nearby_atom.get_parent().id[1]}:{nearby_atom.id}",
                            'distance': round(float(distance), 2)
                        })

        # Salt bridge detection (charged residue pairs within 4Å)
        positive_residues = ['ARG', 'LYS', 'HIS']
        negative_residues = ['ASP', 'GLU']

        for chain in model:
            for residue in chain:
                if is_aa(residue) and residue.get_resname() in positive_residues:
                    if 'CA' not in residue:
                        continue
                    # Find nearby negative residues
                    nearby = ns.search(residue['CA'].coord, 4.0, level='R')
                    for nearby_res in nearby:
                        if is_aa(nearby_res) and nearby_res.get_resname() in negative_residues:
                            if nearby_res != residue:
                                distance = residue['CA'] - nearby_res['CA']
                                salt_bridges.append({
                                    'residue1': f"{residue.get_resname()}{residue.id[1]}",
                                    'residue2': f"{nearby_res.get_resname()}{nearby_res.id[1]}",
                                    'distance': round(float(distance), 2)
                                })

        os.remove(filepath)

        return jsonify({
            'success': True,
            'hbonds': hbonds[:50],  # Limit to 50
            'hbond_count': len(hbonds),
            'salt_bridges': salt_bridges[:20],  # Limit to 20
            'salt_bridge_count': len(salt_bridges)
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/structure/secondary_structure', methods=['POST'])
def secondary_structure_analysis():
    """Enhanced DSSP secondary structure analysis"""
    try:
        file = request.files['file']
        filename = secure_filename(file.filename)
        filepath = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
        file.save(filepath)

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('structure', filepath)
        model = structure[0]

        result = {
            'dssp_available': False,
            'secondary_structure': {},
            'residue_details': []
        }

        # Try DSSP analysis
        if DSSP:
            try:
                dssp = DSSP(model, filepath)
                ss_counts = {}
                ss_mapping = {
                    'H': 'Alpha helix',
                    'B': 'Beta bridge',
                    'E': 'Extended strand',
                    'G': '3-helix',
                    'I': '5-helix',
                    'T': 'Turn',
                    'S': 'Bend',
                    '-': 'Loop/coil'
                }

                for key in dssp.keys():
                    ss = dssp[key][2]
                    acc = dssp[key][3]  # Accessibility
                    phi = dssp[key][4]  # Phi angle
                    psi = dssp[key][5]  # Psi angle

                    ss_counts[ss] = ss_counts.get(ss, 0) + 1

                    # Store first 50 residues
                    if len(result['residue_details']) < 50:
                        result['residue_details'].append({
                            'chain': key[0],
                            'residue': key[1][1],
                            'ss': ss,
                            'ss_name': ss_mapping.get(ss, 'Unknown'),
                            'accessibility': round(acc, 2),
                            'phi': round(phi, 1) if phi else None,
                            'psi': round(psi, 1) if psi else None
                        })

                result['dssp_available'] = True
                result['secondary_structure'] = {
                    'counts': ss_counts,
                    'mapping': ss_mapping,
                    'total_residues': sum(ss_counts.values())
                }
            except Exception as e:
                result['dssp_error'] = str(e)

        os.remove(filepath)

        return jsonify({
            'success': True,
            'result': result
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/structure/ramachandran', methods=['POST'])
def ramachandran_analysis():
    """Calculate phi/psi angles for Ramachandran plot"""
    try:
        file = request.files['file']
        chain_id = request.form.get('chain_id', 'A')

        filename = secure_filename(file.filename)
        filepath = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
        file.save(filepath)

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('structure', filepath)
        model = structure[0]

        if chain_id not in model:
            return jsonify({'success': False, 'error': f'Chain {chain_id} not found'})

        chain = model[chain_id]
        ppb = PPBuilder()
        phi_psi_data = []

        for pp in ppb.build_peptides(chain):
            phi_psi_list = pp.get_phi_psi_list()
            for i, (phi, psi) in enumerate(phi_psi_list):
                if phi and psi:  # Both angles available
                    residue = list(pp)[i]
                    phi_psi_data.append({
                        'residue': f"{residue.get_resname()}{residue.id[1]}",
                        'phi': round(np.degrees(phi), 2),
                        'psi': round(np.degrees(psi), 2),
                        'resname': residue.get_resname()
                    })

        # Classify into regions
        favored = 0
        allowed = 0
        outliers = 0

        for data in phi_psi_data:
            phi, psi = data['phi'], data['psi']
            # Simplified classification
            if (-180 <= phi <= -30 and -180 <= psi <= 50):  # Beta region
                favored += 1
            elif (-90 <= phi <= -30 and -75 <= psi <= -10):  # Alpha region
                favored += 1
            elif (-180 <= phi <= 180 and -180 <= psi <= 180):  # Allowed
                allowed += 1
            else:
                outliers += 1

        os.remove(filepath)

        return jsonify({
            'success': True,
            'phi_psi_data': phi_psi_data[:100],  # Limit to 100 for display
            'total_residues': len(phi_psi_data),
            'classification': {
                'favored': favored,
                'allowed': allowed,
                'outliers': outliers
            }
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/structure/sasa', methods=['POST'])
def sasa_analysis():
    """Calculate Solvent Accessible Surface Area"""
    try:
        file = request.files['file']
        filename = secure_filename(file.filename)
        filepath = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
        file.save(filepath)

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('structure', filepath)
        model = structure[0]

        try:
            from Bio.PDB import SASA
            sr = SASA.ShrakeRupley()
            sr.compute(structure, level="S")
            total_sasa = structure.sasa

            # Per-chain SASA
            chain_sasa = {}
            for chain in model:
                sr.compute(chain, level="C")
                chain_sasa[chain.id] = round(chain.sasa, 2)

            # Per-residue SASA (first 50)
            residue_sasa = []
            for chain in model:
                for residue in list(chain)[:50]:
                    if is_aa(residue):
                        sr.compute(residue, level="R")
                        if hasattr(residue, 'sasa'):
                            residue_sasa.append({
                                'chain': chain.id,
                                'residue': f"{residue.get_resname()}{residue.id[1]}",
                                'sasa': round(residue.sasa, 2)
                            })

            os.remove(filepath)

            return jsonify({
                'success': True,
                'total_sasa': round(total_sasa, 2),
                'chain_sasa': chain_sasa,
                'residue_sasa': residue_sasa
            })
        except ImportError:
            # Fallback: Use DSSP accessibility as proxy
            if DSSP:
                try:
                    dssp = DSSP(model, filepath)
                    accessibility_data = []
                    chain_acc = {}

                    for key in dssp.keys():
                        acc = dssp[key][3]
                        chain = key[0]
                        chain_acc[chain] = chain_acc.get(chain, 0) + acc

                        if len(accessibility_data) < 50:
                            accessibility_data.append({
                                'chain': chain,
                                'residue': key[1][1],
                                'relative_accessibility': round(acc, 2)
                            })

                    os.remove(filepath)

                    return jsonify({
                        'success': True,
                        'method': 'DSSP (relative accessibility)',
                        'chain_accessibility': {k: round(v, 2) for k, v in chain_acc.items()},
                        'residue_accessibility': accessibility_data
                    })
                except Exception as e:
                    os.remove(filepath)
                    return jsonify({'success': False, 'error': f'DSSP failed: {str(e)}'})
            else:
                os.remove(filepath)
                return jsonify({'success': False, 'error': 'SASA module not available, DSSP not installed'})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/structure/residue_depth', methods=['POST'])
def residue_depth_analysis():
    """Calculate residue depth (burial)"""
    try:
        file = request.files['file']
        chain_id = request.form.get('chain_id', 'A')

        filename = secure_filename(file.filename)
        filepath = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
        file.save(filepath)

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('structure', filepath)
        model = structure[0]

        if chain_id not in model:
            return jsonify({'success': False, 'error': f'Chain {chain_id} not found'})

        try:
            from Bio.PDB.ResidueDepth import ResidueDepth
            rd = ResidueDepth(model)
            depth_data = []

            chain = model[chain_id]
            for residue in chain:
                if is_aa(residue) and residue in rd:
                    depth, ca_depth = rd[residue]
                    depth_data.append({
                        'residue': f"{residue.get_resname()}{residue.id[1]}",
                        'depth': round(depth, 2),
                        'ca_depth': round(ca_depth, 2)
                    })

            # Calculate statistics
            depths = [d['depth'] for d in depth_data]
            stats = {}
            if depths:
                stats = {
                    'mean': round(np.mean(depths), 2),
                    'median': round(np.median(depths), 2),
                    'min': round(np.min(depths), 2),
                    'max': round(np.max(depths), 2)
                }

            os.remove(filepath)

            return jsonify({
                'success': True,
                'depth_data': depth_data[:50],
                'total_residues': len(depth_data),
                'statistics': stats
            })
        except ImportError:
            os.remove(filepath)
            return jsonify({'success': False, 'error': 'ResidueDepth module requires MSMS (not available)'})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/structure/hse', methods=['POST'])
def hse_analysis():
    """Calculate Half Sphere Exposure"""
    try:
        file = request.files['file']
        chain_id = request.form.get('chain_id', 'A')

        filename = secure_filename(file.filename)
        filepath = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
        file.save(filepath)

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('structure', filepath)
        model = structure[0]

        if chain_id not in model:
            return jsonify({'success': False, 'error': f'Chain {chain_id} not found'})

        try:
            from Bio.PDB.HSExposure import HSExposureCA, HSExposureCB

            # Calculate HSE using CA atoms
            hse_ca = HSExposureCA(model)
            hse_data = []

            chain = model[chain_id]
            for residue in chain:
                if is_aa(residue) and residue in hse_ca:
                    hse_up, hse_down, angle = hse_ca[residue]
                    hse_data.append({
                        'residue': f"{residue.get_resname()}{residue.id[1]}",
                        'hse_up': round(hse_up, 2),
                        'hse_down': round(hse_down, 2),
                        'angle': round(np.degrees(angle), 2) if angle else None,
                        'cn': hse_up + hse_down  # Contact number
                    })

            # Statistics
            cn_values = [d['cn'] for d in hse_data]
            stats = {}
            if cn_values:
                stats = {
                    'mean_cn': round(np.mean(cn_values), 2),
                    'median_cn': round(np.median(cn_values), 2)
                }

            os.remove(filepath)

            return jsonify({
                'success': True,
                'hse_data': hse_data[:50],
                'total_residues': len(hse_data),
                'statistics': stats
            })
        except ImportError:
            os.remove(filepath)
            return jsonify({'success': False, 'error': 'HSExposure module not available'})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/structure/extract', methods=['POST'])
def extract_structure():
    """Extract specific chains, residues, or atoms"""
    try:
        file = request.files['file']
        selection_type = request.form.get('selection_type', 'chain')
        selection_value = request.form.get('selection_value', 'A')

        filename = secure_filename(file.filename)
        filepath = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
        file.save(filepath)

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('structure', filepath)
        model = structure[0]

        extraction_info = {
            'selection_type': selection_type,
            'selection_value': selection_value,
            'extracted_count': 0,
            'details': []
        }

        if selection_type == 'chain':
            if selection_value in model:
                chain = model[selection_value]
                extraction_info['extracted_count'] = len(list(chain.get_residues()))
                extraction_info['details'] = [{
                    'chain': chain.id,
                    'residues': len(list(chain.get_residues())),
                    'atoms': len(list(chain.get_atoms()))
                }]
            else:
                return jsonify({'success': False, 'error': f'Chain {selection_value} not found'})

        elif selection_type == 'residue_range':
            # Format: "A:1-50" (chain:start-end)
            try:
                parts = selection_value.split(':')
                chain_id = parts[0]
                range_parts = parts[1].split('-')
                start_res = int(range_parts[0])
                end_res = int(range_parts[1])

                if chain_id in model:
                    chain = model[chain_id]
                    extracted_residues = []
                    for residue in chain:
                        if start_res <= residue.id[1] <= end_res:
                            extracted_residues.append({
                                'residue': f"{residue.get_resname()}{residue.id[1]}",
                                'atoms': len(list(residue.get_atoms()))
                            })
                    extraction_info['extracted_count'] = len(extracted_residues)
                    extraction_info['details'] = extracted_residues[:50]
                else:
                    return jsonify({'success': False, 'error': f'Chain {chain_id} not found'})
            except (IndexError, ValueError) as e:
                return jsonify({'success': False, 'error': 'Invalid residue range format. Use "A:1-50"'})

        elif selection_type == 'atom_type':
            # Extract specific atom types (e.g., "CA")
            atom_count = 0
            for chain in model:
                for residue in chain:
                    if selection_value in residue:
                        atom = residue[selection_value]
                        atom_count += 1
                        if len(extraction_info['details']) < 50:
                            extraction_info['details'].append({
                                'chain': chain.id,
                                'residue': f"{residue.get_resname()}{residue.id[1]}",
                                'atom': atom.id
                            })
            extraction_info['extracted_count'] = atom_count

        os.remove(filepath)

        return jsonify({
            'success': True,
            'extraction': extraction_info
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})
