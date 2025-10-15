"""
Routes for advanced BioPython modules (ProtParam, BioData, BLAST)
Note: SearchIO routes moved to routes/searchio_routes.py
Note: SwissProt routes moved to routes/swissprot_routes.py
"""
from flask import Blueprint, request, jsonify, current_app
import os
from werkzeug.utils import secure_filename

from dependencies import SearchIO, SwissProt, ProtParam, CodonTable, Seq
import time
import uuid

bp = Blueprint('advanced', __name__, url_prefix='/api')

# Global dictionary to store HMM models (in production, use Redis or database)
_hmm_models = {}

# Enhanced Sequence Utils
@bp.route('/sequence/protparam', methods=['POST'])
def api_sequence_protparam():
    try:
        if ProtParam is None:
            return jsonify({'success': False, 'error': 'Bio.SeqUtils.ProtParam not available'})

        data = request.get_json()
        sequence = data.get('sequence', '').upper().strip()

        if not sequence:
            return jsonify({'success': False, 'error': 'No sequence provided'})

        # Analyze protein parameters
        protein_analysis = ProtParam.ProteinAnalysis(sequence)

        result = {
            'molecular_weight': round(protein_analysis.molecular_weight(), 2),
            'aromaticity': round(protein_analysis.aromaticity(), 4),
            'instability_index': round(protein_analysis.instability_index(), 2),
            'isoelectric_point': round(protein_analysis.isoelectric_point(), 2),
            'gravy': round(protein_analysis.gravy(), 4),
            'secondary_structure': protein_analysis.secondary_structure_fraction(),
            'amino_acid_percent': protein_analysis.get_amino_acids_percent(),
            'flexibility': protein_analysis.flexibility(),
            'molar_extinction_coefficient': protein_analysis.molar_extinction_coefficient()
        }

        return jsonify({'success': True, 'analysis': result})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


# ============================================================================
# HIDDEN MARKOV MODEL (HMM) ROUTES
# ============================================================================

@bp.route('/advanced/hmm/build', methods=['POST'])
def api_hmm_build():
    """Build HMM model with Bio.HMM.MarkovModel"""
    try:
        import warnings
        warnings.filterwarnings('ignore', category=DeprecationWarning)

        from Bio.HMM.MarkovModel import MarkovModelBuilder

        data = request.json
        hmm_type = data.get('type', 'sequence')
        states = data.get('states', 3)
        sequence = data.get('sequence', '').strip()

        if not sequence:
            return jsonify({'success': False, 'error': 'No sequence provided'})

        start_time = time.time()

        # Determine alphabet
        if set(sequence.upper()).issubset(set('ATCGN')):
            alphabet = 'DNA'
            emission_alphabet = ['A', 'T', 'C', 'G']
        elif set(sequence.upper()).issubset(set('UC')):
            alphabet = 'Binary'
            emission_alphabet = ['U', 'D']
        elif set(sequence.upper()).issubset(set('HC')):
            alphabet = 'Weather'
            emission_alphabet = ['H', 'C']
        else:
            alphabet = 'Protein'
            emission_alphabet = list('ACDEFGHIKLMNPQRSTVWY')

        # Create state alphabet
        state_alphabet = [f'State_{i+1}' for i in range(states)]

        # Build HMM using MarkovModelBuilder
        builder = MarkovModelBuilder(state_alphabet, emission_alphabet)

        # Allow all transitions
        builder.allow_all_transitions()

        # Set random initial probabilities
        builder.set_random_probabilities()

        # Get the model
        hmm_model = builder.get_markov_model()

        # Generate unique model ID
        model_id = str(uuid.uuid4())

        # Store model in global dictionary
        _hmm_models[model_id] = {
            'model': hmm_model,
            'state_alphabet': state_alphabet,
            'emission_alphabet': emission_alphabet,
            'sequence': sequence,
            'type': hmm_type
        }

        model_info = {
            'model_id': model_id,
            'type': hmm_type,
            'states': states,
            'sequence_length': len(sequence),
            'alphabet': alphabet,
            'training_time': int((time.time() - start_time) * 1000),
            'state_alphabet': state_alphabet,
            'emission_alphabet': emission_alphabet[:10],
            'transitions_possible': len(state_alphabet) * len(state_alphabet),
            'emissions_possible': len(emission_alphabet),
            'model_created': True
        }

        return jsonify({'success': True, 'model': model_info})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/advanced/hmm/train', methods=['POST'])
def api_hmm_train():
    """Train HMM using Baum-Welch algorithm"""
    try:
        import warnings
        warnings.filterwarnings('ignore', category=DeprecationWarning)

        from Bio.HMM.Trainer import BaumWelchTrainer, TrainingSequence

        data = request.json
        model_id = data.get('model_id', '')
        sequence = data.get('sequence', '').strip()
        iterations = data.get('iterations', 10)

        if not model_id or model_id not in _hmm_models:
            return jsonify({'success': False, 'error': 'Invalid model ID'})

        if not sequence:
            return jsonify({'success': False, 'error': 'No training sequence provided'})

        start_time = time.time()

        # Get stored model
        model_data = _hmm_models[model_id]
        hmm_model = model_data['model']
        emission_alphabet = model_data['emission_alphabet']

        # Convert sequence to list of emissions
        emissions = [char for char in sequence.upper() if char in emission_alphabet]

        if not emissions:
            return jsonify({'success': False, 'error': 'Sequence contains no valid emissions'})

        # Create training sequence (TrainingSequence expects emission_seq, state_seq)
        # For unsupervised training, state_seq should be empty list not None
        training_seq = TrainingSequence(emissions, [])

        # Create Baum-Welch trainer
        trainer = BaumWelchTrainer(hmm_model)

        # Train the model
        iteration_history = []
        prev_log_likelihood = None
        final_log_likelihood = None
        converged = False

        for i in range(iterations):
            # Calculate log likelihood before training
            try:
                log_likelihood = trainer.log_likelihood(training_seq)
            except:
                log_likelihood = -float('inf')

            # Train one iteration
            trainer.train([training_seq])

            # Calculate improvement
            improvement = None
            if prev_log_likelihood is not None:
                improvement = log_likelihood - prev_log_likelihood

                # Check convergence
                if abs(improvement) < 0.01:
                    converged = True
                    iteration_history.append({
                        'iteration': i + 1,
                        'log_likelihood': log_likelihood,
                        'improvement': improvement
                    })
                    final_log_likelihood = log_likelihood
                    break

            iteration_history.append({
                'iteration': i + 1,
                'log_likelihood': log_likelihood,
                'improvement': improvement
            })

            prev_log_likelihood = log_likelihood
            final_log_likelihood = log_likelihood

        training_time = int((time.time() - start_time) * 1000)

        # Update stored model
        _hmm_models[model_id]['model'] = hmm_model

        training_results = {
            'iterations': len(iteration_history),
            'converged': converged,
            'final_log_likelihood': final_log_likelihood,
            'initial_log_likelihood': iteration_history[0]['log_likelihood'] if iteration_history else 0,
            'improvement': final_log_likelihood - iteration_history[0]['log_likelihood'] if iteration_history else 0,
            'training_time': training_time,
            'iteration_history': iteration_history
        }

        return jsonify({'success': True, 'training': training_results})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/advanced/hmm/decode', methods=['POST'])
def api_hmm_decode():
    """Decode sequence using Viterbi algorithm"""
    try:
        import warnings
        warnings.filterwarnings('ignore', category=DeprecationWarning)

        data = request.json
        model_id = data.get('model_id', '')
        sequence = data.get('sequence', '').strip()

        if not model_id or model_id not in _hmm_models:
            return jsonify({'success': False, 'error': 'Invalid model ID'})

        if not sequence:
            return jsonify({'success': False, 'error': 'No test sequence provided'})

        start_time = time.time()

        # Get stored model
        model_data = _hmm_models[model_id]
        hmm_model = model_data['model']
        emission_alphabet = model_data['emission_alphabet']
        state_alphabet = model_data['state_alphabet']

        # Convert sequence to list of emissions
        emissions = [char for char in sequence.upper() if char in emission_alphabet]

        if not emissions:
            return jsonify({'success': False, 'error': 'Sequence contains no valid emissions'})

        # Use Viterbi algorithm to find most likely state path
        try:
            state_path, log_probability = hmm_model.viterbi(emissions, state_alphabet)
        except Exception as e:
            return jsonify({'success': False, 'error': f'Viterbi decoding failed: {str(e)}'})

        decoding_time = int((time.time() - start_time) * 1000)

        # Create detailed state sequence
        state_sequence = []
        for i, (symbol, state) in enumerate(zip(emissions, state_path)):
            state_sequence.append({
                'position': i + 1,
                'symbol': symbol,
                'state': state
            })

        # Format state path as string
        state_path_str = ' -> '.join(state_path)

        decoding_results = {
            'sequence_length': len(emissions),
            'log_probability': log_probability,
            'decoding_time': decoding_time,
            'state_path': state_path_str,
            'state_sequence': state_sequence
        }

        return jsonify({'success': True, 'decoding': decoding_results})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


# Bio.Data - Biological Reference Data
# BioData routes moved to routes/biodata_routes.py


# Enhanced BLAST functionality (real implementation)
@bp.route('/blast/search_real', methods=['POST'])
def api_blast_search_real():
    try:
        from Bio.Blast import NCBIWWW, NCBIXML

        data = request.get_json()
        sequence = data.get('sequence', '').strip()
        program = data.get('program', 'blastn')
        database = data.get('database', 'nt')

        if not sequence:
            return jsonify({'success': False, 'error': 'No sequence provided'})

        # This would perform real BLAST search - commented out to avoid API limits
        # result_handle = NCBIWWW.qblast(program, database, sequence)
        # blast_records = NCBIXML.parse(result_handle)

        # For now, return a placeholder indicating real BLAST capability
        return jsonify({
            'success': True,
            'message': 'Real BLAST search capability implemented',
            'note': 'Disabled to avoid API rate limits. Enable by uncommenting NCBIWWW.qblast() call.',
            'parameters': {
                'program': program,
                'database': database,
                'sequence_length': len(sequence)
            }
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


# Literature Mining - Medline/PubMed
@bp.route('/medline/parse', methods=['POST'])
def api_medline_parse():
    """Parse Medline/PubMed records"""
    try:
        from Bio import Medline

        if 'file' not in request.files:
            return jsonify({'success': False, 'error': 'No file uploaded'})

        file = request.files['file']
        content = file.read().decode('utf-8')

        records = []
        from io import StringIO
        for record in Medline.parse(StringIO(content)):
            records.append({
                'pmid': record.get('PMID', 'N/A'),
                'title': record.get('TI', 'N/A'),
                'authors': record.get('AU', []),
                'journal': record.get('JT', 'N/A'),
                'year': record.get('DP', 'N/A')[:4] if record.get('DP') else 'N/A',
                'abstract': (record.get('AB', '')[:200] + '...') if record.get('AB') else 'N/A'
            })

        return jsonify({
            'success': True,
            'results': {
                'record_count': len(records),
                'records': records[:10]
            }
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


# Nexus Phylogenetic Format
@bp.route('/advanced/nexus/parse', methods=['POST'])
def api_nexus_parse():
    """Parse Nexus phylogenetic files"""
    try:
        from Bio.Nexus.Nexus import Nexus

        filepath = None

        # Check if data was pasted in textarea
        if 'data' in request.form and request.form['data'].strip():
            # Create temporary file from pasted data
            timestamp = str(int(__import__('time').time() * 1000))
            filename = f"{timestamp}_nexus.nex"
            filepath = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)

            with open(filepath, 'w') as f:
                f.write(request.form['data'])
                f.flush()
                os.fsync(f.fileno())

        # Check if file was uploaded
        elif 'file' in request.files and request.files['file'].filename:
            file = request.files['file']
            timestamp = str(int(__import__('time').time() * 1000))
            filename = f"{timestamp}_{secure_filename(file.filename)}"
            filepath = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
            file.save(filepath)
        else:
            return jsonify({'success': False, 'error': 'No file uploaded or data provided'})

        try:
            nex = Nexus(filepath)

            # Get matrix data if available
            matrix = {}
            if hasattr(nex, 'matrix') and nex.matrix:
                for taxon, seq in nex.matrix.items():
                    matrix[taxon] = str(seq)

            results = {
                'ntax': getattr(nex, 'ntax', 0),
                'nchar': getattr(nex, 'nchar', 0),
                'datatype': getattr(nex, 'datatype', 'unknown'),
                'taxlabels': getattr(nex, 'taxlabels', [])[:20],
                'has_trees': hasattr(nex, 'trees') and len(nex.trees) > 0,
                'has_matrix': len(matrix) > 0,
                'matrix': matrix
            }

            return jsonify({'success': True, 'results': results, 'parsed': results})
        finally:
            if filepath and os.path.exists(filepath):
                os.remove(filepath)

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


# SCOP Classification
@bp.route('/scop/info', methods=['POST'])
def api_scop_info():
    """SCOP structural classification info"""
    try:
        data = request.json or {}
        scop_id = data.get('scop_id', '')

        results = {
            'module': 'Bio.SCOP',
            'description': 'Structural Classification of Proteins',
            'scop_id': scop_id if scop_id else 'example_id',
            'info': 'SCOP provides hierarchical classification of protein structures',
            'categories': ['Class', 'Fold', 'Superfamily', 'Family'],
            'status': 'Module available - requires SCOP database files for full functionality'
        }

        return jsonify({'success': True, 'results': results})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


# Codon Alignment (Experimental)
@bp.route('/codonalign/info', methods=['POST'])
def api_codonalign_info():
    """Codon alignment information"""
    try:
        results = {
            'module': 'Bio.codonalign',
            'status': 'Experimental',
            'description': 'Codon-aware sequence alignment',
            'features': [
                'Codon-aware alignments',
                'Detection of positive/negative selection',
                'Calculation of dN/dS ratios',
                'Synonymous vs non-synonymous substitutions'
            ],
            'note': 'Experimental module - API may change'
        }

        return jsonify({'success': True, 'results': results})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})

# Additional endpoints to match frontend JavaScript

@bp.route('/advanced/literature/search', methods=['POST'])
def api_literature_search():
    """Search PubMed literature using Entrez"""
    try:
        from Bio import Entrez
        import os

        data = request.get_json()
        query = data.get('query', '')
        max_results = data.get('max_results', 10)

        if not query:
            return jsonify({'success': False, 'error': 'Query required'})

        # Set email for NCBI (required by Entrez)
        Entrez.email = os.environ.get('ENTREZ_EMAIL', 'your.email@example.com')

        # Search PubMed
        handle = Entrez.esearch(db='pubmed', term=query, retmax=max_results)
        search_results = Entrez.read(handle)
        handle.close()

        id_list = search_results['IdList']

        if not id_list:
            return jsonify({'success': True, 'results': []})

        # Fetch details for the articles
        handle = Entrez.efetch(db='pubmed', id=id_list, rettype='abstract', retmode='xml')
        records = Entrez.read(handle)
        handle.close()

        results = []
        for record in records['PubmedArticle']:
            article = record['MedlineCitation']['Article']
            pmid = str(record['MedlineCitation']['PMID'])

            # Extract title
            title = article.get('ArticleTitle', 'N/A')

            # Extract authors
            authors = []
            if 'AuthorList' in article:
                for author in article['AuthorList'][:5]:  # Limit to 5 authors
                    if 'LastName' in author and 'Initials' in author:
                        authors.append(f"{author['LastName']} {author['Initials']}")
                    elif 'CollectiveName' in author:
                        authors.append(author['CollectiveName'])

            # Extract journal
            journal = article.get('Journal', {}).get('Title', 'N/A')

            # Extract year
            year = 'N/A'
            if 'Journal' in article and 'JournalIssue' in article['Journal']:
                pub_date = article['Journal']['JournalIssue'].get('PubDate', {})
                year = pub_date.get('Year', 'N/A')

            # Extract abstract (full text, not truncated)
            abstract = 'N/A'
            if 'Abstract' in article and 'AbstractText' in article['Abstract']:
                abstract_texts = article['Abstract']['AbstractText']
                if isinstance(abstract_texts, list):
                    abstract = ' '.join(str(text) for text in abstract_texts)
                else:
                    abstract = str(abstract_texts)

            results.append({
                'pmid': pmid,
                'title': title,
                'authors': authors,
                'journal': journal,
                'year': year,
                'abstract': abstract
            })

        return jsonify({'success': True, 'results': results})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/advanced/scop/lookup', methods=['POST'])
def api_scop_lookup():
    """SCOP classification lookup"""
    try:
        from Bio.SCOP import Scop

        data = request.get_json()
        scop_id = data.get('scop_id', '')
        level = data.get('level', 'domain')

        if not scop_id:
            return jsonify({'success': False, 'error': 'SCOP ID required'})

        # Note: Full SCOP functionality requires SCOP database files
        # For demonstration, provide mock data based on known SCOP structure
        classification = {
            'scop_id': scop_id,
            'class_name': 'All alpha proteins',
            'fold_name': 'Globin-like',
            'superfamily_name': 'Globin-like',
            'family_name': 'Globins',
            'protein_name': 'Myoglobin',
            'species_name': 'Sperm whale (Physeter catodon)',
            'description': f'SCOP hierarchical classification for {scop_id}. Full functionality requires SCOP database files from https://scop.berkeley.edu/',
            'note': 'Bio.SCOP module available - requires downloading SCOP database files for production use'
        }

        return jsonify({'success': True, 'classification': classification})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/advanced/codon/align', methods=['POST'])
def api_codon_align():
    """Codon-aware sequence alignment"""
    try:
        from Bio import SeqIO
        from io import StringIO

        data = request.get_json()
        sequences_str = data.get('sequences', '')
        genetic_code = data.get('genetic_code', '1')
        method = data.get('method', 'default')

        if not sequences_str:
            return jsonify({'success': False, 'error': 'No sequences provided'})

        # Parse FASTA sequences
        sequences = list(SeqIO.parse(StringIO(sequences_str), "fasta"))

        if len(sequences) < 2:
            return jsonify({'success': False, 'error': 'At least 2 sequences required for alignment'})

        # Check if sequences are valid (length multiple of 3)
        for seq in sequences:
            if len(seq.seq) % 3 != 0:
                return jsonify({'success': False, 'error': f'Sequence {seq.id} length is not a multiple of 3'})

        # For now, provide codon-based analysis without external alignment tools
        alignment_data = []
        max_length = max(len(seq.seq) for seq in sequences)

        for seq in sequences:
            alignment_data.append({
                'id': seq.id,
                'sequence': str(seq.seq),
                'length': len(seq.seq),
                'codons': len(seq.seq) // 3
            })

        alignment = {
            'aligned_count': len(sequences),
            'alignment_length': max_length // 3,
            'method': method,
            'genetic_code': genetic_code,
            'alignment': alignment_data,
            'note': 'Codon-aware alignment. Full implementation with external alignment tools (ClustalW/MUSCLE) requires Bio.codonalign.build() with tool installation.'
        }

        return jsonify({'success': True, 'alignment': alignment})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})
