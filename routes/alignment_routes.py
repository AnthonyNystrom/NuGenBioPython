"""
Routes for sequence alignment operations
"""
from flask import Blueprint, request, jsonify, current_app
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio import AlignIO
from io import StringIO
import os
from werkzeug.utils import secure_filename
from collections import Counter

from dependencies import Seq

bp = Blueprint('alignment', __name__, url_prefix='/api')


@bp.route('/alignment/pairwise', methods=['POST'])
def pairwise_alignment():
    try:
        data = request.json
        seq1 = Seq.Seq(data.get('sequence1', ''))
        seq2 = Seq.Seq(data.get('sequence2', ''))

        aligner = PairwiseAligner()

        # Get alignment mode (default: global)
        mode = data.get('mode', 'global')
        if mode == 'local':
            aligner.mode = 'local'
        elif mode == 'global':
            aligner.mode = 'global'

        # Handle semi-global alignment
        if data.get('end_gap_score_query') is not None:
            aligner.query_end_gap_score = float(data.get('end_gap_score_query', 0))
        if data.get('end_gap_score_target') is not None:
            aligner.target_end_gap_score = float(data.get('end_gap_score_target', 0))

        # Get substitution matrix if specified
        matrix_name = data.get('substitution_matrix', None)
        if matrix_name and matrix_name != 'none':
            try:
                matrix = substitution_matrices.load(matrix_name)
                aligner.substitution_matrix = matrix
            except:
                # If matrix loading fails, use manual scoring
                aligner.match_score = float(data.get('match_score', 2))
                aligner.mismatch_score = float(data.get('mismatch_score', -1))
        else:
            # Use manual match/mismatch scores
            aligner.match_score = float(data.get('match_score', 2))
            aligner.mismatch_score = float(data.get('mismatch_score', -1))

        # Set gap penalties
        aligner.open_gap_score = float(data.get('gap_open', -2))
        aligner.extend_gap_score = float(data.get('gap_extend', -0.5))

        alignments = aligner.align(seq1, seq2)
        top_alignment = alignments[0]

        # Calculate alignment statistics
        alignment_str = str(top_alignment)
        lines = alignment_str.strip().split('\n')

        # Extract aligned sequences from the alignment string
        # BioPython format: "target  0 ACGT-CG 7"
        # We need to extract just the sequence part
        matches = 0
        gaps = 0
        identity_percent = 0
        gap_percent = 0
        aligned_seq1 = ""
        aligned_seq2 = ""

        if len(lines) >= 3:
            # Parse first line (target sequence)
            parts1 = lines[0].split()
            if len(parts1) >= 3:
                aligned_seq1 = parts1[-2]

            # Parse third line (query sequence)
            parts2 = lines[2].split()
            if len(parts2) >= 3:
                aligned_seq2 = parts2[-2]

            # Calculate identity, similarity, and gaps
            length = min(len(aligned_seq1), len(aligned_seq2))

            for i in range(length):
                if aligned_seq1[i] == aligned_seq2[i]:
                    matches += 1
                if aligned_seq1[i] == '-' or aligned_seq2[i] == '-':
                    gaps += 1

            identity_percent = (matches / length * 100) if length > 0 else 0
            gap_percent = (gaps / length * 100) if length > 0 else 0

        return jsonify({
            'success': True,
            'alignment': alignment_str,
            'score': float(top_alignment.score),
            'statistics': {
                'identity_percent': round(identity_percent, 2),
                'gap_percent': round(gap_percent, 2),
                'matches': matches,
                'gaps': gaps,
                'alignment_length': len(aligned_seq1) if len(lines) >= 2 else 0
            }
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/alignment/matrices', methods=['GET'])
def get_substitution_matrices():
    """Return list of available substitution matrices"""
    try:
        # Common matrices available in BioPython
        protein_matrices = [
            'BLOSUM45', 'BLOSUM50', 'BLOSUM62', 'BLOSUM80', 'BLOSUM90',
            'PAM30', 'PAM70', 'PAM250',
            'BENNER6', 'BENNER22', 'BENNER74',
            'GONNET1992', 'HOXD70', 'JOHNSON', 'JONES', 'LEVIN',
            'MCLACHLAN', 'MDM78', 'RAO', 'RISLER', 'SCHNEIDER', 'STR'
        ]

        dna_matrices = [
            'NUC.4.4', 'MATCH'
        ]

        all_matrices = protein_matrices + dna_matrices

        return jsonify({
            'success': True,
            'matrices': all_matrices,
            'protein_matrices': protein_matrices,
            'dna_matrices': dna_matrices
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/alignment/multiple', methods=['POST'])
def multiple_sequence_alignment():
    """Perform multiple sequence alignment using progressive alignment"""
    try:
        data = request.json
        sequences = data.get('sequences', [])

        if len(sequences) < 2:
            return jsonify({'success': False, 'error': 'At least 2 sequences required'})

        # Convert to Seq objects
        seq_objects = [Seq.Seq(seq) for seq in sequences]

        # Use pairwise alignment to build progressive MSA
        from Bio.Align import PairwiseAligner

        aligner = PairwiseAligner()
        aligner.mode = data.get('mode', 'global')
        aligner.match_score = float(data.get('match_score', 2))
        aligner.mismatch_score = float(data.get('mismatch_score', -1))
        aligner.open_gap_score = float(data.get('gap_open', -2))
        aligner.extend_gap_score = float(data.get('gap_extend', -0.5))

        # Simple progressive alignment: align first two, then add others
        alignments_list = []
        current_alignment = None

        for i, seq in enumerate(seq_objects):
            if i == 0:
                continue
            elif i == 1:
                # Align first two sequences
                alignments = aligner.align(seq_objects[0], seq)
                current_alignment = str(alignments[0])
                alignments_list.append({
                    'index': i,
                    'alignment': current_alignment,
                    'score': float(alignments[0].score)
                })
            else:
                # Align subsequent sequences to first sequence
                alignments = aligner.align(seq_objects[0], seq)
                alignments_list.append({
                    'index': i,
                    'alignment': str(alignments[0]),
                    'score': float(alignments[0].score)
                })

        # Build consensus and statistics
        total_score = sum(a['score'] for a in alignments_list)
        avg_score = total_score / len(alignments_list) if alignments_list else 0

        return jsonify({
            'success': True,
            'alignments': alignments_list,
            'total_score': round(total_score, 2),
            'average_score': round(avg_score, 2),
            'num_sequences': len(sequences)
        })

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})

# ============================================================================
# ADVANCED ALIGNMENT FEATURES
# ============================================================================

@bp.route('/alignment/consensus', methods=['POST'])
def generate_consensus():
    """Generate consensus sequence from alignment"""
    try:
        data = request.json
        sequences = data.get('sequences', [])

        if len(sequences) < 2:
            return jsonify({'success': False, 'error': 'At least 2 sequences required'})

        max_len = max(len(seq) for seq in sequences)
        padded_seqs = [seq + '-' * (max_len - len(seq)) for seq in sequences]

        consensus = []
        conservation_scores = []

        for pos in range(max_len):
            residues = [seq[pos] for seq in padded_seqs if pos < len(seq)]
            residues = [r for r in residues if r != '-']

            if residues:
                counter = Counter(residues)
                most_common = counter.most_common(1)[0]
                consensus_residue = most_common[0]
                frequency = most_common[1] / len(residues)

                consensus.append(consensus_residue)
                conservation_scores.append(round(frequency * 100, 1))
            else:
                consensus.append('-')
                conservation_scores.append(0)

        consensus_seq = ''.join(consensus)
        avg_conservation = round(sum(conservation_scores) / len(conservation_scores), 2) if conservation_scores else 0

        return jsonify({
            'success': True,
            'consensus': consensus_seq,
            'conservation_scores': conservation_scores,
            'average_conservation': avg_conservation
        })

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/alignment/conservation', methods=['POST'])
def analyze_conservation():
    """Analyze conservation across aligned sequences"""
    try:
        data = request.json
        sequences = data.get('sequences', [])

        if len(sequences) < 2:
            return jsonify({'success': False, 'error': 'At least 2 sequences required'})

        max_len = max(len(seq) for seq in sequences)
        padded_seqs = [seq + '-' * (max_len - len(seq)) for seq in sequences]

        conservation_data = []

        for pos in range(max_len):
            residues = [seq[pos] for seq in padded_seqs if pos < len(seq)]
            residues_no_gaps = [r for r in residues if r != '-']

            if residues_no_gaps:
                counter = Counter(residues_no_gaps)
                total = len(residues_no_gaps)

                max_freq = counter.most_common(1)[0][1]
                conservation_score = max_freq / total

                import math
                entropy = 0
                for count in counter.values():
                    p = count / total
                    entropy -= p * math.log2(p)

                conservation_data.append({
                    'position': pos + 1,
                    'conservation': round(conservation_score * 100, 1),
                    'entropy': round(entropy, 3),
                    'most_common': counter.most_common(1)[0][0],
                    'diversity': len(counter)
                })
            else:
                conservation_data.append({
                    'position': pos + 1,
                    'conservation': 0,
                    'entropy': 0,
                    'most_common': '-',
                    'diversity': 0
                })

        avg_conservation = round(sum(d['conservation'] for d in conservation_data) / len(conservation_data), 2)

        return jsonify({
            'success': True,
            'conservation_data': conservation_data,
            'average_conservation': avg_conservation,
            'length': len(conservation_data)
        })

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/alignment/parse_file', methods=['POST'])
def parse_alignment_file():
    """Parse alignment file in various formats"""
    try:
        file = request.files['file']
        file_format = request.form.get('format', 'fasta')

        if file:
            filename = secure_filename(file.filename)
            filepath = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
            file.save(filepath)

            alignment = AlignIO.read(filepath, file_format)

            sequences = []
            for record in alignment:
                sequences.append({
                    'id': record.id,
                    'description': record.description,
                    'sequence': str(record.seq),
                    'length': len(record.seq)
                })

            os.remove(filepath)

            return jsonify({
                'success': True,
                'sequences': sequences,
                'num_sequences': len(sequences),
                'alignment_length': len(alignment[0]) if len(alignment) > 0 else 0
            })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/alignment/export', methods=['POST'])
def export_alignment():
    """Export alignment in various formats"""
    try:
        data = request.json
        sequences = data.get('sequences', [])
        output_format = data.get('format', 'fasta')

        if not sequences:
            return jsonify({'success': False, 'error': 'No sequences provided'})

        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.Align import MultipleSeqAlignment

        records = []
        for i, seq_data in enumerate(sequences):
            seq_id = seq_data.get('id', f'seq{i+1}')
            sequence = seq_data.get('sequence', '')
            description = seq_data.get('description', '')

            record = SeqRecord(Seq(sequence), id=seq_id, description=description)
            records.append(record)

        alignment = MultipleSeqAlignment(records)

        output = StringIO()
        AlignIO.write(alignment, output, output_format)
        result = output.getvalue()

        return jsonify({
            'success': True,
            'exported': result,
            'format': output_format
        })

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/alignment/trim', methods=['POST'])
def trim_alignment():
    """Trim alignment by removing poorly aligned regions"""
    try:
        data = request.json
        sequences = data.get('sequences', [])
        min_conservation = data.get('min_conservation', 50)

        if len(sequences) < 2:
            return jsonify({'success': False, 'error': 'At least 2 sequences required'})

        max_len = max(len(seq) for seq in sequences)
        padded_seqs = [seq + '-' * (max_len - len(seq)) for seq in sequences]

        positions_to_keep = []

        for pos in range(max_len):
            residues = [seq[pos] for seq in padded_seqs if pos < len(seq)]
            residues_no_gaps = [r for r in residues if r != '-']

            if residues_no_gaps:
                counter = Counter(residues_no_gaps)
                max_freq = counter.most_common(1)[0][1]
                conservation = (max_freq / len(residues_no_gaps)) * 100

                if conservation >= min_conservation:
                    positions_to_keep.append(pos)

        trimmed_sequences = []
        for seq in padded_seqs:
            trimmed = ''.join([seq[pos] for pos in positions_to_keep])
            trimmed_sequences.append(trimmed)

        return jsonify({
            'success': True,
            'trimmed_sequences': trimmed_sequences,
            'original_length': max_len,
            'trimmed_length': len(positions_to_keep),
            'positions_removed': max_len - len(positions_to_keep)
        })

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})

# ============================================================================
# NEW ADVANCED FEATURES  
# ============================================================================

@bp.route('/alignment/all_alignments', methods=['POST'])
def get_all_alignments():
    """Generate all possible alignments (not just top 1)"""
    try:
        data = request.json
        seq1 = Seq.Seq(data.get('sequence1', ''))
        seq2 = Seq.Seq(data.get('sequence2', ''))
        max_alignments = data.get('max_alignments', 10)

        aligner = PairwiseAligner()
        aligner.mode = data.get('mode', 'global')
        aligner.match_score = float(data.get('match_score', 2))
        aligner.mismatch_score = float(data.get('mismatch_score', -1))
        aligner.open_gap_score = float(data.get('gap_open', -2))
        aligner.extend_gap_score = float(data.get('gap_extend', -0.5))

        alignments = aligner.align(seq1, seq2)

        alignment_list = []
        count = 0
        for alignment in alignments:
            if count >= max_alignments:
                break
            alignment_list.append({
                'alignment': str(alignment),
                'score': float(alignment.score),
                'index': count + 1
            })
            count += 1

        return jsonify({
            'success': True,
            'alignments': alignment_list,
            'total_found': count,
            'requested_max': max_alignments
        })

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/alignment/identity_matrix', methods=['POST'])
def pairwise_identity_matrix():
    """Calculate pairwise identity matrix for multiple sequences"""
    try:
        data = request.json
        sequences = data.get('sequences', [])

        if len(sequences) < 2:
            return jsonify({'success': False, 'error': 'At least 2 sequences required'})

        seq_objects = [Seq.Seq(seq) for seq in sequences]
        n = len(seq_objects)

        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 1
        aligner.mismatch_score = 0

        matrix = []
        labels = [f"Seq{i+1}" for i in range(n)]

        for i in range(n):
            row = []
            for j in range(n):
                if i == j:
                    row.append(100.0)
                elif i < j:
                    alignments = aligner.align(seq_objects[i], seq_objects[j])
                    alignment = alignments[0]
                    alignment_str = str(alignment)
                    lines = alignment_str.strip().split('\n')

                    if len(lines) >= 3:
                        aligned_seq1 = lines[0].split()[-2] if len(lines[0].split()) >= 3 else ""
                        aligned_seq2 = lines[2].split()[-2] if len(lines[2].split()) >= 3 else ""

                        matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != '-')
                        length = max(len(aligned_seq1), len(aligned_seq2))
                        identity = (matches / length * 100) if length > 0 else 0
                        row.append(round(identity, 2))
                    else:
                        row.append(0.0)
                else:
                    row.append(matrix[j][i])
            matrix.append(row)

        return jsonify({
            'success': True,
            'matrix': matrix,
            'labels': labels,
            'num_sequences': n
        })

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/alignment/coordinates', methods=['POST'])
def alignment_coordinates():
    """Get alignment coordinates and path information"""
    try:
        data = request.json
        seq1 = Seq.Seq(data.get('sequence1', ''))
        seq2 = Seq.Seq(data.get('sequence2', ''))

        aligner = PairwiseAligner()
        aligner.mode = data.get('mode', 'global')
        aligner.match_score = float(data.get('match_score', 2))
        aligner.mismatch_score = float(data.get('mismatch_score', -1))
        aligner.open_gap_score = float(data.get('gap_open', -2))
        aligner.extend_gap_score = float(data.get('gap_extend', -0.5))

        alignments = aligner.align(seq1, seq2)
        alignment = alignments[0]

        # Get alignment path/coordinates
        alignment_str = str(alignment)
        lines = alignment_str.strip().split('\n')

        coordinates = {
            'target_start': 0,
            'target_end': len(seq1),
            'query_start': 0,
            'query_end': len(seq2),
            'score': float(alignment.score)
        }

        # Extract coordinates from alignment string if available
        if len(lines) > 0:
            # Format: "target  start SEQUENCE end"
            parts = lines[0].split()
            if len(parts) >= 4:
                coordinates['target_start'] = int(parts[1])
                coordinates['target_end'] = int(parts[-1])

            if len(lines) >= 3:
                parts = lines[2].split()
                if len(parts) >= 4:
                    coordinates['query_start'] = int(parts[1])
                    coordinates['query_end'] = int(parts[-1])

        return jsonify({
            'success': True,
            'coordinates': coordinates,
            'alignment': alignment_str,
            'mode': aligner.mode
        })

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/alignment/codon_aware', methods=['POST'])
def codon_aware_alignment():
    """Perform codon-aware alignment for coding sequences"""
    try:
        data = request.json
        seq1_str = data.get('sequence1', '')
        seq2_str = data.get('sequence2', '')

        # Validate sequences are multiples of 3
        if len(seq1_str) % 3 != 0 or len(seq2_str) % 3 != 0:
            return jsonify({
                'success': False,
                'error': 'Sequences must be multiples of 3 for codon alignment'
            })

        seq1 = Seq.Seq(seq1_str)
        seq2 = Seq.Seq(seq2_str)

        # Align at codon level
        aligner = PairwiseAligner()
        aligner.mode = 'global'

        # Use higher penalties for gaps to preserve codon structure
        aligner.match_score = float(data.get('match_score', 3))
        aligner.mismatch_score = float(data.get('mismatch_score', -1))
        aligner.open_gap_score = float(data.get('gap_open', -5))  # Higher penalty
        aligner.extend_gap_score = float(data.get('gap_extend', -2))  # Higher penalty

        alignments = aligner.align(seq1, seq2)
        alignment = alignments[0]

        alignment_str = str(alignment)
        lines = alignment_str.strip().split('\n')

        # Calculate codon statistics
        if len(lines) >= 3:
            aligned_seq1 = lines[0].split()[-2] if len(lines[0].split()) >= 3 else ""
            aligned_seq2 = lines[2].split()[-2] if len(lines[2].split()) >= 3 else ""

            # Count synonymous vs non-synonymous changes
            codon_matches = 0
            codon_mismatches = 0

            # Group into codons
            for i in range(0, min(len(aligned_seq1), len(aligned_seq2)) - 2, 3):
                codon1 = aligned_seq1[i:i+3]
                codon2 = aligned_seq2[i:i+3]

                if '-' not in codon1 and '-' not in codon2:
                    if codon1 == codon2:
                        codon_matches += 1
                    else:
                        codon_mismatches += 1

            total_codons = codon_matches + codon_mismatches
            codon_identity = (codon_matches / total_codons * 100) if total_codons > 0 else 0

            return jsonify({
                'success': True,
                'alignment': alignment_str,
                'score': float(alignment.score),
                'codon_statistics': {
                    'codon_matches': codon_matches,
                    'codon_mismatches': codon_mismatches,
                    'codon_identity_percent': round(codon_identity, 2),
                    'total_codons': total_codons
                }
            })

        return jsonify({
            'success': True,
            'alignment': alignment_str,
            'score': float(alignment.score)
        })

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/alignment/detailed_stats', methods=['POST'])
def detailed_alignment_stats():
    """Get detailed statistics for an alignment"""
    try:
        data = request.json
        seq1 = Seq.Seq(data.get('sequence1', ''))
        seq2 = Seq.Seq(data.get('sequence2', ''))

        aligner = PairwiseAligner()
        aligner.mode = data.get('mode', 'global')
        aligner.match_score = float(data.get('match_score', 2))
        aligner.mismatch_score = float(data.get('mismatch_score', -1))
        aligner.open_gap_score = float(data.get('gap_open', -2))
        aligner.extend_gap_score = float(data.get('gap_extend', -0.5))

        alignments = aligner.align(seq1, seq2)
        alignment = alignments[0]

        alignment_str = str(alignment)
        lines = alignment_str.strip().split('\n')

        stats = {
            'score': float(alignment.score),
            'mode': aligner.mode,
            'seq1_length': len(seq1),
            'seq2_length': len(seq2)
        }

        if len(lines) >= 3:
            aligned_seq1 = lines[0].split()[-2] if len(lines[0].split()) >= 3 else ""
            aligned_seq2 = lines[2].split()[-2] if len(lines[2].split()) >= 3 else ""

            matches = 0
            mismatches = 0
            gaps_seq1 = 0
            gaps_seq2 = 0
            gap_opens = 0
            gap_extends = 0
            in_gap = False

            for i in range(min(len(aligned_seq1), len(aligned_seq2))):
                a, b = aligned_seq1[i], aligned_seq2[i]

                if a == '-':
                    gaps_seq1 += 1
                    if not in_gap:
                        gap_opens += 1
                        in_gap = True
                    else:
                        gap_extends += 1
                elif b == '-':
                    gaps_seq2 += 1
                    if not in_gap:
                        gap_opens += 1
                        in_gap = True
                    else:
                        gap_extends += 1
                else:
                    in_gap = False
                    if a == b:
                        matches += 1
                    else:
                        mismatches += 1

            alignment_length = max(len(aligned_seq1), len(aligned_seq2))
            total_gaps = gaps_seq1 + gaps_seq2

            stats.update({
                'alignment_length': alignment_length,
                'matches': matches,
                'mismatches': mismatches,
                'gaps_seq1': gaps_seq1,
                'gaps_seq2': gaps_seq2,
                'total_gaps': total_gaps,
                'gap_opens': gap_opens,
                'gap_extends': gap_extends,
                'identity_percent': round((matches / alignment_length * 100), 2) if alignment_length > 0 else 0,
                'similarity_percent': round(((matches) / (matches + mismatches) * 100), 2) if (matches + mismatches) > 0 else 0,
                'gap_percent': round((total_gaps / alignment_length * 100), 2) if alignment_length > 0 else 0
            })

        return jsonify({
            'success': True,
            'statistics': stats,
            'alignment': alignment_str
        })

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})
