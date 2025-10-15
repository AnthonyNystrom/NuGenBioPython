"""
Routes for sequence analysis and SeqIO operations
"""
from flask import Blueprint, request, jsonify, send_file, current_app
import os
from werkzeug.utils import secure_filename
from io import StringIO

from dependencies import Seq, SeqUtils, SeqIO, SeqRecord
import re

bp = Blueprint('sequence', __name__, url_prefix='/api')


@bp.route('/sequence/analyze', methods=['POST'])
def analyze_sequence():
    try:
        data = request.json
        sequence_str = data.get('sequence', '').upper()
        seq_type = data.get('type', 'dna')

        # Create sequence object based on type
        if seq_type == 'dna':
            seq = Seq.Seq(sequence_str)
            composition = dict(zip(['A', 'T', 'G', 'C'], [sequence_str.count(base) for base in 'ATGC']))
            gc_content = SeqUtils.gc_fraction(seq) * 100
            complement = str(seq.complement())
            reverse_complement = str(seq.reverse_complement())
            translation = str(seq.translate())
        elif seq_type == 'rna':
            # For RNA, convert U to T for BioPython operations, then convert back
            dna_sequence = sequence_str.replace('U', 'T')
            seq = Seq.Seq(dna_sequence)
            composition = dict(zip(['A', 'U', 'G', 'C'], [sequence_str.count(base) for base in 'AUGC']))
            # Calculate GC content manually for RNA
            gc_count = sequence_str.count('G') + sequence_str.count('C')
            gc_content = (gc_count / len(sequence_str)) * 100 if len(sequence_str) > 0 else 0
            # Convert complement back to RNA (T -> U)
            complement = str(seq.complement()).replace('T', 'U')
            reverse_complement = str(seq.reverse_complement()).replace('T', 'U')
            translation = str(seq.translate())
        else:  # protein
            # For proteins, we need to handle amino acids differently
            # Create a DNA sequence for BioPython operations (using a dummy sequence)
            dummy_dna = 'ATG' * (len(sequence_str) // 3) + 'ATG'[:len(sequence_str) % 3]
            seq = Seq.Seq(dummy_dna)
            # For proteins, count amino acids
            amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
            composition = dict(zip(list(amino_acids), [sequence_str.count(aa) for aa in amino_acids]))
            gc_content = None  # Not applicable for proteins
            complement = None  # Not applicable for proteins
            reverse_complement = None  # Not applicable for proteins
            translation = None  # Not applicable for proteins

        # Calculate molecular weight based on sequence type
        if seq_type == 'protein':
            # For proteins, calculate molecular weight manually
            # Amino acid molecular weights (Da)
            aa_weights = {
                'A': 89.1, 'C': 121.2, 'D': 133.1, 'E': 147.1, 'F': 165.2,
                'G': 75.1, 'H': 155.2, 'I': 131.2, 'K': 146.2, 'L': 131.2,
                'M': 149.2, 'N': 132.1, 'P': 115.1, 'Q': 146.2, 'R': 174.2,
                'S': 105.1, 'T': 119.1, 'V': 117.1, 'W': 204.2, 'Y': 181.2
            }
            mol_weight = sum(aa_weights.get(aa, 0) for aa in sequence_str)
        else:
            mol_weight = SeqUtils.molecular_weight(seq)

        analysis = {
            'length': len(sequence_str),  # Use original sequence length
            'composition': composition,
            'gc_content': gc_content,
            'molecular_weight': mol_weight,
            'complement': complement,
            'reverse_complement': reverse_complement,
            'translation': translation
        }

        return jsonify({'success': True, 'analysis': analysis})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/sequence/melting_temp', methods=['POST'])
def melting_temp():
    """Calculate melting temperature of DNA sequence"""
    try:
        data = request.json
        sequence = data.get('sequence', '').strip().upper()

        if not sequence:
            return jsonify({'success': False, 'error': 'Please provide a DNA sequence'})

        from Bio.Seq import Seq
        from Bio.SeqUtils import MeltingTemp as mt

        seq = Seq(sequence)

        # Calculate using different methods
        tm_wallace = mt.Tm_Wallace(seq)  # For short sequences
        tm_gc = mt.Tm_GC(seq)  # GC-based method

        result = {
            'tm_wallace': round(tm_wallace, 2),
            'tm_gc': round(tm_gc, 2),
            'length': len(seq),
            'gc_content': round((seq.count('G') + seq.count('C')) / len(seq) * 100, 2)
        }

        # Add NN method for longer sequences
        if len(seq) > 13:
            tm_nn = mt.Tm_NN(seq)
            result['tm_nn'] = round(tm_nn, 2)
        else:
            result['tm_nn'] = None

        return jsonify({'success': True, **result})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})

@bp.route('/sequence/find_orfs', methods=['POST'])
def find_orfs():
    """Find Open Reading Frames in sequence"""
    try:
        data = request.json
        sequence = data.get('sequence', '').strip().upper()
        min_length = data.get('min_length', 75)

        if not sequence:
            return jsonify({'success': False, 'error': 'Please provide a DNA sequence'})

        from Bio.Seq import Seq

        seq = Seq(sequence)
        orfs = []

        # Find ORFs in all 6 frames
        for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
            for frame in range(3):
                trans = str(nuc[frame:].translate())
                trans_len = len(trans)
                aa_start = 0

                while aa_start < trans_len:
                    aa_end = trans.find('*', aa_start)
                    if aa_end == -1:
                        aa_end = trans_len

                    orf_length_bp = (aa_end - aa_start) * 3
                    if orf_length_bp >= min_length:
                        orf_seq = trans[aa_start:aa_end]
                        if orf_seq:
                            nt_start = aa_start * 3 + frame
                            nt_end = aa_end * 3 + frame
                            orfs.append({
                                'frame': frame + 1 if strand == 1 else -(frame + 1),
                                'strand': strand,
                                'start': nt_start,
                                'end': nt_end,
                                'length': orf_length_bp,
                                'protein_length': len(orf_seq),
                                'protein_sequence': orf_seq
                            })

                    aa_start = aa_end + 1

        # Sort by length
        orfs.sort(key=lambda x: x['length'], reverse=True)

        return jsonify({'success': True, 'orfs': orfs[:20], 'total_found': len(orfs)})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})

@bp.route('/sequence/codon_usage_analysis', methods=['POST'])
def codon_usage():
    """Analyze codon usage in sequence"""
    try:
        data = request.json
        sequence = data.get('sequence', '').strip().upper()

        if not sequence:
            return jsonify({'success': False, 'error': 'Please provide a sequence'})

        from Bio.Seq import Seq
        from collections import Counter

        seq = Seq(sequence)

        # Get codons
        codons = [str(seq[i:i+3]) for i in range(0, len(seq) - 2, 3) if len(seq[i:i+3]) == 3]

        # Count codon occurrences
        codon_counts = Counter(codons)
        total_codons = len(codons)

        # Calculate frequencies and prepare top codons
        top_codons = []
        for codon, count in codon_counts.most_common(20):
            top_codons.append({
                'codon': codon,
                'count': count,
                'frequency': round(count / total_codons * 100, 2)
            })

        return jsonify({
            'success': True,
            'total_codons': total_codons,
            'unique_codons': len(codon_counts),
            'top_codons': top_codons
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})

@bp.route('/sequence/checksums', methods=['POST'])
def checksums():
    """Calculate sequence checksums"""
    try:
        data = request.json
        sequence = data.get('sequence', '').strip().upper()

        if not sequence:
            return jsonify({'success': False, 'error': 'Please provide a sequence'})

        from Bio.Seq import Seq
        from Bio.SeqUtils.CheckSum import seguid, crc32, crc64

        seq = Seq(sequence)
        
        return jsonify({
            'success': True,
            'seguid': seguid(seq),
            'crc32': crc32(seq),
            'crc64': crc64(seq)
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})

@bp.route('/sequence/protparam_advanced', methods=['POST'])
def protparam_advanced():
    """Advanced ProtParam analysis"""
    try:
        data = request.json
        sequence = data.get('sequence', '').strip().upper()

        if not sequence:
            return jsonify({'success': False, 'error': 'Please provide a protein sequence'})

        from Bio.SeqUtils.ProtParam import ProteinAnalysis

        analyzed_seq = ProteinAnalysis(sequence)
        
        return jsonify({
            'success': True,
            'molar_extinction': {
                'reduced': analyzed_seq.molar_extinction_coefficient()[0],
                'oxidized': analyzed_seq.molar_extinction_coefficient()[1]
            },
            'charge_at_pH7': round(analyzed_seq.charge_at_pH(7.0), 4),
            'flexibility': analyzed_seq.flexibility()[:50] if len(analyzed_seq.flexibility()) > 0 else []
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})

@bp.route('/sequence/molecular_weight_advanced', methods=['POST'])
def molecular_weight_advanced():
    """Calculate molecular weight with variations"""
    try:
        data = request.json
        sequence = data.get('sequence', '').strip().upper()
        seq_type = data.get('type', 'dna')

        if not sequence:
            return jsonify({'success': False, 'error': 'Please provide a sequence'})

        from Bio.Seq import Seq
        from Bio.SeqUtils import molecular_weight

        seq = Seq(sequence)

        molecular_weights = {}

        if seq_type == 'dna':
            molecular_weights['single_stranded'] = round(molecular_weight(seq, seq_type='DNA'), 2)
            molecular_weights['double_stranded'] = round(molecular_weight(seq, seq_type='DNA', double_stranded=True), 2)
        elif seq_type == 'rna':
            molecular_weights['rna'] = round(molecular_weight(seq, seq_type='RNA'), 2)
        else:
            molecular_weights['protein'] = round(molecular_weight(seq, seq_type='protein'), 2)

        return jsonify({'success': True, 'molecular_weights': molecular_weights})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})

@bp.route('/sequence/reference_data', methods=['GET'])
def reference_data():
    """Get reference data (codon tables, IUPAC codes)"""
    try:
        from Bio.Data import CodonTable, IUPACData

        # Get standard genetic code
        standard_table = CodonTable.unambiguous_dna_by_name["Standard"]

        # Format codon table for UI
        codon_table = {
            'forward_table': dict(standard_table.forward_table),
            'stop_codons': list(standard_table.stop_codons),
            'start_codons': list(standard_table.start_codons)
        }

        # IUPAC codes
        iupac_dna = dict(IUPACData.ambiguous_dna_values)
        iupac_rna = dict(IUPACData.ambiguous_rna_values)
        iupac_protein = dict(IUPACData.protein_letters_1to3)

        return jsonify({
            'success': True,
            'codon_table': codon_table,
            'iupac_dna': iupac_dna,
            'iupac_rna': iupac_rna,
            'iupac_protein': iupac_protein
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})

@bp.route('/sequence/transcribe', methods=['POST'])
def transcribe():
    """Transcribe DNA to RNA or back-transcribe RNA to DNA"""
    try:
        data = request.json
        sequence = data.get('sequence', '').strip().upper()
        operation = data.get('operation', 'transcribe')

        if not sequence:
            return jsonify({'success': False, 'error': 'Please provide a sequence'})

        from Bio.Seq import Seq

        seq = Seq(sequence)

        if operation == 'transcribe':
            # DNA to RNA
            result = str(seq.transcribe())
            operation_name = 'DNA → RNA Transcription'
        elif operation == 'back_transcribe':
            # RNA to DNA
            result = str(seq.back_transcribe())
            operation_name = 'RNA → DNA Back-transcription'
        else:
            return jsonify({'success': False, 'error': 'Invalid operation'})

        return jsonify({
            'success': True,
            'operation': operation_name,
            'original': sequence,
            'result': result
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})

@bp.route('/sequence/gc_analysis', methods=['POST'])
def gc_analysis():
    """Perform advanced GC analysis including GC skew and GC123"""
    try:
        data = request.json
        sequence = data.get('sequence', '').strip().upper()
        window = data.get('window', 100)

        if not sequence:
            return jsonify({'success': False, 'error': 'Please provide a sequence'})

        from Bio.Seq import Seq
        from Bio.SeqUtils import gc_fraction, GC123, GC_skew

        seq = Seq(sequence)

        gc_content = round(gc_fraction(seq) * 100, 2)

        # GC123 - GC content by codon position
        gc123 = None
        if len(sequence) >= 3:
            try:
                gc123 = GC123(seq)
            except:
                pass

        # GC Skew
        gc_skew_values = []
        if len(sequence) >= window:
            try:
                gc_skew_values = GC_skew(seq, window=window)
            except:
                pass

        return jsonify({
            'success': True,
            'gc_content': gc_content,
            'gc123': gc123,
            'gc_skew': gc_skew_values,
            'window_size': window
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})

@bp.route('/sequence/protein_convert', methods=['POST'])
def protein_convert():
    """Convert protein sequences between 1-letter and 3-letter codes"""
    try:
        data = request.json
        sequence = data.get('sequence', '').strip().upper()
        conversion = data.get('conversion', 'to_three')

        if not sequence:
            return jsonify({'success': False, 'error': 'Please provide a sequence'})

        from Bio.SeqUtils import seq1, seq3

        if conversion == 'to_three':
            result = seq3(sequence)
            conversion_name = '1-Letter → 3-Letter Conversion'
        elif conversion == 'to_one':
            result = seq1(sequence)
            conversion_name = '3-Letter → 1-Letter Conversion'
        else:
            return jsonify({'success': False, 'error': 'Invalid conversion type'})

        return jsonify({
            'success': True,
            'conversion': conversion_name,
            'original': sequence,
            'result': result
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})

@bp.route('/sequence/manipulate', methods=['POST'])
def manipulate():
    """Manipulate sequences (ungap, uppercase, lowercase)"""
    try:
        data = request.json
        sequence = data.get('sequence', '').strip()
        operation = data.get('operation', 'ungap')

        if not sequence:
            return jsonify({'success': False, 'error': 'Please provide a sequence'})

        from Bio.Seq import Seq

        seq = Seq(sequence)
        original_length = len(sequence)

        if operation == 'ungap':
            result = str(seq.ungap())
            operation_name = 'Remove Gaps'
        elif operation == 'upper':
            result = sequence.upper()
            operation_name = 'Convert to Uppercase'
        elif operation == 'lower':
            result = sequence.lower()
            operation_name = 'Convert to Lowercase'
        else:
            return jsonify({'success': False, 'error': 'Invalid operation'})

        return jsonify({
            'success': True,
            'operation': operation_name,
            'original': sequence,
            'result': result,
            'original_length': original_length,
            'result_length': len(result)
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})

@bp.route('/sequence/search', methods=['POST'])
def search():
    """Search for patterns in sequences"""
    try:
        data = request.json
        sequence = data.get('sequence', '').strip().upper()
        pattern = data.get('pattern', '').strip().upper()
        search_type = data.get('search_type', 'find')

        if not sequence or not pattern:
            return jsonify({'success': False, 'error': 'Please provide both sequence and pattern'})

        from Bio.Seq import Seq
        from Bio.SeqUtils import nt_search

        seq = Seq(sequence)

        result = {}
        result['pattern'] = pattern
        result['type'] = search_type

        if search_type == 'find':
            position = sequence.find(pattern)
            result['position'] = position
            result['found'] = position >= 0
        elif search_type == 'count':
            count = sequence.count(pattern)
            result['count'] = count
        elif search_type == 'count_overlap':
            # Count overlapping occurrences
            count = 0
            start = 0
            while True:
                pos = sequence.find(pattern, start)
                if pos >= 0:
                    count += 1
                    start = pos + 1
                else:
                    break
            result['count'] = count
        elif search_type == 'nt_search':
            search_result = nt_search(sequence, pattern)
            result['count'] = len(search_result) - 1  # First element is the pattern
            result['positions'] = search_result[1:] if len(search_result) > 1 else []
        else:
            return jsonify({'success': False, 'error': 'Invalid search type'})

        return jsonify({
            'success': True,
            'result': result
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})
