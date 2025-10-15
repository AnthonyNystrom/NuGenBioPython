"""
Routes for sequence feature operations (ORF, CDS, annotations)
"""
from flask import Blueprint, request, jsonify, current_app
import os
from werkzeug.utils import secure_filename
import re

from dependencies import Seq, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation, ExactPosition, BeforePosition, AfterPosition

bp = Blueprint('features', __name__, url_prefix='/api')


@bp.route('/features/orf_find', methods=['POST'])
def find_orfs():
    """Find Open Reading Frames in a sequence"""
    try:
        data = request.get_json()
        sequence = data.get('sequence', '').upper().strip()
        min_length = int(data.get('min_length', 100))  # Minimum ORF length in nucleotides
        strand = data.get('strand', 'both')  # 'forward', 'reverse', 'both'

        if not sequence:
            return jsonify({'success': False, 'error': 'No sequence provided'})

        seq_obj = Seq.Seq(sequence)
        orfs = []

        # Start and stop codons
        start_codons = ['ATG']
        stop_codons = ['TAA', 'TAG', 'TGA']

        def find_orfs_in_sequence(seq, frame_offset, strand_name):
            """Find ORFs in a given sequence and frame"""
            seq_str = str(seq)
            orf_list = []

            for frame in range(3):
                start_pos = frame + frame_offset
                i = start_pos

                while i < len(seq_str) - 2:
                    codon = seq_str[i:i+3]

                    if codon in start_codons:
                        # Found start codon, look for stop
                        start = i
                        j = i + 3

                        while j < len(seq_str) - 2:
                            stop_codon = seq_str[j:j+3]
                            if stop_codon in stop_codons:
                                # Found stop codon
                                end = j + 3
                                orf_length = end - start

                                if orf_length >= min_length:
                                    orf_seq = seq_str[start:end]
                                    protein = str(Seq.Seq(orf_seq).translate())

                                    orf_list.append({
                                        'start': start + 1,  # 1-indexed
                                        'end': end,
                                        'length': orf_length,
                                        'frame': frame + 1,
                                        'strand': strand_name,
                                        'sequence': orf_seq,
                                        'protein': protein,
                                        'protein_length': len(protein) - 1  # Exclude stop codon
                                    })

                                i = j + 3  # Move past this ORF
                                break
                            j += 3
                        else:
                            # No stop codon found, reached end of sequence
                            i += 3
                    else:
                        i += 3

            return orf_list

        # Find ORFs on forward strand
        if strand in ['forward', 'both']:
            orfs.extend(find_orfs_in_sequence(seq_obj, 0, 'forward'))

        # Find ORFs on reverse strand
        if strand in ['reverse', 'both']:
            rev_comp = seq_obj.reverse_complement()
            reverse_orfs = find_orfs_in_sequence(rev_comp, 0, 'reverse')
            # Adjust coordinates for reverse strand
            for orf in reverse_orfs:
                orig_start = len(sequence) - orf['end'] + 1
                orig_end = len(sequence) - orf['start'] + 1
                orf['start'] = orig_start
                orf['end'] = orig_end
            orfs.extend(reverse_orfs)

        # Sort by length (longest first)
        orfs.sort(key=lambda x: x['length'], reverse=True)

        return jsonify({
            'success': True,
            'orf_count': len(orfs),
            'orfs': orfs[:50],  # Limit to 50 for display
            'min_length': min_length,
            'sequence_length': len(sequence)
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/features/create', methods=['POST'])
def create_feature():
    """Create a SeqFeature and add it to a sequence"""
    try:
        data = request.get_json()
        sequence = data.get('sequence', '').upper().strip()
        feature_type = data.get('feature_type', 'CDS')
        start = int(data.get('start', 0))
        end = int(data.get('end', 0))
        strand = int(data.get('strand', 1))  # 1 for forward, -1 for reverse
        qualifiers = data.get('qualifiers', {})

        if not sequence:
            return jsonify({'success': False, 'error': 'No sequence provided'})

        # Create SeqRecord
        seq_record = SeqRecord(
            Seq.Seq(sequence),
            id=data.get('seq_id', 'sequence'),
            description=data.get('description', '')
        )

        # Create feature location
        location = FeatureLocation(start, end, strand=strand)

        # Create feature
        feature = SeqFeature(location, type=feature_type, qualifiers=qualifiers)

        # Add feature to record
        seq_record.features.append(feature)

        # Extract feature sequence
        feature_seq = str(feature.extract(seq_record.seq))

        return jsonify({
            'success': True,
            'feature': {
                'type': feature_type,
                'start': start,
                'end': end,
                'strand': '+' if strand == 1 else '-',
                'length': end - start,
                'sequence': feature_seq,
                'qualifiers': qualifiers
            }
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/features/parse_genbank', methods=['POST'])
def parse_genbank_features():
    """Parse features from GenBank or EMBL file"""
    try:
        file = request.files.get('file')
        file_format = request.form.get('format', 'genbank')

        if not file:
            return jsonify({'success': False, 'error': 'No file provided'})

        filename = secure_filename(file.filename)
        filepath = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
        file.save(filepath)

        # Parse file
        record = SeqIO.read(filepath, file_format)

        # Extract features
        features_list = []
        for feature in record.features:
            feature_data = {
                'type': feature.type,
                'start': int(feature.location.start),
                'end': int(feature.location.end),
                'strand': feature.location.strand,
                'strand_symbol': '+' if feature.location.strand == 1 else '-' if feature.location.strand == -1 else '?',
                'length': len(feature),
                'qualifiers': {}
            }

            # Extract qualifiers
            for key, value in feature.qualifiers.items():
                if isinstance(value, list):
                    feature_data['qualifiers'][key] = value[0] if len(value) == 1 else value
                else:
                    feature_data['qualifiers'][key] = value

            # Extract feature sequence if possible
            try:
                feature_seq = str(feature.extract(record.seq))
                feature_data['sequence'] = feature_seq[:100] + ('...' if len(feature_seq) > 100 else '')
                feature_data['full_length'] = len(feature_seq)
            except:
                feature_data['sequence'] = None

            features_list.append(feature_data)

        os.remove(filepath)

        return jsonify({
            'success': True,
            'record_id': record.id,
            'record_description': record.description,
            'sequence_length': len(record.seq),
            'feature_count': len(features_list),
            'features': features_list
        })
    except Exception as e:
        if os.path.exists(filepath):
            os.remove(filepath)
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/features/extract', methods=['POST'])
def extract_feature():
    """Extract and optionally translate a feature from a sequence"""
    try:
        data = request.get_json()
        sequence = data.get('sequence', '').upper().strip()
        start = int(data.get('start', 0))
        end = int(data.get('end', 0))
        strand = int(data.get('strand', 1))
        translate = data.get('translate', False)
        feature_type = data.get('feature_type', 'gene')

        if not sequence:
            return jsonify({'success': False, 'error': 'No sequence provided'})

        # Create SeqRecord
        seq_record = SeqRecord(Seq.Seq(sequence), id='sequence')

        # Create feature location
        location = FeatureLocation(start, end, strand=strand)
        feature = SeqFeature(location, type=feature_type)

        # Extract feature sequence
        feature_seq = feature.extract(seq_record.seq)
        feature_str = str(feature_seq)

        result = {
            'start': start,
            'end': end,
            'strand': '+' if strand == 1 else '-',
            'length': end - start,
            'sequence': feature_str,
            'type': feature_type
        }

        # Translate if requested
        if translate and feature_type in ['CDS', 'gene']:
            try:
                protein = str(feature_seq.translate())
                result['protein'] = protein
                result['protein_length'] = len(protein)
            except Exception as e:
                result['translation_error'] = str(e)

        return jsonify({
            'success': True,
            'feature': result
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/features/compound_location', methods=['POST'])
def create_compound_location():
    """Create a compound location (for split features like introns/exons)"""
    try:
        data = request.get_json()
        sequence = data.get('sequence', '').upper().strip()
        locations = data.get('locations', [])  # List of [start, end, strand]

        if not sequence or not locations:
            return jsonify({'success': False, 'error': 'Sequence and locations required'})

        # Create individual FeatureLocations
        feature_locations = []
        for loc in locations:
            start = int(loc['start'])
            end = int(loc['end'])
            strand = int(loc.get('strand', 1))
            feature_locations.append(FeatureLocation(start, end, strand=strand))

        # Create compound location
        compound_loc = CompoundLocation(feature_locations)

        # Create SeqRecord and feature
        seq_record = SeqRecord(Seq.Seq(sequence), id='sequence')
        feature = SeqFeature(compound_loc, type='mRNA')

        # Extract sequence (will join the parts)
        feature_seq = str(feature.extract(seq_record.seq))

        # Calculate total length
        total_length = sum(int(loc['end']) - int(loc['start']) for loc in locations)

        return jsonify({
            'success': True,
            'compound_feature': {
                'parts': len(locations),
                'total_length': total_length,
                'sequence': feature_seq,
                'locations': [{'start': loc['start'], 'end': loc['end'], 'strand': loc.get('strand', 1)} for loc in locations]
            }
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/features/annotate', methods=['POST'])
def annotate_sequence():
    """Add annotations to a sequence and export as GenBank"""
    try:
        data = request.get_json()
        sequence = data.get('sequence', '').upper().strip()
        seq_id = data.get('seq_id', 'sequence')
        description = data.get('description', 'Annotated sequence')
        features_data = data.get('features', [])

        if not sequence:
            return jsonify({'success': False, 'error': 'No sequence provided'})

        # Create SeqRecord
        seq_record = SeqRecord(
            Seq.Seq(sequence),
            id=seq_id,
            description=description
        )
        seq_record.annotations['molecule_type'] = 'DNA'

        # Add features
        for feat_data in features_data:
            start = int(feat_data.get('start', 0))
            end = int(feat_data.get('end', 0))
            strand = int(feat_data.get('strand', 1))
            feature_type = feat_data.get('type', 'misc_feature')
            qualifiers = feat_data.get('qualifiers', {})

            location = FeatureLocation(start, end, strand=strand)
            feature = SeqFeature(location, type=feature_type, qualifiers=qualifiers)
            seq_record.features.append(feature)

        # Convert to GenBank format
        from io import StringIO
        output = StringIO()
        SeqIO.write(seq_record, output, 'genbank')
        genbank_str = output.getvalue()

        return jsonify({
            'success': True,
            'feature_count': len(seq_record.features),
            'genbank': genbank_str,
            'summary': {
                'id': seq_id,
                'length': len(sequence),
                'features': len(seq_record.features)
            }
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/features/feature_types', methods=['GET'])
def get_feature_types():
    """Return list of standard feature types"""
    feature_types = [
        'gene', 'CDS', 'mRNA', 'tRNA', 'rRNA',
        'exon', 'intron', 'promoter', 'terminator',
        'regulatory', 'enhancer', 'RBS', 'primer_bind',
        'misc_feature', 'misc_RNA', 'repeat_region',
        'source', 'variation', 'stem_loop', 'UTR',
        '5\'UTR', '3\'UTR', 'polyA_signal', 'polyA_site'
    ]

    return jsonify({
        'success': True,
        'feature_types': feature_types
    })
