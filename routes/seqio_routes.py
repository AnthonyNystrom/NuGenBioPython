from flask import Blueprint, request, jsonify, send_file, current_app
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqUtils
from werkzeug.utils import secure_filename
from io import StringIO, BytesIO
import os
import re
import gzip
import tempfile

bp = Blueprint('seqio', __name__, url_prefix='/api')

# ============================================================================
# BASIC SEQIO OPERATIONS
# ============================================================================

@bp.route('/seqio/parse', methods=['POST'])
def parse_sequence_file():
    """Parse sequence files in multiple formats with compression support"""
    try:
        file = request.files['file']
        file_format = request.form.get('format', 'fasta')

        if not file:
            return jsonify({'success': False, 'error': 'No file provided'})

        filename = secure_filename(file.filename)
        filepath = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
        file.save(filepath)

        # Check if file is gzipped
        is_gzipped = filename.endswith('.gz')

        sequences = []
        try:
            if is_gzipped:
                with gzip.open(filepath, 'rt') as handle:
                    for record in SeqIO.parse(handle, file_format):
                        sequences.append({
                            'id': record.id,
                            'description': record.description,
                            'sequence': str(record.seq),
                            'length': len(record.seq)
                        })
            else:
                for record in SeqIO.parse(filepath, file_format):
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
                'format': file_format,
                'compressed': is_gzipped
            })
        except Exception as e:
            os.remove(filepath)
            return jsonify({'success': False, 'error': f'Failed to parse file: {str(e)}'})

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/seqio/convert', methods=['POST'])
def convert_sequence_format():
    """Convert sequences between different formats"""
    try:
        data = request.json
        sequences = data.get('sequences', [])
        output_format = data.get('output_format', 'fasta')

        if not sequences:
            return jsonify({'success': False, 'error': 'No sequences provided'})

        output = StringIO()
        seq_records = []

        for seq_data in sequences:
            record = SeqRecord(
                Seq(seq_data.get('sequence', '')),
                id=seq_data.get('id', 'seq'),
                description=seq_data.get('description', '')
            )
            # Add molecule_type annotation for formats that require it
            if output_format.lower() in ['genbank', 'gb', 'embl']:
                record.annotations['molecule_type'] = 'DNA'
            seq_records.append(record)

        SeqIO.write(seq_records, output, output_format)

        return jsonify({
            'success': True,
            'converted': output.getvalue(),
            'count': len(seq_records)
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/seqio/write', methods=['POST'])
def write_sequences():
    """Write sequences to downloadable file"""
    try:
        data = request.json
        sequences = data.get('sequences', [])
        output_format = data.get('output_format', 'fasta')
        filename = data.get('filename', f'sequences.{output_format}')

        if not sequences:
            return jsonify({'success': False, 'error': 'No sequences provided'})

        seq_records = []
        for seq_data in sequences:
            record = SeqRecord(
                Seq(seq_data.get('sequence', '')),
                id=seq_data.get('id', 'seq'),
                description=seq_data.get('description', '')
            )
            if output_format.lower() in ['genbank', 'gb', 'embl']:
                record.annotations['molecule_type'] = 'DNA'
            seq_records.append(record)

        # Create temporary file
        temp_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=f'.{output_format}')
        SeqIO.write(seq_records, temp_file.name, output_format)
        temp_file.close()

        return send_file(temp_file.name, as_attachment=True, download_name=filename)

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/seqio/to_dict', methods=['POST'])
def sequences_to_dict():
    """Convert sequences to dictionary indexed by ID"""
    try:
        data = request.json
        sequences = data.get('sequences', [])

        seq_dict = {}
        for seq_data in sequences:
            seq_id = seq_data.get('id', '')
            seq_dict[seq_id] = seq_data

        return jsonify({
            'success': True,
            'sequence_dict': seq_dict,
            'count': len(seq_dict),
            'ids': list(seq_dict.keys())
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


# ============================================================================
# FASTQ-SPECIFIC OPERATIONS
# ============================================================================

@bp.route('/seqio/parse_fastq', methods=['POST'])
def parse_fastq():
    """Parse FASTQ files with quality scores"""
    try:
        file = request.files['file']
        file_format = request.form.get('format', 'fastq')

        if not file:
            return jsonify({'success': False, 'error': 'No file provided'})

        filename = secure_filename(file.filename)
        filepath = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
        file.save(filepath)

        # Check if gzipped
        is_gzipped = filename.endswith('.gz')

        sequences = []
        try:
            if is_gzipped:
                with gzip.open(filepath, 'rt') as handle:
                    for record in SeqIO.parse(handle, file_format):
                        seq_dict = _extract_fastq_record(record)
                        sequences.append(seq_dict)
            else:
                for record in SeqIO.parse(filepath, file_format):
                    seq_dict = _extract_fastq_record(record)
                    sequences.append(seq_dict)

            os.remove(filepath)
            return jsonify({
                'success': True,
                'sequences': sequences,
                'format': file_format,
                'compressed': is_gzipped
            })
        except Exception as e:
            os.remove(filepath)
            return jsonify({'success': False, 'error': f'Failed to parse FASTQ: {str(e)}'})

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


def _extract_fastq_record(record):
    """Helper function to extract FASTQ record data"""
    seq_dict = {
        'id': record.id,
        'description': record.description,
        'sequence': str(record.seq),
        'length': len(record.seq)
    }

    # Add quality scores if available
    if hasattr(record, 'letter_annotations') and 'phred_quality' in record.letter_annotations:
        seq_dict['quality_scores'] = record.letter_annotations['phred_quality']
        seq_dict['avg_quality'] = round(sum(record.letter_annotations['phred_quality']) / len(record.letter_annotations['phred_quality']), 2)

    return seq_dict


@bp.route('/seqio/fastq_filter_quality', methods=['POST'])
def fastq_filter_quality():
    """Filter FASTQ sequences by quality thresholds"""
    try:
        # Check if file upload or JSON data
        if 'file' in request.files:
            file = request.files['file']
            file_format = request.form.get('format', 'fastq')
            min_quality = int(request.form.get('min_quality', 20))
            min_percent = float(request.form.get('min_percent', 90))

            filename = secure_filename(file.filename)
            filepath = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
            file.save(filepath)

            # Parse FASTQ and filter
            filtered = []
            total = 0
            try:
                for record in SeqIO.parse(filepath, file_format):
                    total += 1
                    quality_scores = record.letter_annotations.get('phred_quality', [])
                    if quality_scores:
                        high_quality_bases = sum(1 for q in quality_scores if q >= min_quality)
                        percent_high_quality = (high_quality_bases / len(quality_scores) * 100)
                        if percent_high_quality >= min_percent:
                            filtered.append(_extract_fastq_record(record))
                os.remove(filepath)
            except:
                if os.path.exists(filepath):
                    os.remove(filepath)
                raise
        else:
            # JSON data
            data = request.json
            sequences = data.get('sequences', [])
            min_quality = data.get('min_quality', 20)
            min_percent = data.get('min_percent', 90)
            total = len(sequences)

            filtered = []
            for seq_data in sequences:
                quality_scores = seq_data.get('quality_scores', [])
                if not quality_scores:
                    continue
                high_quality_bases = sum(1 for q in quality_scores if q >= min_quality)
                percent_high_quality = (high_quality_bases / len(quality_scores) * 100) if quality_scores else 0
                if percent_high_quality >= min_percent:
                    seq_data['percent_high_quality'] = round(percent_high_quality, 2)
                    filtered.append(seq_data)

        return jsonify({
            'success': True,
            'filtered': filtered,
            'total_sequences': total,
            'passed_filter': len(filtered),
            'threshold': f'Q{min_quality}',
            'min_percent': min_percent
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/seqio/fastq_trim_quality', methods=['POST'])
def fastq_trim_quality():
    """Trim sequences based on quality scores"""
    try:
        # Check if file upload or JSON data
        if 'file' in request.files:
            file = request.files['file']
            file_format = request.form.get('format', 'fastq')
            quality_cutoff = int(request.form.get('quality_cutoff', 20))
            trim_from = request.form.get('trim_from', 'both')

            filename = secure_filename(file.filename)
            filepath = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
            file.save(filepath)

            trimmed = []
            try:
                for record in SeqIO.parse(filepath, file_format):
                    quality_scores = record.letter_annotations.get('phred_quality', [])
                    sequence = str(record.seq)

                    if not quality_scores or not sequence:
                        continue

                    start_trim = 0
                    end_trim = len(sequence)

                    if trim_from in ['start', 'both']:
                        for i, q in enumerate(quality_scores):
                            if q >= quality_cutoff:
                                start_trim = i
                                break

                    if trim_from in ['end', 'both']:
                        for i in range(len(quality_scores) - 1, -1, -1):
                            if quality_scores[i] >= quality_cutoff:
                                end_trim = i + 1
                                break

                    trimmed_seq = sequence[start_trim:end_trim]
                    trimmed_qual = quality_scores[start_trim:end_trim]

                    if trimmed_seq:
                        trimmed.append({
                            'id': record.id,
                            'description': record.description,
                            'sequence': trimmed_seq,
                            'length': len(trimmed_seq),
                            'original_length': len(sequence),
                            'trimmed_start': start_trim,
                            'trimmed_end': len(sequence) - end_trim,
                            'quality_scores': trimmed_qual,
                            'avg_quality': round(sum(trimmed_qual) / len(trimmed_qual), 2) if trimmed_qual else 0
                        })
                os.remove(filepath)
            except:
                if os.path.exists(filepath):
                    os.remove(filepath)
                raise
        else:
            # JSON data
            data = request.json
            sequences = data.get('sequences', [])
            quality_cutoff = data.get('quality_cutoff', 20)
            trim_from = data.get('trim_from', 'both')

            trimmed = []
            for seq_data in sequences:
                quality_scores = seq_data.get('quality_scores', [])
                sequence = seq_data.get('sequence', '')

                if not quality_scores or not sequence:
                    continue

                start_trim = 0
                end_trim = len(sequence)

                if trim_from in ['start', 'both']:
                    for i, q in enumerate(quality_scores):
                        if q >= quality_cutoff:
                            start_trim = i
                            break

                if trim_from in ['end', 'both']:
                    for i in range(len(quality_scores) - 1, -1, -1):
                        if quality_scores[i] >= quality_cutoff:
                            end_trim = i + 1
                            break

                trimmed_seq = sequence[start_trim:end_trim]
                trimmed_qual = quality_scores[start_trim:end_trim]

                if trimmed_seq:
                    trimmed.append({
                        'id': seq_data.get('id', ''),
                        'description': seq_data.get('description', ''),
                        'sequence': trimmed_seq,
                        'length': len(trimmed_seq),
                        'original_length': len(sequence),
                        'trimmed_start': start_trim,
                        'trimmed_end': len(sequence) - end_trim,
                        'quality_scores': trimmed_qual,
                        'avg_quality': round(sum(trimmed_qual) / len(trimmed_qual), 2) if trimmed_qual else 0
                    })

        return jsonify({
            'success': True,
            'sequences': trimmed,
            'original_count': len(sequences) if 'file' not in request.files else len(trimmed),
            'trimmed_count': len(trimmed),
            'quality_cutoff': quality_cutoff
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


# ============================================================================
# FEATURE EXTRACTION
# ============================================================================

@bp.route('/seqio/extract_features', methods=['POST'])
def extract_features():
    """Extract features from GenBank/EMBL files"""
    try:
        file = request.files['file']
        file_format = request.form.get('format', 'genbank')
        feature_type = request.form.get('feature_type', 'all')  # 'all', 'CDS', 'gene', etc.

        if not file:
            return jsonify({'success': False, 'error': 'No file provided'})

        filename = secure_filename(file.filename)
        filepath = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
        file.save(filepath)

        extracted_features = []
        try:
            for record in SeqIO.parse(filepath, file_format):
                for feature in record.features:
                    if feature_type == 'all' or feature.type == feature_type:
                        feature_data = {
                            'type': feature.type,
                            'location': str(feature.location),
                            'strand': feature.location.strand,
                            'qualifiers': dict(feature.qualifiers)
                        }

                        # Extract sequence if possible
                        try:
                            feature_seq = feature.extract(record.seq)
                            feature_data['sequence'] = str(feature_seq)
                            feature_data['length'] = len(feature_seq)
                        except:
                            feature_data['sequence'] = None

                        extracted_features.append(feature_data)

            os.remove(filepath)
            return jsonify({
                'success': True,
                'features': extracted_features,
                'count': len(extracted_features),
                'format': file_format
            })
        except Exception as e:
            os.remove(filepath)
            return jsonify({'success': False, 'error': f'Failed to extract features: {str(e)}'})

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


# ============================================================================
# FILTERING & SORTING
# ============================================================================

@bp.route('/seqio/filter', methods=['POST'])
def filter_sequences():
    """Filter sequences by length, ID pattern, or GC content"""
    try:
        data = request.json
        sequences = data.get('sequences', [])
        filter_type = data.get('filter_type', 'length')
        filter_params = data.get('filter_params', {})

        # Support both nested filter_params and direct parameters
        min_length = data.get('min_length', filter_params.get('min_length', 0))
        max_length = data.get('max_length', filter_params.get('max_length', float('inf')))
        pattern = data.get('pattern', filter_params.get('pattern', ''))
        min_gc = data.get('min_gc', filter_params.get('min_gc', 0))
        max_gc = data.get('max_gc', filter_params.get('max_gc', 100))

        filtered = []

        for seq_data in sequences:
            seq_id = seq_data.get('id', '')
            sequence = seq_data.get('sequence', '')
            length = len(sequence)

            include = False

            if filter_type == 'length':
                include = min_length <= length <= max_length

            elif filter_type == 'id_pattern':
                if pattern:
                    include = bool(re.search(pattern, seq_id, re.IGNORECASE))
                else:
                    include = True

            elif filter_type == 'gc_content':
                seq_obj = Seq(sequence)
                gc = SeqUtils.gc_fraction(seq_obj) * 100
                include = min_gc <= gc <= max_gc

            if include:
                filtered.append(seq_data)

        return jsonify({
            'success': True,
            'filtered': filtered,
            'original_count': len(sequences),
            'filtered_count': len(filtered)
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/seqio/sort', methods=['POST'])
def sort_sequences():
    """Sort sequences by length, ID, or GC content"""
    try:
        data = request.json
        sequences = data.get('sequences', [])
        sort_by = data.get('sort_by', 'length')
        reverse = data.get('reverse', False)

        if sort_by == 'length':
            sorted_seqs = sorted(sequences, key=lambda x: len(x.get('sequence', '')), reverse=reverse)
        elif sort_by == 'id':
            sorted_seqs = sorted(sequences, key=lambda x: x.get('id', ''), reverse=reverse)
        elif sort_by == 'gc_content':
            def get_gc(seq_data):
                seq = seq_data.get('sequence', '')
                if seq:
                    seq_obj = Seq(seq)
                    return SeqUtils.gc_fraction(seq_obj) * 100
                return 0
            sorted_seqs = sorted(sequences, key=get_gc, reverse=reverse)
        else:
            sorted_seqs = sequences

        return jsonify({
            'success': True,
            'sorted': sorted_seqs,
            'sort_by': sort_by,
            'reverse': reverse
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


# ============================================================================
# SEQUENCE MANIPULATION
# ============================================================================

@bp.route('/seqio/slice', methods=['POST'])
def slice_sequences():
    """Extract subsequences or slice sequence collection"""
    try:
        data = request.json
        sequences = data.get('sequences', [])
        start = data.get('start', 0)
        end = data.get('end', None)
        slice_type = data.get('slice_type', 'sequence')  # 'sequence' or 'collection'

        if slice_type == 'collection':
            # Slice the collection of sequences
            if end is None:
                sliced = sequences[start:]
            else:
                sliced = sequences[start:end]
        else:
            # Slice each individual sequence (default behavior)
            sliced = []
            for seq_data in sequences:
                original_seq = seq_data.get('sequence', '')
                if end is None:
                    sliced_seq = original_seq[start:]
                else:
                    sliced_seq = original_seq[start:end]

                sliced.append({
                    'id': seq_data.get('id', ''),
                    'description': seq_data.get('description', ''),
                    'sequence': sliced_seq,
                    'length': len(sliced_seq),
                    'original_length': len(original_seq)
                })

        return jsonify({
            'success': True,
            'sliced': sliced,
            'start': start,
            'end': end,
            'slice_type': slice_type
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


# ============================================================================
# STATISTICS & ANALYSIS
# ============================================================================

@bp.route('/seqio/statistics', methods=['POST'])
def sequence_statistics():
    """Calculate comprehensive statistics for sequences"""
    try:
        data = request.json
        sequences = data.get('sequences', [])

        stats = {
            'total_sequences': len(sequences),
            'total_length': 0,
            'avg_length': 0,
            'min_length': float('inf'),
            'max_length': 0,
            'gc_content_avg': 0,
            'composition': {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        }

        gc_values = []

        for seq_data in sequences:
            seq = seq_data.get('sequence', '')
            length = len(seq)

            stats['total_length'] += length
            stats['min_length'] = min(stats['min_length'], length)
            stats['max_length'] = max(stats['max_length'], length)

            # Composition
            stats['composition']['A'] += seq.count('A')
            stats['composition']['T'] += seq.count('T')
            stats['composition']['G'] += seq.count('G')
            stats['composition']['C'] += seq.count('C')

            # GC content
            if seq:
                seq_obj = Seq(seq)
                gc = SeqUtils.gc_fraction(seq_obj) * 100
                gc_values.append(gc)

        if sequences:
            stats['avg_length'] = round(stats['total_length'] / len(sequences), 2)
            stats['gc_content_avg'] = round(sum(gc_values) / len(gc_values), 2) if gc_values else 0

        if stats['min_length'] == float('inf'):
            stats['min_length'] = 0

        return jsonify({'success': True, 'statistics': stats})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


# ============================================================================
# BATCH PROCESSING
# ============================================================================

@bp.route('/seqio/batch_convert', methods=['POST'])
def batch_convert():
    """Convert multiple files at once"""
    try:
        files = request.files.getlist('files')
        input_format = request.form.get('input_format', 'fasta')
        output_format = request.form.get('output_format', 'genbank')

        if not files:
            return jsonify({'success': False, 'error': 'No files provided'})

        results = []
        for file in files:
            filename = secure_filename(file.filename)
            filepath = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
            file.save(filepath)

            try:
                # Read sequences
                records = list(SeqIO.parse(filepath, input_format))

                # Add annotations if needed
                for record in records:
                    if output_format.lower() in ['genbank', 'gb', 'embl']:
                        record.annotations['molecule_type'] = 'DNA'

                # Convert
                output = StringIO()
                SeqIO.write(records, output, output_format)

                results.append({
                    'filename': filename,
                    'success': True,
                    'record_count': len(records),
                    'converted': output.getvalue()
                })

                os.remove(filepath)

            except Exception as e:
                os.remove(filepath)
                results.append({
                    'filename': filename,
                    'success': False,
                    'error': str(e)
                })

        return jsonify({
            'success': True,
            'results': results,
            'total_files': len(files),
            'successful': sum(1 for r in results if r.get('success'))
        })

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})
