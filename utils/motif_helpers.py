"""
Helper functions for motif analysis using Bio.motifs
"""
from io import StringIO
import numpy as np
from Bio import motifs
from Bio.Seq import Seq
from scipy.stats import pearsonr


def create_motif_from_sequences(sequences):
    """
    Create motif from list of sequences

    Args:
        sequences: List of sequence strings

    Returns:
        Motif object
    """
    if not sequences:
        raise ValueError("No sequences provided")

    # Convert to Seq objects
    seq_objects = [Seq(seq.upper()) for seq in sequences]

    # Create motif
    m = motifs.create(seq_objects)
    return m


def parse_motif_from_file(file_content, file_format='jaspar'):
    """
    Parse motif from file content

    Args:
        file_content: String content of motif file
        file_format: Format (jaspar, meme, transfac, pfm)

    Returns:
        Motif object or list of motifs
    """
    handle = StringIO(file_content)

    if file_format == 'jaspar':
        m = motifs.read(handle, 'jaspar')
        return m
    elif file_format == 'meme':
        # MEME can contain multiple motifs
        motif_list = motifs.parse(handle, 'meme')
        return list(motif_list)
    elif file_format == 'transfac':
        motif_list = motifs.parse(handle, 'transfac')
        return list(motif_list)
    elif file_format == 'pfm':
        m = motifs.read(handle, 'pfm')
        return m
    else:
        raise ValueError(f"Unsupported format: {file_format}")


def get_motif_info(m, pseudocounts=0.5):
    """
    Extract comprehensive motif information

    Args:
        m: Motif object
        pseudocounts: Pseudocount value for normalization

    Returns:
        Dictionary with motif information
    """
    # Basic info
    info = {
        'consensus': str(m.consensus),
        'anticonsensus': str(m.anticonsensus) if hasattr(m, 'anticonsensus') else None,
        'degenerate_consensus': str(m.degenerate_consensus),
        'length': len(m),
        'name': m.name if hasattr(m, 'name') else None
    }

    # Get counts
    counts = m.counts
    info['num_sequences'] = counts.get('A', [0])[0] + counts.get('C', [0])[0] + \
                           counts.get('G', [0])[0] + counts.get('T', [0])[0]

    # PWM (Position Weight Matrix)
    pwm = counts.normalize(pseudocounts=pseudocounts)
    pwm_data = []
    for i in range(len(m)):
        row = {'position': i + 1}
        for base in ['A', 'C', 'G', 'T']:
            row[base] = float(pwm[base][i])
        pwm_data.append(row)
    info['pwm'] = pwm_data

    # PSSM (Position-Specific Scoring Matrix) with background
    try:
        pssm = pwm.log_odds()
        pssm_data = []
        for i in range(len(m)):
            row = {'position': i + 1}
            for base in ['A', 'C', 'G', 'T']:
                row[base] = float(pssm[base][i])
            pssm_data.append(row)
        info['pssm'] = pssm_data
    except:
        info['pssm'] = None

    # Information content
    try:
        ic = pwm.information_content()
        info['information_content'] = float(ic)
    except:
        info['information_content'] = None

    # Relative entropy at each position
    try:
        re_per_pos = []
        for i in range(len(m)):
            re = 0
            for base in ['A', 'C', 'G', 'T']:
                p = pwm[base][i]
                if p > 0:
                    re += p * np.log2(p / 0.25)  # Assuming uniform background
            re_per_pos.append({'position': i + 1, 'entropy': float(re)})
        info['relative_entropy'] = re_per_pos
    except:
        info['relative_entropy'] = None

    return info


def search_motif_advanced(m, sequence, threshold_type='abs', threshold_value=0, pseudocounts=0.5):
    """
    Advanced motif search using PSSM scoring

    Args:
        m: Motif object
        sequence: Target sequence string
        threshold_type: 'abs' for absolute score, 'rel' for relative (0-1)
        threshold_value: Threshold value
        pseudocounts: Pseudocount value

    Returns:
        List of match dictionaries
    """
    seq = Seq(sequence.upper())
    pwm = m.counts.normalize(pseudocounts=pseudocounts)

    try:
        pssm = pwm.log_odds()
    except:
        # Fallback to PWM if PSSM fails
        pssm = pwm

    matches = []

    # Calculate max and min possible scores
    max_score = sum(max(pssm[base][i] for base in ['A', 'C', 'G', 'T']) for i in range(len(m)))
    min_score = sum(min(pssm[base][i] for base in ['A', 'C', 'G', 'T']) for i in range(len(m)))

    # Search forward strand
    for pos, score in pssm.search(seq, threshold=min_score):
        if threshold_type == 'rel':
            # Normalize score to 0-1
            normalized_score = (score - min_score) / (max_score - min_score) if max_score != min_score else 0
            if normalized_score >= threshold_value:
                matches.append({
                    'position': int(pos + 1),
                    'score': float(score),
                    'normalized_score': float(normalized_score),
                    'sequence': str(seq[pos:pos + len(m)]),
                    'strand': '+'
                })
        else:
            if score >= threshold_value:
                matches.append({
                    'position': int(pos + 1),
                    'score': float(score),
                    'normalized_score': float((score - min_score) / (max_score - min_score)) if max_score != min_score else 0,
                    'sequence': str(seq[pos:pos + len(m)]),
                    'strand': '+'
                })

    # Search reverse complement
    rc_seq = seq.reverse_complement()
    for pos, score in pssm.search(rc_seq, threshold=min_score):
        if threshold_type == 'rel':
            normalized_score = (score - min_score) / (max_score - min_score) if max_score != min_score else 0
            if normalized_score >= threshold_value:
                matches.append({
                    'position': int(len(seq) - pos - len(m) + 1),
                    'score': float(score),
                    'normalized_score': float(normalized_score),
                    'sequence': str(rc_seq[pos:pos + len(m)]),
                    'strand': '-'
                })
        else:
            if score >= threshold_value:
                matches.append({
                    'position': int(len(seq) - pos - len(m) + 1),
                    'score': float(score),
                    'normalized_score': float((score - min_score) / (max_score - min_score)) if max_score != min_score else 0,
                    'sequence': str(rc_seq[pos:pos + len(m)]),
                    'strand': '-'
                })

    # Sort by position
    matches.sort(key=lambda x: x['position'])

    return matches


def compare_motifs(motif1, motif2, pseudocounts=0.5):
    """
    Compare two motifs using Pearson correlation

    Args:
        motif1: First motif object
        motif2: Second motif object
        pseudocounts: Pseudocount value

    Returns:
        Dictionary with comparison results
    """
    pwm1 = motif1.counts.normalize(pseudocounts=pseudocounts)
    pwm2 = motif2.counts.normalize(pseudocounts=pseudocounts)

    # If different lengths, compare only overlapping region
    min_len = min(len(motif1), len(motif2))

    # Flatten PWMs for correlation
    flat1 = []
    flat2 = []

    for i in range(min_len):
        for base in ['A', 'C', 'G', 'T']:
            flat1.append(pwm1[base][i])
            flat2.append(pwm2[base][i])

    # Calculate Pearson correlation
    correlation, p_value = pearsonr(flat1, flat2)

    # Calculate distance metrics
    euclidean_dist = np.sqrt(sum((a - b) ** 2 for a, b in zip(flat1, flat2)))

    return {
        'pearson_correlation': float(correlation),
        'p_value': float(p_value),
        'euclidean_distance': float(euclidean_dist),
        'compared_length': min_len,
        'motif1_length': len(motif1),
        'motif2_length': len(motif2)
    }


def export_motif(m, format_type='jaspar'):
    """
    Export motif to specified format

    Args:
        m: Motif object
        format_type: Export format (jaspar, meme, transfac, pfm)

    Returns:
        String representation of motif in specified format
    """
    output = StringIO()

    if format_type == 'jaspar':
        output.write(m.format('jaspar'))
    elif format_type == 'meme':
        # MEME format requires specific structure
        output.write(f"MEME version 4\n\n")
        output.write(f"ALPHABET= ACGT\n\n")
        output.write(f"strands: + -\n\n")
        output.write(f"MOTIF {m.name if hasattr(m, 'name') else 'motif'}\n")
        output.write(f"letter-probability matrix: alength= 4 w= {len(m)}\n")

        pwm = m.counts.normalize(pseudocounts=0.5)
        for i in range(len(m)):
            output.write(f"  {pwm['A'][i]:.6f}  {pwm['C'][i]:.6f}  {pwm['G'][i]:.6f}  {pwm['T'][i]:.6f}\n")
    elif format_type == 'transfac':
        output.write(m.format('transfac'))
    elif format_type == 'pfm':
        output.write(m.format('pfm'))
    else:
        raise ValueError(f"Unsupported export format: {format_type}")

    return output.getvalue()


def calculate_motif_statistics(m, sequence, pseudocounts=0.5):
    """
    Calculate statistical significance of motif occurrences

    Args:
        m: Motif object
        sequence: Target sequence
        pseudocounts: Pseudocount value

    Returns:
        Dictionary with statistical information
    """
    seq = Seq(sequence.upper())
    pwm = m.counts.normalize(pseudocounts=pseudocounts)

    try:
        pssm = pwm.log_odds()

        # Get score distribution
        scores = []
        for pos in range(len(seq) - len(m) + 1):
            subseq = seq[pos:pos + len(m)]
            try:
                score = pssm.calculate(subseq)
                scores.append(score)
            except:
                continue

        if not scores:
            return None

        return {
            'mean_score': float(np.mean(scores)),
            'std_score': float(np.std(scores)),
            'max_score': float(np.max(scores)),
            'min_score': float(np.min(scores)),
            'median_score': float(np.median(scores)),
            'num_positions': len(scores)
        }
    except:
        return None
