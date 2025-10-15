"""
Helper functions for BLAST operations
"""
from io import StringIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def extract_sequence_from_file(file_content, file_format='fasta'):
    """
    Extract sequence from uploaded file

    Args:
        file_content: File content as string
        file_format: Format (fasta, genbank, or txt)

    Returns:
        Extracted sequence string
    """
    try:
        if file_format == 'txt':
            # Plain text - clean and return
            sequence = file_content.strip().replace(' ', '').replace('\n', '')
            return sequence

        elif file_format in ['fasta', 'fa', 'fna', 'faa']:
            # Parse FASTA format
            handle = StringIO(file_content)
            sequences = []
            for record in SeqIO.parse(handle, 'fasta'):
                sequences.append(str(record.seq))
            return ''.join(sequences)

        elif file_format in ['genbank', 'gb', 'gbk']:
            # Parse GenBank format
            handle = StringIO(file_content)
            sequences = []
            for record in SeqIO.parse(handle, 'genbank'):
                sequences.append(str(record.seq))
            return ''.join(sequences)

        else:
            raise ValueError(f"Unsupported file format: {file_format}")

    except Exception as e:
        raise ValueError(f"Error parsing file: {str(e)}")


def run_ncbi_blast(sequence, program='blastn', database='nt', **kwargs):
    """
    Run BLAST search against NCBI servers

    Args:
        sequence: Query sequence string
        program: BLAST program (blastn, blastp, blastx, tblastn, tblastx)
        database: Target database
        **kwargs: Additional BLAST parameters

    Returns:
        BLAST results in XML format
    """
    # Clean sequence
    sequence = sequence.strip().replace(' ', '').replace('\n', '')

    # Run BLAST
    result_handle = NCBIWWW.qblast(
        program=program,
        database=database,
        sequence=sequence,
        **kwargs
    )

    return result_handle.read()


def parse_blast_xml(blast_xml):
    """
    Parse BLAST XML results

    Args:
        blast_xml: BLAST XML string

    Returns:
        List of parsed BLAST records
    """
    blast_records = []

    with StringIO(blast_xml) as result_handle:
        for blast_record in NCBIXML.parse(result_handle):
            blast_records.append(blast_record)

    return blast_records


def format_blast_results(blast_records, max_hits=50):
    """
    Format BLAST records into a structured dictionary

    Args:
        blast_records: List of BLAST record objects
        max_hits: Maximum number of hits to return

    Returns:
        Dictionary with formatted results
    """
    results = []

    for record in blast_records:
        for alignment in record.alignments[:max_hits]:
            for hsp in alignment.hsps:
                hit = {
                    'accession': alignment.accession,
                    'title': alignment.title,
                    'length': alignment.length,
                    'score': hsp.score,
                    'bits': hsp.bits,
                    'e_value': hsp.expect,
                    'identities': hsp.identities,
                    'positives': hsp.positives,
                    'gaps': hsp.gaps,
                    'align_length': hsp.align_length,
                    'query_start': hsp.query_start,
                    'query_end': hsp.query_end,
                    'sbjct_start': hsp.sbjct_start,
                    'sbjct_end': hsp.sbjct_end,
                    'query': hsp.query,
                    'match': hsp.match,
                    'sbjct': hsp.sbjct,
                    'frame': hsp.frame if hasattr(hsp, 'frame') else None,
                    'identity_percent': round((hsp.identities / hsp.align_length) * 100, 2),
                    'coverage_percent': round((hsp.align_length / alignment.length) * 100, 2)
                }
                results.append(hit)

    return results


def get_alignment_view(hit):
    """
    Generate formatted alignment view for a hit

    Args:
        hit: Dictionary with hit information

    Returns:
        Formatted alignment string
    """
    query = hit['query']
    match = hit['match']
    sbjct = hit['sbjct']

    # Format in blocks of 60 characters
    block_size = 60
    alignment_lines = []

    for i in range(0, len(query), block_size):
        query_block = query[i:i+block_size]
        match_block = match[i:i+block_size]
        sbjct_block = sbjct[i:i+block_size]

        query_pos = hit['query_start'] + i
        sbjct_pos = hit['sbjct_start'] + i

        alignment_lines.append(f"Query {query_pos:6d} {query_block} {query_pos + len(query_block):6d}")
        alignment_lines.append(f"            {match_block}")
        alignment_lines.append(f"Sbjct {sbjct_pos:6d} {sbjct_block} {sbjct_pos + len(sbjct_block):6d}")
        alignment_lines.append("")

    return '\n'.join(alignment_lines)


def validate_sequence(sequence, seq_type='nucleotide'):
    """
    Validate sequence format

    Args:
        sequence: Sequence string
        seq_type: 'nucleotide' or 'protein'

    Returns:
        Tuple (is_valid, error_message)
    """
    sequence = sequence.strip().replace(' ', '').replace('\n', '').upper()

    if not sequence:
        return False, "Sequence is empty"

    if seq_type == 'nucleotide':
        valid_chars = set('ATCGUN')
        invalid = set(sequence) - valid_chars
        if invalid:
            return False, f"Invalid nucleotide characters: {', '.join(sorted(invalid))}"
    elif seq_type == 'protein':
        valid_chars = set('ACDEFGHIKLMNPQRSTVWYX*')
        invalid = set(sequence) - valid_chars
        if invalid:
            return False, f"Invalid amino acid characters: {', '.join(sorted(invalid))}"

    return True, None


def get_blast_program_info(program):
    """
    Get information about BLAST program

    Args:
        program: BLAST program name

    Returns:
        Dictionary with program information
    """
    info = {
        'blastn': {
            'name': 'blastn',
            'query_type': 'nucleotide',
            'database_type': 'nucleotide',
            'description': 'Nucleotide query vs nucleotide database'
        },
        'blastp': {
            'name': 'blastp',
            'query_type': 'protein',
            'database_type': 'protein',
            'description': 'Protein query vs protein database'
        },
        'blastx': {
            'name': 'blastx',
            'query_type': 'nucleotide',
            'database_type': 'protein',
            'description': 'Translated nucleotide query vs protein database'
        },
        'tblastn': {
            'name': 'tblastn',
            'query_type': 'protein',
            'database_type': 'nucleotide',
            'description': 'Protein query vs translated nucleotide database'
        },
        'tblastx': {
            'name': 'tblastx',
            'query_type': 'nucleotide',
            'database_type': 'nucleotide',
            'description': 'Translated nucleotide query vs translated nucleotide database'
        }
    }

    return info.get(program, {})


def filter_results_by_evalue(results, max_evalue):
    """
    Filter BLAST results by E-value threshold

    Args:
        results: List of result dictionaries
        max_evalue: Maximum E-value threshold

    Returns:
        Filtered list of results
    """
    return [r for r in results if r['e_value'] <= max_evalue]


def filter_results_by_identity(results, min_identity):
    """
    Filter BLAST results by identity percentage

    Args:
        results: List of result dictionaries
        min_identity: Minimum identity percentage

    Returns:
        Filtered list of results
    """
    return [r for r in results if r['identity_percent'] >= min_identity]


def get_result_statistics(results):
    """
    Calculate statistics for BLAST results

    Args:
        results: List of result dictionaries

    Returns:
        Dictionary with statistics
    """
    if not results:
        return {}

    identities = [r['identity_percent'] for r in results]
    evalues = [r['e_value'] for r in results]
    scores = [r['score'] for r in results]

    stats = {
        'total_hits': len(results),
        'avg_identity': round(sum(identities) / len(identities), 2),
        'max_identity': max(identities),
        'min_identity': min(identities),
        'avg_evalue': sum(evalues) / len(evalues),
        'min_evalue': min(evalues),
        'max_score': max(scores),
        'avg_score': round(sum(scores) / len(scores), 2)
    }

    return stats
