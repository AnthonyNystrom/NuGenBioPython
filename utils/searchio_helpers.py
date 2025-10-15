"""
Helper functions for Bio.SearchIO operations
Handles parsing, reading, indexing, converting, filtering, and writing search results
"""
from Bio import SearchIO
from io import StringIO
import os
import tempfile


def searchio_parse(file_path, format_name):
    """
    Parse search results file containing one or more queries

    Args:
        file_path: Path to search results file
        format_name: Format identifier (blast-xml, blast-tab, hmmer3-tab, etc.)

    Returns:
        List of query results with hits
    """
    results = []
    try:
        for query_result in SearchIO.parse(file_path, format_name):
            query_data = {
                'id': str(getattr(query_result, 'id', 'Unknown')),
                'description': str(getattr(query_result, 'description', '')),
                'seq_len': getattr(query_result, 'seq_len', 'N/A'),
                'hits': []
            }

            # Get top 20 hits
            for hit in query_result[:20]:
                hit_data = {
                    'id': str(getattr(hit, 'id', 'Unknown')),
                    'description': str(getattr(hit, 'description', '')),
                    'evalue': 'N/A',
                    'bitscore': 'N/A',
                    'length': getattr(hit, 'seq_len', 'N/A')
                }

                # Extract evalue and bitscore from first HSP
                if hasattr(hit, 'hsps') and len(hit.hsps) > 0:
                    hsp = hit.hsps[0]
                    if hasattr(hsp, 'evalue'):
                        hit_data['evalue'] = f"{hsp.evalue:.2e}" if isinstance(hsp.evalue, (int, float)) else str(hsp.evalue)
                    if hasattr(hsp, 'bitscore'):
                        hit_data['bitscore'] = f"{hsp.bitscore:.1f}" if isinstance(hsp.bitscore, (int, float)) else str(hsp.bitscore)

                query_data['hits'].append(hit_data)

            results.append(query_data)

        return results
    except Exception as e:
        raise Exception(f"Parse error: {str(e)}")


def searchio_read(file_path, format_name):
    """
    Read a search results file containing exactly one query

    Args:
        file_path: Path to search results file
        format_name: Format identifier

    Returns:
        Single query result with detailed information
    """
    try:
        query_result = SearchIO.read(file_path, format_name)

        result = {
            'id': str(getattr(query_result, 'id', 'Unknown')),
            'description': str(getattr(query_result, 'description', '')),
            'seq_len': getattr(query_result, 'seq_len', 'N/A'),
            'num_hits': len(query_result),
            'hits': []
        }

        # Get all hits with detailed HSP information
        for hit in query_result[:50]:
            hit_data = {
                'id': str(getattr(hit, 'id', 'Unknown')),
                'description': str(getattr(hit, 'description', '')),
                'length': getattr(hit, 'seq_len', 'N/A'),
                'num_hsps': len(hit.hsps) if hasattr(hit, 'hsps') else 0,
                'hsps': []
            }

            # Extract HSP details
            if hasattr(hit, 'hsps'):
                for hsp in hit.hsps[:5]:  # Top 5 HSPs per hit
                    hsp_data = {
                        'evalue': f"{hsp.evalue:.2e}" if hasattr(hsp, 'evalue') and isinstance(hsp.evalue, (int, float)) else 'N/A',
                        'bitscore': f"{hsp.bitscore:.1f}" if hasattr(hsp, 'bitscore') and isinstance(hsp.bitscore, (int, float)) else 'N/A',
                        'query_start': getattr(hsp, 'query_start', 'N/A'),
                        'query_end': getattr(hsp, 'query_end', 'N/A'),
                        'hit_start': getattr(hsp, 'hit_start', 'N/A'),
                        'hit_end': getattr(hsp, 'hit_end', 'N/A')
                    }
                    hit_data['hsps'].append(hsp_data)

            result['hits'].append(hit_data)

        return result
    except Exception as e:
        raise Exception(f"Read error: {str(e)}")


def searchio_index(file_path, format_name, limit=100):
    """
    Create indexed dictionary-like access to search results
    Useful for large files with many queries

    Args:
        file_path: Path to search results file
        format_name: Format identifier
        limit: Maximum number of indexed queries to return

    Returns:
        Dictionary mapping query IDs to basic information
    """
    try:
        indexed = SearchIO.index(file_path, format_name)

        results = {
            'total_queries': len(indexed),
            'queries': []
        }

        # Get first N query IDs and basic info
        for i, query_id in enumerate(indexed.keys()):
            if i >= limit:
                break

            query = indexed[query_id]
            query_info = {
                'id': str(query_id),
                'description': str(getattr(query, 'description', '')),
                'seq_len': getattr(query, 'seq_len', 'N/A'),
                'num_hits': len(query)
            }
            results['queries'].append(query_info)

        indexed.close()
        return results
    except Exception as e:
        raise Exception(f"Index error: {str(e)}")


def searchio_convert(input_path, input_format, output_format):
    """
    Convert search results from one format to another

    Args:
        input_path: Path to input file
        input_format: Input format identifier
        output_format: Output format identifier

    Returns:
        Path to converted file
    """
    try:
        # Create temporary output file
        output_fd, output_path = tempfile.mkstemp(suffix=f'.{output_format}', text=True)
        os.close(output_fd)

        # Convert using SearchIO - Note: Not all format combinations are supported
        try:
            count = SearchIO.convert(input_path, input_format, output_path, output_format)
        except (AttributeError, ValueError, NotImplementedError) as conv_error:
            # Some conversions aren't supported - return helpful message
            raise Exception(f"Conversion from {input_format} to {output_format} not fully supported by BioPython. Try using blast-xml to blast-xml or hmmer formats.")

        # Read converted content
        with open(output_path, 'r') as f:
            content = f.read()

        return {
            'count': count,
            'output_path': output_path,
            'content_preview': content[:2000] if len(content) > 2000 else content,
            'content_size': len(content)
        }
    except Exception as e:
        if "not fully supported" in str(e):
            raise
        raise Exception(f"Convert error: {str(e)}")


def searchio_filter(file_path, format_name, evalue_threshold=None, bitscore_threshold=None, min_identity=None):
    """
    Filter search results by E-value, bit score, or identity

    Args:
        file_path: Path to search results file
        format_name: Format identifier
        evalue_threshold: Maximum E-value (hits with E-value > threshold are excluded)
        bitscore_threshold: Minimum bit score
        min_identity: Minimum identity percentage

    Returns:
        Filtered results
    """
    results = []
    try:
        for query_result in SearchIO.parse(file_path, format_name):
            query_data = {
                'id': str(getattr(query_result, 'id', 'Unknown')),
                'description': str(getattr(query_result, 'description', '')),
                'seq_len': getattr(query_result, 'seq_len', 'N/A'),
                'hits': [],
                'filtered_count': 0,
                'original_count': len(query_result)
            }

            for hit in query_result:
                # Check if hit passes filters
                passes_filter = True

                if hasattr(hit, 'hsps') and len(hit.hsps) > 0:
                    hsp = hit.hsps[0]

                    # E-value filter
                    if evalue_threshold is not None and hasattr(hsp, 'evalue'):
                        try:
                            if float(hsp.evalue) > float(evalue_threshold):
                                passes_filter = False
                        except:
                            pass

                    # Bit score filter
                    if bitscore_threshold is not None and hasattr(hsp, 'bitscore'):
                        try:
                            if float(hsp.bitscore) < float(bitscore_threshold):
                                passes_filter = False
                        except:
                            pass

                    # Identity filter
                    if min_identity is not None and hasattr(hsp, 'ident_pct'):
                        try:
                            if float(hsp.ident_pct) < float(min_identity):
                                passes_filter = False
                        except:
                            pass

                    if passes_filter:
                        hit_data = {
                            'id': str(getattr(hit, 'id', 'Unknown')),
                            'description': str(getattr(hit, 'description', '')),
                            'evalue': f"{hsp.evalue:.2e}" if hasattr(hsp, 'evalue') and isinstance(hsp.evalue, (int, float)) else 'N/A',
                            'bitscore': f"{hsp.bitscore:.1f}" if hasattr(hsp, 'bitscore') and isinstance(hsp.bitscore, (int, float)) else 'N/A',
                            'identity': f"{hsp.ident_pct:.1f}%" if hasattr(hsp, 'ident_pct') else 'N/A'
                        }
                        query_data['hits'].append(hit_data)

            query_data['filtered_count'] = len(query_data['hits'])
            results.append(query_data)

        return results
    except Exception as e:
        raise Exception(f"Filter error: {str(e)}")


def searchio_write(file_path, format_name, output_format, max_queries=10):
    """
    Write search results to a different format

    Args:
        file_path: Path to input search results file
        format_name: Input format identifier
        output_format: Output format identifier
        max_queries: Maximum number of queries to write

    Returns:
        Written content and statistics
    """
    try:
        # Read queries from input file
        queries = []
        for i, qresult in enumerate(SearchIO.parse(file_path, format_name)):
            if i >= max_queries:
                break
            queries.append(qresult)

        # Write to temporary file
        output_fd, output_path = tempfile.mkstemp(suffix=f'.{output_format}', text=True)
        os.close(output_fd)

        count = SearchIO.write(queries, output_path, output_format)

        # Read written content
        with open(output_path, 'r') as f:
            content = f.read()

        return {
            'count': count,
            'output_path': output_path,
            'content': content[:3000] if len(content) > 3000 else content,
            'content_size': len(content)
        }
    except Exception as e:
        raise Exception(f"Write error: {str(e)}")


def get_supported_formats():
    """
    Get list of all supported SearchIO formats

    Returns:
        Dictionary of format categories with format details
    """
    formats = {
        'BLAST': [
            {'id': 'blast-xml', 'name': 'BLAST XML', 'description': 'NCBI BLAST XML output', 'can_write': True},
            {'id': 'blast-tab', 'name': 'BLAST Tabular', 'description': 'BLAST tabular format (-outfmt 6/7)', 'can_write': True},
            {'id': 'blast-text', 'name': 'BLAST Text', 'description': 'BLAST pairwise text output', 'can_write': False}
        ],
        'HMMER': [
            {'id': 'hmmer3-text', 'name': 'HMMER3 Text', 'description': 'HMMER3 standard output', 'can_write': False},
            {'id': 'hmmer3-tab', 'name': 'HMMER3 Table', 'description': 'HMMER3 table output (--tblout)', 'can_write': True},
            {'id': 'hmmer3-domtab', 'name': 'HMMER3 Domain Table', 'description': 'HMMER3 domain table (--domtblout)', 'can_write': True}
        ],
        'Other': [
            {'id': 'blat-psl', 'name': 'BLAT PSL', 'description': 'BLAT PSL output format', 'can_write': True},
            {'id': 'fasta-m10', 'name': 'FASTA -m 10', 'description': 'FASTA output format (-m 10)', 'can_write': False},
            {'id': 'exonerate-text', 'name': 'Exonerate Text', 'description': 'Exonerate text output', 'can_write': False},
            {'id': 'exonerate-cigar', 'name': 'Exonerate CIGAR', 'description': 'Exonerate CIGAR output', 'can_write': False}
        ]
    }

    return formats
