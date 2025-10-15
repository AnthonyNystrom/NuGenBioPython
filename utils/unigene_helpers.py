"""
Helper functions for Bio.UniGene operations
Handles parsing and reading UniGene cluster data
"""
from Bio import UniGene


def unigene_parse(file_path, max_records=10):
    """
    Parse multiple UniGene records from a file

    Args:
        file_path: Path to UniGene .data file
        max_records: Maximum number of records to return

    Returns:
        List of parsed UniGene records with comprehensive data
    """
    records = []
    try:
        with open(file_path, 'r') as handle:
            for record in UniGene.parse(handle):
                sequences = getattr(record, 'sequence', [])
                protsims = getattr(record, 'protsim', [])
                sts_entries = getattr(record, 'sts', [])
                express_data = getattr(record, 'express', [])

                # Extract tissue names from expression data
                tissues = []
                for expr in express_data:
                    tissue = getattr(expr, 'tissue', None)
                    if tissue and tissue not in tissues:
                        tissues.append(tissue)

                record_data = {
                    'cluster_id': str(getattr(record, 'ID', 'Unknown')),
                    'species': str(getattr(record, 'species', 'N/A')),
                    'gene_symbol': str(getattr(record, 'symbol', 'N/A')) if getattr(record, 'symbol', None) else 'N/A',
                    'title': str(getattr(record, 'title', 'Unknown')),
                    'chromosome': str(getattr(record, 'chromosome', 'N/A')) if getattr(record, 'chromosome', None) else 'N/A',
                    'sequence_count': len(sequences),
                    'tissues': tissues[:10],  # Limit to 10 tissues
                    'tissue_count': len(tissues)
                }

                records.append(record_data)

                if len(records) >= max_records:
                    break

        return records
    except Exception as e:
        raise Exception(f"Parse error: {str(e)}")


def unigene_read(file_path):
    """
    Read a single UniGene record from a file with full details

    Args:
        file_path: Path to UniGene .data file containing one record

    Returns:
        Detailed UniGene record information
    """
    try:
        with open(file_path, 'r') as handle:
            record = UniGene.read(handle)

            sequences = getattr(record, 'sequence', [])
            protsims = getattr(record, 'protsim', [])
            sts_entries = getattr(record, 'sts', [])
            express_data = getattr(record, 'express', [])

            # Extract tissue expression details
            tissues = []
            tissue_details = []
            for expr in express_data:
                tissue = getattr(expr, 'tissue', None)
                freq = getattr(expr, 'freq', None)
                if tissue:
                    if tissue not in tissues:
                        tissues.append(tissue)
                    tissue_details.append({
                        'tissue': str(tissue),
                        'frequency': int(freq) if freq else 0
                    })

            record_data = {
                'cluster_id': str(getattr(record, 'ID', 'Unknown')),
                'species': str(getattr(record, 'species', 'N/A')),
                'gene_symbol': str(getattr(record, 'symbol', 'N/A')) if getattr(record, 'symbol', None) else 'N/A',
                'title': str(getattr(record, 'title', 'Unknown')),
                'chromosome': str(getattr(record, 'chromosome', 'N/A')) if getattr(record, 'chromosome', None) else 'N/A',
                'cytoband': str(getattr(record, 'cytoband', 'N/A')) if getattr(record, 'cytoband', None) else 'N/A',
                'gene_id': str(getattr(record, 'gene_id', 'N/A')) if getattr(record, 'gene_id', None) else 'N/A',
                'locuslink': str(getattr(record, 'locuslink', 'N/A')) if getattr(record, 'locuslink', None) else 'N/A',
                'homol': str(getattr(record, 'homol', 'N/A')) if getattr(record, 'homol', None) else 'N/A',
                'restr_expr': str(getattr(record, 'restr_expr', 'N/A')) if getattr(record, 'restr_expr', None) else 'N/A',
                'gnm_terminus': str(getattr(record, 'gnm_terminus', 'N/A')) if getattr(record, 'gnm_terminus', None) else 'N/A',
                'txmap': str(getattr(record, 'txmap', 'N/A')) if getattr(record, 'txmap', None) else 'N/A',
                'sequence_count': len(sequences),
                'tissue_count': len(tissues),
                'tissues': tissues,
                'tissue_details': tissue_details,
                'sequences': [],
                'protein_similarities': [],
                'sts': []
            }

            # Add sequence details (all sequences)
            for seq in sequences:
                seq_data = {
                    'acc': str(getattr(seq, 'acc', 'N/A')),
                    'seqtype': str(getattr(seq, 'seqtype', 'N/A')),
                    'clone': str(getattr(seq, 'clone', 'N/A')) if getattr(seq, 'clone', None) else 'N/A',
                    'end': str(getattr(seq, 'end', 'N/A')) if getattr(seq, 'end', None) else 'N/A',
                    'lid': str(getattr(seq, 'lid', 'N/A')) if getattr(seq, 'lid', None) else 'N/A'
                }
                record_data['sequences'].append(seq_data)

            # Add protein similarity details (all)
            for prot in protsims:
                prot_data = {
                    'organism': str(getattr(prot, 'org', 'N/A')),
                    'protein_id': str(getattr(prot, 'protid', 'N/A')),
                    'percent': str(getattr(prot, 'pct', 'N/A')) if getattr(prot, 'pct', None) else 'N/A',
                    'alignment_length': str(getattr(prot, 'aln', 'N/A')) if getattr(prot, 'aln', None) else 'N/A'
                }
                record_data['protein_similarities'].append(prot_data)

            # Add STS details (all)
            for sts in sts_entries:
                sts_data = {
                    'acc': str(getattr(sts, 'acc', 'N/A')),
                    'unists': str(getattr(sts, 'unists', 'N/A')) if getattr(sts, 'unists', None) else 'N/A'
                }
                record_data['sts'].append(sts_data)

            return record_data
    except Exception as e:
        raise Exception(f"Read error: {str(e)}")
