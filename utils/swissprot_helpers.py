"""
Helper functions for Bio.SwissProt operations
Handles parsing and reading SwissProt/UniProt protein database records
"""
from Bio import SwissProt


def swissprot_parse(file_path, max_records=10):
    """
    Parse multiple SwissProt records from a file

    Args:
        file_path: Path to SwissProt .dat file
        max_records: Maximum number of records to return

    Returns:
        List of parsed SwissProt records
    """
    records = []
    try:
        with open(file_path, 'r') as handle:
            for record in SwissProt.parse(handle):
                record_data = {
                    'accessions': list(record.accessions),
                    'entry_name': str(record.entry_name),
                    'description': str(record.description),
                    'gene_name': format_gene_name(record.gene_name),
                    'organism': str(record.organism),
                    'organism_classification': list(record.organism_classification) if hasattr(record, 'organism_classification') else [],
                    'sequence_length': record.sequence_length,
                    'keywords': list(record.keywords),
                    'features': [],
                    'cross_references': [],
                    'comments': []
                }

                # Add features (limit to 10)
                for feature in record.features[:10]:
                    try:
                        # SwissProt features are tuples: (type, from, to, description)
                        if isinstance(feature, tuple) and len(feature) >= 3:
                            feature_data = {
                                'type': str(feature[0]),
                                'location': f"{feature[1]}-{feature[2]}",
                                'description': str(feature[3]) if len(feature) > 3 else ''
                            }
                        else:
                            # Fallback for other formats
                            feature_data = {
                                'type': str(getattr(feature, 'type', 'Unknown')),
                                'location': str(getattr(feature, 'location', '')),
                                'description': str(getattr(feature, 'description', ''))
                            }
                        record_data['features'].append(feature_data)
                    except Exception as fe:
                        # Skip malformed features
                        continue

                # Add cross-references (limit to 10)
                for ref in record.cross_references[:10]:
                    record_data['cross_references'].append({
                        'database': str(ref[0]) if len(ref) > 0 else '',
                        'id': str(ref[1]) if len(ref) > 1 else '',
                        'info': str(ref[2]) if len(ref) > 2 else ''
                    })

                # Add comments (limit to 5)
                for comment in record.comments[:5]:
                    record_data['comments'].append(str(comment))

                records.append(record_data)

                if len(records) >= max_records:
                    break

        return records
    except Exception as e:
        raise Exception(f"Parse error: {str(e)}")


def swissprot_read(file_path):
    """
    Read a single SwissProt record from a file with full details

    Args:
        file_path: Path to SwissProt .dat file containing one record

    Returns:
        Detailed SwissProt record information
    """
    try:
        with open(file_path, 'r') as handle:
            record = SwissProt.read(handle)

            record_data = {
                'accessions': list(record.accessions),
                'entry_name': str(record.entry_name),
                'data_class': str(record.data_class) if hasattr(record, 'data_class') else 'N/A',
                'molecule_type': str(record.molecule_type) if hasattr(record, 'molecule_type') else 'N/A',
                'description': str(record.description),
                'gene_name': format_gene_name(record.gene_name),
                'organism': str(record.organism),
                'organism_classification': list(record.organism_classification) if hasattr(record, 'organism_classification') else [],
                'taxonomy_id': list(record.taxonomy_id) if hasattr(record, 'taxonomy_id') else [],
                'host_organism': list(record.host_organism) if hasattr(record, 'host_organism') else [],
                'sequence': str(record.sequence),
                'sequence_length': record.sequence_length,
                'molecular_weight': getattr(record, 'seqinfo', (0, 0))[0] if hasattr(record, 'seqinfo') else 0,
                'keywords': list(record.keywords),
                'protein_existence': str(record.protein_existence) if hasattr(record, 'protein_existence') else 'N/A',
                'seqinfo': str(record.seqinfo) if hasattr(record, 'seqinfo') else 'N/A',
                'features': [],
                'cross_references': [],
                'comments': [],
                'references': []
            }

            # Add all features
            for feature in record.features:
                try:
                    # SwissProt features are tuples: (type, from, to, description)
                    if isinstance(feature, tuple) and len(feature) >= 3:
                        feature_data = {
                            'type': str(feature[0]),
                            'location': f"{feature[1]}-{feature[2]}",
                            'description': str(feature[3]) if len(feature) > 3 else ''
                        }
                    else:
                        # Fallback for other formats
                        feature_data = {
                            'type': str(getattr(feature, 'type', 'Unknown')),
                            'location': str(getattr(feature, 'location', '')),
                            'description': str(getattr(feature, 'description', ''))
                        }
                    record_data['features'].append(feature_data)
                except Exception as fe:
                    # Skip malformed features
                    continue

            # Add all cross-references
            for ref in record.cross_references:
                record_data['cross_references'].append({
                    'database': str(ref[0]) if len(ref) > 0 else '',
                    'id': str(ref[1]) if len(ref) > 1 else '',
                    'info': str(ref[2]) if len(ref) > 2 else ''
                })

            # Add all comments
            for comment in record.comments:
                record_data['comments'].append(str(comment))

            # Add references
            for reference in record.references[:10]:  # Limit to 10 references
                ref_data = {
                    'number': str(reference.number) if hasattr(reference, 'number') else '',
                    'authors': str(reference.authors) if hasattr(reference, 'authors') else '',
                    'title': str(reference.title) if hasattr(reference, 'title') else '',
                    'location': str(reference.location) if hasattr(reference, 'location') else ''
                }
                record_data['references'].append(ref_data)

            return record_data
    except Exception as e:
        raise Exception(f"Read error: {str(e)}")


def format_gene_name(gene_name):
    """
    Format gene name from SwissProt record

    Args:
        gene_name: Gene name data from SwissProt record

    Returns:
        Formatted gene name string
    """
    if not gene_name or gene_name == 'N/A':
        return 'N/A'

    # Handle if it's already a string
    if isinstance(gene_name, str):
        return gene_name

    # Handle array of gene objects
    if isinstance(gene_name, list):
        return ', '.join([
            f"{gene.get('Name', '')} ({', '.join(gene.get('Synonyms', []))})" if gene.get('Synonyms')
            else gene.get('Name', '')
            for gene in gene_name
            if isinstance(gene, dict)
        ]) or 'N/A'

    # Handle single gene object
    if isinstance(gene_name, dict):
        result = gene_name.get('Name', '')
        synonyms = gene_name.get('Synonyms', [])
        if synonyms:
            result += f" ({', '.join(synonyms)})"
        return result or 'N/A'

    return 'N/A'
