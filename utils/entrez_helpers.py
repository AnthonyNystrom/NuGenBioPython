"""
Helper functions for NCBI Entrez database operations using Bio.Entrez
"""
from Bio import Entrez
from io import StringIO
import xml.etree.ElementTree as ET


def search_entrez(database, term, email, retmax=20, retstart=0,
                  sort='relevance', date_from=None, date_to=None, field=None):
    """
    Search NCBI Entrez database with advanced options

    Args:
        database: NCBI database name
        term: Search term
        email: User email (required by NCBI)
        retmax: Maximum results to return
        retstart: Starting index for pagination
        sort: Sort order (relevance, pub_date, etc.)
        date_from: Start date filter (YYYY/MM/DD)
        date_to: End date filter (YYYY/MM/DD)
        field: Specific field to search

    Returns:
        Dictionary with search results
    """
    Entrez.email = email

    # Build search term with date range if provided
    search_term = term
    if date_from and date_to:
        search_term += f" AND {date_from}:{date_to}[PDAT]"
    elif date_from:
        search_term += f" AND {date_from}:3000[PDAT]"

    # Add field restriction if provided
    if field:
        search_term += f"[{field}]"

    # Perform search
    handle = Entrez.esearch(
        db=database,
        term=search_term,
        retmax=retmax,
        retstart=retstart,
        sort=sort
    )
    results = Entrez.read(handle)
    handle.close()

    return results


def fetch_summaries(database, id_list, email):
    """
    Fetch summaries for list of IDs

    Args:
        database: NCBI database name
        id_list: List of record IDs
        email: User email

    Returns:
        List of summary dictionaries
    """
    Entrez.email = email

    if not id_list:
        return []

    ids = ','.join([str(i) for i in id_list])
    handle = Entrez.esummary(db=database, id=ids)
    summaries = Entrez.read(handle)
    handle.close()

    return summaries


def global_query(term, email):
    """
    Perform global query across all NCBI databases

    Args:
        term: Search term
        email: User email

    Returns:
        Dictionary with hit counts per database
    """
    Entrez.email = email

    handle = Entrez.egquery(term=term)
    results = Entrez.read(handle)
    handle.close()

    # Parse results
    db_stats = []
    for result in results['eGQueryResult']:
        db_stats.append({
            'database': result['DbName'],
            'count': int(result['Count']),
            'status': result['Status']
        })

    # Sort by count descending
    db_stats.sort(key=lambda x: x['count'], reverse=True)

    return db_stats


def fetch_records(database, id_list, email, rettype='fasta', retmode='text'):
    """
    Fetch full records from NCBI

    Args:
        database: NCBI database name
        id_list: List of record IDs
        email: User email
        rettype: Return type (fasta, gb, gp, xml, etc.)
        retmode: Return mode (text, xml)

    Returns:
        String containing fetched records
    """
    Entrez.email = email

    if not id_list:
        return ""

    ids = ','.join([str(i) for i in id_list])

    handle = Entrez.efetch(
        db=database,
        id=ids,
        rettype=rettype,
        retmode=retmode
    )

    records = handle.read()
    handle.close()

    # Convert bytes to string if needed
    if isinstance(records, bytes):
        records = records.decode('utf-8')

    return records


def find_related_records(record_id, from_db, to_db, email):
    """
    Find records in other databases related to given record

    Args:
        record_id: Source record ID
        from_db: Source database
        to_db: Target database
        email: User email

    Returns:
        List of related record IDs
    """
    Entrez.email = email

    handle = Entrez.elink(
        dbfrom=from_db,
        db=to_db,
        id=record_id
    )
    results = Entrez.read(handle)
    handle.close()

    # Extract linked IDs
    linked_ids = []
    if results and len(results) > 0:
        link_set_dbs = results[0].get('LinkSetDb', [])
        for link_set in link_set_dbs:
            for link in link_set.get('Link', []):
                linked_ids.append(link['Id'])

    return linked_ids


def get_database_info(database, email):
    """
    Get information about NCBI database

    Args:
        database: NCBI database name (empty string for all databases list)
        email: User email

    Returns:
        Dictionary with database information
    """
    Entrez.email = email

    if database:
        # Get info for specific database
        handle = Entrez.einfo(db=database)
    else:
        # Get list of all databases
        handle = Entrez.einfo()

    results = Entrez.read(handle)
    handle.close()

    return results


def parse_pubmed_summary(summary):
    """
    Parse PubMed summary into clean dictionary

    Args:
        summary: PubMed summary from Entrez

    Returns:
        Dictionary with parsed fields
    """
    return {
        'id': summary.get('Id', ''),
        'title': summary.get('Title', ''),
        'authors': ', '.join(summary.get('AuthorList', [])[:3]) + ('...' if len(summary.get('AuthorList', [])) > 3 else ''),
        'journal': summary.get('Source', ''),
        'date': summary.get('PubDate', ''),
        'doi': summary.get('DOI', ''),
        'pmid': summary.get('Id', '')
    }


def parse_nucleotide_summary(summary):
    """
    Parse Nucleotide summary into clean dictionary

    Args:
        summary: Nucleotide summary from Entrez

    Returns:
        Dictionary with parsed fields
    """
    return {
        'id': summary.get('Id', ''),
        'accession': summary.get('AccessionVersion', summary.get('Caption', '')),
        'title': summary.get('Title', ''),
        'organism': summary.get('Organism', ''),
        'length': summary.get('Length', ''),
        'date': summary.get('UpdateDate', '')
    }


def parse_protein_summary(summary):
    """
    Parse Protein summary into clean dictionary

    Args:
        summary: Protein summary from Entrez

    Returns:
        Dictionary with parsed fields
    """
    return {
        'id': summary.get('Id', ''),
        'accession': summary.get('AccessionVersion', summary.get('Caption', '')),
        'title': summary.get('Title', ''),
        'organism': summary.get('Organism', ''),
        'length': summary.get('Length', ''),
        'date': summary.get('UpdateDate', '')
    }


def parse_gene_summary(summary):
    """
    Parse Gene summary into clean dictionary

    Args:
        summary: Gene summary from Entrez

    Returns:
        Dictionary with parsed fields
    """
    return {
        'id': summary.get('Id', ''),
        'name': summary.get('Name', ''),
        'description': summary.get('Description', ''),
        'organism': summary.get('Organism', summary.get('OrganismName', '')),
        'chromosome': summary.get('Chromosome', ''),
        'map_location': summary.get('MapLocation', '')
    }


def format_search_results(database, summaries):
    """
    Format search results based on database type

    Args:
        database: Database name
        summaries: List of summary records

    Returns:
        List of formatted result dictionaries
    """
    formatted = []

    for summary in summaries:
        if database == 'pubmed':
            formatted.append(parse_pubmed_summary(summary))
        elif database == 'nucleotide':
            formatted.append(parse_nucleotide_summary(summary))
        elif database == 'protein':
            formatted.append(parse_protein_summary(summary))
        elif database == 'gene':
            formatted.append(parse_gene_summary(summary))
        else:
            # Generic format
            formatted.append({
                'id': summary.get('Id', ''),
                'title': str(summary.get('Title', summary.get('Name', 'N/A'))),
                'info': str(summary)[:200] + '...'
            })

    return formatted
