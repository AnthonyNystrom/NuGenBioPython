"""
Helper functions for NCBI Entrez database operations using Bio.Entrez
"""
import logging
import os

from Bio import Entrez
from io import StringIO
import xml.etree.ElementTree as ET

from utils.request_helpers import remote_timeout

log = logging.getLogger(__name__)

# Entrez ceiling — NCBI usually responds in <5s but large queries can run longer.
_ENTREZ_TIMEOUT = 30

# Addresses that identify nobody. Callers used to default to one of these when
# the client sent no email, which is worse than sending nothing: it labels every
# request from every deployment as the same fake contact, and that is what gets
# a host IP throttled or blocked by NCBI.
_PLACEHOLDER_EMAILS = {
    'user@example.com',
    'your.email@example.com',
    'your@email.com',
    'email@example.com',
}


def _looks_like_email(value):
    value = (value or '').strip()
    if '@' not in value or value.lower() in _PLACEHOLDER_EMAILS:
        return False
    local, _, domain = value.partition('@')
    return bool(local) and '.' in domain and not domain.endswith('.')


def configure_entrez(email=None):
    """Set the NCBI identity used by the next Entrez call, and return it.

    NCBI's usage policy asks every request to carry a real contact address,
    and grants a higher rate ceiling (10 req/s rather than 3) when an API key
    is supplied. Resolution order:

      1. a valid address supplied by the caller,
      2. the ENTREZ_EMAIL environment variable,
      3. nothing — the call still proceeds, but is logged as unidentified so
         an operator can see why NCBI may be throttling them.

    ENTREZ_API_KEY is applied whenever it is set. Entrez.tool is always set so
    NCBI can attribute traffic to this application.
    """
    Entrez.tool = 'NuGenBioPython'

    resolved = email if _looks_like_email(email) else None
    if resolved is None:
        env_email = os.environ.get('ENTREZ_EMAIL', '')
        resolved = env_email if _looks_like_email(env_email) else None

    if resolved:
        Entrez.email = resolved
    else:
        log.warning(
            'Entrez request has no valid contact email. NCBI asks every '
            'request to identify a real address; set ENTREZ_EMAIL (and '
            'optionally ENTREZ_API_KEY for a 10 req/s ceiling) to avoid '
            'being throttled or blocked.'
        )

    api_key = os.environ.get('ENTREZ_API_KEY', '').strip()
    if api_key:
        Entrez.api_key = api_key

    return resolved


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
    configure_entrez(email)

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
    with remote_timeout(_ENTREZ_TIMEOUT, service='ncbi'):
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
    configure_entrez(email)

    if not id_list:
        return []

    ids = ','.join([str(i) for i in id_list])
    with remote_timeout(_ENTREZ_TIMEOUT, service='ncbi'):
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
    configure_entrez(email)

    with remote_timeout(_ENTREZ_TIMEOUT, service='ncbi'):
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
    configure_entrez(email)

    if not id_list:
        return ""

    ids = ','.join([str(i) for i in id_list])

    with remote_timeout(_ENTREZ_TIMEOUT, service='ncbi'):
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
    configure_entrez(email)

    with remote_timeout(_ENTREZ_TIMEOUT, service='ncbi'):
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
    configure_entrez(email)

    with remote_timeout(_ENTREZ_TIMEOUT, service='ncbi'):
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
