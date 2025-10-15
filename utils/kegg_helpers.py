"""
Helper functions for KEGG database operations using Bio.KEGG.REST
"""
from Bio.KEGG import REST


def kegg_find_search(database, query, organism=None):
    """
    Search KEGG database using kegg_find

    Args:
        database: KEGG database (pathway, genes, compound, etc.)
        query: Search query term
        organism: Organism code (optional, for gene searches)

    Returns:
        List of dictionaries with search results
    """
    results = []

    try:
        # Construct search query
        if database == 'genes' and organism:
            search_result = REST.kegg_find(organism, query)
        else:
            search_result = REST.kegg_find(database, query)

        # Parse results
        search_data = search_result.read().strip()

        if isinstance(search_data, bytes):
            search_data = search_data.decode('utf-8')

        if not search_data:
            return []

        for line in search_data.split('\n'):
            if line and '\t' in line:
                parts = line.split('\t', 1)
                if len(parts) == 2:
                    entry_id, definition = parts
                    results.append({
                        'id': entry_id.strip(),
                        'definition': definition.strip(),
                        'database': database
                    })

        return results

    except Exception as e:
        raise Exception(f"KEGG find error: {str(e)}")


def kegg_list_entries(database, organism=None):
    """
    List all entries in a KEGG database using kegg_list

    Args:
        database: KEGG database name or 'organism' for organism list
        organism: Organism code (optional)

    Returns:
        List of dictionaries with entry information
    """
    results = []

    try:
        # Get list of entries
        if organism:
            list_result = REST.kegg_list(f"{organism}:{database}")
        else:
            list_result = REST.kegg_list(database)

        list_data = list_result.read().strip()

        if isinstance(list_data, bytes):
            list_data = list_data.decode('utf-8')

        if not list_data:
            return []

        for line in list_data.split('\n'):
            if line and '\t' in line:
                parts = line.split('\t', 1)
                if len(parts) == 2:
                    entry_id, definition = parts
                    results.append({
                        'id': entry_id.strip(),
                        'definition': definition.strip()
                    })
            elif line:
                # Some entries might not have tab separation
                results.append({
                    'id': line.strip(),
                    'definition': ''
                })

        return results

    except Exception as e:
        raise Exception(f"KEGG list error: {str(e)}")


def kegg_link_entries(target_db, source_db, source_id=None):
    """
    Find related entries using kegg_link

    Args:
        target_db: Target database to link to
        source_db: Source database
        source_id: Specific entry ID (optional)

    Returns:
        List of dictionaries with linked entries
    """
    results = []

    try:
        # Get links
        if source_id:
            link_result = REST.kegg_link(target_db, source_id)
        else:
            link_result = REST.kegg_link(target_db, source_db)

        link_data = link_result.read().strip()

        if isinstance(link_data, bytes):
            link_data = link_data.decode('utf-8')

        if not link_data:
            return []

        for line in link_data.split('\n'):
            if line and '\t' in line:
                parts = line.split('\t')
                if len(parts) == 2:
                    source, target = parts
                    results.append({
                        'source': source.strip(),
                        'target': target.strip(),
                        'source_db': source_db,
                        'target_db': target_db
                    })

        return results

    except Exception as e:
        raise Exception(f"KEGG link error: {str(e)}")


def kegg_convert_ids(target_db, source_db, ids=None):
    """
    Convert identifiers between databases using kegg_conv

    Args:
        target_db: Target database ID system
        source_db: Source database ID system
        ids: List of IDs to convert (optional)

    Returns:
        List of dictionaries with ID conversions
    """
    results = []

    try:
        # Perform conversion
        if ids:
            # Convert specific IDs
            if isinstance(ids, list):
                ids = ' '.join(ids)
            conv_result = REST.kegg_conv(target_db, ids)
        else:
            # Get all conversions between databases
            conv_result = REST.kegg_conv(target_db, source_db)

        conv_data = conv_result.read().strip()

        if isinstance(conv_data, bytes):
            conv_data = conv_data.decode('utf-8')

        if not conv_data:
            return []

        for line in conv_data.split('\n'):
            if line and '\t' in line:
                parts = line.split('\t')
                if len(parts) == 2:
                    source_id, target_id = parts
                    results.append({
                        'source_id': source_id.strip(),
                        'target_id': target_id.strip(),
                        'source_db': source_db,
                        'target_db': target_db
                    })

        return results

    except Exception as e:
        raise Exception(f"KEGG conv error: {str(e)}")


def kegg_get_info(database):
    """
    Get database information using kegg_info

    Args:
        database: Database name (or None for list of all databases)

    Returns:
        Dictionary with database information
    """
    try:
        # Get info
        if database:
            info_result = REST.kegg_info(database)
        else:
            info_result = REST.kegg_info('kegg')

        info_data = info_result.read().strip()

        if isinstance(info_data, bytes):
            info_data = info_data.decode('utf-8')

        # Parse info data
        info = {
            'raw_data': info_data,
            'fields': {}
        }

        for line in info_data.split('\n'):
            if line and line.strip():
                if ':' in line or '  ' in line:
                    parts = line.split(None, 1)
                    if len(parts) == 2:
                        key, value = parts
                        info['fields'][key.strip()] = value.strip()

        return info

    except Exception as e:
        raise Exception(f"KEGG info error: {str(e)}")


def kegg_get_entry(entry_id):
    """
    Get detailed entry information using kegg_get

    Args:
        entry_id: KEGG entry identifier

    Returns:
        String with entry data
    """
    try:
        entry_result = REST.kegg_get(entry_id)
        entry_data = entry_result.read()

        if isinstance(entry_data, bytes):
            entry_data = entry_data.decode('utf-8')

        return entry_data

    except Exception as e:
        raise Exception(f"KEGG get error: {str(e)}")


def parse_kegg_entry(entry_data):
    """
    Parse KEGG entry format into structured data

    Args:
        entry_data: Raw KEGG entry string

    Returns:
        Dictionary with parsed fields
    """
    parsed = {}
    current_field = None
    current_value = []

    for line in entry_data.split('\n'):
        if not line.strip() or line.startswith('///'):
            continue

        # Check if line starts a new field
        if line and not line.startswith(' '):
            # Save previous field
            if current_field:
                parsed[current_field] = '\n'.join(current_value)

            # Start new field
            parts = line.split(None, 1)
            if len(parts) == 2:
                current_field = parts[0]
                current_value = [parts[1]]
            else:
                current_field = parts[0]
                current_value = []
        else:
            # Continuation of current field
            if current_field:
                current_value.append(line.strip())

    # Save last field
    if current_field:
        parsed[current_field] = '\n'.join(current_value)

    return parsed
