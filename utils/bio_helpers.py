"""
Helper functions for BioPython operations
"""

def extract_pubmed_record(root):
    """Extract structured data from PubMed XML"""
    try:
        record = {}

        # Find the article element
        article = root.find('.//PubmedArticle')
        if article is not None:
            # Title
            title_elem = article.find('.//ArticleTitle')
            record['title'] = title_elem.text if title_elem is not None else 'No title'

            # Authors
            authors = []
            author_list = article.find('.//AuthorList')
            if author_list is not None:
                for author in author_list.findall('Author'):
                    last_name = author.find('LastName')
                    first_name = author.find('ForeName')
                    if last_name is not None:
                        name = last_name.text
                        if first_name is not None:
                            name = f"{first_name.text} {name}"
                        authors.append(name)
            record['authors'] = authors

            # Journal
            journal_elem = article.find('.//Journal/Title')
            record['journal'] = journal_elem.text if journal_elem is not None else 'Unknown journal'

            # Publication date
            pub_date = article.find('.//PubDate')
            if pub_date is not None:
                year = pub_date.find('Year')
                month = pub_date.find('Month')
                day = pub_date.find('Day')
                date_parts = []
                if year is not None:
                    date_parts.append(year.text)
                if month is not None:
                    date_parts.append(month.text)
                if day is not None:
                    date_parts.append(day.text)
                record['date'] = ' '.join(date_parts) if date_parts else 'Unknown date'

            # Abstract
            abstract_elem = article.find('.//AbstractText')
            record['abstract'] = abstract_elem.text if abstract_elem is not None else 'No abstract available'

            # PMID
            pmid_elem = article.find('.//PMID')
            record['pmid'] = pmid_elem.text if pmid_elem is not None else 'Unknown'

        return record
    except:
        return {'error': 'Failed to parse PubMed record'}


def generate_sample_kegg_data(database, organism, query):
    """Generate sample KEGG data for demonstration when API is not available"""
    sample_data = {
        'genes': {
            'MAPK': [
                {
                    'id': f'{organism}:5594',
                    'definition': 'mitogen-activated protein kinase 1',
                    'pathway': 'MAPK signaling pathway - Homo sapiens (human)',
                    'organism': 'Homo sapiens',
                    'class': 'Protein kinase',
                    'kegg_id': 'K04371'
                },
                {
                    'id': f'{organism}:5595',
                    'definition': 'mitogen-activated protein kinase 3',
                    'pathway': 'MAPK signaling pathway - Homo sapiens (human)',
                    'organism': 'Homo sapiens',
                    'class': 'Protein kinase',
                    'kegg_id': 'K04371'
                },
                {
                    'id': f'{organism}:1432',
                    'definition': 'mitogen-activated protein kinase kinase 1',
                    'pathway': 'MAPK signaling pathway - Homo sapiens (human)',
                    'organism': 'Homo sapiens',
                    'class': 'Protein kinase',
                    'kegg_id': 'K04368'
                }
            ],
            'insulin': [
                {
                    'id': f'{organism}:3630',
                    'definition': 'insulin',
                    'pathway': 'Insulin signaling pathway - Homo sapiens (human)',
                    'organism': 'Homo sapiens',
                    'class': 'Hormone',
                    'kegg_id': 'K04526'
                },
                {
                    'id': f'{organism}:3643',
                    'definition': 'insulin receptor',
                    'pathway': 'Insulin signaling pathway - Homo sapiens (human)',
                    'organism': 'Homo sapiens',
                    'class': 'Receptor tyrosine kinase',
                    'kegg_id': 'K04527'
                }
            ],
            'glycolysis': [
                {
                    'id': f'{organism}:2597',
                    'definition': 'glucose-6-phosphate isomerase',
                    'pathway': 'Glycolysis / Gluconeogenesis - Homo sapiens (human)',
                    'organism': 'Homo sapiens',
                    'class': 'Isomerase',
                    'kegg_id': 'K01810'
                },
                {
                    'id': f'{organism}:5211',
                    'definition': 'phosphofructokinase, liver type',
                    'pathway': 'Glycolysis / Gluconeogenesis - Homo sapiens (human)',
                    'organism': 'Homo sapiens',
                    'class': 'Kinase',
                    'kegg_id': 'K00850'
                }
            ]
        },
        'pathway': {
            'MAPK': [
                {
                    'id': 'path:hsa04010',
                    'definition': 'MAPK signaling pathway - Homo sapiens (human)',
                    'pathway': 'Signal transduction',
                    'organism': 'Homo sapiens',
                    'class': 'Signaling pathway',
                    'kegg_id': 'hsa04010'
                }
            ],
            'insulin': [
                {
                    'id': 'path:hsa04910',
                    'definition': 'Insulin signaling pathway - Homo sapiens (human)',
                    'pathway': 'Endocrine system',
                    'organism': 'Homo sapiens',
                    'class': 'Metabolic pathway',
                    'kegg_id': 'hsa04910'
                }
            ],
            'glycolysis': [
                {
                    'id': 'path:hsa00010',
                    'definition': 'Glycolysis / Gluconeogenesis - Homo sapiens (human)',
                    'pathway': 'Carbohydrate metabolism',
                    'organism': 'Homo sapiens',
                    'class': 'Metabolic pathway',
                    'kegg_id': 'hsa00010'
                }
            ]
        }
    }

    # Return sample data based on query and database
    query_lower = query.lower()
    if database in sample_data:
        for key, data_list in sample_data[database].items():
            if key.lower() in query_lower:
                return data_list

    # Add more sample data for other databases
    if database == 'enzyme':
        return [
            {
                'id': 'ec:2.7.1.1',
                'definition': 'hexokinase',
                'pathway': 'Glycolysis / Gluconeogenesis',
                'organism': 'KEGG Database',
                'class': 'Transferase',
                'kegg_id': 'K00844'
            },
            {
                'id': 'ec:2.7.1.2',
                'definition': 'glucokinase',
                'pathway': 'Glycolysis / Gluconeogenesis',
                'organism': 'KEGG Database',
                'class': 'Transferase',
                'kegg_id': 'K12407'
            }
        ]
    elif database == 'compound':
        return [
            {
                'id': 'cpd:C00031',
                'definition': 'D-Glucose',
                'pathway': 'Glycolysis / Gluconeogenesis',
                'organism': 'KEGG Database',
                'class': 'Monosaccharide',
                'kegg_id': 'C00031'
            },
            {
                'id': 'cpd:C00024',
                'definition': 'Acetyl-CoA',
                'pathway': 'Citrate cycle (TCA cycle)',
                'organism': 'KEGG Database',
                'class': 'Coenzyme',
                'kegg_id': 'C00024'
            }
        ]
    elif database == 'disease':
        return [
            {
                'id': 'ds:H00031',
                'definition': 'Diabetes mellitus, type 2',
                'pathway': 'Endocrine and metabolic disease',
                'organism': 'KEGG Database',
                'class': 'Disease',
                'kegg_id': 'H00031'
            },
            {
                'id': 'ds:H00001',
                'definition': 'Alzheimer disease',
                'pathway': 'Neurodegenerative disease',
                'organism': 'KEGG Database',
                'class': 'Disease',
                'kegg_id': 'H00001'
            }
        ]

    # Default sample data
    return sample_data.get('genes', {}).get('MAPK', [])
