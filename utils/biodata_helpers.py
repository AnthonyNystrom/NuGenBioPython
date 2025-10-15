"""
Helper functions for Bio.Data operations
Handles CodonTable, IUPACData, and PDBData functionality
"""
from Bio.Data import CodonTable, IUPACData, PDBData
from Bio.Seq import Seq


def get_codon_tables():
    """Get all available codon tables"""
    tables = {}
    for table_id in range(1, 34):
        try:
            table = CodonTable.unambiguous_dna_by_id[table_id]
            tables[table_id] = {
                'id': table_id,
                'name': table.names[0] if table.names else f"Table {table_id}",
                'start_codons': list(table.start_codons),
                'stop_codons': list(table.stop_codons)
            }
        except KeyError:
            continue
    return tables


def translate_sequence(sequence, table_id=1):
    """Translate DNA sequence using specified codon table"""
    try:
        seq_obj = Seq(sequence)
        translation = seq_obj.translate(table=table_id)
        table = CodonTable.unambiguous_dna_by_id[table_id]

        return {
            'success': True,
            'original_sequence': sequence,
            'translated_sequence': str(translation),
            'table_id': table_id,
            'table_name': table.names[0] if table.names else f"Table {table_id}",
            'start_codons': list(table.start_codons),
            'stop_codons': list(table.stop_codons),
            'length_original': len(sequence),
            'length_translated': len(translation)
        }
    except Exception as e:
        return {
            'success': False,
            'error': str(e)
        }


def get_iupac_codes(code_type='dna'):
    """Get IUPAC ambiguity codes for DNA or RNA"""
    if code_type == 'dna':
        return {
            'letters': IUPACData.ambiguous_dna_letters,
            'values': IUPACData.ambiguous_dna_values,
            'complement': IUPACData.ambiguous_dna_complement
        }
    elif code_type == 'rna':
        return {
            'letters': IUPACData.ambiguous_rna_letters,
            'values': IUPACData.ambiguous_rna_values,
            'complement': IUPACData.ambiguous_rna_complement
        }
    elif code_type == 'protein':
        return {
            'letters': IUPACData.protein_letters,
            'extended_letters': IUPACData.extended_protein_letters,
            'one_to_three': IUPACData.protein_letters_1to3,
            'three_to_one': IUPACData.protein_letters_3to1,
            'extended_one_to_three': IUPACData.protein_letters_1to3_extended,
            'extended_three_to_one': IUPACData.protein_letters_3to1_extended
        }
    return {}


def convert_protein_letters(input_str, conversion_type='1to3'):
    """Convert protein letters between 1-letter and 3-letter codes"""
    try:
        if conversion_type == '1to3':
            # Convert 1-letter to 3-letter
            result = []
            for char in input_str.upper():
                if char in IUPACData.protein_letters_1to3:
                    result.append(IUPACData.protein_letters_1to3[char])
                elif char in IUPACData.protein_letters_1to3_extended:
                    result.append(IUPACData.protein_letters_1to3_extended[char])
                else:
                    result.append('???')
            return {'success': True, 'result': '-'.join(result), 'count': len(result)}
        else:  # 3to1
            # Convert 3-letter to 1-letter
            codes = [c.strip() for c in input_str.upper().replace('-', ' ').split()]
            result = []
            for code in codes:
                if code in IUPACData.protein_letters_3to1:
                    result.append(IUPACData.protein_letters_3to1[code])
                elif code in IUPACData.protein_letters_3to1_extended:
                    result.append(IUPACData.protein_letters_3to1_extended[code])
                else:
                    result.append('X')
            return {'success': True, 'result': ''.join(result), 'count': len(result)}
    except Exception as e:
        return {'success': False, 'error': str(e)}


def calculate_molecular_weight(sequence, seq_type='protein', weight_type='average'):
    """Calculate molecular weight of a sequence"""
    try:
        sequence = sequence.upper().replace(' ', '').replace('\n', '')
        total_weight = 0.0
        composition = {}

        if seq_type == 'protein':
            if weight_type == 'monoisotopic':
                weights = IUPACData.monoisotopic_protein_weights
            else:
                weights = IUPACData.protein_weights

            for aa in sequence:
                if aa in weights:
                    total_weight += weights[aa]
                    composition[aa] = composition.get(aa, 0) + 1

        elif seq_type == 'dna':
            if weight_type == 'monoisotopic':
                weights = IUPACData.monoisotopic_unambiguous_dna_weights
            else:
                weights = IUPACData.unambiguous_dna_weights

            for base in sequence:
                if base in weights:
                    total_weight += weights[base]
                    composition[base] = composition.get(base, 0) + 1

        elif seq_type == 'rna':
            if weight_type == 'monoisotopic':
                weights = IUPACData.monoisotopic_unambiguous_rna_weights
            else:
                weights = IUPACData.unambiguous_rna_weights

            for base in sequence:
                if base in weights:
                    total_weight += weights[base]
                    composition[base] = composition.get(base, 0) + 1

        return {
            'success': True,
            'weight': round(total_weight, 2),
            'length': len(sequence),
            'composition': composition,
            'weight_type': weight_type
        }
    except Exception as e:
        return {'success': False, 'error': str(e)}


def get_pdb_conversions(conversion_type='protein_1to3'):
    """Get PDB residue name conversions"""
    if conversion_type == 'protein_1to3':
        return dict(PDBData.protein_letters_1to3)
    elif conversion_type == 'protein_3to1':
        return dict(PDBData.protein_letters_3to1)
    elif conversion_type == 'protein_3to1_extended':
        return dict(PDBData.protein_letters_3to1_extended)
    elif conversion_type == 'nucleic_3to1':
        return dict(PDBData.nucleic_letters_3to1)
    elif conversion_type == 'nucleic_3to1_extended':
        return dict(PDBData.nucleic_letters_3to1_extended)
    return {}


def get_atom_weights():
    """Get atomic weights for common elements"""
    return dict(IUPACData.atom_weights)


def lookup_iupac_code(code, code_type='dna'):
    """Look up what bases a specific IUPAC code represents"""
    try:
        code = code.upper()
        if code_type == 'dna':
            if code in IUPACData.ambiguous_dna_values:
                return {
                    'success': True,
                    'code': code,
                    'bases': IUPACData.ambiguous_dna_values[code],
                    'complement': IUPACData.ambiguous_dna_complement.get(code, 'N')
                }
        elif code_type == 'rna':
            if code in IUPACData.ambiguous_rna_values:
                return {
                    'success': True,
                    'code': code,
                    'bases': IUPACData.ambiguous_rna_values[code],
                    'complement': IUPACData.ambiguous_rna_complement.get(code, 'N')
                }
        return {'success': False, 'error': f'Code {code} not found in {code_type.upper()} IUPAC codes'}
    except Exception as e:
        return {'success': False, 'error': str(e)}
