from Bio import Restriction
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import IUPACData

aa_set = set(IUPACData.protein_letters)
dna_set = set(IUPACData.unambiguous_dna_letters)
rna_set = set(IUPACData.unambiguous_rna_letters)

INDETERMINATE = 'Indeterminate'
INVALID = 'Invalid'
RNA = 'rna'
DNA = 'dna'
AA = 'Amino Acid'


def determine_sequence_type(sequence: str):
    sequence_set = set(sequence)
    is_dna = sequence_set.issubset(dna_set)
    is_rna = sequence_set.issubset(rna_set)
    is_aa = sequence_set.issubset(aa_set)

    true_count = 0
    for val in [is_dna, is_aa, is_rna]:
        if val:
            true_count += 1
    if true_count == 0:
        return INVALID
    if true_count == 1:
        if is_rna:
            return RNA
        if is_aa:
            return AA
        if is_dna:
            return DNA
    return INDETERMINATE


def clean_sequence(sequence: str):
    return ("".join(sequence.split())).upper()


def turn_string_sequence_into_amino(sequence: str, sequence_type: str = ''):
    sequence = clean_sequence(sequence)
    if sequence_type == '':
        sequence_type = determine_sequence_type(sequence)

    # handle invalid sequence type stuff
    if sequence_type == INDETERMINATE or sequence_type == INVALID:
        raise Exception("Can't figure out the sequence type")

    # transform to a protein
    if sequence_type == RNA:
        sequence = Seq(sequence, IUPAC.unambiguous_rna).translate()
    elif sequence_type == DNA:
        sequence = Seq(sequence, IUPAC.unambiguous_dna).translate()
    elif sequence_type == AA:
        sequence = Seq(sequence, IUPAC.protein)

    return sequence


def get_rest_enzymes_from_list(restriction_enzymes: list):
    return Restriction.RestrictionBatch(
        [Restriction.AllEnzymes.get(enz) for enz in restriction_enzymes]
    )


def get_rest_enzymes_from_string(restriction_enzymes: str):
    return get_rest_enzymes_from_list(restriction_enzymes.split())
