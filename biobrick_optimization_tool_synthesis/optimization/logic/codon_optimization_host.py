
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from biobrick_optimization_tool_synthesis.optimization.logic import test_sequences
from biobrick_optimization_tool_synthesis.optimization.logic.set_codon_table import fetch_codon_table

def determine_ideal_codon_optimized_sequence(sequence: Seq, codon_usage: dict) -> Seq:
    # requires amino acid sequence and a dict mapping amino acids to codons

    if sequence.alphabet is not IUPAC.protein:
        raise Exception("It ain't right alphabet, mus tbe protein")
    for amino_acid in sequence:



    return sequence


if __name__ == '__main__':
    sequence = Seq(test_sequences.valid_protein, IUPAC.protein)
    codon_usage = fetch_codon_table()
    determine_ideal_codon_optimized_sequence(sequence, codon_usage)