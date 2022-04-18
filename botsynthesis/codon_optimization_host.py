from Bio.Seq import Seq


def determine_ideal_codon_optimized_sequence(
    sequence: Seq, codon_usage: dict
) -> Seq:
    """
        Determine the codon optimized sequence that would be ideal in increase
        expression in the host
    :param sequence: (Bio.Seq.Seq) ~ (.alphabet == IUPAC.protein) sequence
        representing the amino acids to be synthesized
    :param codon_usage: (dict) represents the codon frequencies for each amino
        acid for the given host
    :return: (Bio.Seq.Seq) sequence representing the optimal codon sequence for
        expression
    """
    # requires amino acid sequence and a dict mapping amino acids to codons
    codon_optimized_sequence = ""
    for amino_acid in sequence:
        codon_optimized_sequence += max(
            codon_usage[amino_acid],
            key=lambda key: codon_usage[amino_acid][key],
        )
    return Seq(codon_optimized_sequence)


if __name__ == "__main__":
    pass
    # sequence = Seq(cleaning.clean_sequence(test_parameters.valid_protein),
    # IUPAC.protein)
    # codon_usage = fetch_codon_table()
    # codon_optimized = determine_ideal_codon_optimized_sequence(sequence,
    # codon_usage)
    # print(codon_usage)
    # print(codon_optimized)
