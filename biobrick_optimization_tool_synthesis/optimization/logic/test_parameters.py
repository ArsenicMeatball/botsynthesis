from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

from biobrick_optimization_tool_synthesis.optimization.logic import interpret_sequence
from biobrick_optimization_tool_synthesis.optimization.logic.codon_optimization_host import \
    determine_ideal_codon_optimized_sequence
from biobrick_optimization_tool_synthesis.optimization.logic.set_codon_table import fetch_codon_table

valid_protein = 'AL*RLEGIGLN ATQACSSTVG DPKDGSVKVT FLVMISKEYP ASEMIFEKSS KDARIKVFGT ' \
                'HMVILDLHNA GDGPGVIDAR RLNLPTLLQL IAFVKGYKIT QMCQAYKREE QNTWPSIFDA ' \
                'TLEAGSRWSK AGEQLVYREA IIKRCAARPE RLDGPHRNKN TEYEHTFDHL SLQVPELGVM ' \
                'FPNGVVFDLQ ANSSVD*ITLY '
valid_protein2 = 'AL*RLEGIGLN '
valid_dna = 'CGGCTTGCCCTTTCTGCCTGTAGATCCATTGGACTGGTGCCAACGCGCAGGCATAGTTCGAGGAGAATTATCCGGGGGCAATGACAACCAGCATCTCGGGTCTTGCCCAACCCGTCTACACGCTGTTATAGCGAATCAGCGGGAACCCGGTGCCACGCGATGGAACGTCCTTAACTCTGGCAGGCAATTAAAGGGAACGTATATATAACGCAAAGAAGCTGGAAAATTGGCGAGAGAATCTTCTCTCTGTCTATCGAAGAATGGCCACGCGGTGGCAACCGTCATGCTAGCGTGCGGGGTACACTTGCTAACCATTTGGGACACGGGACACTCGCTGTTTTCGAAATTACCCTTTATGCGCGGGTATTGAACCACGCTTATGCCCAGCATCGTTACAAGCAGACTCATACTAGATGTATTATGCCCGCCATGCAGACGAAACCAGTCGGAGATTACCGAGCATTCTATCACGTCGGCGACCACTAGTGAGCTACTGGAGCCGAGGGGTAACGTTGATGCCCCTAAGAACCTCTCGGTCGACGCAAGCGATTACACTCCTGTCACATCATAATCGTTTGCTATTCAGGGCTTGACCAACACTGGATTGCTTTTCACTTAAAGTATTATGCACGACAGGGTGCGTGTACCATGTAAACCTGTTATAACTTACCTCAGACTAGTTGGAAGTGTGGCTAGATCTTAGCTTACGTCACTAGAGGGTCCACGTTTAGTTTTTAAGATCCATTGATCTCCTAAACGCTGCAAGATTCGCAACCTGGTATACTTAGCGCTAGGTCCTAGTGCAGCGGGACTTTTTTTCTAAAGTCGTTGAGAGGAGGAGTCGTCAGACCAGATAGCTTTGATGTCCTGATCGGAAGGATCGTTGGCCCCCGACCCTTAGACTCTGTACTCAGTTCTATAAACGAGCCATTGGATACGAGATCCGTAGATTGATAAGGGACACGGAATATCCCCGGACGCAATAGACGGACAGCTTGGT'
valid_dna2 = 'CGGCTT'
valid_rna = valid_dna.replace('T', 'U')

invalid_protein = valid_protein.replace('A', '=', 3)
invalid_dna = valid_dna.replace('T', '=', 3)
invalid_rna = valid_rna.replace('U', '=', 3)

unknown = valid_dna.replace('T', 'C')

seq = Seq(interpret_sequence.clean_sequence(valid_protein), IUPAC.protein)
codon_usage = fetch_codon_table()
codon_optimized = determine_ideal_codon_optimized_sequence(seq, codon_usage)
algorithm_params = {
    'population_size': 5,
    'num_param_considered': 8,
    'prot_seq': seq,
    'codon_opt_seq': codon_optimized,
    'codon_usage': codon_usage,
    'mutation_%': 5,
    'archive_size': 2,
    'crossover_%': 15,
    'gc_parameters': {
        'min': 0.1,
        'max': 0.9,
        'window_size': 20
    },
    'restriction_sites': 'EcoRI'
}