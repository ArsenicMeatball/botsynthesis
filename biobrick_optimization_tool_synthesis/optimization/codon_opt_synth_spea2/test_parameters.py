from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

from biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2 import cleaning
from biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.cleaning import \
    get_rest_enzymes_from_string
from biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.codon_optimization_host import \
    determine_ideal_codon_optimized_sequence
from biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.set_codon_table import fetch_codon_table

valid_protein = 'AL*RLEGIGLN ATQACSSTVG DPKDGSVKVT FLVMISKEYP ASEMIFEKSS KDARIKVFGT ' \
                'HMVILDLHNA GDGPGVIDAR RLNLPTLLQL IAFVKGYKIT QMCQAYKREE QNTWPSIFDA ' \
                'TLEAGSRWSK AGEQLVYREA IIKRCAARPE RLDGPHRNKN TEYEHTFDHL SLQVPELGVM ' \
                'FPNGVVFDLQ ANSSVD*ITLY'
valid_protein2 = 'AL*RLEGIGLN '
valid_dna = 'CGGCTTGCCCTTTCTGCCTGTAGATCCATTGGACTGGTGCCAACGCGCAGGCATAGTTCGAGGAGAATTATCCGGGGGCAATGACAACCAGCATCTCGGGTCTTGCCCAACCCGTCTACACGCTGTTATAGCGAATCAGCGGGAACCCGGTGCCACGCGATGGAACGTCCTTAACTCTGGCAGGCAATTAAAGGGAACGTATATATAACGCAAAGAAGCTGGAAAATTGGCGAGAGAATCTTCTCTCTGTCTATCGAAGAATGGCCACGCGGTGGCAACCGTCATGCTAGCGTGCGGGGTACACTTGCTAACCATTTGGGACACGGGACACTCGCTGTTTTCGAAATTACCCTTTATGCGCGGGTATTGAACCACGCTTATGCCCAGCATCGTTACAAGCAGACTCATACTAGATGTATTATGCCCGCCATGCAGACGAAACCAGTCGGAGATTACCGAGCATTCTATCACGTCGGCGACCACTAGTGAGCTACTGGAGCCGAGGGGTAACGTTGATGCCCCTAAGAACCTCTCGGTCGACGCAAGCGATTACACTCCTGTCACATCATAATCGTTTGCTATTCAGGGCTTGACCAACACTGGATTGCTTTTCACTTAAAGTATTATGCACGACAGGGTGCGTGTACCATGTAAACCTGTTATAACTTACCTCAGACTAGTTGGAAGTGTGGCTAGATCTTAGCTTACGTCACTAGAGGGTCCACGTTTAGTTTTTAAGATCCATTGATCTCCTAAACGCTGCAAGATTCGCAACCTGGTATACTTAGCGCTAGGTCCTAGTGCAGCGGGACTTTTTTTCTAAAGTCGTTGAGAGGAGGAGTCGTCAGACCAGATAGCTTTGATGTCCTGATCGGAAGGATCGTTGGCCCCCGACCCTTAGACTCTGTACTCAGTTCTATAAACGAGCCATTGGATACGAGATCCGTAGATTGATAAGGGACACGGAATATCCCCGGACGCAATAGACGGACAGCTTGGTAGCTTGGTTCGACA'
valid_dna2 = 'CGGCTT'
valid_rna = valid_dna.replace('T', 'U')

invalid_protein = valid_protein.replace('A', '=', 3)
invalid_dna = valid_dna.replace('T', '=', 3)
invalid_rna = valid_rna.replace('U', '=', 3)

unknown = valid_dna.replace('T', 'C')

seq = Seq(cleaning.clean_sequence(valid_protein), IUPAC.protein)
codon_usage = fetch_codon_table()
codon_optimized = determine_ideal_codon_optimized_sequence(seq, codon_usage)
restriction_sites = 'NdeI XhoI HpaI PstI EcoRV NcoI BamHI'
algorithm_params = {
    'population_size': 5,
    'max init population attempts': 25,
    'prot_seq': seq,
    'codon_opt_seq': codon_optimized,
    'codon_usage': codon_usage,
    'mutation_%': 5,
    'archive size': 2,
    'crossover_%': 15,
    'linear': True,
    'gc_parameters': {
        'min': 0.1,
        'max': 0.68,
        'window_size': 20
    },
    'generations': 1,
    'restriction_sites': get_rest_enzymes_from_string(restriction_sites),
    'overlapping': True,
    'locations': True,
    'repeat_size': 10,
    'homopolymer_size': 4,
    'shortest_loop_length': 3,
    'longest_loop_length': 9,
    'stem_length': 10,

}
