from django.test import TestCase
from biobrick_optimization_tool_synthesis.optimization.logic import interpret_sequence


# Create your tests here.
class InterpretSequenceCase(TestCase):

    def test_determine_sequence_type(self):
        valid_protein = 'AL*RLEGIGLN ATQACSSTVG DPKDGSVKVT FLVMISKEYP ASEMIFEKSS KDARIKVFGT ' \
                        'HMVILDLHNA GDGPGVIDAR RLNLPTLLQL IAFVKGYKIT QMCQAYKREE QNTWPSIFDA ' \
                        'TLEAGSRWSK AGEQLVYREA IIKRCAARPE RLDGPHRNKN TEYEHTFDHL SLQVPELGVM ' \
                        'FPNGVVFDLQ ANSSVD*ITLY '
        valid_dna = 'CGGCTTGCCCTTTCTGCCTGTAGATCCATTGGACTGGTGCCAACGCGCAGGCATAGTTCGAGGAGAATTATCCGGGGGCAATGACAACCAGCATCTCGGGTCTTGCCCAACCCGTCTACACGCTGTTATAGCGAATCAGCGGGAACCCGGTGCCACGCGATGGAACGTCCTTAACTCTGGCAGGCAATTAAAGGGAACGTATATATAACGCAAAGAAGCTGGAAAATTGGCGAGAGAATCTTCTCTCTGTCTATCGAAGAATGGCCACGCGGTGGCAACCGTCATGCTAGCGTGCGGGGTACACTTGCTAACCATTTGGGACACGGGACACTCGCTGTTTTCGAAATTACCCTTTATGCGCGGGTATTGAACCACGCTTATGCCCAGCATCGTTACAAGCAGACTCATACTAGATGTATTATGCCCGCCATGCAGACGAAACCAGTCGGAGATTACCGAGCATTCTATCACGTCGGCGACCACTAGTGAGCTACTGGAGCCGAGGGGTAACGTTGATGCCCCTAAGAACCTCTCGGTCGACGCAAGCGATTACACTCCTGTCACATCATAATCGTTTGCTATTCAGGGCTTGACCAACACTGGATTGCTTTTCACTTAAAGTATTATGCACGACAGGGTGCGTGTACCATGTAAACCTGTTATAACTTACCTCAGACTAGTTGGAAGTGTGGCTAGATCTTAGCTTACGTCACTAGAGGGTCCACGTTTAGTTTTTAAGATCCATTGATCTCCTAAACGCTGCAAGATTCGCAACCTGGTATACTTAGCGCTAGGTCCTAGTGCAGCGGGACTTTTTTTCTAAAGTCGTTGAGAGGAGGAGTCGTCAGACCAGATAGCTTTGATGTCCTGATCGGAAGGATCGTTGGCCCCCGACCCTTAGACTCTGTACTCAGTTCTATAAACGAGCCATTGGATACGAGATCCGTAGATTGATAAGGGACACGGAATATCCCCGGACGCAATAGACGGACAGCTTGGT'
        valid_rna = valid_dna.replace('T', 'U')

        invalid_protein = valid_protein.replace('A', '=', 3)
        invalid_dna = valid_dna.replace('T', '=', 3)
        invalid_rna = valid_rna.replace('U', '=', 3)

        unknown = valid_dna.replace('T', 'C')

        self.assertEqual(interpret_sequence.AA, interpret_sequence.determine_sequence_type(valid_protein),
                         'Determine valid protein sequence')
        self.assertEqual(interpret_sequence.DNA, interpret_sequence.determine_sequence_type(valid_dna),
                         'Determine valid DNA sequence')
        self.assertEqual(interpret_sequence.RNA, interpret_sequence.determine_sequence_type(valid_rna),
                         'Determine valid RNA sequence')
        self.assertEqual(interpret_sequence.INVALID, interpret_sequence.determine_sequence_type(invalid_protein),
                         'Determine invalid Protein sequence')
        self.assertEqual(interpret_sequence.INVALID, interpret_sequence.determine_sequence_type(invalid_dna),
                         'Determine invalid DNA sequence')
        self.assertEqual(interpret_sequence.INVALID, interpret_sequence.determine_sequence_type(invalid_rna),
                         'Determine invalid RNA sequence')
        self.assertEqual(interpret_sequence.INDETERMINATE, interpret_sequence.determine_sequence_type(unknown),
                         'Determine unknown sequence')
