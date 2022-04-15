import unittest

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

from src.codon_optimization_host import *


class COH(unittest.TestCase):
    def test_determine_ideal_codon_optimized_sequence_pass(self):
        sequence = Seq('ATRPKKTT', IUPAC.protein)
        usage = {'A': {'aaa': 0.9, 'bbb': 0.1}, 'T': {'ttt': 0.2, 'bbb': 0.1}, 'R': {'rrr': 0.5, 'bbb': 0.5},
                 'P': {'ppp': 0.9, 'bbb': 0.1}, 'K': {'kkk': 0.9, 'bbb': 0.1}}
        expected = 'aaatttrrrpppkkkkkktttttt'
        actual = determine_ideal_codon_optimized_sequence(sequence, usage)
        self.assertEqual(expected, actual, 'failed to determine the correct optimized sequence')

    def test_determine_ideal_codon_optimized_sequence_fail(self):
        sequence = Seq('ATRPKKTT', IUPAC.protein)
        # wrong keys
        usage = {}
        self.assertRaises(KeyError, determine_ideal_codon_optimized_sequence, sequence, usage)
        # no value
        usage = {'A': {}, 'T': {'ttt': 0.2, 'bbb': 0.1}, 'R': {'rrr': 0.5, 'bbb': 0.5},
                 'P': {'ppp': 0.9, 'bbb': 0.1}, 'K': {'kkk': 0.9, 'bbb': 0.1}}
        self.assertRaises(ValueError, determine_ideal_codon_optimized_sequence, sequence, usage)
        # bad numerical value
        # usage = {'A': {'aaa': 'REEE', 'bbb': 'ARRR'}, 'T': {'ttt': 0.2, 'bbb': 0.1}, 'R': {'rrr': 0.5, 'bbb': 0.5},
        #          'P': {'ppp': 0.9, 'bbb': 0.1}, 'K': {'kkk': 0.9, 'bbb': 0.1}}
        # self.assertRaises(ValueError, determine_ideal_codon_optimized_sequence, sequence, usage)
        # bad alphabet
        sequence = Seq('ACTGTACGTACG', IUPAC.unambiguous_dna)
        usage = {'A': {'aaa': 0.9, 'bbb': 0.1}, 'T': {'ttt': 0.2, 'bbb': 0.1}, 'R': {'rrr': 0.5, 'bbb': 0.5},
                 'P': {'ppp': 0.9, 'bbb': 0.1}, 'K': {'kkk': 0.9, 'bbb': 0.1}}
        self.assertRaises(ValueError, determine_ideal_codon_optimized_sequence, sequence, usage)


if __name__ == '__main__':
    unittest.main()
