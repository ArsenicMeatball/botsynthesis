import unittest

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

from biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.mutations import *


class TestMutations(unittest.TestCase):
    def test_mutate_codon(self):
        self.assertEqual(True, False)

    def test_mutate_seq(self):
        self.assertEqual(True, False)

    def test_initialize_population(self):
        self.assertEqual(True, False)

    def test_tournament_selection_without_replacement(self):
        key = 'key'
        population = {
            1: {key: 9099990.393}, 2: {key: 9099990.392}, 3: {key: 9099990.394}
        }
        self.assertRaises(ValueError, tournament_selection_without_replacement, population, 1)
        winner = tournament_selection_without_replacement(population, fitness_key_name=key)
        self.assertIn(winner, population.keys(), 'winner should be in the population')
        winner = tournament_selection_without_replacement(population, n_ary=100, fitness_key_name=key)
        expected = 2
        self.assertEqual(expected, winner, 'winner should be 2, (2/3)^200 chance of not occurring')
        winner = tournament_selection_without_replacement(population, n_ary=100, minima=False, fitness_key_name=key)
        expected = 3
        self.assertEqual(expected, winner, 'winner should be 3, (2/3)^200 chance of not occurring')

    def test_generate_mating_pool_from_archive(self):
        self.assertEqual(True, False)

    def test_recombine(self):
        # test that from 2 strings, the new one is created
        # new one is same aa
        s1 = 'actagctacggaattagcagaagcaatgctagc'
        s2 = 'acaagttatgggatctcgcgatcgaacgcaagt'
        aa1 = Seq(s1, IUPAC.unambiguous_dna).translate()
        aa2 = Seq(s2, IUPAC.unambiguous_dna).translate()
        self.assertEqual(aa2, aa1, 'The setup for recombine is bad, the test cannot be ran')
        for num in range(1, 20):
            r = recombine_dna_sequence(s1, s2, num)
            raa = Seq(r, IUPAC.unambiguous_dna).translate()
            self.assertEqual(raa, aa1, 'Recombined sequence as AA is not the same AA sequence')
            for pos in range(len(r)):
                self.assertTrue(r[pos] == s1[pos] or r[pos] == s2[pos], 'the sequence is not composed of its parents')
            if num > 10:
                self.assertTrue(r != s1, 'at large numbers of sites, the string should theoretically nearly never be '
                                         'the same as the parent chance=1/11^10')
        self.assertRaises(ValueError, recombine_dna_sequence, s1, s2, 0)
        self.assertRaises(ValueError, recombine_dna_sequence, s1, s2[:-2], 5)

    def test_generate_population_from_archive(self):
        self.assertEqual(True, False)


if __name__ == '__main__':
    unittest.main()
