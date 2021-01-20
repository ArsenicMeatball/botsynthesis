import unittest

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.Data import CodonTable

import biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.mutations as mut
import biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.dict_functions as dictf


# TODO finish tests
class TestMutations(unittest.TestCase):
    def test_mutate_codon_pass(self):
        s1 = 'actagctacggaattagcagaagcaatgctagc'
        s2 = 'acaagttatgggatctcgcgatcgaacgcaagt'
        # regular case
        initial_codon = 'aac'.upper()
        codon_table = CodonTable.generic_by_id[1]
        new_codon = mut.mutate_codon(initial_codon, codon_table)
        self.assertNotEqual(initial_codon, new_codon, 'mutate codon failed to change the codon when given the '
                                                      'opportunity')
        possible_codons = dictf.invert_dict(codon_table.forward_table)[codon_table.forward_table[initial_codon]]
        self.assertIn(new_codon, possible_codons, 'codon did not become a possible codon - wtf happened')
        # works with another translation dict
        codon_table = CodonTable.generic_by_id[2]
        new_codon = mut.mutate_codon(initial_codon, codon_table)
        self.assertNotEqual(initial_codon, new_codon, 'mutate codon failed to change the codon when given the '
                                                      'opportunity - using different table')
        possible_codons = dictf.invert_dict(codon_table.forward_table)[codon_table.forward_table[initial_codon]]
        self.assertIn(new_codon, possible_codons, 'codon did not become a possible codon - wtf happened - '
                                                  'using different table')
        # no other codons case
        initial_codon = 'ttg'.upper()
        codon_table = CodonTable.generic_by_id[1]
        new_codon = mut.mutate_codon(initial_codon, codon_table)
        self.assertEqual(initial_codon, new_codon, 'mutate codon failed to keep the codon when given only 1 option')

    def test_mutate_codon_fail(self):
        # smol codon
        initial_codon = 'ac'
        codon_table = 'not a codon tbale object'
        self.assertRaises(ValueError, mut.mutate_codon, initial_codon, codon_table)
        # bad codon usage
        initial_codon = 'acg'
        self.assertRaises(AttributeError, mut.mutate_codon, initial_codon, codon_table)


    def test_mutate_seq(self):
        self.assertEqual(True, False)

    def test_initialize_population(self):
        self.assertEqual(True, False)

    def test_tournament_selection_without_replacement(self):
        key = 'key'
        population = {
            1: {key: 9099990.393}, 2: {key: 9099990.392}, 3: {key: 9099990.394}
        }
        self.assertRaises(ValueError, mut.tournament_selection_without_replacement, population, 1)
        winner = mut.tournament_selection_without_replacement(population, fitness_key_name=key)
        self.assertIn(winner, population.keys(), 'winner should be in the population')
        winner = mut.tournament_selection_without_replacement(population, n_ary=100, fitness_key_name=key)
        expected = 2
        self.assertEqual(expected, winner, 'winner should be 2, (2/3)^200 chance of not occurring')
        winner = mut.tournament_selection_without_replacement(population, n_ary=100, minima=False, fitness_key_name=key)
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
            r = mut.recombine_dna_sequence(s1, s2, num)
            raa = Seq(r, IUPAC.unambiguous_dna).translate()
            self.assertEqual(raa, aa1, 'Recombined sequence as AA is not the same AA sequence')
            for pos in range(len(r)):
                self.assertTrue(r[pos] == s1[pos] or r[pos] == s2[pos], 'the sequence is not composed of its parents')
            if num > 10:
                self.assertTrue(r != s1, 'at large numbers of sites, the string should theoretically nearly never be '
                                         'the same as the parent chance=1/11^10')
        self.assertRaises(ValueError, mut.recombine_dna_sequence, s1, s2, 0)
        self.assertRaises(ValueError, mut.recombine_dna_sequence, s1, s2[:-2], 5)

    def test_generate_population_from_archive(self):
        self.assertEqual(True, False)


if __name__ == '__main__':
    unittest.main()
